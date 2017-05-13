#!/usr/bin/env python

import gc
import os
import re
import sys
import json
import time
import numpy
import struct
import bisect
import tarfile
import tempfile
import argparse
import itertools
import subprocess
import collections
import cffi
ffi = cffi.FFI()

parser = argparse.ArgumentParser(     description="This tool converys FASTQ files to microq (uQ) files and back.")
parser.add_argument("-i", "--input",  required=True,                      help='Required. Input FASTQ/uQ file path.')
parser.add_argument("-o", "--output",                                     help='Optional. FASTQ->uQ only. Default is to append .uQ to input filename.')
parser.add_argument("--compressor",   action='store',                     help='Optional. Path to compression program that accepts data on stdin and prints to stdout/pipe)')
parser.add_argument("--sort",         action='store',                     help='Optional. [DNA/QUAL/QNAME/None] Resort FASTQ. See output of --test for optimium method.')
parser.add_argument("--raw",          nargs='+', metavar='file name',     help='Optional. [DNA/QUAL/QNAME/None] Default is to store DNA/QUAL/QNAMEs in sorted order with a key. Specifying Saves data as a raw table rather than unique & sorted +key combo. [DNA/QUAL/QNAME]')
parser.add_argument("--pattern",      nargs='+', metavar='pattern',       help='Optional. [0.1/0.2/1.1/1.2/2.1/2.2/3.1/3.2] x2 (DNA|QUAL) See output of --test for optimium.')
parser.add_argument("--temp",                                             help='Optional. Directory to write temporary files to. Default is OS dependant.')
parser.add_argument("--test",         action='store_true', default=False, help='Optional. Try all possible sort/raw/pattern combinations (unless user specifies a fixed sort/raw/pattern as well). If --compressor not used, will optimise for decompressed filesize.')
parser.add_argument("--notricks",     action='store_true', default=False, help='Optional. Prevents conversion of N to the most popular base. See documentation for details.')
parser.add_argument("--pad",          action='store_true', default=False, help='Optional. Pads DNA/QUAL to the nearest 2/4/8 bits. Some compressors do a better job when data is padded.')
parser.add_argument("--peek",         action='store_true', default=False, help='Optional. No output files are created, the input is just scanned and report printed to terminal.')
parser.add_argument("--decode",       action='store_true', default=False, help='Requred if you want to convert a .uQ back to .fastq')
args = parser.parse_args()

if args.compressor:
    ## Unfortunately, it seems subprocessing in Linux has an unhelpful requirement that all parent process memory is copied and made avalible to the subprocess - at least as a 
    ## possibility, which in turn means starting a subprocess means doubling your memory space. This is not an issue if doubling memory is OK for you, but for most people operating
    ## at over half system memory, this is unacceptible. The work around is to spawn a bunch of subprocesses (656 to be precise) BEFORE any real memory is used, and thus the subprocesses
    ## will use little system memory. This is an incredibly dumb hack, but unfortunately it's built into Linux and the typical solutions (use vfork or posix_spawn) do not have easy implimentations
    last_subprocess_used = 0
    subprocesses = []
    DEVNULL = open(os.devnull,'w')
    for x in range(656):
        p = subprocess.Popen( args.compressor + ' | wc -c', stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=DEVNULL, shell=True)
        subprocesses.append(p)

def error(message):
    print message
    exit()

if args.pattern:
    if len(args.pattern) != 2: error('ERROR: There must be 2 values for --pattern!')
    elif not all([pattern in ['0.1','1.1','2.1','3.1','0.2','1.2','2.2','3.2'] for pattern in args.pattern]): error('ERROR: Pattern values are incorrect!')
    ## User must supply two patterns. No real reason for this, other than user is unlikely to know one table's pattern without knowing the others...

if args.sort:
    if not args.sort.lower() in ['dna','qual','qname','none']: error('ERROR: --sort value is incorrect!')
    if args.sort.lower() == 'none': args.sort = (None,)
    ## If the user specifies --sort "None", this means do not sort during --test. If user does not specify --sort, will not sort when making uQ without --test, but will try all sorts when run with --test. i.e., "None" enforces no sort.

if args.raw:
    if not all([raw.lower() in ['dna','qual','qname','none'] for raw in args.raw]): error('ERROR: --raw values are incorrect!')
    args.raw = set(args.raw)
    if 'none' in args.raw:
        args.raw.add(None)
        args.raw.discard('none')
    ## As above, if --raw is not set then when making uQs default is to store all tables as sorted/unique and not raw. When running --test default is to try everything.
    ## If user secifies "--raw None" however, then --test will not try different sorts. The difference here is args.raw == None, and args.raw == set([None,])

if not os.path.isfile(args.input):  error('ERROR: Sorry, the input path you have specified is not a file!')

if not args.decode:
    ## File output name:
    if args.output == None: args.output = args.input + '.uQ'

    ## Make a temp directory to work in:
    temp_directory = args.temp if args.temp else tempfile.mkdtemp()
    print 'Temporary directory:',temp_directory

    ## This class prints the status messages
    class init_status:
        def __init__(self,inFile):
            print 'Warming up...'
            self.total = int(subprocess.check_output('wc -l ' + inFile,shell=True).split()[0])
            if self.total % 4 == 0: self.total /= 4
            else: error('ERROR: The FASTQ file provided contains' + str(self.total) +'rows, which is not divisible by 4!')
            self.current = 0
            self.time = time.time()
            self.split = time.time()
            self.message = ''
        def update(self):
            self.current += 1
            if self.current & 65535 == 65535: # an efficient way of saying "do somthing every couple of FASTQ entrys"
                sys.stdout.write('\033[A')
                sys.stdout.flush()
                print self.message, '(' + str( (self.current/float(self.total))*100)[:5] + '%)'
        def split_time(self):
            split = time.time() - self.split
            self.split = time.time()
            return '(' + str(split/60) + ' minutes)'
        def total_time(self): return time.time() - self.time
    status = init_status(args.input)

    ## This function does the actual encoding of the DNA and QUAL values of fixed length. In an idea world the two would be seperated to reduce memory 
    ## consumption and simplify the code, however as large performance benefits exist if "N" bases are stored as a regular base, but with a unique 
    ## quality value, encoding of DNA and QUAL is best done simultaniously for now. Could also write both (+QNAME) straight to disk, then load it back in...?
    def encoder_fixed(total_reads,bases,qualities,N_qual,dna_bytes_per_row,quality_bytes_per_row,status,file_path,bits_per_base,bits_per_quality):
        with open(file_path,'rb') as f:
            status.current = 0
            ffi.cdef("""
                uint8_t *malloc(size_t size);
                uint8_t free(void *ptr);
            """)
            lib = ffi.dlopen(None)

            ## Old way - simple but I can't manually free the memory later, which is a pity as that would be nessecary for other projects.
            #dna_array = ffi.new("uint8_t["+str(total_reads)+"]["+str(dna_bytes_per_row)+"]")
            #qual_array = ffi.new("uint8_t["+str(total_reads)+"]["+str(quality_bytes_per_row)+"]")

            ## New way - more complicated/manual, but allows me to delete the array whenever i want, and not wait around for the Garbage Collector to identify it's not needed any more and remove it.
            dna_pointer = lib.malloc(total_reads*dna_bytes_per_row)
            qual_pointer = lib.malloc(total_reads*quality_bytes_per_row)
            if dna_pointer is None or qual_pointer is None: print "ERROR: Your system does not have the required amount of memory to store all of the DNA/QUALITY scores :("; exit()
            dna_array = ffi.cast("uint8_t["+str(total_reads)+"]["+str(dna_bytes_per_row)+"]", dna_pointer) # Cast the flat array to an array of arrays, also defining its size, type, etc.
            qual_array = ffi.cast("uint8_t["+str(total_reads)+"]["+str(quality_bytes_per_row)+"]", qual_pointer) # Cast the flat array to an array of arrays, also defining its size, type, etc.

            bases = ''.join(bases)
            qualities = ''.join(qualities)
            N_base = 0

            for row in xrange(total_reads):
                try:
                    next(f)
                    dna = next(f)[-2::-1]   # Store sequence backwards, omitting last newline character
                    next(f)
                    quals = next(f)[-2::-1] # Store sequence backwards, omitting last newline character
                    #if 'N' in dna: continue
                except StopIteration: break
                column             = 0
                temp_dna           = 0
                temp_qual          = 0
                dna_bits_done      = 0
                qual_bits_done     = 0
                dna_byte_position  = dna_bytes_per_row -1
                qual_byte_position = quality_bytes_per_row -1
                for base in range(len(dna)):
                    try:
                        temp_dna  += bases.index(dna[base])       << dna_bits_done
                        temp_qual += qualities.index(quals[base]) << qual_bits_done
                    except ValueError:
                        temp_dna  += N_base                       << dna_bits_done
                        temp_qual += N_qual[dna[base]]            << qual_bits_done  # This line is why qualities and DNA are analysed at the same time.
                    dna_bits_done += bits_per_base
                    qual_bits_done += bits_per_quality

                    while dna_bits_done > 8:
                        dna_bits_done -= 8
                        dna_array[row][dna_byte_position] = temp_dna & 255
                        temp_dna >>= 8
                        dna_byte_position -= 1

                    while qual_bits_done > 8:
                        qual_bits_done -= 8
                        qual_array[row][qual_byte_position] = temp_qual & 255
                        temp_qual >>= 8
                        qual_byte_position -= 1

                # Add the last byte of data:
                dna_array[row][dna_byte_position]   = temp_dna  ; dna_byte_position -= 1
                qual_array[row][qual_byte_position] = temp_qual ; qual_byte_position -= 1
                # As dna_byte_position/qual_byte_position may not be zero, memory not written to needs to be zero'd out (I use malloc not calloc):
                while dna_byte_position  >= 1: dna_byte_position  -= 1; dna_array[row][dna_byte_position]   = 0
                while qual_byte_position >= 1: qual_byte_position -= 1; qual_array[row][qual_byte_position] = 0
                status.update()

        ## We send back both the numpy representation of the array, and the raw pointer (to delete the array manually if we wish), and the functions to free it.
        dna_array = numpy.asarray( ffi.buffer(dna_array))
        dna_array.shape = (total_reads,dna_bytes_per_row)
        qual_array = numpy.asarray( ffi.buffer(qual_array))
        qual_array.shape = (total_reads,quality_bytes_per_row)
        return ((dna_array,dna_pointer),(qual_array,qual_pointer),lib)

    ## This function is exactly the same as the above, except comments have been removed and a small block of code at the end of every row has been added.
    ## This was done as it results in a small speed improvement for fixed_length encodings not having that block around wasting cycles.
    ## This adds a '1' to the most significant bit of the binary encoding, telling the decoder to ignore all bits before the 1, and the 1 itself. But ofcourse
    ## It takes 1 extra bit of binary to store variable length things in this way.
    def encoder_variable(total_reads,bases,qualities,N_qual,dna_bytes_per_row,quality_bytes_per_row,status,file_path,bits_per_base,bits_per_quality):
        with open(file_path,'rb') as f:
            status.current = 0
            ffi = cffi.FFI()
            ffi.cdef("""
                uint8_t *malloc(size_t size);
                uint8_t free(void *ptr);
            """)
            lib=ffi.dlopen(None)
            dna_pointer = lib.malloc(total_reads*dna_bytes_per_row)
            qual_pointer = lib.malloc(total_reads*quality_bytes_per_row)
            if dna_pointer is None or qual_pointer is None: print "ERROR: Your system does not have the required amount of memory to store all of the DNA/QUALITY scores :("; exit()
            dna_array = ffi.cast("uint8_t["+str(total_reads)+"]["+str(dna_bytes_per_row)+"]", dna_pointer)
            qual_array = ffi.cast("uint8_t["+str(total_reads)+"]["+str(quality_bytes_per_row)+"]", qual_pointer)
            bases = ''.join(bases)
            qualities = ''.join(qualities)
            N_base = 0
            for row in xrange(total_reads):
                try:
                    next(f)
                    dna   = next(f)[-2::-1]
                    next(f)
                    quals = next(f)[-2::-1]
                except StopIteration: break
                column             = 0
                temp_dna           = 0
                temp_qual          = 0
                dna_bits_done      = 0
                qual_bits_done     = 0
                dna_byte_position  = dna_bytes_per_row -1
                qual_byte_position = quality_bytes_per_row -1
                for base in range(len(dna)):
                    try:
                        temp_dna  += bases.index(dna[base])       << dna_bits_done
                        temp_qual += qualities.index(quals[base]) << qual_bits_done
                    except ValueError:
                        temp_dna  += N_base                       << dna_bits_done
                        temp_qual += N_qual[dna[base]]            << qual_bits_done  # This line is why qualities and DNA are analysed at the same time.
                    dna_bits_done += bits_per_base
                    qual_bits_done += bits_per_quality

                    while dna_bits_done  > 8:
                        dna_bits_done -= 8
                        dna_array[row][dna_byte_position] = temp_dna & 255
                        temp_dna >>= 8
                        dna_byte_position -= 1

                    while qual_bits_done  > 8:
                        qual_bits_done -= 8
                        qual_array[row][qual_byte_position] = temp_qual & 255
                        temp_qual >>= 8
                        qual_byte_position -= 1

                # Add the '1' as the most significant bit:
                dna_array[row][dna_byte_position]   = temp_dna  + (variable_read_lengths << (dna_bits_done))  ; dna_byte_position -= 1
                qual_array[row][qual_byte_position] = temp_qual + (variable_read_lengths << (qual_bits_done)) ; qual_byte_position -= 1
                # As dna_byte_position/qual_byte_position may not be zero, memory not written to needs to be zero'd out (I use malloc not calloc):
                while dna_byte_position  >= 1: dna_byte_position  -= 1; dna_array[row][dna_byte_position]   = 0
                while qual_byte_position >= 1: qual_byte_position -= 1; qual_array[row][qual_byte_position] = 0
                status.update()

        ## We send back both the numpy representation of the array, and the raw pointer (to delete the array manually if we wish), and the functions to free it.
        dna_array = numpy.asarray( ffi.buffer(dna_array))
        dna_array.shape = (total_reads,dna_bytes_per_row)
        qual_array = numpy.asarray( ffi.buffer(qual_array))
        qual_array.shape = (total_reads,quality_bytes_per_row)
        return ((dna_array,dna_pointer),(qual_array,qual_pointer),lib)

    ## Depending on the user defined or optimal --pattern for the DNA and QUAL tables, how the bytes of those tables are arranged on disk (the pattern) is determined below:
    def write_pattern(table,filename):
        if args.pattern is None: pattern = '0.1'; args.pattern = ['0.1','0.1']
        elif filename.startswith('DNA'):  pattern = args.pattern[0]
        elif filename.startswith('QUAL'): pattern = args.pattern[1]
        else: error("ERROR: This should never happen!")
        with open(os.path.join(temp_directory,filename),'wb') as f:
            if pattern == '0.1': numpy.save( f, numpy.ascontiguousarray(             table    ))
            if pattern == '1.1': numpy.save( f, numpy.ascontiguousarray( numpy.rot90(table,1) ))
            if pattern == '2.1': numpy.save( f, numpy.ascontiguousarray( numpy.rot90(table,2) ))
            if pattern == '3.1': numpy.save( f, numpy.ascontiguousarray( numpy.rot90(table,3) ))
            if pattern == '0.2': numpy.save( f, numpy.asfortranarray(                table    ))
            if pattern == '1.2': numpy.save( f, numpy.asfortranarray(    numpy.rot90(table,1) ))
            if pattern == '2.2': numpy.save( f, numpy.asfortranarray(    numpy.rot90(table,2) ))
            if pattern == '3.2': numpy.save( f, numpy.asfortranarray(    numpy.rot90(table,3) ))

    ## Write out without the ".npy" file name extension at the end.
    def write_out(table,filename):
        with open(os.path.join(temp_directory,filename),'wb') as f: numpy.save( f, table )

    # Tests to see if sending only the first X Mb of data doesn't work. The results are way off compared to sending the whole file for test compression.
    def compressed_size(table):
        if not args.compressor: return table.nbytes
        global last_subprocess_used
        p = subprocesses[last_subprocess_used]
        try: numpy.save(p.stdin,table)
        except IOError: pass
        last_subprocess_used += 1
        print last_subprocess_used
        return int(p.communicate()[0])

    ## How uq.py tests different --patterns to see which compresses the best. 
    ## Note that numpy tables are not always freed from memory when python has stopped using them, as technically they exist outside of python-land. As a result, the below can result in
    ## 8x the memory usage per test. This can be fixed by making/cleaning up C arrays manually, which is something i intend to do as it will make encoding literally 1000x faster.
    def test_patterns(table,filename):
        if args.compressor is None: return [table.nbytes, '0.1'] # If the user isn't using compression, patterning won't make a difference, but keys might, so return size. Default pattern is 0.1 for no particular reason. 
        if filename.startswith('DNA'):    pattern = None if args.pattern is None else args.pattern[0] 
        elif filename.startswith('QUAL'): pattern = None if args.pattern is None else args.pattern[1] 
        else: error("ERROR: This should never happen!")

        ## This is the same as the below, but without all the gc.collect()s. Perhaps without the gc.collect is better?
        #results = []
        #if pattern in ('0.1',None): results.append([ compressed_size(numpy.ascontiguousarray(table)               ), '0.1' ])
        #if pattern in ('2.1',None): results.append([ compressed_size(numpy.ascontiguousarray(numpy.rot90(table,2))), '2.1' ])
        #if pattern in ('0.2',None): results.append([ compressed_size(numpy.asfortranarray(table)                  ), '0.2' ])
        #if pattern in ('2.2',None): results.append([ compressed_size(numpy.asfortranarray(numpy.rot90(table,2))   ), '2.2' ])
        #if pattern in ('1.1',None): results.append([ compressed_size(numpy.ascontiguousarray(numpy.rot90(table,1))), '1.1' ]) # This
        #if pattern in ('3.1',None): results.append([ compressed_size(numpy.ascontiguousarray(numpy.rot90(table,3))), '3.1' ]) # This
        #if pattern in ('1.2',None): results.append([ compressed_size(numpy.asfortranarray(numpy.rot90(table,1))   ), '1.2' ]) # This
        #if pattern in ('3.2',None): results.append([ compressed_size(numpy.asfortranarray(numpy.rot90(table,3))   ), '3.2' ]) # This
        #return sorted(results)[0]

        results = []
        if pattern in ('0.1','0.2',None):
            if pattern in ('0.1',None): table = numpy.ascontiguousarray(table) ; gc.collect() ; results.append([compressed_size(table),'0.1'])
            if pattern in ('0.2',None): table = numpy.asfortranarray(table)    ; gc.collect() ; results.append([compressed_size(table),'0.2'])

        if pattern in ('1.1','1.2',None):
            table = numpy.rot90(table,1) ; gc.collect()

            if pattern in ('1.2',None): table = numpy.asfortranarray(table)    ; gc.collect() ; results.append([compressed_size(table),'1.2'])
            if pattern in ('1.1',None): table = numpy.ascontiguousarray(table) ; gc.collect() ; results.append([compressed_size(table),'1.1'])

        if pattern in ('2.1','2.2',None):
            if pattern is None: table = numpy.rot90(table,1) ; gc.collect()
            else:               table = numpy.rot90(table,2) ; gc.collect()

            if pattern in ('2.1',None): table = numpy.ascontiguousarray(table) ; gc.collect() ; results.append([compressed_size(table),'2.1'])
            if pattern in ('2.2',None): table = numpy.asfortranarray(table)    ; gc.collect() ; results.append([compressed_size(table),'2.2'])

        if pattern in ('3.1','3.2',None):
            if pattern is None: table = numpy.rot90(table,1) ; gc.collect()
            else:               table = numpy.rot90(table,3) ; gc.collect()

            if pattern in ('3.2',None): table = numpy.asfortranarray(table)    ; gc.collect() ; results.append([compressed_size(table),'3.2'])
            if pattern in ('3.1',None): table = numpy.ascontiguousarray(table) ; gc.collect() ; results.append([compressed_size(table),'3.1'])

        #for result in sorted(results): print filename,result
        return sorted(results)[0]


    ## OK Now we really get into it:
    status.message = 'Pass 1 of 4: First pass through file to find optimal/safe encoding strategy.'
    with open(args.input,'rb') as f:

        ## Read just the first 4 lines:
        line1,line2,line3,line4 = next(f)[:-1],next(f)[:-1],next(f),next(f)[:-1]
        status.update()

        ## Some very basic FASTQ checks. Used to be more checks, but FASTQ is such a silly format most checks end up being overzealous or true for FASTA too.
        if not line1.startswith('@'): error('ERROR: This does not look like a FASTA/FASTQ file! (first line does not start with @)')

        # QNAME
        prefix = line1
        suffix = line1[:]
        separators = collections.defaultdict(int)
        not_separators = set()

        # DNA
        dna_bases = set(line2)
        dna_min = len(line2)
        dna_max = len(line2)

        # +        This should probably also check that the stuff after the +, if present, is the same as the SEQID...
        if not line3.startswith('+'): error('ERROR: This does not look like a FASTA/FASTQ file! (third line does not start with +)')

        # QUAL
        quals = set(line4)

        ## Check length of DNA == QUAL
        if len(line2) != len(line4): error('ERROR: This does not look like a FASTA/FASTQ file! (SEQ and QUAL lines are not the same length)')

        ## The below generates a count of the qualities per base. This tells us both the counts for each base, but also if any qualities are always associated with a specific base (say, an "N")
        static_qualities = {}
        for base,qual in itertools.izip(line2,line4):
            try:
                static_qualities[base][qual] += 1
            except KeyError:
                static_qualities[base] = collections.defaultdict(int)
                static_qualities[base][qual] += 1

        ## The below essentially feeds the rest of the data in the file into the above constructs. This isn't particularly simple stuff, but working with FASTQ never is...
        while True:
            try:
                qname = next(f)[:-1]
                dna = next(f)[:-1]
                if next(f)[0] != '+': error('ERROR: For entry' + str(status.current) + 'the third line does not start with +')
                qualities = next(f)[:-1]
                status.update()
            except StopIteration: break

            # Quick Check:
            if len(dna) != len(qualities):
                print 'ERROR: Length of DNA does not match the length of the quality scores for entry',status.current
                print '       DNA: ',dna
                print '       QUAL:',qualities
                exit()

            # QNAME Checks:
            if not qname.startswith(prefix):
                for idx,character in enumerate(prefix):
                    if character != qname[idx]:
                        for sep in prefix[idx:]:
                            if sep not in not_separators: separators[sep] += 1
                        prefix = prefix[:idx]
                        break

            if not qname.endswith(suffix):
                for idx,character in enumerate(reversed(suffix)):
                    if character != qname[-1-idx]:
                        if idx == 0: suffix = ''
                        else: suffix = suffix[-idx:]
                        break

            for sep in separators.copy():
                if qname[len(prefix):].count(sep) != separators[sep]:
                    del separators[sep]; gc.collect()
                    not_separators.add(sep)

            # Get DNA max/min:
            if dna_max < len(dna): dna_max = len(dna)
            if dna_min > len(dna): dna_min = len(dna)

            # Get counts of each quality symbol per base:
            for base,qual in itertools.izip(dna,qualities):
                try:
                    static_qualities[base][qual] += 1
                except KeyError:
                    static_qualities[base] = collections.defaultdict(int)
                    static_qualities[base][qual] += 1

    # Tidy up separator/suffix results:
    for sep in separators.copy():
        if suffix.count(sep) != 0:
            separators[sep] -= suffix.count(sep)
            if separators[sep] == 0: del separators[sep]; gc.collect()

    def order_seps(qname,prefix,suffix,separators):
        unordered_seps = ''.join(separators) # concatinates separators
        ordered_seps = re.findall('([' + unordered_seps + ']+)',qname[len(prefix):-1-len(suffix)])
        return ''.join(ordered_seps)

    if order_seps(qname,prefix,suffix,separators) == order_seps(line1,prefix,suffix,separators):
        separators = order_seps(qname,prefix,suffix,separators)
    else:
        print "ERROR: Sorry, the separators used in this file's QNAME/headers are so unusual/improbable that"
        print "       I didn't think it was worth the time to write the code on how to deal with it, only identify it."
        print "       To continue I need to finish encoding QNAMEs as strings, which is on the to-do list, but if you e-mail"
        print "       me about it now i'll get to it sooner rather than later. Sorry about the inconvenience! - john"; exit()

    ## Get data for Base/Quality distribution charts.
    ## Dev note: Could also print a quality distribution per base? Could also order by ASCII code and not most frequent quality? Could scale to largest count rather all counts?) Really this is all a job for a SeQC view.
    base_graph = collections.defaultdict(int)
    qual_graph = collections.defaultdict(int)
    for base,counts_per_qual in static_qualities.items():
        base_graph[base] = sum(counts_per_qual.values())
        for qual,count in counts_per_qual.items(): qual_graph[qual] += count
    total_bases = float(sum(base_graph.values()))
    #dna_bases = map(lambda x: x[0],sorted(base_graph.items(), key=lambda v: v[1], reverse=True)) # Sorting the bases
    #quals     = map(lambda x: x[0],sorted(qual_graph.items(), key=lambda v: v[1], reverse=True)) # and quals by frequency
    dna_bases = sorted(base_graph.keys())
    quals     = sorted(qual_graph.keys())

    ## Print statistics:
    print 'Finished analysing file!',status.split_time()
    print '    - total reads:',status.current
    print '    - total bases:',int(total_bases)
    print ''
    print 'QNAME Analysis:'
    if len(prefix):     print '    - all QNAMEs prefixed with:',prefix
    if len(suffix):     print '    - all QNAMEs suffixed with:',suffix
    if len(separators): print '    - QNAME field separators:'  ,separators,'(',len(separators),'symbols )'
    print ''
    print 'DNA Analysis:'
    print '    - DNA sequences contained the following characters:',' '.join(sorted(dna_bases)),'(' + str(len(dna_bases)) + ' in total)'
    print '    - There distribution is:'
    for base in dna_bases:
        percentage = (base_graph[base]/total_bases)*100
        print '     ',base,'|'+('#'*int(percentage/2)).ljust(50)+'|', ('%.3f'%percentage).rjust(7)+'%', str(base_graph[base]).rjust(12)

    # Decide if special N encoding will be used.
    N_qual = {}
    total_quals = len(quals)
    if args.notricks is False:
        for base,result in static_qualities.items():
            if len(dna_bases) == 1: continue # It could be possible that all bases have unique qualities, and thus no DNA is encoded at all as qualities are enough. But for now, will require at least 1 base of DNA, else code needs refactoring.
            if len(result) == 1:
                dna_bases.remove(base)
                for quality,count in result.items():
                    if count == qual_graph[quality]:
                        print '    - The base',base,'consistently had the same quality value of',quality,'(it was also the only base to receive this value)'
                        print '      Therefore, this DNA base can be encoded as',dna_bases[0],'with its unique quality score distinguishing it.'
                        N_qual[base] = quals.index(quality)
                    else:
                        print '    - The DNA base',base,'consistently had the same quality value of',quality
                        print '    - Therefore, this DNA base can be encoded as',dna_bases[0],'but with a unique/new quality distinguishing it from',dna_bases[0]
                        print '      (dont worry, it will be changed back upon decoding, this is only important if you use .uq files directly )'
                        total_quals += 1
                        N_qual[base] = total_quals # the value of an encoded qual is its index in quals, which is total quals - 1 due to 0-based indexing.

    ## Determine the minimum number of bits to store the DNA (the whole sequence is padded to the nearest byte, but not the individual bases (unless --pad used).
    if   len(dna_bases) <= 4:                    bits_per_base = 2
    elif len(dna_bases) <= 8   and not args.pad: bits_per_base = 3
    elif len(dna_bases) <= 16:                   bits_per_base = 4
    elif len(dna_bases) <= 32  and not args.pad: bits_per_base = 5
    elif len(dna_bases) <= 64  and not args.pad: bits_per_base = 6
    elif len(dna_bases) <= 128 and not args.pad: bits_per_base = 7
    else:                                        bits_per_base = 8 # Since we take a byte at a time, there can't possibly be more than this.

    ## Print some more stats
    print '    - this means we will store each letter of DNA in',bits_per_base,'bits.'
    if dna_min == dna_max:
        print '    - all DNA sequences are', dna_min, 'bases long'
    else: 
        print '    - the largest DNA sequence is',dna_max,'bases long'
        print '    - the smallest DNA sequence is',dna_min,'bases long'
    if dna_min == dna_max: variable_read_lengths = False
    else:                  variable_read_lengths = True
    dna_bits = bits_per_base * (dna_max + variable_read_lengths)
    dna_columns_needed = -(-dna_bits//8)
    print '    - therefore, as',dna_bits,'bits must be used to store the DNA for each sequence, we will use',dna_columns_needed,'bytes per unique DNA sequence.\n'

    ## Rinse and repeat for QUALs
    print 'QUAL Analysis:'
    for base,result in static_qualities.items():
        print '  [Breakdown for "'+str(base)+'"]'
        for qual in sorted(static_qualities[base]):
            percentage = (static_qualities[base][qual]/float(base_graph[base]))*100
            print '     ',qual,'|'+('#'*int(percentage/2)).ljust(50)+'|', ('%.3f'%percentage).rjust(7)+'%', str(static_qualities[base][qual]).rjust(12)
        print ''
    print '  [Total distribution]'
    print '    - the following values were seen as quality scores:',' '.join(sorted(quals)),'(' + str(len(quals)) + ' in total)'
    print '    - There distribution is:'

    for qual in quals:
        percentage = (qual_graph[qual]/total_bases)*100
        print '     ', qual, '|'+('#'*int(percentage/2)).ljust(50)+'|', ('%.3f'%percentage).rjust(7)+'%', str(qual_graph[qual]).rjust(12)

    if   total_quals <= 4:                    bits_per_quality = 2
    elif total_quals <= 8   and not args.pad: bits_per_quality = 3
    elif total_quals <= 16:                   bits_per_quality = 4
    elif total_quals <= 32  and not args.pad: bits_per_quality = 5
    elif total_quals <= 64  and not args.pad: bits_per_quality = 6
    elif total_quals <= 128 and not args.pad: bits_per_quality = 7
    else:                                     bits_per_quality = 8
    if total_quals != len(quals): print '    - as mentioned above,',total_quals-len(quals),'unique quality scores will be added to the',len(quals),'above.'
    print '    - this means we will store each quality symbol in',bits_per_quality,'bits.'
    qual_bits = bits_per_quality * (dna_max + variable_read_lengths)
    qual_columns_needed = -(-qual_bits//8)
    print '    - therefore, as',qual_bits,'bits must be used to store the qualities for each sequence, we will use',qual_columns_needed,'bytes per unique quality sequence.\n\n'

    # Turn the list into a string to reduce time taken to index() it later.
    quals     = ''.join(quals)
    dna_bases = ''.join(dna_bases)



    ## QNAMEs in FASTQ are a total mess. For many reasons. In the previous pass through the file the delimiters (characters that always appear, in the same order, a set number of times per row) 
    ## were determined. Now we can read each "column" of data out of the QNAME and encode it into binary.
    status.message = 'Pass 2 of 4: Now QNAME delimiters are known, QNAMEs are being analysed.'
    if len(separators):
        def qname_reader(inFile,prefix,suffix):
            with open(inFile,'rb') as f:
                status.current = 0
                start = len(prefix)
                end = -1 - len(suffix)
                regex = re.compile('(.*)'.join(separators))
                while True:
                    try:
                        yield re.split(regex, next(f)[start:end])
                        next(f)
                        next(f)
                        next(f)
                    except StopIteration: break
                    status.update()
        cols = False
        target = 10000
        fallback = False
        columns = []

        '''
        The following is fairly complicated, so please allow me to explain:
        All columns start off as a 'mapping' which means a dictionary/map of all known values will be used to assign values to numbers, just like the DNA/QUAL.
        Columns may get bumped down from mappings to 'integers' if the mapping gets too large - specifically if the size of the mapping contains 
        more than 1000 items, and whatever that number above 1000 is, it is more than 1/10th of the number of sequences read so far.
        Next down the type list is integers, which are always numbers. Then below that is strings, which isn't implimented yet. It's important to give integers
        the possibility of being mappings because very large numbers (say 1 billion, which 4 bytes to store) may appear rarely, and thus a mapping is much smaller. 
        However, it is also important to give mappings the chance to be integers, particularly if that doesn't require more bytes to store (but will be faster to decode).
        '''

        def check_format(columns,entries_read):
            for column in columns:
                # This is where mappings get bumped down to integers or strings if they become too large (relatively).
                if column['format'] == 'mapping':
                    if len(column['map']) > entries_read/10:
                        try:
                            _ = map(int,column['map'])
                            column['min'] = min(_)
                            column['max'] = max(_)
                            column['format'] = 'integers'
                            del column['map']
                            del _
                        except ValueError:
                            column['format'] = 'strings'
                            column['longest'] = max(map(len,column['map']))
                            del column['map']
                        gc.collect()
        qnames = qname_reader(args.input,prefix,suffix)
        for entries_read,qname in enumerate(qnames):
            if len(qname) != cols:
                if cols == False:
                    cols = len(qname)
                    for col in range(cols): columns.append({ 'name':'QNAME_'+str(col+1), 'format':'mapping', 'map':set() })
                else:
                    print 'WARNING: I was unable to correctly guess the delimiters of the QNAME in this FASTQ file.'
                    print '         Will fall back to assuming there is no delimiting of data.'
                    fallback = True
                    break
            for idx,column in enumerate(qname):
                # By default all unique items are put into a set, where every item of this set will later be encoded with a binary uint:
                if columns[idx]['format'] == 'mapping': columns[idx]['map'].add(column)
                # If that set becomes very large (relative to the qnames checked) and their all ints, we store then as u/ints directly.
                # For now we'll just store the min and max values:
                elif columns[idx]['format'] == 'integers':
                    try:
                        column = int(column)
                        if   column < columns[idx]['min']: columns[idx]['min'] = column
                        elif column > columns[idx]['max']: columns[idx]['max'] = column
                    # But if they're not even numbers to begin with, tripping this error, we'll have to store them in strings (least efficient way).
                    except ValueError: 
                        # We just store the length of the longest string for now:
                        columns[idx]['format'] = 'strings'
                        columns[idx]['longest'] = len(str(columns[idx]['max'])) # The length of the string representation of the largest number.
                        del columns[idx]['min']
                        del columns[idx]['max']
                        gc.collect()
                elif columns[idx]['format'] == 'strings':
                    if len(column) > columns[idx]['longest']: columns[idx]['longest'] = len(column)
            if entries_read == target:
                check_format(columns,entries_read)
                target *= 2
        if fallback: error('Encoding QNAMEs as strings has not been implimented yet. Please send John a message, ideally with some sample data that generates this message.')
        else: check_format(columns,entries_read) # Check the last bit of data that hadn't made it through yet.

        print 'Optimal encoding method for delimited data in QNAME determined!',status.split_time()
        for column in columns:
            # Now we give maps the chance to become integers:
            if column['format'] == 'mapping':
                map_len = len(column['map'])
                if   map_len <= 255:                  map_len = 255                  ; column['dtype'] = 'uint8'
                elif map_len <= 65535:                map_len = 65535                ; column['dtype'] = 'uint16'
                elif map_len <= 4294967295:           map_len = 4294967295           ; column['dtype'] = 'uint32'
                elif map_len <= 18446744073709551615: map_len = 18446744073709551615 ; column['dtype'] = 'uint64'
                try:
                    # Here we convert a map to an integer if it doesn't effect the dtype size (but would make encoding/decoding quicker)
                    _ = map(int,column['map'])
                    if max(_) - min(_) <= map_len:
                        column['format'] = 'integers'
                        column['max'] = max(_)
                        column['min'] = min(_)
                        if min(_) < 0 or max(_) > map_len: column['offset'] = True
                        else: column['offset'] = False
                        del column['map']; gc.collect()
                    else: column['map'] = sorted(column['map'])
                except:
                    column['map'] = sorted(column['map']) # ASCII sort the map (and change set to list for JSON encoding)

            elif column['format'] == 'integers':
                int_len = column['max']-column['min']
                if   int_len <= 255:                  int_len = 255                  ; column['dtype'] = 'uint8'
                elif int_len <= 65535:                int_len = 65535                ; column['dtype'] = 'uint16'
                elif int_len <= 4294967295:           int_len = 4294967295           ; column['dtype'] = 'uint32'
                elif int_len <= 18446744073709551615: int_len = 18446744073709551615 ; column['dtype'] = 'uint64'
                if column['min'] < 0 or column['max'] > int_len: column['offset'] = True
                else: column['offset'] = False

            elif column['format'] == 'strings':
                print 'I havent implimented this yet'; exit()

            print '    - Column',idx+1,'is type',column['format'],'stored as',column['dtype']
            gc.collect()
    else: print 'Pass 2 of 4: Skipped as there are no delimiters common to all QNAMEs'
    print '\n'

    ## The beginnings of the config.json file:
    config = {
        'base_distribution':        base_graph,
        'qual_distribution':        qual_graph,
        'reads':                    status.total,
        'bases':                    dna_bases,
        'qualities':                quals,
        'variable_read_lengths':    variable_read_lengths,
        'bits_per_base':            bits_per_base,
        'bits_per_quality':         bits_per_quality,
        'N_qual':                   N_qual,
        'dna_max':                  dna_max,
        'QNAME_prefix':             prefix,
        'QNAME_suffix':             suffix,
        'QNAME_separators':         separators,
        'QNAME_columns':            columns
    }

    if args.peek:
        print 'The config.json would look like:'
        print json.dumps(config,indent=4,sort_keys=True)
        print 'However, the sort:, raw:, and pattern: items have not been added. To get values for these, you will need to run --test with or without a --compressor and read the actual config for a file.'
        exit()

    ## Here the DNA/QUALs are actually encoded and written to disk. In an ideal world, encoded data would be immediately written to disk, and not stored in memory first.
    status.message = 'Pass 3 of 4: Encoding DNA and Quality scores...'

    if variable_read_lengths: arrays = encoder_variable(status.total,dna_bases,quals,N_qual,dna_columns_needed,qual_columns_needed,status,args.input,bits_per_base,bits_per_quality)
    else:                     arrays = encoder_fixed(status.total,dna_bases,quals,N_qual,dna_columns_needed,qual_columns_needed,status,args.input,bits_per_base,bits_per_quality)

    write_out(arrays[0][0],'DNA.temp')  # Free
    write_out(arrays[1][0],'QUAL.temp') # up
    arrays[2].free(arrays[0][1])        # system
    arrays[2].free(arrays[1][1])        # memory
    del arrays; gc.collect()            # manually
    print ''

    ## Here the QNAMEs are parsed and encoded. What a time to be alive.
    status.message = 'Pass 4 of 4: Encoding QNAMEs...'
    columns_list = []

    for column in columns: columns_list.append( ffi.new(column['dtype']+"_t["+str(status.total)+"]") )
    for row,qname in enumerate(qname_reader(args.input,prefix,suffix)):
        for column,data in enumerate(qname):
            if   columns[column]['format'] == 'mapping' : columns_list[column][row] = bisect.bisect_left(columns[column]['map'],data)
            elif columns[column]['format'] == 'integers' and columns[column]['offset']: columns_list[column][row] = int(data) - columns[column]['min']  # slows things down if not needed? if columns[column]['offset']: ...
            elif columns[column]['format'] == 'integers': columns_list[column][row] = int(data)
            elif columns[column]['format'] == 'strings' : columns_list[column][row] = data

    for idx,column in enumerate(columns):
        columns_list[idx] = numpy.asarray(  ffi.buffer(columns_list[idx]))
        # The following messing around is because numpy doesnt understand ffi.new properly, always assumed dtype is uint8, so we have to go via void to restore.
        columns_list[idx].dtype = 'V'+str(numpy.dtype(columns[idx]['dtype']).itemsize)
        columns_list[idx].dtype = columns[idx]['dtype']
        write_out(columns_list[idx],column['name']+'.temp')
    del columns_list; gc.collect()
    print ''

    ## This function runs a "mix", which is a specific combination of --sort and --raw, and perhaps others one day like delta encodings or --pad. It's reused when doing --test with different parameters.
    def run_mix(sorted_on,raw_tables,test):
        if test: test = {'sorted_on':sorted_on,'raw_tables':raw_tables}
        if sorted_on in ['DNA','QUAL']:
            not_sorted_on = 'DNA' if sorted_on == 'QUAL' else 'QUAL'
            sort_order = encode_dna_qual(False, sorted_on, True if sorted_on in raw_tables else False,test) ; gc.collect()
            encode_dna_qual(        sort_order, not_sorted_on, True if not_sorted_on in raw_tables else False,test) ; gc.collect()
            encode_qname(sort_order, True if 'QNAME' in raw_tables else False,test) ; gc.collect()
        elif sorted_on == 'QNAME':
            sort_order = encode_qname(False,True if 'QNAME' in raw_tables else False,test)   ; gc.collect()
            encode_dna_qual(sort_order, 'DNA', True if 'DNA'  in raw_tables else False,test) ; gc.collect()
            encode_dna_qual(sort_order, 'QUAL', True if 'QUAL' in raw_tables else False,test); gc.collect()
        else:
            encode_qname(None, True if 'QNAME' in raw_tables else False,test)          ; gc.collect()
            encode_dna_qual(None, 'DNA', True if 'DNA'  in raw_tables else False,test) ; gc.collect()
            encode_dna_qual(None, 'QUAL', True if 'QUAL' in raw_tables else False,test); gc.collect()

        if test:
            total_size = 0
            for key,value in test.items():
                if key.startswith('QNAME') or key.endswith('key'): total_size += value
                elif key.startswith('DNA') or key.startswith('QUAL'): total_size += value[0]
            test['total_size'] = total_size
        return test
        # Add totals to test_results here.

    ## The function used by run_mix to test/write the mix for DNA/QUALs
    def encode_dna_qual(sort_order,table_name,raw,test): # If sort_order is a numpy.ndarray, will sort output with it. If False, will sort and return a sort_order. If None, no sorting done.
        table = numpy.load(os.path.join(temp_directory,table_name+'.temp'))
        if raw:
            out_name = table_name + '.raw'
            if sort_order is None:
                if test: test[out_name] = test_patterns(table,out_name)
                else:                     write_pattern(table,out_name)
            else:
                original_dtype = table.dtype            # pypy needs this
                table.dtype = 'V' + str(len(table[0]))  # to be quick :(
                if sort_order is False: sort_order = numpy.argsort(table,axis=0)
                sort_order.shape = sort_order.shape[0]
                table = table[sort_order]
                gc.collect()
                table.dtype = original_dtype
                if test: test[out_name] = test_patterns(table,out_name)
                else:                     write_pattern(table,out_name)
        else:
            out_name = table_name + '.key'
            original_dtype = table.dtype                         # pypy currently
            table.dtype = 'V' + str(len(table[0]))               # needs this
            table,key = numpy.unique(table, return_inverse=True) # to be
            gc.collect()
            table.shape = (len(table),1)
            table.dtype = original_dtype                         # quick :(
            key = key.astype(numpy.min_scalar_type(max(key))) ######################################################## This adds considerable time to the running of the program, because pypy hasn't impliemnted it properly yet.
            gc.collect()
            if sort_order is None:
                if test: test[out_name] = compressed_size(key)
                else: write_out(key,out_name)
            else:
                if sort_order is False: sort_order = numpy.argsort(key)
                if test: test[out_name] = compressed_size(key[sort_order])
                else: write_out(key[sort_order],out_name)
            del key; gc.collect()
            table.dtype = original_dtype
            if test: test[table_name] = test_patterns(table,table_name)
            else:                       write_pattern(table,table_name)
        del table
        gc.collect()
        return sort_order

    ## The function used by run_mix to test/write the mix for QNAMEs
    def encode_qname(sort_order,raw,test): # If sort_order is False, will generate a sort_order. If sort_order is something truthy, will be used to sort. If None, no sorting.
        columns_data = []
        for column in columns:
            columns_data.append(numpy.load(os.path.join(temp_directory,column['name']+'.temp')))
        if raw:
            if sort_order is False:
                common_dtype_columns_data = numpy.dstack(columns_data)[0] # Stack makes a table where the dtype of every column is the same as the largest column of the inputs. This will be fixed later when writing to disk.
                common_dtype_columns_data.dtype = ','.join(common_dtype_columns_data.shape[1]*[str(common_dtype_columns_data.dtype)]) # Structifys the list so rows are sorted, not each columns individually...
                sort_order = numpy.argsort(common_dtype_columns_data,axis=0)
                sort_order.shape = sort_order.shape[0] 
                del common_dtype_columns_data; gc.collect()
            if type(sort_order) is numpy.ndarray:
                for idx,column in enumerate(columns):
                    if test: test[column['name']+'.raw'] = compressed_size(columns_data[idx][sort_order])
                    else: write_out(columns_data[idx][sort_order],column['name']+'.raw')
            else:
                for idx,column in enumerate(columns):
                    if test: test[column['name']+'.raw'] = compressed_size(columns_data[idx])
                    else: write_out(columns_data[idx],column['name']+'.raw')
        else:
            common_dtype_columns_data = numpy.dstack(columns_data)[0]
            common_dtype_columns_data.dtype = ','.join(common_dtype_columns_data.shape[1]*[str(common_dtype_columns_data.dtype)])
            common_dtype_columns_data, columns_key = numpy.unique(common_dtype_columns_data, return_inverse=True)
            gc.collect()
            columns_key = columns_key.astype(numpy.min_scalar_type(max(columns_key))) ######################################################## This adds considerable time to the running of the program in pypy, because pypy hasn't impliemnted it properly yet.
            if sort_order is False: sort_order = numpy.argsort(columns_key)
            if type(sort_order) is numpy.ndarray:
                if test: test['QNAME.key'] = compressed_size(columns_key[sort_order])
                else: write_out(columns_key[sort_order],'QNAME.key')
                del columns_key; gc.collect()
            else:
                if test: test['QNAME.key'] = compressed_size(columns_key)
                else: write_out(columns_key,'QNAME.key')
                del columns_key; gc.collect()
            old_shape = len(common_dtype_columns_data), len(common_dtype_columns_data[0])
            common_dtype_columns_data.dtype = common_dtype_columns_data.dtype[0]
            common_dtype_columns_data.shape = old_shape
            for idx,column in enumerate(columns):
                if test: test[column['name']] = compressed_size(common_dtype_columns_data[:,idx].astype(column['dtype']))
                else: write_out(common_dtype_columns_data[:,idx].astype(column['dtype']),column['name'])
                gc.collect()
            del common_dtype_columns_data; gc.collect()
        del columns_data; gc.collect()
        if type(sort_order) is numpy.ndarray: return sort_order


    ## If we are --test ing, try all combinations unless the user has locked in a property like a --sort or --raw DNA, etc. Patterns handled in the previous test_patterns function.
    if args.test:
        all_results = []
        print 'Starting tests...'
        print 'Time:                       Sort:      Raw Tables:                      Patterns:'
        status.split_time()
        for raw_tables in [('DNA','QUAL','QNAME'),('DNA','QUAL'),('QUAL','QNAME'),('DNA','QNAME'),('DNA',),('QUAL',),('QNAME',),(None,)] if args.raw is None else [args.raw]:
            if args.compressor is None: args.sort = (None,)
            for to_sort in ['DNA','QUAL','QNAME',None] if args.sort is None else [args.sort]:
                all_results.append(run_mix(to_sort,raw_tables,True))
                gc.collect()
                print status.split_time().ljust(27),str(to_sort).ljust(10),str(tuple(raw_tables)).ljust(32),'All' if args.raw is None else str(args.pattern)

        ## Determine the best run_mix for this file:
        best_result = sorted(all_results, key=lambda k: k['total_size'])[0]
        args.sort = best_result['sorted_on']
        args.raw = best_result['raw_tables']
        args.pattern = (best_result['DNA.raw'][1] if 'DNA.raw' in best_result else best_result['DNA'][1],best_result['QUAL.raw'][1] if 'QUAL.raw' in best_result else best_result['QUAL'][1])

        print '\nAll done!'
        if args.compressor is None: print '             Size    Sort:    Raw Tables:                 Raw stats:' 
        else:                       print 'Size (compressed)    Sort:    Raw Tables:                 Raw stats:' 
        ## Print out some output so user can see what else was tried:
        for result in sorted(all_results, key=lambda k: k['total_size']):
            print str(result['total_size']).rjust(17)+'   ',
            print str(result['sorted_on']).ljust(8),
            print str(tuple(result['raw_tables'])).ljust(27),
            del result['total_size'],result['sorted_on'],result['raw_tables']; gc.collect()
            for key,value in result.items():
                print str(key)+':'+str(value),
            print ''
        print 'Parameters found to be the best for this data type:\n  ',
        if args.sort is not None and args.sort != (None,): print ' --sort',args.sort,
        if args.raw is not None: print ' --raw',' '.join(map(str,args.raw)),
        if args.pattern is not None: print ' --pattern',' '.join(map(str,args.pattern)),
        print ''

    ## Whether the user --test'd or not, we run another mix without testing using either the user supplied settings, or the best settings from the --test.
    ## Note encoding takes so much more time than rotating/patterning/sorting data that although this is not the most efficient way of doing things, it will do just fine for now.
    if args.sort is None: args.sort = (None,)
    if args.raw is None: args.raw = (None,)
    run_mix(args.sort,args.raw,False)

    ## Add final sort/raw/pattern mix to the config:
    config['sort'] = args.sort
    config['raw'] = list(args.raw)
    config['pattern'] = args.pattern
    print '\nWriting final config...'
    path = os.path.join(temp_directory,'config.json')
    with open(path,'wb') as f: f.write(json.dumps(config,indent=4,sort_keys=True))

    ## Delete temp files, put everything into a tar archive, and move to user's --output location if specified.
    print 'Archiving results and cleaning up temp directory...'
    try:
        for f in os.listdir(temp_directory):
            if f.endswith('.temp'): os.remove(os.path.join(temp_directory,f))
        with tarfile.open(os.path.join(temp_directory,'temp.uq'), mode='w') as temp_out:
            for f in os.listdir(temp_directory): temp_out.add(os.path.join(temp_directory,f),arcname=f)
        os.rename(os.path.join(temp_directory,'temp.uq'), args.output)
        for f in os.listdir(temp_directory): os.remove(os.path.join(temp_directory,f))
    except Exception as e:
        print 'ERROR: There was an error taking the encoded data and putting it all in a single tar file:'
        print '      ',e
        print '       You might be able to still do it manually. Take a look inside',temp_directory

    print "All Done! :) "





## Below is all the code for decoding a uQ file back to FASTQ:
else:
    '''
    There are many conceivable ways to turn a uQ into a FASTQ file. The code below uses a simple approach, which is to load the whole uQ file into memory first,
    convert non-raw tables back to raw, and then iterate everything from first to last row, printing out the FASTQ formatting as it goes.
    This is definitely not the best nor fastest way to do it. The non-raw tables do not need to be converted to raw. There's no reason for everything to be loaded into memory first,
    particularly the raw tables. But at the end of the day, I don't suspect anyone will actually use the uQ format, and so sending another week optimising something no one will use
    is not the best use of the time I have on this planet :) As such, there's no need for the below to be overly optimised, only to work.
    '''

    ## Check user sanity
    if os.path.isfile(args.input):
        if not tarfile.is_tarfile(args.input): error('ERROR: Sorry, the path you have provided as input is a file, but not a tar file, and therefore cannot be a .uq file!')
    else: error('ERROR: Sorry, the path you have provided for input is not a file - please check it and try again!')
    uq = tarfile.open(args.input)

    ## Define a few very basic functions:
    def check_if_in_tar(uq,file_name): return file_name in uq.getnames()
    def load_from_tar(uq,file_name,pattern='0.1'):
        if pattern.startswith('0.'): return              numpy.load(uq.extractfile(file_name))
        else:                        return numpy.rot90( numpy.load(uq.extractfile(file_name)) , -int(pattern[0]))

    ## Check all the required stuff is there:
    if check_if_in_tar(uq,'config.json'): config = json.loads(uq.extractfile('config.json').read())
    else: error('ERROR: No config.json file was found in your input path! I cannot decode data without it!')

    if check_if_in_tar(uq,'DNA.raw'): DNA = load_from_tar(uq,'DNA.raw',config['pattern'][0])
    elif check_if_in_tar(uq,'DNA') and check_if_in_tar(uq,'DNA.key'):
        DNA = load_from_tar(uq,'DNA',config['pattern'][0])[load_from_tar(uq,'DNA.key')]
    else: error('ERROR: No DNA data was found in this uQ file?!')

    if check_if_in_tar(uq,'QUAL.raw'): QUAL = load_from_tar(uq,'QUAL.raw',config['pattern'][1])
    elif check_if_in_tar(uq,'QUAL') and check_if_in_tar(uq,'QUAL.key'): QUAL = load_from_tar(uq,'QUAL',config['pattern'][1])[load_from_tar(uq,'QUAL.key')]
    else: error('ERROR: No Quality Score data was found in this uQ file?!')

    QNAMES = []
    QNAME_KEY = False
    if check_if_in_tar(uq,'QNAME.key'):
        QNAME_KEY = load_from_tar(uq,'QNAME.key')
        for file_name in uq.getnames():
            if file_name.startswith('QNAME') and not file_name.endswith('.raw') and not file_name.endswith('.key'):
                QNAMES.append(load_from_tar(uq,file_name))
        if len(QNAMES) == 0: error('ERROR: A QNAME.key file exists but no QNAME data is present?')
    else:
        for file_name in uq.getnames():
            if file_name.startswith('QNAME') and file_name.endswith('.raw'): QNAMES.append(load_from_tar(uq,file_name))
        if len(QNAMES) == 0: error('ERROR: No QNAME data exists in this uQ file?')
    QNAMES = numpy.dstack(QNAMES)[0]
    if QNAME_KEY is not False: QNAMES = QNAMES[QNAME_KEY]

    '''
    There was once the idea to allow certain data tables to be deleted, but --decode would fill in the blanks as best it could, regarding QNAMEs.
    However, this idea was eventually scrapped, because I couldn't see anyone using it. Would save a ton of space if you didn't need QNAMEs though...
    if check_if_in_tar(uq,'QNAME.raw'):                                   QNAME = convert_qname(load_from_tar(uq,'QNAME.raw'), prefix, suffix, columnType)
    elif check_if_in_tar(uq,'QNAME') and check_if_in_tar(uq,'QNAME.key'): QNAME = convert_qname(load_from_tar(uq,'QNAME')[load_from_tar(uq,'QNAME.key')], prefix, suffix, columnType)
    #elif check_if_in_tar(uq,'QNAME.key'):                                QNAME = load_from_tar(uq,'QNAME.key') # Pairs still mated
    #else:                                                                QNAME = xrange(1,len(DNA)+1)          # Every read is unique.
    else: error('ERROR: No QNAME data was found in this uQ file?!')
    '''

    # Check everything we need is in the config file by reassigning it. This also speeds things up a little for CPython users.
    bases = config['bases']
    N_qual = config['N_qual']
    prefix = config['QNAME_prefix']
    suffix = config['QNAME_suffix']
    dna_max = config['dna_max']
    qualities = config['qualities']
    columnType = config['QNAME_columns']
    separators = config['QNAME_separators']
    bits_per_base = config['bits_per_base']
    bits_per_quality = config['bits_per_quality']
    variable_read_lengths = config['variable_read_lengths']
    total_dna_bits = (variable_read_lengths+dna_max)*bits_per_base
    total_qual_bits = (variable_read_lengths+dna_max)*bits_per_quality
    qual_N = dict( (v,k) for k,v in N_qual.iteritems()) # reverse key/val pairs

    ## DNA/QUAL decoder:
    def split_bits(numpy_array,total_bits,bits_per_x):
        bitmask = int('1'*bits_per_x,2)
        static = range(total_bits-bits_per_x,-bits_per_x,-bits_per_x)
        for row in numpy_array:
            the_number = sum([int(y) << (8*x) for x,y in enumerate(reversed(row))])
            yield [(the_number >> x) & bitmask for x in static]

    # QNAME decoder. Note that it is compiled at runtime. Check out convert_qname_fuction before the exec() to see what it's up too :)
    convert_qname_fuction = '''
def convert_qname(numpy_array,prefix,suffix,columnType):
    for row in numpy_array:
        yield "'''+prefix+'" + '
    for idx,column in enumerate(columnType):
        if column['format'] == 'mapping':
            convert_qname_fuction += 'columnType['+str(idx)+']["map"][row['+str(idx)+']] + '
        if column['format'] == 'integers':
            if column['offset'] == False: convert_qname_fuction += 'str(row['+str(idx)+']) + '
            else:                         convert_qname_fuction += 'str(row['+str(idx)+']+'+str(columnType[idx]['min'])+') + '
        if column['format'] == 'string':
            error('ERROR: I dont support string encoding yet.')
        try: convert_qname_fuction += '"'+separators[idx]+'" + '
        except: pass # The last column has no separator
    convert_qname_fuction += '"'+suffix+'"'
    #print convert_qname_fuction
    exec(convert_qname_fuction) # This is much more efficient than having to parse columnType for every QNAME entry. We compile a decoder once and reuse it.

    ## The above decoders are just generator functions. The below zips all those generators together to give decoded DNA/QUAL/QNAME in ASCII:
    try:
        gen = itertools.izip(split_bits(DNA,total_dna_bits,bits_per_base), split_bits(QUAL,total_qual_bits,bits_per_quality), convert_qname(QNAMES, prefix, suffix, columnType))
        if len(N_qual) != 0:
            for dna,qual,qname in gen:
                # Here we decode the encoded value into the base/quality value, making sure to convert the N_qual values too.
                for idx,q in enumerate(qual):
                    qual[idx] = qualities[q]
                    if q in qual_N: dna[idx] = qual_N[q]
                    else:           dna[idx] = bases[dna[idx]]
                ## remove padding if variable length.
                if variable_read_lengths:
                    dna = dna[1+dna.index(bases[1]):]
                    qual = qual[1+qual.index(qualities[1]):]
                print qname
                print ''.join(dna)
                print '+'            # I'm seeing the qname here in a lot of publicly avalible files instead of +. Perhaps could add a flag to config.json?
                print ''.join(qual)
        else:
            # We can decode everything a little bit faster if theres nothing in N_qual/qual_N to convert.
            for dna,qual,qname in gen:
                for idx,q in enumerate(qual):
                    dna[idx] = bases[dna[idx]]
                    qual[idx] = qualities[q]
                if variable_read_lengths:
                    dna = dna[1+dna.index(bases[1]):]
                    qual = qual[1+qual.index(qualities[1]):]
                print qname
                print ''.join(dna)
                print '+'
                print ''.join(qual)
    except IOError:
        pass

'''
Developer To Dos / Notes:

1) Try encoding QNAMEs in strings (and support strings if mapping fails)

2) To reduce the final output filesize further, i would like to write all reads containing Ns to a separate file, then once the unique/sorted
   list is created (for DNA/QUAL), try replacing the N with any other letter to get a fragment that is already in the unique dataset. Since the list
   is sorted, we can use a binary intersect and figure that out near-instantaniously, so it's reasonable. If no match is found, replace Ns with the most
   common base, which is what is currently done. This is part of the LOSSLESS part of the program, so really it would benefit all.

3) DNA could actually be stored in ACGTrie format, rather than a table. Idx points to the end row in the trie (and working backwards gets you the DNA)
   Same could work for the QUAL too.

6) There are no repeated QNAMEs (or few) so key method isn't working so well. Should make a key per column after sort/unique'ing columns individually, then even sort/unique the stack of those keys.

'''
