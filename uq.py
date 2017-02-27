#!/usr/bin/env python

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

'''
Hello,
This program takes FASTQ data and re-encodes it to be as small as possible without compression,
(for large in-memory calculations with hundreds/thousands of FASTQ files being used simultaniously).
It also lays out the data on the disk to be optimised for compressors to pack it into a really tiny
file, for long term storage or mass distrubution - often several times smaller than what those same
compression programs would have acheived on a standard FASTQ file.

Note, it does not do any compression itself - it is simply a repackaging of the existing data which still
allows for instant random reading/writing (using the right functions) just like a typical FASTQ file. This 
repacking/encoding procedure is completely lossless, safe, and reversible. Unlike other encoders, this 
program does not have any concept of DNA, phred scores, CASAVA naming conventions, etc etc. There is no 
"interpretation" of the data - no expectations of what your data "should" look like - thus a lot of the issues 
which could lead to mistakes or loss of information when compressing are avoided.

On the topic of safety, uQ.py is one of the only FASTQ encoding/compression programs that reads the whole
file first before deciding on how to encode the data optimally/safely. This makes it a little slower than
other tools, but arguably much safer.

But it's important to remember that this isn't just about storing ASCII bytes in the smallest number of bits. 
It's about storing those bits in such a way that compressors can very efficiently compress the data - whatever 
compressor you choose to use. I use the very effient lpac9 algorithm for my own data, although other algorithms
from the same fantastic developer (Matt Mahoney) can get the files even smaller if you have more cpu cycles to spare.
For distribution however, I would use something more prevelent like gzip or 7zip, as that's what your users will
already have installed.

The decoding of .uq files to .fastq can also be done with this tool, and this tool is included in the .uq file itself,
which is really just a tar file holding some binary blobs. To decode the data, run this tool with --decode, and 
point the --input to either the .uq tar file, or the untarred .uq directory (should contain a config.json).

Dependancies for encoding/decoding are python 2 with numpy, and unix tool "wc". If you would like me to remove the 
numpy/wc dependancies, and/or port the code to python 3, i'm happy to do that for a beer :)

In the future this tool can support some lossy parameters, most being exactly what you'd expect (renaming read IDs, 
binning quality scores, etc), one is worth looking at, --reorder, because reordering the reads in the file is not 
seen as 'lossy' to many people, and that one currently works.

The three future lossy methods are
--qname   Change the read names from whatever they are off the sequencer to just a constantly increasing number. 
          You must specify either "--qname unique" if you want every read to have a unique number/name and you dont
          care about the order, and "--qname duplicates" if you want read names in the original file that are duplicated
          (for example, in paired-end sequencing), to also have the same number. Regarding space savings, the former is
          like deleting the .QNAME and .QNAME.key tables from the final output, while the latter is like deleting just
          the .QNAME table but keeping the .QNAME.key table.

--qual    Bins quality values so there are less possible quality scores to have to remember. Again you must specify one
          of two options: "--qual nearest N" which will bin similar-values together, into N bins. The new value of binned
          values is the value from the original distribution with the highest count. For example, if you had 10 quality 
          scores, "--qual nearest 5" would bin 1&2, 3&4, 5&6, 7&8, 9&10, to give 1,3,6,7,10 - for example - because
          quality score 1 had more counts than quality score 2, 3 had more than 4, 6 more than 5, etc. If all bins cannot
          contain the same number of quality scores in them, for example 10 original scores and --qual 3 is used, then
          the lowest scores (1,2,3,4) will be binned together, not the highest (8,9,10). You can see the distributions by
          running --test.
          "--qual weighted N" however takes the distribution into consideration. It tries to bin the data into N bins 
          where the bin with the lowest score gets merged with a neighbour first, until there are N bins remaining. 
          This often results in more 'natural' looking distributions - however, the binning method you should use really 
          depends on the tools you use downstream. In terms of space-savings, it depends on how 'complex' the sequencing
          was. Genomic sequencing, etc, will have so much complexity in the DNA that what's going on in the qualities is 
          sort of irrelevent. Transcription Factor ChIP-Seq however, would benefit greatly, as quality score complexity 
          makes up a larger proportion of the total data complexity/entropy.

--dna     Two more options for you again, "--dna drop", which drops any read with an N in it (N's are identified as bases 
          with only 1 quality score). The other is "--dna hide", which converts any N base into whatever the most common 
          base in the file was, but give it the lowest possible quality score seen in the file. Mappers, will treat these 
          bases like an N, but other programs (like GATKs BaseRecalibrator) will not like it at all, and you're better off
          dropping the whole fragment. The space-savings here depend on the number of Ns in your sequencing.

--reorder The order of the reads in the input FASTQ will not be stored - or rather, the data is re-sorted for optimum compression. 
          The resulting FASTQ file will have all the same information as the original, but the reads will be in a different order 
          than in the actual file.

Note: I'm starting to think I should scrap all lossy methods except --order and --qname, because messing with the 
      sequencing/quality data leaves a bad precident. 

'''

parser = argparse.ArgumentParser(     description="Put in a FASTQ file, get out a microq (uq) file.")
parser.add_argument("-i", "--input",  required=True,                      help='Required. Input FASTQ file path.')
parser.add_argument("-o", "--output",                                     help='Optional. Default is to append .uq to input filename.')
parser.add_argument("--raw",          nargs='+', metavar='file name',     help='Optional. Saves data as a raw table rather than unique & sorted +key combo. [DNA/QUAL/QNAME]')
parser.add_argument("--temp",                                             help='Optional. Directory to write temporary files to.')
#parser.add_argument("--dna",                                              help='Optional. [drop/hide] - see readme/code for details.')
#parser.add_argument("--qual",                                             help='Optional. [nearest/weighted] [N] - see readme/code for details.')
#parser.add_argument("--qname",                                            help='Optional. [unique/duplicate] - see readme/code for details.')
parser.add_argument("--compressor",   action='store',                     help='Optional. Path to compression program.')
parser.add_argument("--pattern",      nargs='+', metavar='pattern',       help='Optional. See output of --test for details.')
#parser.add_argument("--reorder",      action='store_true', default=False, help='Optional. No parameters - see readme/code for details.')
parser.add_argument("--test",         action='store_true', default=False, help='Optional. Writes a uq tar file with every possible data output orientation. Compress them all with your compressor of choice to find the appropriate values for --raw.')
parser.add_argument("--pad",          action='store_true', default=False, help='Optional. Pads DNA/QUAL to the nearest 2/4/8 bits. Some compressors do a better job when data is padded.')
parser.add_argument("--peek",         action='store_true', default=False, help='Optional. No output files are created, the input is just scanned for the parameters that *would* be used.')
#parser.add_argument("--paranoid",     action='store_true', default=False, help='Optional. After encoding, will attempt to decode the output and compare its MD5 to the original input FASTQ (only works for lossless encoding)')
#parser.add_argument("--delete",       action='store_true', default=False, help='Optional. Input file is deleted if no issues occur during encoding. If --paranoid used, only deletes if this returns successfully.')
parser.add_argument("--decode",       action='store_true', default=False, help='Requred if you want to decode a .uq file back to .fastq :)')
args = parser.parse_args()

if not args.pattern: args.pattern = ['1.1','1.1'] # default pattern
elif len(args.pattern) != 2: print 'ERROR: There must be 2 values for --pattern!'; exit()
elif not all([pattern in ['1.1','1.2','1.1','1.2','1.1','1.2','1.1','1.2'] for pattern in args.pattern]): print 'ERROR: Pattern values are incorrect!'; exit()

#if args.qname or args.qual or args.dna or args.paranoid or args.delete: print 'Not implimented yet.'; exit()

if args.raw:
    for x in args.raw:
        if x not in ['DNA','QUAL','QNAME']:
            print 'ERROR: When using --raw, you must specificy if which tables, "DNA", "QUAL" or "QNAME", you wish to keep as raw.'
            exit()
    args.raw = set(args.raw)
else: args.raw = set()

if args.decode:
    if os.path.isfile(args.input):
        if not tarfile.is_tarfile(args.input): print 'ERROR: Sorry, the path you have provided as input is a file, but not a tar file, and therefore cannot be a .uq file!'; exit()
        else: is_tar = True
    elif os.path.isdir(args.input):
        if not os.path.isfile(os.path.join(args.input,'config.json')):
            print 'ERROR: Sorry, the DIRECTORY you have used as a decoding input does not include the config.json file.'; exit()
        else: is_tar = False
    else: print 'ERROR: Sorry, the path you have provided for input does not exist - please check it and try again :)'; exit() 
else:
    if not os.path.isfile(args.input):  print 'ERROR: Sorry, the input path you have specified is not a file path. Please check it and try again! :)'; exit()
    temp_directory = args.temp if args.temp else tempfile.mkdtemp()
    if args.output == None: args.output = args.input + '.uq'
    print temp_directory

    # A small class we can use to query the current status. On initlization it wc -l's the file. Will extend in the future to unzip gzipped fastq.
    class init_status:
        def __init__(self,inFile):
            print 'Warming up...'
            self.total = int(subprocess.check_output('wc -l ' + inFile,shell=True).split()[0])
            if self.total % 4 == 0: self.total /= 4
            else: print 'ERROR: The FASTQ file provided contains',self.total,'rows, which is not divisible by 4!'; exit()
            self.current = 0
            self.time = time.time()
            self.split = time.time()
            self.message = ''
        def update(self):
            self.current += 1
            if self.current & 8191 == 8191: # an efficient way of saying "do something every couple FASTQ entires"
                sys.stdout.write('\033[A')
                sys.stdout.flush()
                print self.message, '(' + str( (self.current/float(self.total))*100)[:5] + '%)'
        def split_time(self):
            split = time.time() - self.split
            self.split = time.time()
            return '(' + str(split/60) + ' minutes)'
        def total_time(self): return time.time() - self.time
    status = init_status(args.input)


    status.message = 'Step 1 of 4: Analysing file to find an optimal/safe encoding strategy.'
    ## Read first 4 lines:
    with open(args.input,'rb') as f:
        line1,line2,line3,line4 = next(f)[:-1],next(f)[:-1],next(f),next(f)[:-1]
        status.update()

        # Some very basic FASTQ checks. Used to be more checks, but FASTQ is such a silly format most checks end up being overzealous or true for FASTA too.
        if not line1.startswith('@'): print 'ERROR: This does not look like a FASTA/FASTQ file! (first line does not start with @)'; exit()

        # QNAME
        prefix = line1
        suffix = line1[:]
        separators = collections.defaultdict(int)
        not_separators = set()

        # DNA
        dna_bases = set(line2)
        dna_min = len(line2)
        dna_max = len(line2)

        # +
        if not line3.startswith('+'): print 'ERROR: This does not look like a FASTA/FASTQ file! (third line does not start with +)'; exit()
        '(no data collected, always decodes to a +)'

        # QUAL
        quals = set(line4)

        if len(line2) != len(line4): print 'ERROR: This does not look like a FASTA/FASTQ file! (SEQ and QUAL lines are not the same length)'; exit()

        # N encoded in quality score:
        static_qualities = {}
        for base,qual in itertools.izip(line2,line4):
            try:
                static_qualities[base][qual] += 1
            except KeyError:
                static_qualities[base] = collections.defaultdict(int)
                static_qualities[base][qual] += 1

        # Process the rest of file:
        while True:
            try:
                qname = next(f)[:-1]
                dna = next(f)[:-1]
                if next(f)[0] != '+': print "ERROR: For entry",status.current,"the third line does not start with +"; exit()
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
                for idx,character in enumerate(suffix[::-1]):
                    if character != qname[-1-idx]:
                        suffix = suffix[-idx:]
                        break

            for sep in separators.copy():
                if qname[len(prefix):].count(sep) != separators[sep]:
                    del separators[sep]
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

    # Tidy up separator suffix results:
    for sep in separators.copy():
        if suffix.count(sep) != 0:
            separators[sep] -= suffix.count(sep)
            if separators[sep] == 0: del separators[sep]

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

    # Get data for Base/Quality distribution charts. (could print a chart per base plus a combined? Could scale to largest count rather than 100%? Meh.)
    base_graph = collections.defaultdict(int)
    qual_graph = collections.defaultdict(int)
    for base,counts_per_qual in static_qualities.items():
        base_graph[base] = sum(counts_per_qual.values())
        for qual,count in counts_per_qual.items(): qual_graph[qual] += count
    total_bases = float(sum(base_graph.values()))
    dna_bases = map(lambda x: x[0],sorted(base_graph.items(), key=lambda v: v[1], reverse=True))
    quals     = map(lambda x: x[0],sorted(qual_graph.items(), key=lambda v: v[1], reverse=True))
    #dna_bases = sorted(base_graph.keys())  # Jury is out if the above actually helps and this isnt
    #quals     = sorted(qual_graph.keys())  # smaller - need to run this on more data.

    # Print output:
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
    for x,base in enumerate(dna_bases):
        percentage = (base_graph[base]/total_bases)*100
        print '     ',base,'|'+('#'*int(percentage)).ljust(100)+'|', str(x).rjust(3), ('%.3f'%percentage).rjust(7)+'%', str(base_graph[base]).rjust(12)

    # Decide if special N encoding will be used.
    N_qual = {}
    total_quals = len(quals)
    for base,result in static_qualities.items():
        if len(dna_bases) == 1: continue # We need at least 1 base after all!
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

    if   len(dna_bases) <= 4:                    bits_per_base = 2
    elif len(dna_bases) <= 8   and not args.pad: bits_per_base = 3
    elif len(dna_bases) <= 16:                   bits_per_base = 4
    elif len(dna_bases) <= 32  and not args.pad: bits_per_base = 5
    elif len(dna_bases) <= 64  and not args.pad: bits_per_base = 6
    elif len(dna_bases) <= 128 and not args.pad: bits_per_base = 7
    else:                                        bits_per_base = 8 # Since we take a byte at a time, there can't possibly be more than this.
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

    print 'QUAL Analysis'
    print '    - the following values were seen as quality scores:',' '.join(sorted(quals)),'(' + str(len(quals)) + ' in total)'
    print '    - There distribution is:'
    for x,qual in enumerate(quals):
        percentage = (qual_graph[qual]/total_bases)*100
        print '     ',qual,'|'+('#'*int(percentage)).ljust(100)+'|', str(x).rjust(3), ('%.3f'%percentage).rjust(7)+'%', str(qual_graph[qual]).rjust(12)
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









    status.message = 'Step 2 of 4: Encoding DNA and Quality Scores using optimized settings.'
    if not args.peek:
        class encoder:
            def this(self):
                self.status.current = 0
                with open(self.file_path,'rb') as f:
                    for row in xrange(self.total_reads):
                        try:
                            next(f)
                            dna = next(f)[:-1]
                            next(f)
                            quals = next(f)[:-1]
                            #if 'N' in dna: continue
                        except StopIteration: break
                        temp_dna = self.variable_read_lengths  #   either 1 or 0
                        temp_qual = self.variable_read_lengths # for True or False
                        for base,quality in itertools.izip(dna,quals):
                            try:
                                temp_dna = (temp_dna << self.bits_per_base) + self.bases.index(base)
                                temp_qual = (temp_qual << self.bits_per_quality) + self.qualities.index(quality)
                            except ValueError:
                                temp_dna = (temp_dna << self.bits_per_base) + self.N_base
                                temp_qual = (temp_qual << self.bits_per_quality) + self.N_qual[base]
                        self.dna_array[row] = tuple((temp_dna >> x) & 255 for x in self.dna_magic)
                        self.qual_array[row] = tuple((temp_qual >> x) & 255 for x in self.quality_magic)
                        self.status.update()

                return (self.dna_array,self.qual_array)

            def __init__(self,total_reads,bases,qualities,N_qual,dna_bytes_per_row,quality_bytes_per_row,variable_read_lengths,status,file_path):
                self.status = status
                self.file_path = file_path
                self.dna_magic = range( (dna_bytes_per_row-1) * 8,-8,-8)
                self.quality_magic = range( (quality_bytes_per_row-1) * 8,-8,-8)
                self.total_reads = total_reads
                self.bases = ''.join(bases)
                self.qualities = ''.join(qualities)
                self.variable_read_lengths = int(variable_read_lengths)
                self.bits_per_base = (len(bases)-1).bit_length()
                self.bits_per_quality = len(qualities).bit_length()
                self.N_base = 0 # The most common base (with presumably the smallest post-compression code)
                self.N_qual = N_qual
                self.dna_array  = numpy.zeros(total_reads,numpy.dtype(','.join(    dna_bytes_per_row *['uint8'])))
                self.qual_array = numpy.zeros(total_reads,numpy.dtype(','.join(quality_bytes_per_row *['uint8'])))
                self.more_dna_magic = range( dna_bits - self.bits_per_base ,-self.bits_per_base,-self.bits_per_base)

        encode = encoder(status.total,dna_bases,quals,N_qual,dna_columns_needed,qual_columns_needed,variable_read_lengths,status,args.input)

        dna_array,qual_array = encode.this()
        print 'Finished encoding',status.split_time(),'sorting qualities...'

        def struct_to_std(ar):
            length = len(ar)
            columns = len(ar[0])
            ar.dtype = 'uint8'
            ar.shape = (length,columns)

        def std_to_struct(ar):
            length,columns = ar.shape
            ar.dtype = ','.join(columns*['uint8'])
            ar.shape = length

        def write_pattern(table,filename):
            struct_to_std(table)
            if filename.startswith('DNA'):    pattern = args.pattern[0]
            elif filename.startswith('QUAL'): pattern = args.pattern[1]
            with open(os.path.join(temp_directory,filename),'wb') as f:
                if pattern == '1.1': numpy.save( f, numpy.ascontiguousarray(             table    ))
                if pattern == '2.1': numpy.save( f, numpy.ascontiguousarray( numpy.rot90(table,1) ))
                if pattern == '3.1': numpy.save( f, numpy.ascontiguousarray( numpy.rot90(table,2) ))
                if pattern == '4.1': numpy.save( f, numpy.ascontiguousarray( numpy.rot90(table,3) ))
                if pattern == '2.2': numpy.save( f, numpy.asfortranarray(                table    ))
                if pattern == '3.2': numpy.save( f, numpy.asfortranarray(    numpy.rot90(table,1) ))
                if pattern == '4.2': numpy.save( f, numpy.asfortranarray(    numpy.rot90(table,2) ))
                if pattern == '1.2': numpy.save( f, numpy.asfortranarray(    numpy.rot90(table,3) ))

            std_to_struct(table)

        def write_out(table,filename):
            #struct_to_std(table)
            with open(os.path.join(temp_directory,filename),'wb') as f: numpy.save( f, table )

        def read_in(filename):
            return numpy.load(os.path.join(temp_directory,filename))

        def compressed_size(table):
            p = subprocess.Popen(args.compressor + ' | wc -c', stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
            numpy.save(p.stdin,table)
            return int(p.communicate()[0])

        def test_patterns(table):
            if not args.compressor: return [table.size, '1.1'] # If the user isn't using compression, patterning won't make a difference, but keys might, so return size
            results = []
            struct_to_std(table)
            # Rotate & Flip -- Note, flipping with flipud/fliplr does nothing the below doesn't catch. Endian might make a difference...
            results.append([ compressed_size(numpy.ascontiguousarray(table)               ), '1.1' ])
            results.append([ compressed_size(numpy.ascontiguousarray(numpy.rot90(table,1))), '2.1' ])
            results.append([ compressed_size(numpy.ascontiguousarray(numpy.rot90(table,2))), '3.1' ])
            results.append([ compressed_size(numpy.ascontiguousarray(numpy.rot90(table,3))), '4.1' ])
            results.append([ compressed_size(numpy.asfortranarray(table)                  ), '2.2' ])
            results.append([ compressed_size(numpy.asfortranarray(numpy.rot90(table,1))   ), '3.2' ])
            results.append([ compressed_size(numpy.asfortranarray(numpy.rot90(table,2))   ), '4.2' ])
            results.append([ compressed_size(numpy.asfortranarray(numpy.rot90(table,3))   ), '1.2' ])
            for result,name in sorted(results): print name.rjust(10) + ':' + str(result).rjust(10)
            std_to_struct(table)
            return sorted(results)[0]

        test_sizes = {}

        write_out(dna_array,'DNA.temp'); del dna_array # Free up some RAM for the impending sort. Will process it later. Probably a waste of time for big compute servers.

        if 'QUAL' in args.raw or args.test:
            if args.test:
                print 'QUAL.raw:'
                test_sizes['QUAL.raw'] = test_patterns(qual_array)
                args.pattern[1] = test_sizes['QUAL.raw'][1]
                print ''
            write_pattern(qual_array,'QUAL.raw')

        if 'QUAL' not in args.raw or args.test:
            original_dtype = qual_array.dtype
            qual_array.dtype = 'V' + str(len(qual_array[0])) # pypy is much faster if the row type is just a single void, and not a whole bunch of columns.
            qual_array,qual_key = numpy.unique(qual_array, return_inverse=True)
            if   len(qual_array) <= 256:                  qual_key_dtype = 'uint8'
            elif len(qual_array) <= 65536:                qual_key_dtype = 'uint16'
            elif len(qual_array) <= 4294967296:           qual_key_dtype = 'uint32'
            elif len(qual_array) <= 18446744073709551616: qual_key_dtype = 'uint64'
            write_out(qual_key,'QUAL.key')
            if args.test:
                if args.compressor: test_sizes['QUAL.key'] = compressed_size(qual_key)
                else:               test_sizes['QUAL.key'] = qual_key.size
            del qual_key
            qual_array.dtype = original_dtype
            if args.test:
                print 'QUAL:'
                test_sizes['QUAL'] = test_patterns(qual_array)
                args.pattern[1] = test_sizes['QUAL'][1]
                print ''
            write_pattern(qual_array,'QUAL')
        del qual_array

        print 'Finished sorting qualities',status.split_time(),'sorting DNA...'

        dna_array = read_in('DNA.temp')
        os.remove(os.path.join(temp_directory,'DNA.temp'))

        if 'DNA' in args.raw or args.test:
            if args.test:
                print 'DNA.raw:'
                test_sizes['DNA.raw'] = test_patterns(dna_array)
                args.pattern[0] = test_sizes['DNA.raw'][1]
                print ''
            write_pattern(dna_array,'DNA.raw')

        if 'DNA' not in args.raw or args.test:
            original_dtype = dna_array.dtype
            dna_array.dtype = 'V' + str(len(dna_array[0])) # pypy is much faster if the row type is just a single void, and not a whole bunch of columns.
            dna_array,dna_key = numpy.unique(dna_array, return_inverse=True)
            if   len(dna_array) <= 256:                  dna_key_dtype = 'uint8'
            elif len(dna_array) <= 65536:                dna_key_dtype = 'uint16'
            elif len(dna_array) <= 4294967296:           dna_key_dtype = 'uint32'
            elif len(dna_array) <= 18446744073709551616: dna_key_dtype = 'uint64'
            write_out(dna_key,'DNA.key')
            if args.test:
                if args.compressor: test_sizes['DNA.key'] = compressed_size(dna_key)
                else: test_sizes['DNA.key'] = dna_key.size
            del dna_key
            dna_array.dtype = original_dtype
            if args.test:
                print 'DNA:'
                test_sizes['DNA'] = test_patterns(dna_array)
                args.pattern[0] = test_sizes['DNA'][1]
                print ''
            write_pattern(dna_array,'DNA')
        del dna_array
        print 'Encoded DNA and Qualities written to disk!\n',status.split_time()

    status.message = 'Step 3 of 4: Analysing QNAME data between delimiters.'
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
        columns explination:
        All columns start off as 'coded' which means a dictionary/mapping of all known values will be used to assign values to numbers, just like the DNA/QUAL.
        Columns may get bumped down from mappings to 'integers' or 'strings' format if this mapping gets very large - specifically if the size of the mapping contains 
        more than 1000 items, and whatever that number above 1000 is, it is more than 1/10th of the number of sequences read so far.
        Next down the type list is integers, which are always numbers. Then below that is strings. It's important to give integers the possibility of being
        mappings because very large numbers may appear rarely, and thus a mapping is much smaller.
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
                elif columns[idx]['format'] == 'strings':
                    if len(column) > columns[idx]['longest']: columns[idx]['longest'] = len(column)

            if entries_read == target:
                check_format(columns,entries_read)
                target *= 2
        if fallback:
            print 'I havent implimented this yet!'
            exit()
        else: check_format(columns,entries_read) # Check the last bit of data that hadn't made it through yet.

        print 'Optimal encoding method for delimited data in QNAME determined!',status.split_time()

        for column in columns:
            # Now we fiddle with things a bit:
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
                        del column['map']
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
    else: print 'Step 3 of 4: Skipped as there are no delimiters common to all QNAMEs'

    decoder_ring = {
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
        for x,y in decoder_ring.items():
            print type(x),type(y)
        print json.dumps(decoder_ring,indent=4)

    else:
        print 'Writing config to disk...',
        path = os.path.join(temp_directory,'config.json')
        with open(path,'wb') as f:
            f.write(json.dumps(decoder_ring,indent=4))
        print 'success!\n'

        status.message = 'Step 4 of 4: Encoding QNAME with optimized settings.'
        columns_data = []
        for column in columns: columns_data.append( numpy.zeros(decoder_ring['reads'],column['dtype']) )
        for entries_read,qname in enumerate(qname_reader(args.input,prefix,suffix)):
            for column_idx,column_data in enumerate(qname):
                if columns[column_idx]['format'] == 'mapping':                                      columns_data[column_idx][entries_read] = bisect.bisect_left(columns[column_idx]['map'],column_data)
                elif columns[column_idx]['format'] == 'integers' and columns[column_idx]['offset']: columns_data[column_idx][entries_read] = int(column_data) - columns[column_idx]['min']
                elif columns[column_idx]['format'] == 'integers':                                   columns_data[column_idx][entries_read] = int(column_data)
                elif columns[column_idx]['format'] == 'strings':                                    columns_data[column_idx][entries_read] = column_data

        if args.test or 'QNAME' in args.raw:
            if args.test: test_sizes['QNAME.raw'] = 0; print ''
            for idx,column in enumerate(columns):
                write_out(columns_data[idx],column['name']+'.raw')
                if args.test:
                    if args.compressor: this_size = compressed_size(columns_data[idx])
                    else:               this_size = columns_data[idx].size
                    print column['name']+'.raw :',this_size
                    test_sizes['QNAME.raw'] += this_size
            if args.test: print 'Total QNAME.raw :',test_sizes['QNAME.raw'],'\n'

        if args.test or 'QNAME' not in args.raw:
            # Now we just check how many QNAMEs were repeated
            common_dtype_columns_data = numpy.stack(columns_data, axis=1)
            common_dtype_columns_data.dtype = ','.join(common_dtype_columns_data.shape[1]*[str(common_dtype_columns_data.dtype)]) # Given the dtype of the largest array automatically during stack()
            print 'Sorting QNAMEs.'
            common_dtype_columns_data, columns_key = numpy.unique(common_dtype_columns_data, return_inverse=True)
            old_shape = len(common_dtype_columns_data), len(common_dtype_columns_data[0])
            common_dtype_columns_data.dtype = common_dtype_columns_data.dtype[0]
            common_dtype_columns_data.shape = old_shape
            for column_idx in range(len(columns)):
                columns_data[column_idx] = common_dtype_columns_data[:,column_idx].astype(columns[column_idx]['dtype'])

            if args.test: test_sizes['QNAME'] = 0; print ''
            for idx,column in enumerate(columns):
                write_out(columns_data[idx],column['name'])
                if args.test:
                    if args.compressor: this_size = compressed_size(columns_data[idx])
                    else:               this_size = columns_data[idx].size
                    print column['name']+' :',this_size
                    test_sizes['QNAME'] += this_size
            if args.test: print 'Total QNAME :',test_sizes['QNAME'],

            write_out(columns_key,'QNAME.key')
            if args.test:
                if args.compressor: test_sizes['QNAME.key'] = compressed_size(columns_key)
                else:               test_sizes['QNAME.key'] = columns_key.size
            if args.test: print '( plus',test_sizes['QNAME.key'],'for QNAME.key!)','\n'
            del columns_key

            if len(columns_data[0]) != status.total:
                if len(columns_data[0])*2 != status.total:
                    print 'INFO: Of the',status.total,'reads in the file,',status.total-len(columns_data[0]),'contained duplicated read names/headers/QNAMEs.'
                    print '      What is unusual is that this number is not half the file (which might make sense for paired-end sequencing).'
                    print '      Thus, there may be a mistake in your QNAME formatting, however, it is not my place to judge - encoding will be totally uneffected.'
                else: 
                    print 'INFO: This file looks like it contains properly-formatted paired end data.'
            else:
                print 'INFO: This file looks like it contains single-end sequencing data, or is 1 of 2 paired-end files.'
                print '      If that is NOT the case, then there may be a mistake in your QNAME formatting. Otherwise everything is fine :)'

        print '\nEncoded QNAMEs written to disk!',status.split_time()

        test_raw = []
        test_pattern = []

        if args.test:
            #print 'Bonus Step: Reordering reads to optimize compression.'
            #master = numpy.stack([ read_in('DNA.key'), read_in('QNAME.key'), read_in('QUAL.key') ], axis=1)
            #master = master[master[:,2].argsort(kind='quicksort')]
            #master = master[master[:,1].argsort(kind='mergesort')]
            #master = master[master[:,0].argsort(kind='mergesort')]
            #write_out(master,'master.key') # This should be rotated or individual columns again.

            if test_sizes['QUAL.raw'][0] - (test_sizes['QUAL'][0] + test_sizes['QUAL.key']) > 0:
                #os.remove(os.path.join(temp_directory,'QUAL'))
                #os.remove(os.path.join(temp_directory,'QUAL.key'))
                test_raw.append('QUAL')
                test_pattern.append(test_sizes['QUAL.raw'][1])
            else:
                #os.remove(os.path.join(temp_directory,'QUAL.raw'))
                test_pattern.append(test_sizes['QUAL'][1])

            if test_sizes['DNA.raw'][0] - (test_sizes['DNA'][0] + test_sizes['DNA.key']) > 0:
                #os.remove(os.path.join(temp_directory,'DNA'))
                #os.remove(os.path.join(temp_directory,'DNA.key'))
                test_raw.append('DNA')
                test_pattern.append(test_sizes['DNA.raw'][1])
            else:
                #os.remove(os.path.join(temp_directory,'DNA.raw'))
                test_pattern.append(test_sizes['DNA'][1])

            if test_sizes['QNAME.raw'] - (test_sizes['QNAME'] + test_sizes['QNAME.key']) > 0:
                for idx,column in enumerate(columns): os.remove(os.path.join(temp_directory,column['name']))
                #os.remove(os.path.join(temp_directory,'QNAME.key'))
                test_raw.append('QNAME')
            else:
                for idx,column in enumerate(columns): os.remove(os.path.join(temp_directory,column['name']+'.raw'))

            print 'TEST COMPLETE: The following options were found to be optimal:'
            print '    --raw',' '.join(test_raw),' --pattern',' '.join(test_pattern)
            print '\n',test_sizes

        try:
            with tarfile.open(os.path.join(temp_directory,'temp.uq'), mode='w') as temp_out:
                for f in os.listdir(temp_directory): temp_out.add(os.path.join(temp_directory,f),arcname=f)
            os.rename(os.path.join(temp_directory,'temp.uq'), args.output)
            for f in os.listdir(temp_directory): os.remove(os.path.join(temp_directory,f))

        except Exception as e:
            print 'ERROR: There was an error taking the encoded data and putting it all in a single tar file:'
            print '      ',e
            print '       You might be able to still do it manually. Take a look inside',temp_directory



if args.decode: # or args.paranoid
    if args.decode:
        decode_input = tarfile.open(args.input) if is_tar else args.input    # for decode
    else:
        decode_input = tarfile.open(args.output) if is_tar else args.output    # for paranoid

    # Define a few very basic functions:
    def check_dir(decode_input,file_name): return os.path.isfile(os.path.join(decode_input,file_name))
    def check_tar(decode_input,file_name): return file_name in decode_input.getnames()
    def load_tar(decode_input,file_name):
        return numpy.load(decode_input.extractfile(file_name))
        #if data['rotations'] == 0: return data['table']
        #else:                      return numpy.rot90(data['table'],-data['rotations']) # undo the rotations. Could undo a diff here too.
    def load_dir(decode_input,file_name):
        return numpy.load(os.path.join(decode_input,file_name))
        #if data['rotations'] == 0: return data['table']
        #else:                      return numpy.rot90(data['table'],-data['rotations']) # undo the rotations. Could undo a diff here too.
    check_file = check_tar if is_tar else check_dir
    load_file = load_tar if is_tar else load_dir

    # Check for a config file, and load it. We reassign things in the dict for speed reasons.
    if check_file(decode_input,'config.json'):
        if is_tar:
            decoder_ring = json.loads(decode_input.extractfile('config.json'))
        else:
            with open(os.path.join(decode_input,'config.json'),'rb') as f: decoder_ring = json.loads(f.read())
    else:
        print 'ERROR: No config.json file was found in your input path! I cant decode data without it!'; exit()
    bases = decoder_ring['bases']
    N_qual = decoder_ring['N_qual']
    prefix = decoder_ring['prefix']
    suffix = decoder_ring['suffix']
    dna_max = decoder_ring['dna_max']
    qualities = decoder_ring['qualities']
    columnType = decoder_ring['columnType']
    separators = decoder_ring['separators']
    bits_per_base = decoder_ring['bits_per_base']
    bits_per_quality = decoder_ring['bits_per_quality']
    total_dna_bits = dna_max*bits_per_base
    total_qual_bits = dna_max*bits_per_quality
    qual_N = dict( (v,k) for k,v in N_qual.iteritems()) # reverse key/val pairs

    # DNA/QUAL decoder:
    # This function takes either a DNA or QUAL array, and for each row sums all the uint8 columns together in such a way that 
    # its as if we had stored a single number in uintX - where X is 8*[the number of columns]. There is no uint152, but with this
    # there can be using 19 uint8 columns (for example). Then, on that single number, we can bit shift around to get out values. 
    # Note that much like the encoder, this would be waaaaaay faster if it was written in a language that could read/write/bit 
    # shift arbitary bytes of binary without this dumb python "trick". Maybe numpy's void type can do it? Who knows.
    # For each row, it returns a list of numbers, where each number is the encoded base/quality score.
    def split_bits(numpy_array,total_bits,bits_per_x,decoder):
        bitmask = int('1'*bits_per_x,2)
        static = range(total_bits-bits_per_x,-bits_per_x,-bits_per_x)
        for row in numpy_array:
            the_number = sum([int(y) << (8*x) for x,y in enumerate(reversed(row))])
            yield [(the_number >> x) & bitmask for x in static]

    # QNAME decoder:
    # This function is a little special, because unlike the above 'static' function, this function, for speed reasons,
    # will be compiled at runtime using exec(). You give it your QNAME prefix/suffix/coding table, and the qname numpy array,
    # and it returns the decoded QNAMEs.
    convert_qname_fuction = '''
def convert_qname(numpy_array,prefix,suffix,columnType):
    for row in numpy_array:
        yield "'''+prefix+'" + '
    for idx,column in enumerate(columnType):
        if column[0] == 'coded':
            convert_qname_fuction += 'columnType['+str(idx)+'][1][row['+str(idx)+']] + '
        if column[0] == 'ints':
            convert_qname_fuction += 'str(row['+str(idx)+']+columnType['+str(idx)+'][1]) + '
        if column[0] == 'string':
            print 'ERROR: I dont support string encoding yet.'
        try: convert_qname_fuction += '"'+separators[idx]+'" + '
        except: pass # The last column has no separator
    convert_qname_fuction += '"'+suffix+'"'
    exec(convert_qname_fuction) # This is much more efficient than having to parse columnType for every QNAME entry. We compile a decoder once and reuse it.




    if check_file(decode_input,'DNA.raw'):                                      DNA = load_file(decode_input,'DNA.raw')
    elif check_file(decode_input,'DNA') and check_file(decode_input,'DNA.key'): DNA = load_file(decode_input,'DNA')[load_file(decode_input,'DNA.key')]

    if check_file(decode_input,'QUAL.raw'):                                       QUAL = load_file(decode_input,'QUAL.raw')
    elif check_file(decode_input,'QUAL') and check_file(decode_input,'QUAL.key'): QUAL = load_file(decode_input,'QUAL')[load_file(decode_input,'QUAL.key')]

    if check_file(decode_input,'QNAME.raw'):                                        QNAME = convert_qname(load_file(decode_input,'QNAME.raw')                                 ,prefix,suffix,columnType)
    elif check_file(decode_input,'QNAME') and check_file(decode_input,'QNAME.key'): QNAME = convert_qname(load_file(decode_input,'QNAME')[load_file(decode_input,'QNAME.key')],prefix,suffix,columnType)
    elif check_file(decode_input,'QNAME.key'):                                      QNAME = load_file(decode_input,'QNAME.key') # Pairs still mated
    else:                                                                           QNAME = xrange(1,len(DNA)+1) # Every read is unique.

    try:
        gen = itertools.izip(split_bits(DNA,total_dna_bits,bits_per_base,bases),split_bits(QUAL,total_qual_bits,bits_per_quality,qualities),QNAME)
        if len(N_qual) != 0:
            for dna,qual,qname in gen:
                # Here we decode the encoded value into the base/quality value, making sure to convert the N_qual values too.
                for idx,q in enumerate(qual):
                    qual[idx] = qualities[q]
                    if q in qual_N: dna[idx] = qual_N[q]
                    else:           dna[idx] = bases[dna[idx]]

                print qname
                print ''.join(dna)
                print '+'            # I'm seeing the qname here in a lot of publicly avalible files instead of +. Perhaps could add a flag to config.json?
                print ''.join(qual)

        else:
            # We can decode everything a little bit faster if theres nothing in N_qual/qual_N to convert.
            for dna,qual,qname in gen:
                for idx,q in enumerate(qual):
                    qual[idx] = qualities[q]
                    dna[idx] = bases[dna[idx]]
                print qname
                print ''.join(dna)
                print '+'
                print ''.join(qual)
    except IOError:
        pass




'''
Developer To Dos / Notes:

1) Variable-length encoded into the qual scores.

2) get --order working

3) To reduce the final output filesize further, i would like to write all reads containing Ns to a separate file, then once the unique/sorted
   list is created (for DNA/QUAL), try replacing the N with any other letter to get a fragment that is already in the unique dataset. Since the list
   is sorted, we can use a binary intersect and figure that out near-instantaniously, so it's reasonable. If no match is found, replace Ns with the most
   common base, which is what is currently done. This is part of the LOSSLESS part of the program, so really it would benefit all.

4) The analysing/encoding of the data is 3-10x faster in pypy than regular python, however,
   pypy does not support numpy like CPython does. Specifically the numpy.unique and numpy.fromfile / numpy.load
   functions are not implimented fully. numpy.save is, however.
   Having said that, ideally we would want to drop the numpy dependancy all together. Although numpy is a huge project
   I don't know if save/load/unique/rot90 will all work in 10-20 years from now. Python 2 and 3 intepreters, however
   definitely will. For that reason, I would like to impliment in pure python the unique and rot90 functions.
   If the user uses pypy, this script will encode/decode much faster - and of course can still be decoded by standard python.
   But really, if people like the project, would be re-written in C, parallelized, etc. Then again. C rarely compiles first time.

5) Could ACGTrie all this DNA, then the only change would be the idx points to the end row in the trie (and working backwards gets you the DNA)
   Same could work for the QUAL too.

6) You can't rotated the QNAMES John - it'll be a mix of uints and text...?

7) pypy is really really slow at doing rotations (hours) compared to cpython (seconds). Perhaps rotation is unessecary and only writing/appending data via index is needed to acheive same result.
   If using C struct this is irrelevent anyway.  

8) There are no repeated QNAMEs (or few) so key method isn't working so well. Should make a key per column after sort/unique'ing columns individually, then even sort/unique the stack of those keys. Or all keys. Or something.

'''
