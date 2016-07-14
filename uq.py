#!/usr/bin/env python

import os
import re
import sys
import json
import time
import numpy
import struct
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

On the topic of safety, microq.py is one of the only FASTQ encoding/compression programs that reads the whole
file first before deciding on how to encode the data optimally/safely. This makes it a little slower than
other tools, but argubly much safer.

But it's important to remember that this isn't just about storing ASCII bytes in the smallest number of bits. 
It's about storing those bits in such a way that compressors can very efficiently compress the data - whatever 
compressor you choose to use. I use the very effient lpac9 algorithum for my own data, although other algorithums
from the same fantastic developer (Matt Mahoney) can get the files even smaller if you have more cpu cycles to spare.
For distribution however, I would use something more prevelent like gzip or 7zip, as that's what your users will
already have installed.

The decoding of .uq files to .fastq can also be done with this tool, and this tool is included in the .uq file itself,
which is really just a tar file holding some binary blobs. To decode the data, run this tool will --decode, and 
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
parser.add_argument("--rotate",       nargs='+', metavar='file name',     help='Optional. Rotate an output struct 90 degrees before saving. Helps with some compressors. Will be un-rotated automatically. Typing the name twice rotates twice. [DNA/QUAL/QNAME]')
parser.add_argument("--raw",          nargs='+', metavar='file name',     help='Optional. Saves data as a raw table rather than unique/key combo. Helps with some compressors. [DNA/QUAL/QNAME]')
parser.add_argument("--temp",                                             help='Optional. Directory to write small temporary files to.')
parser.add_argument("--dna",                                              help='Optional. [drop/hide] - see readme/code for details.')
parser.add_argument("--qual",                                             help='Optional. [nearest/weighted] [N] - see readme/code for details.')
parser.add_argument("--qname",                                            help='Optional. [unique/duplicate] - see readme/code for details.')
parser.add_argument("--reorder",      action='store_true', default=False, help='Optional. No parameters - see readme/code for details.')
parser.add_argument("--test",         action='store_true', default=False, help='Optional. Writes a uq tar file with every possible data output orientation. Compress them all with your compressor of choice to find the appropriate values for --raw and --rotate.')
parser.add_argument("--pad",          action='store_true', default=False, help='Optional. Pads DNA/QUAL to the nearest 2/4/8 bits. Some compressors do a better job when data is padded.')
parser.add_argument("--peek",         action='store_true', default=False, help='Optional. No output files are created, the input is just scanned for the parameters that *would* be used.')
parser.add_argument("--paranoid",     action='store_true', default=False, help='Optional. After encoding, will attempt to decode the output and compare its MD5 to the original input FASTQ (only works for lossless encoding)')
parser.add_argument("--delete",       action='store_true', default=False, help='Optional. Input file is deleted if no issues occur during encoding. If --paranoid used, only deletes if this returns successfully.')
parser.add_argument("--decode",       action='store_true', default=False, help='Requred if you want to decode a .uq file back to .fastq :)')
args = parser.parse_args()

if args.rotate:
    for x in args.rotate:
        if x not in ['DNA','QUAL','QNAME']:
            print 'ERROR: When using --rotate, you must specificy if you wish to rotate "DNA", "QUAL" or "QNAME" tables.'
            print '       If you specify the same table more than once, it will be rotated 90 degrees more than once.'
            exit()
    rotations = collections.Counter(args.rotate)
else: rotations = collections.Counter([])
if args.raw:
    for x in args.raw:
        if x not in ['DNA','QUAL','QNAME']:
            print 'ERROR: When using --raw, you must specificy if which tables, "DNA", "QUAL" or "QNAME", you wish to keep as raw.'
            exit()
    args.raw = set(args.raw)
else: args.raw = set()


if args.qname or args.qual or args.dna or args.paranoid or args.delete: print 'Not implimented yet.'; exit()

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
        if not line1.startswith('@'):                 print 'ERROR: This does not look like a FASTA/FASTQ file! (first line does not start with @)'; exit()

        # QNAME
        prefix = line1
        suffix = line1[:]
        separators = set()
        for sep in set(line1): separators.add(( sep, line1.count(sep) ))

        # DNA
        dna_bases = set(line2)
        dna_min = len(line2)
        dna_max = len(line2)

        # +
        '(no checks currently done, no data collected, always decodes to a +)'

        # QUAL
        quals = set(line4) 

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
                next(f) # +
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
                    if character != qname[idx]: prefix = prefix[:idx]; break
            if not qname.endswith(suffix):
                for idx,character in enumerate(suffix[::-1]):
                    if character != qname[-1-idx]:
                        suffix = suffix[-idx:]; break
            for sep,sepcount in separators.copy():
                if qname.count(sep) != sepcount: separators.remove( (sep,sepcount) )

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

    # Tidy up separator/prefix/suffix results:
    separators = list(separators)
    for idx,sepdata in enumerate(separators):
        sep,sepcount = sepdata
        if prefix.count(sep) != 0: separators[idx] = ( sep, sepcount - prefix.count(sep) )
        if suffix.count(sep) != 0: separators[idx] = ( sep, sepcount - suffix.count(sep) )
    for idx in range(len(separators)-1,-1,-1): # iterate the list backwards so we can delete things without upsetting the order.
        if separators[idx][1] == 0: separators.pop(idx)

    def order_seps(string,separators):
        sep_list = []
        for sep in separators:
            last = 0
            for _ in range(sep[1]):
                sep_list.append( (qname.index(sep[0],last),sep[0]) )
                last = qname.index(sep[0],last)+1
        return ''.join([y for x,y in sorted(sep_list)])

    if order_seps(qname,separators) == order_seps(line1,separators):
        separators = order_seps(qname,separators)
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
    if len(separators): print '    - QNAME field separators:'  ,separators
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

        # Some numpy functions work on all elements of a multi-dimentional array by first flattening it (for example, diff)
        # Others require that rows be 'grouped' and treated as single elements, like unique.
        # We switch between the two formats (near-instantaniously) with the following two functions:
        def struct_to_std(ar):
            length = len(ar)
            columns = len(ar[0])
            ar.dtype = 'uint8'
            ar.shape = (length,columns)
        def std_to_struct(ar):
            length,columns = ar.shape
            ar.dtype = ','.join(columns*['uint8'])
            ar.shape = length

        def write_out(table,filename,rotations):
            with open(os.path.join(temp_directory,filename),'wb') as out_file:
                numpy.savez(out_file,table=table,rotations=rotations)
        def read_in(filename): return numpy.load(os.path.join(temp_directory,filename))

        write_out(dna_array,'DNA.raw',0); del dna_array # Free up some RAM for the impending sort. Will process it later.
        if args.test:
            struct_to_std(qual_array)
            write_out(numpy.rot90(qual_array,0),'QUAL.raw',0)
            write_out(numpy.rot90(qual_array,1),'QUAL.raw.rotate1',1)
            write_out(numpy.rot90(qual_array,2),'QUAL.raw.rotate2',2)
            write_out(numpy.rot90(qual_array,3),'QUAL.raw.rotate3',3)
            std_to_struct(qual_array)

        if 'QUAL' in args.raw and not args.test:
            struct_to_std(qual_array)
            write_out(numpy.rot90(qual_array,rotations['QUAL']),'QUAL.raw',rotations['QUAL'])
            del qual_array
        else:
            qual_array,qual_key = numpy.unique(qual_array, return_inverse=True)
            if   len(qual_array) <= 256:                  qual_key_dtype = 'uint8'   #; qual_key = qual_key.astype(qual_key_dtype)
            elif len(qual_array) <= 65536:                qual_key_dtype = 'uint16'  #; qual_key = qual_key.astype(qual_key_dtype)
            elif len(qual_array) <= 4294967296:           qual_key_dtype = 'uint32'  #; qual_key = qual_key.astype(qual_key_dtype)
            elif len(qual_array) <= 18446744073709551616: qual_key_dtype = 'uint64'  #; qual_key = qual_key.astype(qual_key_dtype)
            write_out(qual_key,'QUAL.key',0); del qual_key
            struct_to_std(qual_array)
            if args.test:
                write_out(numpy.rot90(qual_array,0),'QUAL',0)
                write_out(numpy.rot90(qual_array,1),'QUAL.rotate1',1)
                write_out(numpy.rot90(qual_array,2),'QUAL.rotate2',2)
                write_out(numpy.rot90(qual_array,3),'QUAL.rotate3',3)
            else:
                write_out(numpy.rot90(qual_array,rotations['QUAL']),'QUAL',rotations['QUAL'])
            del qual_array

        print 'Finished sorting qualities',status.split_time(),'sorting DNA...'

        dna_array = read_in('DNA.raw')['table']
        os.remove(os.path.join(temp_directory,'DNA.raw'))  # I apprechiate this is dumb if we're going to write it out again, but need to allow for rotation.
        if args.test:
            struct_to_std(dna_array)
            write_out(numpy.rot90(dna_array,0),'DNA.raw',0)
            write_out(numpy.rot90(dna_array,1),'DNA.raw.rotate1',1)
            write_out(numpy.rot90(dna_array,2),'DNA.raw.rotate2',2)
            write_out(numpy.rot90(dna_array,3),'DNA.raw.rotate3',3)
            std_to_struct(dna_array)

        if 'DNA' in args.raw and not args.test:
            struct_to_std(dna_array)
            write_out(numpy.rot90(dna_array,rotations['DNA']),'DNA.raw',rotations['DNA'])
            del dna_array
        else:
            dna_array,dna_key = numpy.unique(dna_array, return_inverse=True)
            if   len(dna_array) <= 256:                  dna_key_dtype = 'uint8'   #; dna_key = dna_key.astype(dna_key_dtype)
            elif len(dna_array) <= 65536:                dna_key_dtype = 'uint16'  #; dna_key = dna_key.astype(dna_key_dtype)
            elif len(dna_array) <= 4294967296:           dna_key_dtype = 'uint32'  #; dna_key = dna_key.astype(dna_key_dtype)
            elif len(dna_array) <= 18446744073709551616: dna_key_dtype = 'uint64'  #; dna_key = dna_key.astype(dna_key_dtype)
            write_out(dna_key,'DNA.key',0); del dna_key
            struct_to_std(dna_array)
            if args.test:
                write_out(numpy.rot90(dna_array,0),'DNA',0)
                write_out(numpy.rot90(dna_array,1),'DNA.rotate1',1)
                write_out(numpy.rot90(dna_array,2),'DNA.rotate2',2)
                write_out(numpy.rot90(dna_array,3),'DNA.rotate3',3)
            else:
                write_out(numpy.rot90(dna_array,rotations['DNA']),'DNA',rotations['DNA'])
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
        columnType = []
        '''
        columnType explination:
        Columns start off as 'coded' which means a tuple of all known values will be used to assign "things" to numbers, just like the DNA/QUAL.
        Columns may get bumped down to 'ints' or 'strings' type if this tuple gets very large - specifically if the size of the tuple contains 
        more than 1000 items, and whatever that number above 1000 is, it is more than 1/10th of the number of sequences read so far.
        Next down the type list is ints, which are always numbers. Then below that is strings. It's important to give ints the possibility of being
        tuple-coded, rather than find a struct big enough to hold them as numbers, because very large numbers may appear rarely, and thus coding is better.
        '''
        qnames = qname_reader(args.input,prefix,suffix)
        for counter,qname in enumerate(qnames):
            if len(qname) != cols:
                if cols == False:
                    cols = len(qname)
                    for _ in range(cols): columnType.append(['coded',set()])
                else:
                    print 'WARNING: I was unable to correctly guess the delimiters of the QNAME in this FASTA file.'
                    print '         Will fall back to assuming there is no delimiting of data.'
                    fallback = True
                    break
            for idx,column in enumerate(qname):
                if columnType[idx][0] == 'coded': columnType[idx][1].add(column)
                elif columnType[idx][0] == 'ints':
                    try:
                        column = int(column)
                        if column < columnType[idx][1]: columnType[idx][1] = column   # min
                        elif column > columnType[idx][2]: columnType[idx][2] = column # max
                    except ValueError: 
                        # Get bumped down to a string
                        columnType[idx][0] = 'strings'
                        columnType[idx][1] = len(str(columnType[idx][2]))
                        del columnType[idx][2]
                # To be here we have to be a string:
                elif len(column) > columnType[idx][1]:
                    columnType[idx][1] = len(column)

            if counter == target:
                for column_description in columnType:
                    if column_description[0] == 'coded':
                        if len(column_description[1]) > target/10:
                            try:
                                column_description[1] = map(int,column_description[1])
                                column_description[0] = 'ints'
                                column_description.append(max(column_description[1]))
                                column_description[1] = min(column_description[1])
                            except ValueError:
                                column_description[0] = 'strings'
                                column_description[1] = max(map(len,yolo))
                target *= 10

        if fallback:
            print 'I havent implimented this yet!'
            exit()
        else:
            # Check the last bit:
            for column_description in columnType:
                if column_description[0] == 'coded':
                    if len(column_description[1]) > target/10:
                        try:
                            column_description[1] = map(int,column_description[1])
                            column_description[0] = 'ints'
                            column_description.append(max(column_description[1]))
                            column_description[1] = min(column_description[1])
                        except ValueError:
                            column_description[0] = 'strings'
                            column_description[1] = max(map(len,yolo))

        print 'Optimal encoding method for delimited data in QNAME determined!',status.split_time()

        qname_dtype = []
        for idx,column_description in enumerate(columnType):
            # Now we fiddle with things a bit:
            if column_description[0] == 'coded':
                i = len(column_description[1])
                if   i <= 255:                  column_description.append('uint8')
                elif i <= 65535:                column_description.append('uint16')
                elif i <= 4294967295:           column_description.append('uint32')
                elif i <= 18446744073709551615: column_description.append('uint64')
                column_description[1] = tuple(sorted(column_description[1])) # sets have no order, so we tuple to fix things in place.
            elif column_description[0] == 'ints':
                i = column_description[2]-column_description[1] +1
                if   i <= 255:                  column_description[2] = 'uint8'
                elif i <= 65535:                column_description[2] = 'uint16'
                elif i <= 4294967295:           column_description[2] = 'uint32'
                elif i <= 18446744073709551615: column_description[2] = 'uint64'
            elif column_description[0] == 'strings':
                print 'I havent implimented this yet'; exit()
            print '    - Column',idx+1,'is type',column_description[0],'stored as',column_description[2]
            qname_dtype.append(column_description[2])
        max_dtype = 'uint' + str(max([ int(x[4:]) for x in qname_dtype])) # The largest dtype in all the QNAME columns.
        qname_dtype = ','.join(qname_dtype)
    else: print 'Step 3 of 4: Skipped as there are no delimiters common to all QNAMEs'

    decoder_ring = {
        'base_distribution':        base_graph,
        'qual_distribution':        qual_graph,
        'prefix':                   prefix,
        'suffix':                   suffix,
        'separators':               separators,
        'bases':                    dna_bases,
        'qualities':                quals,
        'variable_read_lengths':    variable_read_lengths,
        'bits_per_base':            bits_per_base,
        'bits_per_quality':         bits_per_quality,
        'N_qual':                   N_qual,
        'dna_max':                  dna_max,
        'columnType':               columnType
    }

    if args.peek:
        print 'The config.json would look like:'
        print json.dumps(decoder_ring,indent=4)

    else:
        print 'Writing config to disk...',
        path = os.path.join(temp_directory,'config.json')
        with open(path,'wb') as f:
            f.write(json.dumps(decoder_ring,indent=4))
        print 'success!\n'

        status.message = 'Step 4 of 4: Encoding QNAME with optimized settings.'
        qname_array = numpy.zeros(status.total,numpy.dtype(qname_dtype))
        rows = len(qname_array)
        columns = len(qname_array.dtype)
        qnames = qname_reader(args.input,prefix,suffix)
        encoder = 'for row,qname in enumerate(qnames):\n    qname_array[row] = ('
        code = []
        for idx,column_description in enumerate(columnType):
            if column_description[0] == 'coded':     code.append(' columnType['+str(idx)+'][1].index(qname['+str(idx)+'])')
            elif column_description[0] == 'ints':    code.append(' int(qname['+str(idx)+'])-'+str(column_description[1])  )
            elif column_description[0] == 'strings': code.append(' qname['+str(idx)+']'                                   )
        encoder += ','.join(code) + ' )'
        exec(encoder) # dirty but it works well.


        if args.test:
            write_out(qname_array,'QNAME.raw',0) # The raw array with the smallest dtype
            temp_array = numpy.asarray(qname_array,','.join([max_dtype]*columns)) # Now all columns have the dtype to the largest. Required for rotation.
            temp_array.dtype = max_dtype
            temp_array.shape = len(qname_array),len(qname_array[0])
            write_out(numpy.rot90(temp_array,1),'QNAME.raw.rotate1',1)
            write_out(numpy.rot90(temp_array,2),'QNAME.raw.rotate2',2)
            write_out(numpy.rot90(temp_array,3),'QNAME.raw.rotate3',3)
            del temp_array

        if 'QNAME' in args.raw and not args.test:
            if 'QNAME' in rotations:
                qname_array = numpy.asarray(qname_array,','.join([max_dtype]*columns)) # Now all columns have the dtype of the largest column. Required for rotation.
                qname_array.dtype = max_dtype
                qname_array.shape = columns,rows
                write_out(numpy.rot90(qname_array,rotations['QNAME']),'QNAME.raw',rotations['QNAME'])
            else:
                write_out(qname_array,'QNAME.raw',0)
        else:
            beforeLength = len(qname_array)
            qname_array,qname_key = numpy.unique(qname_array, return_inverse=True)
            if   len(qname_array) <= 256:                  qname_key_dtype = 'uint8';  qname_key = qname_key.astype(qname_key_dtype)
            elif len(qname_array) <= 65536:                qname_key_dtype = 'uint16'; qname_key = qname_key.astype(qname_key_dtype)
            elif len(qname_array) <= 4294967296:           qname_key_dtype = 'uint32'; qname_key = qname_key.astype(qname_key_dtype)
            elif len(qname_array) <= 18446744073709551616: qname_key_dtype = 'uint64'; qname_key = qname_key.astype(qname_key_dtype)
            if len(qname_array) != beforeLength:
                if len(qname_array)*2 != beforeLength:
                    print 'INFO: Of the',beforeLength,'reads in the file,',beforeLength-len(qname_array),'contained duplicated read names/headers/QNAMEs.'
                    print '      What is unusual is that this number is not half the file (which might make sense for paired-end sequencing).'
                    print '      Thus, there may be a mistake in your QNAME formatting, however, it is not my place to judge - encoding will be totally uneffected.'
                else: 
                    print 'INFO: This file looks like it contains properly-formatted paired end data.'
            else:
                print 'INFO: This file looks like it contains single-end sequencing data.'
                print '      If that is NOT the case, then there may be a mistake in your QNAME formatting. Otherwise everything is fine :)'
            write_out(qname_key,'QNAME.key',0); del qname_key
            if args.test:
                write_out(qname_array,'QNAME',0) # This array will have the smallest dtype for each column possible
                qname_array = numpy.asarray(qname_array,','.join([max_dtype]*columns)) # Now all columns have the dtype of the largest column. Required for rotation.
                qname_array.dtype = max_dtype
                qname_array.shape = columns,rows
                write_out(numpy.rot90(qname_array,1),'QNAME.rotate1',1)
                write_out(numpy.rot90(qname_array,2),'QNAME.rotate2',2)
                write_out(numpy.rot90(qname_array,3),'QNAME.rotate3',3)
            else:
                if 'QNAME' in rotations:
                    qname_array = numpy.asarray(qname_array,','.join([max_dtype]*columns)) # Now all columns have the dtype to the largest. Required for rotation.
                    qname_array.dtype = max_dtype
                    qname_array.shape = columns,rows
                    write_out(numpy.rot90(qname_array,rotations['QNAME']),'QNAME',rotations['QNAME'])
                else:
                    write_out(qname_array,'QNAME',0)
        del qname_array
        print 'Encoded QNAMEs written to disk!',status.split_time()

        if args.reorder:
            print 'Bonus Step: Reordering reads to optimize compression.'
            master = numpy.stack([ read_in('DNA.key'), read_in('QNAME.key'), read_in('QUAL.key') ], axis=1)
            master = master[master[:,2].argsort(kind='quicksort')]
            master = master[master[:,1].argsort(kind='mergesort')]
            master = master[master[:,0].argsort(kind='mergesort')]
            write_out(numpy.rot90(master),'master.key',1)

        try:
            with tarfile.open(os.path.join(temp_directory,'temp.uq'), mode='w') as temp_out:
                for f in os.listdir(temp_directory): temp_out.add(os.path.join(temp_directory,f),arcname='.')
            os.rename(os.path.join(temp_directory,'temp.uq'), args.output)
            for f in os.listdir(temp_directory): os.remove(os.path.join(temp_directory,f))
        except Exception as e:
            print 'ERROR: There was an error taking the encoded data and putting it all in a single tar file:'
            print '      ',e
            print '       You might be able to still do it manually. Take a look inside',temp_directory

if args.decode or args.paranoid:
    if args.decode:
        decode_input = tarfile.open(args.input) if is_tar else args.input    # for decode
    else:
        decode_input = tarfile.open(args.output) if is_tar else args.output    # for paranoid

    # Define a few very basic functions:
    def check_dir(decode_input,file_name): return os.path.isfile(os.path.join(decode_input,file_name))
    def check_tar(decode_input,file_name): return file_name in decode_input.getnames()
    def load_tar(decode_input,file_name):
        data = numpy.load(decode_input.extractfile(file_name))
        if data['rotations'] == 0: return data['table']
        else:                      return numpy.rot90(data['table'],-data['rotations']) # undo the rotations. Could undo a diff here too.
    def load_dir(decode_input,file_name):
        data = numpy.load(os.path.join(decode_input,file_name))
        if data['rotations'] == 0: return data['table']
        else:                      return numpy.rot90(data['table'],-data['rotations']) # undo the rotations. Could undo a diff here too.
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

'''
