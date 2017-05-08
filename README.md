# uQ - binary FASTQ
  uQ is a lossless FASTQ binary encoder.
  
  uQ does not perform compression itself (uQ files are just numpy arrays),
  but if provided the path to a compression program, uQ can test different 
  settings to optimise the binary encoding for your specific compressor and
  your specific data. This means one can use whatever compression algorithm
  they like, and be confident they are getting a very small output.

  LZMA+uQ will typically give the smallest filesize of any FASTQ compressor,
  particularly for files with long/complex QNAMEs (as this data is delimited,
  encoded, and stored in a table, which aids access and compression).

# Why use uQ?

1. Being a binary format, not only does a uQ file take up considerably less space than a FASTQ file before compression (meaning one can store many of them in memory at once), but a uQ file will often compress down to half the size of it's FASTQ equivilent **whatever compression program you use**. This is ideal if you wish to use a popular compression tool for long-term storage.

2. Being a binary format, data in a uQ file can be read and operated upon much faster than a FASTQ file can be. The data is really nothing more than a collection of Numpy arrays (in NPY format) all tar'd together into a single file (with a json config file). This means one could store, for example, all the DNA sequences from the ENCODE project in memory, on modest hardware.

# How to use uQ:

Good compression depends on the compressor, the amount of information stored in the file, and how that information is physically arranged on disk before compression. As such, there are plenty of non-linear variables at play, and the only way to know what is best is to try everything (or have a good guess based on some similar data). The variables uQ currently accepts include:

### --sort:
```
  Re-sorts the FASTQ file to reduce file size. Values can be "DNA", "QUAL", "QNAME" and "None", 
  case insensitive. This will change the final output (it will be sorted), however if sorting 
  by "DNA" it is likely downstream mapping will be much faster, although "QUAL" usually gives 
  the smallest files.
```

### --raw:
```
  By default tables in uQ format are sorted and unique, with an index used to pull out the 
  values in the correct order.

  [ A, C, C, C, G, T, C, C, A ]   ----\   [ A, C, C, C, G, T, C, C, A ] 
  [ A, C, C, C, G, T, C, C, A ]        )  [ G, G, G, G, G, T, C, C, A ]  +  [1, 1, 3, 2]
  [ T, T, C, C, G, T, C, C, A ]   ----/   [ T, T, C, C, G, T, C, C, A ]
  [ G, G, G, G, G, T, C, C, A ]

  Thus, --raw prevents this behaviour, and instead the data is stored as a single table with no index.
  Sorting still depends on --sort. How desirable this is depends on the nature of your data 
  (how repetitive it is), and what in-memory analyses you are planning to do. Valid values are 
  "DNA", "QUAL", "QNAME", and "None", and more than one value can be specified.
```

### --pattern: 
  ```
  While --raw changes the logical layout of the data, --pattern changes the physical layout,
  i.e. the way in which bytes are physically written to disk, such as writing/reading all 
  rows backwards, writing data in column order and not row order, etc. Because there are 
  so many combinations of transformations, these are simply reffered to by ids: "0.1", "0.2", 
  "1.1", "1.2", "2.1", "2.2", "3.1", and "3.2". No one is expected to remember what these 
  numbers mean, because...
  ```

### --test:
```
  ...when --test is specified, all possible values for the above settings are tried 
  (excluding any that you have locked-in by specifying them explicitly) and the
  best possible combination for your compressor is written out, with details on the less-good 
  arrangments printed to the terminal. This can be done on a small sample of the data 
  if you wish (for speed), but getting a random-enough sample is left to the user for now. 
  To be totally honest, I usually just run --test on the whole file for one type of data set,
  say an H3K9me3 ChIP, which after some time returns the best possible file encoding, and 
  then I use those same --pattern/--raw/--sort values for all other H3K9me3 datasets.
  Do make sure to test every sample kind individually though, as results generated for an 
  RNA-Seq will be different to ChIP-Seq which is different to WGBS, etc.
```

### --compressor:
```
  How does uQ know which arrangements of the above options compresses the best? It compresses the data 
  with a compressor of your choosing with the --compressor option. This is the path or alias, i.e. "gzip".
  Compressor is simply subprocessed, data is fed to the stdin and it's stdout is piped to "wc -c" to 
  count bytes. The code can be modified to work with compression programs that have to read/write files 
  and not from stdin by uncommenting a few lines. Note that if you do not specify a --compressor, uQ 
  will try to optimise for smallest uncompressed filesize, and the results will be very different to
  what a compressor would likely give.
```

### --temp:
```
  At present, data has to be written to a temporary directory after encoding to reduce memory load. 
  This is handled by the python tempfile module, which uses the system preference (usually a random 
  directory in /tmp), however you can override that here.
```

### --pad:
```
  Some compressors prefer it when "symbols" are byte-aligned, i.e. some factor of 8 (2/4/8/16/32). 
  --pad does this. This will definitely increase the size of the uncompressed data, but may or may 
  not decrease the size of the file after compression.
```

### --notricks:
```
  In FASTQ, "N" bases also have to have quality scores. While this does not make much sense, some 
  FASTQs have all Ns with a special, unique quality score. This is enough to identify Ns, and thus
  Ns in the sequence can be temporarily replaced with As (or rather, the most common base), for good 
  compression. This is significant, as it takes 2 bits to store A/C/G/T but 3bits to store A/C/G/T/N 
  - a 50% increase in file size! Typically the quality encodings have bits to spare, as there are so
  many quality scores. For this reason its beneficial to create a new, unique quality score for Ns,
  where possible (has to be only 1 quality score for the N, otherwise it doesn't look like an N to uQ,
  as uQ will encode anything FASTQish reliably). Long story short, --notricks will prevent this trickery.
  This will benefit users not looking for optimum file size reduction, but 'cleaner' datasets that can 
  be operated on plainly. Note, the tricks are used by default, but it will not change the file after 
  decoding.
```

### --peek:
```
  This flag essentially makes uQ abort after analysing the file but before encoding starts.
```

### --input:
```
  FASTQ files in when encoding, uQ files in when decoding. Unfortunately, this will have to be a 
  file and not a pipe, as the file has to be read multiple times to tune the encoder to the data
  and prevent issues that might result in an undecodable file 10 years from now.
```

### --output: 
```
  Path to output. If no output specified, default when making uQ files is to append ".uQ" to the end
  of the original file name, and when decompressing uQ files back to FASTQ, default is to print to
  the standard out/terminal.
```

### --decode:
```
  This simply reverses the encoding operation. No parameters need to be specified, it's all stored in the
  config file in the uQ tar file. I have somehow manged to write a decoder that is slower than the encoder,
  which is quite impressive. I'm sure if there is any interest in uQ files, the speed will be improved over time.
```

# How to *(actually)* use uQ:

If you wish to encode a very large file and have no idea what sort/raw/pattern to use, you can --test on a smaller sample:
head -400000 ./rep3.fastq > ./sample.fastq
uq.py --input sample.fastq --test

Outputs:
```
(pypy) Johns-MacBook:Desktop John$ python uq.py --input /Users/John/Desktop/rep3.fastq --test --compressor pigz
Temporary directory: /var/folders/3t/qctkm9xn6zv871cvwy6pw31m0000gn/T/tmpfEggaT
Pass 1 of 4: First pass through file to find optimal/safe encoding strategy. (99.93%)
Finished analysing file! (0.743970485528 minutes)
    - total reads: 9852944
    - total bases: 354705984

QNAME Analysis:
    - all QNAMEs prefixed with: @SL-XBG_1_FC30CT9AAXX:8:
    - QNAME field separators: :: ( 2 symbols )

DNA Analysis:
    - DNA sequences contained the following characters: A C G N T (5 in total)
    - There distribution is:
      A |############                                      |  25.406%     90116798
      C |#############                                     |  26.154%     92768093
      G |#########                                         |  19.033%     67512877
      N |#                                                 |   2.535%      8990390
      T |#############                                     |  26.872%     95317826
    - this means we will store each letter of DNA in 3 bits.
    - all DNA sequences are 36 bases long
    - therefore, as 108 bits must be used to store the DNA for each sequence, we will use 14 bytes per unique DNA sequence.

QUAL Analysis:
  [Breakdown for "G"]
      $ |                                                  |   0.456%       308036
      % |#                                                 |   3.385%      2285609
      & |##                                                |   5.083%      3431647
      ( |#                                                 |   2.932%      1979456
      ) |####                                              |   8.007%      5405472
      * |                                                  |   0.203%       136927
      + |#####                                             |  10.314%      6963235
      , |                                                  |   1.457%       983831
      - |                                                  |   0.127%        85592
      . |###                                               |   7.962%      5375579
      / |#                                                 |   3.517%      2374618
      0 |                                                  |   0.000%          291
      1 |#                                                 |   3.842%      2594021
      2 |                                                  |   0.478%       322432
      3 |####                                              |   8.669%      5852521
      4 |                                                  |   1.359%       917338
      5 |#                                                 |   2.075%      1401087
      6 |##                                                |   4.676%      3156908
      7 |###                                               |   6.382%      4308856
      8 |                                                  |   1.875%      1265747
      9 |#                                                 |   2.552%      1722737
      : |###                                               |   6.490%      4381867
      ; |                                                  |   1.755%      1184962
      < |########                                          |  16.308%     11010049
      = |                                                  |   0.095%        64030
      A |                                                  |   0.000%           29

  [Breakdown for "A"]
      $ |                                                  |   0.408%       367726
      % |#                                                 |   3.083%      2778268
      & |##                                                |   5.094%      4590233
      ( |#                                                 |   3.022%      2723077
      ) |####                                              |   8.291%      7471607
      * |                                                  |   0.225%       203200
      + |#######                                           |  15.059%     13570635
      , |                                                  |   1.178%      1061225
      - |                                                  |   0.097%        87742
      . |####                                              |   9.422%      8490983
      / |###                                               |   7.844%      7068486
      0 |                                                  |   0.000%          447
      1 |#                                                 |   2.915%      2627267
      2 |                                                  |   0.226%       203328
      3 |######                                            |  12.108%     10910930
      4 |                                                  |   0.501%       451869
      5 |#                                                 |   2.597%      2339894
      6 |#                                                 |   3.897%      3512185
      7 |####                                              |   8.511%      7670236
      8 |                                                  |   0.729%       657133
      9 |                                                  |   1.883%      1697209
      : |##                                                |   4.383%      3950221
      ; |                                                  |   1.072%       965744
      < |###                                               |   7.382%      6652453
      = |                                                  |   0.071%        63793
      A |                                                  |   0.001%          907

  [Breakdown for "T"]
      $ |                                                  |   0.415%       395230
      % |#                                                 |   2.907%      2770617
      & |##                                                |   4.588%      4373481
      ( |#                                                 |   2.695%      2569069
      ) |###                                               |   7.507%      7155051
      * |                                                  |   0.220%       209925
      + |#####                                             |  10.418%      9930407
      , |                                                  |   1.395%      1329887
      - |                                                  |   0.227%       216063
      . |###                                               |   7.507%      7155625
      / |##                                                |   4.359%      4154636
      0 |                                                  |   0.002%         1474
      1 |#                                                 |   3.526%      3361251
      2 |                                                  |   0.517%       493228
      3 |####                                              |   9.858%      9396297
      4 |                                                  |   1.277%      1217337
      5 |#                                                 |   2.652%      2527413
      6 |#                                                 |   3.659%      3487553
      7 |####                                              |   8.204%      7819513
      8 |                                                  |   0.982%       935838
      9 |                                                  |   1.985%      1891951
      : |###                                               |   7.766%      7402697
      ; |#                                                 |   2.828%      2695321
      < |#######                                           |  14.384%     13710745
      = |                                                  |   0.123%       116930
      A |                                                  |   0.000%          287

  [Breakdown for "C"]
      $ |                                                  |   0.418%       387846
      % |#                                                 |   3.597%      3336522
      & |##                                                |   5.812%      5392143
      ( |#                                                 |   3.538%      3282477
      ) |####                                              |   8.387%      7780511
      * |                                                  |   0.288%       267131
      + |######                                            |  13.802%     12803864
      , |                                                  |   1.584%      1469195
      - |                                                  |   0.204%       189194
      . |####                                              |   8.548%      7930263
      / |###                                               |   6.604%      6126651
      0 |                                                  |   0.001%          756
      1 |#                                                 |   2.782%      2580984
      2 |                                                  |   0.392%       363830
      3 |####                                              |   9.969%      9247691
      4 |                                                  |   0.740%       686197
      5 |#                                                 |   2.492%      2312161
      6 |#                                                 |   2.470%      2290990
      7 |###                                               |   7.746%      7185716
      8 |                                                  |   0.518%       480712
      9 |                                                  |   1.392%      1291054
      : |###                                               |   7.222%      6699892
      ; |#                                                 |   2.207%      2047474
      < |####                                              |   9.229%      8561396
      = |                                                  |   0.057%        52857
      A |                                                  |   0.001%          586

  [Breakdown for "N"]
      $ |#                                                 |   2.220%       199560
      % |################################################  |  97.780%      8790830

  [Total distribution]
    - the following values were seen as quality scores: $ % & ( ) * + , - . / 0 1 2 3 4 5 6 7 8 9 : ; < = A (26 in total)
    - There distribution is:
      $ |                                                  |   0.468%      1658398
      % |##                                                |   5.628%     19961846
      & |##                                                |   5.015%     17787504
      ( |#                                                 |   2.975%     10554079
      ) |###                                               |   7.841%     27812641
      * |                                                  |   0.230%       817183
      + |######                                            |  12.198%     43268141
      , |                                                  |   1.366%      4844138
      - |                                                  |   0.163%       578591
      . |####                                              |   8.162%     28952450
      / |##                                                |   5.561%     19724391
      0 |                                                  |   0.001%         2968
      1 |#                                                 |   3.147%     11163523
      2 |                                                  |   0.390%      1382818
      3 |####                                              |   9.982%     35407439
      4 |                                                  |   0.923%      3272741
      5 |#                                                 |   2.419%      8580555
      6 |#                                                 |   3.509%     12447636
      7 |###                                               |   7.608%     26984321
      8 |                                                  |   0.941%      3339430
      9 |                                                  |   1.862%      6602951
      : |###                                               |   6.325%     22434677
      ; |                                                  |   1.943%      6893501
      < |#####                                             |  11.259%     39934643
      = |                                                  |   0.084%       297610
      A |                                                  |   0.001%         1809
    - this means we will store each quality symbol in 5 bits.
    - therefore, as 180 bits must be used to store the qualities for each sequence, we will use 23 bytes per unique quality sequence.

Pass 2 of 4: Now QNAME delimiters are known, QNAMEs are being analysed. (99.93%)
Optimal encoding method for delimited data in QNAME determined! (0.370465902487 minutes)
    - Column 3 is type integers stored as uint8
    - Column 3 is type integers stored as uint16
    - Column 3 is type integers stored as uint16

Pass 3 of 4: Encoding DNA and Quality scores... (99.93%)
Pass 4 of 4: Encoding QNAMEs... (99.93%)

Starting tests...
Time:                       Sort:      Raw Tables:                      Patterns:
(2.3345758001 minutes)      DNA        ('DNA', 'QUAL', 'QNAME')         All
(2.91747886737 minutes)     QUAL       ('DNA', 'QUAL', 'QNAME')         All
(3.85382165114 minutes)     QNAME      ('DNA', 'QUAL', 'QNAME')         All
(2.83744401932 minutes)     None       ('DNA', 'QUAL', 'QNAME')         All
(5.93935928345 minutes)     DNA        ('DNA', 'QUAL')                  All
(5.48512871663 minutes)     QUAL       ('DNA', 'QUAL')                  All
(4.82888338168 minutes)     QNAME      ('DNA', 'QUAL')                  All
(5.01118438641 minutes)     None       ('DNA', 'QUAL')                  All
(3.38339544932 minutes)     DNA        ('QUAL', 'QNAME')                All
(3.45257088343 minutes)     QUAL       ('QUAL', 'QNAME')                All
(5.16640766462 minutes)     QNAME      ('QUAL', 'QNAME')                All
(4.31718771458 minutes)     None       ('QUAL', 'QNAME')                All
(4.29210780064 minutes)     DNA        ('DNA', 'QNAME')                 All
(4.27810624838 minutes)     QUAL       ('DNA', 'QNAME')                 All
(6.30826818546 minutes)     QNAME      ('DNA', 'QNAME')                 All
(4.31755655209 minutes)     None       ('DNA', 'QNAME')                 All
(7.60646115144 minutes)     DNA        ('DNA',)                         All
(10.3870936831 minutes)     QUAL       ('DNA',)                         All
(5.85630671581 minutes)     QNAME      ('DNA',)                         All
(6.39971796672 minutes)     None       ('DNA',)                         All
(6.96401426792 minutes)     DNA        ('QUAL',)                        All
(6.8364989837 minutes)      QUAL       ('QUAL',)                        All
(6.95274269978 minutes)     QNAME      ('QUAL',)                        All
(7.8917345643 minutes)      None       ('QUAL',)                        All
(6.10804513295 minutes)     DNA        ('QNAME',)                       All
(6.31539324919 minutes)     QUAL       ('QNAME',)                       All
(7.17378214995 minutes)     QNAME      ('QNAME',)                       All
(4.65198340019 minutes)     None       ('QNAME',)                       All
(7.19780898492 minutes)     DNA        (None,)                          All
(7.36412746906 minutes)     QUAL       (None,)                          All
(7.51234973272 minutes)     QNAME      (None,)                          All
(7.21773494879 minutes)     None       (None,)                          All

All done!
Size (compressed)    Sort:    Raw Tables:                 Raw stats:
        271821249    QUAL     ('DNA', 'QUAL', 'QNAME')    QUAL.raw : [137272420, '2.2'] DNA.raw : [95416946, '2.2'] QNAME_1.raw : 7617093 QNAME_2.raw : 15645298 QNAME_3.raw : 15869492
        274112292    QNAME    ('DNA', 'QUAL', 'QNAME')    QNAME_1.raw : 11523 QNAME_2.raw : 396943 QNAME_3.raw : 14676332 DNA.raw : [93924895, '2.2'] QUAL.raw : [165102599, '1.1']
        278399077    QUAL     ('QUAL', 'QNAME')           QUAL.raw : [137272420, '2.2'] DNA.key : 32923615 DNA : [69071159, '0.2'] QNAME_1.raw : 7617093 QNAME_2.raw : 15645298 QNAME_3.raw : 15869492
        278755833    DNA      ('DNA', 'QNAME')            DNA.raw : [69145778, '0.2'] QUAL.key : 32929386 QUAL : [137254691, '3.1'] QNAME_1.raw : 7826712 QNAME_2.raw : 15727921 QNAME_3.raw : 15871345
        278923935    DNA      ('DNA', 'QUAL', 'QNAME')    DNA.raw : [69145778, '0.2'] QUAL.raw : [170352179, '1.1'] QNAME_1.raw : 7826712 QNAME_2.raw : 15727921 QNAME_3.raw : 15871345
        279267421    QNAME    ('DNA', 'QNAME')            QNAME_1.raw : 11523 QNAME_2.raw : 396943 QNAME_3.raw : 14676332 DNA.raw : [93924895, '2.2'] QUAL.key : 33003037 QUAL : [137254691, '3.1']
        280817584    QUAL     ('DNA', 'QUAL')             QUAL.raw : [137272420, '2.2'] DNA.raw : [95416946, '2.2'] QNAME.key : 33043420 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332
        282419721    QNAME    ('QUAL', 'QNAME')           QNAME_1.raw : 11523 QNAME_2.raw : 396943 QNAME_3.raw : 14676332 DNA.key : 33161165 DNA : [69071159, '0.2'] QUAL.raw : [165102599, '1.1']
        285310498    QUAL     ('DNA', 'QNAME')            QUAL.key : 13506978 QUAL : [137254691, '3.1'] DNA.raw : [95416946, '2.2'] QNAME_1.raw : 7617093 QNAME_2.raw : 15645298 QNAME_3.raw : 15869492
        287395412    QUAL     ('QUAL',)                   QUAL.raw : [137272420, '2.2'] DNA.key : 32923615 DNA : [69071159, '0.2'] QNAME.key : 33043420 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332
        287574850    QNAME    ('QNAME',)                  QNAME_1.raw : 11523 QNAME_2.raw : 396943 QNAME_3.raw : 14676332 DNA.key : 33161165 DNA : [69071159, '0.2'] QUAL.key : 33003037 QUAL : [137254691, '3.1']
        287615721    DNA      ('DNA',)                    DNA.raw : [69145778, '0.2'] QUAL.key : 32929386 QUAL : [137254691, '3.1'] QNAME.key : 33201068 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332
        287741941    QNAME    ('DNA', 'QUAL')             QNAME.key : 13629649 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332 DNA.raw : [93924895, '2.2'] QUAL.raw : [165102599, '1.1']
        287783823    DNA      ('DNA', 'QUAL')             DNA.raw : [69145778, '0.2'] QUAL.raw : [170352179, '1.1'] QNAME.key : 33201068 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332
        291888326    QUAL     ('QNAME',)                  QUAL.key : 13506978 QUAL : [137254691, '3.1'] DNA.key : 32923615 DNA : [69071159, '0.2'] QNAME_1.raw : 7617093 QNAME_2.raw : 15645298 QNAME_3.raw : 15869492
        292143545    None     ('DNA', 'QUAL', 'QNAME')    QNAME_1.raw : 11523 QNAME_2.raw : 15675788 QNAME_3.raw : 15942048 DNA.raw : [93845134, '0.2'] QUAL.raw : [166669052, '1.1']
        292183269    DNA      ('QNAME',)                  DNA.key : 13502055 DNA : [69071159, '0.2'] QUAL.key : 32929386 QUAL : [137254691, '3.1'] QNAME_1.raw : 7826712 QNAME_2.raw : 15727921 QNAME_3.raw : 15871345
        292351371    DNA      ('QUAL', 'QNAME')           DNA.key : 13502055 DNA : [69071159, '0.2'] QUAL.raw : [170352179, '1.1'] QNAME_1.raw : 7826712 QNAME_2.raw : 15727921 QNAME_3.raw : 15871345
        292897070    QNAME    ('DNA',)                    QNAME.key : 13629649 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332 DNA.raw : [93924895, '2.2'] QUAL.key : 33003037 QUAL : [137254691, '3.1']
        294306833    QUAL     ('DNA',)                    QUAL.key : 13506978 QUAL : [137254691, '3.1'] DNA.raw : [95416946, '2.2'] QNAME.key : 33043420 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332
        295743894    None     ('DNA', 'QNAME')            QNAME_1.raw : 11523 QNAME_2.raw : 15675788 QNAME_3.raw : 15942048 DNA.raw : [93845134, '0.2'] QUAL.key : 33014710 QUAL : [137254691, '3.1']
        296049370    QNAME    ('QUAL',)                   QNAME.key : 13629649 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332 DNA.key : 33161165 DNA : [69071159, '0.2'] QUAL.raw : [165102599, '1.1']
        300334196    None     ('QUAL', 'QNAME')           QNAME_1.raw : 11523 QNAME_2.raw : 15675788 QNAME_3.raw : 15942048 DNA.key : 32964626 DNA : [69071159, '0.2'] QUAL.raw : [166669052, '1.1']
        300884661    QUAL     (None,)                     QUAL.key : 13506978 QUAL : [137254691, '3.1'] DNA.key : 32923615 DNA : [69071159, '0.2'] QNAME.key : 33043420 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332
        301043157    DNA      (None,)                     DNA.key : 13502055 DNA : [69071159, '0.2'] QUAL.key : 32929386 QUAL : [137254691, '3.1'] QNAME.key : 33201068 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332
        301204499    QNAME    (None,)                     QNAME.key : 13629649 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332 DNA.key : 33161165 DNA : [69071159, '0.2'] QUAL.key : 33003037 QUAL : [137254691, '3.1']
        301211259    DNA      ('QUAL',)                   DNA.key : 13502055 DNA : [69071159, '0.2'] QUAL.raw : [170352179, '1.1'] QNAME.key : 33201068 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332
        302886612    None     ('DNA', 'QUAL')             QNAME.key : 27287628 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332 DNA.raw : [93845134, '0.2'] QUAL.raw : [166669052, '1.1']
        303934545    None     ('QNAME',)                  QNAME_1.raw : 11523 QNAME_2.raw : 15675788 QNAME_3.raw : 15942048 DNA.key : 32964626 DNA : [69071159, '0.2'] QUAL.key : 33014710 QUAL : [137254691, '3.1']
        306486961    None     ('DNA',)                    QNAME.key : 27287628 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332 DNA.raw : [93845134, '0.2'] QUAL.key : 33014710 QUAL : [137254691, '3.1']
        311077263    None     ('QUAL',)                   QNAME.key : 27287628 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332 DNA.key : 32964626 DNA : [69071159, '0.2'] QUAL.raw : [166669052, '1.1']
        314677612    None     (None,)                     QNAME.key : 27287628 QNAME_1 : 11523 QNAME_2 : 396943 QNAME_3 : 14676332 DNA.key : 32964626 DNA : [69071159, '0.2'] QUAL.key : 33014710 QUAL : [137254691, '3.1']
Parameters found to be the best for this data type:
    --sort QUAL  --raw DNA QUAL QNAME  --pattern 2.2 2.2

Writing final config...
Archiving results and cleaning up temp directory...
All Done! :)
```

This will output sample.fastq.uQ, which is the data encoded with the best combination of parameters found, however as this was only a sample of the data we are really far more interested in what parameter combination was found to be best, output on the final line. Note, as you become familiar with sort/raw/pattern, you will know what you want (for example, definitely sort on QUAL, definitely DNA as a raw table, etc), then --test will just try the parameters you didn't specify, and output the best combination, and thus in many ways running on a sample is not really worth it. You might as well try it on the full dataset.
  
Regardless, one would then run:
`python uq.py --input ./rep3.fastq --sort QUAL --raw DNA QUAL QNAME --pattern 0.2 2.2`

Which outputs much the same as above (but more data for the metrics) and a file "rep3.fastq.uQ" To really get the smallest files however, you have to use a slower compressor, and try with/without the --pad option, and with/without the --notricks option, to cover all the bases. As --pad and --notricks change the actual encoding and not just patterning/sorting, and encoding is currently really rather slow compared to what it could be in C, these options are not checked by --test.

#Results

Adding the results from `python uq.py --input ./rep3.fastq --test --compressor lzma`, gives us the following for this rep3 dataset (ENCODE ENCFF000BUQ.fastq - rep3.fastq.gz is the file as downloaded from the ENCODE).
Standard LZMA applied with -e. For fqzcomp and DSRC, all the maximum lossless compression modes were used. DSRC decompressed to be byte-for-byte identical to input fastq. fqzcomp was a little lossy in that quality scores for N were lost. uQ kept the quality scores for N, but did changed the order of entries to be sorted on QUAL.

![alt text](http://i.imgur.com/CPKb1Li.png "File size of the same dataset in FASTQ and uQ encoding, and gzip/lzma/fqzcomp compression")
