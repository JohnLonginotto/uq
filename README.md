# uQ - binary FASTQ

```
  uQ is a tool to convert FASTQ files into a smaller, easily compressable binary format, and back again.
  To see results, skip to the bottom.
```

# Why use uQ?

1. Being a binary format, not only does a uQ file take up considerably less space than a FASTQ file before compression, but a uQ file will often compress down to half the size of a FASTQ file **after compression**. Win win - and you don't have to rely on some no-longer-supported, one-off compression program for long term storage. Use whatever compressor you're comfortable with, as uQ will work with your compressor to define a encoding scheme specific to your data and your compressor.

2. Being a *structured* binary format, data in a uQ file can be read and operated on much quicker than a FASTQ file. The data is really nothing more than a collection of Numpy arrays (in NPY format) all tar'd together into a single file (with a json config file). This means one could store a very large number of uQ files, all in-memory, and perform certain analyses very quickly.

# How to use uQ:

Good compression depends on the compressor, the amount of information stored in the file, and how that information is physically arranged on disk before compression. As such, there are plenty of non-linear variables at play, and the only way to know what is best is to try everything (or have a good guess based on some similar data). The variables uQ currently accepts include:

###--sort:
```
  Re-sorts the FASTQ file to reduce file size. Values can be "DNA", "QUAL", "QNAME" and "None", 
  case insensitive. This will change the final output (it will be sorted), however if sorting 
  by "DNA" it is likely downstream mapping will be much faster, although "QUAL" usually gives 
  the smallest files.
```

###--raw:
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

###--pattern: 
  ```
  While --raw changes the logical layout of the data, --pattern changes the physical layout,
  i.e. the way in which bytes are physically written to disk, such as writing/reading all 
  rows backwards, writing data in column order and not row order, etc. Because there are 
  so many combinations of transformations, these are simply reffered to by ids: "0.1", "0.2", 
  "1.1", "1.2", "2.1", "2.2", "3.1", and "3.2". No one is expected to remember what these 
  numbers mean, because...
  ```

###--test:
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

###--compressor:
```
  How does uQ know which arrangements of the above options compresses the best? It compresses the data 
  with a compressor of your choosing with the --compressor option. This is the path or alias, i.e. "gzip".
  Compressor is simply subprocessed, data is fed to the stdin and it's stdout is piped to "wc -c" to 
  count bytes. The code can be modified to work with compression programs that have to read/write files 
  and not from stdin by uncommenting a few lines. Note that if you do not specify a --compressor, uQ 
  will try to optimised for smallest uncompressed filesize, and the results will be very different to
  what a compressor would likely give.
```

###--temp:
```
  At present, data has to be written to a temporary directory after encoding to reduce memory load. 
  This is handled by the python tempfile module, which uses the system preference (usually a random 
  directory in /tmp), however you can override that here.
```

###--pad:
```
  Some compressors prefer it when "symbols" are byte-aligned, i.e. some factor of 8 (2/4/8/16/32). 
  --pad does this. This will definitely increase the size of the uncompressed data, but may or may 
  not decrease the size of the file after compression.
```

###--notricks:
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

###--peek:
```
  This flag essentially makes uQ abort after analysing the file but before encoding starts.
```

###--input:
```
  FASTQ files in when encoding, uQ files in when decoding. Unfortunately, this will have to be a 
  file and not a pipe, as the file has to be read multiple times to tune the encoder to the data
  and prevent issues that might result in an undecodable file 10 years from now.
```

###--output: 
```
  Path to output. If no output specified, default when making uQ files is to append ".uQ" to the end
  of the original file name, and when decompressing uQ files back to FASTQ, default is to print to
  the standard out/terminal.
```

###--decode:
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
(pypy) John Laptop:Desktop John$ python uq.py --input sample.fastq --test --compressor gzip
Temporary directory: /var/folders/3t/qctkm9xn6zv871cvwy6pw31m0000gn/T/tmpAo529j
Pass 1 of 4: First pass through file to find optimal/safe encoding strategy. (98.30%)
Finished analysing file! (0.00661060015361 minutes)
    - total reads: 100000
    - total bases: 3600000

QNAME Analysis:
    - all QNAMEs prefixed with: @SL-XBG_1_FC30CT9AAXX:8:
    - QNAME field separators: :: ( 2 symbols )

DNA Analysis:
    - DNA sequences contained the following characters: A C G N T (5 in total)
    - There distribution is:
      T |##############################                                                                      |  30.064%      1082294
      A |###########################                                                                         |  27.341%       984265
      G |####################                                                                                |  20.578%       740814
      C |####################                                                                                |  20.547%       739681
      N |#                                                                                                   |   1.471%        52946
    - this means we will store each letter of DNA in 3 bits.
    - all DNA sequences are 36 bases long
    - therefore, as 108 bits must be used to store the DNA for each sequence, we will use 14 bytes per unique DNA sequence.

QUAL Analysis:
  [Breakdown for "G"]
      $ |                                                                                                    |   0.595%         4408
      % |###                                                                                                 |   3.373%        24991
      & |#####                                                                                               |   5.177%        38354
      ( |##                                                                                                  |   2.744%        20330
      ) |######                                                                                              |   6.942%        51428
      * |                                                                                                    |   0.214%         1584
      + |#########                                                                                           |   9.773%        72403
      , |#                                                                                                   |   1.043%         7730
      - |                                                                                                    |   0.113%          838
      . |########                                                                                            |   8.143%        60324
      / |####                                                                                                |   4.448%        32954
      0 |                                                                                                    |   0.001%           10
      1 |###                                                                                                 |   3.939%        29177
      2 |                                                                                                    |   0.240%         1781
      3 |########                                                                                            |   8.454%        62626
      4 |#                                                                                                   |   1.388%        10284
      5 |#                                                                                                   |   1.688%        12505
      6 |####                                                                                                |   4.924%        36477
      7 |#######                                                                                             |   7.407%        54872
      8 |#                                                                                                   |   1.350%         9999
      9 |###                                                                                                 |   3.522%        26094
      : |######                                                                                              |   6.002%        44465
      ; |##                                                                                                  |   2.642%        19575
      < |###############                                                                                     |  15.681%       116169
      = |                                                                                                    |   0.194%         1435
      A |                                                                                                    |   0.000%            1

  [Breakdown for "A"]
      $ |                                                                                                    |   0.480%         4729
      % |###                                                                                                 |   3.199%        31491
      & |#####                                                                                               |   5.459%        53734
      ( |###                                                                                                 |   3.353%        33007
      ) |#######                                                                                             |   7.207%        70934
      * |                                                                                                    |   0.183%         1804
      + |##############                                                                                      |  14.229%       140047
      , |                                                                                                    |   0.572%         5626
      - |                                                                                                    |   0.096%          943
      . |#######                                                                                             |   7.300%        71847
      / |###########                                                                                         |  11.924%       117368
      0 |                                                                                                    |   0.001%            6
      1 |#                                                                                                   |   1.836%        18076
      2 |                                                                                                    |   0.080%          789
      3 |###############                                                                                     |  15.467%       152241
      4 |                                                                                                    |   0.336%         3304
      5 |#                                                                                                   |   1.881%        18516
      6 |###                                                                                                 |   3.310%        32575
      7 |##########                                                                                          |  10.586%       104197
      8 |                                                                                                    |   0.247%         2430
      9 |###                                                                                                 |   3.341%        32883
      : |###                                                                                                 |   3.708%        36496
      ; |#                                                                                                   |   1.477%        14542
      < |###                                                                                                 |   3.604%        35472
      = |                                                                                                    |   0.119%         1169
      A |                                                                                                    |   0.004%           39

  [Breakdown for "T"]
      $ |                                                                                                    |   0.499%         5403
      % |##                                                                                                  |   2.930%        31707
      & |####                                                                                                |   4.459%        48261
      ( |##                                                                                                  |   2.381%        25765
      ) |######                                                                                              |   6.185%        66940
      * |                                                                                                    |   0.205%         2219
      + |#########                                                                                           |   9.602%       103923
      , |#                                                                                                   |   1.058%        11447
      - |                                                                                                    |   0.245%         2647
      . |#######                                                                                             |   7.229%        78235
      / |#####                                                                                               |   5.350%        57900
      0 |                                                                                                    |   0.005%           58
      1 |###                                                                                                 |   3.487%        37742
      2 |                                                                                                    |   0.382%         4129
      3 |##########                                                                                          |  10.117%       109497
      4 |#                                                                                                   |   1.545%        16720
      5 |##                                                                                                  |   2.343%        25355
      6 |###                                                                                                 |   3.613%        39106
      7 |#########                                                                                           |   9.978%       107994
      8 |                                                                                                    |   0.571%         6183
      9 |###                                                                                                 |   3.536%        38271
      : |#######                                                                                             |   7.500%        81177
      ; |####                                                                                                |   4.589%        49671
      < |###########                                                                                         |  11.935%       129172
      = |                                                                                                    |   0.255%         2765
      A |                                                                                                    |   0.001%            7

  [Breakdown for "C"]
      $ |                                                                                                    |   0.540%         3997
      % |####                                                                                                |   4.071%        30112
      & |######                                                                                              |   6.663%        49282
      ( |###                                                                                                 |   3.902%        28861
      ) |########                                                                                            |   8.338%        61671
      * |                                                                                                    |   0.227%         1682
      + |##############                                                                                      |  14.084%       104175
      , |                                                                                                    |   0.731%         5406
      - |                                                                                                    |   0.143%         1059
      . |######                                                                                              |   6.885%        50927
      / |##########                                                                                          |  10.246%        75789
      0 |                                                                                                    |   0.001%            6
      1 |#                                                                                                   |   1.904%        14087
      2 |                                                                                                    |   0.152%         1128
      3 |############                                                                                        |  12.332%        91219
      4 |                                                                                                    |   0.477%         3527
      5 |#                                                                                                   |   1.656%        12248
      6 |##                                                                                                  |   2.335%        17275
      7 |##########                                                                                          |  10.423%        77095
      8 |                                                                                                    |   0.173%         1279
      9 |###                                                                                                 |   3.292%        24351
      : |####                                                                                                |   4.757%        35188
      ; |#                                                                                                   |   1.983%        14665
      < |####                                                                                                |   4.616%        34143
      = |                                                                                                    |   0.068%          500
      A |                                                                                                    |   0.001%            9

  [Breakdown for "N"]
      $ |###                                                                                                 |   3.651%         1933
      % |################################################################################################    |  96.349%        51013

  [Total distribution]
    - the following values were seen as quality scores: $ % & ( ) * + , - . / 0 1 2 3 4 5 6 7 8 9 : ; < = A (26 in total)
    - There distribution is:
      + |###########                                                                                         |  11.682%       420548
      3 |###########                                                                                         |  11.544%       415583
      7 |#########                                                                                           |   9.560%       344158
      < |########                                                                                            |   8.749%       314956
      / |#######                                                                                             |   7.889%       284011
      . |#######                                                                                             |   7.259%       261333
      ) |######                                                                                              |   6.971%       250973
      : |#####                                                                                               |   5.481%       197326
      & |#####                                                                                               |   5.268%       189631
      % |####                                                                                                |   4.703%       169314
      6 |###                                                                                                 |   3.484%       125433
      9 |###                                                                                                 |   3.378%       121599
      ( |##                                                                                                  |   2.999%       107963
      1 |##                                                                                                  |   2.752%        99082
      ; |##                                                                                                  |   2.735%        98453
      5 |#                                                                                                   |   1.906%        68624
      4 |                                                                                                    |   0.940%        33835
      , |                                                                                                    |   0.839%        30209
      $ |                                                                                                    |   0.569%        20470
      8 |                                                                                                    |   0.553%        19891
      2 |                                                                                                    |   0.217%         7827
      * |                                                                                                    |   0.202%         7289
      = |                                                                                                    |   0.163%         5869
      - |                                                                                                    |   0.152%         5487
      0 |                                                                                                    |   0.002%           80
      A |                                                                                                    |   0.002%           56
    - this means we will store each quality symbol in 5 bits.
    - therefore, as 180 bits must be used to store the qualities for each sequence, we will use 23 bytes per unique quality sequence.

Pass 2 of 4: Now QNAME delimiters are known, QNAMEs are being analysed. (98.30%)
Optimal encoding method for delimited data in QNAME determined! (0.00494136810303 minutes)
    - Column 3 is type integers stored as uint8
    - Column 3 is type integers stored as uint16
    - Column 3 is type integers stored as uint16

Pass 3 of 4: Encoding DNA and Quality scores... (98.30%)
Pass 4 of 4: Encoding QNAMEs... (98.30%)

Starting tests...
Time:                       Sort:      Raw Tables:                      Patterns:
(0.0623526493708 minutes)   DNA        ('DNA', 'QUAL', 'QNAME')         All
(0.0601833820343 minutes)   QUAL       ('DNA', 'QUAL', 'QNAME')         All
(0.072601834933 minutes)    QNAME      ('DNA', 'QUAL', 'QNAME')         All
(0.0634381532669 minutes)   None       ('DNA', 'QUAL', 'QNAME')         All
(0.0788209994634 minutes)   DNA        ('DNA', 'QUAL')                  All
(0.0796946843465 minutes)   QUAL       ('DNA', 'QUAL')                  All
(0.082397600015 minutes)    QNAME      ('DNA', 'QUAL')                  All
(0.0825541973114 minutes)   None       ('DNA', 'QUAL')                  All
(0.0675637483597 minutes)   DNA        ('QUAL', 'QNAME')                All
(0.0637503345807 minutes)   QUAL       ('QUAL', 'QNAME')                All
(0.0769990324974 minutes)   QNAME      ('QUAL', 'QNAME')                All
(0.0671203494072 minutes)   None       ('QUAL', 'QNAME')                All
(0.0639218171438 minutes)   DNA        ('DNA', 'QNAME')                 All
(0.0672367493312 minutes)   QUAL       ('DNA', 'QNAME')                 All
(0.0754310687383 minutes)   QNAME      ('DNA', 'QNAME')                 All
(0.068313852946 minutes)    None       ('DNA', 'QNAME')                 All
(0.0819273670514 minutes)   DNA        ('DNA',)                         All
(0.085627647241 minutes)    QUAL       ('DNA',)                         All
(0.0848837812742 minutes)   QNAME      ('DNA',)                         All
(0.0848152001699 minutes)   None       ('DNA',)                         All
(0.0869405666987 minutes)   DNA        ('QUAL',)                        All
(0.0826579332352 minutes)   QUAL       ('QUAL',)                        All
(0.0877449154854 minutes)   QNAME      ('QUAL',)                        All
(0.0856060345968 minutes)   None       ('QUAL',)                        All
(0.0703301469485 minutes)   DNA        ('QNAME',)                       All
(0.0703383008639 minutes)   QUAL       ('QNAME',)                       All
(0.0802314162254 minutes)   QNAME      ('QNAME',)                       All
(0.0705535173416 minutes)   None       ('QNAME',)                       All
(0.0906253139178 minutes)   DNA        (None,)                          All
(0.0890094319979 minutes)   QUAL       (None,)                          All
(0.0900878508886 minutes)   QNAME      (None,)                          All
(0.0892335017522 minutes)   None       (None,)                          All

All done!
Size (compressed)    Sort:    Raw Tables:                 Raw stats:
          2745473    QUAL     ('DNA', 'QUAL', 'QNAME')    QUAL.raw : [1487597, '2.2'] DNA.raw : [935959, '0.2'] QNAME_1.raw : 510 QNAME_2.raw : 159528 QNAME_3.raw : 161879
          2768732    DNA      ('DNA', 'QUAL', 'QNAME')    DNA.raw : [756426, '1.1'] QUAL.raw : [1690312, '2.2'] QNAME_1.raw : 508 QNAME_2.raw : 159808 QNAME_3.raw : 161678
          2802126    QNAME    ('DNA', 'QUAL', 'QNAME')    QNAME_1.raw : 214 QNAME_2.raw : 4192 QNAME_3.raw : 149519 DNA.raw : [946462, '0.2'] QUAL.raw : [1701739, '1.1']
          2837749    DNA      ('DNA', 'QNAME')            DNA.raw : [756426, '1.1'] QUAL.key : 273826 QUAL : [1485503, '2.2'] QNAME_1.raw : 508 QNAME_2.raw : 159808 QNAME_3.raw : 161678
          2838636    QUAL     ('QUAL', 'QNAME')           QUAL.raw : [1487597, '2.2'] DNA.key : 273861 DNA : [755261, '0.2'] QNAME_1.raw : 510 QNAME_2.raw : 159528 QNAME_3.raw : 161879
          2854580    QUAL     ('DNA', 'QUAL')             QUAL.raw : [1487597, '2.2'] DNA.raw : [935959, '0.2'] QNAME.key : 277099 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519
          2862409    QNAME    ('DNA', 'QNAME')            QNAME_1.raw : 214 QNAME_2.raw : 4192 QNAME_3.raw : 149519 DNA.raw : [946462, '0.2'] QUAL.key : 276519 QUAL : [1485503, '2.2']
          2877502    DNA      ('DNA', 'QUAL')             DNA.raw : [756426, '1.1'] QUAL.raw : [1690312, '2.2'] QNAME.key : 276839 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519
          2880615    QUAL     ('DNA', 'QNAME')            QUAL.key : 137236 QUAL : [1485503, '2.2'] DNA.raw : [935959, '0.2'] QNAME_1.raw : 510 QNAME_2.raw : 159528 QNAME_3.raw : 161879
          2887448    QNAME    ('QUAL', 'QNAME')           QNAME_1.raw : 214 QNAME_2.raw : 4192 QNAME_3.raw : 149519 DNA.key : 276523 DNA : [755261, '0.2'] QUAL.raw : [1701739, '1.1']
          2904682    DNA      ('QUAL', 'QNAME')           DNA.key : 137115 DNA : [755261, '0.2'] QUAL.raw : [1690312, '2.2'] QNAME_1.raw : 508 QNAME_2.raw : 159808 QNAME_3.raw : 161678
          2940580    QNAME    ('DNA', 'QUAL')             QNAME.key : 138454 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519 DNA.raw : [946462, '0.2'] QUAL.raw : [1701739, '1.1']
          2946519    DNA      ('DNA',)                    DNA.raw : [756426, '1.1'] QUAL.key : 273826 QUAL : [1485503, '2.2'] QNAME.key : 276839 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519
          2947731    QNAME    ('QNAME',)                  QNAME_1.raw : 214 QNAME_2.raw : 4192 QNAME_3.raw : 149519 DNA.key : 276523 DNA : [755261, '0.2'] QUAL.key : 276519 QUAL : [1485503, '2.2']
          2947743    QUAL     ('QUAL',)                   QUAL.raw : [1487597, '2.2'] DNA.key : 273861 DNA : [755261, '0.2'] QNAME.key : 277099 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519
          2965989    None     ('DNA', 'QUAL', 'QNAME')    QNAME_1.raw : 214 QNAME_2.raw : 159554 QNAME_3.raw : 161787 DNA.raw : [942589, '3.1'] QUAL.raw : [1701845, '1.1']
          2973699    DNA      ('QNAME',)                  DNA.key : 137115 DNA : [755261, '0.2'] QUAL.key : 273826 QUAL : [1485503, '2.2'] QNAME_1.raw : 508 QNAME_2.raw : 159808 QNAME_3.raw : 161678
          2973778    QUAL     ('QNAME',)                  QUAL.key : 137236 QUAL : [1485503, '2.2'] DNA.key : 273861 DNA : [755261, '0.2'] QNAME_1.raw : 510 QNAME_2.raw : 159528 QNAME_3.raw : 161879
          2989722    QUAL     ('DNA',)                    QUAL.key : 137236 QUAL : [1485503, '2.2'] DNA.raw : [935959, '0.2'] QNAME.key : 277099 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519
          3000863    QNAME    ('DNA',)                    QNAME.key : 138454 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519 DNA.raw : [946462, '0.2'] QUAL.key : 276519 QUAL : [1485503, '2.2']
          3013452    DNA      ('QUAL',)                   DNA.key : 137115 DNA : [755261, '0.2'] QUAL.raw : [1690312, '2.2'] QNAME.key : 276839 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519
          3024924    None     ('DNA', 'QNAME')            QNAME_1.raw : 214 QNAME_2.raw : 159554 QNAME_3.raw : 161787 DNA.raw : [942589, '3.1'] QUAL.key : 275277 QUAL : [1485503, '2.2']
          3025902    QNAME    ('QUAL',)                   QNAME.key : 138454 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519 DNA.key : 276523 DNA : [755261, '0.2'] QUAL.raw : [1701739, '1.1']
          3052158    None     ('QUAL', 'QNAME')           QNAME_1.raw : 214 QNAME_2.raw : 159554 QNAME_3.raw : 161787 DNA.key : 273497 DNA : [755261, '0.2'] QUAL.raw : [1701845, '1.1']
          3075112    None     ('DNA', 'QUAL')             QNAME.key : 276753 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519 DNA.raw : [942589, '3.1'] QUAL.raw : [1701845, '1.1']
          3082469    DNA      (None,)                     DNA.key : 137115 DNA : [755261, '0.2'] QUAL.key : 273826 QUAL : [1485503, '2.2'] QNAME.key : 276839 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519
          3082885    QUAL     (None,)                     QUAL.key : 137236 QUAL : [1485503, '2.2'] DNA.key : 273861 DNA : [755261, '0.2'] QNAME.key : 277099 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519
          3086185    QNAME    (None,)                     QNAME.key : 138454 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519 DNA.key : 276523 DNA : [755261, '0.2'] QUAL.key : 276519 QUAL : [1485503, '2.2']
          3111093    None     ('QNAME',)                  QNAME_1.raw : 214 QNAME_2.raw : 159554 QNAME_3.raw : 161787 DNA.key : 273497 DNA : [755261, '0.2'] QUAL.key : 275277 QUAL : [1485503, '2.2']
          3134047    None     ('DNA',)                    QNAME.key : 276753 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519 DNA.raw : [942589, '3.1'] QUAL.key : 275277 QUAL : [1485503, '2.2']
          3161281    None     ('QUAL',)                   QNAME.key : 276753 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519 DNA.key : 273497 DNA : [755261, '0.2'] QUAL.raw : [1701845, '1.1']
          3220216    None     (None,)                     QNAME.key : 276753 QNAME_1 : 214 QNAME_2 : 4192 QNAME_3 : 149519 DNA.key : 273497 DNA : [755261, '0.2'] QUAL.key : 275277 QUAL : [1485503, '2.2']

  Parameters found to be the best for this data type:
      --sort QUAL --raw DNA QUAL QNAME --pattern 0.2 2.2

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
