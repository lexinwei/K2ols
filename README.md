# K2ols
 Tools for Kraken2 DB build and taxonomy classification.

## kraken2M.py
```
$ python kraken2M.py --help
usage: kraken2M [-h] -i INPUT -s SUFFIX -d DB -k KRAKEN -kt KRAKEN_TOOLS [-o OUTPUT] [-c CONFIDENCE] [-t THREADS]
                [--gzip-compressed]

Species classification for multiple (paired/single-end) fastq files one time with Kraken2, in order not to load the
super index repeatly. Most of the following arguments are originated from Kraken2, not supporting the modification
of other Kraken2 arguments, just using their default value setted by Kraken2. Please clone KrakenTools by
jenniferlu717 from https://github.com/lexinwei/KrakenTools.git before running.

optional arguments:
  -h, --help            show this help message and exit

REQUIRED PARAMETERS:
  -i INPUT, --input INPUT
                        A directory of multiple fastq files, only the files with names ending in specific suffix
                        (given by -s arguments) will be loaded. (default: None)
  -s SUFFIX, --suffix SUFFIX
                        Specify the suffix of reads files, e.g. 'R1.fastq,R2.fastq' for paired-end files, '.fq' for
                        single-end files. (default: None)
  -d DB, --db DB        Name for Kraken2 DB. (default: None)
  -k KRAKEN, --kraken KRAKEN
                        Path of Kraken2 program. (default: None)
  -kt KRAKEN_TOOLS, --kraken-tools KRAKEN_TOOLS
                        Path of KrakenTools by jenniferlu717. (default: None)

OPTIONAL PARAMETERS:
  -o OUTPUT, --output OUTPUT
                        A directory for saving output files. (default: ./kraken2_output)
  -c CONFIDENCE, --confidence CONFIDENCE
                        Confidence score threshold, must be in [0, 1]. (default: 0)
  -t THREADS, --threads THREADS
                        Number of threads to running Kraken2. (default: 1)
  --gzip-compressed     Input files are compressed with gzip. (default: False)
```
