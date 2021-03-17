# README for kmerRefFilter
The huge amount of data generated by high throughput sequencing technologies (NGS) is a limiting factor for much of computer analysis. Performing these analyzes requires access to infrastructures with sufficient storage space and appropriate computing power. But for some analyzes, only small parts of the initial data are of interest. Reducing the data set, in this case, makes it possible to optimize the storage space, the processor calculation time and the analysis accuracy.
We developed a dedicated algorithm (kmerRefFilter) to reduce, either locally or directly on a download stream, any NGS dataset. This tool produce an exhaustive library of k-mers present in the set of reference sequences, from which duplicates and low-complexity k-mers are removed. This library is next use to filter and reduce the NGS dataset. test

## Dependencies
kmerRefFilter is written with Python 2.7.5 and requires biopython (https://github.com/biopython/biopython) package installed before running.

## Running kmerRefFilter
```
usage: kmerRefFilter [-h] [-v] [-prog] [-nc] [-a] [-y] [-knf] [-mem MAXMEMORY]
                     [-P MULTIPROC] [-o OUTPUTDIR] [-nolog] [-nov]
                     [-k KMERSIZE] [-ks [KMERSTATS]] [-m MINMATCHKEEPSEQ]
                     [-z MINSHENTROPY] [-q MAXRATIOAMBIGOUS]
                     [-j MAXREPEATSFREQ] [-b] [-e [EXTREMITY5P3P]]
                     [-r [REFERENCESFASTA [REFERENCESFASTA ...]]]
                     [-i KMERINFILE] [-p KMEROUTPUT]
                     [-x [EXCLUDEFASTA [EXCLUDEFASTA ...]]] [-maxRD MAXREADS]
                     [-f FASTQFILE] [-ct] [-1 FWDFASTQ] [-2 REVFASTQ]
                     [-u FASTQURL] [-u1 FWDFASTQURL] [-u2 REVFASTQURL]
                     [-ugzip] [-s STREAMINPUT] [-c FWDREVCHOICE]
                     [-kc KMERFWDREV] [-morf MINPORF] [-Morf MAXPORF]
                     [-mkc MINKMERCOUNT] [-Mkc MAXKMERCOUNT] [-sorf STRANDORF]
                     [-pfreq MAXPROTFREQ] [-sbc SINGBOTHCODING]

Tool for raw reads kmer-based filtering.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -prog, --printprogress
                        print progession if set
  -nc, --nokmerclean    Disable kmers cleaning step
  -a, --append          if fastq filtered exists, add new reads at the end of
                        file
  -y, --yamlout         generate yaml output on stdout
  -knf, --keepNotFiltered
                        keep not filtered reads in a file _KEEPED.fastq
  -mem MAXMEMORY, --maxmemory MAXMEMORY
                        maximum memory used in Go (default = 5 Go) - set to 0
                        for unlimited memory usage
  -P MULTIPROC, --multiproc MULTIPROC
                        use multiprocessing for kmer dictionary generation --
                        number of CPU
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        output directory for filtered fastq files
  -nolog, --nologfile   Disable log file output
  -nov, --noverbose     Disable output log on stderr
  -k KMERSIZE, --kmersize KMERSIZE
                        size for kmers - default=20
  -ks [KMERSTATS], --kmerstats [KMERSTATS]
                        generate kmer statistics output table USED or ALL
                        (USED by default)
  -m MINMATCHKEEPSEQ, --minMatchKeepSeq MINMATCHKEEPSEQ
                        minimum number of ref reads matched to keep fastq read
                        -- default=1
  -z MINSHENTROPY, --minShEntropy MINSHENTROPY
                        minimum Shannon Entropy for kmers (0 to 2.0)--
                        defaut=0.8
  -q MAXRATIOAMBIGOUS, --maxratioAmbigous MAXRATIOAMBIGOUS
                        maximum frequency of ambigous bases accepted in kmers
                        -- defaut=0.2
  -j MAXREPEATSFREQ, --maxRepeatsFreq MAXREPEATSFREQ
                        maximum frequency of 2-4 mers repetitions accepted in
                        kmers -- defaut=0.8
  -b, --mutedkmers      generates all single bases replacement possibilities
                        for each dictionary kmers (huge memory consumption)
  -e [EXTREMITY5P3P], --extremity5p3p [EXTREMITY5P3P]
                        compare 5' and 3' read extremity only -- nbr of bases
                        to test from extremities default=1
  -r [REFERENCESFASTA [REFERENCESFASTA ...]], --referencesfasta [REFERENCESFASTA [REFERENCESFASTA ...]]
                        fasta files with references sequences
  -i KMERINFILE, --kmerinfile KMERINFILE
                        load kmer dictionary from a saved dictionary file
  -p KMEROUTPUT, --kmeroutput KMEROUTPUT
                        only save kmer dictionary in a file and exit without
                        filtering
  -x [EXCLUDEFASTA [EXCLUDEFASTA ...]], --excludefasta [EXCLUDEFASTA [EXCLUDEFASTA ...]]
                        fasta files with sequences to exclude
  -maxRD MAXREADS, --maxreads MAXREADS
                        maximum number of read to analyze
  -f FASTQFILE, --fastqfile FASTQFILE
                        Fastq to filter (fastq, gz, bz2)
  -ct, --concatenatedfastq
                        input fastq files are concatenated (FWD+REV in one
                        read). Split half reads and quality values. Should be
                        use with -f or -l or -u or -s options only
  -1 FWDFASTQ, --fwdfastq FWDFASTQ
                        forward Fastq to filter (fastq, gz, bz2) respect file
                        order with reverse files
  -2 REVFASTQ, --revfastq REVFASTQ
                        reverse Fastq to filter (fastq, gz, bz2) respect file
                        order with forward files
  -u FASTQURL, --fastqurl FASTQURL
                        url to Fastq to filter (fastq, gz)
  -u1 FWDFASTQURL, --fwdfastqurl FWDFASTQURL
                        url to forward Fastq to filter (fastq, gz)
  -u2 REVFASTQURL, --revfastqurl REVFASTQURL
                        url to reverse Fastq to filter (fastq, gz)
  -ugzip, --urlgzip     indicate that url's given in -u or in -u1 and -u2 are
                        gz files
  -s STREAMINPUT, --streamInput STREAMINPUT
                        use stdind as input - specify output filename
  -c FWDREVCHOICE, --fwdRevChoice FWDREVCHOICE
                        for paired filtering: filtering on forward reads only
                        (FWD), reverse reads only (REV) or both (BOTH - by
                        default)
  -kc KMERFWDREV, --kmerFwdRev KMERFWDREV
                        kmders dictionary generation: use forward kmers only
                        (FWD), reverse kmers only (REV) or both (BOTH - by
                        default)
  -morf MINPORF, --minporf MINPORF
                        minimum percent of ORF in read
  -Morf MAXPORF, --maxporf MAXPORF
                        maximum percent of ORF in read
  -mkc MINKMERCOUNT, --minkmercount MINKMERCOUNT
                        minimum kmer count for recruited reads
  -Mkc MAXKMERCOUNT, --maxkmercount MAXKMERCOUNT
                        maximum kmer count for recruited reads
  -sorf STRANDORF, --strandorf STRANDORF
                        force strand for ORF filtering (both=0, +=1, -=2) --
                        default=0
  -pfreq MAXPROTFREQ, --maxprotfreq MAXPROTFREQ
                        maximum frequency of AA repetitions in proteics
                        sequences from reads (default = 0.6)
  -sbc SINGBOTHCODING, --singbothcoding SINGBOTHCODING
                        apply ORF threshold for single (SINGLE); both (BOTH)
                        forward and reverse reads to keep the pair or single
                        exclude both (SINGEXBOTH) (SINGLE - by default)
```

### Basic command line usage

using fasta file as reference:

`python kmerRefFilter.py -r reference_sequences.fasta -1 fastq1_R1 fastq2_R1 -2 fastq1_R2 fastq2_R2 -o output_directory`

or using a saved dictionary (previously saved with -p option):

`python kmerRefFilter.py -i saved_dictionary -1 fastq1_R1 fastq2_R1 -2 fastq1_R2 fastq2_R2 -o output_directory`

filtered fastq files are created in output_directory.

### Basic Command line usage on download stream

Fastq files are directly filtered on dowload stream. If the files are compressed in gzip format, use -ugzip option.

example for gz compressed fastq files:

`python kmerRefFilter.py -r reference_sequences.fasta -u1 fastq_url_R1 -u2 fastq_url_R2 -ugzip -o output_directory`

### exclude option (-x)

Need a fasta file as input. Generates a set of kmers, using the fasta file, and next exclude them from the kmer reference dictionary.

### Contact Information
Mathieu Genete

Email: mathieu.genete@univ-lille.fr

### Licence Agreement
This software is covered by GNU General Public License, version 3 (GPL-3.0).
