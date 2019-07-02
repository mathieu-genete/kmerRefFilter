# README for kmerRefFilter
The huge amount of data generated by high throughput sequencing technologies (NGS) is a limiting factor for much of computer analysis. Performing these analyzes requires access to infrastructures with sufficient storage space and appropriate computing power. But for some analyzes, only small parts of the initial data are of interest. Reducing the data set, in this case, makes it possible to optimize the storage space, the processor calculation time and the analysis accuracy.
We developed a dedicated algorithm (kmerRefFilter) to reduce, either locally or directly on a download stream, any NGS dataset. This tool produce an exhaustive library of k-mers present in the set of reference sequences, from which duplicates and low-complexity k-mers are removed. This library is next use to filter and reduce the NGS dataset. 

## Dependencies
kmerRefFilter is written with Python 2.7.5 and requires biopython (https://github.com/biopython/biopython) package installed before running.

## Running kmerRefFilter
----------------------------
```
usage: kmerRefFilter [-h] [-v] [-prog] [-a] [-d] [-y] [-knf] [-o OUTPUTDIR]
                     [-k KMERSIZE] [-m MINMATCHKEEPSEQ] [-z MINSHENTROPY]
                     [-q MAXRATIOAMBIGOUS] [-e [EXTREMITY5P3P]]
                     [-r [REFERENCESFASTA [REFERENCESFASTA ...]]]
                     [-i KMERINFILE] [-p KMEROUTPUT]
                     [-x [EXCLUDEFASTA [EXCLUDEFASTA ...]]] [-maxRD MAXREADS]
                     [-f [FASTQFILE [FASTQFILE ...]]] [-l FASTQLIST]
                     [-1 [FWDFASTQ [FWDFASTQ ...]]]
                     [-2 [REVFASTQ [REVFASTQ ...]]]
                     [-u [FASTQURL [FASTQURL ...]]] [-u1 FWDFASTQURL]
                     [-u2 REVFASTQURL] [-ugzip] [-s STREAMINPUT]
                     [-c FWDREVCHOICE] [-kc KMERFWDREV]

Tool for raw reads kmer-based filtering.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -prog, --printprogress
                        print progession if set
  -a, --append          if fastq filtered exists, add new reads at the end of
                        file
  -d, --ambigousDNA     generate all possible kmers from an ambigous kmer
  -y, --yamlout         generate yaml output on stdout
  -knf, --keepNotFiltered
                        keep not filtered reads in a file _KEEPED.fastq
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        output directory for filtered fastq files
  -k KMERSIZE, --kmersize KMERSIZE
                        size for kmers - default=20
  -m MINMATCHKEEPSEQ, --minMatchKeepSeq MINMATCHKEEPSEQ
                        minimum number of ref reads matched to keep fastq read
                        -- default=1
  -z MINSHENTROPY, --minShEntropy MINSHENTROPY
                        minimum Shannon Entropy for kmers (0 to 2.0)--
                        defaut=0.8
  -q MAXRATIOAMBIGOUS, --maxratioAmbigous MAXRATIOAMBIGOUS
                        maximum ambigous bases accepted in kmers in % --
                        defaut=0.2
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
  -f [FASTQFILE [FASTQFILE ...]], --fastqfile [FASTQFILE [FASTQFILE ...]]
                        Fastq to filter (fastq, gz, bz2)
  -l FASTQLIST, --fastqlist FASTQLIST
                        text file with all fastq (fastq, gz, bz2) path -- 1
                        per line
  -1 [FWDFASTQ [FWDFASTQ ...]], --fwdfastq [FWDFASTQ [FWDFASTQ ...]]
                        forward Fastq to filter (fastq, gz, bz2) respect file
                        order with reverse files
  -2 [REVFASTQ [REVFASTQ ...]], --revfastq [REVFASTQ [REVFASTQ ...]]
                        reverse Fastq to filter (fastq, gz, bz2) respect file
                        order with forward files
  -u [FASTQURL [FASTQURL ...]], --fastqurl [FASTQURL [FASTQURL ...]]
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
```

### Basic command line usage:

`python kmerRefFilter.py -r reference_sequences.fasta -1 fastq1_R1 fastq2_R1 -2 fastq1_R2 fastq2_R2 -o output_directory`

### Command line usage for URL files path

- if the files are compressed in gzip format, use -ugzip option.
      
example for gz compressed fastq files:
`python kmerRefFilter.py -r reference_sequences.fasta -u1 fastq_url_R1 -u2 fastq_url_R2 -ugzip -o output_directory`

## Contact Information
Mathieu Genete

Email: mathieu.genete@univ-lille.fr

## Licence Agreement
This software is covered by GNU General Public License, version 3 (GPL-3.0).
