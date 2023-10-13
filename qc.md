# Quality Control and Trimming

## Let's get started!

In this practical you will learn to import, view and check the quality of raw high thoughput sequencing sequencing data.

The first dataset you will be working with is from an Illumina MiSeq dataset.

The sequenced organism is an enterohaemorrhagic E. coli (EHEC) of the serotype O157, a potentially fatal gastrointestinal pathogen.

The sequenced bacterium was part of an outbreak investigation in the St. Louis area, USA in 2011.

The sequencing was done as paired-end 2x150bp.

Click the **Start** button to move to the next step.

## File Formats

### The fasta format

The fasta format was invented in 1988 and designed to represent nucleotide or peptide sequences. It originates from the [FASTA](https://en.wikipedia.org/wiki/FASTA) software package, but is now a standard in the world of bioinformatics.

The first line in a FASTA file starts with a ">" (greater-than) symbol followed by the description or identifier of the sequence. Following the initial line (used for a unique description of the sequence) is the actual sequence itself in standard one-letter code.

A few sample sequences:

```
>KX580312.1 Homo sapiens truncated breast cancer 1 (BRCA1) gene, exon 15 and partial cds
GTCATCCCCTTCTAAATGCCCATCATTAGATGATAGGTGGTACATGCACAGTTGCTCTGGGAGTCTTCAG
AATAGAAACTACCCATCTCAAGAGGAGCTCATTAAGGTTGTTGATGTGGAGGAGTAACAGCTGGAAGAGT
CTGGGCCACACGATTTGACGGAAACATCTTACTTGCCAAGGCAAGATCTAG
```

```
>KRN06561.1 heat shock [Lactobacillus sucicola DSM 21376 = JCM 15457]
MSLVMANELTNRFNNWMKQDDFFGNLGRSFFDLDNSVNRALKTDVKETDKAYEVRIDVPGIDKKDITVDY
HDGVLSVNAKRDSFNDESDSEGNVIASERSYGRFARQYSLPNVDESGIKAKCEDGVLKLTLPKLAEEKIN
GNHIEIE
```

A fasta file can contain multiple sequence. Each sequence will be separated by their "header" line, starting by ">".

Example:

```
>KRN06561.1 heat shock [Lactobacillus sucicola DSM 21376 = JCM 15457]
MSLVMANELTNRFNNWMKQDDFFGNLGRSFFDLDNSVNRALKTDVKETDKAYEVRIDVPGIDKKDITVDY
HDGVLSVNAKRDSFNDESDSEGNVIASERSYGRFARQYSLPNVDESGIKAKCEDGVLKLTLPKLAEEKIN
GNHIEIE
>3HHU_A Chain A, Human Heat-Shock Protein 90 (Hsp90)
MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNKEIFLRELISNSSDALDKIRYESLTDPSKL
DSGKELHINLIPNKQDRTLTIVDTGIGMTKADLINNLGTIAKSGTKAFMEALQAGADISMIGQFGVGFYS
AYLVAEKVTVITKHNDDEQYAWESSAGGSFTVRTDTGEPMGRGTKVILHLKEDQTEYLEERRIKEIVKKH
SQFIGYPITLFVEK
```

### The fastq format

The fastq format is also a text based format to represent nucleotide sequences, but also contains the corresponding quality of each nucleotide. It is the standard for storing the output of high-throughput sequencing instruments such as the Illumina machines.

A fastq file uses four lines per sequence:

* Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
* Line 2 is the raw sequence letters.
* Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
* Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

An example sequence in fastq format:

```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

#### Quality

The quality, also called phred score, is the probability that the corresponding basecall is incorrect.

Phred scores use a logarithmic scale, and are represented by ASCII characters, mapping to a quality usually going from 0 to 40.

| Phred Quality Score | Probability of incorrect base call | Base call accuracy
| --- | --- | ---
| 10 | 1 in 10 | 90%
| 20 | 1 in 100 | 99%
| 30 | 1 in 1000 | 99.9%
| 40 | 1 in 10,000 | 99.99%
| 50 | 1 in 100,000 | 99.999%
| 60 | 1 in 1,000,000 | 99.9999%

## Downloading the data

The raw data were deposited at the European Nucleotide Archive, under the accession number SRR957824.
You could go to the ENA [website](http://www.ebi.ac.uk/ena) and search for the run with the accession SRR957824.

However these files contain about 3 million reads and are therefore quite big.
We are only going to use a subset of the original dataset for this tutorial.

First create a `data/` directory in your home folder

```bash
mkdir ~/data
```

now let's download the subset

```bash
cd ~/data
curl -O -J -L https://osf.io/shqpv/download
curl -O -J -L https://osf.io/9m3ch/download
```

Let’s make sure we downloaded all of our data using md5sum.

```bash
md5sum SRR957824_500K_R1.fastq.gz SRR957824_500K_R2.fastq.gz
```

you should see this

```
1e8cf249e3217a5a0bcc0d8a654585fb  SRR957824_500K_R1.fastq.gz
70c726a31f05f856fe942d727613adb7  SRR957824_500K_R2.fastq.gz
```

and now look at the file names and their size

```bash
ls -l
```

```
total 97M
-rw-r--r-- 1 hadrien 48M Nov 19 18:44 SRR957824_500K_R1.fastq.gz
-rw-r--r-- 1 hadrien 50M Nov 19 18:53 SRR957824_500K_R2.fastq.gz
```

There are 500 000 paired-end reads taken randomly from the original data

One last thing before we get to the quality control: those files are writeable.
By default, UNIX makes things writeable by the file owner.
This poses an issue with creating typos or errors in raw data.
We fix that before going further

```bash
chmod u-w *
```

## Working Directory

First we make a work directory: a directory where we can play around with a copy of the data without messing with the original

```bash
mkdir ~/work
cd ~/work
```

Now we make a link of the data in our working directory

```bash
ln -s ~/data/* .
```

The files that we've downloaded are FASTQ files. Take a look at one of them with

```bash
zless SRR957824_500K_R1.fastq.gz
```

:bulb:
Use the spacebar to scroll down, and type ‘q’ to exit ‘less’

You can read more on the FASTQ format in the [File Formats](file_formats.md) lesson.

:question:
Where does the filename come from?

:question:
Why are there 1 and 2 in the file names?

## FastQC

To check the quality of the sequence data we will use a tool called FastQC.

FastQC has a graphical interface and can be downloaded and run on a Windows or Linux computer without installation.
It is available [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

However, FastQC is also available as a command line utility on the training server you are using.
To run FastQC on our two files

```bash
fastqc SRR957824_500K_R1.fastq.gz SRR957824_500K_R2.fastq.gz
```

and look what FastQC has produced

```
ls *fastqc*
```

For each file, FastQC has produced both a .zip archive containing all the plots, and a html report.

Download and open the html files with your favourite web browser.

Alternatively you can look a these copies of them:

- [SRR957824_500K_R1_fastqc.html](data/fastqc/SRR957824_500K_R1_fastqc.html)
- [SRR957824_500K_R2_fastqc.html](data/fastqc/SRR957824_500K_R2_fastqc.html)

:question:
What should you pay attention to in the FastQC report?

:question:
Which file is of better quality?

Pay special attention to the per base sequence quality and sequence length distribution.
Explanations for the various quality modules can be found [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).
Also, have a look at examples of a [good](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) and a [bad](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) illumina read set for comparison.

You will note that the reads in your uploaded dataset have fairly poor quality (<20) towards the end. There are also outlier reads that have very poor quality for most of the second half of the reads.

## Scythe

Now we'll do some trimming!

Scythe uses a Naive Bayesian approach to classify contaminant substrings in sequence reads.
It considers quality information, which can make it robust in picking out 3'-end adapters, which often include poor quality bases.

The first thing we need is the adapters to trim off

```bash
curl -O -J -L https://osf.io/v24pt/download
```

Now we run scythe on both our read files

```bash
scythe -a adapters.fasta SRR957824_500K_R1.fastq.gz -o SRR957824_adapt_R1.fastq 
scythe -a adapters.fasta SRR957824_500K_R2.fastq.gz -o SRR957824_adapt_R2.fastq 
```

:question:
What adapters do you use?

## Sickle

Most modern sequencing technologies produce reads that have deteriorating quality towards the 3'-end and some towards the 5'-end as well.
Incorrectly called bases in both regions negatively impact assembles, mapping, and downstream bioinformatics analyses.

We will trim each read individually down to the good quality part to keep the bad part from interfering with downstream applications.

To do so, we will use sickle. Sickle is a tool that uses sliding windows along with quality and length thresholds to determine when quality is sufficiently low to trim the 3'-end of reads and also determines when the quality is sufficiently high enough to trim the 5'-end of reads. It will also discard reads based upon a length threshold.

To run sickle

```bash
sickle pe -f SRR957824_adapt_R1.fastq -r SRR957824_adapt_R2.fastq \
    -t sanger -o SRR957824_trimmed_R1.fastq -p SRR957824_trimmed_R2.fastq \
    -s /dev/null -q 25
```

which should output something like

```
PE forward file: SRR957824_trimmed_R1.fastq
PE reverse file: SRR957824_trimmed_R2.fastq

Total input FastQ records: 1000000 (500000 pairs)

FastQ paired records kept: 834570 (417285 pairs)
FastQ single records kept: 13263 (from PE1: 11094, from PE2: 2169)
FastQ paired records discarded: 138904 (69452 pairs)
FastQ single records discarded: 13263 (from PE1: 2169, from PE2: 11094)
```

## FastQC again

Run fastqc again on the filtered reads

```bash
fastqc SRR957824_trimmed_R1.fastq SRR957824_trimmed_R2.fastq
```

and look at the reports

- [SRR957824_trimmed_R1_fastqc.html](data/fastqc/SRR957824_trimmed_R1_fastqc.html)
- [SRR957824_trimmed_R2_fastqc.html](data/fastqc/SRR957824_trimmed_R2_fastqc.html)

## MultiQC

[MultiQC](http://multiqc.info) is a tool that aggreagtes results from several popular QC bioinformatics software into one html report.

Let's run MultiQC in our current directory

```bash
multiqc .
```

You can download the report or view it by clicking on the link below

- [multiqc_report.html](data/fastqc/multiqc_report.html)

:question:
What did the trimming do to the per-base sequence quality, the per sequence quality scores and the sequence length distribution?