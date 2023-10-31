# Genome Assembly

This tutorial is based on tutorials distributed by **Hadrien Gourlé** (<https://www.hadriengourle.com/tutorials/>) and **Michael Schatz - appliedgenomics course** (<http://schatz-lab.org/appliedgenomics2023/assignments/assignment2/>).

## Let's get started!

In this practical we will perform the assembly of *Mycoplasma genitalium*, a bacterium published in 1995 by Fraser et al in Science.

Click the **Start** button to move to the next step.

## Preparing the environment

Let's install `Miniconda` to manage easily install tools in your virtual machine.

`Miniconda` is a free minimal installer for **conda** that includes only `conda`, Python and a small number of other useful packages. `Miniconda` allows you to create a minimal self contained Python installation, and then use the `conda` command to install additional packages.

Run the following commands to install `Miniconda`:
```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
```

Reload your shell to enable the `conda` command
```bash
bash
```

Add `bioconda` and `conda-forge` to your channels
```bash 
conda config --add channels bioconda
conda config --add channels conda-forge
```

Install `mamba` to your `base` environment
```bash
conda install -y mamba
```

All of the tools needed for this assignment (except `SPAdes` and `BUSCO`) can be installed in the same environment using:
```bash
mamba create -n genome_assembly_env samtools bowtie2 quast megahit pilon
```

Let's now create a separate environment for `SPAdes`
```bash
mamba create -n spades_env spades
```

Finally another environment for `BUSCO`
```bash
mamba create -n busco_env python=3.7 busco=5.2.2
conda activate busco_env
```

## Downloading the data

The raw data were deposited at the European Nucleotide Archive (ENA), under the accession number `ERR486840`.
You can go to the ENA [website](http://www.ebi.ac.uk/ena) and search for the run with the accession `ERR486840`.

First create a `data/` directory in your home folder
```bash
mkdir ~/data
cd ~/data
```

now let's download the subset
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR486/ERR486840/ERR486840_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR486/ERR486840/ERR486840_2.fastq.gz
```

The first step in any sequencing project is to check the quality of the sequence data. 
This is important because it will determine the type of analysis that can be performed and the confidence that can be placed in the results.

Luckly, the files deposited in ENA were already trimmed, so we do not have to trim ourselves!

**Optional:** You can follow the instructions in the [Fastq Quality Control (QC) tutorial](https://github.com/RushAlz/IAMSPE-CS31-Genomica_Computacional/blob/main/tutorials/qc.md) to check/perform the quality control of the data.

## de novo assembly

We will use two tools to assemble our bacterium: `Megahit` and `SPAdes`

First create another folder to work on processed data
```bash
mkdir ~/work
cd ~/work
```

Create links to the `fastq` files
```bash
ln -s ../data/*.fastq.gz .
```

Starting with `SPAdes`. First activate the `conda` environment.
```bash
conda activate spades_env
```

Normally spades would try several values of k-mer and merge the results together, but here we will force it to just use `k=21` to save time. The assembly should take a few minutes.
```bash
spades.py --isolate -1 ERR486840_1.fastq.gz -2 ERR486840_2.fastq.gz -o m_genitalium_spades -t 2 -k 21
```

Check the results in the folder
```bash
ls -lh ~/work/m_genitalium_spades
```

How many contigs were produced?
```bash
grep -c '>' ~/work/m_genitalium_spades/contigs.fasta
```

Now `Megahit`. Activate the conda environment first.
```bash
conda activate genome_assembly_env
```

Run `Megahit`
```bash
megahit -1 ERR486840_1.fastq.gz -2 ERR486840_2.fastq.gz -o m_genitalium_megahit
```

Check the results in the folder
```bash
ls -lh ~/work/m_genitalium_megahit
```

Check how many contigs were produced
```bash
grep -c '>' ~/work/m_genitalium_megahit/final.contigs.fa
```

**Question:** What tool resulted in the better assembly?

## Quality of the Assembly

A single metric that evaluates the quality of an assembly does not exist, especially in the absence of a reference genome and in the case of a metagenomic assembly. However, we can examine some general aspects, such as `N50` (the length for which the collection of all contigs of that length or longer covers at least half of the total assembly size) or largest contig, or fraction of reads that map to your assembly, or number of genes detected, and more. But: 1) these metrics lack any context to indicate if they are "good" or not unless we compare different assemblies of the same data; and 2) it is possible that an assembly with "worse" overall summary statistics than another assembly may actually perform better for our specific purpose than the assembly with more favorable summary statistics. Therefore, we have to consider what our objectives are, and understand that selecting the "best" assembly from our data is not necessarily a simple or clear-cut task. Having a reference genome like we do in this case simplifies things as we will see next, and then we will discuss some options if we did not have a reference.

`QUAST` is a software evaluating the quality of genome assemblies by computing various metrics.

Run `QUAST` on both assemblies for comparison:
```bash
cd ~/work
quast.py -o quast_report \
  -l "SPAdes,Megahit" \
  m_genitalium_spades/contigs.fasta \
  m_genitalium_megahit/final.contigs.fa
```

and take a look at the text report
```bash
cat quast_report/report.txt
```

You also download the html(`quast_report/report.html`) report and open it in your own web browser.

**Question:** N50: length for which the collection of all contigs of that length or longer covers at least 50% of assembly length

**Question:** How well does the assembly total consensus size and coverage correspond to your earlier estimation?

**Question:** How many contigs in total did the assembly produce?

**Question:** What is the N50 of the assembly? What does this mean?

## Fixing misassemblies

`Pilon` is a software tool which can be used to automatically improve draft assemblies. It attempts to make improvements to the input genome, including:

* Single base differences
* Small Indels
* Larger Indels or block substitution events
* Gap filling
* Identification of local misassemblies, including optional opening of new gaps

`Pilon` then outputs a `FASTA` file containing an improved representation of the genome from the read data and an optional `VCF` file detailing variation seen between the read data and the input genome.

The result of the assembly is in the directory `m_genitalium` under the name `final.contigs.fa`

Let's first create another folder and make a copy of the assembly results using `Megahit`.
```bash
mkdir -p ~/work/pilon
cd ~/work/pilon
ln -s ~/work/m_genitalium_megahit/final.contigs.fa m_genitalium.fasta
```

Before running `Pilon` itself, we have to align our reads against the assembly. This will take some minutes...
```bash
bowtie2-build m_genitalium.fasta m_genitalium
bowtie2 -x m_genitalium -1 ../ERR486840_1.fastq.gz -2 ../ERR486840_2.fastq.gz | \
    samtools view -bS -o m_genitalium.bam
samtools sort m_genitalium.bam -o m_genitalium.sorted.bam
samtools index m_genitalium.sorted.bam
```

then run `Pilon`
```bash
pilon --genome m_genitalium.fasta --frags m_genitalium.sorted.bam --output m_genitalium_improved
```

Compare results using `QUAST`
```bash
quast.py -o quast_report \
  -l "Megahit,Megahit_Improved" \
  m_genitalium.fasta \
  m_genitalium_improved.fasta
```

and take a look at the text report
```bash
cat quast_report/report.txt
```

## Assembly Completeness

Although `QUAST` output a range of metric to assess how contiguous our assembly is, having a long `N50` does not guarantee a good assembly: it could be riddled by misassemblies!

We will run `BUSCO` to try to find marker genes in our assembly. Marker genes are conserved across a range of species and finding intact conserved genes in our assembly would be a good indication of its quality

First let's create another folder, create a link to the assembly and activate the `conda` environment
```bash
mkdir -p ~/work/busco
cd ~/work/busco
ln ~/work/m_genitalium_megahit/final.contigs.fa m_genitalium.fasta
conda activate busco_env
```

then we can run `BUSCO` with
```bash
busco -i m_genitalium.fasta -l bacteria_odb10 -o busco_genitalium -m genome
```

**Question:** How many marker genes has `BUSCO` found?

## Conclusion

<walkthrough-conclusion-trophy></walkthrough-conclusion-trophy>

Congratulations!

You've performed *de novo* genome assembly from sequencing data!
