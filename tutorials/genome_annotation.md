# Genome Annotation

This tutorial is based on tutorials distributed by **Hadrien Gourlé** (<https://www.hadriengourle.com/tutorials/>).

## Let's get started!

After you have *de novo* assembled your genome sequencing reads into `contigs`, it is useful to know what genomic features are on those `contigs`.

The process of identifying and labeling those features is called `genome annotation`.

Click the **Start** button to move to the next step.

## Preparing the environment

Let's install `Miniconda` to manage and easily install tools in your virtual machine.

`Miniconda` is a free minimal installer for **conda** that includes only `conda`, Python, and a small number of other useful packages. `Miniconda` allows you to create a minimal, self-contained Python installation and then use the `conda` command to install additional packages.

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

Accept the `Terms of Service` for `conda`:
```bash
conda tos accept
```

Create a new environment and install `Prokka` using:
```bash
conda create -n prokka_env -c conda-forge -c bioconda -c defaults prokka 
```

## Input data

Prokka requires assembled contigs. You can prepare your working directory for this annotation tutorial.
```bash
mkdir ~/annotation
cd ~/annotation
```

You will download an improved assembly of *Mycoplasma genitalium* into your data directory:
```bash
curl -O -J -L https://osf.io/7eaky/download
```

You will also need a set of proteins specific to *Mycoplasma* for the annotation. Here is a file containing the *Mycoplasma* proteins retrieved from `Swiss-Prot` database (3041 sequences)
```bash
curl -O -J -L https://osf.io/xjm3n/download
```

## Running Prokka

`Prokka` is a `"wrapper"`; it collects together several pieces of software (from various authors), and so avoids "re-inventing the wheel".

`Prokka` finds and annotates features (both protein-coding regions and RNA genes, i.e., tRNA, rRNA) present in a sequence. `Prokka` uses a two-step process for the annotation of protein-coding regions: first, protein-coding regions on the genome are identified using [Prodigal](http://compbio.ornl.gov/prodigal/); second, the function of the encoded protein is predicted by similarity to proteins in one of many protein or protein domain databases.

`Prokka` is a software tool that can be used to annotate bacterial, archaeal, and viral genomes quickly, generating standard output files in `GenBank`, `EMBL` and `gff` formats.

More information about `Prokka` can be found [here](https://github.com/tseemann/prokka).

Activate your conda environment
```bash
conda activate prokka_env
```

```bash
prokka --outdir annotation --kingdom Bacteria \
--proteins uniprot_mycoplasma_reviewed.faa m_genetalium_improved.fasta
```

Once Prokka has finished, examine each of its output files.

-   The GFF and GBK files contain all of the information about the features annotated (in different formats.)
-   The .txt file contains a summary of the number of features annotated.
-   The .faa file contains the protein sequences of the genes annotated.
-   The .ffn file contains the nucleotide sequences of the genes annotated.

## Visualizing the annotation

`IGV` is a graphical Java program to browse genomes. It supports many different data formats, including `BAM`, `BED`, `GFF`, `VCF` and `WIG`. It was developed by the Broad Institute and is available for free download at <http://www.broadinstitute.org/software/igv/>.

`IGV` requires the genome fasta to be indexed. Let's do this using `sammtols`:

```bash
conda install -y samtools
```

```bash
samtools faidx m_genetalium_improved.fasta
```

Now, copy the genome assembly `.fasta`, the index `.fasta.fai` and the `.gff` file produced by `Prokka`.

Here we will use the `IGV Web App`, but you can also download the standalone version and run locally.

Got to this [link](https://igv.org/app/) and click on `Genome` then `Local File` and select the `.fasta` and `.fasta.fai` files you just copied. This will load the genome into the browser.

Now load the gene annotations. Click on `Tracks` then `Local File` and select the `.gff` file you just copied. This will load the annotations into the browser.

Now explore the genome annotations. You can zoom in and out using the mouse wheel. You can also search for specific genes using the search bar.

## Conclusion

<walkthrough-conclusion-trophy></walkthrough-conclusion-trophy>

Congratulations!

You've performed bacterial genome annotation from sequencing data!
