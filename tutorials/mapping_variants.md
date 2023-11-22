## Mapping genetic variants with GATK

In this tutorial we will learn how to map germline genetic variants of a human DNA sample compared to the reference genome. 
We will follow the [GATK best practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels)

For purposes of demonstration, in this tutorial sequencing files were reduced to a small subset of reads. 
Also the reference genome has been restricted to chromosome 22.

References from which some of the tutorial was based on: 
* [Bioinformagician](https://github.com/kpatel427/YouTubeTutorials/blob/main/variant_calling.sh)
* [sib-swiss](https://sib-swiss.github.io/NGS-variants-training/2023.2/day2/annotation/)
* [Joe McGirr](https://joemcgirr.github.io/files/code_tutorials/my_genome/SnpEFF.html)

Click the **Start** button to move to the next step.

## Preparing the environment

Let's install `Miniconda` to manage and easily install tools in your virtual machine.

`Miniconda` is a free minimal installer for **conda** that includes only `conda`, Python, and a small number of other useful packages. 
`Miniconda` allows you to create a minimal, self-contained Python installation and then use the `conda` command to install additional packages.

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

Create a new environment and install `GATK` and other requirements using:
```bash
mamba create -y -n gatk_env -c conda-forge -c bioconda -c defaults gatk4 samtools sambamba bwa R snpeff snpsift
```

Activate the environment with:
```bash
conda activate gatk_env
```

## Let's get started!

The Github repository for this tutorial contains a folder called `files` with the data required to run the tutorial.

Let's check the contents of the `files` folder:
```bash
ls -alh tutorials/files/*
```

Now, let's set some variables to make it easier to run the commands in this tutorial:

Go the the `files` folder:
```bash
cd tutorials/files
```

Create the following variables:
```bash
FASTQ_FOLDER="${PWD}/samples_reads/"
RESOURCES_FOLDER="${PWD}/resources/"
REF="${PWD}/reference_genome/hg38.chr22.fasta"
```

To start let's decompress the reference genome file:
```bash
cd reference_genome
gunzip hg38.chr22.fasta.gz
```

Finally, let's set a folder to work on the analysis and save the results:
```bash
mkdir -p ~/analysis
cd ~/analysis
```

## Align reads to the reference genome

The first step is to map the reads to the reference genome. But before we do that, we need to index the reference genome.
```bash
bwa index ${REF}
```

We also need `samtools` index:
```bash
samtools faidx ${REF}
```

You can check the contents of the folder to see the new files created by `bwa`:
```bash
ls -alh ${REF}*
```

Now, create a folder to store alignment results (note that we are using a variable to save the folder path):
```bash
ALIGN_FOLDER=${PWD}/"align"
mkdir -p ${ALIGN_FOLDER}
```

Mapping reads to the reference using `BWA-MEM`
```bash
bwa mem -M \
  -t 2 \
  -R "@RG\tID:tiny\tSM:tiny\tPL:ILLUMINA" \
  ${REF} \
  ${FASTQ_FOLDER}/tiny_R1.fastq.gz \
  ${FASTQ_FOLDER}/tiny_R2.fastq.gz > ${ALIGN_FOLDER}/tiny.sam
```

Convert `sam` to `bam`
```bash
samtools view -bS ${ALIGN_FOLDER}/tiny.sam > ${ALIGN_FOLDER}/tiny.bam
```

Sort and index `bam` file
```bash
sambamba sort -t 2 -o ${ALIGN_FOLDER}/tiny.sort.bam ${ALIGN_FOLDER}/tiny.bam
```  

## GATK MarkDuplicates

This tool locates and tags duplicate reads in a `BAM` or `SAM` file, where duplicate reads are defined as originating from a single fragment of DNA. 
Duplicates can arise during sample preparation e.g. library construction using PCR. 
Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. 
These duplication artifacts are referred to as optical duplicates.

First, create another folder to store metrics and technical data:
```bash
DATA_FOLDER=${PWD}/"data"
mkdir -p ${DATA_FOLDER}
```

Before starting with `GATK`. The workflow requires a dictionary file need to be created. 
This contains information about the reference genome.
```bash
gatk CreateSequenceDictionary -R ${REF} 
```

Now, you are ready to run `MarkDuplicates` to flag duplicated reads:
```bash
gatk MarkDuplicates \
  -I ${ALIGN_FOLDER}/tiny.sort.bam \
  -O ${ALIGN_FOLDER}/tiny.sort.dedup.bam \
  -M ${DATA_FOLDER}/tiny.mkd.metric
```

## BaseRecalibrator

`BaseRecalibrator` generates a recalibration table for Base Quality Score Recalibration (BQSR).

The base recalibration process involves two key steps: first the `BaseRecalibrator` tool builds a model of covariation based on the input data and a set of known variants, producing a recalibration file; then the `ApplyBQSR` tool adjusts the base quality scores in the data based on the model, producing a new BAM file.

Tables are based on specified covariates, calculated by-locus, operating only at sites that are in the known sites VCF. 
ExAc, gnomAD, or dbSNP resources can be used as known sites of variation. 
We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. 
Since there is a large amount of data one can then calculate an empirical probability of error given the particular covariates seen at this site, where `p(error) = num mismatches / num observations`.
The output file is a table (of the several covariate values, num observations, num mismatches, empirical quality score).

Note that we are using the `dbsnp_138.hg38.chr22.vcf.gz` file as a known sites.

First, we build the model:
```bash
gatk BaseRecalibrator \
  -I ${ALIGN_FOLDER}/tiny.sort.dedup.bam \
  -R ${REF} \
  --known-sites ${RESOURCES_FOLDER}/dbsnp_138.hg38.chr22.vcf.gz \
  -O ${DATA_FOLDER}/recal_data.table
```

Then we apply the model to adjust the base quality scores:
```bash
gatk ApplyBQSR \
  -I ${ALIGN_FOLDER}/tiny.sort.dedup.bam \
  -R ${REF} \
  --bqsr-recal-file ${DATA_FOLDER}/recal_data.table \
  -O ${ALIGN_FOLDER}/tiny.sort.dedup.bqsr.bam 
```

## Collect Alignment & Insert Size Metrics

`CollectAlignmentSummaryMetrics` (Picard) produces a summary of alignment metrics from a `SAM` or `BAM` file. 
This tool takes a `SAM/BAM` file input and produces metrics detailing the quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters.
```bash
gatk CollectAlignmentSummaryMetrics \
  -R ${REF} \
  -I ${ALIGN_FOLDER}/tiny.sort.dedup.bqsr.bam \
  -O ${DATA_FOLDER}/alignment_metrics.txt
```

`CollectInsertSizeMetrics` does not estimate but (as by the name) collects insert size metrics, which is nothing different than parsing the TLEN field from the `BAM` file. 
It is intended to only be used for paired-end data, not single-end, as single-end data do not have a value in that field in the `BAM` file.
```bash
gatk CollectInsertSizeMetrics \
  -I ${ALIGN_FOLDER}/tiny.sort.dedup.bqsr.bam \
  -O ${DATA_FOLDER}/insert_size_metrics.txt \
  -H ${DATA_FOLDER}/insert_size_histogram.pdf
```

## Call Variants - gatk haplotype caller

The `HaplotypeCaller` is capable of calling `SNPs` and `indels` simultaneously via local *de-novo* assembly of haplotypes in an active region. 
In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region.

First lets create a folder to store the results:
```bash
RESULTS_FOLDER="${PWD}/results"
mkdir -p ${RESULTS_FOLDER}
```

Now, run the `HaplotypeCaller`:
```bash
gatk HaplotypeCaller \
  -R ${REF} \
  -I ${ALIGN_FOLDER}/tiny.sort.dedup.bqsr.bam \
  -O ${RESULTS_FOLDER}/raw_variants.vcf
```

Results are salved in the file `raw_variants.vcf`, we can then extract `SNPs` only:
```bash
gatk SelectVariants -R ${REF} \
  -V ${RESULTS_FOLDER}/raw_variants.vcf \
  --select-type SNP \
  -O ${RESULTS_FOLDER}/raw_snps.vcf
```

Same thing for `indels`:
```bash
gatk SelectVariants -R ${REF} \
  -V ${RESULTS_FOLDER}/raw_variants.vcf \
  --select-type INDEL \
  -O ${RESULTS_FOLDER}/raw_indels.vcf
```

## Variant Filtering & evaluation (Hard filtering)

The developers of `gatk` strongly advise to do Variant Quality Score Recalibration (VQSR) for filtering SNPs and INDELs. 
However, this is not always possible. For example, in the case of limited data availability and/or in the case you are working with non-model organisms and/or in the case you are a bit lazy and okay with a number of false positives. Our dataset is too small to apply VQSR. We will therefore do hard filtering instead.

`VariantFiltration` tool is designed for hard-filtering variant calls based on certain criteria. 
Records are hard-filtered by changing the value in the FILTER field to something other than PASS. 
Filtered records will be preserved in the output unless their removal is requested in the command line.

Add flags to filter `SNPs`:
```bash
gatk VariantFiltration \
	-R ${REF} \
	-V ${RESULTS_FOLDER}/raw_snps.vcf \
	-O ${RESULTS_FOLDER}/filtered_snps.vcf \
  --filter-expression "QD < 2.0"              --filter-name "QD2" \
  --filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \
  --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
  --filter-expression "FS > 60.0"             --filter-name "FS60" \
  --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  --genotype-filter-expression "DP < 10"      --genotype-filter-name "DP_filter" \
  --genotype-filter-expression "GQ < 10"      --genotype-filter-name "GQ_filter"
```

Check the number of variants that passed the filters:
```bash
grep -v "^#" ${RESULTS_FOLDER}/filtered_snps.vcf | cut -f 7 | sort | uniq -c
```

Now, let's add filter for `INDELS`:
```bash
gatk VariantFiltration \
	-R ${REF} \
	-V ${RESULTS_FOLDER}/raw_indels.vcf \
	-O ${RESULTS_FOLDER}/filtered_indels.vcf \
  --filter-expression "QD < 2.0"                  --filter-name "QD2" \
  --filter-expression "QUAL < 30.0"               --filter-name "QUAL30" \
  --filter-expression "FS > 200.0"                --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0"    --filter-name "ReadPosRankSum-20" \
	--genotype-filter-expression "DP < 10" 	        --genotype-filter-name "DP_filter" \
	--genotype-filter-expression "GQ < 10"        	--genotype-filter-name "GQ_filter"
```

Now, we can easily select only variants that `PASS` filters. First `SNPs`:
```bash
gatk SelectVariants \
	--exclude-filtered \
	-V ${RESULTS_FOLDER}/filtered_snps.vcf \
	-O ${RESULTS_FOLDER}/analysis-ready-snps.vcf
```

Then, `indels`:
```bash
gatk SelectVariants \
	--exclude-filtered \
	-V ${RESULTS_FOLDER}/filtered_indels.vcf \
	-O ${RESULTS_FOLDER}/analysis-ready-indels.vcf
```

To exclude variants that failed the `genotype` filters, we can run:
```bash
cat ${RESULTS_FOLDER}/analysis-ready-snps.vcf | grep -v -E "DP_filter|GQ_filter" > ${RESULTS_FOLDER}/analysis-ready-snps-filteredGT.vcf
cat ${RESULTS_FOLDER}/analysis-ready-indels.vcf | grep -v -E "DP_filter|GQ_filter" > ${RESULTS_FOLDER}/analysis-ready-indels-filteredGT.vcf
```

Double-check the number of variants that passed the filters:
```bash
grep -v "^#" ${RESULTS_FOLDER}/analysis-ready-snps-filteredGT.vcf | cut -f 7 | sort | uniq -c
```

Finally, let's merge these two files into one:
```bash
gatk MergeVcfs \
	-I ${RESULTS_FOLDER}/analysis-ready-snps-filteredGT.vcf \
	-I ${RESULTS_FOLDER}/analysis-ready-indels-filteredGT.vcf \
	-O ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.vcf
```

## Annotate Variants (part 1)

Now, the last step is to annotate the variant mapped with relevant information.
We will use `snpEff`, a variant annotation and effect prediction tool. 
It annotates and predicts the effects of genetic variants.

First we need to download the data sources:
```bash
snpEff download -v GRCh38.99
```

You can run snpEff like so:
```bash
snpEff -Xmx4g -v GRCh38.99 ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.vcf > ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.ann.vcf
```

Now, we can check the missense variants, for example:
```bash
grep missense ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.ann.vcf
```

We can convert the vcf into a table, to easily check the variants:
```bash
gatk VariantsToTable \
	-V ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.ann.vcf \
	-F CHROM -F POS -F TYPE -F ID -F ANN -F LOF -F NMD -GF AD -GF DP -GF GQ -GF GT \
	-O ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.ann.table
```

Check the first lines:
```bash
head ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.ann.table
```

## Annotate Variants (part 2)

The annotated vcf output by SnpEFF has lots of information about how a variant influences molecular phenotypes (not necessarily disease phenotypes, which are explored below with SnpSift and Clinvar). 
Molecular effects are described by a sequence ontology term and associated with an estimate the magnitude of the functional impact.

SnpSift is automatically installed along with SnpEFF. SnpSift takes annotations from one vcf and transfers them to matching variants in another vcf.
I annotate my Nebula vcf using a vcf curated by Clinvar. The Clinvar vcf reports the clinical significance of each variant based on supporting literature. I use this to prioritize possible variants of concern (annotated as pathogenic). Download the Clinvar vcf.

```bash
cd ${RESOURCES_FOLDER}
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
tabix -p vcf clinvar.vcf.gz
```

Now go back to the analysis folder
```bash
cd ~/analysis
```

Annotate my vcf with Clinvar vcf and convert to tab delimited table with GATK.
```bash
SnpSift annotate ${RESOURCES_FOLDER}/clinvar.vcf.gz ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.vcf > ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.clinvar.vcf
```

Again, we can convert the vcf into a table:
```bash
gatk VariantsToTable \
	-V ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.clinvar.vcf \
	-F CHROM -F POS -F TYPE -F ID -F ALLELEID -F CLNDN -F CLNSIG -F CLNSIGCONF -F CLNSIGINCL -F CLNVC -F GENEINFO -GF AD -GF GQ -GF GT \
	-O ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.clinvar.table
```

Finally, we can check for pathogenic variants, by searching for the word `Pathogenic` or `Likely_pathogenic`. 
For this example data, there are no pathogenic variants mapped, but we can find some variants annotated as `Benign`:
```bash
grep "Benign" ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.clinvar.table
```

## Annotate Variants (part 3)

Alternatively, we can use `GATK Funcotator`. First we need to download the data sources:
```bash
cd ${RESOURCES_FOLDER}
gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download
```

A folder named `funcotator_dataSources.v1.7.20200521g` with several files will be created. 

Now go back to the analysis folder
```bash
cd ~/analysis
```

Then Annotate SNPs and INDELS using `Funcotator`:
```bash
gatk Funcotator \
	--variant ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT.vcf \
	--reference ${REF} \
	--ref-version hg38 \
	--data-sources-path ${RESOURCES_FOLDER}/funcotator_dataSources.v1.7.20200521g \
	--output ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT-functotated.vcf \
	--output-file-format VCF 
```

## Visualizing the VCF using IGV

`IGV` is a graphical Java program to browse genomes.
It supports many different data formats, including `BAM`, `BED`, `GFF`, `VCF` and `WIG`. 
It was developed by the Broad Institute and is available for free download at <http://www.broadinstitute.org/software/igv/>.

Now, copy the resulting `VCF` file with the variants and annotations to your local computer. 

Here we will use the `IGV Web App`, but you can also download the standalone version and run locally.

Got to this [link](https://igv.org/app/). Make the genome is set to `hg38`.

Now load variants file. Click on `Tracks` then `Local File` and select the `analysis-ready-variants-filteredGT-functotated.vcf.gz` file you just copied. 

Now explore the variants in the chromosome 22. You can zoom in and out using the mouse wheel. You can also search for specific genes using the search bar.

## Conclusion

<walkthrough-conclusion-trophy></walkthrough-conclusion-trophy>

Congratulations!

You've performed genetic variant detection from sequencing data!
