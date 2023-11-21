# Mapping genetic variants with GATK

In this tutorial we will learn how to map germline genetic variants to a reference genome. 
We will follow the [GATK best practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels)

For purposes of demonstration, in this tutorial sequencing files were reduced to a small subset of reads. 
Also the reference genome has been restricted to a small region of chromosome 22.

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
mamba create -y -n gatk_env -c conda-forge -c bioconda -c defaults gatk4 samtools bcftools bwa  
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

To start let's decompress the reference genome file:
```bash
cd tutorials/files/reference_genome
gunzip hg38.chr22.fasta.gz
```

Now, let's set some variables to make it easier to run the commands in this tutorial:
```bash
FASTQ_FOLDER="~/cloudshell_open/IAMSPE-CS31-Genomica_Computacional/tutorials/files/samples_reads/"
RESOURCES_FOLDER="~/cloudshell_open/IAMSPE-CS31-Genomica_Computacional/tutorials/files/resources/"
REF="~/cloudshell_open/IAMSPE-CS31-Genomica_Computacional/tutorials/files/reference_genome/hg38.chr22.fasta"
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

You can check the contents of the folder to see the new files created by `bwa`:
```bash
ls -alh ${REF}*
```

Now, create a folder to store alignment results (note that we are using a variable to save the folder path):
```bash
ALIGN_FOLDER=${PWD}/"align"
mkdir -p ${ALIGN_FOLDER}
```

Mapping reads to the reference using BWA-MEM
```bash
## ALIGN RAW READS
mkdir -p ${ALIGN_FOLDER}
bwa mem -M -t 2 \
  -R "@RG\tID:tiny\tSM:tiny\tPL:ILLUMINA" \
  ${REF} ${FASTQ_FOLDER}/tiny_R1.fastq.gz ${FASTQ_FOLDER}/tiny_R2.fastq.gz > ${ALIGN_FOLDER}/tiny.sam
```

Convert sam to bam
```bash
samtools view -bS ${ALIGN_FOLDER}/tiny.sam > ${ALIGN_FOLDER}/tiny.bam
```

Sort and index bam file
```bash
sambamba sort -t 2 -o ${ALIGN_FOLDER}/tiny.sort.bam ${ALIGN_FOLDER}/tiny.bam
```  

## GATK MarkDuplicates

Create another folder to store metrics and technical data:
```bash
DATA_FOLDER=${PWD}/"data"
mkdir -p ${DATA_FOLDER}
```

Before starting with `GATK`. The workflow requires a dictionary file need to be created. 
This contains information about the reference genome.
```bash
gatk CreateSequenceDictionary -R ${REF} 
```

Now, you are ready to run `MarkDuplicates` to flag PCR duplicates. 
```bash
gatk MarkDuplicates \
  -I ${ALIGN_FOLDER}/tiny.sort.bam \
  -O ${ALIGN_FOLDER}/tiny.sort.dedup.bam \
  -M ${DATA_FOLDER}/tiny.mkd.metric
```

## BaseRecalibrator

Now let's run `BaseRecalibrator` to generate a recalibration table.
Note that we are using the `dbsnp_138.hg38.chr22.vcf.gz` file as a known sites.

Fire we build the model:
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

```bash
gatk CollectAlignmentSummaryMetrics \
  -R ${REF} \
  -I ${ALIGN_FOLDER}/tiny.sort.dedup.bqsr.bam \
  -O ${DATA_FOLDER}/alignment_metrics.txt
  
gatk CollectInsertSizeMetrics \
  -I ${ALIGN_FOLDER}/tiny.sort.dedup.bqsr.bam \
  -O ${DATA_FOLDER}/insert_size_metrics.txt \
  -H ${DATA_FOLDER}/insert_size_histogram.pdf
```

## Call Variants - gatk haplotype caller

```bash
mkdir -p ${RESULTS_FOLDER}

gatk HaplotypeCaller \
  -R ${REF} \
  -I ${ALIGN_FOLDER}/tiny.sort.dedup.bqsr.bam \
  -O ${RESULTS_FOLDER}/raw_variants.vcf
```

Extract SNPs:
```bash
gatk SelectVariants -R ${REF} \
  -V ${RESULTS_FOLDER}/raw_variants.vcf \
  --select-type SNP \
  -O ${RESULTS_FOLDER}/raw_snps.vcf
```

Extract INDELS:
```bash
gatk SelectVariants -R ${REF} \
  -V ${RESULTS_FOLDER}/raw_variants.vcf \
  --select-type INDEL \
  -O ${RESULTS_FOLDER}/raw_indels.vcf
```

## Filter Variants

Filter SNPs
```bash
gatk VariantFiltration \
	-R ${REF} \
	-V ${RESULTS_FOLDER}/raw_snps.vcf \
	-O ${RESULTS_FOLDER}/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"
```

Filter INDELS
```bash
gatk VariantFiltration \
	-R ${REF} \
	-V ${RESULTS_FOLDER}/raw_indels.vcf \
	-O ${RESULTS_FOLDER}/filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"
```

Select Variants that PASS filters:
```bash
gatk SelectVariants \
	--exclude-filtered \
	-V ${RESULTS_FOLDER}/filtered_snps.vcf \
	-O ${RESULTS_FOLDER}/analysis-ready-snps.vcf

gatk SelectVariants \
	--exclude-filtered \
	-V ${RESULTS_FOLDER}/filtered_indels.vcf \
	-O ${RESULTS_FOLDER}/analysis-ready-indels.vcf
```

To exclude variants that failed genotype filters:
```bash
cat ${RESULTS_FOLDER}/analysis-ready-snps.vcf |grep -v -E "DP_filter|GQ_filter" > ${RESULTS_FOLDER}/analysis-ready-snps-filteredGT.vcf
cat ${RESULTS_FOLDER}/analysis-ready-indels.vcf | grep -v -E "DP_filter|GQ_filter" > ${RESULTS_FOLDER}/analysis-ready-indels-filteredGT.vcf
```

## Annotate Variants with GATK4 Funcotator

We can annotate the variant mapped with relevant information using `GATK Funcotator`. First we need to download the data sources:
```bash
cd ${RESOURCES_FOLDER}
gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download
```

A folder named `funcotator_dataSources.v1.7.20200521g` with several files will be created. 
Let's also decompress `gnomad` files to enable querying to their database:
```bash
cd ${RESOURCES_FOLDER}/"funcotator_dataSources.v1.7.20200521g"
tar -zxf gnomAD_exome.tar.gz
tar -zxf gnomAD_genome.tar.gz
```

Now go back to the analysis folder
```bash
cd ~/analysis
```

Then Annotate SNPs using `Funcotator`:
```bash
gatk Funcotator \
	--variant ${RESULTS_FOLDER}/analysis-ready-snps-filteredGT.vcf \
	--reference ${REF} \
	--ref-version hg38 \
	--data-sources-path ${RESOURCES_FOLDER}/funcotator_dataSources.v1.7.20200521g \
	--output ${RESULTS_FOLDER}/analysis-ready-snps-filteredGT-functotated.vcf \
	--output-file-format VCF 
```

Same thing for INDELS:
```bash
gatk Funcotator \
	--variant ${RESULTS_FOLDER}/analysis-ready-indels-filteredGT.vcf \
	--reference ${REF} \
	--ref-version hg38 \
	--data-sources-path ${RESOURCES_FOLDER}/funcotator_dataSources.v1.7.20200521g \
	--output ${RESULTS_FOLDER}/analysis-ready-indels-filteredGT-functotated.vcf \
	--output-file-format VCF
```

Finally, let's merge these two files into one:
```bash
gatk MergeVcfs \
	-I ${RESULTS_FOLDER}/analysis-ready-snps-filteredGT-functotated.vcf \
	-I ${RESULTS_FOLDER}/analysis-ready-indels-filteredGT-functotated.vcf \
	-O ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT-functotated.vcf
```

Compress and index the file:
```bash
bgzip ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT-functotated.vcf
tabix -p vcf ${RESULTS_FOLDER}/analysis-ready-variants-filteredGT-functotated.vcf.gz
```

## Visualizing the annotation using IGV

`IGV` is a graphical Java program to browse genomes. It supports many different data formats, including `BAM`, `BED`, `GFF`, `VCF` and `WIG`. It was developed by the Broad Institute and is available for free download at <http://www.broadinstitute.org/software/igv/>.

Now, copy the resulting `VCF` file with the variants and annotations to your local computer. 

Here we will use the `IGV Web App`, but you can also download the standalone version and run locally.

Got to this [link](https://igv.org/app/). Make the genome is set to `hg38`.

Now load variants file. Click on `Tracks` then `Local File` and select the `analysis-ready-variants-filteredGT-functotated.vcf.gz` file you just copied. 

Now explore the variants in the chromossome 22. You can zoom in and out using the mouse wheel. You can also search for specific genes using the search bar.

## Conclusion

<walkthrough-conclusion-trophy></walkthrough-conclusion-trophy>

Congratulations!

You've performed genetic variant detection from sequencing data!
