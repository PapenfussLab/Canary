# Canary

Canary is a self contained amplicon NGS pipeline that takes a pair of zipped Fastq files and generates an annotated VCF file of variants. Canary correctly describes variants in HGVS nomenclature by 
making a number of calls to the Mutalyzer.nl website to render them with a Refseq transcript as an HGVSg, HGVSc and HGVSp and in their most 3' form. The preferred Refseq transcript for a gene
is taken from a gene to transcript mapping file which may edited by the user if required. Complex multi-nucleotide variants are correctly rendered as delins variants in their most parsimonious form
to save error-prone manual interpretation of these variants - a common cause of variant description error.

Clinical diagnostics is being transformed by the technology capable of analysing patient DNA at the nucleotide level.

Amplicon sequencing is a cost effective way of deeply sequencing a panel of genes of interest. Canary has been implemented for use with a number of clinical cancer diagnostic panels
at the Peter MacCallum Cancer Centre, the largest cancer only hospital in Australia.

The accuracy, turn around time and reproducibility of clinical diagnostic sequencing relies heavily on bioinformatics pipelines to convert raw sequencing data into meaningful variants
that can be curated for patient reporting.
Canary is a clinically tested tool which performs the key pipeline tasks of:
- alignment of amplicon reads, 
- calling of variants, 
- 3' shifting, 
- coalescing of MNPs (multi-nucleotide polymorphisms), 
- Refseq transcript selection and 
- rendering of variants into HGVS nomenclature. 

## Technology Platform
Canary takes advantage of many open-source and public Java libraries to implement an enterprise-grade application suitable for clinical use.

It is implemented in Groovy (a Java byte compatible language) and runs on any computer with Java installed.

## Installation
This repository contains a pre-built uber-jar bundling all dependencies for Canary. Cloning the repository will download the [uber-jar](https://github.com/PapenfussLab/Canary/blob/master/lib/Canary-all-1.0.0.jar) ready for running the [example data](https://github.com/PapenfussLab/Canary#example-data).

If you are more adventurous, the repository also contains all the dependencies (listed below) to build Canary. From the install directory, run:

	% export CANARY_HOME=/canary/target/dir   # replace with where you want to install Canary - must be absolute path - starts with '/'
	% gradle uploadArchives
	% PATH=$CANARY_HOME/bin:$PATH

There are a number of dependencies including the following
- Java JDK 1.7 from [Oracle](http://www.oracle.com/technetwork/java/javase/downloads/java-archive-downloads-javase7-521261.html)
- Groovy 2.1.9 from [here](http://groovy-lang.org/download.html)
- Gradle for building Canary (we use 1.10 at time of writing) [Gradle 1.10](https://services.gradle.org/distributions/gradle-1.10-bin.zip)
- Genome Analysis Toolkit (Currently GATK 3.3) and the Sting utility JAR (currently 2.1.8) available from [here](https://software.broadinstitute.org/gatk/download/)
- The PathOS Core library available from [PathosCore-all-1.3.jar](https://github.com/PapenfussLab/Canary/blob/master/repos/PathosCore-all-1.3.jar) and maintained [here](https://github.com/PapenfussLab/PathOS). 
- JNI wrapper to the striped Smith-Waterman alignment library SSW see [here](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)

## Example Data
The `Test` directory contains an example shell script `runCanary.sh` for running Canary against sample FASTQ reads files in the `Fastq` directory. These reads were generated on 
an Illumina MiSeq platform using the TruSeq assay, a 48 gene, targeted cancer amplicon panel from Illumina. The seqencing yielded 771,606 reads for this sample of which 93.5% were mapped to the amplicons at a coverage of ~1000X.

To run test data:

	% export CANARY_HOME=/canary/install/dir   # where you cloned the repository to
	% cd $CAN_DIR/Test
	% runCanary.sh

    2017-02-02 11:07:42,696 [main] INFO  org.petermac.pathos.pipeline.Canary - Canary [--mutalyzer, https://mutalyzer.nl, --amplicon, /usr/local/dev/Canary/Amplicon/amplicon.fa, --primers, /usr/local/dev/Canary/Amplicon/amplicon.primers.tsv, --transcript, /usr/local/dev/Canary/etc/transcript.tsv, --columns, /usr/local/dev/Canary/etc/cols, --reads, 10, --complex, --output, out.canary.tsv, --vcf, out.canary.vcf, --bam, out.canary.bam, --normalise, out.norm.vcf, --tsv, out.norm.tsv, /usr/local/dev/Canary/Fastq/14M6168_AACCCCTC-TAGACCTA_L001_R1_001.fastq.gz, /usr/local/dev/Canary/Fastq/14M6168_AACCCCTC-TAGACCTA_L001_R2_001.fastq.gz]
    2017-02-02 11:07:43,036 [main] INFO  org.petermac.pathos.pipeline.Canary - Loaded 14808 gene/transcripts from /usr/local/dev/Canary/etc/transcript.tsv
    Using PathOS Configuration File [/usr/local/dev/Canary/etc/canary.properties] PathOS Home [.]
    2017-02-02 11:07:43,068 [main] INFO  org.petermac.util.SmithWaterman - Loading SmithWaterman JNI Library from /usr/local/dev/Canary/lib/libsswjni.jnilib
    2017-02-02 11:07:43,110 [main] INFO  org.petermac.pathos.pipeline.Canary - Found Reads file: readlen=150 filesize=30495012 estimated reads=392429
    2017-02-02 11:07:43,119 [main] INFO  org.petermac.pathos.pipeline.Canary - Found 221 Amplicons
    2017-02-02 11:07:43,335 [main] INFO  org.petermac.pathos.pipeline.Canary - PIK3CA4_11.chr3.178936074.178936095_tile_1 has fwd rev matches to Off_target_1_PIK3CA4_11.chr3.178936074.178936095_tile_1-PIK3CA4_11.chr3.178936074.178936095_tile_1
    2017-02-02 11:07:43,336 [main] INFO  org.petermac.pathos.pipeline.Canary - PIK3CA12.chr3.178938860.178938860_tile_1 has fwd rev matches to Off_target_2_PIK3CA12.chr3.178938860.178938860_tile_1-PIK3CA12.chr3.178938860.178938860_tile_1
    2017-02-02 11:07:43,340 [main] INFO  org.petermac.pathos.pipeline.Canary - PTEN13.chr10.89720716.89720852_tile_1 has fwd rev matches to Off_target_4_PTEN13.chr10.89720716.89720852_tile_1-PTEN13.chr10.89720716.89720852_tile_1
    2017-02-02 11:07:43,342 [main] INFO  org.petermac.pathos.pipeline.Canary - Using 212 Amplicons
    2017-02-02 11:08:00,800 [main] WARN  org.petermac.pathos.pipeline.Canary - Maximum estimated read limit reached: 39242 read pairs
    2017-02-02 11:08:00,867 [main] INFO  org.petermac.pathos.pipeline.HGVS - Complex indel: chr17:g.7578373_7578395delinsGCTGCTCACCATCGCT
    2017-02-02 11:08:01,143 [main] INFO  org.petermac.pathos.pipeline.Canary - Found 46 events
    2017-02-02 11:08:01,170 [main] INFO  org.petermac.pathos.pipeline.Canary - Found 44 variants
    2017-02-02 11:08:01,807 [main] INFO  org.petermac.pathos.pipeline.Mutalyzer - Mutalyzer: No Proxy set
    2017-02-02 11:08:01,807 [main] INFO  org.petermac.pathos.pipeline.Mutalyzer - Mutalyzer: No Proxy set
    2017-02-02 11:08:01,829 [main] INFO  org.petermac.pathos.pipeline.HGVS - Complex indel: chr17:g.7578373_7578395delinsGCTGCTCACCATCGCT
    2017-02-02 11:08:01,833 [main] INFO  org.petermac.pathos.pipeline.MutalyzerUtil - Loaded 44 variants from out.canary.vcf
    2017-02-02 11:08:05,581 [main] INFO  org.petermac.pathos.pipeline.Mutalyzer - Job running: retrieve with: https://mutalyzer.nl/batch-job-result/batch-job-d8086e97-9a8d-49c2-9be3-59646a9e3210.txt
    2017-02-02 11:08:20,916 [main] INFO  org.petermac.pathos.pipeline.Mutalyzer - Job running: retrieve with: https://mutalyzer.nl/batch-job-result/batch-job-9bf359c8-a8cf-4656-97ee-83260d3c1dae.txt
    2017-02-02 11:08:34,345 [main] INFO  org.petermac.pathos.pipeline.MutalyzerUtil - Found 0 redesignated variants
    2017-02-02 11:08:34,401 [main] INFO  org.petermac.pathos.pipeline.HGVS - Complex indel: chr17:g.7578373_7578395delinsGCTGCTCACCATCGCT
    2017-02-02 11:08:34,521 [main] INFO  org.petermac.pathos.pipeline.HGVS - Complex indel: chr17:g.7578373_7578395delinsGCTGCTCACCATCGCT
    2017-02-02 11:08:34,622 [main] INFO  org.petermac.pathos.pipeline.MutalyzerUtil - convertVcf(out.canary.vcf): In 44 Out 44
    2017-02-02 11:08:34,656 [main] INFO  org.petermac.pathos.pipeline.Canary - Done: processed 313944 lines 39243 read pairs into out.canary.vcf in 52 seconds

Typical parameters for Canary:

	% Canary	\
		--mutalyzer 'https://mutalyzer.nl' \
		--amplicon   $CANARY_HOME/Amplicon/amplicon.fa \
		--primers    $CANARY_HOME/Amplicon/amplicon.primers.tsv \
		--transcript $CANARY_HOME/etc/transcript.tsv \
		--columns    $CANARY_HOME/etc/cols \
		--reads      10 \
		--complex	    \
		--output     out.canary.tsv \
		--vcf        out.canary.vcf \
		--bam        out.canary.bam \
		--normalise  out.norm.vcf   \
		--tsv        out.norm.tsv   \
		$CANARY_HOME/Fastq/*R1_001.fastq.gz   \
		$CANARY_HOME/Fastq/*R2_001.fastq.gz

## Configuration

There are a number of files that make Canary work. These are referenced in the `bin/Canary` shell script and the properties file at `etc/canary.properties`.

The location of logging can be configured by editing the properties file `lib/log4j.properties`.

When using the `--normalise` flag, Canary will need access to a full reference genome FASTA file for correctly annotation 3' shifted variants. The location of this file is set by the properties file located in `etc/canary.properties`. Update the `genome.path` property to point to the genome reference. These can be obtained from [1000g](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/) for various genome builds.

The amplicons generating the reads need to be described to Canary with two files. A fasta file (see `Amplicon/amplicon.fa`) containing the genomic sequences of the amplicons including primers. A second tab-delimited file describing the genomic position of the amplicons is required with the following format. (see `Amplicon/amplicon.primers.tsv`)

| Column | Description | Example |
| --- | --- | --- |
|1|Genomic coordinates of amplicon (including primers)|1:43814982-43815163|
|2|Length of start primer (bp)|24 |
|3|Length of end primer (bp)|26|
|4|Name of amplicon (Off target amplicons must start with "Off")|MPL1_2.chr1.43815008.43815009_tile_1|

## Contact
Ken Doig, Bioinformatics, Cancer Research Department, Data Scientist, Molecular Pathology Department
Peter MacCallum Cancer Centre, Victorian Comprehensive Cancer Centre Building
305 Grattan Street, Melbourne Victoria 3000 Australia
Ph: +61 411 225 178 Mail: ken.doig@petermac.org

