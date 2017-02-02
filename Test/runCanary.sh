#!/bin/bash
#
#	Run Canary for testing
#

##
##	Set to Canary install directory
##
CANARY_HOME=/usr/local/dev/Canary

$CANARY_HOME/bin/Canary	\
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
