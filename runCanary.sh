#!/bin/bash
#
#	Run Canary for testing
#

bin/Canary	\
		--mutalyzer 'https://mutalyzer.nl' \
		--amplicon   Amplicon/amplicon.fa \
		--primers    Amplicon/amplicon.primers.tsv \
		--transcript etc/transcript.tsv \
		--columns    etc/cols \
		--reads      10 \
		--complex	 \
		--output     out.canary.tsv \
		--vcf        out.canary.vcf \
		--bam        out.canary.bam \
		--normalise  out.norm.vcf   \
		--tsv        out.norm.tsv   \
		Fastq/*R1_001.fastq.gz      \
		Fastq/*R2_001.fastq.gz
