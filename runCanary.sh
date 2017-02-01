#!/bin/bash
#
#	Run Canary for testing
#

#bin/Canary	--mutalyzer 'https://mutalyzer.nl' \
bin/Canary	\
		--mutalyzer  'https://vmts-mutalyzer1.unix.petermac.org.au' \
		--amplicon   Amplicon/amplicon.fa \
		--primers    Amplicon/amplicon.primers.tsv \
		--transcript etc/transcript.tsv \
		--columns    etc/cols \
		--reads      1 \
		--output     out.canary.tsv \
		--vcf        out.canary.vcf \
		--bam        out.canary.bam \
		--normalise  out.norm.vcf   \
		--tsv        out.norm.tsv   \
		Fastq/*R1_001.fastq.gz      \
		Fastq/*R2_001.fastq.gz
