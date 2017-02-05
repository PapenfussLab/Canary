#!/bin/bash
#
#	Run Canary for testing
#

RUNDIR=`dirname $0`

Canary	\
		--mutalyzer 'https://mutalyzer.nl' \
		--amplicon   $RUNDIR/../Amplicon/amplicon.fa \
		--primers    $RUNDIR/../Amplicon/amplicon.primers.tsv \
		--transcript $CANARY_HOME/etc/transcript.tsv \
		--columns    $CANARY_HOME/etc/cols \
		--reads      10 \
		--complex	    \
		--output     out.canary.tsv \
		--vcf        out.canary.vcf \
		--bam        out.canary.bam \
		--normalise  out.norm.vcf   \
		--tsv        out.norm.tsv   \
		$RUNDIR/../Fastq/*R1_001.fastq.gz   \
		$RUNDIR/../Fastq/*R2_001.fastq.gz
