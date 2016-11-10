# Canary

Canary is a self contained amplicon NGS pipeline that takes a pair of zipped Fastq files and generates an annotated VCF file of variants. Canary correctly describes variants in HGVS nomenclature by 
making a number of calls to the Mutalyzer.nl website to render them with a Refseq transcript as an HGVSg, HGVSc and HGVSp and in their most 3' form. The preferred Refseq transcript for a gene
is taken from a gene to transcript mapping file which may edited by the user if required. Complex multi-nucleotide variants are correctly rendered as delins variants in their most parsimonious form
to save error-prone manual interpretation of these variants - a common cause of variant description error.

Clinical diagnostics is being transformed by the technology capable of analysing patient DNA at the nucleotide level.

Amplicon sequencing is a cost effective way of deeply sequencing a panel of genes of interest. Canary has been implemented for use with a number of clinical cancer diagnostic panels
at the Peter MacCallum Cancer Centre, the largest cancer only hospital in Australia.

The accuracy, turn around time and reproducability of clinical diagnostic sequencing relies heavily on bioinformatics pipelines to convert raw sequencing data into meaningful variants
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
There are a number of dependencies including the following
- Java JDK 1.7 from http://www.oracle.com/technetwork/java/javase/downloads/jdk7-downloads-1880260.html
- Groovy 2.1.9 from http://groovy-lang.org/download.html
- Gradle for building Canary (we use 1.10 at time of writing)
- Genome Analysis Toolkit (Currently GATK 3.3) and the Sting utility JAR (currently 2.1.8) available from here https://software.broadinstitute.org/gatk/download/
- The PathOS Core library available from here https://github.com/PapenfussLab/PathOS
- JNI wrapper to the striped Smith-Waterman alignment library SSW see https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library


## Contact
Ken Doig, Bioinformatics, Cancer Research Department, Data Scientist, Molecular Pathology Department
Peter MacCallum Cancer Centre, Victorian Comprehensive Cancer Centre Building
305 Grattan Street, Melbourne Victoria 3000 Australia
Ph: +61 411 225 178 Mail: ken.doig@petermac.org

