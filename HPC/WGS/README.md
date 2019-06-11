<!-- ABOUT THE PROJECT -->
## Portable Pipeline for Whole Exome and Genome Sequencing

This is a portable workflow for converting the raw sequencing data in FASTQ format to results usable in a clinical environment. 
In combination with the Forome Anfisa Variant Curation Tool, it can yield a clinically actionable report.

### Prerequisites

This workflow uses the following tools:

 * [GATK4](https://software.broadinstitute.org/gatk/) 
 * [Samtools](http://www.htslib.org/)
 * [Sambamba](http://lomereiter.github.io/sambamba/)
 * [bcftools](https://samtools.github.io/bcftools/bcftools.html)
 * [bgzip](http://www.htslib.org/doc/bgzip.html)
 * [tabix](http://www.htslib.org/doc/tabix.html)
 * [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)
 * Python 2.7.9
 * [BWA](http://bio-bwa.sourceforge.net/)
 * [NovoCaller](https://github.com/bgm-cwg/novoCaller)
 
## Running

This workflow can be run with [Cromwell](https://github.com/broadinstitute/cromwell). You can run this workflow with the following command:

```sh
java -Dconfig.file=cromwell.config -jar cromwell.jar run wgs_main.wdl --inputs inputs.json
```
 
