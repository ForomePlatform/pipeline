<!-- ABOUT THE PROJECT -->
## Portable Pipeline for Whole Exome and Genome Sequencing

This is a portable workflow for converting the raw sequencing data in FASTQ format to results usable in a clinical environment. 
In combination with the Forome Anfisa Variant Curation Tool, it can yield a clinically actionable report.

This workflow consists of 3 main components:
  * Upstream
  * Downstream
  * QC report

Upstream takes the raw .fastq files as an input and produces the aligned .bam files (by utilizing BWA) and g.vcf files. This component follows the GATk Best Practices.

Downstream takes all available g.vcf files (including files from all of previous runs), performs Genotyping and applies set of custom callers on the produced vcf. This component outputs a .vcf file that is ready to be upload in the curation tool.

Qc report gathers metrics and produces reports and plots that are helpful in the case analysis.

The overall algorithm (wgs_main.wdl) can be described as following:
  ```sh
  for each sample in case:
      call upstream
  call downstream
  call qc-report
  ```

### Prerequisites

This workflow uses the following tools:

 * [GATK4](https://software.broadinstitute.org/gatk/)
 * [GATK3.8](https://software.broadinstitute.org/gatk/download/archive)
 * [Samtools](http://www.htslib.org/)
 * [Sambamba](http://lomereiter.github.io/sambamba/)
 * [bcftools](https://samtools.github.io/bcftools/bcftools.html)
 * [bgzip](http://www.htslib.org/doc/bgzip.html)
 * [tabix](http://www.htslib.org/doc/tabix.html)
 * [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)
 * Python 2.7.9
   * [pandas 0.18.0](https://pandas.pydata.org/pandas-docs/version/0.18.0/)
   * [matplotlib](https://matplotlib.org/1.5.1/contents.html?fbclid=IwAR2fohH2ja4E6-jWl3_YaPxLckdv5OsWXVButxB8bhmU46nwSshnB5cMVtI)
 * [BWA](http://bio-bwa.sourceforge.net/)
 * [NovoCaller](https://github.com/bgm-cwg/novoCaller)
 * [vcf2geno](http://csg.sph.umich.edu/chaolong/LASER/index.html)
 * [trace](http://csg.sph.umich.edu/chaolong/LASER/index.html)
 * [plink](https://www.cog-genomics.org/plink2)
 * [bedtools](https://bedtools.readthedocs.io/en/latest/)
 

Due to the legacy pipeline code structure and the cluster configuration, there are ssome constraints on the directoris structure and storage format:
 * This workflow was created assuming that g.vcf files are split by chromosome for each sample.
 * All g.vcf files from all of priveous runs are stored under the same root directory. This root directory is passed as an input to get a list of all available g.vcf files. For example, here is the tree for the directory /data/gvcf_storage/:
    ```sh
    |-- case01_sample1
    |   |-- chr1.gvcf.gz
    |   |-- chr1.gvcf.gz.tbi
    |   |-- chr2.gvcf.gz
    |   |-- chr2.gvcf.gz.tbi
    |   |-- ..................
    |   |-- chrY.gvcf.gz
    |   |-- chrY.gvcf.gz.tbi
    |-- case02_sample1
    |   |-- chr1.gvcf.gz
    |   |-- chr1.gvcf.gz.tbi
    |   |-- chr2.gvcf.gz
    |   |-- chr2.gvcf.gz.tbi
    |   |-- ..................
    |   |-- chrY.gvcf.gz
    |   |-- chrY.gvcf.gz.tbi
    ```
    It contains 2 subdirectories with the g.vcf files for each chromosome from 2 samples. In this example, "/data/gvcf_storage/" will be passed as an input to the workflow.
 * The alligned .bam files from all of priveous runs are also stored under the same root directory. This root directory is passed as an input to get a list of all available .bam files (at the moment, this is required for the de-novo caller). For example, here is the tree for the directory /data/bam_storage/:
    ```sh
    |-- case01_sample1.bam
    |-- case01_sample1.bam.bai
    |-- case02_sample1.bam
    |-- case02_sample1.bam.bai
    ```
 * All executables and .jar files are stored under the same directory. This directory should be passed as "tools" input to the workflow. Here is the tree for the sample /home/tools/ directory:
    ```sh
    |-- bcftools
    |-- bedtools
    |-- bgzip
    |-- bwa
    |-- compound_het.out
    |-- denovo.out
    |-- gatk-4.1.2.0
    |   |-- gatk
    |   |-- gatk-completion.sh
    |   |-- gatkcondaenv.yml
    |   |-- GATKConfig.EXAMPLE.properties
    |   |-- gatk-package-4.1.2.0-local.jar
    |   |-- gatk-package-4.1.2.0-spark.jar
    |   |-- gatkPythonPackageArchive.zip
    |   |-- README.md
    |-- homrec.out
    |-- LASER-2.02
    |   |-- laser
    |   |-- license
    |   |-- pileup2seq
    |   |-- trace
    |   |-- vcf2geno
    |-- parallel
    |-- plink
    |-- pysam_prac.py
    |-- sambamba
    |-- samtools
    |-- split
    |-- tabix
    ```
 
## Running

This workflow can be run with [Cromwell](https://github.com/broadinstitute/cromwell). You can run this workflow in the Cromwells command line mode with the following command:

```sh
java -Dconfig.file=cromwell.config -jar cromwell.jar run wgs_main.wdl --inputs inputs.json
```
 
