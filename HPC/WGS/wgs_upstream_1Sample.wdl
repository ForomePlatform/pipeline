workflow wgs_upstream
{
      #RESOURCES SECTION
      File dbSNP_vcf
      File dbSNP_vcf_idx
	
      File known_indels_sites_VCF
      File known_indels_sites_idx
      
      File ref_fasta
      File ref_dict
      File ref_fasta_index
      File ref_bwt
      File ref_sa
      File ref_amb
      File ref_ann
      File ref_pac
	
      #INPUT SECTION
      Array[Pair[File, File]] fastqs
      String sampleName
	  
      String base_name
      #String rg_string = "'@RG\\tID:@HS2000\\tPL:ILLUMINA\\tLB:bgm_lib\\tSM:" + base_name + sampleName + "'"	
      
      String bwa_commandline="bwa mem -R '@RG\\tID:" + base_name + sampleName + "\\tPL:ILLUMINA\\tLB:bgm_lib\\tSM:" + base_name + sampleName + "' -K 100000000 -v 3 -t $bwa_threads -Y $bash_ref_fasta"
      
      Array[String] scattered_calling_intervals = ["chrM", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
	
      #Optimization flags
      Int bwa_threads
      Int samtools_threads

      String tools
      String gatk_version = "4.1.0.0"

      String gvcfs_dir 

      call CreateSequenceGroupingTSV
      {
         input:
           ref_dict = ref_dict
      }

      call GetBwaVersion
      {
         input:
           tools = tools
      }
   
      #String key_sampleName         =sampleName
      #Pair[File, File] value_fastqs = (fastq1, fastq2)

      scatter(index in range(length(fastqs)))
      {
         #Int i = index
         Pair[File, File] fq_pair = fastqs[index]

         call SplitFQ_parallel as SplitFQ1
         {
            input:
              fq = fq_pair.left,
              suffix = "1",
              case_name = base_name + sampleName + "_" + index,
              cpu = 2,
              memory = "1 GB"
         }

         call SplitFQ_parallel as SplitFQ2
         {
            input:
              fq = fq_pair.right,
              suffix = "2",
              case_name = base_name + sampleName + "_" + index,
              cpu = 2,
              memory = "1 GB"
         }
      }

      call GatherFiles as GatherFQ1
      {
         input:
           files = SplitFQ1.fqs,
           suffix = "fq1"
      }
  
      call GatherFiles as GatherFQ2
      {
         input:
           files = SplitFQ2.fqs,
           suffix = "fq2"
      }
      
      scatter(index in range(length(GatherFQ1.o1)))
      {
          Int i = index
          Pair[File, File] fastq_pair = (GatherFQ1.o1[i], GatherFQ2.o1[i])
      
          call Fastq_to_uBAM
      	  {
      	    input:
      	      fastq1               = fastq_pair.left,
              fastq2               = fastq_pair.right,
              sampleName           = base_name + sampleName,
              output_uBAM_baseName = base_name + sampleName + "_" + i + ".u",
              cpu                  = 1,
              tools                = tools,
              memory               = "4096 MB",
              java_heap            = "-Xmx4096m",
              gatk_version         = gatk_version
	      }
	  
	      call SortSam_byQuery as Fastq_to_uBAM_sort
	      {
	         input:
	           input_bam           = Fastq_to_uBAM.uBAM,
	           output_bam_basename = "sorted_uBAM_" + base_name + sampleName + "." + i,
	           compressionLvl      = 5,
	           cpus                = 8,
                   sambamba            = tools + "/my_sambamba/sambamba",
	           flag                = "-n"
	      }
	  
	      call BwaMem
	      {
	         input:
	           fastq1               = fastq_pair.left,
	           fastq2               = fastq_pair.right,
	           bwa_commandline      = bwa_commandline,
	           output_bam_basename  = base_name + "." + sampleName + ".unmerged",
	           ref_fasta            = ref_fasta,
	           ref_dict             = ref_dict,
                   ref_fasta_index      = ref_fasta_index,
                   ref_bwt              = ref_bwt,
                   ref_sa               = ref_sa,
                   ref_amb              = ref_amb,
                   ref_ann              = ref_ann,
                   ref_pac              = ref_pac,
	           cpu                  = bwa_threads,
                   memory               = "16 GB",
                   tools                = tools
	      }
	  
	      call MergeBamAlignment
	      {
	         input:
	           alignedBAM          = BwaMem.output_bam,
	           uBAM                = Fastq_to_uBAM_sort.output_bam,
	           uBAM_index          = Fastq_to_uBAM_sort.output_bam_index,
	           output_bam_basename = base_name + "." + sampleName + "_" + i + "merged.aligned",
	           ref_fasta           = ref_fasta,
	           ref_dict            = ref_dict,
                   ref_fasta_index     = ref_fasta_index,
                   ref_bwt             = ref_bwt,
                   ref_sa              = ref_sa,
                   ref_amb             = ref_amb,
                   ref_ann             = ref_ann,
                   ref_pac             = ref_pac,
	           bwa_version         = GetBwaVersion.version,
	           bwa_mem_commandline = bwa_commandline,
                   tools               = tools,
                   gatk_version        = gatk_version,
                   tmp_dir             = "/net/home/isaevt/tmp_dir"
	      }

              #call SortSam_byQuery as MergeAlignmentSort
              #{
              #   input:
              #     input_bam           = MergeBamAlignment.output_bam,
              #     output_bam_basename = "sorted_aligned_" + base_name + sampleName + "." + i,
              #     compressionLvl      = 5,
              #     cpus                = 8,
              #     sambamba            = tools + "/my_sambamba/sambamba",
              #     flag                = "-n"
              #}
	  
      }
	
      #call MarkDuplicates
      #{
      #   input:
      #     input_bams           = MergeAlignmentSort.output_bam,
      #     output_bam_basename  = base_name + "." + sampleName + ".aligned.unsorted.duplicates_marked",
      #     metrics_filename     = "duplicate_metrics",
      #     tools                = tools,
      #     gatk_version         = gatk_version
      #}

      call MarkDuplicatesSpark
      {
         input:
           input_bams           = MergeBamAlignment.output_bam,
           output_bam_basename  = base_name + "." + sampleName + ".aligned.unsorted.duplicates_marked",
           metrics_filename     = "duplicate_metrics",
           tools                = tools,
           gatk_version         = gatk_version,
           memory               = "128 GB",
           cpu                  = 4,
           tmp_dir              = "/net/home/isaevt/tmp_dir"
      }

      call SortSam_byQuery as sortBeforeFix
      {
         input:
           input_bam           = MarkDuplicatesSpark.output_bam,
           output_bam_basename = base_name + "." + sampleName + ".aligned.duplicate_marked.sorted",
           compressionLvl      = 5,
           cpus                = 8,
           sambamba            = tools + "/my_sambamba/sambamba",
           flag                = ""
      }

      call FixTags as FixSampleBam
      {
         input:
           input_bam                 = sortBeforeFix.output_bam,
           input_bam_index           = sortBeforeFix.output_bam_index,
           output_bam_basename       = base_name + "." + sampleName + ".aligned.duplicate_marked.sorted.fixed",
           ref_fasta                 = ref_fasta,
           ref_dict                  = ref_dict,
           ref_fasta_index           = ref_fasta_index,
           ref_bwt                   = ref_bwt,
           ref_sa                    = ref_sa,
           ref_amb                   = ref_amb,
           ref_ann                   = ref_ann,
           ref_pac                   = ref_pac,
           fix_tags_java_heap_memory = "8000m",
           tools                     = tools,
           gatk_version              = gatk_version,
           memory                    = "8 GB",
           cpu                       = 1
      }

      scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping)
      {
         call BaseRecalibrator
         {
            input:
              input_bam                     = FixSampleBam.output_bam,
              input_bam_index               = FixSampleBam.output_bam_index,
              recalibration_report_filename = base_name + "." + sampleName + ".recal_data.csv",
              sequence_group_interval       = subgroup,
              dbSNP_vcf                     = dbSNP_vcf,
              dbSNP_vcf_idx                 = dbSNP_vcf_idx,
              known_indels_sites_VCF        = known_indels_sites_VCF,
              known_indels_sites_idx        = known_indels_sites_idx,
              ref_fasta                     = ref_fasta,
              ref_dict                      = ref_dict,
              ref_fasta_index               = ref_fasta_index,
              ref_bwt                       = ref_bwt,
              ref_sa                        = ref_sa,
              ref_amb                       = ref_amb,
              ref_ann                       = ref_ann,
              ref_pac                       = ref_pac,
              tools                         = tools,
              gatk_version                  = gatk_version
         }
      }

      call GatherBqsrReports
      {
         input:
           input_bqsr_reports     = BaseRecalibrator.recalibration_report,
           output_report_filename = base_name + "." + sampleName + ".recal_data.csv",
           tools                  = tools,
           gatk_version           = gatk_version
      }

      scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped)
      {
         call ApplyBQSR
         {
           input:
             input_bam               = FixSampleBam.output_bam,
             input_bam_index         = FixSampleBam.output_bam_index,
             output_bam_basename     = base_name + "." + sampleName + ".aligned.duplicates_marked.recalibrated",
             recalibration_report    = GatherBqsrReports.output_bqsr_report,
             sequence_group_interval = subgroup,
             tools                   = tools,
             gatk_version            = gatk_version,
             ref_fasta               = ref_fasta,
             ref_dict                = ref_dict,
             ref_fasta_index         = ref_fasta_index,
             ref_bwt                 = ref_bwt,
             ref_sa                  = ref_sa,
             ref_amb                 = ref_amb,
             ref_ann                 = ref_ann,
             ref_pac                 = ref_pac
         }
      }

      call GatherBamFiles
      {
         input:
           input_bams          = ApplyBQSR.recalibrated_bam,
           input_bams_inds     = ApplyBQSR.recalibrated_bam_index,
           output_bam_basename = base_name + "." + sampleName + ".recal.realign.dedup",
           cpu                 = 1,
           tools               = tools,
           gatk_version        = gatk_version
      }

      scatter (subInterval in scattered_calling_intervals)
      {
         call HaplotypeCaller
         {
            input:
              input_bam        = GatherBamFiles.output_bam,
              input_bam_index  = GatherBamFiles.output_bam_index,
              interval         = subInterval,
              interval_padding = 500,
              dbsnp            = dbSNP_vcf,
              dbsnp_index      = dbSNP_vcf_idx,
              gvcf_basename    = subInterval + ".g",
              ref_fasta        = ref_fasta,
              ref_dict         = ref_dict,
              ref_fasta_index  = ref_fasta_index,
              ref_bwt          = ref_bwt,
              ref_sa           = ref_sa,
              ref_amb          = ref_amb,
              ref_ann          = ref_ann,
              ref_pac          = ref_pac,
              memory           = "8GB",
              cpu              = 1,
              tools            = tools,
              gatk_version     = gatk_version
         }

         call MoveGvcf
         {
            input:
              gvcf     = HaplotypeCaller.output_gvcf,
              gvcf_idx = HaplotypeCaller.output_gvcf_index,
              dir      = gvcfs_dir
         }
      }

      #call GatherResults
      #{
      #   input:
      #     bam = GatherBamFiles.output_bam,
      #     bai = GatherBamFiles.output_bam_index,
      #     gvcfs = HaplotypeCaller.output_gvcf,
      #     gvcfs_idx = HaplotypeCaller.output_gvcf_index,
      #     dir_name = base_name + sampleName
      #}

      output
      {
         Array[File] gvcfs = HaplotypeCaller.output_gvcf
         Array[File] gvcfs_idx = HaplotypeCaller.output_gvcf_index
         Pair[String, Pair[File, File]] sample_bam = (sampleName, (GatherBamFiles.output_bam, GatherBamFiles.output_bam_index))
         Int done = 1
      }
}

task MoveGvcf
{
   File gvcf
   File gvcf_idx
   String dir

   command
   <<<
     cp ${gvcf} ${dir}
     cp ${gvcf_idx} ${dir}
   >>>

   runtime
   {
      cpu: 1
      memory: "100 MB"
   }
}

task HaplotypeCaller
{ 
  File input_bam
  File input_bam_index
 
  String interval
  File dbsnp
  File dbsnp_index
  Int interval_padding
  String gvcf_basename
  
  File ref_fasta
  File ref_dict
  File ref_fasta_index
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac
      
  String memory
  Int cpu
  
  String tools
  String gatk_version 

  command
  <<<
    ${tools}/gatk-${gatk_version}/gatk --java-options "-Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" HaplotypeCaller \
      -R ${ref_fasta} \
      -O ${gvcf_basename}vcf.gz \
      -I ${input_bam} \
      -L ${interval} \
      -D ${dbsnp} \
      -ip ${interval_padding} \
      -ERC GVCF \
      --max-alternate-alleles 6 \
      -RF OverclippedReadFilter \
      --smith-waterman FASTEST_AVAILABLE \
      --native-pair-hmm-use-double-precision true
  >>>

  runtime
  {
     memory: memory
     cpu: cpu
     docker:"broadinstitute/gatk:latest"
  }

  output
  {
     File output_gvcf = "${gvcf_basename}vcf.gz"
     File output_gvcf_index = "${gvcf_basename}vcf.gz.tbi"
  }
}

task GatherBamFiles
{
  Array[File] input_bams
  Array[File] input_bams_inds
  String output_bam_basename
  Int cpu
  String tools
  String gatk_version 

  command
  {
    ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms2000m -Xmx4000m" GatherBamFiles \
      --INPUT ${sep=' --INPUT ' input_bams} \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
  }

  runtime
  {
    memory: "4096 MB"
    cpu: cpu
    docker: "broadinstitute/gatk:latest"
  }

  output
  {
    File output_bam       = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
  }
}

task Fastq_to_uBAM
{
   File fastq1
   File fastq2

   String sampleName
   String output_uBAM_baseName

   String tools 
   String gatk_version

   Int cpu
   String memory
   String java_heap

   command
   {
     ${tools}/gatk-${gatk_version}/gatk --java-options "${java_heap}" FastqToSam \
     --FASTQ ${fastq1} \
     --FASTQ2 ${fastq2} \
     --OUTPUT ${output_uBAM_baseName}.bam \
     --READ_GROUP_NAME ${sampleName} \
     --SAMPLE_NAME ${sampleName} \
     --LIBRARY_NAME bgm_lib \
     --PLATFORM ILLUMINA
   }

   runtime
   {
      memory: memory
      cpu: cpu
   }

   output
   {
      File uBAM = "${output_uBAM_baseName}.bam" 
   }
}



# Get version of BWA
task GetBwaVersion
{
  String tools
  command 
  {
    # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed 
    # because the sed may also fail with that error and that is something we actually want to fail on.
    ${tools}/bwa/bwa-0.7.17/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime 
  {
    memory: "100 MB"
    cpu: 1
  }

  output 
  {
    String version = read_string(stdout())
  }
}


task BwaMem
{
  File fastq1
  File fastq2

  String bwa_commandline
  String output_bam_basename

  File ref_fasta
  File ref_dict
  File ref_fasta_index
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac
      
  String tools
  String memory
  Int cpu

  command
  {
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    bwa_threads=${cpu}
    
    ${tools}/bwa/bwa-0.7.17/${bwa_commandline} ${fastq1} ${fastq2} | \
    samtools view -1 -@ ${cpu} - > ${output_bam_basename}.bam 2>> bwa_err.log
  }

  runtime
  {
    memory: memory
    cpu: cpu
    docker: "timuris/bwa:latest"
  }

  output 
  {
    File output_bam     = "${output_bam_basename}.bam"
    File bwa_stderr_log = "bwa_err.log"
  }
}

task SortSam_byQuery
{
  File input_bam
  String output_bam_basename
  Int compressionLvl
  Int cpus
  String flag
  String sambamba

  command
  {
     ${sambamba} sort -m 8G -o ${output_bam_basename}.bam ${flag} -l ${compressionLvl} -t ${cpus} ${input_bam}
     ${sambamba} index ${output_bam_basename}.bam
  }

  runtime 
  {
    cpu: cpus
    memory: "8000 MB"
    docker: "timuris/bam_sort:latest"
  }

  output
  {
    File output_bam       = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bam.bai"
  }
}

task MergeBamAlignment
{
   File alignedBAM
   File uBAM
   File uBAM_index
   
   String output_bam_basename
   
   File ref_fasta
   File ref_dict
   File ref_fasta_index
   File ref_bwt
   File ref_sa
   File ref_amb
   File ref_ann
   File ref_pac
      
   String tmp_dir
 
   String bwa_version
   String bwa_mem_commandline   

   String tools
   String gatk_version

   command
   {   
     ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms3G" MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ATTRIBUTES_TO_REMOVE NM \
      --ATTRIBUTES_TO_REMOVE MD \
      --ALIGNED_BAM ${alignedBAM} \
      --UNMAPPED_BAM ${uBAM} \
      --OUTPUT ${output_bam_basename}.bam \
      --REFERENCE_SEQUENCE ${ref_fasta} \
      --PAIRED_RUN true \
      --SORT_ORDER "queryname" \
      --CREATE_INDEX=true \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 2000000 \
      --TMP_DIR ${tmp_dir} \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --PROGRAM_RECORD_ID "bwamem" \
      --PROGRAM_GROUP_VERSION "${bwa_version}" \
      --PROGRAM_GROUP_COMMAND_LINE "${bwa_mem_commandline}" \
      --PROGRAM_GROUP_NAME "bwamem" \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --UNMAP_CONTAMINANT_READS true \
      --ADD_PG_TAG_TO_READS false
   }

   runtime
   {
      cpu: 2
      memory: "16 GB"
      docker: "broadinstitute/gatk:latest"
   }

   output
   {
      File output_bam = "${output_bam_basename}.bam"
   }
  
}

task MarkDuplicatesSpark
{
   Array[File] input_bams
   #Array[File] input_bams_inds
   String output_bam_basename
   String metrics_filename
   String tools
   String gatk_version
   String memory
   String tmp_dir
   Int cpu

   command
   {
      ${tools}/gatk-${gatk_version}/gatk  MarkDuplicatesSpark \
      -I ${sep=' -I ' input_bams} \
      -O ${output_bam_basename}.bam \
      -M ${metrics_filename} \
      -VS SILENT \
      --optical-duplicate-pixel-distance 2500 \
      --spark-master local[${cpu}] \
      --tmp-dir ${tmp_dir}
   }

   runtime
   {
      cpu: cpu
      memory: memory
   }

   output
   {
      File output_bam        = "${output_bam_basename}.bam"
      File output_bam_index  = "${output_bam_basename}.bam.bai"
      File duplicate_metrics = "${metrics_filename}"
   }
}

task MarkDuplicates
{
   Array[String] input_bams
   String output_bam_basename
   String metrics_filename
   String tools
   String gatk_version

   command
   { 
      ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms16g" MarkDuplicates \
      --INPUT ${sep=' --INPUT ' input_bams} \
      --OUTPUT ${output_bam_basename}.bam \
      --METRICS_FILE ${metrics_filename} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER queryname \
      --CLEAR_DT false \
      --ADD_PG_TAG_TO_READS false \
      --CREATE_MD5_FILE true \
      --COMPRESSION_LEVEL 5 \
      --MAX_RECORDS_IN_RAM 1000000
   }

   runtime
   {
      docker: "broadinstitute/gatk:latest"
      cpu: 1
      memory: "64 GB"
   }

  output
  {
     File output_bam        = "${output_bam_basename}.bam"
     File duplicate_metrics = "${metrics_filename}"
  }
}

task CreateSequenceGroupingTSV 
{
  File ref_dict

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    # hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence

    hg38_protection_tag = "" #hg19 not need it?    

    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>

  runtime
  {
    memory: "3750 MB"
    cpu: 1
    docker: "timuris/bam_sort:latest"
  }

  output 
  {
    Array[Array[String]] sequence_grouping               = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

task BaseRecalibrator 
{
  File input_bam
  File input_bam_index
  
  String recalibration_report_filename
  Array[String] sequence_group_interval

  File dbSNP_vcf
  File dbSNP_vcf_idx

  File known_indels_sites_VCF
  File known_indels_sites_idx

  File ref_fasta
  File ref_dict
  File ref_fasta_index
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac
  
  String tools
  String gatk_version

  command 
  {
    ${tools}/gatk-${gatk_version}/gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -OQ \
      -O ${recalibration_report_filename} \
      -known-sites ${dbSNP_vcf} \
      -known-sites ${known_indels_sites_VCF} \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime 
  {
    memory: "8192 MB"
    cpu: 1
    docker: "broadinstitute/gatk:latest"
  }

  output 
  {
    File recalibration_report = "${recalibration_report_filename}"
  }
}

task FixTags
{
  File input_bam
  File input_bam_index
  
  String output_bam_basename
  
  File ref_fasta
  File ref_dict
  File ref_fasta_index
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac
      
  String tools
  String gatk_version
  String fix_tags_java_heap_memory
  String memory
  Int cpu
  
  command
  {
    ${tools}/gatk-${gatk_version}/gatk --java-options "-Xmx${fix_tags_java_heap_memory}" SetNmAndUqTags \
    --INPUT ${input_bam} \
    --OUTPUT ${output_bam_basename}.bam \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true \
    --REFERENCE_SEQUENCE ${ref_fasta}
  }

  runtime
  {
    docker: "broadinstitute/gatk:latest"
    memory: memory
    cpu: cpu
  }

  output
  {
    File output_bam       = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5   = "${output_bam_basename}.bam.md5"
  }
}


task GatherBqsrReports 
{
  Array[File] input_bqsr_reports
  String output_report_filename
  String tools
  String gatk_version

  command
  {
      ${tools}/gatk-${gatk_version}/gatk --java-options "-Xmx2000m" GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename} 
  }

  runtime 
  {
    memory: "2048 MB"
    cpu: 1
    docker: "broadinstitute/gatk:latest"
  }

  output
  {
    File output_bqsr_report = "${output_report_filename}"
  }
}

task ApplyBQSR
{
  File input_bam
  File input_bam_index
  
  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  
  File ref_fasta
  File ref_dict
  File ref_fasta_index
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac
  
  String tools
  String gatk_version

  command 
  {
    ${tools}/gatk-${gatk_version}/gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
      ApplyBQSR \
      -OBM \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -OQ \
      -O ${output_bam_basename}.bam \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      -L ${sep=" -L " sequence_group_interval}
  }

  runtime
  {
    memory: "4096 MB"
    cpu: 1
    docker: "broadinstitute/gatk:latest"
  }

  output
  {
    File recalibrated_bam          = "${output_bam_basename}.bam"
    File recalibrated_bam_index    = "${output_bam_basename}.bai"
    File recalibrated_bam_checksum = "${output_bam_basename}.bam.md5"
  }
}

task SplitFQ_parallel
{
  File fq
  String case_name
  String suffix

  Int cpu
  String memory

  command
  <<<
     set -e
     zcat ${fq} | /net/home/isaevt/bin/bin/split --additional-suffix=.fq${suffix} -l 60000000 - ${case_name}
     parallel -j${cpu} gzip ::: *.fq${suffix}
  >>>

  runtime
  {
     cpu: cpu
     memory: memory
  }
  output
  {
     Array[File] fqs = glob("*.fq${suffix}.gz")
  }
}

task GatherFiles
{
   Array[Array[File]] files
   String suffix
	
   command
   <<<
     find ../ -type f -name "*${suffix}.gz" -exec mv {} "$(pwd)" \;
   >>>

   runtime
   {
      cpu: 1
      memory: "1024 MB"
      docker: "ubuntu:latest"
   }

   output
   {
      Array[File] o1 = glob("*${suffix}.gz")
   }
}
