workflow downstream
{
   Array[Int] finished
   String ref_fasta = "/net/bgm/resources/hg19.fasta"
   String gvcfs_dir = "/net/bgm/resources/wgs/" #all g.vcfs live here
   Array[String] types_to_include = ["SNP", "INDEL"]

   String tools     = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/tools"
   #String db_path   = "/net/home/isaevt/db_test_2/all_dbs/"

   File chr_intervals_list     = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/intervals_new.txt"
   File chr_list = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/chr_list.txt"
   Array[String] chrs = read_lines(chr_list)

   #File chr_list    = "/net/bgm/resources/variant_caller_intervals.file"
   #Array[String] chr_intervals = read_lines(chr_intervals_list)

   File dbSNP_vcf = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/dbsnp_147.hg19.vcf"
   File dbSNP_vcf_idx = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/dbsnp_147.hg19.vcf.idx"

   #Array[String] chromosomes = ["chrM", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
   Array[String] chromosomes = ["chrM", "chr21", "chr22", "chrX", "chrY"]


   String gvcfsFolder = "/net/bgm/resources/wes/mem"
   #File intervals_file = "/net/home/isaevt/python/intervals/new.intervals.txt"
   #new.intervals.chr2122.txt

   File intervals_file = "/net/home/isaevt/python/intervals/count.txt"
   #File intervals_file = "/net/bgm/resources/variant_caller_intervals.file"

   Array[String] chr_intervals = read_lines(intervals_file)

   File fam = "/net/bgm/cases/bgm0102_wes/bgm0102.fam"

   #Array[String] recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"]
   Array[String] recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0"]
   Array[String] recalibration_annotation_values = ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]

   #Array[String] SNPsVariantRecalibrator_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0"]
   Array[String] SNPsVariantRecalibrator_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.0"]
   Array[String] SNPsVariantRecalibrator_recalibration_annotation_values = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]
      
   File dbsnp_resource_vcf = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/dbsnp_147.hg19.vcf"
   File dbsnp_resource_vcf_index = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/dbsnp_147.hg19.vcf.idx"

   File mills_resource_vcf = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
   File mills_resource_vcf_index = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx"

   File hapmap_resource_vcf = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/resources/hapmap_3.3.hg19.sites.vcf"
   File hapmap_resource_vcf_index = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/resources/hapmap_3.3.hg19.sites.vcf.idx"
   File omni_resource_vcf = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/resources/1000G_omni2.5.hg19.sites.vcf"
   File omni_resource_vcf_index = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/resources/1000G_omni2.5.hg19.sites.vcf.idx"
   File one_thousand_genomes_resource_vcf = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/resources/1000G_phase1.indels.hg19.vcf"
   File one_thousand_genomes_resource_vcf_index = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/resources/1000G_phase1.indels.hg19.vcf.idx"

   File case_ids_for_test = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/case_ids_for_test.txt"
   File case_bams_for_test = "/net/home/isaevt/cromwell/wgs_workflow_gatk4/wes/case_bams_for_test.txt"

   Array[String] case_sample_names = ["bgm0102a1", "bgm0102u1", "bgm0102u2"]

  # scatter (chromosome in chromosomes)
  # {
  #    call GetChrGVCFs
  #    {
  #       input:
  #         workflow_id = "3",
  #         chr = chromosome + ".",
  #         folder_to_search = gvcfsFolder,
  #         tools = tools
  #    }
  #
  #    call GetChrIntervals
  #    {
  #       input:
  #         workflow_id = "1",
  #         chr = chromosome + ":",
  #         intervals_file = intervals_file
  #    }
  #
  #    call GenomicsDB_Import.importGVCFs
  #    {
  #       input:
  #         gvcfs = GetChrGVCFs.gvcfs,
  #         intervals = GetChrIntervals.intervals,
  #         tools = tools,
  #         db_path = "/net/home/isaevt/db_test_final/",
  #         workflow_id = "3",
  #         ref_fasta = ref_fasta,
  #         dbSNP_vcf = dbSNP_vcf,
  #         dbSNP_vcf_idx = dbSNP_vcf_idx
  #    }
  # }

  # Array[Array[File]] genotyped_vcfs = importGVCFs.vcfs

  # call CollectFilenames
  # {
  #    input:
  #      filename_arrays = genotyped_vcfs,
  #      python_flatten  = "/net/home/isaevt/python/flatten_arrays.py"
  # }


   #Array[String] res = importGVCFs.finished
   String db_path = "/net/home/isaevt/genomics_db/"

   scatter (index in range(length(chr_intervals)))
   {
      Int i = index
      String chr = chr_intervals[i]
      String sub_interval = sub(chr, ":", "_")

      call GetChrGVCFs
      {
         input:
           workflow_id = "3",
           chr = sub(chr, ":.+", "") + ".",
           folder_to_search = gvcfsFolder,
           tools = tools
      }

      call DB_Import
      {
         input:
           workflow_id = "4",
           gatk_launch = tools + "/gatk-4.0.10.1/gatk",
           gvcfs = GetChrGVCFs.gvcfs,
           interval = chr,
           ref      = ref_fasta,
           #db_path = db_path + sub(sub_interval, "-", "_")
           db_path = db_path + i
      }
 
      #call ConsolidateGVCFs
      #{
      #   input:
      #     gatk_launch = tools + "/gatk-4.0.9.0/gatk",
      #     gvcfs = GetChrGVCFs.gvcfs,
      #     interval = chr,
      #     workflow_id = 6,
      #     ref = ref_fasta,
      #     dbSNP_vcf = dbSNP_vcf,
      #     out_basename = sub(chr, ":.+", "") + "_" + i
      #}
 
      #call GetChrIntervals
      #{
      #   input:
      #     workflow_id = "1",
      #     chr = sub(chr, ":.+", "") + ":",
      #     intervals_file = intervals_file
      #}

      call GenotypeGVCFs
      {
         input:
           db_path = DB_Import.path,
           ref_fasta = ref_fasta,
           gatk_launch = tools + "/gatk-4.0.10.0/gatk",
        
           chr = chr,
           dbSNP_vcf = dbSNP_vcf,
           dbSNP_vcf_idx = dbSNP_vcf_idx,
           out_base_name = sub(chr, ":.+", "") + "_" + i
      }

      #call GenotypeGVCFs_no_db
      #{
      #   input:
      #     ref_fasta = ref_fasta,
      #     gatk_launch = tools + "/gatk-4.0.9.0/gatk",
      #     gvcf = ConsolidateGVCFs.gvcf_out,
      #     interval = chr,
      #     dbSNP_vcf = dbSNP_vcf,
      #     dbSNP_vcf_idx = dbSNP_vcf_idx,
      #     out_base_name = sub(chr, ":.+", "") + "_" + i
      #}
   }

  #Array[File] vcfs = GenotypeGVCFs.output_vcf
  #Array[File] vcfs_idx = GenotypeGVCFs.output_vcf_idx

  #call GatherVcfs
  #{
  #   input:
  #     vcfs = GenotypeGVCFs.output_vcf,
  #     vcfs_idx = GenotypeGVCFs.output_vcf_idx,
  #     picard = tools + "/picard.jar",
  #     out_base_name = "combined"
  #}

  call BcfToolsMergeVcf as MergeGenotypeRes
  {
     input:
       vcfs = GenotypeGVCFs.output_vcf,
       vcfs_ind = GenotypeGVCFs.output_vcf_idx,
       tools = tools,
       out_base_name = "genotyped",
       out_extension = "vcf.gz",
       out_type = "z",
       num_threads = 4
  }
 
  
  #call SelectVariantsByType as SelectIndels
  #{
  #   input:
  #     type_to_select = "INDEL",
  #     ref_fasta = ref_fasta,
  #     gatk_launch = tools + "/gatk-4.0.5.1/gatk",
  #     vcf = GatherVcfs.combined_vcf,
  #     output_base_name = "selectedVariants_indels" 
  #}

  #call SelectVariantsByType as SelectSnps
  #{
  #   input:
  #     type_to_select = "SNP",
  #     ref_fasta = ref_fasta,
  #     gatk_launch = tools + "/gatk-4.0.5.1/gatk",
  #     vcf = GatherVcfs.combined_vcf,
  #     output_base_name = "selectedVariants_snps"
  #}

  call IndelsVariantRecalibrator
  {
     input:
       recalibration_filename = "indels.recalibration.file",
       tranches_filename = "indels.tranches.file",
       recalibration_tranche_values = recalibration_tranche_values,
       recalibration_annotation_values = recalibration_annotation_values,
       vcf = MergeGenotypeRes.out_vcf,
       vcf_index = MergeGenotypeRes.out_vcf_index,
       dbsnp_resource_vcf = dbsnp_resource_vcf,
       dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
       mills_resource_vcf = mills_resource_vcf,
       mills_resource_vcf_index = mills_resource_vcf_index,
       gatk_launch = tools + "/gatk-4.0.10.0/gatk"
  }

  call SNPsVariantRecalibrator
  {
     input:
       recalibration_filename = "indels.recalibration.file",
       tranches_filename = "indels.tranches.file",
       recalibration_tranche_values = SNPsVariantRecalibrator_recalibration_tranche_values,
       recalibration_annotation_values = SNPsVariantRecalibrator_recalibration_annotation_values,
       vcf = MergeGenotypeRes.out_vcf,
       vcf_index = MergeGenotypeRes.out_vcf_index,
       dbsnp_resource_vcf = dbsnp_resource_vcf,
       dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
       hapmap_resource_vcf = hapmap_resource_vcf,
       hapmap_resource_vcf_index = hapmap_resource_vcf_index,
       omni_resource_vcf = omni_resource_vcf,
       omni_resource_vcf_index = omni_resource_vcf_index,
       one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
       one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
       gatk_launch = tools + "/gatk-4.0.10.0/gatk"
  }

  call ApplyRecalibration
  {
     input:
       recalibrated_vcf_filename = "recalibrated.vcf",
       input_vcf =  MergeGenotypeRes.out_vcf,
       input_vcf_index = MergeGenotypeRes.out_vcf_index,

       indels_recalibration = IndelsVariantRecalibrator.recalibration,
       indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
       indels_tranches = IndelsVariantRecalibrator.tranches,

       snps_recalibration = SNPsVariantRecalibrator.recalibration,
       snps_recalibration_index = SNPsVariantRecalibrator.recalibration_index,
       snps_tranches = SNPsVariantRecalibrator.tranches,

       indel_filter_level = 95.0,
       snp_filter_level   = 99.0,

       gatk_path = tools + "/gatk-4.0.10.0/gatk"  
  }

  call VepAnnotate
  {
     input:
       number_of_threads = 32,
       vep_path = "/net/bgm/tools/vep-85/variant_effect_predictor.pl",
       db_path  = "/net/bgm/resources/vep-85-cache",
       plugin_path = "/net/bgm/tools/vep-85/plugins",
       exac_path = "/net/bgm/resources/ExAC.r0.3.sites.vep.b37.vcf.gz",
       input_vcf = ApplyRecalibration.recalibrated_vcf,
       input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
       output_basename = "final.vep.annotated"
  }

  call FindAllBams
  {
     input:
       dir_to_search  = "/net/bgm/cases/",
       platform       = "wes",
       tools          = tools
  }

  scatter (chr in chrs)
  {
     call SplitVcfByChr
     {
        input:
          in_vcf = VepAnnotate.recalibrated_vcf,
          in_vcf_index = VepAnnotate.recalibrated_vcf_index,
          chr = chr
     }

     call BgmBayesDeNovo
     {
        input:
          case_ids  = case_ids_for_test,
          input_vcf = SplitVcfByChr.vcf,
          unrelated_bams = FindAllBams.all_bams,
          case_bams = case_bams_for_test, 
          x = 1.0,
          p = 0.005,
          e = 0.05,
          tools = tools
     }

     call ProcessBayesDeNovoSt1Res
     {
        input:
          res = BgmBayesDeNovo.stage_one_out,
          tools = tools
     }

     call BgmCompoundHetCaller
     {
        input:
          case_ids = case_ids_for_test,
          input_vcf = SplitVcfByChr.vcf,
          t = 0.9,
          e = 0.1,
          a = 0.2,
          tools = tools
     }

     call ProcessBgmCompHetRes
     {
        input:
          res = BgmCompoundHetCaller.out,
          tools = tools,
          vcf = SplitVcfByChr.vcf 
     }

     call BgmHomRec
     {
        input:
          case_ids = case_ids_for_test,
          input_vcf = SplitVcfByChr.vcf,
          x = 1.0,
          p = 0.9,
          e = 0.05,
          a = 0.2,
          tools = tools
     }

     call ProcessHomRecRes
     {
        input:
          res = BgmHomRec.out,
          tools = tools
     }

     call ABCaller as AutosomalDominantCaller
     {
        input:
          in_vcf = SplitVcfByChr.vcf,
          fam    = fam,
          id     = 0,
          py_script = tools + "/ab_caller.py",
          output_name = "autosomalDominant.calls" 
     }

     call ABCaller as HomozygousRecessiveCaller
     {
        input:
          in_vcf = SplitVcfByChr.vcf,
          fam    = fam,
          id     = 1,
          py_script = tools + "/ab_caller.py",
          output_name = "homozygousRecessive.calls"
     }

     call ABCaller as DeNovoCaller
     {
        input:
          in_vcf = SplitVcfByChr.vcf,
          fam    = fam,
          id     = 2,
          py_script = tools + "/ab_caller.py",
          output_name = "deNovoCaller.calls"
     }

     call ABCaller as CompoundHeterozygousCaller
     {
        input:
          in_vcf = SplitVcfByChr.vcf,
          fam    = fam,
          id     = 3,
          py_script = tools + "/ab_caller.py",
          output_name = "compoundHeterozygous.calls"
     }

     call GatherCalls
     {
        input:
          tools = tools,
          bgmHomRec_res = ProcessHomRecRes.hom_rec_calls,
          compoundHeterozygous_res = CompoundHeterozygousCaller.calls,
          deNovo_res = DeNovoCaller.calls,
          homozygousRecessive_res = HomozygousRecessiveCaller.calls,
          autosomalDominant_res = AutosomalDominantCaller.calls,
          bgmCompoundHet_res = ProcessBgmCompHetRes.comp_het_calls,
          bgmDeNovo_res = ProcessBayesDeNovoSt1Res.bayes_deNovo_calls_st1,
          vep_vcf = VepAnnotate.recalibrated_vcf,
          workflow_id = "1"
     }

     call Bgzip
     {
        input:
          tools = tools,
          vcf = GatherCalls.vcf_out,
          out_base_name = "calls." + chr
     }
  }

  #call GatherVcfs as GatherCallsVcf
  #{
  #   input:
  #     vcfs = GatherCalls.vcf_out,
  #     picard = tools + "/picard.jar",
  #     out_base_name = "calls"
  #}

  call BcfToolsMergeVcf
  {
     input:
       vcfs = Bgzip.vcf_out,
       vcfs_ind = Bgzip.vcf_out_index,
       tools = tools,
       out_base_name = "calls",
       out_extension = "vcf",
       out_type = "v",
       num_threads = 4
  }

  call MergeVcfCalls
  {
     input:
       calls = BcfToolsMergeVcf.out_vcf,
       vcf = VepAnnotate.recalibrated_vcf,
       vcf_index = VepAnnotate.recalibrated_vcf_index,
       tools = tools 
  }

  call SelectCaseVariants
  {
     input:
       tools = tools,
       gatk_version = "4.0.8.1",
       ref = ref_fasta,
       in_vcf = MergeVcfCalls.vcf_out,
       samples = case_sample_names,
       caller_names = ["BGM_DE_NOVO", "BGM_CMPD_HET", "BGM_HOM_REC", "BGM_AUTO_DOM", "BGM_BAYES_HOM_REC", "BGM_BAYES_DE_NOVO", "BGM_BAYES_CMPD_HET"],
       out_base_name = "xbrowse.vep"
  }


  output
  {
     File vcf = VepAnnotate.recalibrated_vcf
     File vcf_idx = VepAnnotate.recalibrated_vcf_index
  }
}

#task DB_Import
#{
#  String gatk_launch
#  Array[String] gvcfs
#  String interval
#  String db_path
#  String workflow_id
#
#  command
#  {
#     ${gatk_launch} --java-options "-Xms4g -XX:ParallelGCThreads=1" GenomicsDBImport \
#       --batch-size 50 \
#       -ip 250 \
#       --reader-threads 5 \
#       --overwrite-existing-genomicsdb-workspace true \
#       -L ${interval} \
#       --genomicsdb-workspace-path ${db_path} \
#       -V ${sep=' -V ' gvcfs}
#
#     echo ${db_path}
#     echo ${workflow_id}
#  }
#
#  runtime
#  {
#     cpu: 2
#     memory: "6 GB"
#  }
#
#  output
#  {
#     String path = "${db_path}"
#  }
#
#}

task ConsolidateGVCFs
{
   String gatk_launch
   Array[String] gvcfs
   String interval
   String workflow_id
   String ref
   String out_basename
   String dbSNP_vcf

   command
   {
      ${gatk_launch} --java-options "-Xms1g -XX:ParallelGCThreads=1" CombineGVCFs \
        -R ${ref} \
        -V ${sep=' -V ' gvcfs} \
        -L ${interval} \
        -D ${dbSNP_vcf} \
        -ip 200 \
        -O ${out_basename}.g.vcf.gz

      #${gatk_launch} --java-options "-Xms1g -XX:ParallelGCThreads=1" CombineGVCFs \
      echo ${workflow_id}
   }

   runtime
   {
      cpu: 1
      memory: "1500 MB"
   }

   output
   {
      File gvcf_out = "${out_basename}.g.vcf.gz"
      File gvcf_idx_out = "${out_basename}.g.vcf.gz.tbi"
   }
}

task DB_Import
{
  String gatk_launch
  Array[String] gvcfs
  String interval
  String db_path
  String workflow_id
  String ref

  command
  {
     ${gatk_launch} --java-options "-Xms4g -XX:ParallelGCThreads=1" GenomicsDBImport \
       --batch-size 50 \
       -ip 100 \
       -R ${ref} \
       --reader-threads 5 \
       --overwrite-existing-genomicsdb-workspace true \
       -L ${interval} \
       --genomicsdb-workspace-path ${db_path} \
       --tmp-dir /net/home/isaevt/tmp_db_import \
       -V ${sep=' -V ' gvcfs}

     echo ${db_path}
     echo ${workflow_id}
  }

  runtime
  {
     cpu: 2
     memory: "4 GB"
  }

  output
  {
     String path = "${db_path}"
  }

}



task SelectCaseVariants
{
   String tools
   String gatk_version
   String ref
   File in_vcf
   Array[String] samples
   Array[String] caller_names
   String out_base_name

   command
   {
      ${tools}/gatk-${gatk_version}/gatk --java-options "-XX:ParallelGCThreads=1 -Xms1000m" \
        SelectVariants \
        -R ${ref} \
        -sn ${sep=" -sn " samples} \
        -select "(!vc.isMonomorphicInSamples()) || ${sep=" || " caller_names}" \
        -V ${in_vcf} \
        -O ${out_base_name}.vcf  
   }

   runtime
   {
      memory: "1500 MB"
      cpu: 1
   }

   output
   {
      File out_vcf = "${out_base_name}.vcf"
   }
}

task CollectFilenames
{
    Array[Array[String]] filename_arrays
    File python_flatten

    command
    {
        echo "${sep='\n' filename_arrays}" > raw_array
        python ${python_flatten} < raw_array > file_of_filenames
    }

    output
    {
        Array[String] vcfs = read_lines("./file_of_filenames")
    }

    runtime
    {
       memory: "100MB"
       cpu: 1
    }
}



task GatherCalls
{
   File bgmHomRec_res
   File compoundHeterozygous_res
   File deNovo_res 
   File homozygousRecessive_res
   File autosomalDominant_res
   File bgmCompoundHet_res
   File bgmDeNovo_res
   File vep_vcf

   String tools

   String workflow_id

   command
   {
      echo ${workflow_id}
      python ${tools}/callsGathering.py ${deNovo_res} ${bgmHomRec_res} ${compoundHeterozygous_res} ${autosomalDominant_res} ${bgmDeNovo_res} ${homozygousRecessive_res} ${autosomalDominant_res} ${vep_vcf}
   }

   runtime
   {
      memory: "500 MB"
      cpu: 1
   }

   output
   {
      File vcf_out = "calls.vcf"
   } 
}

task BcfToolsMergeVcf
{
   Array[String] vcfs
   Array[String] vcfs_ind
   String tools
   String out_base_name
   String out_extension
   String out_type
   Int num_threads

   command
   {
      ${tools}/bcftools-1.9/bcftools concat --threads ${num_threads} -a -D -O ${out_type} -o ${out_base_name}.${out_extension} ${sep=' ' vcfs}
      if [ "${out_type}" = "z" ]; then
         tabix -p vcf ${out_base_name}.${out_extension}
      fi 
   }

   runtime
   {
      memory: "1 GB"
      cpu: num_threads
   }

   output
   {
      File out_vcf = "${out_base_name}.${out_extension}"
      File? out_vcf_index = "${out_base_name}.vcf.gz.tbi"
   }
}

task Bgzip
{
  File vcf
  String tools
  String out_base_name  

  command
  {
     ${tools}/htslib-1.9/bgzip -@ 4 -c ${vcf} > ${out_base_name}.vcf.gz 
     ${tools}/htslib-1.9/tabix -p vcf ${out_base_name}.vcf.gz
  }

  runtime
  {
     memory: "1 GB"
     cpu: 4
  }

  output
  {
     File vcf_out = "${out_base_name}.vcf.gz"
     File vcf_out_index = "${out_base_name}.vcf.gz.tbi"
  }
}

task MergeVcfCalls
{
   File calls
   File vcf
   File vcf_index
   String tools

   command
   {
      ${tools}/htslib-1.9/bgzip -@ 4 ${calls}
      ${tools}/htslib-1.9/tabix -p vcf ${calls}.gz
      ${tools}/bcftools-1.9/bcftools annotate -a ${calls}.gz -c INFO ${vcf} > original.with.calls.vcf
   }

   runtime
   {
      memory: "1 GB"
      cpu: 4
   }

   output
   {
      File vcf_out = "original.with.calls.vcf"
   }
} 

task ABCaller
{
   File in_vcf
   File fam
   File py_script
   Int id
   String output_name

   command
   {
      python ${py_script} ${in_vcf} ${fam} ${id} > ${output_name}
   }

   runtime
   {
      memory: "1 GB"
      cpu: 1
   }

   output
   {
      File calls = "${output_name}"
   }
}

task ProcessHomRecRes
{
   File res
   String tools

   command
   <<<
     python ${tools}/processHomRec.py ${res} > hom_rec_calls.txt
   >>>

   runtime
   {
      memory: "500 MB"
      cpu: 1
   }

   output
   {
      File hom_rec_calls = "hom_rec_calls.txt"
   }
}

task ProcessBgmCompHetRes
{
   File res
   File vcf
   String tools

   command
   <<<
     python ${tools}/processBgmCompHet.py ${res} ${vcf} > bgm_comp_het_calls.txt
   >>>

   runtime
   {
      memory: "500 MB"
      cpu: 1
   }

   output
   {
      File comp_het_calls = "bgm_comp_het_calls.txt"
   }
}

task ProcessBayesDeNovoSt1Res
{
   File res
   String tools

   command
   <<<
     python ${tools}/processBayesDeNovo.py ${res} > bayes_de_novo1_calls.txt 
   >>>

   runtime
   {
     memory: "500 MB"
     cpu: 1
   }

   output
   {
     File bayes_deNovo_calls_st1 = "bayes_de_novo1_calls.txt"
   }
}

task BgmBayesDeNovo
{
  File case_ids
  File input_vcf

  File unrelated_bams
  File case_bams

  Float x
  Float p
  Float e  

  String tools

  command
  {
     ${tools}/denovo.out -I ${input_vcf} -O bayes_de_novo_st1.out.vcf -T ${case_ids} -X ${x} -P ${p} -E ${e}
     #python ${tools}/pysam_prac.py -I bayes_de_novo_st1.out.vcf -U ${unrelated_bams} -T ${case_bams} -O bayes_de_novo_st2.out
  }

  runtime
  {
     memory: "1 GB"
     cpu: 1
  }

  output
  {
     File stage_one_out = "bayes_de_novo_st1.out.vcf"
     #File stage_two_out = "bayes_de_novo_st2.out"
  }
}

task BgmCompoundHetCaller
{
  File case_ids
  File input_vcf
  Float t
  Float e
  Float a
  
  String tools

  command
  {
     ${tools}/compound_het.out -I ${input_vcf} -O bgm_cmp_het.out -TR ${case_ids} -T ${t} -E ${e} -A ${a}
  }

  runtime
  {
     memory: "1 GB"
     cpu: 1
  }

  output
  {
     File out = "bgm_cmp_het.out"
  }
}

task BgmHomRec
{
   File case_ids
   File input_vcf
   Float x
   Float p
   Float e
   Float a  

   String tools

   command
   {
      ${tools}/homrec.out -I ${input_vcf} -O bgm_hom_rec.out -T ${case_ids} -X ${x} -P ${p} -E ${e} -A ${a}
   }

   runtime
   {
      memory: "1 GB"
      cpu: 1 
   }

   output
   {
      File out = "bgm_hom_rec.out"
   }
}

task GetChrIntervals
{
   String chr
   String workflow_id
   File intervals_file

   command
   <<<
      echo ${workflow_id} > out.txt
      awk '/${chr}/ {print $1}' ${intervals_file}
   >>>

   output
   {
      Array[String] intervals = read_lines(stdout())
   } 
}

task GetChrGVCFs
{
   String workflow_id
   String tools
   String chr
   String folder_to_search

   command
   {
      echo ${workflow_id} > out.txt
      python ${tools}/get_gvcfs.py ${chr} ${folder_to_search}
   }

   output
   {
      Array[String] gvcfs = read_lines(stdout())
   }
}

task FindAllBams
{
   String dir_to_search
   String platform
   String tools
      
   command
   {
      python ${tools}/findBams.py ${dir_to_search} ${platform}
   }

   runtime
   {
      memory: "100 MB"
      cpu: 1
   }

   output
   {
      File all_bams = "allBams.txt"
   }
}

task SplitVcfByChr
{
   File in_vcf
   File in_vcf_index
   String chr

   command
   {
      tabix -h ${in_vcf} ${chr} > ${chr}.vcf 
   }

   runtime
   {
      memory: "1 GB"
      cpu: 1
   }

   output
   {
      File vcf = "${chr}.vcf"
   }
}

task VepAnnotate
{
   String vep_path
   String db_path
   String plugin_path
   String exac_path
   File input_vcf
   File input_vcf_index
   String output_basename
   Int number_of_threads

   command
   {
      perl ${vep_path} --format vcf --force_overwrite --dir_cache ${db_path} --dir_plugin ${plugin_path} --offline --no_stats --vcf --everything \
           --plugin ExAC,${exac_path} --allele_number --input_file ${input_vcf} --output_file ${output_basename}.vcf.gz --tabix --fork ${number_of_threads} --buffer_size 250000
   }

   runtime
   {
     memory: "1 GB"
     cpu: number_of_threads
   }

   output
   {
     File recalibrated_vcf = "${output_basename}.vcf.gz"
     File recalibrated_vcf_index = "${output_basename}.vcf.gz.tbi"
   }

}

task ApplyRecalibration
{
  String recalibrated_vcf_filename
  File input_vcf
  File input_vcf_index
  File indels_recalibration
  File indels_recalibration_index
  File indels_tranches
  File snps_recalibration
  File snps_recalibration_index
  File snps_tranches

  Float indel_filter_level
  Float snp_filter_level

  String gatk_path

  command 
  {
    set -e
    ${gatk_path} --java-options "-Xms5g" \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ${input_vcf} \
      --recal-file ${indels_recalibration} \
      --tranches-file ${indels_tranches} \
      --truth-sensitivity-filter-level ${indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    ${gatk_path} --java-options "-Xmx5g" \
      ApplyVQSR \
      -O ${recalibrated_vcf_filename} \
      -V tmp.indel.recalibrated.vcf \
      --recal-file ${snps_recalibration} \
      --tranches-file ${snps_tranches} \
      --truth-sensitivity-filter-level ${snp_filter_level} \
      --create-output-variant-index true \
      -mode SNP
  }

  runtime 
  {
    memory: "5 GB"
    cpu: "1"
  }

  output 
  {
    File recalibrated_vcf = "${recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "${recalibrated_vcf_filename}.idx"
  }
}

task SNPsVariantRecalibrator
{
  String recalibration_filename
  String tranches_filename

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File vcf
  File vcf_index

  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index

  String gatk_launch

  command
  {
    ${gatk_launch} --java-options "-Xms3g" \
      VariantRecalibrator \
      -V ${vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode SNP \
      --max-gaussians 6 \
      -resource hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
      -resource omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
      -resource 1000G,known=false,training=true,truth=false,prior=10:${one_thousand_genomes_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_resource_vcf}
  }

  runtime
  {
    memory: "3 GB"
    cpu: "2"
  }

  output
  {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
  }
}


task IndelsVariantRecalibrator 
{
   String recalibration_filename
   String tranches_filename

   Array[String] recalibration_tranche_values
   Array[String] recalibration_annotation_values

   File vcf
   File vcf_index

   File mills_resource_vcf
   File dbsnp_resource_vcf

   File mills_resource_vcf_index
   File dbsnp_resource_vcf_index

   String gatk_launch

   command
   {
      ${gatk_launch} --java-options "-Xmx20g -Xms5g" \
                     VariantRecalibrator \
                     -V ${vcf} \
                     -O ${recalibration_filename} \
                     --tranches-file ${tranches_filename} \
                     --trust-all-polymorphic \
                     -tranche ${sep=' -tranche ' recalibration_tranche_values} \
                     -an ${sep=' -an ' recalibration_annotation_values} \
                     -mode INDEL \
                     --max-gaussians 4 \
                     -resource mills,known=false,training=true,truth=true,prior=12:${mills_resource_vcf} \
                     -resource dbsnp,known=true,training=false,truth=false,prior=2:${dbsnp_resource_vcf}
   }

   output
   {
      File recalibration = "${recalibration_filename}"
      File recalibration_index = "${recalibration_filename}.idx"
      File tranches = "${tranches_filename}"
   }

   runtime 
   {
      memory: "5 GB"
      cpu: "2"
   }
}

task GatherVcfs
{
   Array[String] vcfs
   Array[String] vcfs_idx
   String picard
   String out_base_name
  
   command
   {
      java -Xms2g -XX:ParallelGCThreads=1 -jar ${picard} GatherVcfs INPUT=${sep=' INPUT=' vcfs} OUTPUT=${out_base_name}.vcf.gz
   }

   output
   {
      File combined_vcf = "${out_base_name}.vcf.gz"
      File combined_vcf_index = "${out_base_name}.vcf.gz.tbi"
   }

   runtime
   {
      memory: "2 GB"
      cpu: 1
   }
}

task SelectVariantsByType
{
   String type_to_select
   String gatk_launch
   String ref_fasta
   File vcf
   String output_base_name

   command
   {
      ${gatk_launch} --java-options "-Xms2g -Xmx4g -XX:ParallelGCThreads=1" \
                     SelectVariants -R ${ref_fasta} \
                                    --select-type-to-include ${type_to_select} \
                                    -V ${vcf} \
                                    -O ${output_base_name}.vcf.gz 
   }

   output
   {
      File out_vcf = "${output_base_name}.vcf.gz"
   }

   runtime
   {
      memory: "2 GB"
      cpu: 1
   }
}  

task GenotypeGVCFs_no_db
{
   String gvcf
   String ref_fasta
   String gatk_launch

   String interval
   String out_base_name

   String dbSNP_vcf
   String dbSNP_vcf_idx

   command
   {
      ${gatk_launch} --java-options "-Xms2g -XX:ParallelGCThreads=2" GenotypeGVCFs \
                     -R ${ref_fasta} \
                     -V ${gvcf} \
                     -D ${dbSNP_vcf} \
                     -O ${out_base_name}.raw.vcf \
                     -L ${interval} \
                     -ip 250 \
                     --tmp-dir=/net/home/isaevt/tmp_for_gatk4
   }

   output
   {
      File output_vcf = "${out_base_name}.raw.vcf"
      File output_vcf_idx = "${out_base_name}.raw.vcf.idx"
   }

   runtime
   {
      memory: "1 GB"
      cpu: 2
   }
}


task GenotypeGVCFs
{
   String db_path
   String ref_fasta
   String gatk_launch

   String chr
   String out_base_name   

   String dbSNP_vcf
   String dbSNP_vcf_idx

   command
   {
      echo ${db_path}
      ${gatk_launch} --java-options "-Xms4g -XX:ParallelGCThreads=2" GenotypeGVCFs \
                     -R ${ref_fasta} \
                     -V gendb://${db_path} \
                     -D ${dbSNP_vcf} \
                     -O ${out_base_name}.raw.vcf.gz \
                     -L ${chr} \
                     -new-qual true \
                     --tmp-dir /net/home/isaevt/tmp_for_gatk4 \
   }

   output
   {
      File output_vcf = "${out_base_name}.raw.vcf.gz"
      File output_vcf_idx = "${out_base_name}.raw.vcf.gz.tbi"
   }

   runtime
   {
      memory: "8 GB"
      cpu: 1
   }
}

