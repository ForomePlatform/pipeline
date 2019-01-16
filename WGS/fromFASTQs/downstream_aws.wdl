workflow wes_downstream
{
   File ref_fasta = "s3://bgm-cromwell-results/resources/hg19.fasta"
   File ref_dict = "s3://bgm-cromwell-results/resources/hg19.dict"
   File ref_fasta_index = "s3://bgm-cromwell-results/resources/hg19.fasta.fai"
   File ref_bwt = "s3://bgm-cromwell-results/resources/hg19.fasta.bwt"
   File ref_sa = "s3://bgm-cromwell-results/resources/hg19.fasta.sa"
   File ref_amb = "s3://bgm-cromwell-results/resources/hg19.fasta.amb"
   File ref_ann = "s3://bgm-cromwell-results/resources/hg19.fasta.ann"
   File ref_pac = "s3://bgm-cromwell-results/resources/hg19.fasta.pac"

   String data_base = "s3://bgm-cromwell-results/data/"

   File dbSNP_vcf = "s3://bgm-cromwell-results/resources/dbsnp_147.hg19.vcf.gz"
   File dbSNP_vcf_idx = "s3://bgm-cromwell-results/resources/dbsnp_147.hg19.vcf.gz.tbi"

   File hapmap_vcf = "s3://bgm-cromwell-results/resources/hapmap_3.3.hg19.sites.vcf"
   File hapmap_vcf_idx = "s3://bgm-cromwell-results/resources/hapmap_3.3.hg19.sites.vcf.idx"

   File omni_vcf = "s3://bgm-cromwell-results/resources/1000G_omni2.5.hg19.sites.vcf"
   File omni_vcf_idx = "s3://bgm-cromwell-results/resources/1000G_omni2.5.hg19.sites.vcf.idx"

   File one_thousand_genomes_vcf = "s3://bgm-cromwell-results/resources/1000G_phase1.indels.hg19.vcf"
   File one_thousand_genomes_vcf_idx = "s3://bgm-cromwell-results/resources/1000G_phase1.indels.hg19.vcf.idx"

   File mills_vcf = "s3://bgm-cromwell-results/resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
   File mills_vcf_index = "s3://bgm-cromwell-results/resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx"

   #Array[String] chromosomes = ["chrM", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
   Array[String] chromosomes = ["chr20", "chr21", "chr22"]

   #File intervals_file = "s3://bgm-cromwell-results/resources/for.aws.long.txt"
   #File intervals_file = "/net/bgm/resources/variant_caller_intervals.file"

   File intervals_file = "s3://bgm-cromwell-results/resources/chrs.txt"

   Array[String] chr_intervals = read_lines(intervals_file)

   File fam = "s3://bgm-cromwell-results/data/bgm0102.fam"

   #Array[String] recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"]
   Array[String] recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0"]
   Array[String] recalibration_annotation_values = ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]

   #Array[String] SNPsVariantRecalibrator_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0"]
   Array[String] SNPsVariantRecalibrator_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.0"]
   Array[String] SNPsVariantRecalibrator_recalibration_annotation_values = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]

   File case_ids_for_test = "s3://bgm-cromwell-results/resources/case_ids_for_test.txt"

   Array[String] case_sample_names = ["bgm0102a1", "bgm0102u1", "bgm0102u2"]

   File ab_caller_py = "s3://bgm-cromwell-results/tools/ab_caller.py"
   File call_gathering_py = "s3://bgm-cromwell-results/tools/callsGathering.py"

   File exac_vcf = "s3://bgm-cromwell-results/resources/ExAC.r0.3.sites.vep.b37.vcf.gz"
   File exac_vcf_index = "s3://bgm-cromwell-results/resources/ExAC.r0.3.sites.vep.b37.vcf.gz.tbi"
   File cache = "s3://bgm-cromwell-results/resources/homo_sapiens_vep_94_GRCh37.tar.gz"

   scatter (index in range(length(chr_intervals)))
   {
      Int i = index
      String chr = chr_intervals[i]
      String sub_interval = sub(chr, ":", "_")

      #call GetChrGVCFs
      #{
      #   input:
      #     workflow_id = "3",
      #     chr = sub(chr, ":.+", "") + ".",
      #     folder_to_search = gvcfsFolder,
      #     tools = tools
      #}

      call DB_Import
      {
         input:
           gvcfs_archive = data_base + "wes." +  sub(chr_intervals[i], ":.+", "") + ".tar.gz",
           #dir = "/gatk/net/home/isaevt/mem/" + sub(chr_intervals[i], ":.+", ""),
           interval = chr_intervals[i],
           db_path = sub(chr_intervals[i], ":.+", "") + "_" + i + "db",
           padding = 100,
           #ref_fasta = ref_fasta,
           #ref_dict = ref_dict,
           #ref_fasta_index = ref_fasta_index,
           #ref_bwt = ref_bwt,
           #ref_sa = ref_sa,
           #ref_amb = ref_amb,
           #ref_ann = ref_ann,
           #ref_pac = ref_pac,
           docker_in = "broadinstitute/gatk:latest",
           gatk_launch = "/gatk/gatk",
           cpu = 2,
           memory = "6000 MB"
      }
  
      call GenotypeGVCFs
      {
         input:
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_fasta_index = ref_fasta_index,
           ref_bwt = ref_bwt,
           ref_sa = ref_sa,
           ref_amb = ref_amb,
           ref_ann = ref_ann,
           ref_pac = ref_pac,
           docker_in = "broadinstitute/gatk:latest",
           gatk_launch = "/gatk/gatk",
           cpu = 1,
           memory = "4096 MB",
           chr = chr_intervals[i],
           out_base_name = sub(chr_intervals[i], ":.+", "") + "_" + i,
           dbSNP_vcf = dbSNP_vcf,
           dbSNP_vcf_idx = dbSNP_vcf_idx,
           workspace_tar = DB_Import.db 
      }
   }

  call MergeGenotypeRes
  {
     input:
       vcfs = GenotypeGVCFs.output_vcf,
       vcfs_ind = GenotypeGVCFs.output_vcf_idx,
       out_base_name = "genotyped",
       out_extension = "vcf.gz",
       num_threads = 2,
       memory = "4096 MB",
       docker = "timuris/bcftools:bgzip"
  }

  #File out_vcf = "s3://bgm-cromwell-results/cromwell-execution/wes_downstream/5ee7e602-aecc-45ca-89c5-e45e7d515ec9/call-MergeGenotypeRes/genotyped.vcf.gz"
  #File out_vcf_index = "s3://bgm-cromwell-results/cromwell-execution/wes_downstream/5ee7e602-aecc-45ca-89c5-e45e7d515ec9/call-MergeGenotypeRes/genotyped.vcf.gz.tbi"

  call IndelsVariantRecalibrator
  {
     input:
       recalibration_filename = "indels.recalibration.file",
       tranches_filename = "indels.tranches.file",
       recalibration_tranche_values = recalibration_tranche_values,
       recalibration_annotation_values = recalibration_annotation_values,
       vcf = MergeGenotypeRes.out_vcf,
       vcf_index = MergeGenotypeRes.out_vcf_index,
       #vcf = out_vcf,
       #vcf_index = out_vcf_index,
       dbsnp_resource_vcf = dbSNP_vcf,
       dbsnp_resource_vcf_index = dbSNP_vcf_idx,
       mills_resource_vcf = mills_vcf,
       mills_resource_vcf_index = mills_vcf_index,
       gatk_launch = "/gatk/gatk",
       cpu = 2,
       docker = "broadinstitute/gatk:latest",
       memory = "16 GB"
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
       #vcf = out_vcf,
       #vcf_index = out_vcf_index,
       dbsnp_resource_vcf = dbSNP_vcf,
       dbsnp_resource_vcf_index = dbSNP_vcf_idx,
       hapmap_resource_vcf = hapmap_vcf,
       hapmap_resource_vcf_index = hapmap_vcf_idx,
       omni_resource_vcf = omni_vcf,
       omni_resource_vcf_index = omni_vcf_idx,
       one_thousand_genomes_resource_vcf = one_thousand_genomes_vcf,
       one_thousand_genomes_resource_vcf_index = one_thousand_genomes_vcf_idx,
       gatk_launch = "/gatk/gatk",
       cpu = 2,
       docker = "broadinstitute/gatk:latest",
       memory = "4 GB"
  }

  call ApplyRecalibration
  {
     input:
       recalibrated_vcf_filename = "recalibrated.vcf.gz",
       input_vcf =  MergeGenotypeRes.out_vcf,
       input_vcf_index = MergeGenotypeRes.out_vcf_index,
  
       #input_vcf = out_vcf,
       #input_vcf_index = out_vcf_index,
  
       indels_recalibration = IndelsVariantRecalibrator.recalibration,
       indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
       indels_tranches = IndelsVariantRecalibrator.tranches,
   
       snps_recalibration = SNPsVariantRecalibrator.recalibration,
       snps_recalibration_index = SNPsVariantRecalibrator.recalibration_index,
       snps_tranches = SNPsVariantRecalibrator.tranches,
   
       indel_filter_level = 95.0,
       snp_filter_level   = 99.0,
   
       gatk_launch = "/gatk/gatk",
       docker = "broadinstitute/gatk:latest",
       cpu = 1,
       memory = "4 GB"   
  }

  #call VepAnnotate
  #{
  #   input:
  #     number_of_threads = 8,
  #     cache_path  = "/opt/vep/.vep",
  #     plugin_path = "/opt/vep/.vep/Plugins",
  #     exac_path = "/opt/vep/src/ensembl-vep/ExAC.r0.3.sites.vep.b37.vcf.gz",
  #     input_vcf = ApplyRecalibration.recalibrated_vcf,
  #     input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
  #     output_basename = "final.vep.annotated"
  #}

  #call FindAllBams
  #{
  #   input:
  #     dir_to_search  = "/net/bgm/cases/",
  #     platform       = "wes",
  #     tools          = tools
  #}

  #File test = "s3://bgm-cromwell-results/cromwell-execution/wes_downstream/1c8bc563-6562-40d4-bfe3-6036a2e1d3cc/call-VepAnnotate/shard-0/chr22.vep.vcf"
  #File test = "s3://bgm-cromwell-results/data/chr22.vep.vcf"  

  scatter (chr in chromosomes)
  {
     call SplitVcfByChr
     {
        input:
          in_vcf = ApplyRecalibration.recalibrated_vcf,
          in_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          #in_vcf = VepAnnotate.recalibrated_vcf,
          #in_vcf_index = VepAnnotate.recalibrated_vcf_index,
          chr = chr,
          memory = "1 GB",
          docker = "timuris/bwa"
     }
  
     call VepAnnotate
     {
        input:
          exac_vcf = exac_vcf,
          exac_vcf_index = exac_vcf_index,
          cache = cache,
          input_vcf = SplitVcfByChr.vcf,
          input_vcf_idx = SplitVcfByChr.vcf_idx,
          out_name = chr + ".vep.vcf.gz",
          buffer = 50000,
          cpu = 4,
          memory = "16 GB",
          docker = "timuris/vep:empty"
     }     

     call TabixIndexing
     {
        input:
          vcf_in = VepAnnotate.out,
          chr = chr,
          #vcf_in = test,
          memory = "512 MB",
          docker = "timuris/bwa"
     }

     call BgmBayesDeNovo_stage1
     {
        input:
          case_ids  = case_ids_for_test,
          input_vcf = TabixIndexing.vcf, 
          x = 1.0,
          p = 0.005,
          e = 0.05,
          docker = "timuris/novo_caller:full",
          cpu = 1,
          memory = "1 GB",
          command = "novoCaller"
     }
  
     call ProcessBayesDeNovoSt1Res
     {
        input:
          res = BgmBayesDeNovo_stage1.stage_one_out,
          docker = "timuris/python2",
          memory = "512 MB",
          script = "s3://bgm-cromwell-results/tools/processBayesDeNovo.py"
     }
  
     call BgmCompoundHetCaller
     {
        input:
          case_ids = case_ids_for_test,
          input_vcf = TabixIndexing.vcf,
          t = 0.9,
          e = 0.1,
          a = 0.2,
          docker = "timuris/novo_caller:full",
          cpu = 1,
          memory = "1 GB"
     }
  
     call ProcessBgmCompHetRes
     {
        input:
          res = BgmCompoundHetCaller.out,
          docker = "timuris/python2",
          memory = "512 MB",
          script = "s3://bgm-cromwell-results/tools/processBgmCompHet.py",
          cpu = 1,
          vcf = TabixIndexing.vcf 
     }
  
     call BgmHomRec
     {
        input:
          case_ids = case_ids_for_test,
          input_vcf = TabixIndexing.vcf,
          x = 1.0,
          p = 0.9,
          e = 0.05,
          a = 0.2,
          docker = "timuris/novo_caller:full",
          cpu = 1,
          memory = "1 GB"
     }
  
     call ProcessHomRecRes
     {
        input:
          res = BgmHomRec.out,
          docker = "timuris/python2",
          script = "s3://bgm-cromwell-results/tools/processHomRec.py",
          memory = "512 MB",
          cpu = 1
     }
  
     call ABCaller as AutosomalDominantCaller
     {
        input:
          in_vcf = TabixIndexing.vcf,
          fam    = fam,
          id     = 0,
          py_script = ab_caller_py,
          output_name = "autosomalDominant.calls",
          docker = "timuris/python2" 
     }
  
     call ABCaller as HomozygousRecessiveCaller
     {
        input:
          in_vcf = TabixIndexing.vcf,
          fam    = fam,
          id     = 1,
          py_script = ab_caller_py,
          output_name = "homozygousRecessive.calls",
          docker = "timuris/python2"
     }
  
     call ABCaller as DeNovoCaller
     {
        input:
          in_vcf = TabixIndexing.vcf,
          fam    = fam,
          id     = 2,
          py_script = ab_caller_py,
          output_name = "deNovoCaller.calls",
          docker = "timuris/python2"
     }
  
     call ABCaller as CompoundHeterozygousCaller
     {
        input:
          in_vcf = TabixIndexing.vcf,
          fam    = fam,
          id     = 3,
          py_script = ab_caller_py,
          output_name = "compoundHeterozygous.calls",
          docker = "timuris/python2"
     } 
  
     call GatherCalls
     {
        input:
          bgmHomRec_res = ProcessHomRecRes.hom_rec_calls,
          compoundHeterozygous_res = CompoundHeterozygousCaller.calls,
          deNovo_res = DeNovoCaller.calls,
          homozygousRecessive_res = HomozygousRecessiveCaller.calls,
          autosomalDominant_res = AutosomalDominantCaller.calls,
          bgmCompoundHet_res = ProcessBgmCompHetRes.comp_het_calls,
          bgmDeNovo_res = ProcessBayesDeNovoSt1Res.bayes_deNovo_calls_st1,
          vcf = ApplyRecalibration.recalibrated_vcf,
          vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          docker = "timuris/python2",
          py_script = call_gathering_py,
          memory = "512 MB",
          cpu = 1
     }
  
     call Bgzip
     {
        input:
          docker = "timuris/bwa",
          memory = "1 GB",
          cpu = 1,
          vcf = GatherCalls.vcf_out,
          out_base_name = "calls." + chr
     }

     call MergeVcfCalls
     {
        input:
          calls = Bgzip.vcf_out,
          calls_idx = Bgzip.vcf_out_index,
          vcf = TabixIndexing.vcf,
          vcf_index = TabixIndexing.vcf_idx,
          out_base_name = chr + ".calls",
          docker = "timuris/bcftools:bgzip",
          memory = "1024 MB",
          cpu = 1,
     }

  }

  call BcfToolsMergeVcf
  {
     input:
       vcfs = MergeVcfCalls.vcf_out,
       vcfs_ind = MergeVcfCalls.vcf_out_index,
       out_base_name = "calls",
       out_extension = "vcf.gz",
       out_type = "z",
       num_threads = 4,
       memory = "2 GB",
       docker = "timuris/bcftools:bgzip"
  }

  ##call MergeVcfCalls
  ##{
  ##   input:
  ##     calls = BcfToolsMergeVcf.out_vcf,
  ##     vcf = VepAnnotate.recalibrated_vcf,
  ##     vcf_index = VepAnnotate.recalibrated_vcf_index,
  ##     tools = tools 
  ##}

  #File pre_final_vcf = "s3://bgm-cromwell-results/cromwell-execution/wes_downstream/f0fdb229-5b56-435f-aa04-82cc8d3f692f/call-BcfToolsMergeVcf/calls.vcf"

  call SelectCaseVariants
  {
     input:
       ref_fasta = ref_fasta,
       ref_dict = ref_dict,
       ref_fasta_index = ref_fasta_index,
       ref_bwt = ref_bwt,
       ref_sa = ref_sa,
       ref_amb = ref_amb,
       ref_ann = ref_ann,
       ref_pac = ref_pac,
       docker = "broadinstitute/gatk:latest",
       memory = "2048 MB",
       cpu = 1,
       gatk_launch = "/gatk/gatk",
       in_vcf = BcfToolsMergeVcf.out_vcf,
       in_vcf_idx = BcfToolsMergeVcf.out_vcf_idx,
       #in_vcf = pre_final_vcf,
       samples = case_sample_names,
       caller_names = ["BGM_DE_NOVO", "BGM_CMPD_HET", "BGM_HOM_REC", "BGM_AUTO_DOM", "BGM_BAYES_HOM_REC", "BGM_BAYES_DE_NOVO", "BGM_BAYES_CMPD_HET"],
       out_base_name = "xbrowse.vep"
  }

  #output
  #{
  #   File vcf = VepAnnotate.recalibrated_vcf
  #   File vcf_idx = VepAnnotate.recalibrated_vcf_index
  #}
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

#task DB_Import_Big
#{
#  File gvcfs_archive
#
#  File intervals
#
#  Int padding
#
#  String docker_in
#  String gatk_launch
#
#  Int cpu
#  String memory
#
#  String db_path
#
#  command
#  <<<
#     tar -xzf ${gvcfs_archive} && \
#     find "$(pwd)" -name "*gvcf.gz" -type f -exec echo -V {} \; > variants.txt && \
#     value=`cat variants.txt`
#
#     l=`awk '{ print "-L " $0 }' < ${intervals}`
#
#     ${gatk_launch} --java-options "-Xms2g -Xmx4g" GenomicsDBImport \
#       --batch-size 50 \
#       -ip ${padding} \
#       -R ${ref_fasta} \
#       --max-num-intervals-to-import-in-parallel 4
#       --reader-threads 30 \
#       --overwrite-existing-genomicsdb-workspace true \
#       $l \
#       --genomicsdb-workspace-path ${db_path} \
#       $value
#
#     tar -cf ${db_path}.tar ${db_path}
#  >>>
#
#  runtime
#  {
#     cpu: cpu
#     memory: memory
#     docker: docker_in
#  }
#
#  output
#  {
#     File db = "${db_path}.tar"
#  }
#}

task DB_Import
{
  File gvcfs_archive
  #String dir

  String interval

  Int padding
 
  #File ref_fasta
  #File ref_dict
  #File ref_fasta_index
  #File ref_bwt
  #File ref_sa
  #File ref_amb
  #File ref_ann
  #File ref_pac

  String docker_in
  String gatk_launch

  Int cpu
  String memory

  String db_path 

  command
  <<<
     tar -xzf ${gvcfs_archive} && \
     find "$(pwd)" -name "*gvcf.gz" -type f -exec echo -V {} \; > variants.txt && \
     value=`cat variants.txt`

     ${gatk_launch} --java-options "-Xms2g -Xmx4g" GenomicsDBImport \
       --batch-size 40 \
       -ip ${padding} \
       --reader-threads 5 \
       --overwrite-existing-genomicsdb-workspace true \
       -L ${interval} \
       --genomicsdb-workspace-path ${db_path} \
       $value

     tar -cf ${db_path}.tar ${db_path}
  >>>

  runtime
  {
     cpu: cpu
     memory: memory
     docker: docker_in
  }

  output
  {
     File db = "${db_path}.tar"
  }

}



task SelectCaseVariants
{
   String gatk_launch
   
   File ref_fasta
   File ref_dict
   File ref_fasta_index
   File ref_bwt
   File ref_sa
   File ref_amb
   File ref_ann
   File ref_pac

   File in_vcf
   File in_vcf_idx
   Array[String] samples
   Array[String] caller_names
   String out_base_name

   String docker
   String memory
   Int cpu

   command
   {
      ${gatk_launch} --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap -XX:MaxRAMFraction=1" \
        SelectVariants \
        -R ${ref_fasta} \
        -sn ${sep=" -sn " samples} \
        -select "(!vc.isMonomorphicInSamples()) || ${sep=" || " caller_names}" \
        -V ${in_vcf} \
        -O ${out_base_name}.vcf.gz  
   }

   runtime
   {
      memory: memory
      cpu: cpu
      docker: docker
   }

   output
   {
      File out_vcf = "${out_base_name}.vcf.gz"
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
   File vcf
   File vcf_index
   File py_script
   String docker
   String memory
   Int cpu

   command
   {
      python2.7 ${py_script} ${deNovo_res} ${bgmHomRec_res} ${compoundHeterozygous_res} ${autosomalDominant_res} ${bgmDeNovo_res} ${homozygousRecessive_res} ${autosomalDominant_res} ${vcf}
   }

   runtime
   {
      memory: memory
      cpu: cpu
      docker: docker
   }

   output
   {
      File vcf_out = "calls.vcf"
   } 
}

task MergeGenotypeRes
{
   Array[File] vcfs
   Array[File] vcfs_ind
   String out_base_name
   String out_extension

   Int num_threads
   String memory
   String docker

   command
   {
      set -e
      bcftools concat --threads ${num_threads} -a -D -O z -o ${out_base_name}.${out_extension} ${sep=' ' vcfs}
      tabix -p vcf ${out_base_name}.${out_extension}
   }

   runtime
   {
      memory: memory
      cpu: num_threads
      docker: docker
   }

   output
   {
      File out_vcf = "${out_base_name}.${out_extension}"
      File out_vcf_index = "${out_base_name}.vcf.gz.tbi"
   }
}

task BcfToolsMergeVcf
{
   Array[File] vcfs
   Array[File] vcfs_ind
   String out_base_name
   String out_extension
   String out_type

   Int num_threads
   String memory
   String docker

   command
   {
      bcftools concat --threads ${num_threads} -a -D -O ${out_type} -o ${out_base_name}.${out_extension} ${sep=' ' vcfs}
      if [ "${out_type}" = "z" ]; then
         tabix -p vcf ${out_base_name}.${out_extension}
      fi 
   }

   runtime
   {
      memory: memory
      cpu: num_threads
      docker: docker
   }

   output
   {
      File out_vcf = "${out_base_name}.${out_extension}"
      File? out_vcf_idx = "${out_base_name}.vcf.gz.tbi"
   }
}

##################
##################

task TabixIndexing
{
   File vcf_in
   String docker
   String memory
   String chr

   command
   {
      set -e
      zcat ${vcf_in} > decompressed.vcf
      bgzip -c decompressed.vcf > ${chr}.decompressed.vcf.gz
      tabix -p vcf ${chr}.decompressed.vcf.gz
   }

   runtime
   {
      memory: memory
      cpu: 1
      docker: docker
   }

   output
   {
      File vcf = "${chr}.decompressed.vcf.gz"
      File vcf_idx = "${chr}.decompressed.vcf.gz.tbi"
   }
}

task Bgzip
{
  File vcf
  String out_base_name

  String docker
  String memory
  Int cpu  

  command
  {
     bgzip -c ${vcf} > ${out_base_name}.vcf.gz 
     tabix -p vcf ${out_base_name}.vcf.gz
  }

  runtime
  {
     memory: memory
     cpu: cpu
     docker: docker
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
   File calls_idx
   File vcf
   File vcf_index

   String out_base_name

   String docker
   String memory
   Int cpu

   command
   {
      set -e
      bcftools annotate -a ${calls} -c INFO ${vcf} > ${out_base_name}.vcf
      bgzip -c ${out_base_name}.vcf > ${out_base_name}.vcf.gz
      tabix -p vcf ${out_base_name}.vcf.gz      
   }

   runtime
   {
      memory: memory
      cpu: cpu
      docker: docker
   }

   output
   {
      File vcf_out = "${out_base_name}.vcf.gz"
      File vcf_out_index = "${out_base_name}.vcf.gz.tbi"
   }
} 

task ABCaller
{
   File in_vcf
   File fam
   File py_script
   Int id
   String output_name
   String docker

   command
   {
      set -e
      zcat ${in_vcf} > decompressed.vcf
      python2.7 ${py_script} decompressed.vcf ${fam} ${id} > ${output_name}
   }

   runtime
   {
      memory: "1 GB"
      cpu: 1
      docker: docker
   }

   output
   {
      File calls = "${output_name}"
   }
}

task ProcessHomRecRes
{
   File res
   File script

   String docker
   String memory
   Int cpu

   command
   <<<
     python2.7 ${script} ${res} > hom_rec_calls.txt
   >>>

   runtime
   {
      memory: memory
      cpu: cpu
      docker: docker
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
   File script

   String docker
   String memory
   Int cpu

   command
   <<<
     set -e
     zcat ${vcf} > decompressed.vcf
     python2.7 ${script} ${res} decompressed.vcf > bgm_comp_het_calls.txt
   >>>

   runtime
   {
      memory: memory
      cpu: cpu
      docker: docker
   }

   output
   {
      File comp_het_calls = "bgm_comp_het_calls.txt"
   }
}

task ProcessBayesDeNovoSt1Res
{
   File res
   File script
   String docker
   String memory

   command
   <<<
     python2.7 ${script} ${res} > bayes_de_novo1_calls.txt 
   >>>

   runtime
   {
     memory: memory
     cpu: 1
     docker: docker
   }

   output
   {
     File bayes_deNovo_calls_st1 = "bayes_de_novo1_calls.txt"
   }
}

task BgmBayesDeNovo_stage1
{
  File case_ids
  File input_vcf

  Float x
  Float p
  Float e  

  String docker
  Int cpu
  String memory

  String command

  command
  {
     zcat ${input_vcf} > decompressed.vcf
     ${command} -I decompressed.vcf -O bayes_de_novo_st1.out.vcf -T ${case_ids} -X ${x} -P ${p} -E ${e}
  }

  runtime
  {
     memory: memory
     cpu: cpu
     docker: docker
  }

  output
  {
     File stage_one_out = "bayes_de_novo_st1.out.vcf"
  }
}

task BgmCompoundHetCaller
{
  File case_ids
  File input_vcf
  Float t
  Float e
  Float a
  
  String docker
  String memory
  Int cpu

  command
  {
     set -e
     zcat ${input_vcf} > decompressed.vcf
     cmpHet -I decompressed.vcf -O bgm_cmp_het.out -TR ${case_ids} -T ${t} -E ${e} -A ${a}
  }

  runtime
  {
     memory: memory
     cpu: cpu
     docker: docker
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

   String docker
   String memory
   Int cpu

   command
   {
      set -e 
      zcat ${input_vcf} > decompressed.vcf
      homrec -I decompressed.vcf -O bgm_hom_rec.out -T ${case_ids} -X ${x} -P ${p} -E ${e} -A ${a}
   }

   runtime
   {
      memory: memory
      cpu: cpu
      docker: docker 
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
   
   String docker
   String memory

   command
   {
      tabix -h ${in_vcf} ${chr} > ${chr}.vcf
      bgzip -c ${chr}.vcf > ${chr}.vcf.gz
      tabix -p vcf ${chr}.vcf.gz
   }

   runtime
   {
      memory: memory
      cpu: 1
      docker: docker
   }

   output
   {
      File vcf = "${chr}.vcf.gz"
      File vcf_idx = "${chr}.vcf.gz.tbi"
   }
}

#task VepAnnotate
#{
#   String cache_path
#   String exac_path
#   String plugin_path
#   File input_vcf
#   File input_vcf_index
#   String output_basename
#   Int number_of_threads
#
#   command
#   {
#      ./vep --format vcf --force_overwrite --dir_cache ${cache_path} --dir_plugin ${plugin_path} --offline --no_stats --vcf --everything \
#            --plugin ExAC,${exac_path} --allele_number --input_file ${input_vcf} --output_file ${output_basename}.vcf.gz --tabix --fork ${number_of_threads} --buffer_size 50000
#   }
#
#   runtime
#   {
#     memory: "16 GB"
#     cpu: number_of_threads
#     docker: "timuris/vep"
#   }
#
#   output
#   {
#     File recalibrated_vcf = "${output_basename}.vcf.gz"
#     File recalibrated_vcf_index = "${output_basename}.vcf.gz.tbi"
#   }
#
#}

task VepAnnotate
{
   File exac_vcf
   File exac_vcf_index
   File cache
   File input_vcf
   File input_vcf_idx
   String out_name
   Int buffer

   Int cpu
   String memory
   String docker

   command
   {
      tar -xzf ${cache}
      cd /opt/vep/src/ensembl-vep
      perl INSTALL.pl -a p --PLUGINS ExAC

      vep --format vcf --force_overwrite --dir_cache /cromwell_root/ --offline --no_stats --vcf --everything \
          --plugin ExAC,${exac_vcf} --allele_number --input_file ${input_vcf} --compress_output gzip --output_file /cromwell_root/${out_name} -v --fork ${cpu} --buffer_size ${buffer}
   }

   runtime
   {
     memory: memory
     cpu: cpu
     docker: docker
   }

   output
   {
     File out = "${out_name}"
     File out_index = "${out_name}.tbi"
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

  String gatk_launch
  String docker
  Int cpu
  String memory

  command 
  {
    set -e
    ${gatk_launch} --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap -XX:MaxRAMFraction=1" \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ${input_vcf} \
      --recal-file ${indels_recalibration} \
      --tranches-file ${indels_tranches} \
      --truth-sensitivity-filter-level ${indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    ${gatk_launch} --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap -XX:MaxRAMFraction=1" \
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
    memory: memory
    cpu: cpu
    docker: docker
  }

  output 
  {
    File recalibrated_vcf = "${recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "${recalibrated_vcf_filename}.tbi"
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

  Int cpu
  String memory
  String docker

  command
  {
    ${gatk_launch} --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap -XX:MaxRAMFraction=1" \
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
    memory: memory
    cpu: cpu
    docker: docker
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
   String docker
   String memory
   Int cpu

   command
   {
      ${gatk_launch} --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap -XX:MaxRAMFraction=1" \
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
      memory: memory
      cpu: cpu
      docker: docker
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
   File ref_fasta
   File ref_dict
   File ref_fasta_index
   File ref_bwt
   File ref_sa
   File ref_amb
   File ref_ann
   File ref_pac

   String docker_in
   String gatk_launch

   Int cpu
   String memory

   String chr
   String out_base_name   

   File dbSNP_vcf
   File dbSNP_vcf_idx

   File workspace_tar

   command
   {
      set -e
      tar -xf ${workspace_tar}
      WORKSPACE=$( basename ${workspace_tar} .tar)

      ${gatk_launch} --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap -XX:MaxRAMFraction=1" GenotypeGVCFs \
                     -R ${ref_fasta} \
                     -V gendb://$WORKSPACE \
                     -D ${dbSNP_vcf} \
                     -O ${out_base_name}.raw.vcf.gz \
                     -L ${chr} \
                     -G StandardAnnotation \
                     -new-qual true
   }

   output
   {
      File output_vcf = "${out_base_name}.raw.vcf.gz"
      File output_vcf_idx = "${out_base_name}.raw.vcf.gz.tbi"
   }

   runtime
   {
      memory: memory
      cpu: cpu
      docker: docker_in
   }
}

