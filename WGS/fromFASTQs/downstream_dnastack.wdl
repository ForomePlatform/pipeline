workflow wes_downstream
{
   File ref_fasta
   File ref_dict
   File ref_fasta_index
   File ref_bwt
   File ref_sa
   File ref_amb
   File ref_ann
   File ref_pac

   File data_base

   File dbSNP_vcf
   File dbSNP_vcf_idx

   File hapmap_vcf
   File hapmap_vcf_idx

   File omni_vcf
   File omni_vcf_idx

   File one_thousand_genomes_vcf
   File one_thousand_genomes_vcf_idx

   File mills_vcf
   File mills_vcf_index

   Array[String] chromosomes = ["chr22"]

   File intervals_file

   Array[String] chr_intervals = read_lines(intervals_file)

   File fam

   Array[String] recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0"]
   Array[String] recalibration_annotation_values = ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]

   Array[String] SNPsVariantRecalibrator_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.0"]
   Array[String] SNPsVariantRecalibrator_recalibration_annotation_values = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]

   File case_ids_for_test

   Array[String] case_sample_names

   File ab_caller_py
   File call_gathering_py
   File processBayesDenovo_py
   File processBgmCompHet_py
   File processHomRec_py

   File exac_vcf
   File exac_vcf_index

   scatter (index in range(length(chr_intervals)))
   {
      Int i = index
      String chr = chr_intervals[i]
      String sub_interval = sub(chr, ":", "_")

      call DB_Import
      {
         input:
           gvcfs_archive = data_base,
           interval = chr_intervals[i],
           db_path = sub(chr_intervals[i], ":.+", "") + "_" + i + "db",
           padding = 100,
           docker_in = "broadinstitute/gatk:latest",
           gatk_launch = "/gatk/gatk",
           cpu = 2,
           memory = "7.5 GB"
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
           cpu = 2,
           memory = "7.5 GB",
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
       memory = "7.5 GB",
       docker = "timuris/bcftools:bgzip"
  }

  call IndelsVariantRecalibrator
  {
     input:
       recalibration_filename = "indels.recalibration.file",
       tranches_filename = "indels.tranches.file",
       recalibration_tranche_values = recalibration_tranche_values,
       recalibration_annotation_values = recalibration_annotation_values,
       vcf = MergeGenotypeRes.out_vcf,
       vcf_index = MergeGenotypeRes.out_vcf_index,
       dbsnp_resource_vcf = dbSNP_vcf,
       dbsnp_resource_vcf_index = dbSNP_vcf_idx,
       mills_resource_vcf = mills_vcf,
       mills_resource_vcf_index = mills_vcf_index,
       gatk_launch = "/gatk/gatk",
       cpu = 4,
       docker = "broadinstitute/gatk:4.1.0.0",
       memory = "26 GB"
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
       docker = "broadinstitute/gatk:4.1.0.0",
       memory = "7.5 GB"
  }

  call ApplyRecalibration
  {
     input:
       recalibrated_vcf_filename = "recalibrated.vcf.gz",
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
       gatk_launch = "/gatk/gatk",
       docker = "broadinstitute/gatk:latest",
       cpu = 2,
       memory = "7.5 GB"   
  }

  scatter (chr in chromosomes)
  {
     call SplitVcfByChr
     {
        input:
          in_vcf = ApplyRecalibration.recalibrated_vcf,
          in_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          chr = chr,
          memory = "3.75 GB",
          docker = "timuris/bwa"
     }
  
     call VepAnnotate
     {
        input:
          exac_vcf = exac_vcf,
          exac_vcf_index = exac_vcf_index,
          input_vcf = SplitVcfByChr.vcf,
          input_vcf_idx = SplitVcfByChr.vcf_idx,
          out_name = chr + ".vep.vcf.gz",
          buffer = 50000,
          cpu = 4,
          memory = "26 GB",
          docker = "timuris/vep:cache"
     }     

     call BgmBayesDeNovo_stage1
     {
        input:
          case_ids  = case_ids_for_test,
          input_vcf = VepAnnotate.out, 
          x = 1.0,
          p = 0.005,
          e = 0.05,
          docker = "timuris/novo_caller:full",
          cpu = 1,
          memory = "3.75 GB",
          command = "novoCaller"
     }
  
     call ProcessBayesDeNovoSt1Res
     {
        input:
          res = BgmBayesDeNovo_stage1.stage_one_out,
          docker = "timuris/python2",
          memory = "3.75 GB",
          script = processBayesDenovo_py
     }
  
     call BgmCompoundHetCaller
     {
        input:
          case_ids = case_ids_for_test,
          input_vcf = VepAnnotate.out,
          t = 0.9,
          e = 0.1,
          a = 0.2,
          docker = "timuris/novo_caller:full",
          cpu = 1,
          memory = "3.75 GB"
     }

     call ProcessBgmCompHetRes
     {
        input:
          res = BgmCompoundHetCaller.out,
          docker = "timuris/python2",
          memory = "3.75 GB",
          script = processBgmCompHet_py,
          cpu = 1,
          vcf = VepAnnotate.out 
     }
  
     call BgmHomRec
     {
        input:
          case_ids = case_ids_for_test,
          input_vcf = VepAnnotate.out,
          x = 1.0,
          p = 0.9,
          e = 0.05,
          a = 0.2,
          docker = "timuris/novo_caller:full",
          cpu = 1,
          memory = "3.75 GB"
     }
  
     call ProcessHomRecRes
     {
        input:
          res = BgmHomRec.out,
          docker = "timuris/python2",
          script = processHomRec_py,
          memory = "3.75 GB",
          cpu = 1
     }
  
     call ABCaller as AutosomalDominantCaller
     {
        input:
          in_vcf = VepAnnotate.out,
          fam    = fam,
          id     = 0,
          py_script = ab_caller_py,
          output_name = "autosomalDominant.calls",
          docker = "timuris/python2" 
     }
  
     call ABCaller as HomozygousRecessiveCaller
     {
        input:
          in_vcf = VepAnnotate.out,
          fam    = fam,
          id     = 1,
          py_script = ab_caller_py,
          output_name = "homozygousRecessive.calls",
          docker = "timuris/python2"
     }
  
     call ABCaller as DeNovoCaller
     {
        input:
          in_vcf = VepAnnotate.out,
          fam    = fam,
          id     = 2,
          py_script = ab_caller_py,
          output_name = "deNovoCaller.calls",
          docker = "timuris/python2"
     }
  
     call ABCaller as CompoundHeterozygousCaller
     {
        input:
          in_vcf = VepAnnotate.out,
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
          memory = "3.75 GB",
          cpu = 1
     }
  
     call Bgzip
     {
        input:
          docker = "timuris/bwa",
          memory = "3.75 GB",
          cpu = 1,
          vcf = GatherCalls.vcf_out,
          out_base_name = "calls." + chr
     }

     call MergeVcfCalls
     {
        input:
          calls = Bgzip.vcf_out,
          calls_idx = Bgzip.vcf_out_index,
          vcf = VepAnnotate.out,
          vcf_index = VepAnnotate.out_index,
          out_base_name = chr + ".calls",
          docker = "timuris/bcftools:bgzip",
          memory = "3.75 GB",
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
       memory = "3.6 GB",
       docker = "timuris/bcftools:bgzip"
  }

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
       memory = "3.75 GB",
       cpu = 1,
       gatk_launch = "/gatk/gatk",
       in_vcf = BcfToolsMergeVcf.out_vcf,
       in_vcf_idx = BcfToolsMergeVcf.out_vcf_idx,
       samples = case_sample_names,
       caller_names = ["BGM_DE_NOVO", "BGM_CMPD_HET", "BGM_HOM_REC", "BGM_AUTO_DOM", "BGM_BAYES_HOM_REC", "BGM_BAYES_DE_NOVO", "BGM_BAYES_CMPD_HET"],
       out_base_name = "xbrowse.vep"
  }

  output {
    File out_vcf = SelectCaseVariants.out_vcf
  }

}


task DB_Import
{
  File gvcfs_archive

  String interval

  Int padding

  String docker_in
  String gatk_launch

  Int cpu
  String memory

  String db_path 

  Int disk_size = ceil(size(gvcfs_archive, "GB")*4 + 50)

  command
  <<<
     tar -zxvf ${gvcfs_archive} && \
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

     tar -cvf ${db_path}.tar ${db_path}
  >>>

  runtime
  {
     cpu: cpu
     memory: memory
     docker: docker_in
     disks: "local-disk " + disk_size + " HDD"
     preemptible: 4
  }

  output
  {
     File db = "${db_path}.tar"
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

   Int disk_size = ceil(size(ref_fasta, "GB")*8 + size(dbSNP_vcf, "GB")*2 + size(workspace_tar, "GB")*4 + 100)

   command
   {
      mv ${ref_fasta} ${ref_dict} ${ref_fasta_index} .
      set -e
      tar -xf ${workspace_tar}
      WORKSPACE=$( basename ${workspace_tar} .tar)

      ${gatk_launch} --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap -XX:MaxRAMFraction=1" GenotypeGVCFs \
                     -R $(basename ${ref_fasta}) \
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
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
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

   Int disk_size = ceil(size(vcfs[0], "GB")*length(vcfs)*3 + 50)

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
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
   }

   output
   {
      File out_vcf = "${out_base_name}.${out_extension}"
      File out_vcf_index = "${out_base_name}.vcf.gz.tbi"
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

   Int disk_size = ceil((size(vcf, "GB") + size(mills_resource_vcf, "GB") + size(dbsnp_resource_vcf, "GB"))*3 + 50)

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
                     --resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf} \
                     --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp_resource_vcf}
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
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
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

  Int disk_size = ceil((size(vcf, "GB")*2 + size(hapmap_resource_vcf, "GB") + size(omni_resource_vcf, "GB") + size(one_thousand_genomes_resource_vcf, "GB") + size(dbsnp_resource_vcf, "GB"))*2 + 50)

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
      --resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \
      --resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \
      --resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp_resource_vcf}
  }

  runtime
  {
    memory: memory
    cpu: cpu
    docker: docker
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 4
  }

  output
  {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
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

  Int disk_size = ceil((size(input_vcf, "GB")*2 + size(indels_recalibration, "GB") + size(snps_recalibration, "GB"))*2 + 50)

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
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 4
  }

  output 
  {
    File recalibrated_vcf = "${recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "${recalibrated_vcf_filename}.tbi"
  }
}


task SplitVcfByChr
{
   File in_vcf
   File in_vcf_index
   String chr
   
   String docker
   String memory

   Int disk_size = ceil(size(in_vcf, "GB")*4 + 20)

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
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
   }

   output
   {
      File vcf = "${chr}.vcf.gz"
      File vcf_idx = "${chr}.vcf.gz.tbi"
   }
}


task VepAnnotate
{
   File exac_vcf
   File exac_vcf_index
   File input_vcf
   File input_vcf_idx
   String out_name
   Int buffer

   Int cpu
   String memory
   String docker

   Int disk_size = ceil((size(exac_vcf, "GB") + size(input_vcf, "GB"))*2 + 50)

   command
   {
    perl /opt/vep/src/ensembl-vep/INSTALL.pl -n -a p --PLUGINS ExAC
    export PERL5LIB=$PERL5LIB:/root/.vep/Plugins/

    vep --format vcf --force_overwrite --dir_cache /opt/vep/.vep/ --offline --no_stats --vcf --everything \
        --plugin ExAC,${exac_vcf} --allele_number --input_file ${input_vcf} --compress_output bgzip --output_file ${out_name} -v --fork ${cpu} --buffer_size ${buffer}

    tabix ${out_name}

   }

   runtime
   {
     memory: memory
     cpu: cpu
     docker: docker
     disks: "local-disk " + disk_size + " HDD"
     preemptible: 4
     bootDiskSizeGb: 50
   }

   output
   {
     File out = "${out_name}"
     File out_index = "${out_name}.tbi"
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

  Int disk_size = ceil(size(input_vcf, "GB")*6 + 50)

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
     disks: "local-disk " + disk_size + " HDD"
     preemptible: 4
  }

  output
  {
     File stage_one_out = "bayes_de_novo_st1.out.vcf"
  }
}


task ProcessBayesDeNovoSt1Res
{
   File res
   File script
   String docker
   String memory

   Int disk_size = ceil(size(res, "GB")*2 + 20)

   command
   <<<
     python2.7 ${script} ${res} > bayes_de_novo1_calls.txt 
   >>>

   runtime
   {
     memory: memory
     cpu: 1
     docker: docker
     disks: "local-disk " + disk_size + " HDD"
     preemptible: 4
   }

   output
   {
     File bayes_deNovo_calls_st1 = "bayes_de_novo1_calls.txt"
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

  Int disk_size = ceil(size(input_vcf, "GB")*6 + 50)

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
     disks: "local-disk " + disk_size + " HDD"
     preemptible: 4
  }

  output
  {
     File out = "bgm_cmp_het.out"
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

   Int disk_size = ceil((size(res, "GB") + size(vcf, "GB"))*6 + 50)

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
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
   }

   output
   {
      File comp_het_calls = "bgm_comp_het_calls.txt"
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

   Int disk_size = ceil(size(input_vcf, "GB")*6 + 50)

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
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
   }

   output
   {
      File out = "bgm_hom_rec.out"
   }
}


task ProcessHomRecRes
{
   File res
   File script

   String docker
   String memory
   Int cpu

   Int disk_size = ceil(size(res, "GB")*2 + 20)

   command
   <<<
     python2.7 ${script} ${res} > hom_rec_calls.txt
   >>>

   runtime
   {
      memory: memory
      cpu: cpu
      docker: docker
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
   }

   output
   {
      File hom_rec_calls = "hom_rec_calls.txt"
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

   Int disk_size = ceil((size(in_vcf, "GB") + size(fam, "GB"))*6 + 50)

   command
   {
      set -e
      zcat ${in_vcf} > decompressed.vcf
      python2.7 ${py_script} decompressed.vcf ${fam} ${id} > ${output_name}
   }

   runtime
   {
      memory: "3.75 GB"
      cpu: 1
      docker: docker
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
   }

   output
   {
      File calls = "${output_name}"
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

   Int disk_size = ceil((size(vcf, "GB") + size(bgmHomRec_res, "GB") + size(compoundHeterozygous_res, "GB") + size(deNovo_res, "GB") + size(homozygousRecessive_res, "GB") + size(autosomalDominant_res, "GB") + size(bgmCompoundHet_res, "GB") + size(bgmDeNovo_res, "GB"))*2 + 50)

   command
   {
      python2.7 ${py_script} ${deNovo_res} ${bgmHomRec_res} ${compoundHeterozygous_res} ${autosomalDominant_res} ${bgmDeNovo_res} ${homozygousRecessive_res} ${autosomalDominant_res} ${vcf}
   }

   runtime
   {
      memory: memory
      cpu: cpu
      docker: docker
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
   }

   output
   {
      File vcf_out = "calls.vcf"
   } 
}


task Bgzip
{
  File vcf
  String out_base_name

  String docker
  String memory
  Int cpu

  Int disk_size = ceil(size(vcf, "GB")*3 + 20)

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
     disks: "local-disk " + disk_size + " HDD"
     preemptible: 4
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

   Int disk_size = ceil((size(calls, "GB") + size(vcf, "GB"))*4 + 50)

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
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
   }

   output
   {
      File vcf_out = "${out_base_name}.vcf.gz"
      File vcf_out_index = "${out_base_name}.vcf.gz.tbi"
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

   Int disk_size = ceil((size(vcfs[0], "GB")*length(vcfs))*3 + 50 )

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
      disks: "local-disk " + disk_size + " HDD"
      preemptible: 4
   }

   output
   {
      File out_vcf = "${out_base_name}.${out_extension}"
      File out_vcf_idx = "${out_base_name}.vcf.gz.tbi"
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
   
   String ref_fasta_basename = basename(ref_fasta)

   File in_vcf
   File in_vcf_idx
   Array[String] samples
   Array[String] caller_names
   String out_base_name

   String docker
   String memory
   Int cpu

   Int disk_size = ceil(size(ref_fasta, "GB")*8 + size(in_vcf, "GB")*3 + 50)

   command
   {
      mv ${ref_fasta} ${ref_dict} ${ref_fasta_index} ${ref_bwt} ${ref_sa} ${ref_amb} ${ref_ann} ${ref_pac} .
      ${gatk_launch} --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap -XX:MaxRAMFraction=1" \
        SelectVariants \
        -R ${ref_fasta_basename} \
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
      disks: "local-disk " + disk_size + " HDD"
      preempible: 4
   }

   output
   {
      File out_vcf = "${out_base_name}.vcf.gz"
   }
}