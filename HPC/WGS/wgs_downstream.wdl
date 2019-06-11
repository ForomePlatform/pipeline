workflow wgs_downstream 
{
   String run_id

   String gatk_version
   String gatk_version_old
   File ref_fasta
   File ref_dict       
   File ref_fasta_index
   File ref_bwt         
   File ref_sa         
   File ref_amb         
   File ref_ann         
   File ref_pac        
   
   String directory_to_search
   String tmp_dir
   String tools
   
   File chr_intervals_file
   Array[String] chr_intervals = read_lines(chr_intervals_file)

   Array[String] chromosomes = ["chrM", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

   Array[String] recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0"]
   Array[String] recalibration_annotation_values = ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]

   Array[String] SNPsVariantRecalibrator_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.0"]
   Array[String] SNPsVariantRecalibrator_recalibration_annotation_values = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]

   Array[File] family_bams
   Array[File] family_bams_idx
   String novo_dir_to_search      

   File dbsnp_vcf
   File dbsnp_vcf_index

   File mills_vcf
   File mills_vcf_index

   File hapmap_vcf
   File hapmap_vcf_idx
   File omni_vcf
   File omni_vcf_idx
   File one_thousand_genomes_vcf
   File one_thousand_genomes_vcf_idx
   
   Array[String] case_sample_names
   
   File call_gathering_py
   
   String cache_path
   String exac_path
   String plugin_path
   
   File fam 
   File case_ids_file

   File stub_file #for select_first function to cast Type? to Type

   Int family_size

   String v_env_path_activation
   String gvcf_naming_pattern
   
   String vep_exe

   scatter (index in range(length(chr_intervals)))
   {
      Int i = index
      String chr = chr_intervals[i]
      String sub_interval = sub(chr, ":", "_")

      call CombineGVCFs
      {
         input:
           run_id = run_id,
           directory_to_search = directory_to_search,
           tmp_dir = tmp_dir,
           chr_to_get = sub(chr, ":.+", ""),
           ref_fasta = ref_fasta,
           ref_dict = ref_dict,
           ref_fasta_index = ref_fasta_index,
           ref_bwt = ref_bwt,
           ref_sa = ref_sa,
           ref_amb = ref_amb,
           ref_ann = ref_ann,
           ref_pac = ref_pac,
           output_base_name = sub(chr_intervals[i], ":.+", "") + "_" + i,
           tools = tools,
           gatk_version = gatk_version,
           cpu = 2,
           memory = "16000 MB",
           interval = chr_intervals[i],
           dbsnp = dbsnp_vcf,
           dbsnp_ind = dbsnp_vcf_index,
           gvcf_naming_pattern = gvcf_naming_pattern
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
           tools = tools,
           gatk_version = gatk_version,	
           cpu = 1,
           memory = "18 GB",	
           chr = chr_intervals[i],
           out_base_name = sub(chr_intervals[i], ":.+", "") + "_" + i,   		
           dbSNP_vcf = dbsnp_vcf,
           dbSNP_vcf_idx = dbsnp_vcf_index,
           gvcf = CombineGVCFs.gvcf,
           gvcf_idx = CombineGVCFs.gvcf_ind,
           tmp_dir = tmp_dir
      }
  }

  call MergeGenotypeRes
  {
     input:
       vcfs = GenotypeGVCFs.output_vcf,
       vcfs_ind = GenotypeGVCFs.output_vcf_idx,
       out_base_name = "genotyped",
       out_extension = "vcf.gz",
       tools = tools,
       num_threads = 2,
       memory = "4096 MB",
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
       dbsnp_resource_vcf = dbsnp_vcf,
       dbsnp_resource_vcf_index = dbsnp_vcf_index,
       mills_resource_vcf = mills_vcf,
       mills_resource_vcf_index = mills_vcf_index,
       tools = tools,
       gatk_version = gatk_version_old,	
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
       dbsnp_resource_vcf = dbsnp_vcf,
       dbsnp_resource_vcf_index = dbsnp_vcf_index,
       hapmap_resource_vcf = hapmap_vcf,
       hapmap_resource_vcf_index = hapmap_vcf_idx,
       omni_resource_vcf = omni_vcf,
       omni_resource_vcf_index = omni_vcf_idx,
       one_thousand_genomes_resource_vcf = one_thousand_genomes_vcf,
       one_thousand_genomes_resource_vcf_index = one_thousand_genomes_vcf_idx,
       tools = tools,
       gatk_version = gatk_version_old,
       cpu = 2,
       docker = "broadinstitute/gatk:latest",
       memory = "16 GB"
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
   
       tools = tools,
       gatk_version = gatk_version_old,
       docker = "broadinstitute/gatk:latest",
       cpu = 1,
       memory = "16 GB"   
  } 

  scatter (chr in chromosomes)
  {
     call SplitVcfByChr
     {
        input:
          in_vcf = ApplyRecalibration.recalibrated_vcf,
          in_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          chr = chr,
          tools = tools,
          memory = "2 GB",
          docker = "timuris/bwa"
     }

     call VepAnnotate
     {
        input:
          cache_path = cache_path,
          exac_path = exac_path,
          plugin_path = plugin_path,
          input_vcf = SplitVcfByChr.vcf,
          input_vcf_index = SplitVcfByChr.vcf_idx,
          output_basename = chr + ".vep",
          buffer = 50000,
          number_of_threads = 4,
          memory = "2 GB",
          vep_exe = vep_exe
     }

     if (family_size == 1)
     {
        call SelectCaseVariants as SelectSingleProbandVariants
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
             memory = "4096 MB",
             cpu = 1,
             tools = tools,
             gatk_version = gatk_version,
             in_vcf = VepAnnotate.out_vcf,
             in_vcf_idx = VepAnnotate.out_vcf_index,
             samples = case_sample_names,
             caller_names = ["BGM_DE_NOVO", "BGM_CMPD_HET", "BGM_HOM_REC", "BGM_AUTO_DOM", "BGM_BAYES_HOM_REC", "BGM_BAYES_DE_NOVO", "BGM_BAYES_CMPD_HET"],
             out_base_name = chr + ".xbrowse.vep"
        }
     }
     
     File SelectOutputVcf = select_first([SelectSingleProbandVariants.out_vcf, stub_file])
     File SelectOutputVcfIdx = select_first([SelectSingleProbandVariants.out_vcf_idx, stub_file])

     if (family_size > 1)
     {
        call BgmBayesDeNovo_stage1
        {
           input:
             case_ids  = case_ids_file,
             input_vcf = VepAnnotate.out_vcf, 
             x = 1.0,
             p = 0.005,
             e = 0.05,
             docker = "timuris/novo_caller:full",
             cpu = 1,
             memory = "2 GB",
             command = tools + "/denovo.out"
        }

        call ProcessBayesDeNovoSt1Res
        {
           input:
             res = BgmBayesDeNovo_stage1.stage_one_out,
             docker = "timuris/python2",
             memory = "1024 MB"
        }

        call BgmBayesDeNovo_stage2
        {
           input:
             st1_res = BgmBayesDeNovo_stage1.stage_one_out,
             family_bams = family_bams,
             family_bams_idx = family_bams_idx,
             dir_to_search = novo_dir_to_search,
             include_in_path = "wgs",
             exclude_from_path = run_id,
             cpu = 1,
             memory = "1024 MB", 
             chr = chr,
             script = tools + "/pysam_prac.py",
             v_env_path_activation = v_env_path_activation
        }

        call ProcessBayesDeNovoSt2Res
        {
           input:
             res = BgmBayesDeNovo_stage2.st2_res,
             docker = "timuris/python2",
             memory = "1024 MB"
        }

        call ABCaller as AutosomalDominantCaller
        {
           input:
             in_vcf = SplitVcfByChr.vcf,
             fam    = fam,
             id     = 0,
             output_name = "autosomalDominant.calls"
        }

        call ABCaller as HomozygousRecessiveCaller
        {
           input:
             in_vcf = SplitVcfByChr.vcf,
             fam    = fam,
             id     = 1,
             output_name = "homozygousRecessive.calls"
        }

        call ABCaller as DeNovoCaller
        {
           input:
             in_vcf      = SplitVcfByChr.vcf,
             fam         = fam,
             id          = 2,
             output_name = "deNovoCaller.calls"
        }

        call ABCaller as CompoundHeterozygousCaller
        {
           input:
             in_vcf      = SplitVcfByChr.vcf,
             fam         = fam,
             id          = 3,
             output_name = "compoundHeterozygous.calls"
        }


        call BgmCompoundHetCaller
        {
           input:
             case_ids  = case_ids_file,
             input_vcf = VepAnnotate.out_vcf,
             t         = 0.9,
             e         = 0.1,
             a         = 0.2,
             docker    = "timuris/novo_caller:full",
             cpu       = 1,
             memory    = "2 GB",
             command   = tools + "/compound_het.out"
        }

        call ProcessBgmCompHetRes
        {
           input:
             res = BgmCompoundHetCaller.out,
             docker = "timuris/python2",
             memory = "1024 MB",
             cpu = 1,
             vcf = VepAnnotate.out_vcf
        }

        call BgmHomRec
        {
           input:
             case_ids = case_ids_file,
             input_vcf = VepAnnotate.out_vcf,
             x = 1.0,
             p = 0.9,
             e = 0.05,
             a = 0.2,
             docker = "timuris/novo_caller:full",
             cpu = 1,
             memory = "2 GB",
             command = tools + "/homrec.out"
        }

        call ProcessHomRecRes
        {
           input:
             res = BgmHomRec.out,
             docker = "timuris/python2",
             memory = "1024 MB",
             cpu = 1
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
             bgmDeNovo_res = ProcessBayesDeNovoSt2Res.bayes_deNovo_calls_st2,
             vcf = ApplyRecalibration.recalibrated_vcf,
             vcf_index = ApplyRecalibration.recalibrated_vcf_index,
             docker = "timuris/python2",
             memory = "1024 MB",
             cpu = 1
        }

        call Bgzip
        {
           input:
             docker = "timuris/bwa",
             memory = "2 GB",
             cpu = 1,
             vcf = GatherCalls.vcf_out,
             tools = tools,
             out_base_name = "calls." + chr
        }

        call MergeVcfCalls
        {
           input:
             calls = Bgzip.vcf_out,
             calls_idx = Bgzip.vcf_out_index,
             vcf = VepAnnotate.out_vcf,
             vcf_index = VepAnnotate.out_vcf_index,
             out_base_name = chr + ".calls",
             docker = "timuris/bcftools:bgzip",
             tools = tools,
             memory = "2048 MB",
             cpu = 1,
        }
     }

     File MergeVcfCallsVcf = select_first([MergeVcfCalls.vcf_out, stub_file])
     File MergeVcfCallsVcfIdx = select_first([MergeVcfCalls.vcf_out_index, stub_file])
  }

  if (family_size == 1)
  {
     call BcfToolsMergeVcf as SingleProbandOutputMerge
     {
        input:
          vcfs = SelectOutputVcf,
          vcfs_ind = SelectOutputVcfIdx,
          out_base_name = "final",
          out_extension = "vcf.gz",
          out_type = "z",
          num_threads = 4,
          memory = "1 GB",
          docker = "timuris/bcftools:bgzip"
     }
  }

  if (family_size > 1)
  {
     call BcfToolsMergeVcf
     {
        input:
          vcfs = MergeVcfCallsVcf,
          vcfs_ind = MergeVcfCallsVcfIdx,
          out_base_name = "calls",
          out_extension = "vcf.gz",
          out_type = "z",
          num_threads = 4,
          memory = "1 GB",
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
          memory = "4096 MB",
          cpu = 1,
          tools = tools,
          gatk_version = gatk_version,
          in_vcf = BcfToolsMergeVcf.out_vcf,
          in_vcf_idx = BcfToolsMergeVcf.out_vcf_idx,
          #in_vcf = pre_final_vcf,
          samples = case_sample_names,
          caller_names = ["BGM_DE_NOVO", "BGM_CMPD_HET", "BGM_HOM_REC", "BGM_AUTO_DOM", "BGM_BAYES_HOM_REC", "BGM_BAYES_DE_NOVO", "BGM_BAYES_CMPD_HET"],
          out_base_name = "xbrowse.vep"
     }
  }

  output
  {
     File Vcf_for_QC         = select_first([BcfToolsMergeVcf.out_vcf, SingleProbandOutputMerge.out_vcf]) 
     File VcfIdx_for_QC      = select_first([BcfToolsMergeVcf.out_vcf_idx, SingleProbandOutputMerge.out_vcf_idx])
     File Vcf_for_xBrowse    = select_first([SelectCaseVariants.out_vcf, SingleProbandOutputMerge.out_vcf])
     File VcfIdx_for_xBrowse = select_first([SelectCaseVariants.out_vcf_idx, SingleProbandOutputMerge.out_vcf_idx])
  }
}

task DB_Import
{
  String directory_to_search
  String chr_to_get
  String tmp_dir
  
  String interval

  Int padding
 
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

  Int cpu
  String memory

  String db_path 

  command
  <<<
     find ${directory_to_search} -not -path '*/\.*' -not -path '*test*' -name "${chr_to_get}.gvcf.gz" -type f -exec echo -V {} \; > variants.txt && \
     value=`cat variants.txt`

     ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms2g -XX:ParallelGCThreads=1" GenomicsDBImport \
       --batch-size 50 \
       -ip ${padding} \
       --reader-threads 5 \
       --overwrite-existing-genomicsdb-workspace true \
       -L ${interval} \
       -R ${ref_fasta} \
       --tmp-dir ${tmp_dir} \
       --genomicsdb-vcf-buffer-size 32768 \
       --genomicsdb-workspace-path ${db_path} \
       $value

     tar -cf ${db_path}.tar ${db_path}
  >>>

  runtime
  {
     cpu: cpu
     memory: memory
  }

  output
  {
     File db = "${db_path}.tar"
  }

}

task SelectCaseVariants
{
   String tools
   String gatk_version
   
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
      ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms2g" \
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
      File out_vcf_idx = "${out_base_name}.vcf.gz.tbi"
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
   String docker
   String memory
   Int cpu

   command
   <<<
      python <<CODE
      import vcf as pyvcf;

      def add_headers(reader):
        info = pyvcf.parser._Info('BGM_AUTO_DOM', 1, 'Flag', 'Autosomal dominant by BGM allele balance caller', None, None);
        reader.infos[info.id] = info;

        info = pyvcf.parser._Info('BGM_HOM_REC', 1, 'Flag', 'Homozygous recessive by BGM allele balance caller', None, None);
        reader.infos[info.id] = info;

        info = pyvcf.parser._Info('BGM_DE_NOVO', 1, 'Flag', 'De novo by BGM allele balance caller', None, None);
        reader.infos[info.id] = info;

        info = pyvcf.parser._Info('BGM_CMPD_HET', 1, 'Integer', 'Compound heterozygous by BGM allele balance caller', None, None);
        reader.infos[info.id] = info;

        info = pyvcf.parser._Info('BGM_BAYES_CMPD_HET', 1, 'Float', 'Probability of compund het by BGM Bayes caller', None, None);
        reader.infos[info.id] = info;

        info = pyvcf.parser._Info('BGM_BAYES_DE_NOVO', 1, 'Float', 'Probability of de novo by BGM Bayes caller', None, None);
        reader.infos[info.id] = info;

        info = pyvcf.parser._Info('BGM_BAYES_HOM_REC', 1, 'Float', 'Probability of hom rec by BGM Bayes caller', None, None);
        reader.infos[info.id] = info;

      def process_site(site, vcf_writer, listOfCallers):
        site_split = site.split("[");
        alts = list(site_split[1].replace("]", "").split(", "));
        
        site_split = site_split[0].split(", ")[:-1];
        record = pyvcf.model._Record('', 0, None, '', [None], None, None, None, None, None);

        record.CHROM = site_split[0];
        record.POS = site_split[1];
        record.REF = site_split[2];
        record.ALT = [pyvcf.model._Substitution(alt) for alt in alts]
        record.INFO = INFO=listOfCallers;

        if "chr" not in record.CHROM:
          record.CHROM = "chr" + record.CHROM;

        vcf_writer.write_record(record)
      #################################
      deNovoCallerRes         = "${deNovo_res}";
      bgmHomRecRes            = "${bgmHomRec_res}";
      compoundHeterozygousRes = "${compoundHeterozygous_res}";
      autosomalDominantRes    = "${autosomalDominant_res}";
      bgmBayesDeNovoRes       = "${bgmDeNovo_res}";
      homRecRes               = "${homozygousRecessive_res}";
      bgmCompoundHetRes       = "${autosomalDominant_res}";
      vep_vcf                 = "${vcf}";

      vcf_reader = pyvcf.Reader(filename=vep_vcf);
      vcf_reader.infos = dict();
      vcf_reader.samples = list();

      add_headers(vcf_reader);
      calls_vcf = open('calls.vcf', 'w+');
      vcf_writer = pyvcf.Writer(calls_vcf, template=vcf_reader);

      site_callers = {};
      listOfResults = [bgmCompoundHetRes, homRecRes, compoundHeterozygousRes, bgmHomRecRes, deNovoCallerRes, autosomalDominantRes, bgmBayesDeNovoRes];

      for result in listOfResults:
        with open(result) as f:
          line = f.readline();
          while line:
            line_split = line[1:-2].replace("'", "").split("), ");
            site = line_split[0].replace("(", "");

            if not site in site_callers:
              site_callers[site] = {};

            callerName = line_split[1].split(", ")[-1].replace("'", "");
            callerValue = line_split[1].split(", ")[-2].replace("'", "");

            if callerValue == 'True':
              callerValue = True;

            site_callers[site][callerName] = callerValue;
            line = f.readline();

      for key in sorted(site_callers.iterkeys(), key=lambda k: int(k.split(", ")[1])):
        process_site(key, vcf_writer, site_callers[key]);
      CODE
   >>>

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

   String tools

   Int num_threads
   String memory
   String docker

   command
   {
      set -e
      ${tools}/bcftools concat --threads ${num_threads} -a -D -O z -o ${out_base_name}.${out_extension} ${sep=' ' vcfs}
      ${tools}/tabix -p vcf ${out_base_name}.${out_extension}
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

task Bgzip
{
  File vcf
  String out_base_name

  String docker
  String memory
  Int cpu  

  String tools

  command
  {
     ${tools}/bgzip -c ${vcf} > ${out_base_name}.vcf.gz 
     ${tools}/tabix -p vcf ${out_base_name}.vcf.gz
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
   String tools

   String docker
   String memory
   Int cpu

   command
   {
      set -e
      ${tools}/bcftools annotate -a ${calls} -c INFO ${vcf} > ${out_base_name}.vcf
      ${tools}/bgzip -c ${out_base_name}.vcf > ${out_base_name}.vcf.gz
      ${tools}/tabix -p vcf ${out_base_name}.vcf.gz      
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

task ProcessHomRecRes
{
   File res

   String docker
   String memory
   Int cpu

   command
   <<<
     python <<CODE
     import re;
     
     HOM_REC_CALL_PATT = re.compile(r'^(?P<CHROM>\S+)\s+(?P<POS>\d+)\s+(\S+)\s+(?P<REF>[ACGT]+)\s+(?P<ALT>[ACGT,]+).*?P_hom_rec=(?P<PROB>[^;\n]+)', re.M);
     res_file = open("hom_rec_calls.txt", "w+");
     
     with open("${res}", "r") as f:
       call_iter = HOM_REC_CALL_PATT.finditer(f.read());

     listOfCalledTuples = [];
     for call in call_iter:
       site = (call.group('CHROM'),
               int(call.group('POS')),
               call.group('REF'),
               call.group('ALT').split(','));
       prob = float(call.group('PROB');
       if prob > 0.5:
         call_tuple = (site, prob, 'BGM_BAYES_HOM_REC');
         listOfCalledTuples.append(call_tuple);
         
     for call in listOfCalledTuples:
       print >> res_file, call;
     res_file.close();
     CODE
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

   String docker
   String memory
   Int cpu

   command
   <<<
     python <<CODE
     import re;
     import sys;
     import vcf as pyvcf;
     import os;
     
     cmd = "zcat ${vcf} > decompressed.vcf";
     os.system(cmd);
     res_file = open("bgm_comp_het_calls.txt", "w+");
     
     CMPD_HET_CALL_PATT = re.compile(r'^Probty=(?P<PROB>\S+)\s+LOC=(?P<CHROM>\S+)\s+(?P<POS>\d+)', re.M);

     with open("${res}", 'r') as f:
       call_iter = CMPD_HET_CALL_PATT.finditer(f.read());

     reader = pyvcf.Reader(filename=sys.argv[2]);
     read_iter = reader.__iter__();

     listOfCalls = [];
     for call in call_iter:
       listOfCalls.append(call);

     listOfCalls.sort(key=lambda x: int(x.group('POS')));     
     lisOfCalledTuples = [];
     
     for call in listOfCalls:
       while True:
         record = read_iter.next();
         if int(record.POS) == int(call.group('POS')):
           break;

       site = ("chr" + call.group('CHROM'),
               int(call.group('POS')),
               record.REF,
               [str(alt) for alt in record.ALT]);
       prob = float(call.group('PROB'));
       if prob > 0.5:
         call_tuple = (site, prob, 'BGM_BAYES_CMPD_HET');
         lisOfCalledTuples.append(call_tuple);
     
     for call in listOfCalledTuples:
       print >> res_file, call;
     res_file.close();
     CODE
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
   String docker
   String memory

   command
   <<<
     python <<CODE
     import re;
     
     result_file = open("bayes_de_novo1_calls.txt", "w+");
     DE_NOVO_VCF_CALL_PATT = re.compile(r'^(?P<CHROM>\S+)\s+(?P<POS>\d+)\s+.+\s+(?P<REF>[ACGT]+)\s+(?P<ALT>[ACGT,]+).*?P_denovo=(?P<PP>\d\.\d+).+', re.M);
     with open("${res}", "r") as f:
       call_iter = DE_NOVO_VCF_CALL_PATT.finditer(f.read());
       
     list_of_tuples = [];
     for call in call_iter:
       site = (call.group('CHROM'),
               int(call.group('POS')),
               call.group('REF'),
               call.group('ALT').split(','));
       prob = float(call.group('PP'));       
       
       if prob >= 0.01:
         call_tuple = (site, prob, 'BGM_BAYES_DE_NOVO');
         list_of_tuples.append(call_tuple);
         
     for result in list_of_tuples:
       print >> result_file, result
     result_file.close();
     CODE
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

task ProcessBayesDeNovoSt2Res
{
   File res
   String docker
   String memory

   command
   <<<
     python <<CODE
     import re;
     
     result_file = open("bayes_de_novo2_calls.txt", "w+");
     DE_NOVO_CALL_PATT = re.compile(r'^\d+\)\s+(?P<CHROM>\S+)\s+(?P<POS>\d+)\s+(?P<REF>[ACGT]+)\s+(?P<ALT>[ACGT,]+).*?PP=(?P<PP>\S+)', re.M);

     list_of_tuples = [];
     with open("${res}", "r") as f:
       call_iter = DE_NOVO_CALL_PATT.finditer(f.read());

     for call in call_iter:
       site = (call.group('CHROM'),
               int(call.group('POS')),
               call.group('REF'),
               call.group('ALT').split(','));
       prob = float(call.group('PP'));

       if prob >= 0.01:
         call_tuple = (site, prob, 'BGM_BAYES_DE_NOVO');
         list_of_tuples.append(call_tuple);
         
     for result in list_of_tuples:
       print >> result_file, result
     result_file.close();
     CODE
   >>>

   runtime
   {
     memory: memory
     cpu: 1
     docker: docker
   }

   output
   {
     File bayes_deNovo_calls_st2 = "bayes_de_novo2_calls.txt"
   }
}

task BgmBayesDeNovo_stage2
{
   File st1_res
   Array[File] family_bams
   Array[File] family_bams_idx
   
   String dir_to_search
   String include_in_path
   String exclude_from_path
   Int cpu
   String memory
   
   String chr
   File script

   #optional
   String v_env_path_activation
  
   command
   <<<
      if [[ ! -s ${st1_res} ]]; then 
        echo "stage 1 res is empty"
        touch ${chr}.st2.res
        exit 0;
      fi;
      ${v_env_path_activation}
      echo -e ' ${sep='\n' family_bams}' > case.bams
      find ${dir_to_search} -type f -name "*bam" -path "${include_in_path}" -not -path "${exclude_from_path}" > allBams.txt
      python ${script} -I ${st1_res} -U allBams.txt -T case.bams -O ${chr}.st2.res
   >>>
   
   runtime
   {
      memory: memory
      cpu: cpu
   }
   
   output
   {
      File st2_res = "${chr}.st2.res"
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

     input="${case_ids}"
     while IFS= read -r var
     do
       echo "$var"
       if [ "$var" = "-1" ]; then
         echo "case_ids has -1. Not running DeNovo caller"
         touch bayes_de_novo_st1.out.vcf
         exit 0
       fi
     done < "$input"

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

  String command

  command
  {
     set -e
     zcat ${input_vcf} > decompressed.vcf
     ${command} -I decompressed.vcf -O bgm_cmp_het.out -TR ${case_ids} -T ${t} -E ${e} -A ${a}
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

   String command

   command
   {
      set -e 
      zcat ${input_vcf} > decompressed.vcf
      ${command} -I decompressed.vcf -O bgm_hom_rec.out -T ${case_ids} -X ${x} -P ${p} -E ${e} -A ${a}
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

task SplitVcfByChr
{
   File in_vcf
   File in_vcf_index
   String chr
   
   String tools
   
   String docker
   String memory

   command
   {
      ${tools}/tabix -h ${in_vcf} ${chr} > ${chr}.vcf
      ${tools}/bgzip -c ${chr}.vcf > ${chr}.vcf.gz
      ${tools}/tabix -p vcf ${chr}.vcf.gz
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

task VepAnnotate
{
   String cache_path
   String exac_path
   String plugin_path
   File input_vcf
   File input_vcf_index
   String output_basename
   Int number_of_threads
   Int buffer
   String memory
   String vep_exe
   
   command
   {
      ${vep_exe} --format vcf --force_overwrite --dir_cache ${cache_path} --dir_plugin ${plugin_path} --offline --no_stats --vcf --everything \
          --plugin ExAC,${exac_path} --allele_number --input_file ${input_vcf} --output_file ${output_basename}.vcf.gz --tabix --fork ${number_of_threads} --buffer_size ${buffer}
   }

   runtime
   {
     memory: memory
     cpu: number_of_threads
     docker: "timuris/vep"
   }

   output
   {
     File out_vcf = "${output_basename}.vcf.gz"
     File out_vcf_index = "${output_basename}.vcf.gz.tbi"
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

  String tools
  String gatk_version
  String docker
  Int cpu
  String memory

  command 
  {
    set -e
    ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms4g -XX:ParallelGCThreads=1" \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ${input_vcf} \
      --recal-file ${indels_recalibration} \
      --tranches-file ${indels_tranches} \
      --truth-sensitivity-filter-level ${indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms4g -XX:ParallelGCThreads=1" \
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

  String tools
  String gatk_version

  Int cpu
  String memory
  String docker

  command
  {
    ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms4g -XX:ParallelGCThreads=1" \
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

   String tools
   String gatk_version
   String docker
   String memory
   Int cpu

   command
   {
      ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms4g -XX:ParallelGCThreads=1" \
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

task ABCaller
{
   File in_vcf
   File fam
   Int id
   String output_name

   command
   {
      python <<CODE
      from __future__ import division
      from collections import namedtuple
      import vcf as pyvcf

      res_file = open("${output_name}", "w+");

      def all_gts_match(sample_gts, *gts):
          return all(gt in gts for gt in sample_gts)

      FAM_TUPLE = namedtuple('fam_structure', ['samples', 'affected', 'unaffected', 'de_novo'])

      def parse_fam_file(fam_file):

          samples = dict()

          for line in fam_file:
              if line.startswith('#'):
                  continue
              try:
                  _, id, fid, mid, _, pheno_code = line.split()
                  mid = None if mid in ('0', '-9')  else mid
                  fid = None if fid in ('0', '-9')  else fid
                  affected = (pheno_code == '2')
                  samples[id] = {'fid': fid,
                                 'mid': mid,
                                 'affected': affected,
                                 'in_vcf': True}
              except:
                  raiseException('Could not parse fam file line: {}'
                                  .format(line.strip()))

          all_mids = {s['mid'] for s in samples.values() if s['mid']}
          all_fids = {s['fid'] for s in samples.values() if s['fid']}
          unknown_ids = (all_mids | all_fids) - set(samples.keys())
          for id in unknown_ids:
              samples[id] = {'mid': None,
                             'fid': None,
                             'affected': False,
                             'in_vcf': False}

          affected_ids = [id for id, s in samples.items() if s['affected']]
          unaffected_ids = [id for id, s in samples.items()
                            if s['in_vcf'] and not s['affected']]
          naffected = len(affected_ids)
          if naffected <= 0:
              raiseException('At least one affected sample should be present')

          de_novo = None
          if naffected == 1:
              id, mid, fid = [(id, s['mid'], s['fid']) for id, s in samples.items()
                                 if s['affected']].pop()
              if (mid is not None and fid is not None and
                  samples[mid]['in_vcf'] and samples[fid]['in_vcf']):
                  de_novo = (id, mid, fid)

          return FAM_TUPLE(samples, affected_ids, unaffected_ids, de_novo)


      def make_ab_call(reads):
          """
          0-5% alt reads: 0/0
          5-85% alt reads: 0/1
          80-100% alt reads: 1/1
          """

          nref = reads[0]
          nalt = sum(reads[1:])

          if nref + nalt == 0:
              return None

          prop_alt = nalt / (nref + nalt)

          if prop_alt <= 0.05:
              return '0/0'
          if 0.05 < prop_alt <= 0.85:
              return '0/1'
          if 0.85 < prop_alt:
              return '1/1'


      def recall_genotypes(record):
          for call in record.samples:
              ad = getattr(call.data, 'AD', None)
              ad = map(int, ad);
              if not ad:
                  continue
              new_gt = make_ab_call(ad)

              if new_gt:
                  call.gt_alleles = new_gt.split('/')
                  call.gt_nums = new_gt
                  call.called = True
                  call.data = call.data._replace(GT=new_gt)


      class ABCaller(object):

          UNRELATED_AF_THRESHOLD = 0.1
          callerName = "";

          def __init__(self, vcf_file, fam_file, callerName = "Abstract"):
              self.reader = pyvcf.Reader(filename=vcf_file)
              self.callerName = callerName;
              with open(fam_file) as f:
                  self.fam_structure = parse_fam_file(f)
              self.unrelated_ids = set(self.reader.samples) - set(self.fam_structure.affected + self.fam_structure.unaffected)

          def case_gts(self, call_dict):
              affected_gts = [call_dict[id]
                              for id in self.fam_structure.affected]
              unaffected_gts = [call_dict[id]
                                for id in self.fam_structure.unaffected]
              return affected_gts, unaffected_gts

          def unrelated_af(self, call_dict):
              unrelated_gts = [call_dict[id]
                               for id in self.unrelated_ids
                               if call_dict[id] is not None]
              if not unrelated_gts:
                  return 0.
              return sum(unrelated_gts) / (2. * len(unrelated_gts))

          def __iter__(self):
              for record in self.reader:
                  recall_genotypes(record)
                  call_dict = {s.sample: s.gt_type for s in record.samples}

                  if self.unrelated_af(call_dict) > self.UNRELATED_AF_THRESHOLD:
                      continue
                  call = self.make_call(call_dict)
                  if not call:
                      continue
                  yield (
                      record.CHROM,
                      record.POS,
                      record.REF,
                      [str(alt) for alt in record.ALT]
                  ), call, self.callerName;

          def make_call(self, call_dict):
                  raise NotImplementedError


      class AutosomalDominantCaller(ABCaller):

          CALLER_NAME = 'BGM_AUTO_DOM'
          CALLER_NUMBER = 1
          CALLER_TYPE = 'Flag'
          CALLER_DESC = 'Autosomal dominant by BGM allele balance caller'

          def __init__(self, vcfFile, famFile):
              ABCaller.__init__(self, vcfFile, famFile, self.CALLER_NAME);

          def make_call(self, call_dict):
              affected_gts, unaffected_gts = self.case_gts(call_dict)
              return (all_gts_match(affected_gts, 1, 2) and
                      all_gts_match(unaffected_gts, 0, None))


      class HomozygousRecessiveCaller(ABCaller):

          CALLER_NAME = 'BGM_HOM_REC'
          CALLER_NUMBER = 1
          CALLER_TYPE = 'Flag'
          CALLER_DESC = 'Homozygous recessive by BGM allele balance caller'

          def __init__(self, vcfFile, famFile):
              ABCaller.__init__(self, vcfFile, famFile, self.CALLER_NAME);

          def make_call(self, call_dict):
              affected_gts, unaffected_gts = self.case_gts(call_dict)
              return (all_gts_match(affected_gts, 2) and
                      all_gts_match(unaffected_gts, 0, 1, None))


      class DeNovoCaller(ABCaller):

          CALLER_NAME = 'BGM_DE_NOVO'
          CALLER_NUMBER = 1
          CALLER_TYPE = 'Flag'
          CALLER_DESC = 'De novo by BGM allele balance caller'

          def __init__(self, vcfFile, famFile):
              ABCaller.__init__(self, vcfFile, famFile, self.CALLER_NAME);

          def make_call(self, call_dict):
              if not self.fam_structure.de_novo:
                  return False
              affected_gts, unaffected_gts = self.case_gts(call_dict)
              return (all_gts_match(affected_gts, 1, 2) and
                      all_gts_match(unaffected_gts, 0, None))


      class CompoundHeterozygousCaller(ABCaller):

          CALLER_NAME = 'BGM_CMPD_HET'
          CALLER_NUMBER = 1
          CALLER_TYPE = 'Integer'
          CALLER_DESC = 'Compound heterozygous by BGM allele balance caller'

          def __init__(self, vcfFile, famFile):
              ABCaller.__init__(self, vcfFile, famFile, self.CALLER_NAME);

          def make_call(self, call_dict):
              affected_gts, unaffected_gts = self.case_gts(call_dict)

              if not all_gts_match(affected_gts, 1, 2):
                  return False

              # bitwise flag to represent which unaffected samples have the variant
              # e.g. if the unaffected_gts are [0, 1, 0, 1] (the 2nd and 4th samples
              # have the variant), value = 1<<2 + 1<<4 = 2^1 + 2^4 = 18.
              # xBrowse will use bitwise AND operator (&) to select pairs of variants
              # after appropriate filtering such that no unaffected sample has both
              # variants (var1_value & var2_value == 0)
              value = sum(1 << i for i, gt in enumerate(unaffected_gts)
                          if gt in (1, 2))
              return value

      vcfFile  = "${in_vcf}";
      famFile  = "${fam}";
      callerId = int("${id}");

      listOfResults = [];

      if(callerId == 0):
        adCaller = AutosomalDominantCaller(vcfFile, famFile);
        for call in adCaller:
          listOfResults.append(call);

      elif(callerId == 1):
        hrCaller = HomozygousRecessiveCaller(vcfFile, famFile);
        for call in hrCaller:
          listOfResults.append(call);

      elif(callerId == 2):
        deNovoCaller = DeNovoCaller(vcfFile, famFile);
        for call in deNovoCaller:
          listOfResults.append(call);

      elif(callerId == 3):
        cmpHetCaller = CompoundHeterozygousCaller(vcfFile, famFile);
        for call in cmpHetCaller:
          listOfResults.append(call);
      
      for result in listOfResults:
        print >> res_file, result
      res_file.close();
      CODE
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

   String tools
   String gatk_version

   Int cpu
   String memory

   String chr
   String out_base_name

   File dbSNP_vcf
   File dbSNP_vcf_idx

   File gvcf
   File gvcf_idx
   String tmp_dir

   command
   {
      set -e
      ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms6g -XX:ParallelGCThreads=1" GenotypeGVCFs \
                     -R ${ref_fasta} \
                     -V ${gvcf} \
                     -D ${dbSNP_vcf} \
                     -O ${out_base_name}.raw.vcf.gz \
                     -L ${chr} \
                     -G StandardAnnotation \
                     --tmp-dir ${tmp_dir}
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
   }
}

task CombineGVCFs
{
   String run_id
   String directory_to_search
   String tmp_dir
   String chr_to_get
   File ref_fasta
   File ref_dict
   File ref_fasta_index
   File ref_bwt
   File ref_sa
   File ref_amb
   File ref_ann
   File ref_pac

   String output_base_name

   String tools
   String gatk_version

   Int cpu
   String memory

   String interval
   File dbsnp
   File dbsnp_ind
   
   String gvcf_naming_pattern

   command
   <<<
        echo ${run_id}
        find ${directory_to_search} -not -path '*/\.*' -not -path '*test*' -name "${chr_to_get}.${gvcf_naming_pattern}" -type f -exec echo -V {} \; > variants.txt && \
        value=`cat variants.txt`

        ${tools}/gatk-${gatk_version}/gatk --java-options "-Xms2g -XX:ParallelGCThreads=1" CombineGVCFs \
            -R ${ref_fasta} \
            -L ${interval} \
            -D ${dbsnp} \
            -O ${output_base_name}.g.vcf.gz \
            -ip 500 \
            --tmp-dir ${tmp_dir} \
            $value
   >>>

   output
   {
      File gvcf = "${output_base_name}.g.vcf.gz"
      File gvcf_ind = "${output_base_name}.g.vcf.gz.tbi"
   }

   runtime
   {
      memory: memory
      cpu: cpu
   }
}
