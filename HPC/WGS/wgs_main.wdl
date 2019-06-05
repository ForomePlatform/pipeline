import "wgs_upstream_1Sample.wdl" as sample_wf
import "wgs_downstream.wdl" as downstream
import "main_qc_report.wdl" as qc_rep

workflow wgs
{
    String run_id

    String gatk_version = "4.1.2.0"
    String gatk_version_old = "4.1.2.0"
    
    String directory_to_search 
    String tmp_dir 
   
    File chr_intervals_file
    
    Array[String] recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0"]
    Array[String] recalibration_annotation_values = ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]
    Array[String] SNPsVariantRecalibrator_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.0"]
    Array[String] SNPsVariantRecalibrator_recalibration_annotation_values = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]

    #RESOURCES SECTION
    File dbSNP_vcf
    File dbSNP_vcf_idx
    #
    File known_indels_sites_VCF
    File known_indels_sites_idx
    #
    File hapmap_vcf
    File hapmap_vcf_idx
    #
    File omni_vcf
    File omni_vcf_idx
    #
    File one_thousand_genomes_vcf
    File one_thousand_genomes_vcf_idx
    
    #PATHS FOR VEP
    String cache_path
    String exac_path 
    String plugin_path 
    
    File ref_fasta
    File ref_dict
    File ref_fasta_index
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac
    
    File case_fam                                               #File in the following format - https://www.cog-genomics.org/plink2/formats#fam
    File call_gathering_py
    
    #INPUT SECTION
    Map[String, Array[Pair[File, File]]] input_fastqs           #Map that links a sample to the pair of FASTQ files

    String tools                                                #Path to the directory that contains all needed tools and scripts
    String base_name                                            #Base name that will be used in every output file.
    String gvcfs_dir                                            #Path to the directory that will contain the gvcfs of the upstream run

    Int    bwa_threads
    Int    samtools_threads

    #QC REPORT INPUTS
    File bedtools_bin
    File clinical_bed
    File genome_coverage
    File exons
    File plink_bin
    File exome_ranges
    File vcf2geno
    File trace
    String hgdp_geno                                            #Should be File in the future
    String hgdp_coord                                           #Should be File in the future
    String hgdp_coords                                          #Should be File in the future
    File make_PCA_plot_script
    File report_script
    String gatk_launch_3_8

    String novo_dir_to_search                                   #root dir to search for .bam files
    String v_env_path_activation                                #python v_env for stage2 novo_caller
    File stub_file                                              #stub_file for select_first function
    
    String gvcf_naming_pattern                                  #e.g. "g.vcf.gz"

    call GenerateCaseIdsFile
    {
       input:
        fam = case_fam
    }

    call GetSampleIDs
    {
       input:
        fam = case_fam
    }

    Int family_size = length(GetSampleIDs.sample_ids)

    #For every sample call the upstream sub-workflow
    scatter(key_value in input_fastqs)
    {
       String key_sampleName                = key_value.left
       Array[Pair[File, File]] value_fastqs = key_value.right

       call sample_wf.wgs_upstream {
                                      input:
                                        sampleName                       = key_sampleName,
                                        fastqs                           = value_fastqs,
                                        bwa_threads                      = bwa_threads,
                                        samtools_threads                 = samtools_threads,
                                        base_name                        = base_name,
                                        dbSNP_vcf                        = dbSNP_vcf,
                                        dbSNP_vcf_idx                    = dbSNP_vcf_idx,
                                        known_indels_sites_VCF           = known_indels_sites_VCF,
                                        known_indels_sites_idx           = known_indels_sites_idx,
                                        ref_fasta                        = ref_fasta,
                                        ref_dict                         = ref_dict,
                                        ref_fasta_index                  = ref_fasta_index,
                                        ref_bwt                          = ref_bwt,
                                        ref_sa                           = ref_sa,
                                        ref_amb                          = ref_amb,
                                        ref_ann                          = ref_ann,
                                        ref_pac                          = ref_pac,
                                        tmp_dir                          = tmp_dir,
                                        tools                            = tools,
                                        gvcfs_dir                        = gvcfs_dir + "/" + base_name + key_sampleName,
                                        gatk_version                     = gatk_version
                                   }

       Pair[String, Pair[File, File]] res_sample_bam = wgs_upstream.sample_bam
       String res_sample = base_name + res_sample_bam.left
       File res_bam = res_sample_bam.right.left
       File res_bai = res_sample_bam.right.right
    }

    Array[String] res_samples = res_sample
    Array[File] res_bams = res_bam
    Array[File] res_bais = res_bai

    call GetTrioBams
    {
       input:
         fam = case_fam,
         bams = wgs_upstream.family_bam,
         bam_idx = wgs_upstream.family_bam_idx
    }
    
    call downstream.wgs_downstream {
                                      input:
                                        run_id                                                  = run_id,
                                        family_bams                                             = GetTrioBams.family_trio_bams,
                                        family_bams_idx                                         = GetTrioBams.family_trio_bams_idx,
                                        novo_dir_to_search                                      = novo_dir_to_search,
                                        family_size                                             = family_size,
                                        gatk_version                                            = gatk_version,
                                        gatk_version_old                                        = gatk_version_old,
                                        ref_fasta                                               = ref_fasta,
                                        ref_dict                                                = ref_dict,
                                        ref_fasta_index                                         = ref_fasta_index,
                                        ref_bwt                                                 = ref_bwt,
                                        ref_sa                                                  = ref_sa,
                                        ref_amb                                                 = ref_amb,
                                        ref_ann                                                 = ref_ann,
                                        ref_pac                                                 = ref_pac,
                                        tools                                                   = tools,
                                        directory_to_search                                     = directory_to_search, 
                                        tmp_dir                                                 = tmp_dir,
                                        chr_intervals_file                                      = chr_intervals_file,
                                        recalibration_tranche_values                            = recalibration_tranche_values,
                                        recalibration_annotation_values                         = recalibration_annotation_values,
                                        SNPsVariantRecalibrator_recalibration_tranche_values    = SNPsVariantRecalibrator_recalibration_tranche_values,
                                        SNPsVariantRecalibrator_recalibration_annotation_values = SNPsVariantRecalibrator_recalibration_annotation_values,
                                        dbsnp_vcf                                               = dbSNP_vcf,
                                        dbsnp_vcf_index                                         = dbSNP_vcf_idx,
                                        mills_vcf                                               = known_indels_sites_VCF,
                                        mills_vcf_index                                         = known_indels_sites_idx,
                                        hapmap_vcf                                              = hapmap_vcf,
                                        hapmap_vcf_idx                                          = hapmap_vcf_idx,
                                        omni_vcf                                                = omni_vcf,
                                        omni_vcf_idx                                            = omni_vcf_idx,
                                        one_thousand_genomes_vcf                                = one_thousand_genomes_vcf,
                                        one_thousand_genomes_vcf_idx                            = one_thousand_genomes_vcf_idx,
                                        case_sample_names                                       = GetSampleIDs.sample_ids,
                                        call_gathering_py                                       = call_gathering_py,
                                        cache_path                                              = cache_path,
                                        exac_path                                               = exac_path,
                                        plugin_path                                             = plugin_path,
                                        fam                                                     = case_fam,
                                        case_ids_file                                           = GenerateCaseIdsFile.case_ids,
                                        v_env_path_activation                                   = v_env_path_activation,
                                        stub_file                                               = stub_file,
                                        gvcf_naming_pattern                                     = gvcf_naming_pattern
                                   }

    call qc_rep.qc_report          {
                                      input:
                                        ref_fasta                                               = ref_fasta,
                                        ref_dict                                                = ref_dict,
                                        ref_fasta_index                                         = ref_fasta_index,
                                        ref_bwt                                                 = ref_bwt,
                                        ref_sa                                                  = ref_sa,
                                        ref_amb                                                 = ref_amb,
                                        ref_ann                                                 = ref_ann,
                                        ref_pac                                                 = ref_pac,
                                        tools                                                   = tools,
                                        fam                                                     = case_fam,
                                        run_id                                                  = run_id,
                                        bedtools_bin                                            = bedtools_bin,
                                        clinical_bed                                            = clinical_bed,
                                        genome_coverage                                         = genome_coverage,
                                        exons                                                   = exons,
                                        plink_bin                                               = plink_bin,
                                        exome_ranges                                            = exome_ranges,
                                        vcf2geno                                                = vcf2geno,
                                        trace                                                   = trace,
                                        hgdp_geno                                               = hgdp_geno,
                                        hgdp_coord                                              = hgdp_coord,
                                        hgdp_coords                                             = hgdp_coords,
                                        make_PCA_plot_script                                    = make_PCA_plot_script,
                                        report_script                                           = report_script,
                                        dbSNP_vcf                                               = dbSNP_vcf,
                                        dbSNP_vcf_idx                                           = dbSNP_vcf_idx,
                                        in_vcf                                                  = wgs_downstream.Vcf_for_QC,
                                        in_vcf_index                                            = wgs_downstream.VcfIdx_for_QC,
                                        prefix                                                  = "au",
                                        gatk_version                                            = gatk_version,
                                        gatk_launch_3_8                                         = gatk_launch_3_8,
                                        sample_ids                                              = res_samples,
                                        sample_bams                                             = res_bams,
                                        sample_idxs                                             = res_bais
                                   }
                                        
    output
    {
       Array[File] case_bams           = res_bams
       Array[File] case_bais           = res_bais
       File Vcf_xBrowse                = wgs_downstream.Vcf_for_xBrowse
       File Vcf_xBrowse_idx            = wgs_downstream.VcfIdx_for_xBrowse
       File pca_plot                   = qc_report.pcaPlot
       Array[File] reports             = qc_report.sample_reports
       Array[File] coverage_histograms = qc_report.coverage_histograms 
    }
}

task GetSampleIDs
{
  File fam

  command
  <<<
      python <<CODE
       with open("${fam}") as f:
        lines = f.readlines();
        for line in lines:
          if not line.startswith('#'):
            print line.split()[1];
      CODE
  >>>

  runtime
   {
      cpu: 1
      memory: "100 MB"
   }

   output
   {
     Array[String] sample_ids = read_lines(stdout())
   }
}

task GetTrioBams
{
   File fam
   Array[File] bams
   Array[File] bam_idx

   command
   <<<
      python <<CODE
      list_of_bams = ["${sep='", "' bams}"];
      list_of_idx  = ["${sep='", "' bam_idx}"];
      fam_lines = [];
      with open("${fam}") as f:
        lines = f.readlines();
        for line in lines:
          if not line.startswith('#'):
            fam_lines.append(line.split());

      for fam_line in fam_lines:
        if (int(fam_line[-1]) == 2 and '1' in fam_line[1]):
          p_id1 = fam_line[2];
          p_id2 = fam_line[3];
          a_id  = fam_line[1];

      with open("trio.bams.txt", 'w+') as f:
        for bam in list_of_bams:
          bam_name_split = bam.split('/')[-1].split('.');
          bam_sample_name = bam_name_split[0] + bam_name_split[1];
          if p_id1 in bam_sample_name or p_id2 in bam_sample_name or a_id in bam_sample_name:
            f.write(bam);
            f.write("\n");

      with open("trio.bams.idx.txt", 'w+') as f:
        for idx in list_of_idx:
          idx_name_split = idx.split('/')[-1].split('.');
          idx_sample_name = idx_name_split[0] + idx_name_split[1];
          if p_id1 in idx_sample_name or p_id2 in idx_sample_name or a_id in idx_sample_name:
            f.write(idx);
            f.write("\n");
      CODE
   >>>

   runtime
   {
      cpu: 1
      memory: "100 MB"
   }

   output
   {
      Array[File] family_trio_bams = read_lines("trio.bams.txt")
      Array[File] family_trio_bams_idx = read_lines("trio.bams.idx.txt")
   }
}

task GenerateCaseIdsFile
{
   File fam

   command
   <<<
      python <<CODE
      fam_lines = [];
      with open("${fam}") as f:
        lines = f.readlines();
        for line in lines:
          if not line.startswith('#'):
            fam_lines.append(line.split());

      for fam_line in fam_lines:
        if (int(fam_line[-1]) == 2 and '1' in fam_line[1]):
          p_id1 = fam_line[2];
          p_id2 = fam_line[3];
          a_id  = fam_line[1];

      if(not p_id1 == '0' and not p_id2 == '0'): 
        with open("case.ids.txt", 'w+') as f:
          f.write(p_id1 + '\t' + p_id2 + '\t' + a_id);
      else:
        with open("case.ids.txt", 'w+') as f:
          f.write('-1'); 
      CODE
   >>>

   runtime
   {
      cpu: 1
      memory: "100 MB"
   }

   output
   {
      File case_ids = "case.ids.txt"
   }
}
