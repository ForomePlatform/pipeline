import "wgs_upstream_1Sample.wdl" as sample_wf
import "wgs_downstream.wdl" as downstream

workflow wgs
{
    String run_id

    String gatk_version = "4.1.0.0"
    String gatk_version_old = "4.0.12.0"
    
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
    
    File case_fam
    File call_gathering_py
    
    File case_ids_file
    
    #INPUT SECTION
    Map[String, Array[Pair[File, File]]] input_fastqs           #Map that links a sample to the pair of FASTQ files

    String tools                                                #Path to the directory that contains all needed tools and scripts
    String base_name                                            #Base name that will be used in every output file.
    String gvcfs_dir                                            #Path to the directory that will contain the gvcfs of the upstream run

    Int    bwa_threads
    Int    samtools_threads

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
                                        tools                            = tools,
                                        gvcfs_dir                        = gvcfs_dir + "/" + base_name + key_sampleName
                                   }
    }
    
    call downstream.wgs_downstream {
                                      input:
                                        run_id                                                  = run_id,
                                        family_size                                             = length(read_lines(write_map(input_fastqs)))
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
                                        case_sample_names                                       = ["bgm9001u1", "bgm9001u2", "bgm9001a1"],
                                        call_gathering_py                                       = call_gathering_py,
                                        cache_path                                              = cache_path,
                                        exac_path                                               = exac_path,
                                        plugin_path                                             = plugin_path,
                                        fam                                                     = case_fam,
                                        case_ids_file                                           = case_ids_file,
                                   }

    #output
    #{
    #   Array[Array[File]] gvcfs = wes.gvcfs
    #   Array[Array[File]] gvcfs_idx = wes.gvcfs_idx
    #   Array[Pair[String, Pair[File, File]]] case_bams = wes.sample_bam
    #}
}
