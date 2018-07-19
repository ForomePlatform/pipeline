import "wes_upstream_1SampleFastqs_gatk4.wdl" as sample_wf

workflow wes_upstream 
{
    #RESOURCES SECTION 
    File dbSNP_vcf
   
    File dbSNP_vcf_idx
    String known_indels_sites_VCF
    String known_indels_sites_idx
    String ref_fasta
	
    File ref_dict
    File scattered_calling_intervals_list

    #INPUT SECTION
    Map[String, Pair[File, File]] input_fastqs		#Map that links a sample to the pair of FASTQ files
    
    String tools					#Path to the directory that contains all needed tools and scripts
    String base_name					#Base name that will be used in every output file. 
    String res_dir				        #Path to the directory that will contain the results of the upstream run
    
    Int    bwa_threads								
    Int    samtools_threads

    #For every sample call the upstream sub-workflow
    scatter(key_value in input_fastqs)
    {
       String key_sampleName         = key_value.left
       Pair[File, File] value_fastqs = key_value.right
       
       call sample_wf.wes { 
                             input: 
                               sampleName                       = key_sampleName,
                               fastq1                           = value_fastqs.left,
                               fastq2                           = value_fastqs.right,
                               bwa_threads                      = bwa_threads,
                               samtools_threads                 = samtools_threads,
                               tools                            = tools,
                               base_name                        = base_name,
                               res_dir                          = res_dir,
                               dbSNP_vcf                        = dbSNP_vcf,
                               dbSNP_vcf_idx                    = dbSNP_vcf_idx,
                               known_indels_sites_VCF           = known_indels_sites_VCF,
                               known_indels_sites_idx           = known_indels_sites_idx,
                               ref_fasta                        = ref_fasta,
                               ref_dict                         = ref_dict, 
                               scattered_calling_intervals_list = scattered_calling_intervals_list
                          }
    }
}

