workflow qc_report
{
   String run_id

   File ref_fasta
   File ref_dict
   File ref_fasta_index
   File ref_bwt
   File ref_sa
   File ref_amb
   File ref_ann
   File ref_pac

   File fam

   String gatk_version

   String tools
   File bedtools_bin
   File clinical_bed
   File genome_coverage
   File exons
   File plink_bin
   File exome_ranges
   File vcf2geno
   File trace
   String hgdp_geno
   String hgdp_coord
   String hgdp_coords
   File make_PCA_plot_script
   File report_script

   File dbSNP_vcf
   File dbSNP_vcf_idx

   File in_vcf
   File in_vcf_index

   String prefix

   #String gatk_version
   String gatk_launch = tools + "/gatk-" + gatk_version + "/gatk"
   String gatk_launch_3_8

   Array[String] sample_ids
   Array[File] sample_bams
   Array[File] sample_idxs

   scatter(index in range(length(sample_ids)))
   {
      String sample = sample_ids[index]
      File bam = sample_bams[index]
      File bai = sample_idxs[index]

      call SplitSamples
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
           vcf = in_vcf,
           vcf_ind = in_vcf_index,
           out_vcf = sample + ".vcf",
           gatk_launch = gatk_launch,
           sample_name = sample     
      }

      call VariantCounts
      {
         input:
           vcf = SplitSamples.sample_vcf,
           out_base_name = sample + ".variant.counts"  
      }

      call ClinicalVariants_Bedtools
      {
         input:
           bedtools_bin = bedtools_bin,
           clinical_bed = clinical_bed,
           vcf = SplitSamples.sample_vcf,
           genome_coverage = genome_coverage,
           out_name = sample + ".intersect"
      }

      call ClinicalVariants_Python
      {
         input:
           intersect = ClinicalVariants_Bedtools.out,
           out_name = sample + ".clinical.variants.dat.json"
      }
       
      #??
      call BedtoolsCoverage
      {
         input:
           bam = bam,
           bai = bai,
           exons = exons,
           coverage_genome = genome_coverage,
           bed_tools = bedtools_bin,
           out_file = "coverage.data.txt"
      }

      call CoveragePython
      {
         input:
           coverage_file = BedtoolsCoverage.coverage,
           sample_name = sample
      }

      call SexCheck
      {
        input:
          sample_name = sample,
          coverage_file = BedtoolsCoverage.coverage,
          fam = fam
      }
   }

   call Autozygosity
   {
      input:
        vcf = in_vcf,
        vcf_idx = in_vcf_index,
        plink_bin = plink_bin,
        samples = sample_ids,
        prefix = prefix,
   }

   call ProcessAutozygosity
   {
      input:
        out_name = "data_autozygosity.json",
        hom_indiv = Autozygosity.hom_indiv,
        hom = Autozygosity.hom
   }

   call VariantEval
   {
      input:
        gatk_launch = gatk_launch_3_8,
        ref_fasta = ref_fasta,
        ref_dict = ref_dict,
        ref_fasta_index = ref_fasta_index,
        ref_bwt = ref_bwt,
        ref_sa = ref_sa,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        dbsnp = dbSNP_vcf,
        vcf = in_vcf,
        vcf_idx = in_vcf_index,
        nt = 16,
        memory = "1024 MB",
        out_name = run_id + ".variant.eval"
   }

   call SplitCountVariantsEval
   {
      input:
        eval = VariantEval.eval
   }

   call SplitTiTvEval
   {
      input:
        eval = VariantEval.eval
   }

   call QcMetrics
   {
      input:
        ti_tv_split = SplitTiTvEval.split,
        count_variants_split = SplitCountVariantsEval.split,
        sample_ids = sample_ids
   }

   call VcfToGeno
   {
      input:
        vcf = in_vcf,
        vcf_idx = in_vcf_index,
        ranges = exome_ranges,
        sample_ids = sample_ids,
        vcf2geno = vcf2geno,
        prefix = run_id
   }

   call PcaTrace
   {
      input:
        trace = trace,
        geno = VcfToGeno.geno,
        site = VcfToGeno.site,
        hgdp_geno = hgdp_geno,
        hgdp_coord = hgdp_coord,
        out_base_name = run_id
   }

   call PcaPlot
   {
      input:
        hgdp_coords = hgdp_coords,
        coord = PcaTrace.coord,
        make_PCA_plot = make_PCA_plot_script,
        out_base_name = run_id + ".pca.plot.png"
   }

   call MakeReport
   {
      input:
        script = report_script,
        samples = sample_ids,
        coverage_files = CoveragePython.json,
        sex_check_files = SexCheck.json,
        variant_counts_files = VariantCounts.data,
        clinical_variants_files = ClinicalVariants_Python.variants,
        qc_metrics_file = QcMetrics.qc_metrics,
        autozygosity_file = ProcessAutozygosity.json 
   }

   call CoverageComparison
   {
      input:
        qc_reports = MakeReport.reports,
        samples = sample_ids,
        coverage_jsons = CoveragePython.json,
   }

   output
   {
      File pcaPlot = PcaPlot.plot
      Array[File] coverage_histograms = CoveragePython.hist
      Array[File] sample_reports = MakeReport.reports
   }
}

task MakeReport
{
   Array[File]   coverage_files
   Array[File]   sex_check_files
   Array[File]   variant_counts_files
   Array[File]   clinical_variants_files
   File          qc_metrics_file
   File          autozygosity_file
   Array[String] samples

   File script

   command
   <<<
     python ${script} ${sep=',' coverage_files} ${autozygosity_file} ${sep=',' sex_check_files} ${qc_metrics_file} ${sep=',' variant_counts_files} ${sep=',' clinical_variants_files} \
            ${sep=',' samples} 
   >>>

   runtime
   {
     cpu: 1
     memory: "512 MB"
   }

   output
   {
      Array[File] reports = glob("report.*.txt")
   }
}

task SplitSamples
{
   File ref_fasta
   File ref_dict
   File ref_fasta_index
   File ref_bwt
   File ref_sa
   File ref_amb
   File ref_ann
   File ref_pac
   File vcf
   File vcf_ind
   String out_vcf
   String gatk_launch
   String sample_name

   command
   {
      ${gatk_launch} --java-options "-XX:ParallelGCThreads=1" SelectVariants \
      -R ${ref_fasta} \
      -V ${vcf} \
      -sn ${sample_name} \
      --exclude-non-variants true \
      --remove-unused-alternates true \
      -O ${out_vcf}
   }

   runtime
   {
      cpu: 1
      memory: "2 GB"
   }

   output
   {
      File sample_vcf = "${out_vcf}"
   }
}

task BedtoolsCoverage
{
   File bam
   File bai
   File exons
   File coverage_genome
   String bed_tools
   String out_file

   command
   {
      ${bed_tools} coverage -b ${bam} -a ${exons} -sorted -g ${coverage_genome} -mean > ${out_file}
   }

   runtime
   {
      cpu: 1
      memory: "4 GB"
   }

   output
   {
      File coverage = "${out_file}"
   }
}

task CoveragePython
{
   File coverage_file
   String sample_name

   command
   <<<
     python << CODE
     import pandas
     import json
     import matplotlib;
     matplotlib.use('Agg');
     import matplotlib.pyplot as plt;

     def calculate_mean(df):
       return (df.length * df.coverage).sum() / df.length.sum();

     dataframe = pandas.read_csv("${coverage_file}", sep='\t', header=None, names=['chr', 'start', 'end', 'exon', 'coverage']);
     dataframe['length'] = dataframe.end - dataframe.start;
   
     mean = calculate_mean(dataframe);
     mean_high_cov = calculate_mean(dataframe[dataframe.coverage >= 5]);

     fig, ax = plt.subplots(figsize=(12, 9));
     dataframe.hist('coverage', range=(0, 300), bins=50, ax=ax);
     ax.set_title('Coverage of ${sample_name}');
     ax.set_xlabel('Exon depth');
     ax.set_ylabel('Frequency');

     fig.savefig("coverage.histogram.${sample_name}.png");

     f = open("mean.txt", "w");
     f.write(str(mean));
     f.close();

     f = open("mean_high_cov.txt", "w");
     f.write(str(mean_high_cov));
     f.close();

     with open("data_coverage.${sample_name}.json", 'w') as f:
       json.dump({'mean': mean, 'mean_high_cov': mean_high_cov, 'histogram_filename': "coverage.histogram.${sample_name}.png"}, f)

     CODE
   >>>

   runtime
   {
      cpu: 1
      memory: "1024 MB"
   }

   output
   {
      File hist = "coverage.histogram.${sample_name}.png"
      Float mean = read_float("mean.txt")
      Float mean_high_cov = read_float("mean_high_cov.txt")
      File json = "data_coverage.${sample_name}.json"
   }
}

task SexCheck
{
   File fam
   String sample_name
   File coverage_file
 
   command
   <<<
     python << CODE
     import pandas;
     import json;
     from collections import namedtuple;

     SAMPLE_TUPLE = namedtuple('Sample', ['sample_name', 'sex'])
     SEX_CODE_MAP = {'1': 'M', '2': 'F'}

     def parse_fam_file(fam_file):
       sample_data = {};
       with open(fam_file) as f:
         for l in f:
           if l.startswith('#') or not l:
             continue;
           fid, iid, pid, mid, sex_code, pheno = l.split()
           sex = SEX_CODE_MAP.get(sex_code)
           sample_data[iid] = sex;
       return sample_data
 
     sample_data = parse_fam_file("${fam}");

     reported_sex = sample_data["${sample_name}"];

     dataframe = pandas.read_csv("${coverage_file}", sep='\t', header=None, names=['chr', 'start', 'end', 'exon', 'coverage']);
     # exclude PAR and low-coverage exons
     high_cov_non_par_y = dataframe[
                (dataframe.coverage >= 5) &
                (dataframe.chr == 'chrY') &
                (dataframe.start > 2649520) &
                (dataframe.end < 59034050)
        ];

     exon_count = len(high_cov_non_par_y.index)
     male = (exon_count >= 200)  # empirically determined

     passing = ((male and reported_sex == 'M') or (not male and reported_sex == 'F'));

     with open("data_sex_check.${sample_name}.json", 'w') as f:
       json.dump({'passing': passing}, f);

     print passing
     CODE
   >>>

   runtime
   {
      cpu: 1
      memory: "512 MB"
   }

   output
   {
      String passing = read_string(stdout())
      File json = "data_sex_check.${sample_name}.json"
   }
}

task VariantEval
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
   File dbsnp
   File vcf
   File vcf_idx
   String memory
   Int nt

   String out_name

   command
   {
      java -XX:ParallelGCThreads=1 -jar ${gatk_launch} \
      -T VariantEval \
      -R ${ref_fasta} \
      -D ${dbsnp} \
      -noST -noEV -ST Sample -EV CountVariants -EV TiTvVariantEvaluator \
      --eval ${vcf} \
      -o ${out_name} \
      -nt ${nt}
   }

   runtime
   {
     cpu: nt
     memory: memory
   }

   output
   {
      File eval = "${out_name}"
   } 
}

task SplitCountVariantsEval
{
   File eval

   command
   {
      grep "^CountVariants" ${eval} > count_variants.split
   }

   runtime
   {
      cpu: 1
      memory: "512 MB"
   }

   output
   {
      File split = "count_variants.split"
   }
}

task SplitTiTvEval
{
   File eval

   command
   {
      grep "^TiTvVariantEvaluator" ${eval} > ti_tv.split
   }

   runtime
   {
      cpu: 1
      memory: "512 MB"
   }

   output
   {
      File split = "ti_tv.split"
   }
}

task QcMetrics
{
   File ti_tv_split
   File count_variants_split
   Array[String] sample_ids

   command
   <<<
     python << CODE
     import pandas;
     import json;

     sample_list = ['${sep="','" sample_ids}']
     
     full_df = pandas.read_table("${count_variants_split}", sep='\s+');
     titv_df = pandas.read_table("${ti_tv_split}", sep='\s+');

     full_df['titv'] = titv_df.tiTvRatio;
     full_metrics = list(full_df.columns.values)[5:];
     full_df = full_df[full_df.Sample != 'all'];
     full_df['unrelated'] = ~full_df.Sample.isin(sample_list);

     unrelated_data = {metric: (full_df[full_df.unrelated][metric].mean(),
                                   full_df[full_df.unrelated][metric].std())
                          for metric in full_metrics};

     output_data = dict()

     for sid in sample_list:
            all_metrics = list()
            failed_metrics = list()
            for metric in full_metrics:
                value = full_df[full_df.Sample == sid].iloc[0][metric]
                mean, sd = unrelated_data[metric]

                if sd == 0:
                    continue
                z = (value - mean) / sd

                data = {'metric': metric, 'observed': value, 'expected': mean,
                        'sd': sd, 'z': z}
                all_metrics.append(data)
                if abs(data['z']) >= 3:
                    failed_metrics.append(data)

            output_data[sid] = {'all_metrics': all_metrics,
                                'failed_metrics': failed_metrics}     
     
     with open("data_qc_metrics.json", 'w') as f:
       json.dump(output_data, f);
     CODE
   >>>

   runtime
   {
      cpu: 1
      memory: "512 MB"
   }

   output
   {
      File qc_metrics = "data_qc_metrics.json"
   }
}

task VcfToGeno
{
   File vcf
   File vcf_idx
   File ranges
   Array[String] sample_ids
   File vcf2geno
   String prefix

   command
   {
      ${vcf2geno} --inVcf ${vcf} \
                  --rangeFile ${ranges} \
                  --peopleIncludeID ${sep="," sample_ids} \
                  --thresholdGQ 20 \
                  --out ${prefix}
      sed -i 's/chr//g' ${prefix}.site
   }

   runtime
   {
      cpu: 1
      memory: "1024 MB"
   }

   output
   {
      File site = "${prefix}.site"
      File geno = "${prefix}.geno"
   }
}

task PcaTrace
{
   File trace
   File geno
   File site
   String hgdp_geno
   String hgdp_coord
   String out_base_name

   command
   {
      ${trace} -s ${geno} -g ${hgdp_geno} -c ${hgdp_coord} -k 10 -o ${out_base_name}
   }

   runtime
   {
      cpu: 1
      memory: "512 MB"
   }

   output
   {
      File coord = "${out_base_name}.ProPC.coord"
   }
}

task PcaPlot
{
   File hgdp_coords
   File coord
   File make_PCA_plot
   String out_base_name

   command
   {
      Rscript ${make_PCA_plot} ${coord} ${out_base_name} ${hgdp_coords}
   }

   runtime
   {
      cpu: 1
      memory: "256 MB"
   }

   output
   {
      File plot = "${out_base_name}"
   }
}

task VariantCounts
{
   File vcf
   String out_base_name 

   command
   <<<
     python << CODE
     import re;
     import json;

     PROTEIN_CODING = [('synonymous', ('synonymous_variant',)),
                  ('missense', ('missense_variant',)),
                  ('nonsense', ('stop_gained',)),
                  ('frameshift', ('frameshift_variant',)),
                  ('inframe', ('inframe_insertion', 'inframe_deletion')),
                  ('splice', ('splice_acceptor_variant',
                              'splice_donor_variant',
                              'splice_region_variant'))];

     canonical = re.compile(r'(?<==|,)[^,;\s]*\|YES\|[^,;\s]*(?=[;\s,])');

     data = {c_class: 0 for c_class, _ in PROTEIN_CODING};

     with open("${vcf}") as f:
       for l in f:
         if l.startswith('#'):
           continue;

         if l.split()[6] != 'PASS':
           continue;

         csqs = canonical.findall(l);
         for c_class, consequences in PROTEIN_CODING:
           if any(consequence in csq
                  for consequence in consequences
                  for csq in csqs):
             data[c_class] += 1;

     with open("${out_base_name}.dat.json", 'w') as f:
       json.dump(data, f);     
     
     CODE
   >>>

   runtime
   {
      cpu: 1
      memory: "512 MB"
   }  

   output
   {
      File data = "${out_base_name}.dat.json"
   }
}

task ClinicalVariants_Bedtools
{
   File bedtools_bin
   File clinical_bed
   File vcf
   File genome_coverage
   String out_name

   command
   {
      ${bedtools_bin} intersect -a ${clinical_bed} -b ${vcf} -sorted -g ${genome_coverage} -wa -wb -f 1 -F 1 > ${out_name}
   }

   runtime
   {
      cpu: 1
      memory: "512 MB"
   }

   output
   {
      File out = "${out_name}"
   }
}

task Autozygosity
{
   File vcf
   File vcf_idx
   File plink_bin
   String prefix
   Array[String] samples

   command
   <<<
     for sample in ${sep=' ' samples}  ; do
       echo -e "$sample\t$sample" >> keep.txt
     done

     ${plink_bin} --vcf ${vcf} --biallelic-only --double-id --vcf-filter --allow-extra-chr --autosome --snps-only \
     --keep keep.txt --homozyg -out ${prefix}
   >>>

   runtime
   {
      cpu: 1
      memory: "2048 MB"
   }

   output
   {
      File hom_indiv = "${prefix}.hom.indiv"
      File hom = "${prefix}.hom"
   }
}

task ProcessAutozygosity
{
   String out_name
   File hom_indiv
   File hom

   command
   <<<
     python << CODE
     from collections import namedtuple
     from itertools import groupby
     import json
     from operator import itemgetter
     import os
     import re
     import shlex
     from StringIO import StringIO
     import sys
     import subprocess
     import tempfile
     with open("${hom_indiv}") as f:
       f.readline()
       hom_indiv_data = {sid: next(line_iter) for sid, line_iter in groupby((l.split() for l in f), itemgetter(0))}  # only one line per sample
     with open("${hom}") as f:
       hom_data = {sid: list(line_iter) for sid, line_iter in groupby((l.split() for l in f), itemgetter(0))}

     save_data = dict()

     for sid in hom_indiv_data:
       data = {'fraction': float(hom_indiv_data[sid][4]) / 2881000,
               'regions': [{'chr': l[3],
               'start': int(l[6]),
               'end': int(l[7]),
               'kb': float(l[8]),
               'nsnps': int(l[9])} for l in hom_data.get(sid, list())]}

       save_data[sid] = data

     with open("${out_name}", 'w') as f:
       json.dump(save_data, f);
     CODE
   >>>

   runtime
   {
      cpu: 1
      memory: "512 MB"
   }

   output
   {
      File json = "${out_name}"
   }
}

task CoverageComparison
{
   Array[File] qc_reports
   Array[String] samples
   Array[File] coverage_jsons

   command
   <<<
     python << CODE
     import matplotlib
     matplotlib.use('Agg')
     import os
     import numpy as np
     import matplotlib.pyplot as plt
     import json

     def get_current_quality(file):
       with open(file) as f:
         json_data = json.load(f);
         contents = float(json_data['mean']);
       return contents

     from collections import namedtuple

     listOfReports = ["${sep='","' qc_reports}"];
     listOfSamples = ["${sep='","' samples}"];
     listOfJsons   = ["${sep='","' coverage_jsons}"];
     mean_qc_list=[];

     for report in listOfReports:
       with open(report) as qc_file:
         mean_qc_list.append(float(qc_file.read().split('Mean coverage over exome: ')[1].split('\nMean coverage')[0]));

     mean_qc_list.sort();

     mean_qc_wanted = [];
     mean_qc_dict = {};

     for index in range(len(listOfSamples)):
       temp_str = get_current_quality(listOfJsons[index]);
       mean_qc_wanted.append(temp_str);
       more = len([i for i in mean_qc_list if i <= temp_str]) / float(len(mean_qc_list));
       mean_qc_dict[temp_str]=[listOfSamples[index], more];

     plt.title("Quality Histogram")
     plt.xlabel("Coverage")
     plt.ylabel("Counts")
     num_bins = 20
     n, bins, patches = plt.hist(mean_qc_list, num_bins, normed=False, facecolor='blue', alpha=0.3)
     #display current sequence quality
     quality = mean_qc_wanted
     data = []
     for sample_mean in mean_qc_wanted:
       try:
         data.append(n[np.digitize(float(sample_mean), bins)])
       except:
         data.append(0)

     plt.scatter(quality, data, s=100, label='There are {} sequencing quality files'.format(len(mean_qc_list)) )
     plt.legend(loc='upper left')
     #label the plot with probe name
     for i,j in enumerate(quality):
       if mean_qc_dict[j][1]>0.05:
           plt.annotate(mean_qc_dict[j][0]+" qc is ok" ,xy=(j, data[i]))
       else:
           plt.annotate(mean_qc_dict[j][0]+" is in "+str(mean_qc_dict[j][1])+" worst", xy=(j, data[i]))

     plt.savefig('coverage_comparison.png')
     CODE
   >>>

   runtime
   {
      cpu: 1
      memory: "512 MB"
   }

   output
   {
      File plot = "coverage_comparison.png"
   }
}

task ClinicalVariants_Python
{
   File intersect
   String out_name

   command
   <<<
     python << CODE
     import json;
     import re;

     COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'};

     def matching_alleles(a, b):
       if any(not allele for allele in a + b):
          return False;
       acomp = [COMPLEMENTS[allele] for allele in a];
       return set(a) == set(b) or set(acomp) == set(b);

     def parse_frequency(info, alt):
       csq_string = re.search('CSQ=([^;\s]+)[;\s]', info);
       if not csq_string:
         return 0.;
       for csq in csq_string.group(1).split(','):
         fields = csq.split('|');
         print "Fields:"
         print len(fields);
         if fields[0] == alt and fields[65]:
           return float(fields[65]);
       return 0.;

     with open("${intersect}") as f:
       lines = f.readlines();

     variants = list();
     for variant in lines:
       fields = variant.strip().split('\t');
       hgvs_allele_1 = fields[3];
       hgvs_allele_2 = fields[4];
       disease = fields[5];
       gene = fields[6];
       link = fields[7];
       inheritance = fields[8];
       chr = fields[9];
       pos = fields[10];
       ref = fields[12];
       alt = fields[13];
       info = fields[16];
       gt = fields[18][:3];

       if ',' in alt:  # ignore multiallelic for now - very rare!
         continue;

       snp = (len(ref) == 1 and len(alt) == 1)
       if not snp:
         vep_alt = alt[1:];
         if not vep_alt:
           vep_alt = '-';
       else:
         vep_alt = alt;

       print alt, vep_alt;

       if snp and not matching_alleles((hgvs_allele_1, hgvs_allele_2), (ref, alt)):
         continue;

       frequency = parse_frequency(info, vep_alt);
       if frequency > 0.025:
         continue;

       genotype = 'het' if gt == '0/1' else 'hom alt';

       variants.append({'chr': chr,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'genotype': genotype,
                        'frequency': frequency,
                        'gene': gene,
                        'disease': disease,
                        'link': link,
                        'inheritance': inheritance}); 
     
     with open("${out_name}", 'w') as f:
       json.dump({'variants': variants}, f);
     CODE
   >>>

   runtime
   {
      cpu: 1
      memory: "256 MB"
   }

   output
   {
      File variants = "${out_name}"
   }
}
