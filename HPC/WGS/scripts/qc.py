import sys;
import json;
from jinja2 import Template

REPORT_TEMPLATE = '''
QC REPORT FOR {{ sample_id }}


{% if autozygosity %}
--- Autozygosity ---

Percent autozygosity: {{ '{:.3%}'.format(autozygosity.fraction) }}
{% if autozygosity.regions %}
Autozygous regions:
CHR START   END KB  NSNPS
{% for region in autozygosity.regions %}
{{ region.chr }}    {{ region.start }}  {{ region.end }}    {{ region.kb }} {{ region.nsnps }}
{% endfor %}
{% endif %}


{% endif %}
--- Coverage ---

Mean coverage over exome: {{ coverage.mean | round(1) }}
Mean coverage over exome excluding low-coverage exons (DP < 5): {{ coverage.mean_high_cov | round(1) }}
Histogram of coverage by exon available at {{ coverage.histogram_filename }}


--- Sex check ---

{% if sex_check.passing %}
Sex check PASSED
{% else %}
Sex check FAILED

(This may be due to difference in capture kits, but further investigation may be warranted.)
{% endif %}


--- QC metrics ---

{% if qc_metrics.failed_metrics %}
The following QC metrics are outside of expected ranges:
{% for failed in qc_metrics.failed_metrics %}
 - {{ failed.metric }} (expected: {{ failed.expected }} +/- {{ failed.sd }}; observed: {{ failed.observed }}; Z-score: {{ failed.z }})
{% endfor %}
{% else %}
No QC metrics are outside of expected ranges
{% endif %}


--- Variant counts ---

Synonymous: {{ variant_counts.synonymous }}
Missense: {{ variant_counts.missense }}
Nonsense: {{ variant_counts.nonsense }}
Frameshift: {{ variant_counts.frameshift }}
Inframe: {{ variant_counts.inframe }}
Splice site: {{ variant_counts.splice }}


--- Clinical variants ---

{% if clinical_variants.variants %}
The following clinical variants are present:
{% for variant in clinical_variants.variants %}
 - {{ variant.chr }}:{{ variant.pos }} {{ variant.ref }}/{{ variant.alt }} ~ {{ variant.genotype }} ~ AF: {{ variant.frequency }} ~ {{ variant.gene }} ~ {{ variant.disease }} ~ {{ variant.inheritance }} ~ {{ variant.link }}
{% endfor %}
{% else %}
No clinical variants are present
{% endif %}


'''.lstrip()

def fill_dict(context, sample_list, task_name, listOfFiles):
  for json_file in listOfFiles.split(","):
    for sample in sample_list:
      context[sample][task_name] = json.load(open(json_file));

def single_fill_dict(context, sample_list, task_name, file):
  with open(file) as f:
    task_context = json.load(f);
    for sample in sample_list:
      context[sample][task_name] = task_context[sample];

coverage_jsons = sys.argv[1];
autozigosity_json = sys.argv[2];
sex_check_jsons = sys.argv[3];
qc_metrics_json = sys.argv[4];
variant_counts_jsons = sys.argv[5];
clinical_variants_jsons = sys.argv[6];

sample_list = sys.argv[7].split(",");

context = {sample: dict() for sample in sample_list}

fill_dict(context, sample_list, "coverage", coverage_jsons);
fill_dict(context, sample_list, "sex_check", sex_check_jsons);
fill_dict(context, sample_list, "variant_counts", variant_counts_jsons);
fill_dict(context, sample_list, "clinical_variants", clinical_variants_jsons);

single_fill_dict(context, sample_list, "autozygosity", autozigosity_json);
single_fill_dict(context, sample_list, "qc_metrics", qc_metrics_json);

template = Template(REPORT_TEMPLATE, trim_blocks=True)
for sid in sample_list:
  context[sid]['sample_id'] = sid
  file_name = "report." + sid + ".txt";
  with open(file_name, 'w') as f:
    f.write(template.render(context[sid]).encode('utf-8'))
