# Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved
from __future__ import division
import sys
sys.path.insert(0, sys.argv[6]+"/code")
from string import Template
from proc_coverage import SamCoverage
from proc_ercc import load_ercc_conc, dose_response, generate_reports
from ercc_counter import write_output_counts_file
from jqplot_param_generator import chart_series_params, generate_trendline_points, generate_color_legend
from create_summary import create_summary_block

PLUGIN_NAME = sys.argv[4]
PLUGIN_DIR = sys.argv[6]
DATA = PLUGIN_DIR+'/data'
SRC = PLUGIN_DIR+'/code'
RAW_SAM_FILE = sys.argv[1] + '/tmap.sam'
FILTERED_SAM_FILE = sys.argv[1] + '/filtered.sam'
OUTPUT_DIR = sys.argv[2] + '/plugin_out/'+sys.argv[4]+'_out/'
COUNTS_FILE = OUTPUT_DIR + 'ercc.counts'
COUNTS_URL = sys.argv[3] + '/plugin_out/'+sys.argv[4]+'_out/ercc.counts'
GENOME   = DATA + '/ercc.genome'
ERCC_CONC = DATA + '/ERCCconcentration.txt'
ERCC_CONC_COLS = {'ERCC ID': str, 'Pool1': float, 'Pool2': float, 'log2pool1': float, 'log2pool2': float}
try:
  MINIMUM_RSQUARED = float(sys.argv[5])
except ValueError: 
  MINIMUM_RSQUARED = 0.9

coverage = SamCoverage(GENOME)

transcript_names = []
transcript_sizes = []
with open(GENOME) as ercc_genome:
    for line in ercc_genome:
        parsed_line = line.split()
        transcript_names.append(parsed_line[0])
        transcript_sizes.append(parsed_line[1])
counts, all_ercc_counts, total_counts, mean_mapqs = write_output_counts_file(RAW_SAM_FILE,FILTERED_SAM_FILE,COUNTS_FILE,transcript_names) # side-effect is to write counts file, filtered sam file
if (total_counts > 0):
  percent_total_counts_ercc = '%.2f' % (100 * (all_ercc_counts / total_counts))
  percent_total_counts_non_ercc = 100 - float(percent_total_counts_ercc)

coverage.parse_sam(FILTERED_SAM_FILE)

ercc_conc = load_ercc_conc(filter = counts.keys())
ercc_conc.sort()


dr = dose_response(coverage,ercc_conc,counts)
trendline_points = generate_trendline_points(dr)

data_to_display = True
msg_to_user = ""

try:
  report_components = generate_reports(coverage, dr, counts)
except ValueError:
  msg_to_user = "There does not appear to be any ERCC data to display.  Was the ERCC spike-in mix added to the sample?"
  data_to_display = False

if (all_ercc_counts < 250):
  msg_to_user = "The total number of counts (with acceptable mapping quality), "+str(all_ercc_counts)+", is not sufficient for a reliable correlation to be calculated."
  data_to_display = False
elif (dr[8]<3): 
  msg_to_user = "The number of transcripts detected, "+str(dr[8])+", is not sufficient for a reliable correlation to be calculated."
  data_to_display = False

if data_to_display:
  transcript_names_js, transcript_images_js, transcript_counts_js, transcript_sizes_js, series_params, ercc_points = chart_series_params(counts,transcript_sizes,transcript_names,  mean_mapqs,ercc_conc)
  color_legend = generate_color_legend()
  msg_to_user, rsquared = create_summary_block(OUTPUT_DIR,dr,MINIMUM_RSQUARED)
  template = open(SRC + '/ercc_template.html')
  page_we_are_making = Template(template.read())
  print page_we_are_making.substitute(dose_response_image = report_components['dose_response_image'],
                                    results_divs = report_components['results_divs'],
                                    ercc_points = ercc_points,
                                    percent_total_counts_non_ercc = percent_total_counts_non_ercc,
                                    percent_total_counts_ercc = percent_total_counts_ercc,
                                    all_ercc_counts = all_ercc_counts,
                                    rsquared = rsquared,
                                    slope = '%.2f' % (dr[2]),
                                    yintercept = '%.2f' % (dr[3]),
                                    N = str(dr[8]),
                                    trendline_points = trendline_points,
                                    color_legend = color_legend,
                                    counts_file = COUNTS_URL,
                                    msg_to_user = msg_to_user,
                                    series = series_params,
                                    plugin_name = PLUGIN_NAME )
else:
  template = open(SRC + '/ercc_error_template.html')
  page_we_are_making = Template(template.read())
  print page_we_are_making.substitute(msg_to_user = msg_to_user)
