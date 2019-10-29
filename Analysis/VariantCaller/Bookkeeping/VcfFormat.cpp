/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     VcfFormat.cpp
//! @ingroup  VariantCaller
//! @brief    Vcf file formatting & info tags

#include "VcfFormat.h"
#include "MiscUtil.h"
#include "ExtendParameters.h"
#include "TVCVersion.h"


// current date string in YYYYMMDD format
string dateStr()
{
  time_t rawtime;
  struct tm* timeinfo;
  char buffer[80];
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer, 80, "%Y%m%d", timeinfo);
  return string(buffer);
}


string get_time_iso_string(time_t time)
{
  char time_buffer[1024];
  strftime(time_buffer, 1024, "%Y-%m-%dT%H:%M:%S", localtime(&time));
  return string(time_buffer);
}


string getVCFHeader(const ExtendParameters *parameters, ReferenceReader& ref_reader, const vector<string>& sample_list, int primary_sample)
{
  stringstream headerss;
  headerss
    << "##fileformat=VCFv4.1" << endl
    << "##fileDate=" << dateStr() << endl
    << "##fileUTCtime=" << get_time_iso_string(time(NULL)) << endl
    << "##source=\"bamwalker " << TVCVersion::GetVersion() << "-" << TVCVersion::GetRelease() << " (" << TVCVersion::GetGitHash() << ") - Torrent Variant Caller - candidate generator\"" << endl
    << "##bams=\"";

  for(size_t i = 0; i < parameters->bams.size(); ++i) {
      if(i != 0) headerss << ",";
      headerss << parameters->bams[i];
  }
  headerss << "\"\n";

  if (not parameters->params_meta_name.empty())
    headerss << "##parametersName=\"" << parameters->params_meta_name << "\"" << endl;
  if (not parameters->params_meta_details.empty())
    headerss << "##parametersDetails=\"" << parameters->params_meta_details << "\"" << endl;

  if (not parameters->basecaller_version.empty())
    headerss << "##basecallerVersion=\"" << parameters->basecaller_version << "\"" << endl;
  if (not parameters->tmap_version.empty())
    headerss << "##tmapVersion=\"" << parameters->tmap_version << "\"" << endl;

  headerss << "##reference=" << parameters->fasta << endl;
  if (parameters->fasta == "GRCh38.p2") {headerss << "##masked_reference=GRCh38.p2.mask1" << endl;}

  string ref_filename = ref_reader.get_filename();
  size_t pos = ref_filename.rfind(".");
  if (pos != string::npos) {ref_filename = ref_filename.substr(0, pos);}
  pos = ref_filename.rfind("/");
  string mask_file = "maskfile_donot_remove.bed";
  if (pos != string::npos) {
    mask_file =  ref_filename.substr(0, pos+1)+mask_file;
    ref_filename = ref_filename.substr(pos + 1);
  }
  FILE *fp = fopen(mask_file.c_str(), "r");
  if (fp) {
    //read the header
    char line[1000];
    if (fgets(line, sizeof line, fp) and line[0] == '#') {
      char first[1000];
      sscanf(line, "%s", first);
      string tmp(first+1);
      headerss << "##maskVersion=" << tmp << endl;
    }
  }
  string chr_name;
  long chr_size;
  for (int index = 0; (index < ref_reader.chr_count()); ++index) {
    chr_name = ref_reader.chr_str(index);
    chr_size = ref_reader.chr_size(index);
    headerss << "##contig=<ID=" << chr_name << ",length=" << chr_size << ",assembly=" << ref_filename << ">" << endl;
  }

  headerss << "##phasing=none" << endl
    << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << endl

    << "##INFO=<ID=HS,Number=0,Type=Flag,Description=\"Indicate it is at a hot spot\">" << endl

    << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">" << endl
    << "##INFO=<ID=DPF,Number=1,Type=Integer,Description=\"Total forward read depth at the locus\">" << endl
    << "##INFO=<ID=DPR,Number=1,Type=Integer,Description=\"Total reverse read depth at the locus\">" << endl
    << "##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observations\">" << endl
    << "##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observations\">" << endl

    << "##INFO=<ID=SRF,Number=1,Type=Integer,Description=\"Number of reference observations on the forward strand\">" << endl
    << "##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Number of reference observations on the reverse strand\">" << endl
    << "##INFO=<ID=SAF,Number=A,Type=Integer,Description=\"Alternate allele observations on the forward strand\">" << endl
    << "##INFO=<ID=SAR,Number=A,Type=Integer,Description=\"Alternate allele observations on the reverse strand\">" << endl

    << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">" << endl

    << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">" << endl;

  if(parameters->my_controls.use_lod_filter){
    headerss << "##INFO=<ID=LOD,Number=A,Type=Float,Description=\"Limit of Detection at genomic location.\">" << endl;
  }

  headerss << "##FILTER=<ID=NOCALL,Description=\"Generic filter. Filtering details stored in FR info tag.\">" << endl;

  headerss
    << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl
    << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">" << endl

    << "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">" << endl

    << "##FORMAT=<ID=SRF,Number=1,Type=Integer,Description=\"Number of reference observations on the forward strand\">" << endl
    << "##FORMAT=<ID=SRR,Number=1,Type=Integer,Description=\"Number of reference observations on the reverse strand\">" << endl
    << "##FORMAT=<ID=SAF,Number=A,Type=Integer,Description=\"Alternate allele observations on the forward strand\">" << endl
    << "##FORMAT=<ID=SAR,Number=A,Type=Integer,Description=\"Alternate allele observations on the reverse strand\">" << endl

    << "##FORMAT=<ID=RQ,Number=1,Type=Float,Description=\"Mean base quality in the reference allele\">" << endl
    << "##FORMAT=<ID=RQF,Number=1,Type=Float,Description=\"Mean base quality in the reference allele on the forward strand\">" << endl
    << "##FORMAT=<ID=RQR,Number=1,Type=Float,Description=\"Mean base quality in the reference allele on the reverse strand\">" << endl

    << "##FORMAT=<ID=AQ,Number=1,Type=Float,Description=\"Mean base quality in the alternate allele\">" << endl
    << "##FORMAT=<ID=AQF,Number=1,Type=Float,Description=\"Mean base quality in the alternate allele on the forward strand\">" << endl
    << "##FORMAT=<ID=AQR,Number=1,Type=Float,Description=\"Mean base quality in the alternate allele on the reverse strand\">" << endl;

  headerss << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  // Ensure primary sample is always in the first column (IR req)
  headerss << "\t" << sample_list[primary_sample];
  for (size_t i = 0; i < sample_list.size(); i++)
    if (i != (size_t)primary_sample)
      headerss << "\t" << sample_list[i];

  return headerss.str();
}

void ClearVal(vcf::Variant &var, const char *clear_me, const string &sample_name){
  map<string, vector<string> >::iterator it;
  it = var.samples[sample_name].find(clear_me);
  if (it != var.samples[sample_name].end())
    var.samples[sample_name][clear_me].clear();
};

void ClearVal(vcf::Variant &var, const char *clear_me){
  map<string, vector<string> >::iterator it;
  it = var.info.find(clear_me);
  if (it != var.info.end())
    var.info[clear_me].clear();
};


//clear all the info tags
void clearInfoTags(vcf::Variant &var) {
  map<string, vector<string> >::iterator it;

  it = var.info.find("RO");
  if (it != var.info.end())
    var.info["RO"].clear();

  it = var.info.find("AO");
  if (it != var.info.end())
    var.info["AO"].clear();

  it = var.info.find("SAF");
  if (it != var.info.end())
    var.info["SAF"].clear();

  it = var.info.find("SAR");
  if (it != var.info.end())
    var.info["SAR"].clear();

  it = var.info.find("SRF");
  if (it != var.info.end())
    var.info["SRF"].clear();

  it = var.info.find("SRR");
  if (it != var.info.end())
    var.info["SRR"].clear();

  it = var.info.find("DP");
  if (it != var.info.end())
    var.info["DP"].clear();


  it = var.info.find("MLLD");
  if (it != var.info.end())
    var.info["MLLD"].clear();

  it = var.info.find("LOD");
  if (it != var.info.end())
    var.info["LOD"].clear();
}

void NullInfoFields(vcf::Variant &var, bool use_position_bias){
  clearInfoTags(var);
  for (vector<string>::iterator I = var.alt.begin(); I != var.alt.end(); ++I) {
    var.info["AO"].push_back(convertToString(0));
    var.info["SAF"].push_back(convertToString(0));
    var.info["SAR"].push_back(convertToString(0));
  }
  var.info["DP"].push_back(convertToString(0));
  var.info["RO"].push_back(convertToString(0));
  var.info["SRF"].push_back(convertToString(0));
  var.info["SRR"].push_back(convertToString(0));
}

// set up format string
void SetUpFormatString(vcf::Variant &var) {
  var.format.clear();
  var.format.push_back("DP");
  var.format.push_back("RO");
  var.format.push_back("AO");
  var.format.push_back("SAR");
  var.format.push_back("SAF");
  var.format.push_back("SRF");
  var.format.push_back("SRR");
  var.format.push_back("RQ");
  var.format.push_back("RQF");
  var.format.push_back("RQR");
  var.format.push_back("AQ");
  var.format.push_back("AQF");
  var.format.push_back("AQR");
}


int CalculateWeightOfVariant(vcf::Variant &current_variant){

  map<string, vector<string> >::iterator it;
  int weight;

  it = current_variant.info.find("DP");
  if (it != current_variant.info.end())
    weight = atoi(current_variant.info["DP"][0].c_str()); // or is this current sample ident?
  else weight = 1;
  return(weight);
}

float RetrieveQualityTagValue(vcf::Variant &current_variant, const string &tag_wanted, int _allele_index, const string& sample_name){

  map<string, vector<string> >::iterator it;
  float weight;

  it = current_variant.samples[sample_name].find(tag_wanted);
  if (it != current_variant.samples[sample_name].end()){
    // if the index is valid...
    if (current_variant.samples[sample_name][tag_wanted].size()> (unsigned int) _allele_index)
      weight = atof(current_variant.samples[sample_name][tag_wanted][_allele_index].c_str()); // or is this current sample ident?
    else
      weight = 0.0f;
  }
  else weight = 0.0f;
  return(weight);
}

float RetrieveQualityTagValue(vcf::Variant &current_variant, const string &tag_wanted, int _allele_index){

  map<string, vector<string> >::iterator it;
  float weight;

  it = current_variant.info.find(tag_wanted);
  if (it != current_variant.info.end()){
    // if the index is valid...
    if (current_variant.info[tag_wanted].size()> (unsigned int) _allele_index)
      weight = atof(current_variant.info[tag_wanted][_allele_index].c_str()); // or is this current sample ident?
    else
      weight = 0.0f;
  }
  else weight = 0.0f;
  return(weight);
}
// XXX

void NullFilterReason(vcf::Variant &candidate_variant, const string &sample_name){
}

void AddFilterReason(vcf::Variant &candidate_variant, string &additional_reason, const string &sample_name){
  int alt_allele_index = 0;
  for (vector<string>::iterator I = candidate_variant.alt.begin(); I != candidate_variant.alt.end(); ++I) {
    AddFilterReason(candidate_variant, additional_reason, sample_name, alt_allele_index);
    alt_allele_index++;
  }
}

void AddFilterReason(vcf::Variant &candidate_variant, string &additional_reason, const string &sample_name, unsigned int alt_allele_index){
}

void AddInfoReason(vcf::Variant &candidate_variant, string &additional_reason, const string &sample_name){
}

