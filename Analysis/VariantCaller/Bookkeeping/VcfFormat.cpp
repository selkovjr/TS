/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     VcfFormat.cpp
//! @ingroup  VariantCaller
//! @brief    Vcf file formatting & info tags

#include "VcfFormat.h"
#include "MiscUtil.h"
#include "ExtendParameters.h"
#include "IonVersion.h"


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
    << "##source=\"bamwalker " << IonVersion::GetVersion() << "-" << IonVersion::GetRelease() << " (" << IonVersion::GetGitHash() << ") - Torrent Variant Caller - candidate generator\"" << endl
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

    << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">" << endl
    << "##INFO=<ID=HRUN,Number=A,Type=Integer,Description=\"Run length: the number of consecutive repeats of the alternate allele in the reference genome\">" << endl;

  if(parameters->my_controls.use_lod_filter){
    headerss << "##INFO=<ID=LOD,Number=A,Type=Float,Description=\"Limit of Detection at genomic location.\">" << endl;
  }

  headerss << "##FILTER=<ID=NOCALL,Description=\"Generic filter. Filtering details stored in FR info tag.\">" << endl;

  headerss << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
    << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl
    << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">" << endl

    << "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">" << endl

    << "##FORMAT=<ID=SRF,Number=1,Type=Integer,Description=\"Number of reference observations on the forward strand\">" << endl
    << "##FORMAT=<ID=SRR,Number=1,Type=Integer,Description=\"Number of reference observations on the reverse strand\">" << endl
    << "##FORMAT=<ID=SAF,Number=A,Type=Integer,Description=\"Alternate allele observations on the forward strand\">" << endl
    << "##FORMAT=<ID=SAR,Number=A,Type=Integer,Description=\"Alternate allele observations on the reverse strand\">" << endl;

  // If we want to output multiple min-allele-freq
  if(parameters->program_flow.is_multi_min_allele_freq) {
    string multi_min_allele_freq_size = "?";
    string snp_multi_min_allele_freq_size = convertToString(parameters->program_flow.snp_multi_min_allele_freq.size());
    string mnp_multi_min_allele_freq_size = convertToString(parameters->program_flow.mnp_multi_min_allele_freq.size());
    string indel_multi_min_allele_freq_size = convertToString(parameters->program_flow.indel_multi_min_allele_freq.size());
    // the union of all var types
    headerss<< "##FORMAT=<ID=MUAF,Number=" + multi_min_allele_freq_size + ",Type=Float,Description=\"Union of multi-min-allele-freq associated with TYPE.\">" << endl
      << "##FORMAT=<ID=MUQUAL,Number=" + multi_min_allele_freq_size + ",Type=Float,Description=\"QUAL scores for MAF.\">" << endl
      << "##FORMAT=<ID=MUGT,Number=" + multi_min_allele_freq_size + ",Type=String,Description=\"Genotypes for MAF.\">" << endl
      << "##FORMAT=<ID=MUGQ,Number=" + multi_min_allele_freq_size + ",Type=Integer,Description=\"Genotype quality scores for MAF.\">" << endl
      // for snp
      << "##FORMAT=<ID=SMAF,Number=" + snp_multi_min_allele_freq_size + ",Type=Float,Description=\"Values of snp-multi-min-allele-freq.\">" << endl
      << "##FORMAT=<ID=SMQUAL,Number=" + snp_multi_min_allele_freq_size + ",Type=Float,Description=\"QUAL scores for snp-multi-min-allele-freq.\">" << endl
      << "##FORMAT=<ID=SMGT,Number=" + snp_multi_min_allele_freq_size + ",Type=String,Description=\"Genotypes for snp-multi-min-allele-freq.\">" << endl
      << "##FORMAT=<ID=SMGQ,Number=" + snp_multi_min_allele_freq_size + ",Type=Integer,Description=\"Genotype quality scores for snp-multi-min-allele-freq.\">" << endl
      // for mnp
      << "##FORMAT=<ID=MMAF,Number=" + mnp_multi_min_allele_freq_size + ",Type=Float,Description=\"Values of mnp-multi-min-allele-freq.\">" << endl
      << "##FORMAT=<ID=MMQUAL,Number=" + mnp_multi_min_allele_freq_size + ",Type=Float,Description=\"QUAL scores for mnp-multi-min-allele-freq.\">" << endl
      << "##FORMAT=<ID=MMGT,Number=" + mnp_multi_min_allele_freq_size + ",Type=String,Description=\"Genotypes for mnp-multi-min-allele-freq.\">" << endl
      << "##FORMAT=<ID=MMGQ,Number=" + mnp_multi_min_allele_freq_size + ",Type=Integer,Description=\"Genotype quality scores for mnp-multi-min-allele-freq.\">" << endl
      // for indel
      << "##FORMAT=<ID=IMAF,Number=" + indel_multi_min_allele_freq_size + ",Type=Float,Description=\"Values of indel-multi-min-allele-freq.\">" << endl
      << "##FORMAT=<ID=IMQUAL,Number=" + indel_multi_min_allele_freq_size + ",Type=Float,Description=\"QUAL scores for indel-multi-min-allele-freq.\">" << endl
      << "##FORMAT=<ID=IMGT,Number=" + indel_multi_min_allele_freq_size + ",Type=String,Description=\"Genotypes for indel-multi-min-allele-freq.\">" << endl
      << "##FORMAT=<ID=IMGQ,Number=" + indel_multi_min_allele_freq_size + ",Type=Integer,Description=\"Genotype quality scores for indel-multi-min-allele-freq.\">" << endl;
  }


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


  it = var.info.find("HRUN");
  if (it != var.info.end())
    var.info["HRUN"].clear();


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
    var.info["HRUN"].push_back(convertToString(0));
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

// if, for example, missing data
void NullGenotypeAllSamples(vcf::Variant & candidate_variant)
{
  vector<string>& sampleNames = candidate_variant.sampleNames;

  for (vector<string>::iterator its = sampleNames.begin(); its != sampleNames.end(); ++its) {
    string& sampleName = *its;
    map<string, vector<string> >& sampleOutput = candidate_variant.samples[sampleName];

    // for multi-min-allele-freq
    vector<string> tags_for_multi_min_allele_freq = {"MUAF", "MUQUAL", "MUGQ", "MUGT",
      "SMAF", "SMQUAL", "SMGQ", "SMGT",
      "MMAF", "MMQUAL", "MMGQ", "MMGT",
      "IMAF", "IMQUAL", "IMGQ", "IMGT",
      "HMAF", "HMQUAL", "HMGQ", "HMGT"};
    map<string, vector<string> >::iterator it;
    for(unsigned int i_tag = 0; i_tag < tags_for_multi_min_allele_freq.size(); ++i_tag){
      it = sampleOutput.find(tags_for_multi_min_allele_freq[i_tag]);
      if (it != sampleOutput.end()){
        it->second.clear();
        it->second.push_back(".");
      }
    }
  }
}

void OverwriteGenotypeForOneSample(vcf::Variant &candidate_variant, const string &my_sample_name, string &my_genotype, float genotype_quality){
  // will create entry if one does not exist

  map<string, vector<string> >& sampleOutput = candidate_variant.samples[my_sample_name];
  // clear existing values
  map<string, vector<string> >::iterator it;
  it = sampleOutput.find("GT");
  if (it != sampleOutput.end())    sampleOutput["GT"].clear();
  it = sampleOutput.find("GQ");
  if (it != sampleOutput.end())     sampleOutput["GQ"].clear();


  sampleOutput["GT"].push_back(my_genotype);
  // genotype quality should be an "int"
  //cout << "Storing Genotype = " << my_genotype << endl;
  sampleOutput["GQ"].push_back(convertToString((int)genotype_quality));

}

void DetectAndSetFilteredGenotype(vcf::Variant &candidate_variant, map<string, float>& variant_quality, const string &sample_name){
  if (candidate_variant.isFiltered){
    string no_call_genotype = "./.";
    float original_quality = variant_quality[sample_name];
    OverwriteGenotypeForOneSample(candidate_variant, sample_name, no_call_genotype, original_quality);
  }
}


void StoreGenotypeForOneSample(vcf::Variant &candidate_variant, const string &sample_name, string &my_genotype, float genotype_quality, bool multisample) {
  vector<string> sampleNames = candidate_variant.sampleNames;

  if (multisample) {
    map<string, vector<string> >& sampleOutput = candidate_variant.samples[sample_name];
    sampleOutput["GT"].clear();
    sampleOutput["GT"].push_back(my_genotype);
    //cout << "Storing Genotype = " << my_genotype << endl;
    sampleOutput["GQ"].clear();
    sampleOutput["GQ"].push_back(convertToString((int)genotype_quality));
  }
  else {
    for (vector<string>::iterator its = sampleNames.begin(); its != sampleNames.end(); ++its) {
      string& sampleName = *its;
      map<string, vector<string> >& sampleOutput = candidate_variant.samples[sampleName];
      map<string, vector<string> >::iterator it;
      it = sampleOutput.find("GT");
      if (it != sampleOutput.end())    sampleOutput["GT"].clear();
      it = sampleOutput.find("GQ");
      if (it != sampleOutput.end())     sampleOutput["GQ"].clear();

      if (sampleName.compare(sample_name) == 0) { //sample of interest
        //cout << "isNocall " << isNoCall << " genotype = " << my_genotype << endl;

        // if no-call, will reset this entry as a final step, but until then, give me my genotype
        sampleOutput["GT"].push_back(my_genotype);
        //cout << "Storing Genotype = " << my_genotype << endl;
        sampleOutput["GQ"].push_back(convertToString((int)genotype_quality));

      } else { //for all other samples in BAM file just make a no-call at this point.
        sampleOutput["GT"].push_back("./.");
        sampleOutput["GQ"].push_back(convertToString(0));
      }
    }
  }
}

void SetFilteredStatus(vcf::Variant &candidate_variant, bool isFiltered) {
  // filtering only sets the column to no-call
  // choice to put  in filtered file is handled by writing it out
  // genotype is modified by genotype
  if (isFiltered) {
    candidate_variant.filter = "NOCALL";
    candidate_variant.isFiltered = true;
  } else {
    candidate_variant.filter = "PASS";
    candidate_variant.isFiltered = false;
  }

}



void AdjustFDPForRemovedAlleles(vcf::Variant &candidate_variant, int filtered_allele_index, string sampleName)
{
  // first do the "info" tag as it is easier to find
  map<string, vector<string> >::iterator it;
  int total_depth=0;

  it = candidate_variant.info.find("FDP");
  if (it != candidate_variant.info.end())
    total_depth = atoi(candidate_variant.info["FDP"][0].c_str()); // or is this current sample ident?

  int allele_depth = 0;
  it = candidate_variant.info.find("FAO");
  if (it != candidate_variant.info.end())
    allele_depth = atoi(candidate_variant.info["FAO"][filtered_allele_index].c_str());

  total_depth -= allele_depth;
  if (total_depth<0)
    total_depth = 0; // how can this happen?

  ClearVal(candidate_variant, "FDP");
  candidate_variant.info["FDP"].push_back(convertToString(total_depth));

  if (!sampleName.empty()) {
    map<string, vector<string> >& sampleOutput = candidate_variant.samples[sampleName];
    sampleOutput["FDP"].clear();
    sampleOutput["FDP"].push_back(convertToString(total_depth));
  }
}


void RemoveFilteredAlleles(vcf::Variant &candidate_variant, vector<int> &filtered_alleles_index, const string &sample_name) {
  //now that all possible alt. alleles are evaluated decide on which allele is most likely and remove any that
  //that does'nt pass score threshold. Determine Genotype based on alleles that have evidence.
  candidate_variant.updateAlleleIndexes();
  string my_healing_glow = "HEALED";
  vector<string> originalAltAlleles = candidate_variant.alt;
  if (originalAltAlleles.size() > 1  &&
      originalAltAlleles.size() > filtered_alleles_index.size()) {  //remove only when number of alleles more than number of filtered alleles
    //remove filtered alleles with no support
    string altStr;
    int index;
    for (size_t i = 0; i <filtered_alleles_index.size(); i++) {

      index = filtered_alleles_index[i];
      //generate allele index before removing alleles
      altStr = originalAltAlleles[index];
      // specify what the alleles removed are
      //my_healing_glow = "HEALED" + altStr;

      //altStr = (*candidate_variant)->alt[index];
      // Note: need to update index for adjustments
      //AdjustFDPForRemovedAlleles(candidate_variant, index, sample_name);
      //cout << "Removed Fitered allele: index = " << index << " allele = " << altStr << endl;
      // @TODO: removeAlt wrecks the genotype as well
      // fix so we don't remove genotype components.

      candidate_variant.removeAlt(altStr);
      candidate_variant.updateAlleleIndexes();
      // if we are deleting alleles, indicate data potentially damaged at this location
      AddInfoReason(candidate_variant, my_healing_glow, sample_name);
    }
  }
}

// this only needs to know candidate variant, nothing else
void AdjustAlleles(vcf::Variant &candidate_variant, int position_upper_bound)
{
  vector<string>& types  = candidate_variant.info["TYPE"];
  string& refAllele      = candidate_variant.ref;
  vector<string>& alts   = candidate_variant.alt;
  int position     = candidate_variant.position;
  int max_trim = refAllele.length();
  if (position_upper_bound)
    max_trim = min(max_trim,position_upper_bound - position);
  string& altAllele = alts[0];
  //nothing to do if there are multiple allels

  if (types.size() != 1)
    return;

  if ((types[0]).compare("snp") == 0 && refAllele.length() > 1 && refAllele.length() == altAllele.length())  {
    //need to adjust position only in cases where SNP is represent as MNV due to haplotyping - REF= TTC ALT = TTT
    for (int i = 0; i < max_trim; ++i) {
      if (refAllele[i] != altAllele[i]) {
        candidate_variant.position = position;
        candidate_variant.ref = refAllele.substr(i, 1);
        candidate_variant.alt[0] = altAllele.substr(i, 1);
        break;
      }
      position++;
    }

  }

}

