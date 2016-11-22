/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */


#ifndef DECISIONTREEDATA_H
#define DECISIONTREEDATA_H


#include "api/BamReader.h"

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <vector>

#include "ClassifyVariant.h"
#include "StackEngine.h"
#include "VcfFormat.h"

class EvaluatedGenotype{
public:
  bool genotype_already_set;
  float evaluated_genotype_quality;
  float evaluated_variant_quality;
  vector<int> genotype_component;

  EvaluatedGenotype(){
    genotype_already_set = false;
    evaluated_genotype_quality = 0.0f;
    evaluated_variant_quality = 0.0f;
    genotype_component.assign(2,0); // 0/0 = reference call

  };
  string GenotypeAsString();
  bool IsReference();
};

// all the data needed to make a decision for filtration
// characterize the variant, and the outcome from whatever evaluator we use
class DecisionTreeData {
  public:

    vcf::Variant * variant;                         //!< VCF record of this variant position
    vector<AlleleIdentity> allele_identity_vector;  //!< Detailed information for each candidate allele
    vector<string>         info_fields;             //!< Additional information to be printed out in vcf FR tag

    MultiBook all_summary_stats;

    vector<VariantOutputInfo> summary_info_vector;

    vector<int> filteredAllelesIndex;

	map<string, float> variant_quality;


    bool best_variant_filtered;

    bool best_allele_set;
    int best_allele_index;
    bool isBestAlleleSNP;
    bool reference_genotype;

    EvaluatedGenotype eval_genotype;

    float tune_xbias; // not tuned, removed from filters
    float tune_sbias;

    DecisionTreeData(vcf::Variant &candidate_variant) /*: multi_allele(candidate_variant)*/ {
      variant = &candidate_variant;
      best_allele_set = false;
      best_allele_index = 0;
      best_variant_filtered=false;
      isBestAlleleSNP = false;
      reference_genotype = false;


      tune_xbias = 0.005f; // tune calculation of chi-square bias = proportioinal variance by frequency
      tune_sbias = 0.5f; // safety factor for small allele counts for transformed strand bias
    };

    void OverrideFilter(vcf::Variant &candidate_variant, string & _filter_reason, int _allele, const string &sample_name);
    void FilterOneAllele(int i_alt,VariantOutputInfo &l_summary_info,
                         AlleleIdentity &l_variant_identity, const ControlCallAndFilters &my_filters,
                         const VariantSpecificParams& variant_specific_params);
    void FilterAlleles(const ControlCallAndFilters &my_filters, const vector<VariantSpecificParams>& variant_specific_params);

    void AccumulateFilteredAlleles();

    void BestSNPsSuppressInDels(bool heal_snps);
    void FindBestAlleleIdentity();
    void FindBestAlleleByScore();

    void GenotypeFromBestAlleleIndex(vcf::Variant &candidate_variant, const ExtendParameters &parameters);
    void GenotypeFromEvaluator(vcf::Variant &candidate_variant, const ExtendParameters &parameters, const string &sample_name);

    void FilterMyCandidate(vcf::Variant &candidate_variant, const ExtendParameters &parameters, const string &sample_name);
    void BestAlleleFilterMyCandidate(vcf::Variant &candidate_variant, const ExtendParameters &parameters);
    void GenotypeAlleleFilterMyCandidate(vcf::Variant &candidate_variant, const ExtendParameters &parameters, const string &sample_name);

    void SimplifySNPsIfNeeded(VariantCandidate &candidate_variant, const ExtendParameters &parameters, const string &sample_name);


    bool SetGenotype(vcf::Variant &candidate_variant, const ExtendParameters &parameters, float gt_quality);
    void DecisionTreeOutputToVariant(VariantCandidate &candidate_variant, const ExtendParameters &parameters, int sample_index);

    void AggregateFilterInformation(vcf::Variant &candidate_variant, const vector<VariantSpecificParams>& variant_specific_params, const ExtendParameters &parameters, const string &sample_name);
    void FillInFiltersAtEnd(VariantCandidate &candidate_variant,const ExtendParameters &parameters, const string &sample_name);



    void SetupFromMultiAllele(const Evaluator &eval);
    void AddStrandBiasTags(vcf::Variant &candidate_variant, const string &sample_name);
    void AddPositionBiasTags(vcf::Variant &candidate_variant, const string &sample_name);
    void AddLodTags(vcf::Variant &candidate_variant, const string &sample_name, const ExtendParameters &parameters);
    void AddCountInformationTags(vcf::Variant &candidate_variant, const string &sampleName, const ExtendParameters &parameters);
    string GenotypeStringFromAlleles(std::vector<int> &allowedGenotypes, bool refAlleleFound);
    bool AllowedGenotypesFromSummary(std::vector<int> &allowedGenotypes);
    string GenotypeFromStatus(vcf::Variant &candidate_variant, const ExtendParameters &parameters);
    void FilterAlleleHypothesisBias(vcf::Variant &candidate_variant, float ref_bias, float var_bias, float threshold_bias, int _allele, const string &sample_name);
    void FilterOnSpecialTags(vcf::Variant &candidate_variant, const ExtendParameters &parameters, const vector<VariantSpecificParams>& variant_specific_params, const string &sample_name);
    void FilterOnStringency(vcf::Variant &candidate_variant, const float data_quality_stringency,  int _check_allele_index, const string &sample_name);
    void FilterOnPositionBias(const string& allele, int i_alt, MultiBook &m_summary_stats, VariantOutputInfo &l_summary_info, const ControlCallAndFilters &my_filters, const VariantSpecificParams& variant_specific_params);
    void FilterOnLod(const string& allele, int i_alt, MultiBook &m_summary_stats, VariantOutputInfo &l_summary_info, const ControlCallAndFilters &my_filters, const BasicFilters &basic_filter, bool is_reference_call);
    void FilterBlackList(const vector<VariantSpecificParams>& variant_specific_params);
    void FilterSSE(vcf::Variant &candidate_variant, const ClassifyFilters &filter_variant, const vector<VariantSpecificParams>& variant_specific_params, const string &sample_name);
};
void FilterByBasicThresholds(stringstream &s, int i_alt, MultiBook &m_summary_stats,
                             VariantOutputInfo &l_summary_info,
                             const BasicFilters &basic_filter, float tune_xbias, float tune_bias, bool is_reference_call);

void AutoFailTheCandidate(vcf::Variant &candidate_variant, bool use_position_bias, const string &sample_name);
float FreqThresholdByType(AlleleIdentity &variant_identity, const ControlCallAndFilters &my_controls, const VariantSpecificParams& variant_specific_params);
void DetectSSEForNoCall(AlleleIdentity &var_identity, float sseProbThreshold, float minRatioReadsOnNonErrorStrand, float relative_safety_level, vcf::Variant &candidate_variant, unsigned _altAlleIndex);
void SetQualityByDepth(vcf::Variant &candidate_variant, map<string, float>& variant_quality, const string &sample_name);
void RemoveVcfInfo(vcf::Variant &variant, const vector<string> &info_fields, const string &sample_name, int sample_index);

#endif // DECISIONTREEDATA_H
