/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#ifndef STACKENGINE_H
#define STACKENGINE_H

#include "api/BamReader.h"

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <vector>

#include "ExtendedReadInfo.h"
#include "CrossHypotheses.h"
#include "ExtendParameters.h"
#include "PosteriorInference.h"
#include "ClassifyVariant.h"

using namespace std;

// what is the set of parameters describing a complete set of latent variables under some state
class LatentSlate{
  public:
    // frequency
    PosteriorInference cur_posterior;
    // information by strand
    //    PosteriorInference fwd_posterior;
    //    PosteriorInference rev_posterior;

    // // track important latent variables
    // BasicBiasGenerator bias_generator;
    // // track filter
    // BiasChecker bias_checker;
    //
    // // and the other important set of latent variables
    // StrandedSigmaGenerator sigma_generator;
    //
    // // and the third set
    // BasicSkewGenerator skew_generator;

    bool detailed_integral;
    int max_iterations;
    int iter_done;
    vector<float> ll_at_stage;
    vector<float> start_freq_of_winner;

    void FastExecuteInference(ShortStack &total_theory, bool update_frequency, vector<float> &start_frequency);
    void LocalExecuteInference(ShortStack &total_theory, bool update_frequency, vector<float> &start_frequency);
    void FastStep(ShortStack &total_theory, bool update_frequency);
    void DetailedStep(ShortStack &total_theory, bool update_frequency);
    void ScanStrandPosterior(ShortStack &total_theory, bool vs_ref, int max_detail_level);
    void PropagateTuningParameters(EvalTuningParameters &my_params, int num_hyp_no_null);
    LatentSlate(){
      max_iterations = 10;
      detailed_integral = true;
      iter_done = 0;
    };
};

class HypothesisStack{
  public:
    // latent variables under states of the world
    LatentSlate cur_state;
    ShortStack total_theory;
    EvalTuningParameters my_params;
    bool try_alternatives;

    vector<float> ll_record;
    vector<vector <float> > try_hyp_freq;

    HypothesisStack(){
      DefaultValues();
    }
    void AllocateFrequencyStarts(int num_hyp_no_null, vector<AlleleIdentity> &allele_identity_vector);
    void DefaultValues();
    void PropagateTuningParameters(int num_hyp_no_null);

    // starting to make inferences
    void RestoreFullInference();
    void SetAlternateFromMain();
    void ExecuteExtremeInferences();
    void TriangulateRestart();
    float ExecuteOneRestart(vector<float> &restart_hyp, int max_detail_level = 0);
    void ExecuteInference();
    void InitForInference(PersistingThreadObjects &thread_objects, vector<const Alignment *>& read_stack, const InputStructures &global_context, int num_hyp_no_null, vector<AlleleIdentity> &allele_identity_vector);

    // tool for posterior density estimation
    bool CallGermline(float hom_safety, int &genotype_call, float &quasi_phred_quality_score, float &reject_status_quality_score);

    float ReturnMaxLL();
};

class PositionInBam;

class Eval {
  public:
    // Raw read information
    vector<const Alignment *> read_stack;    //!< Reads spanning the variant position
    // Raw alleles information
    vcf::Variant *         variant;                 //!< VCF record of this variant position
    vector<AlleleIdentity> allele_identity_vector;  //!< Detailed information for each candidate allele
    LocalReferenceContext  seq_context;             //!< Reference context of this variant position
    int                    multiallele_window_start;
    int                    multiallele_window_end;
    vector<string>         info_fields;
    bool                   doRealignment;
    // Allele evaluation information
    HypothesisStack allele_eval;
    vector<int> diploid_choice;

    Eval(vcf::Variant &candidate_variant) {
      diploid_choice.assign(2,0);
      diploid_choice[1] = 1; // ref = 0, alt = 1
      variant = &candidate_variant;
      info_fields.clear();
      multiallele_window_start = -1;
      multiallele_window_end = -1;
      doRealignment = false;
      is_detect_best_multi_allele_pair_done_ = false;
      is_hard_classification_for_reads_done_ = false;
      is_hard_classification_for_families_done_ = false;
    };

    //! @brief  Create a detailed picture about this variant and all its alleles
    void SetupAllAlleles(const ExtendParameters &parameters, const InputStructures &global_context,
        const ReferenceReader &ref_reader, int chr_idx);
    void FilterAllAlleles(const ClassifyFilters &filter_variant, const vector<VariantSpecificParams>& variant_specific_params);
    void StackUpOneVariant(const ExtendParameters &parameters, const PositionInProgress& bam_position, int sample_index);
    void SpliceAllelesIntoReads(PersistingThreadObjects &thread_objects, const InputStructures &global_context,
        const ExtendParameters &parameters, const ReferenceReader &ref_reader, int chr_idx);
    void ApproximateHardClassifierForReads();
    void ApproximateHardClassifierForFamilies(); // calculate the family id for cfDNA
    void ScanSupportingEvidence(float &mean_ll_delta, int i_allele);
    int DetectBestMultiAllelePair();
    void ComputePosteriorGenotype(int _alt_allele_index,float local_min_allele_freq, int &genotype_call,
        float &gt_quality_score, float &reject_status_quality_score);
    void MultiAlleleGenotype(float local_min_allele_freq, vector<int> &genotype_component,
        float &gt_quality_score, float &reject_status_quality_score,
        int max_detail_level = 0);

    friend void GlueOutputVariant(Eval &my_ensemble, VariantCandidate &candidate_variant, const ExtendParameters &parameters, int _best_allele_index, int sample_index); // I want to access your private members

    // The following private members are used only in the internal steps at
    // a) int DetectBestMultiAllelePair()
    // b) void ApproximateHardClassifierForReads()
    // c) void ApproximateHardClassifierForFamilies()
  private:
    bool is_detect_best_multi_allele_pair_done_ = false;
    // Tvc used to compute read_id etc. twice. This is not computationally efficient.
    // Now we do it just once.
    // The following private members are the results of approximate hard classification for reads
    bool is_hard_classification_for_reads_done_ = false;
    vector<int> read_id_;        // vector of allele ids per read, -1 = outlier, 0 = ref, >0 real allele
    vector<bool> strand_id_;     // vector of forward (true) or reverse (false) per read
    // for each variant, calculate its' position within the soft clipped read distance to left and distance to right
    vector<int> dist_to_left_;   // vector of distances from allele position to left soft clip per read
    vector<int> dist_to_right_;  // vector of distances from allele position to left soft clip per read

    // The following private members are the results of approximate hard classification for families
    bool is_hard_classification_for_families_done_ = false;
    vector<int> family_id_; // similar to read_id, but it is for family
    vector<bool> family_strand_id_; // similar to strand_id, but it is for family
    vector<float> max_family_responsibility_; // maximum family responsibility each the family in family_id
};


#endif // STACKENGINE_H
