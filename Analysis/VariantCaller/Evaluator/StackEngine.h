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
#include "ExtendParameters.h"
#include "ClassifyVariant.h"
#include "LocusData.h"

using namespace std;

class PositionInBam;

class Evaluator {
  public:
    // Raw read information
    vector<const Alignment *> read_stack;    // Reads spanning the variant position
    vector<const Alignment *> read_stack_n;  // Normal sample subset
    vector<const Alignment *> read_stack_t;  // Tumor sample subset
    // Raw alleles information
    vcf::Variant *         variant;                 //!< VCF record of this variant position
    vector<AlleleIdentity> allele_identity_vector;  //!< Detailed information for each candidate allele
    LocalReferenceContext  seq_context;             //!< Reference context of this variant position
    int                    multiallele_window_start;
    int                    multiallele_window_end;
    vector<string>         info_fields;
    bool                   doRealignment;
    // Allele evaluation information
    LocusData allele_eval;
    vector<int> diploid_choice;

    Evaluator(vcf::Variant &candidate_variant) {
      diploid_choice.assign(2,0);
      diploid_choice[1] = 1; // ref = 0, alt = 1
      variant = &candidate_variant;
      info_fields.clear();
      multiallele_window_start = -1;
      multiallele_window_end = -1;
      doRealignment = false;
    };

    //! @brief  Create a detailed picture about this variant and all its alleles
    void SetupAllAlleles(const ExtendParameters &parameters, const InputStructures &global_context,
        const ReferenceReader &ref_reader, int chr_idx);
    void FilterAllAlleles(const ClassifyFilters &filter_variant, const vector<VariantSpecificParams>& variant_specific_params);
    void StackUpOneVariant(const ExtendParameters &parameters, VariantCandidate &candidate_variant, const PositionInProgress& bam_position);
    void SampleLikelihood(PersistingThreadObjects &thread_objects, const InputStructures &global_context,
        const ExtendParameters &parameters, const ReferenceReader &ref_reader, int chr_idx);
    void ScanSupportingEvidence(float &mean_ll_delta, int i_allele);
    void ComputePosteriorGenotype(int _alt_allele_index,float local_min_allele_freq, int &genotype_call,
        float &gt_quality_score, float &reject_status_quality_score);
    void MultiAlleleGenotype(float local_min_allele_freq, vector<int> &genotype_component,
        float &gt_quality_score, float &reject_status_quality_score,
        int max_detail_level = 0);

    friend void GlueOutputVariant(Evaluator &eval, VariantCandidate &candidate_variant, const ExtendParameters &parameters, int _best_allele_index, int sample_index); // I want to access your private members

  private:
    // Tvc used to compute read_id etc. twice. This is not computationally efficient.
    // Now we do it just once.
    // The following private members are the results of approximate hard classification for reads
    vector<int> read_id_;        // vector of allele ids per read, -1 = outlier, 0 = ref, >0 real allele
    vector<bool> strand_id_;     // vector of forward (true) or reverse (false) per read
    // for each variant, calculate its position within the soft clipped read distance to left and distance to right
    vector<int> dist_to_left_;   // vector of distances from allele position to left soft clip per read
    vector<int> dist_to_right_;  // vector of distances from allele position to left soft clip per read
};


#endif // STACKENGINE_H
