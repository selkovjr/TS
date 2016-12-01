/* Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     HandleVariant.cpp
//! @ingroup  VariantCaller
//! @brief    HP Indel detection

#include "HandleVariant.h"
#include "SpliceVariantHypotheses.h"
#include "DecisionTreeData.h"


void Evaluator::SampleLikelihood (
  PersistingThreadObjects &thread_objects,
  const InputStructures &global_context,
  const ExtendParameters &parameters,
  const ReferenceReader &ref_reader,
  int chr_idx,
  VariantCandidate &candidate_variant
) {
  bool changed_alignment;
  unsigned int  num_valid_reads = 0;
  unsigned int  num_realigned = 0;
  int  num_hyp_no_null = allele_identity_vector.size() + 1; // num alleles +1 for ref
  // generate null+ref+nr.alt hypotheses per read in the case of do_multiallele_eval
  allele_eval.alignments.resize(read_stack.size());

  cerr << "evaluating " << allele_eval.alignments.size() << " hypotheses\n";
  for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
    // --- New splicing function ---
    allele_eval.alignments[i_read].success =
      SpliceVariantHypotheses(
        *read_stack[i_read],
        *this,
        seq_context,
        thread_objects,
        allele_eval.alignments[i_read].basecall,
        allele_eval.alignments[i_read].error_prob,
        changed_alignment,
        global_context,
        ref_reader, chr_idx
    );

    if (allele_eval.alignments[i_read].success){
      num_valid_reads++;
      if (changed_alignment)
        num_realigned++;
    }
  }
  cerr << "num_valid_reads: " << num_valid_reads << endl;
  cerr << "num_realigned: " << num_realigned << endl;

  // Check how many reads had their alignment modified
  std::ostringstream my_info;
  my_info.precision(4);
  if (doRealignment and num_valid_reads > 0) {
    cerr << "checking modified reads\n";
    float frac_realigned = (float)num_realigned / (float)num_valid_reads;
    // And re-do splicing without realignment if we exceed the threshold
    if (frac_realigned > parameters.my_controls.filter_variant.realignment_threshold) {
      my_info << "SKIPREALIGNx" << frac_realigned;
      doRealignment = false;
      for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
        allele_eval.alignments[i_read].success =
          SpliceVariantHypotheses(
            *read_stack[i_read],
            *this,
            seq_context,
            thread_objects,
            allele_eval.alignments[i_read].basecall,
            allele_eval.alignments[i_read].error_prob,
            changed_alignment,
            global_context,
            ref_reader, chr_idx
          );
      }
    }
    else {
      my_info << "REALIGNEDx" << frac_realigned;
    }
    info_fields.push_back(my_info.str());
  }
  cerr << my_info.str() << endl;

  // Tabulate frequencies
  allele_eval.TabulateFrequencies(candidate_variant, read_stack);

  // Computed the log likelihood of joint samples
  double snv_likelihood = 0;
  for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
    double p = 0;
    const string basecall = allele_eval.alignments[i_read].basecall;
    double e = allele_eval.alignments[i_read].error_prob;
    cerr << basecall << endl;
    for (const string &true_base: {"A", "C", "G", "T"}) {
      double p_b_given_a = basecall == true_base ? 1 - e : e / 3;
      cerr << "  P(" << basecall << "|" << true_base << ") = " << p_b_given_a << endl;
      cerr << "  f = " << allele_eval.freq(basecall) << endl;
      p += p_b_given_a * allele_eval.freq(basecall);
    }
    snv_likelihood += log10(p);
  }

  cerr << "total: " << allele_eval.alignments.size() << endl;
  for ( const auto &a: allele_eval.tumor_allele_count ) {
    cerr << "tumor " << a.first << " -> " << a.second << " (" << allele_eval.freq(a.first, "tumor") << ")\n";
  }
  for ( const auto &a: allele_eval.normal_allele_count ) {
    cerr << "normal " << a.first << " -> " << a.second << " (" << allele_eval.freq(a.first, "normal") << ")\n";
  }
  for (const string &base: {"A", "C", "G", "T"}) {
    cerr << base << " -> " << allele_eval.freq(base) << endl;
  }

  cerr << "sample likelihood: " << snv_likelihood << endl;
}




void SummarizeInfoFields(Evaluator &eval, vcf::Variant &candidate_variant, int _cur_allele_index, const string &sample_name) {

  float mean_ll_delta;

  eval.ScanSupportingEvidence(mean_ll_delta, _cur_allele_index);

  candidate_variant.info["MLLD"].push_back(convertToString(mean_ll_delta));

  int fwd_strand = 0;
  int rev_strand = 1;

  int var_hyp = 1;
}


// The structure InferenceResult is used for void MultiMinAlleleFreq(...) only.
// qual, gt, gq are the QUAL, GT, GQ computed for min_allele_freq, respectively
struct InferenceResult{
  float min_allele_freq;
  float qual;
  string gt;
  int gq;
  bool operator<(const InferenceResult &rhs) const { return min_allele_freq < rhs.min_allele_freq; } // use the operator "<" for "std::sort"
  bool operator==(const float &rhs) const { return min_allele_freq == rhs; } // use the operator "==" for "std::find"
};

void MultiMinAlleleFreq(Evaluator &eval, VariantCandidate &candidate_variant, int sample_index, ProgramControlSettings &program_flow, int max_detail_level){
  string sample_name = "";
  if(sample_index >= 0) {
    sample_name = candidate_variant.variant.sampleNames[sample_index];
  }

  // info tags needed for multi-min-allele-freq
  vector<string> tags_for_multi_min_allele_freq = {"MUAF", "MUQUAL", "MUGQ", "MUGT",
    "SMAF", "SMQUAL", "SMGQ", "SMGT",
    "MMAF", "MMQUAL", "MMGQ", "MMGT",
    "IMAF", "IMQUAL", "IMGQ", "IMGT",
    "HMAF", "HMQUAL", "HMGQ", "HMGT"};
  // Add the info tags for multi-min-allele-freq if we have not added yet.
  for(unsigned int i_tag = 0; i_tag < tags_for_multi_min_allele_freq.size(); ++i_tag){
    vector<string>::iterator it_format = find(candidate_variant.variant.format.begin(),
        candidate_variant.variant.format.end(),
        tags_for_multi_min_allele_freq[i_tag]);
    if(it_format == candidate_variant.variant.format.end()){
      candidate_variant.variant.format.push_back(tags_for_multi_min_allele_freq[i_tag]);
    }
  }

  // inference_results_union stores the inference results for the union of the min-allele-freq of all available variant types
  vector<InferenceResult> inference_results_union;
  bool is_snp_done = false;
  bool is_mnp_done = false;
  bool is_hotspot_done = false;
  bool is_indel_done = false;

  for(unsigned int alt_allele_index = 0; alt_allele_index < eval.allele_identity_vector.size(); ++alt_allele_index){
    // ptr_maf_vec = the pointer to the multi_min_allele_freq vector of which type of variant for this allele
    vector<float> *ptr_maf_vec = &(program_flow.snp_multi_min_allele_freq);
    string type_prefix = "S";

    // Note that no override here!
    if(eval.allele_identity_vector[alt_allele_index].status.isHotSpot){
      if(is_hotspot_done){
        continue;
      }
      ptr_maf_vec = &(program_flow.hotspot_multi_min_allele_freq);
      type_prefix = "H";
    }
    else if(eval.allele_identity_vector[alt_allele_index].ActAsSNP()){
      if(is_snp_done){
        continue;
      }
      ptr_maf_vec = &(program_flow.snp_multi_min_allele_freq);
      type_prefix = "S";
    }
    else if(eval.allele_identity_vector[alt_allele_index].ActAsMNP()){
      if(is_mnp_done){
        continue;
      }
      ptr_maf_vec = &(program_flow.mnp_multi_min_allele_freq);
      type_prefix = "M";
    }
    else if(eval.allele_identity_vector[alt_allele_index].ActAsHPIndel()){
      if(is_indel_done){
        continue;
      }
      ptr_maf_vec = &(program_flow.indel_multi_min_allele_freq);
      type_prefix = "I";
    }else{
      if(is_snp_done){
        continue;
      }
      ptr_maf_vec = &(program_flow.snp_multi_min_allele_freq);
      string type_prefix = "S";
    }

    for(unsigned int i_freq = 0; i_freq < ptr_maf_vec->size(); ++i_freq){
      float loc_min_allele_freq = ptr_maf_vec->at(i_freq);
      vector<InferenceResult>::iterator it = find(inference_results_union.begin(), inference_results_union.end(), loc_min_allele_freq);
      float loc_qual;
      int loc_gq;
      string loc_gt;

      if(it == inference_results_union.end()){ // This the first time we get loc_min_allele_freq
        int genotype_call;
        float evaluated_genotype_quality;
        // Let's do the inference for the given loc_min_allele_freq
        // eval.allele_eval.CallGermline(loc_min_allele_freq, genotype_call, evaluated_genotype_quality, loc_qual);

        vector<int> genotype_component = {eval.diploid_choice[0], eval.diploid_choice[1]}; // starts with het var

        if(true) { // genotype_call == 2){ //hom var
          genotype_component[0] = eval.diploid_choice[1];
        }
        else { // if(genotype_call == 0){ //hom ref
          genotype_component[1] = eval.diploid_choice[0];
        }

        loc_gt = convertToString(genotype_component[0]) + "/" + convertToString(genotype_component[1]);
        loc_gq = int(round(evaluated_genotype_quality)); // genotype quality is rounded as an integer.
        // append inference_results_union
        inference_results_union.push_back({loc_min_allele_freq, loc_qual, loc_gt, loc_gq});
      }
      else{ // We've seen loc_min_allele_freq before. Don't need to call CallGermline(...) again.
        loc_qual = it->qual;
        loc_gq = it->gq;
        loc_gt = it->gt;
      }

      // write the info tag for the corresponding var type
      candidate_variant.variant.samples[sample_name][type_prefix + "MUAF"].push_back(convertToString(loc_min_allele_freq));
      candidate_variant.variant.samples[sample_name][type_prefix + "MUQUAL"].push_back(convertToString(loc_qual));
      candidate_variant.variant.samples[sample_name][type_prefix + "MUGT"].push_back(loc_gt);
      candidate_variant.variant.samples[sample_name][type_prefix + "MUGQ"].push_back(convertToString(loc_gq));

      switch(type_prefix[0]){
        case 'S':
          is_snp_done = true;
          break;
        case 'M':
          is_mnp_done = true;
          break;
        case 'I':
          is_indel_done = true;
          break;
        case 'H':
          is_hotspot_done = true;
          break;
      }
    }
  }

  // sort inference_results_union according to min_allele_freq in the ascending order
  sort(inference_results_union.begin(), inference_results_union.end());
  // write the info tag for the union of min_allele_freq of the var types
  for(unsigned int i_freq = 0; i_freq < inference_results_union.size(); ++i_freq){
    candidate_variant.variant.samples[sample_name]["MUAF"].push_back(convertToString(inference_results_union[i_freq].min_allele_freq));
    candidate_variant.variant.samples[sample_name]["MUQUAL"].push_back(convertToString(inference_results_union[i_freq].qual));
    candidate_variant.variant.samples[sample_name]["MUGT"].push_back(inference_results_union[i_freq].gt);
    candidate_variant.variant.samples[sample_name]["MUGQ"].push_back(convertToString(inference_results_union[i_freq].gq));
  }
}


void GlueOutputVariant(Evaluator &eval, VariantCandidate &candidate_variant, const ExtendParameters &parameters, int _best_allele_index, int sample_index){
  string sample_name = "";
  if(sample_index >= 0) {
    sample_name = candidate_variant.variant.sampleNames[sample_index];
  }

  DecisionTreeData my_decision(*(eval.variant));
  my_decision.tune_sbias = parameters.my_controls.sbias_tune;
  my_decision.SetupFromMultiAllele(eval);

  float smallest_allele_freq = 1.0f;
  for (unsigned int _alt_allele_index = 0; _alt_allele_index < my_decision.allele_identity_vector.size(); _alt_allele_index++) {
    // for each alt allele, do my best
    // thresholds here can vary by >type< of allele
    float local_min_allele_freq = FreqThresholdByType(eval.allele_identity_vector[_alt_allele_index], parameters.my_controls,
        candidate_variant.variant_specific_params[_alt_allele_index]);

    if (local_min_allele_freq < smallest_allele_freq){
      smallest_allele_freq = local_min_allele_freq;  // choose least-restrictive amongst multialleles
    }

    /* The following piece of code seems redundant. Perhaps due to historical reasons?
       eval.ComputePosteriorGenotype(_alt_allele_index, local_min_allele_freq,
       my_decision.summary_info_vector[_alt_allele_index].genotype_call,
       my_decision.summary_info_vector[_alt_allele_index].gt_quality_score,
       my_decision.summary_info_vector[_alt_allele_index].variant_qual_score);
       */
    SummarizeInfoFields(eval, *(eval.variant), _alt_allele_index, sample_name);
  }

  my_decision.best_allele_index = _best_allele_index;
  my_decision.best_allele_set = true;

  // tell the evaluator to do a genotype
  // choose a diploid genotype
  // return it and set it so that decision tree cannot override

  //@TODO: fix this frequency to be sensible
  float local_min_allele_freq = smallest_allele_freq; // must choose a qual score relative to some frequency

  //eval.MultiAlleleGenotype(local_min_allele_freq,
  //    my_decision.eval_genotype.genotype_component,
  //    my_decision.eval_genotype.evaluated_genotype_quality,
  //    my_decision.eval_genotype.evaluated_variant_quality,
  //    parameters.my_eval_control.max_detail_level);

  my_decision.eval_genotype.genotype_already_set = true; // because we computed it here

  // Tell me what QUAL means (Is QUAL for a ref call or for a var call?) in case we won't show GT in the info tag.
  candidate_variant.variant.samples[sample_name]["QT"].push_back(convertToString(not my_decision.eval_genotype.IsReference()));

  // and I must also set for each allele so that the per-allele filter works
  for(unsigned int i_allele = 0; i_allele < my_decision.allele_identity_vector.size(); ++i_allele){
    my_decision.summary_info_vector[i_allele].variant_qual_score = my_decision.eval_genotype.evaluated_variant_quality;
    my_decision.summary_info_vector[i_allele].gt_quality_score = my_decision.eval_genotype.evaluated_genotype_quality;
  }

  // now that all the data has been gathered describing the variant, combine to produce the output
  my_decision.DecisionTreeOutputToVariant(candidate_variant, parameters, sample_index);
}

// Read and process records appropriate for this variant; positions are zero based
void Evaluator::StackUpOneVariant(const ExtendParameters &parameters, VariantCandidate &candidate_variant, const PositionInProgress& bam_position) {

  // Initialize random number generator for each stack -> ensure reproducibility
  RandSchrange RandGen(parameters.my_controls.RandSeed);

  read_stack.clear();  // reset the stack
  read_stack.reserve(2 * parameters.my_controls.downSampleCoverage);
  read_stack_n.reserve(parameters.my_controls.downSampleCoverage);
  read_stack_t.reserve(parameters.my_controls.downSampleCoverage);
  int read_counter = 0;
  int read_counter_n = 0;
  int read_counter_t = 0;

  for (Alignment* rai = bam_position.begin; rai != bam_position.end; rai = rai->next) {

    // Check global conditions to stop reading in more alignments
    if (rai->original_position > multiallele_window_start)
      break;

    if (rai->alignment.Position > multiallele_window_start)
      continue;

    if (rai->filtered)
      continue;

    if (rai->alignment.GetEndPosition() < multiallele_window_end)
      continue;

    // Get sample name without the numeric suffix
    string sample_name = candidate_variant.variant.sampleNames[rai->sample_index];
    sample_name = sample_name.substr(0, sample_name.find("."));

    // Reservoir Sampling
    if (read_stack.size() < 2 * (unsigned int)parameters.my_controls.downSampleCoverage) {
      read_counter++;
      read_stack.push_back(rai);
    }
    else {
      read_counter++;
      // produces a uniformly distributed test_position between [0, read_counter - 1]
      unsigned int test_position = ((double)RandGen.Rand() / ((double)RandGen.RandMax + 1.0)) * (double)read_counter;
      if (test_position < (unsigned int)parameters.my_controls.downSampleCoverage)
        read_stack[test_position] = rai;
    }

    if (sample_name == "tumor") {
      if (read_stack_t.size() < (unsigned int)parameters.my_controls.downSampleCoverage) {
        read_counter_t++;
        read_stack_t.push_back(rai);
      }
      else {
        read_counter_t++;
        // produces a uniformly distributed test_position between [0, read_counter - 1]
        unsigned int test_position = ((double)RandGen.Rand() / ((double)RandGen.RandMax + 1.0)) * (double)read_counter;
        if (test_position < (unsigned int)parameters.my_controls.downSampleCoverage)
          read_stack_t[test_position] = rai;
      }
    }

    if (sample_name == "normal") {
      if (read_stack_n.size() < (unsigned int)parameters.my_controls.downSampleCoverage) {
        read_counter_n++;
        read_stack_n.push_back(rai);
      }
      else {
        read_counter_n++;
        // produces a uniformly distributed test_position between [0, read_counter - 1]
        unsigned int test_position = ((double)RandGen.Rand() / ((double)RandGen.RandMax + 1.0)) * (double)read_counter;
        if (test_position < (unsigned int)parameters.my_controls.downSampleCoverage)
          read_stack_n[test_position] = rai;
      }
    }
  }
}

bool ProcessOneVariant (
  PersistingThreadObjects &thread_objects,
  VariantCallerContext& vc,
  VariantCandidate &candidate_variant,
  const PositionInProgress& bam_position
) {
  string sample_name = "";
  //if (sample_index >= 0) {sample_name = candidate_variant.variant.sampleNames[sample_index];}

  cerr << "processing " << sample_name << "\n";

  int chr_idx = vc.ref_reader->chr_idx(candidate_variant.variant.sequenceName.c_str());

  Evaluator eval(candidate_variant.variant);
  eval.SetupAllAlleles(*vc.parameters, *vc.global_context, *vc.ref_reader, chr_idx);
  eval.FilterAllAlleles(vc.parameters->my_controls.filter_variant, candidate_variant.variant_specific_params); // put filtering here in case we want to skip below entries

  // We read in one stack per multi-allele variant
  eval.StackUpOneVariant(*vc.parameters, candidate_variant, bam_position);

  cerr << "read stack sizes: " << eval.read_stack_n.size() << ", " << eval.read_stack_t.size() << ", " << eval.read_stack.size() << endl;
  if (eval.read_stack_n.empty() or eval.read_stack_t.empty()) {
    string reason = "NODATA";
    string sample_name = "";
    if (eval.read_stack_n.empty()) {
      cerr << "Nonfatal: No reads found in normal sample for " << candidate_variant.variant.sequenceName << "\t" << eval.multiallele_window_start << endl;
      sample_name = "normal.1";
    }
    if (eval.read_stack_t.empty()) {
      cerr << "Nonfatal: No reads found in tumor sample for " << candidate_variant.variant.sequenceName << "\t" << eval.multiallele_window_start << endl;
      sample_name = "tumor.1";
    }
    NullFilterReason(candidate_variant.variant, sample_name);
    AddFilterReason(candidate_variant.variant, reason, sample_name);
    SetFilteredStatus(candidate_variant.variant, true);
    candidate_variant.variant.samples[sample_name]["GT"].clear();
    candidate_variant.variant.samples[sample_name]["GT"].push_back("./.");
    return false;
  }


  // handle the unfortunate case in which we must try multiple alleles to be happy
  // try only ref vs alt allele here
  // leave ensemble in ref vs alt state

  // Estimate joint variant likelihood (the product of single-sample likelihoods obtained by piling reads from both samples)
  eval.SampleLikelihood(thread_objects, *vc.global_context, *vc.parameters, *vc.ref_reader, chr_idx, candidate_variant);
  cerr << candidate_variant.variant << endl;
  return true;

  // fill in quantities derived from predictions
  int num_hyp_no_null = eval.allele_identity_vector.size() + 1; // num alleles +1 for ref
  //eval.allele_eval.InitForInference(thread_objects, eval.read_stack, *vc.global_context, num_hyp_no_null, eval.allele_identity_vector);

  // output to variant
  // GlueOutputVariant(eval, candidate_variant, *vc.parameters, best_allele, sample_index);

  // output the inference results (MUQUAL, MUGT, MUGQ, etc.) if I turn on multi_min_allele_freq
  if (true or vc.parameters->program_flow.is_multi_min_allele_freq) {
    // MultiMinAlleleFreq(eval, candidate_variant, sample_index, vc.parameters->program_flow, vc.parameters->my_eval_control.max_detail_level);
  }

  cerr << "diagnostic: " << vc.parameters->program_flow.rich_json_diagnostic << endl;
  // test diagnostic output for this ensemble
  // if (vc.parameters->program_flow.rich_json_diagnostic & (!(eval.variant->isFiltered) | eval.variant->isHotSpot)) // look at everything that came through
  if (vc.parameters->program_flow.rich_json_diagnostic) // look at everything
    JustOneDiagnosis(eval, *vc.global_context, vc.parameters->program_flow.json_plot_dir, true);
  if (vc.parameters->program_flow.minimal_diagnostic & (!(eval.variant->isFiltered) | eval.variant->isHotSpot)) // look at everything that came through
    JustOneDiagnosis(eval, *vc.global_context, vc.parameters->program_flow.json_plot_dir, false);

  return true;
}





