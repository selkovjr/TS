/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#include "StackEngine.h"

// utility function for comparisons
bool compare_best_response(pair<int,float> a, pair<int,float> b){
  return (a.second>b.second);
}


void HypothesisStack::DefaultValues()
{
  try_alternatives = true;
}

// if try_alternatives and my_params.try_few_restart_freq, then I skip the pure
// frequencies of snp or mnp alt alleles, also skip the pure frequency of ref
// allele if "all" alt alleles are snp or mnp
void HypothesisStack::AllocateFrequencyStarts(int num_hyp_no_null, vector<AlleleIdentity> &allele_identity_vector){
  // int num_hyp = 2; // ref + alt, called doesn't count as a "start"
  //int num_start = num_hyp_no_null + 1;
  vector<float> try_me(num_hyp_no_null);
  float safety_zero = 0.0f;

  // I try at most (1 + num_hyp_no_null) frequencies: one uniform allele freq and num_hyp_no_null pure freq
  try_hyp_freq.reserve(1 + num_hyp_no_null);

  // Always try uniform hyp_freq
  try_me.assign(num_hyp_no_null, 1.0f / (float) num_hyp_no_null);
  try_hyp_freq.push_back(try_me);

  // -----------------------------------
  // For now, try_few_restart_freq means don't try alternative
  if(my_params.try_few_restart_freq){
    try_alternatives = false;
  }
  // -----------------------------------

  // try pure frequencies for the alleles
  if(try_alternatives){
    // try pure frequencies for alt alleles
    for(int i_hyp = 1; i_hyp < num_hyp_no_null; ++i_hyp){
      int i_alt = i_hyp - 1;
      if(my_params.try_few_restart_freq and (allele_identity_vector[i_alt].ActAsSNP() or allele_identity_vector[i_alt].ActAsMNP())){
        // skip the pure frequency of a snp or mnp alt allele.
        continue;
      }
      try_me.assign(num_hyp_no_null, safety_zero / float(num_hyp_no_null - 1));
      try_me[i_hyp] = 1.0f - safety_zero;
      try_hyp_freq.push_back(try_me);
    }
    // try the pure frequency of the ref allele if we try at least one pure frequency of alt allele
    if(try_hyp_freq.size() > 1){
      try_me.assign(num_hyp_no_null, safety_zero / float(num_hyp_no_null - 1));
      try_me[0] = 1.0f - safety_zero;
      try_hyp_freq.push_back(try_me);
    }
  }

  ll_record.assign(try_hyp_freq.size(), 0.0f);
}

void HypothesisStack::InitForInference(PersistingThreadObjects &thread_objects, vector<const Alignment *>& read_stack, const InputStructures &global_context, int num_hyp_no_null, vector<AlleleIdentity> &allele_identity_vector) {
  PropagateTuningParameters(num_hyp_no_null); // sub-objects need to know
  // predict given hypotheses per read
  total_theory.FillInPredictionsAndTestFlows(thread_objects, read_stack, global_context);
  total_theory.FindValidIndexes();
  // how many alleles?
  AllocateFrequencyStarts(num_hyp_no_null, allele_identity_vector);
}


void HypothesisStack::ExecuteInference() {
  // now with unrestrained inference
  ExecuteExtremeInferences();
  // set up our filter
  // @TODO: Now with fast scan, we can always do it
}

void HypothesisStack::SetAlternateFromMain() {
  // make sure things are populated

  //  ref_state = cur_state;
  //  var_state = cur_state;
}



// try altenatives to the main function to see if we have a better fit somewhere else
void HypothesisStack::ExecuteExtremeInferences() {
  /*
  for (unsigned int i_start=0; i_start < try_hyp_freq.size(); i_start++){
    ll_record[i_start] = ExecuteOneRestart(try_hyp_freq[i_start], my_params.max_detail_level);
  }
  //TriangulateRestart();
  RestoreFullInference(); // put total_theory to be consistent with whoever won
  */
}

// subobjects need to know their tuning
void HypothesisStack::PropagateTuningParameters(int num_hyp_no_null) {
  total_theory.PropagateTuningParameters(my_params);
  // number of pseudo-data points at no bias
}



// tool for combining items at differing log-levels
float log_sum(float a, float b) {
  float max_ab = max(a, b);
  float log_sum_val = max_ab + log(exp(a - max_ab) + exp(b - max_ab));
  return log_sum_val;
}


void GenotypeByIntegral(vector<float> genotype_interval, int &genotype_call, float &quasi_phred_quality_score){
  //@TODO: do as paired and sort, for clean code
  // best zone = call
  unsigned int best_call = 0;
  float best_val = genotype_interval[0];
  unsigned int worst_call = 0;
  float worst_val = genotype_interval[0];
  for (unsigned int i_geno = 1; i_geno < genotype_interval.size(); i_geno++) {
    if (best_val < genotype_interval[i_geno]) {
      best_val = genotype_interval[i_geno];
      best_call = i_geno;
    }
    if (worst_val > genotype_interval[i_geno]) {
      worst_val = genotype_interval[i_geno];
      worst_call = i_geno;
    }
  }
  float middle_val = 0.0f;
  for (unsigned int i_geno = 0; i_geno < genotype_interval.size(); i_geno++) {
    if ((i_geno != worst_call) & (i_geno != best_call)) {
      middle_val = genotype_interval[i_geno];
    }
  }
  // most likely interval
  genotype_call = best_call;

  // quality score
  float log_alternative = middle_val + log(1 + exp(worst_val - middle_val)); // total mass on alternative intervals
  float log_all = best_val + log(1 + exp(log_alternative - best_val)); // total mass on all intervals

  // output
  quasi_phred_quality_score = 10 * (log_all - log_alternative) / log(10); // -10*log(error), base 10

  if (isnan(quasi_phred_quality_score)) {
    cout << "Warning: quality score NAN "  << endl;
    quasi_phred_quality_score = 0.0f;
  }
}

bool RejectionByIntegral(vector<float> dual_interval, float &reject_status_quality_score){
  // reject ref = quality of rejection call
  bool is_ref = false;
  float log_ref = 0.0f;
  int top = 0;
  int bottom = 1;
  if (dual_interval[1]>dual_interval[0]) {
    // if reference, how strongly can we reject the outer interval
    log_ref = dual_interval[0];
    is_ref = false;
    top=1;
    bottom=0;
  }
  else {
    // if var, how strongly can we reject ref
    log_ref = dual_interval[1];
    is_ref = true;
    top=0;
    bottom=1;
  }

  float log_all = dual_interval[top]+ log(1+exp(dual_interval[bottom]-dual_interval[top]));
  reject_status_quality_score = 10 * (log_all - log_ref) / log(10); // how much mass speaks against the pure reference state
  if (isnan(reject_status_quality_score)) {
    cout << "Warning: reject ref score NAN " << endl;
    reject_status_quality_score = 0.0f;
  }
  return is_ref;
}


// evidence for i_allele vs ref
void Evaluator::ScanSupportingEvidence(float &mean_ll_delta,  int i_allele) {

  mean_ll_delta = 0.0f;
  int count = 0;
  int ref_hyp = 1;
  int alt_hyp = i_allele + 2;  // alt_alleles = 0->n not counting ref >or< null, therefore alt-allele 0 = 2

  for (unsigned int i_read = 0; i_read < allele_eval.total_theory.my_hypotheses.size(); i_read++) {
    if (allele_eval.total_theory.my_hypotheses[i_read].success) {
      // measure disruption

      mean_ll_delta += allele_eval.total_theory.my_hypotheses[i_read].ComputeLLDifference(ref_hyp, alt_hyp);
      count++;
    }
  }
  mean_ll_delta /= (count + 0.01f);
  mean_ll_delta = 10.0f * mean_ll_delta / log(10.0f); // phred-scaled
}

// Hard classify each read using its responsibility
// read_id_[i] = -1 means the i-th read is classified as an outlier.
// read_id_[i] = 0 means the i-th read is classified as ref.
// read_id_[i] = 1 means the i-th read is classified as the variant allele 1, and so on.
void Evaluator::ApproximateHardClassifierForReads() {
  cerr << "ApproximateHardClassifierForReads()\n";
  read_id_.clear();
  strand_id_.clear();
  dist_to_left_.clear();
  dist_to_right_.clear();

  read_id_.assign(read_stack.size(), -1);
  strand_id_.assign(read_stack.size(), false);
  dist_to_left_.assign(read_stack.size(), -1);
  dist_to_right_.assign(read_stack.size(), -1);

  int position0 = variant->position - 1; // variant->position 1-base: vcflib/Variant.h

  for (unsigned int i_read = 0; i_read < read_stack.size(); ++i_read) {
    // compute read_id_
    if (allele_eval.total_theory.my_hypotheses[i_read].success){
      read_id_[i_read] = allele_eval.total_theory.my_hypotheses[i_read].MostResponsible() - 1; // -1 = null, 0 = ref , ...
    }
    else{
      cerr << i_read << " is outlier" << endl;
      read_id_[i_read] = -1; // failure = outlier
    }

    // not an outlier
    if(read_id_[i_read] > -1){
      //fprintf(stdout, "position0 =%d, read_stack[i_read]->align_start = %d, read_stack[i_read]->align_end = %d, read_stack[i_read]->left_sc = %d, read_stack[i_read]->right_sc = %d\n", (int)position0, (int)read_stack[i_read]->align_start, (int)read_stack[i_read]->align_end, (int)read_stack[i_read]->left_sc, (int)read_stack[i_read]->right_sc);
      //fprintf(stdout, "dist_to_left[%d] = =%d, dist_to_right[%d] = %d\n", (int)i_read, (int)(position0 - read_stack[i_read]->align_start), (int)i_read, (int)(read_stack[i_read]->align_end - position0));
      dist_to_left_[i_read] = position0 - read_stack[i_read]->align_start;
      assert ( dist_to_left_[i_read] >=0 );
      dist_to_right_[i_read] = read_stack[i_read]->align_end - position0;
      assert ( dist_to_right_[i_read] >=0 );
    }
    //compute strand_id_
    strand_id_[i_read] = not (read_stack[i_read]->is_reverse_strand);
  }

  is_hard_classification_for_reads_done_ = true;
}


int Evaluator::DetectBestMultiAllelePair(){
  int best_alt_ndx = 0; // forced choice with ref
  //@TODO: just get the plane off the ground
  //@TODO: do the top pair by responsibility
  vector< pair<int,float> > best_allele_test;
  int num_hyp_no_null = allele_eval.total_theory.my_hypotheses[0].responsibility.size()-1;
  best_allele_test.resize(num_hyp_no_null); // null can never be a winner in "best allele" sweepstakes

  for (unsigned int i_alt=0; i_alt<best_allele_test.size(); i_alt++){
    best_allele_test[i_alt].first = i_alt;
    best_allele_test[i_alt].second = 0.0f;
  }

  if(not is_hard_classification_for_reads_done_){
    ApproximateHardClassifierForReads();
  }
  // take responsibility
  for(unsigned int i_read = 0; i_read < read_id_.size(); ++i_read){
    int my_alt = read_id_[i_read];
    if (my_alt > -1){
      best_allele_test[my_alt].second += allele_eval.total_theory.my_hypotheses[i_read].responsibility[my_alt + 1];
    } // otherwise count for nothing
  }

  // pick my pair of best alleles
  sort(best_allele_test.begin(), best_allele_test.end(), compare_best_response);
  // not-null choices
  diploid_choice[0] = best_allele_test[0].first; // index of biggest weight
  diploid_choice[1] = best_allele_test[1].first; // index of second-biggest weight
  // problematic cases:
  // 2 alleles & ref, ref + 1 allele zero, want ref as the comparison
  // all ref, don't care about the second allele for genotype?
  // all zero implies what?
  //cout << best_allele_test.at(0).first << "\t" << best_allele_test.at(0).second << "\t" << best_allele_test.at(1).first << "\t" << best_allele_test.at(1).second << endl;
  if(diploid_choice[0]==0)
    best_alt_ndx = diploid_choice[1]-1;
  else
    best_alt_ndx = diploid_choice[0]-1;

  // sort as final step to avoid best_alt_ndx reflecting a worse allele
  sort(diploid_choice.begin(),diploid_choice.end()); // now in increasing allele order as needed

  is_detect_best_multi_allele_pair_done_ = true;
  return best_alt_ndx;
}

