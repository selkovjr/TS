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


// tool for combining items at differing log-levels
float log_sum(float a, float b) {
  float max_ab = max(a, b);
  float log_sum_val = max_ab + log(exp(a - max_ab) + exp(b - max_ab));
  return log_sum_val;
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

