/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#include "StackEngine.h"

// utility function for comparisons
bool compare_best_response(pair<int,float> a, pair<int,float> b){
  return (a.second>b.second);
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

  for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
    if (allele_eval.alignments[i_read].success) {
      // measure disruption

      mean_ll_delta += allele_eval.alignments[i_read].ComputeLLDifference(ref_hyp, alt_hyp);
      count++;
    }
  }
  mean_ll_delta /= (count + 0.01f);
  mean_ll_delta = 10.0f * mean_ll_delta / log(10.0f); // phred-scaled
}

int Evaluator::DetectBestMultiAllelePair(){
  int best_alt_ndx = 0; // forced choice with ref
  //@TODO: just get the plane off the ground
  //@TODO: do the top pair by responsibility
  vector< pair<int,float> > best_allele_test;
  int num_hyp_no_null = allele_eval.alignments[0].responsibility.size()-1;
  best_allele_test.resize(num_hyp_no_null); // null can never be a winner in "best allele" sweepstakes

  for (unsigned int i_alt=0; i_alt<best_allele_test.size(); i_alt++){
    best_allele_test[i_alt].first = i_alt;
    best_allele_test[i_alt].second = 0.0f;
  }

  // take responsibility
  for(unsigned int i_read = 0; i_read < read_id_.size(); ++i_read){
    int my_alt = read_id_[i_read];
    if (my_alt > -1){
      best_allele_test[my_alt].second += allele_eval.alignments[i_read].responsibility[my_alt + 1];
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

