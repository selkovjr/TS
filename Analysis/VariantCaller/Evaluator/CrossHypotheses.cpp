/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#include "CrossHypotheses.h"


//  control degrees of freedom for tradeoff in outlier resistance/sensitivity
float xTDistOddN(float res, float sigma, float skew, int half_n) {
  // skew t-dist one direction or the other
  float l_sigma;
  if (res>0.0f) {
    l_sigma = sigma*skew;
  } else {
    l_sigma = sigma/skew;
  }

  float x = res/l_sigma;
  float xx = x*x;
  float v = 2*half_n-1; // 1,3,5,7,...
  float my_likelihood = 1.0f/(3.14159f*sqrt(v));
  float my_factor = 1.0f/(1.0f+xx/v);

  for (int i_prod=0; i_prod<half_n; i_prod++) {
    my_likelihood *= my_factor;
  }
  for (int i_prod=1; i_prod<half_n; i_prod++) {
    my_likelihood *= (v+1.0f-2.0f*i_prod)/(v-2.0f*i_prod);
  }
  my_likelihood /= l_sigma;
  // account for skew
  float skew_factor = 2.0f*skew/(skew*skew+1.0f);
  my_likelihood *= skew_factor;
  return(my_likelihood);
}


void CrossHypotheses::CleanAllocate(int num_hyp) {
  // allocate my vectors here
  responsibility.assign(num_hyp, 0.0f);
  log_likelihood.assign(num_hyp, 0.0f);
  scaled_likelihood.assign(num_hyp, 0.0f);

  tmp_prob_f.assign(num_hyp, 0.0f);
  tmp_prob_d.assign(num_hyp, 0.0);

  basic_likelihoods.resize(num_hyp);
}

void CrossHypotheses::FillInPrediction(PersistingThreadObjects &thread_objects, const Alignment& my_read, const InputStructures &global_context) {
  // allocate everything here
  CleanAllocate(instance_of_read_by_state.size());
  if (my_read.is_reverse_strand)
    strand_key = 1;
  else
    strand_key = 0;
}


void CrossHypotheses::InitializeDerivedQualities() {

  InitializeResponsibility(); // depends on hypotheses

  // ComputeBasicLikelihoods(); // depends on residuals and sigma
  // compute log-likelihoods
  ComputeLogLikelihoods();  // depends on test flow(s)
}

void CrossHypotheses::InitializeResponsibility() {
  cerr << "InitializeResponsibility\n";
  responsibility[0] = 1.0f;  // everyone is an outlier until we trust you
  for (unsigned int i_hyp = 1; i_hyp < responsibility.size(); i_hyp++)
    responsibility[i_hyp] = 0.0f;
}



// responsibility depends on the relative global probability of the hypotheses and the likelihoods of the observations under each hypothesis
// divide the global probabilities into "typical" data points and outliers
// divide the variant probabilities into each hypothesis (summing to 1)
// treat the 2 hypothesis case to start with
void CrossHypotheses::UpdateResponsibility(const vector<float > &hyp_prob, float typical_prob) {
  cerr << "UpdateResponsibility()\n";

  if (!success){
    cout << "alert: fail to splice still called" << endl;
    InitializeResponsibility();
  }
  else {
    //  vector<double> tmp_prob(3);
    tmp_prob_d[0] = (1.0f - typical_prob) * scaled_likelihood[0];   // i'm an outlier
    for (unsigned int i_hyp = 1; i_hyp < scaled_likelihood.size(); i_hyp++)
      tmp_prob_d[i_hyp] = typical_prob * hyp_prob[i_hyp - 1] * scaled_likelihood[i_hyp];

    double ll_denom = 0.0;
    for (unsigned int i_hyp = 0; i_hyp<scaled_likelihood.size(); i_hyp++){
      ll_denom += tmp_prob_d[i_hyp];
    }

    for (unsigned int i_hyp=0; i_hyp<responsibility.size(); i_hyp++) {
      responsibility[i_hyp] = tmp_prob_d[i_hyp]/ll_denom;
      cerr << "UpdateResponsibility(): responsibility[" << i_hyp << "] = " << responsibility[i_hyp] << endl;
    }
  }
}

float CrossHypotheses::ComputePosteriorLikelihood(const vector<float > &hyp_prob, float typical_prob) {
  cerr << "ComputePosteriorLikelihood()\n";
  //  vector<float> tmp_prob(3);
  tmp_prob_f[0] = (1.0f-typical_prob)*scaled_likelihood[0];   // i'm an outlier
  for (unsigned int i_hyp=1; i_hyp<scaled_likelihood.size(); i_hyp++){
    tmp_prob_f[i_hyp] = typical_prob * hyp_prob[i_hyp-1] * scaled_likelihood[i_hyp];
  }
  float ll_denom = 0.0f;
  for (unsigned int i_hyp=0; i_hyp<scaled_likelihood.size(); i_hyp++) {
    ll_denom += tmp_prob_f[i_hyp];
  }
  return(log(ll_denom)+ll_scale);  // log-likelihood under current distribution, including common value of log-likelihood-scale
}


void CrossHypotheses::UpdateRelevantLikelihoods() {
  // ComputeBasicLikelihoods();
  ComputeLogLikelihoods(); // automatically over relevant likelihoods
}

void CrossHypotheses::ComputeLogLikelihoodsSum() {
  for (unsigned int i_hyp=0; i_hyp<log_likelihood.size(); i_hyp++) {
    log_likelihood[i_hyp] = 0.0f;
    // Do real computation here
  }
}


// and again:  project onto the first component when correlation high

void CrossHypotheses::JointLogLikelihood() {
  for (unsigned int i_hyp =0; i_hyp<log_likelihood.size(); i_hyp++) {
    // Do real computation here
    float b_likelihood = 1.0f; // my_t.TDistOddN(res_projection,sigma_projection,skew_estimate);
    log_likelihood[i_hyp] = log(b_likelihood);
  }
}

void CrossHypotheses::ComputeLogLikelihoods() {
  //@TODO: magic numbers bad for thresholds
  // if ((fabs(delta_state.delta_correlation)<0.8f) || !use_correlated_likelihood) // suppress this feature for now
  if (true)
    ComputeLogLikelihoodsSum();
  else
    JointLogLikelihood();
  ComputeScaledLikelihood();
}

void CrossHypotheses::ComputeScaledLikelihood() {
  //cout << "ComputeScaledLikelihood" << endl;
  //  scaled_likelihood.resize(log_likelihood.size());
  // doesn't matter who I scale to, as long as we scale together
  ll_scale = log_likelihood[0];
  for (unsigned int i_hyp =0; i_hyp<scaled_likelihood.size(); i_hyp++){
    if (log_likelihood[i_hyp] > ll_scale){
      ll_scale = log_likelihood[i_hyp];
    }
  }

  for (unsigned int i_hyp =0; i_hyp<scaled_likelihood.size(); i_hyp++){
    scaled_likelihood[i_hyp] = exp(log_likelihood[i_hyp]-ll_scale);
  }
  // really shouldn't happen, but sometimes does if splicing has gone awry
  // prevent log(0) events from happening in case we evaluate under weird circumstances
  scaled_likelihood[0] = max(scaled_likelihood[0], MINIMUM_RELATIVE_OUTLIER_PROBABILITY);
}

float CrossHypotheses::ComputeLLDifference(int a_hyp, int b_hyp) {
  // difference in likelihoods between hypotheses
  return(fabs(log_likelihood[a_hyp] - log_likelihood[b_hyp]));
}


int CrossHypotheses::MostResponsible(){
  int most_r = 0;
  float best_r = responsibility[0];
  for (unsigned int i_hyp = 1; i_hyp < responsibility.size(); i_hyp++){
    if (responsibility[i_hyp] > best_r){
      best_r = responsibility[i_hyp];
      most_r = i_hyp;
    }
  }
  return(most_r);
}

// Consider all the "binary" hypothesis hyp_0 (outlier) vs. hyp_i (i>0) with the prior P(hyp_i) = typical_prob.
// I claim the read is "not" an outlier, if I find any binary hypothesis pair hyp_0 vs. hyp_i s.t. log APP(hyp_0) < log APP(hyp_i),
// where APP is the posterior probability.
// Suppose that sigma is the default magic sigma (minimum_sigma_prior = 0.085, slope_sigma_prior = 0.0084).
// Example 1:
// hyp_0: 0-mer, hyp_1: 1-mer, and the measurement = 0-mer
// If outlier_probability < 0.00055 then we claim that it is not an outlier.
// => We need outlier_probability < 5.5E-4 to tolerate a 1-mer error at a single flow.
// Example 2:
// Suppose hyp_0: 0-mer, hyp_1: 2-mer, and the measurement = 0-mer
// With the default value of magic sigma, if outlier_probability < 1.1E-5 then we claim that it is not an outlier.
// => We need outlier_probability < 1.1E-5 to tolerate a 2-mer error at a single flow.
// Example 3:
// Suppose hyp_0: {0-mer, 0-mer}, hyp_1: {1-mer, 1-mer}, and the measurement = {0-mer, 0-mer}
// With the default value of magic sigma, if outlier_probability < 3.1E-7 then we claim that it is not an outlier.
// => We need outlier_probability < 3.1E-7 to tolerate 1-mer, 1-mer errors at two flows.
bool CrossHypotheses::LocalOutlierClassifier(float typical_prob){
  //Use CheckParameterLowerUpperBound(...) to prevent extremely low ol_prob instead
  //if(ol_prob < MINIMUM_RELATIVE_OUTLIER_PROBABILITY){
  //    ol_prob = MINIMUM_RELATIVE_OUTLIER_PROBABILITY;
  //    typical_prob = 1.0f - MINIMUM_RELATIVE_OUTLIER_PROBABILITY;
  //}
  bool is_outlier = true;

  // First check same_as_null_hypothesis
  // The read is not an outlier if any hypothesis not null is the same as the null hypothesis
  // This can prevent the case where we classify a hypothesis the same as null to be an outlier if typical_prob < 0.5.
  for(unsigned int i_hyp = 1; i_hyp < same_as_null_hypothesis.size(); ++i_hyp){
    if(same_as_null_hypothesis[i_hyp]){
      is_outlier = false;
      return is_outlier;
    }
  }

  float log_priori_not_ol = log(typical_prob);
  float log_posterior_likelihood_ol = log(1.0f - typical_prob) + log_likelihood[0];

  for(unsigned int i_hyp = 1; i_hyp < log_likelihood.size(); ++i_hyp){
    float log_posterior_likelihood_not_ol = log_priori_not_ol + log_likelihood[i_hyp];
    // If any hypothesis says "I am not an outlier", then the read is not an outlier read.
    if(log_posterior_likelihood_ol < log_posterior_likelihood_not_ol){
      is_outlier = false;
      return is_outlier;
    }
  }
  return is_outlier;
}

