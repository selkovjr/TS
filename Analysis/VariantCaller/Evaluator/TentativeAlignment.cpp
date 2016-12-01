/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#include "TentativeAlignment.h"


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


void TentativeAlignment::CleanAllocate(int num_hyp) {
  // allocate my vectors here
  responsibility.assign(num_hyp, 0.0f);
  log_likelihood.assign(num_hyp, 0.0f);
  scaled_likelihood.assign(num_hyp, 0.0f);

  tmp_prob_f.assign(num_hyp, 0.0f);
  tmp_prob_d.assign(num_hyp, 0.0);

  basic_likelihoods.resize(num_hyp);
}

void TentativeAlignment::FillInPrediction(PersistingThreadObjects &thread_objects, const Alignment& my_read, const InputStructures &global_context) {
  // allocate everything here
  CleanAllocate(instance_of_read_by_state.size());
  if (my_read.is_reverse_strand)
    strand_key = 1;
  else
    strand_key = 0;
}


void TentativeAlignment::InitializeDerivedQualities() {

  InitializeResponsibility(); // depends on hypotheses

  // ComputeBasicLikelihoods(); // depends on residuals and sigma
  // compute log-likelihoods
  ComputeLogLikelihoods();  // depends on test flow(s)
}

void TentativeAlignment::InitializeResponsibility() {
  cerr << "InitializeResponsibility\n";
  responsibility[0] = 1.0f;  // everyone is an outlier until we trust you
  for (unsigned int i_hyp = 1; i_hyp < responsibility.size(); i_hyp++)
    responsibility[i_hyp] = 0.0f;
}



// responsibility depends on the relative global probability of the hypotheses and the likelihoods of the observations under each hypothesis
// divide the global probabilities into "typical" data points and outliers
// divide the variant probabilities into each hypothesis (summing to 1)
// treat the 2 hypothesis case to start with
void TentativeAlignment::UpdateResponsibility(const vector<float > &hyp_prob, float typical_prob) {
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

float TentativeAlignment::ComputePosteriorLikelihood(const vector<float > &hyp_prob, float typical_prob) {
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


void TentativeAlignment::UpdateRelevantLikelihoods() {
  // ComputeBasicLikelihoods();
  ComputeLogLikelihoods(); // automatically over relevant likelihoods
}

void TentativeAlignment::ComputeLogLikelihoodsSum() {
  for (unsigned int i_hyp=0; i_hyp<log_likelihood.size(); i_hyp++) {
    log_likelihood[i_hyp] = 0.0f;
    // Do real computation here
  }
}


// and again:  project onto the first component when correlation high

void TentativeAlignment::JointLogLikelihood() {
  for (unsigned int i_hyp =0; i_hyp<log_likelihood.size(); i_hyp++) {
    // Do real computation here
    float b_likelihood = 1.0f; // my_t.TDistOddN(res_projection,sigma_projection,skew_estimate);
    log_likelihood[i_hyp] = log(b_likelihood);
  }
}

void TentativeAlignment::ComputeLogLikelihoods() {
  //@TODO: magic numbers bad for thresholds
  // if ((fabs(delta_state.delta_correlation)<0.8f) || !use_correlated_likelihood) // suppress this feature for now
  if (true)
    ComputeLogLikelihoodsSum();
  else
    JointLogLikelihood();
  ComputeScaledLikelihood();
}

void TentativeAlignment::ComputeScaledLikelihood() {
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

float TentativeAlignment::ComputeLLDifference(int a_hyp, int b_hyp) {
  // difference in likelihoods between hypotheses
  return(fabs(log_likelihood[a_hyp] - log_likelihood[b_hyp]));
}

