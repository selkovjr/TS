/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */


#ifndef TENTATIVE_ALIGNMENT_H
#define TENTATIVE_ALIGNMENT_H


#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <vector>

#include "ExtendedReadInfo.h"

// use both strands for evaluating likelihood
#define ALL_STRAND_KEY -1

// no matter what, the read as called should be 1 in a million relative to the variant
// prevent some awkward moments if we divide by zero
#define MINIMUM_RELATIVE_OUTLIER_PROBABILITY 0.000001f

// handle auxiliary variables for one read's associated hypothesis evaluation
class TentativeAlignment {
  public:
    vector<string>         instance_of_read_by_state;       // this read, modified by each state of a variant
    vector<int>            state_spread;
    vector<bool>           same_as_null_hypothesis; // indicates whether a ref or alt hypothesis equals the read as called

    vector<vector<float> > basic_likelihoods; // likelihood given residuals at each test flow of the observation at that flow != likelihood of read

    // size number of hypotheses
    vector<float> responsibility; // how responsible this read is for a given hypothesis under the MAP: size number of hypotheses (including null=outlier)
    vector<float> log_likelihood; // sum over our test flows: logged to avoid under-flows
    vector<float> scaled_likelihood; // actual sum likelihood over test flows, rescaled to null hypothesis (as called), derived from log_likelihood
    float ll_scale; // local scaling factor for scaled likelihood as can't trust null hypothesis to be near data

    // intermediate allocations
    vector<float> tmp_prob_f;
    vector<double> tmp_prob_d;

    // useful hidden variables
    int strand_key;

    int heavy_tailed;

    bool success;

    // functions
    TentativeAlignment (){
      heavy_tailed = 3;  // t_5 degrees of freedom
      strand_key = 0;
      success = true;
      ll_scale = 0.0f;
    };
    void  CleanAllocate(int num_hyp);
    void  FillInPrediction(PersistingThreadObjects &thread_objects, const Alignment &my_read, const InputStructures &global_context);
    void  InitializeDerivedQualities();
    void  ComputeBasicLikelihoods();
    void  ComputeLogLikelihoods();
    void  ComputeLogLikelihoodsSum();
    void  JointLogLikelihood();
    void  ComputeScaledLikelihood();
    float ComputePosteriorLikelihood(const vector<float> &hyp_prob, float typical_prob);
    void  InitializeResponsibility();
    void  UpdateResponsibility(const vector<float > &hyp_prob, float typical_prob);
    void  UpdateRelevantLikelihoods();
    float ComputeLLDifference(int a_hyp, int b_hyp);
    int   MostResponsible();
    // HardOutlierClassifier is used to pre-filter out the outliers in a family.
    // Then it is reasonable to say that a functional family contains no outliers.
    bool  LocalOutlierClassifier(float typical_prob);

};


#endif // TENTATIVE_ALIGNMENT_H
