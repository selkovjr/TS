/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#ifndef SHORTSTACK_H
#define SHORTSTACK_H

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <vector>
#include <map>

#include "api/BamReader.h"
#include "ExtendedReadInfo.h"
#include "CrossHypotheses.h"
#include "ExtendParameters.h"

using namespace std;


// induced theories of the world
class ShortStack{
public:
    vector<CrossHypotheses> my_hypotheses;
    vector<int> valid_indexes;

    void FindValidIndexes(); // only loop over reads where we successfully filled in variants
    void FillInPredictionsAndTestFlows(PersistingThreadObjects &thread_objects, vector<const Alignment *>& read_stack, const InputStructures &global_context);
    void ResetQualities();
    void PropagateTuningParameters(EnsembleEvalTuningParameters &my_params);
    void ResetRelevantResiduals();
    void UpdateRelevantLikelihoods();
    void ResetNullBias();
    unsigned int DetailLevel(void);

    float PosteriorFrequencyLogLikelihood(const vector<float> &hyp_freq, const vector<float> &prior_frequency_weight, float prior_log_normalization, float my_reliability, int strand_key);
    void MultiFrequencyFromResponsibility(vector<float> &hyp_freq, vector<float> &prior_frequency_weight, int strand_key);
    void UpdateResponsibility(const vector<float> &hyp_freq, float data_reliability);
};

#endif // SHORTSTACK_H
