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
    string basecall;
    int qscore;
    double error_prob;
    int pos_in_read;

    // useful hidden variables
    int strand_key;
    int sample_index;

    bool success;

    // functions
    TentativeAlignment (){
      sample_index = -1;
      strand_key = 0;
      success = true;
    };
};


#endif // TENTATIVE_ALIGNMENT_H
