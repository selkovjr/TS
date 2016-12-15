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
#include "TentativeAlignment.h"
#include "ExtendParameters.h"

using namespace std;


// induced theories of the world
class LocusData {
public:
    vector<TentativeAlignment> pileup;
    double depth;
    double tumor_depth;
    double normal_depth;
    map <string, long> allele_count;
    map <string, long> tumor_allele_count;
    map <string, long> normal_allele_count;
    map <string, double> allele_freq;
    map <string, double> tumor_allele_freq;
    map <string, double> normal_allele_freq;
    vector<int> valid_indexes;

    void FindValidIndexes(); // only loop over reads where we successfully filled in variants
    void ResetQualities();
    unsigned int DetailLevel(void);

    void TabulateFrequencies(VariantCandidate &candidate_variant, vector<const Alignment *>& read_stack);
    double freq(string allele);
    double freq(string allele, string sample_name);
};

#endif // SHORTSTACK_H
