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
    vector<TentativeAlignment> alignments;
    vector<int> valid_indexes;

    void FindValidIndexes(); // only loop over reads where we successfully filled in variants
    void ResetQualities();
    unsigned int DetailLevel(void);
};

#endif // SHORTSTACK_H
