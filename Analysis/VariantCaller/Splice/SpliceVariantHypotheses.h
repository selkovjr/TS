/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#ifndef SPLICEVARIANTHYPOTHESES_H
#define SPLICEVARIANTHYPOTHESES_H


#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <assert.h>
#include "stdlib.h"
#include "ctype.h"
#include "ClassifyVariant.h"
#include "ExtendedReadInfo.h"
#include "MiscUtil.h"

using namespace std;

class Evaluator;

bool SpliceVariantHypotheses(
  const Alignment &current_read,
  const Evaluator &eval,
  const LocalReferenceContext &local_context,
  PersistingThreadObjects &thread_objects,
  string &basecall,
  int &qscore,
  double &error_prob,
  int &position,
  bool &changed_alignment,
  const InputStructures
  &global_context,
  const ReferenceReader
  &ref_reader,
  int chr_idx
);


bool SpliceAddVariantAlleles(
  const Alignment &current_read,
  const string& pretty_alignment,
  const Evaluator &eval,
  const LocalReferenceContext &local_context,
  vector<string> &alignments,
  vector<string> &qualities,
  unsigned int pretty_idx, int DEBUG
);


void IncrementAlignmentIndices(const char aln_symbol, int &ref_idx, int &read_idx);

void DecrementAlignmentIndices(const char aln_symbol, int &ref_idx, int &read_idx);

string SpliceDoRealignement (PersistingThreadObjects &thread_objects, const Alignment &current_read, long variant_position,
                         bool &changed_alignment, int DEBUG, const ReferenceReader &ref_reader, int chr_idx);

#endif // SPLICEVARIANTHYPOTHESES_H
