/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#include "LocusData.h"

void LocusData::ResetQualities() {
  // ! does not reset test flows or delta (correctly)
  for (unsigned int i_read = 0; i_read < alignments.size(); i_read++) {
    alignments[i_read].InitializeDerivedQualities();
  }
}


void LocusData::FindValidIndexes() {
  valid_indexes.resize(0);
  // only loop over reads where variant construction worked correctly
  for (unsigned int i_read = 0; i_read < alignments.size(); i_read++) {
    if (alignments[i_read].success) {
      valid_indexes.push_back(i_read);
    }
  }
}

unsigned int LocusData::DetailLevel(void){
  return alignments.size();
}

