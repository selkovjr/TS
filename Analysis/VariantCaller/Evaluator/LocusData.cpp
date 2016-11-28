/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#include "LocusData.h"




// fill in predictions for each hypothesis and initialze test flows
void LocusData::FillInPredictionsAndTestFlows(PersistingThreadObjects &thread_objects, vector<const Alignment *>& read_stack,
    const InputStructures &global_context)
{
  //ion::FlowOrder flow_order(my_data.flow_order, my_data.flow_order.length());
  for (unsigned int i_read = 0; i_read < alignments.size(); i_read++) {
    alignments[i_read].FillInPrediction(thread_objects, *read_stack[i_read], global_context);
  }
}

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

