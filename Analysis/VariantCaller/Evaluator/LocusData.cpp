/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#include "LocusData.h"

void LocusData::TabulateFrequencies(VariantCandidate &candidate_variant, vector<const Alignment *>& read_stack) {
  depth = 0;
  normal_depth = 0;
  tumor_depth = 0;
  for (unsigned int i_read = 0; i_read < pileup.size(); i_read++) {
    depth++;
    string sample_name = candidate_variant.variant.sampleNames[read_stack[i_read]->sample_index];
    sample_name = sample_name.substr(0, sample_name.find("."));
    const string basecall = pileup[i_read].basecall;

    if (allele_count.count(basecall)) {
      allele_count[basecall] += 1;
    }
    else {
      allele_count[basecall] = 1;
    }

    if (sample_name == "normal") {
      normal_depth++;
      if (normal_allele_count.count(basecall)) {
        normal_allele_count[basecall] += 1;
      }
      else {
        normal_allele_count[basecall] = 1;
      }
    }

    if (sample_name == "tumor") {
      tumor_depth++;
      if (tumor_allele_count.count(basecall)) {
        tumor_allele_count[basecall] += 1;
      }
      else {
        tumor_allele_count[basecall] = 1;
      }
    }

    pileup[i_read].sample_index = read_stack[i_read]->sample_index;
  }

  for ( const auto &a: allele_count ) {
    allele_freq[a.first] = 1.0 * a.second / depth;
  }

  for ( const auto &a: normal_allele_count ) {
    normal_allele_freq[a.first] = 1.0 * a.second / normal_depth;
  }

  for ( const auto &a: tumor_allele_count ) {
    tumor_allele_freq[a.first] = 1.0 * a.second / tumor_depth;
  }
}

double LocusData::freq(string allele) {
  if (allele_freq.count(allele)) {
    return allele_freq[allele];
  }
  return 0.0;
}

double LocusData::freq(string allele, string sample_name) {
  if (sample_name == "normal" or sample_name == "tumor") {
    if (sample_name == "normal") {
      if (normal_allele_freq.count(allele)) {
        return normal_allele_freq[allele];
      }
      return 0.0;
    }

    if (sample_name == "tumor") {
      if (tumor_allele_freq.count(allele)) {
        return tumor_allele_freq[allele];
      }
      return 0.0;
    }
  }
  cerr << "invalid sample name '" << sample_name << "'\n";
  exit(-1);
}


void LocusData::FindValidIndexes() {
  valid_indexes.resize(0);
  // only loop over reads where variant construction worked correctly
  for (unsigned int i_read = 0; i_read < pileup.size(); i_read++) {
    if (pileup[i_read].success) {
      valid_indexes.push_back(i_read);
    }
  }
}

unsigned int LocusData::DetailLevel(void){
  return pileup.size();
}

