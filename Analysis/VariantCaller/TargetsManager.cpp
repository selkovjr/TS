/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     TargetsManager.cpp
//! @ingroup  VariantCaller
//! @brief    BED loader

#include "TargetsManager.h"

#include <stdlib.h>
#include <algorithm>
#include "BAMWalkerEngine.h"
//#include "ExtendParameters.h"


TargetsManager::TargetsManager()
{
}

TargetsManager::~TargetsManager()
{
}

bool CompareTargets(TargetsManager::UnmergedTarget *i, TargetsManager::UnmergedTarget *j)
{
  if (i->chr < j->chr)
    return true;
  if (i->chr > j->chr)
    return false;
  if (i->begin < j->begin)
    return true;
  if (i->begin > j->begin)
    return false;
  return i->end < j->end;
}

// Goals:
// * BED loading
//   - Support merged and unmerged BEDs (merge if needed, keep both versions)
//   - Support plain and detailed BEDs (track line, extra columns)
//   - Possibly support unsorted BEDs (low priority)
//   - If no BED provided, build a fake one covering entire reference
//   - Validation

void TargetsManager::Initialize(const ReferenceReader& ref_reader, const string& _targets)
{

  //
  // Step 1. Retrieve raw target definitions
  //

  list<UnmergedTarget>  raw_targets;

  if (not _targets.empty()) {
    LoadRawTargets(ref_reader, _targets, raw_targets);

  } else {
    for (int chr = 0; chr < ref_reader.chr_count(); ++chr) {
      raw_targets.push_back(UnmergedTarget());
      UnmergedTarget& target = raw_targets.back();
      target.begin = 0;
      target.end = ref_reader.chr_size(chr);
      target.chr = chr;
    }
  }

  //
  // Step 2. Sort raw targets and transfer to the vector
  //

  int num_unmerged = raw_targets.size();
  vector<UnmergedTarget*> raw_sort;
  raw_sort.reserve(num_unmerged);
  for (list<UnmergedTarget>::iterator I = raw_targets.begin(); I != raw_targets.end(); ++I)
    raw_sort.push_back(&(*I));
  sort(raw_sort.begin(), raw_sort.end(), CompareTargets);

  unmerged.reserve(num_unmerged);
  bool already_sorted = true;
  list<UnmergedTarget>::iterator I = raw_targets.begin();
  for (int idx = 0; idx < num_unmerged; ++idx, ++I) {
    if (raw_sort[idx] != &(*I) and already_sorted) {
      already_sorted = false;
      cerr << "TargetsManager: BED not sorted at position " << idx;
      cerr << " replaced " << I->name << ":" << I->chr << ":" << I->begin << "-" << I->end;
      cerr << " with " << raw_sort[idx]->name << ":" << raw_sort[idx]->chr << ":" << raw_sort[idx]->begin << "-" << raw_sort[idx]->end << endl;
    }
    unmerged.push_back(*raw_sort[idx]);
  }



  //
  // Step 3. Merge targets and link merged/unmerged entries
  //

  merged.reserve(num_unmerged);
  bool already_merged = true;
  for (int idx = 0; idx < num_unmerged; ++idx) {
    if (idx and merged.back().chr == unmerged[idx].chr and merged.back().end >= unmerged[idx].begin) {
      merged.back().end = max(merged.back().end, unmerged[idx].end);
      already_merged = false;
    } else {
      merged.push_back(MergedTarget());
      merged.back().chr = unmerged[idx].chr;
      merged.back().begin = unmerged[idx].begin;
      merged.back().end = unmerged[idx].end;
      merged.back().first_unmerged = idx;
    }
    unmerged[idx].merged = merged.size();
  }

  if (_targets.empty()) {
    cout << "TargetsManager: No targets file specified, processing entire reference" << endl;

  }
  else  {
    cout << "TargetsManager: Loaded targets file " << _targets << endl;

    cout << "TargetsManager: " << num_unmerged << " target(s)";
    if (not already_merged)
      cout << " (" << merged.size() << " after merging)";
    cout << endl;
    if (not already_sorted)
      cout << "TargetsManager: Targets required sorting" << endl;
  }
}




void TargetsManager::LoadRawTargets(const ReferenceReader& ref_reader, const string& bed_filename, list<UnmergedTarget>& raw_targets)
{
  FILE *bed_file = fopen(bed_filename.c_str(), "r");
  if (not bed_file) {
    cerr << "ERROR: Unable to open target file " << bed_filename << " : " << strerror(errno) << endl;
    exit(1);
  }

  char line[4096];
  char chr_name[4096];
  int begin;
  int end;
  char region_name[4096];
  int line_number = 0;

  while (fgets(line, 4096, bed_file)) {
    ++line_number;

    if (strncmp(line,"track",5) == 0) {
      // Parse track line if needed
      continue;
    }


    int num_fields = sscanf(line, "%s\t%d\t%d\t%s", chr_name, &begin, &end, region_name);
    if (num_fields == 0)
      continue;
    if (num_fields < 3) {
      cerr << "ERROR: Failed to parse target file line " << line_number << endl;
      exit(1);
    }

    raw_targets.push_back(UnmergedTarget());
    UnmergedTarget& target = raw_targets.back();
    target.begin = begin;
    target.end = end;
    target.chr = ref_reader.chr_idx(chr_name);
    if (num_fields > 3 and strcmp(region_name,".") != 0)
      target.name = region_name;

    if (target.chr < 0) {
      cerr << "ERROR: Target region " << target.name << " (" << chr_name << ":" << begin << "-" << end << ")"
        << " has unrecognized chromosome name" << endl;
      exit(1);
    }

    if (begin < 0 || end > ref_reader.chr_size(target.chr)) {
      cerr << "ERROR: Target region " << target.name << " (" << chr_name << ":" << begin << "-" << end << ")"
        << " is outside of reference sequence bounds ("
        << chr_name << ":0-" << ref_reader.chr_size(target.chr) << ")" << endl;
      exit(1);
    }
    if (end < begin) {
      cerr << "ERROR: Target region " << target.name << " (" << chr_name << ":" << begin << "-" << end << ")"
        << " has inverted coordinates" << endl;
      exit(1);
    }
  }

  fclose(bed_file);

  if (raw_targets.empty()) {
    cerr << "ERROR: No targets loaded from " << bed_filename
      << " after parsing " << line_number << " lines" << endl;
    exit(1);
  }
}

