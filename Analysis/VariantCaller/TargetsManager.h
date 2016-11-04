/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     TargetsManager.h
//! @ingroup  VariantCaller
//! @brief    BED loader


#ifndef TARGETSMANAGER_H
#define TARGETSMANAGER_H

#include <vector>
#include <list>
#include <string>
#include "ReferenceReader.h"

struct Alignment;

struct MergedTarget {
  int         chr;
  int         begin;
  int         end;
  int         first_unmerged;
};

class TargetsManager {
public:
  TargetsManager();
  ~TargetsManager();

  void Initialize(const ReferenceReader& ref_reader, const string& _targets);


  struct UnmergedTarget {
    int         chr;
    int         begin;
    int         end;
    string      name;
    int         merged;
  };


  void LoadRawTargets(const ReferenceReader& ref_reader, const string& bed_filename, list<UnmergedTarget>& raw_targets);

  vector<UnmergedTarget>  unmerged;
  vector<MergedTarget>    merged;
};



#endif //TARGETSMANAGER_H
