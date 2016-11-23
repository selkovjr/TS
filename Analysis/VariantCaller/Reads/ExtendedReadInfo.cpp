/* Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     ExtendedReadInfo.cpp
//! @ingroup  VariantCaller
//! @brief    HP Indel detection

#include "ExtendedReadInfo.h"
#include "ion_util.h"
#include "ReferenceReader.h"
#include "BAMWalkerEngine.h"
#include "RandSchrange.h"
#include "MiscUtil.h"


// -------------------------------------------------------

// Sets the member variables ref_aln, seq_aln, pretty_aln, startSC, endSC
void UnpackAlignmentInfo(Alignment *rai)
{
  rai->left_sc = 0;
  rai->right_sc = 0;

  unsigned int num_query_bases = 0;
  bool match_found = false;

  for (vector<CigarOp>::const_iterator cigar = rai->alignment.CigarData.begin(); cigar != rai->alignment.CigarData.end(); ++cigar) {
    switch (cigar->Type) {
      case 'M':
      case '=':
      case 'X':
        match_found = true;
        rai->pretty_aln.append(cigar->Length, '|');
        num_query_bases += cigar->Length;
        break;

      case 'I':
        rai->pretty_aln.append(cigar->Length, '+');
        num_query_bases += cigar->Length;
        break;

      case 'S':
        num_query_bases += cigar->Length;
        if (match_found)
          rai->right_sc = cigar->Length;
        else
          rai->left_sc = cigar->Length;
        break;

      case 'D':
      case 'P':
      case 'N':
        rai->pretty_aln.append(cigar->Length, '-');
        break;
    }
  }
  // after possible trimming
  // rai->align_start = rai->alignment.Position;
  // rai->align_end = rai->alignment.GetEndPosition(false, true);

  // Basic alignment sanity check
  if (num_query_bases != rai->alignment.QueryBases.length()) {
    cerr << "WARNING in ExtendedReadInfo::UnpackAlignmentInfo: Invalid Cigar String in Read " << rai->alignment.Name << " Cigar: ";
    for (vector<CigarOp>::const_iterator cigar = rai->alignment.CigarData.begin(); cigar != rai->alignment.CigarData.end(); ++cigar)
      cerr << cigar->Length << cigar->Type;
    cerr << " Length of query string: " << rai->alignment.QueryBases.length() << endl;
    assert(num_query_bases == rai->alignment.QueryBases.length());
  }

}



// -------------------------------------------------------
// Unpacking read meta data and filtering read if this is not possible

void UnpackOnLoad(Alignment *rai, const InputStructures &global_context)
{
  // No need to waste time if the read is filtered
  if (rai->filtered)
    return;

  rai->is_reverse_strand = rai->alignment.IsReverseStrand();

  // Parse read name and run id

  if (not rai->alignment.Name.empty()) {
    rai->well_rowcol.resize(2);
  }

  // Populate read_bases (bases without rev-comp on reverse-mapped reads)

  rai->read_bases = rai->alignment.QueryBases;
  if (rai->is_reverse_strand)
    RevComplementInPlace(rai->read_bases);
  if (rai->read_bases.empty()){
    cerr << "WARNING: Ignoring length zero read " << rai->alignment.Name << endl;
    rai->filtered = true;
    return;
  }

  // Unpack alignment
  rai->pretty_aln.reserve(rai->alignment.QueryBases.length());

  UnpackAlignmentInfo(rai);
  if (rai->is_reverse_strand)
    rai->start_sc = rai->right_sc;
  else
    rai->start_sc = rai->left_sc;

  // Retrieve read group name
  if (not rai->alignment.GetTag("RG", rai->read_group)) {
    cerr << "WARNING: No read group found in read " << rai->alignment.Name << endl;
    // No big problem, we'll just have to solve the prefix like it's 2013!
    rai->read_group.clear();
  }
}

