/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */


#include "SpliceVariantHypotheses.h"

#include "StackEngine.h"

double prob(int qchar) {
  return 1.0 / (1.0 + pow(10.0, (qchar - 33) / 10.0));
}

bool SpliceVariantHypotheses (
  const Alignment &current_read,
  const Evaluator &eval,
  const LocalReferenceContext &local_context,
  PersistingThreadObjects &thread_objects,
  string &basecall,
  int &qscore,
  double &error_prob,
  int &position,
  bool &changed_alignment,
  const InputStructures &global_context,
  const ReferenceReader &ref_reader,
  int chr_idx
) {

  // Hypotheses: 1) Null; read as called 2) Reference Hypothesis 3-?) Variant Hypotheses
  vector<string> alignments;
  vector<string> qualities;

  alignments.resize(eval.allele_identity_vector.size() + 2);
  qualities.resize(eval.allele_identity_vector.size() + 2);

  // 1) Null hypothesis is read as called
  alignments[0] = current_read.read_bases;
  qualities[0] = current_read.read_qual;

  // Initialize hypotheses variables for splicing
  for (unsigned int i_hyp = 1; i_hyp < alignments.size(); i_hyp++) {
    alignments[i_hyp].clear();
    qualities[i_hyp].clear();
    alignments[i_hyp].reserve(current_read.alignment.QueryBases.length() + 20 + local_context.reference_allele.length());
    qualities[i_hyp].reserve(current_read.alignment.QueryBases.length() + 20 + local_context.reference_allele.length());
  }

  int read_idx = current_read.left_sc;
  int ref_idx  = current_read.alignment.Position;
  int read_idx_max = current_read.alignment.QueryBases.length() - current_read.right_sc;
  bool did_splicing = false;
  string pretty_alignment;
  changed_alignment = false;

  // do realignment of a small region around variant if desired
  if (eval.doRealignment) {
    pretty_alignment = SpliceDoRealignement(
      thread_objects,
      current_read,
      local_context.position0,
      changed_alignment,
      global_context.DEBUG,
      ref_reader,
      chr_idx
    );
    if (pretty_alignment.empty() and global_context.DEBUG > 0)
      cerr << "Realignment returned an empty string in read " << current_read.alignment.Name << endl;
  }

  if (pretty_alignment.empty()) {
    pretty_alignment = current_read.pretty_aln;
    changed_alignment = false;
  }

  // // Now fill in 2) and 3)

  // for (unsigned int pretty_idx = 0; pretty_idx < pretty_alignment.length(); pretty_idx++) {

  //   bool outside_of_window = ref_idx < eval.multiallele_window_start or ref_idx >= eval.multiallele_window_end;
  //   bool outside_ref_allele = (long)ref_idx < local_context.position0 or ref_idx >= (int)(local_context.position0 + local_context.reference_allele.length());

  //   // Basic sanity checks
  //   if (read_idx >= read_idx_max
  //       or  ref_idx > ref_reader.chr_size(chr_idx)
  //       or (ref_idx == ref_reader.chr_size(chr_idx) and pretty_alignment[pretty_idx] != '+')) {
  //     did_splicing = false;
  //     break;
  //   }

  //   // --- Splice ---
  //   if (ref_idx == local_context.position0 and !did_splicing and !outside_of_window) {
  //     // Add insertions before variant window
  //     while (pretty_idx < pretty_alignment.length() and pretty_alignment[pretty_idx] == '+') {
  //       for (unsigned int i_hyp = 1; i_hyp < alignments.size(); i_hyp++) {
  //         alignments[i_hyp].push_back(current_read.alignment.QueryBases[read_idx]);
  //         qualities[i_hyp].push_back(current_read.alignment.Qualities[read_idx]);
  //       }
  //       read_idx++;
  //       pretty_idx++;
  //     }
  //     did_splicing = SpliceAddVariantAlleles(current_read, pretty_alignment, eval, local_context, alignments, qualities, pretty_idx, global_context.DEBUG);
  //   } // --- ---

  //   // Have reference bases inside of window but outside of span of reference allele
  //   if (outside_ref_allele and !outside_of_window and pretty_alignment[pretty_idx] != '+') {
  //     for (unsigned int i_hyp = 1; i_hyp < alignments.size(); i_hyp++) {
  //       alignments[i_hyp].push_back(ref_reader.base(chr_idx, ref_idx));
  //       qualities[i_hyp].push_back(' ');
  //     }
  //   }

  //   // Have read bases as called outside of variant window
  //   if (outside_of_window and pretty_alignment[pretty_idx] != '-') {
  //     for (unsigned int i_hyp = 1; i_hyp < alignments.size(); i_hyp++) {
  //       alignments[i_hyp].push_back(current_read.alignment.QueryBases[read_idx]);
  //       qualities[i_hyp].push_back(current_read.alignment.Qualities[read_idx]);
  //     }
  //   }

  //   IncrementAlignmentIndices(pretty_alignment[pretty_idx], ref_idx, read_idx);
  // } // end of for loop over extended pretty alignment

  // // Check whether the whole reference allele fit
  // // It seems that with primer trimming ion TVC, many a read throw this warning
  // if (ref_idx < (int)(local_context.position0 + local_context.reference_allele.length())) {
  //   did_splicing = false;
  //   if (global_context.DEBUG > 0)
  //     cerr << "Warning in Splicing: Reference allele "<< local_context.reference_allele << " did not fit into read " << current_read.alignment.Name << endl;
  // }

  // if (did_splicing) {
  //   // --- Add soft clipped bases to the right of the alignment and reverse complement ---
  //   for (unsigned int i_hyp = 1; i_hyp < alignments.size(); i_hyp++) {
  //     alignments[i_hyp] += current_read.alignment.QueryBases.substr(current_read.alignment.QueryBases.length() - current_read.right_sc, current_read.right_sc);
  //     qualities[i_hyp] += current_read.alignment.Qualities.substr(current_read.alignment.QueryBases.length() - current_read.right_sc, current_read.right_sc);

  //     if (current_read.is_reverse_strand) {
  //       RevComplementInPlace(alignments[i_hyp]);
  //       RevInPlace(qualities[i_hyp]);
  //     }

  //     if (global_context.DEBUG > 0) {
  //       cerr << "hypothesis[" << i_hyp << "]: " << alignments[i_hyp] << endl;
  //       cerr << " qualities[" << i_hyp << "]: " << qualities[i_hyp] << endl;
  //     }
  //   }
  // }

  // // // Check for non-ACGT bases in hypotheses
  // // bool valid_bases = true;
  // // for (unsigned int i_hyp=0; i_hyp<alignments.size(); i_hyp++) {
  // //   unsigned int iBase = 0;
  // //   while (iBase < alignments[i_hyp].length() and valid_bases) {
  // //     if (alignments[i_hyp].at(iBase) == 'A' or alignments[i_hyp].at(iBase) == 'C' or
  // //         alignments[i_hyp].at(iBase) == 'G' or alignments[i_hyp].at(iBase) == 'T')
  // //       iBase++;
  // //     else
  // //       valid_bases = false;
  // //   }
  // // }
  // // if (not valid_bases){
  // //   cerr << "Non-Fatal ERROR in Splicing for " << local_context.contigName << ":" << local_context.position0 + 1
  // //     << ": Read Hypotheses for " << current_read.alignment.Name << " contain non-ACGT characters." << endl;
  // //   did_splicing = false;
  // // }

  // // --- Fail safe for hypotheses and verbose
  // if (!did_splicing) {
  //   for (unsigned int i_hyp = 1; i_hyp < alignments.size(); i_hyp++) {
  //     alignments[i_hyp] = alignments[0];
  //     qualities[i_hyp] = qualities[0];
  //   }
  //   if (global_context.DEBUG > 1) {
  //     cerr << "Failed to splice " << local_context.reference_allele << "->";
  //     for (unsigned int i_alt = 0; i_alt < eval.allele_identity_vector.size(); i_alt++) {
  //       cerr << eval.allele_identity_vector[i_alt].altAllele;
  //       if (i_alt < eval.allele_identity_vector.size() - 1)
  //         cerr << ",";
  //     }
  //     cerr << " into read " << current_read.alignment.Name << endl;
  //   }
  // }

  // else if (global_context.DEBUG > 1 and alignments[0] != alignments[1]) {
  //   cerr << "Spliced " << local_context.reference_allele << " -> ";
  //   for (unsigned int i_alt = 0; i_alt < eval.allele_identity_vector.size(); i_alt++) {
  //     cerr << eval.allele_identity_vector[i_alt].altAllele;
  //     if (i_alt < eval.allele_identity_vector.size() - 1)
  //       cerr << ",";
  //   }
  //   cerr << " into ";
  //   if (current_read.is_reverse_strand) cerr << "reverse ";
  //   else cerr << "forward ";
  //   cerr << "strand read " << current_read.alignment.Name << endl;
  //   cerr << "- Read as called: " << alignments[0] << endl;
  //   cerr << "- Reference Hyp.: " << alignments[1] << endl;
  //   for (unsigned int i_hyp = 2; i_hyp < alignments.size(); i_hyp++)
  //     cerr << "- Variant Hyp. " << (i_hyp - 1) << ": " << alignments[i_hyp] << endl;
  // }

  if (global_context.DEBUG > 1) cerr << "alignments: " << eval.allele_identity_vector.size() << endl;

  char qual;
  if (global_context.DEBUG > 1) {
    cerr << "equal to ref\n";
    cerr << "  allele: " << eval.allele_identity_vector[eval.allele_identity_vector.size() - 1].altAllele << endl;
  }

  string read_buf;
  string qual_buf;
  if (current_read.is_reverse_strand) {
    read_buf = alignments[0].substr(current_read.right_sc, alignments[0].length() - current_read.left_sc);
    qual_buf = qualities[0].substr(current_read.right_sc, alignments[0].length() - current_read.left_sc);
    RevComplementInPlace(read_buf);
    RevInPlace(qual_buf);
  }
  else {
    read_buf = alignments[0].substr(current_read.left_sc, alignments[0].length() - current_read.right_sc);
    qual_buf = qualities[0].substr(current_read.left_sc, alignments[0].length() - current_read.right_sc);
  }

  for (unsigned int pretty_idx = 0; pretty_idx < pretty_alignment.length(); pretty_idx++) {
    if (pretty_alignment[pretty_idx] == '-') {
      read_buf.insert(pretty_idx, "N");
      qual_buf.insert(pretty_idx, "N");
    }
    if (pretty_alignment[pretty_idx] == '+') {
      read_buf.erase(pretty_idx, 1);
      qual_buf.erase(pretty_idx, 1);
    }
  }

  // not sure whether this test is necessary
  if (local_context.position0 >= current_read.align_start and local_context.position0 <= current_read.align_end) {
    position = local_context.position0 - current_read.align_start;
    basecall = read_buf[position];
    qual = qual_buf[position];
  }
  else {
    basecall = "N";
    qual = ' ';
  }

#if 1
  // This is a rough-and-ready pileup viewer. Mind the fixed offset from
  // current_read.alignment.Position -- it can result in extremely long lines
  // if set too far upstream or std::length_error if set downstream of the
  // leftmost read start.
  cerr << (current_read.is_reverse_strand ? "(R)" : "(F)") << " " << current_read.alignment.Position << ": " << string(current_read.alignment.Position - 7723408, ' ') << read_buf << " (" << current_read.left_sc << ", " << current_read.right_sc << ")\n";
  cerr <<  string(current_read.alignment.Position - 7723408 + 13, ' ') << pretty_alignment << endl;
#endif

  qscore = qual - 33;
  error_prob = prob(qual);
  // return did_splicing;
  return true;
};

// -------------------------------------------------------------------

void IncrementAlignmentIndices(const char aln_symbol, int &ref_idx, int &read_idx) {
  switch (aln_symbol) {
    case ('-'):
      ref_idx++;
      break;
    case ('+'):
    case (' '):
    case ('|'):
      read_idx++;
      if (aln_symbol != '+')
        ref_idx++;
      break;
  }
}

void DecrementAlignmentIndices(const char aln_symbol, int &ref_idx, int &read_idx) {
  switch (aln_symbol) {
    case ('-'):
      ref_idx--;
      break;
    case ('+'):
    case (' '):
    case ('|'):
      read_idx--;
      if (aln_symbol != '+')
        ref_idx--;
      break;
  }
}


// -------------------------------------------------------------------

// This function is used when insertions count towards reference index before them.
bool SpliceAddVariantAlleles(
  const Alignment &current_read,
  const string& pretty_alignment,
  const Evaluator &eval,
  const LocalReferenceContext &local_context,
  vector<string> &alignments,
  vector<string> &qualities,
  unsigned int pretty_idx_orig,
  int DEBUG
) {

  // Splice reference Hypothesis
  alignments[1] += local_context.reference_allele;
  // qualities[1] += string(local_context.reference_allele.length(), 0xb0); // to make gaps visible while debugging
  qualities[1] += string(local_context.reference_allele.length(), '0'); // without a good reason

  for (unsigned int i_hyp = 2; i_hyp < alignments.size(); ++i_hyp) {
    int my_allele_idx = i_hyp - 2;

    // Special SNP splicing to not accidentally split HPs in the presence of insertions at start of HP
    if (eval.allele_identity_vector[my_allele_idx].status.isSNP) {
      int shifted_position = 0;
      unsigned int splice_idx = alignments[i_hyp].length();
      unsigned int pretty_idx = pretty_idx_orig;
      alignments[i_hyp] += local_context.reference_allele;
      qualities[i_hyp] += string(local_context.reference_allele.length(), 0xdb);
      // move left if there are insertions of the same base as the reference hypothesis base
      while (pretty_idx > 0 and pretty_alignment[pretty_idx - 1] == '+' and splice_idx > 0
          and current_read.alignment.QueryBases[splice_idx - 1] == local_context.reference_allele[0]) {
        pretty_idx--;
        splice_idx--;
        shifted_position++;
      }
      if (DEBUG > 1 and shifted_position > 0) {
        // printouts
        cerr << "Shifted splice position by " << shifted_position << " in " << current_read.alignment.Name
          << " " << local_context.position0 << local_context.reference_allele
          << "->" << eval.allele_identity_vector[my_allele_idx].altAllele << endl;
        cerr << alignments[i_hyp] << endl;
      }
      alignments[i_hyp][splice_idx] = eval.allele_identity_vector[my_allele_idx].altAllele[0];
      if (current_read.is_reverse_strand) {
        qualities[i_hyp][splice_idx] = current_read.read_qual[current_read.read_qual.length() - splice_idx - 1];
        qualities[i_hyp][splice_idx] = current_read.read_qual[current_read.read_qual.length() - splice_idx - 1];
      }
      else {
        qualities[i_hyp][splice_idx] = current_read.read_qual[splice_idx];
      }
    }
    else { // Default splicing
      alignments[i_hyp] += eval.allele_identity_vector[my_allele_idx].altAllele;
      alignments[i_hyp] += string(eval.allele_identity_vector[my_allele_idx].altAllele.length(), '*');
    }
  } // end looping over hypotheses
  return true;
}



// -------------------------------------------------------------------

string SpliceDoRealignement (
  PersistingThreadObjects &thread_objects,
  const Alignment &current_read,
  long variant_position,
  bool &changed_alignment,
  int DEBUG,
  const ReferenceReader &ref_reader,
  int chr_idx
) {

  // We do not allow any clipping since we align a short substring
  thread_objects.realigner.SetClipping(0, true);
  string new_alignment;


  // --- Get index positions at snp variant position
  int read_idx = current_read.left_sc;
  int ref_idx  = current_read.alignment.Position;
  unsigned int pretty_idx = 0;

  while (pretty_idx < current_read.pretty_aln.length() and ref_idx < variant_position) {
    IncrementAlignmentIndices(current_read.pretty_aln[pretty_idx], ref_idx, read_idx);
    pretty_idx++;
  }
  if (DEBUG > 1)
    cerr << "Computed variant position as (read, ref, pretty) " << read_idx << " " << ref_idx << " " << pretty_idx << endl;

  if (pretty_idx >= current_read.pretty_aln.length()
      or ref_idx  >= ref_reader.chr_size(chr_idx)
      or read_idx >= (int)current_read.alignment.QueryBases.length() - current_read.right_sc)
    return new_alignment;

  // --- Get small sequence context for very local realignment ------------------------
  int min_bases = 5;

  // Looking at alignment to the left of variant position to find right place to cut sequence
  int read_left = read_idx;
  int ref_left  = ref_idx;
  unsigned int pretty_left = pretty_idx;
  bool continue_looking = pretty_idx > 0;

  while (continue_looking) {
    pretty_left--;
    DecrementAlignmentIndices(current_read.pretty_aln[pretty_left], ref_left, read_left);

    // Stopping criterion
    if (pretty_left < 1) {
      continue_looking = false;
      break;
    }
    if (ref_idx - ref_left < min_bases)
      continue_looking = true;
    else {
      // make sure to start with a matching base and don't split large HPs
      if (current_read.pretty_aln[pretty_left] != '|'
          or (ref_reader.base(chr_idx,ref_left+1) == ref_reader.base(chr_idx,ref_left)))
        continue_looking = true;
      else
        continue_looking = false;
    }
  }
  if (DEBUG > 1)
    cerr << "Computed left realignment window as (read, ref, pretty) " << read_left << " " << ref_left << " " << pretty_left << endl;


  // Looking at alignment to the right to find right place to cut sequence
  int read_right = read_idx;
  int ref_right  = ref_idx;
  unsigned int pretty_right = pretty_idx;
  continue_looking = pretty_idx < current_read.pretty_aln.length() - 1;

  while (continue_looking) {
    IncrementAlignmentIndices(current_read.pretty_aln[pretty_right], ref_right, read_right);
    pretty_right++;
    // Stopping criterion (half open interval)
    if (pretty_right >= current_read.pretty_aln.length()
        or ref_right >= ref_reader.chr_size(chr_idx)) {
      continue_looking = false;
      break;
    }
    if (ref_right - ref_idx < min_bases)
      continue_looking = true;
    else {
      // make sure to stop with a matching base and don't split large HPs
      if (current_read.pretty_aln[pretty_right - 1] != '|'
          or (ref_reader.base(chr_idx,ref_right - 1) == ref_reader.base(chr_idx,ref_right)))
        continue_looking = true;
      else
        continue_looking = false;
    }
  }
  if (DEBUG > 1)
    cerr << "Computed right realignment window as (read, ref, pretty) " << read_right << " " << ref_right << " " << pretty_right << endl;
  // Put in some sanity checks for alignment boundaries found...


  // --- Realign -------------------------
  unsigned int start_position_shift;
  vector<CigarOp>    new_cigar_data;
  vector<MDelement>  new_md_data;

  // printouts
  if (DEBUG > 1) {
    thread_objects.realigner.verbose_ = true;
    cerr << "Realigned " << current_read.alignment.Name << " from " << endl;
  }
  if (read_left >= read_right and ref_left >= ref_right) {
    if (DEBUG > 1)
      cerr << "ERROR: realignment window has zero size! " << endl;
    return new_alignment;
  }

  string old_alignment = current_read.pretty_aln.substr(pretty_left, pretty_right - pretty_left);
  thread_objects.realigner.SetSequences(
    current_read.alignment.QueryBases.substr(read_left, read_right - read_left),
    ref_reader.substr(chr_idx, ref_left, ref_right - ref_left),
    old_alignment,
    true
  );

  if (!thread_objects.realigner.computeSWalignment(new_cigar_data, new_md_data, start_position_shift)) {
    if (DEBUG > 1)
      cerr << "ERROR: realignment failed! " << endl;
    return new_alignment;
  }

  // --- Fuse realigned partial sequence back into pretty_aln string
  new_alignment = current_read.pretty_aln;
  if (old_alignment == thread_objects.realigner.pretty_aln()) {
    changed_alignment = false;
  }
  else {
    new_alignment.replace(pretty_left, (pretty_right - pretty_left), thread_objects.realigner.pretty_aln());
    changed_alignment = true;
  }
  if (DEBUG > 1 and changed_alignment) {
    cerr << "new alignment: " << new_alignment << endl;
    cerr << "new_cigar_data: ";
    for (vector<CigarOp>::iterator iter = new_cigar_data.begin(); (iter != new_cigar_data.end()); ++iter) {
      cerr << iter->Length << string(1, iter->Type);
    }
    cerr << endl;
  }
  return new_alignment;
}

