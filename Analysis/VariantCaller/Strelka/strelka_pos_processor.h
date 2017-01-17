#ifndef STRELKA_POS_PROCESSOR
#define STRELKA_POS_PROCESSOR

#include "HandleVariant.h"
#include "SpliceVariantHypotheses.h"
#include "DecisionTreeData.h"
#include "Evaluator/TentativeAlignment.h"

#include "starling_common/gvcf_aggregator.hh"
#include "strelka/strelka_shared.hh" // for strelka_options

static auto_ptr<gvcf_aggregator> _gvcfer;

void pileup_read_segment (
  vector<TentativeAlignment> pileup,
  vector<const Alignment *> read_stack,
  const ExtendParameters &parameters,
  strelka_options &opt,
  const strelka_deriv_options &dopt,
  const ReferenceReader &ref_reader,
  const reference_contig_segment& ref,
  int chr_idx,
  int sample_no,
  VariantCandidate &candidate_variant
);

#endif
