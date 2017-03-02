/* Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     InputStructures.h
//! @ingroup  VariantCaller
//! @brief    HP Indel detection


#ifndef INPUTSTRUCTURES_H
#define INPUTSTRUCTURES_H

#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <pthread.h>
#include <Variant.h>
#include <VariantCaller/OrderedBAMWriter.h>

#include "LinearCalibrationModel.h"
#include "../ErrorMotifs.h"
#include "ExtendParameters.h"
#include "Realigner.h"

using namespace std;
using namespace BamTools;
using namespace ion;


class InputStructures;
class ReferenceReader;
class TargetsManager;
class BAMWalkerEngine;
class AlleleParser;
class OrderedVCFWriter;
class SampleManager;
class MetricsManager;
class IndelAssembly;

// ==============================================================================

struct VariantCallerContext {

  ExtendParameters *    parameters;                   //! Raw parameters, retrieved from command line and json
  InputStructures *     global_context;               //! Miscellaneous data structures
  ReferenceReader *     ref_reader;                   //! Memory-mapped reference reader
  TargetsManager *      targets_manager;              //! Manages target regions from BED file
  BAMWalkerEngine *     bam_walker;                   //! Manages traversing reference and retrieving reads covering each position
  AlleleParser *        candidate_generator;          //! Candidate variant generator
  OrderedBAMWriter *    bam_writer;                   //! Container for reads ensuring ordering
  OrderedVCFWriter *    vcf_writer;                   //! Sorting, threading friendly VCF writer
  MetricsManager *      metrics_manager;              //! Keeps track of metrics to output in tvc_metrics.json
  SampleManager*        sample_manager;               //! Tracks the sample used in multi-sample analysis
  IndelAssembly*        indel_assembly;               //! Outputs the indel_assembly.vcf for long indels

  pthread_mutex_t       bam_walker_mutex;             //! Mutex for state-altering bam_walker operations
  pthread_mutex_t       read_loading_mutex;           //! Mutex for raw read retrieval
  pthread_mutex_t       candidate_generation_mutex;   //! Mutex for candidate generation
  pthread_mutex_t       read_removal_mutex;           //! Mutex for removal/write of used reads
  pthread_cond_t        memory_contention_cond;       //! Conditional variable for preventing memory contention
  pthread_cond_t        alignment_tail_cond;          //! Conditional variable for blocking until the final alignments are processed

  int                   candidate_counter;            //! Number of candidates generated so far
  // int                   candidate_dot;                //! Number of candidates that will trigger printing next "."
  time_t                dot_time;                     //! Time that will trigger printing next '.'
};

// ==============================================================================

struct VariantSpecificParams {
  VariantSpecificParams() :
      min_allele_freq_override(false), min_allele_freq(0),
      strand_bias_override(false), strand_bias(0),
      strand_bias_pval_override(false), strand_bias_pval(0),

      min_coverage_override(false), min_coverage(0),
      min_coverage_each_strand_override(false), min_coverage_each_strand(0),
      min_variant_score_override(false), min_variant_score(0),
      data_quality_stringency_override(false), data_quality_stringency(0),

      position_bias_override(false), position_bias(0),
      position_bias_pval_override(false), position_bias_pval(0),

      hp_max_length_override(false), hp_max_length(0),
      filter_unusual_predictions_override(false), filter_unusual_predictions(0),
      filter_insertion_predictions_override(false), filter_insertion_predictions(0),
      filter_deletion_predictions_override(false), filter_deletion_predictions(0),
      sse_prob_threshold_override(false), sse_prob_threshold(0),
      black_strand('.') {}

  bool  min_allele_freq_override;
  float min_allele_freq;
  bool  strand_bias_override;
  float strand_bias;
  bool  strand_bias_pval_override;
  float strand_bias_pval;
  bool  min_coverage_override;
  int   min_coverage;
  bool  min_coverage_each_strand_override;
  int   min_coverage_each_strand;
  bool  min_variant_score_override;
  float min_variant_score;
  bool  data_quality_stringency_override;
  float data_quality_stringency;

  bool  position_bias_override;
  float position_bias;
  bool position_bias_pval_override;
  float position_bias_pval;

  bool  hp_max_length_override;
  int   hp_max_length;
  bool  filter_unusual_predictions_override;
  float filter_unusual_predictions;
  bool  filter_insertion_predictions_override;
  float filter_insertion_predictions;
  bool  filter_deletion_predictions_override;
  float filter_deletion_predictions;
  bool  sse_prob_threshold_override;
  float sse_prob_threshold;
  char  black_strand;
};

// ==============================================================================

struct VariantCandidate {
  VariantCandidate(vcf::VariantCallFile& initializer) : variant(initializer), position_upper_bound(0) {}
  vcf::Variant variant;
  vector<VariantSpecificParams> variant_specific_params;
  int position_upper_bound;   // Prevents SNP healing from accidentally reordering vcf records
};


// ==============================================================================
//Input Structures

class InputStructures {
public:
  int                       DEBUG;

  // Reusable objects
  TIonMotifSet              ErrorMotifs;

  InputStructures();
  void Initialize(ExtendParameters &parameters, const ReferenceReader& ref_reader, const SamHeader &bam_header);
  void read_error_motifs(string & fname){ErrorMotifs.load_from_file(fname.c_str());};
};

// ==============================================================================

// A collections of objects that are shared and reused throughout the execution of one tread
class PersistingThreadObjects {
public:

    PersistingThreadObjects(const InputStructures &global_context);
    ~PersistingThreadObjects() {};

    Realigner              realigner;             // realignment tool
};

#endif //INPUTSTRUCTURES_H
