/* Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     parameter.h
//! @ingroup  VariantCaller
//! @brief    Indel detection

#ifndef EXTENDPARAMETERS_H
#define EXTENDPARAMETERS_H

#include <string>
#include <vector>
#include "OptArgs.h"
#include "json/json.h"

using namespace std;


class BasicFilters {
  public:
    float min_allele_freq;
    int min_allele_count;

    float strand_bias_threshold;
    float strand_bias_pval_threshold;

    float min_quality_score;

    int min_cov;
    int min_cov_each_strand;
    int min_var_cov;

    BasicFilters() {
      min_allele_count =2;
      min_allele_freq = 0.2f;
      strand_bias_threshold = 0.8f;
      strand_bias_pval_threshold = 1.0f;

      min_cov = 3;
      min_cov_each_strand = 3;
      min_quality_score = 2.5f;
      min_var_cov = 0;
    };
};


class ControlCallAndFilters {
  public:
    int downSampleCoverage;
    int RandSeed;                  //!< Seed for random number generator to reservoir sample reads.

    // position bias probably should not be variant specific
    bool use_position_bias;
    float position_bias_ref_fraction;
    float position_bias;
    float position_bias_pval;

    // Dima's LOD filter
    bool use_lod_filter;
    float lod_multiplier;

    // tuning parameter for xbias
    //  float xbias_tune;
    float sbias_tune;

    BasicFilters filter_snps;
    BasicFilters filter_mnp;
    BasicFilters filter_hp_indel;

    ControlCallAndFilters();
    void SetOpts(OptArgs &opts, Json::Value& tvc_params);
    void CheckParameterLimits();
};

class ProgramControlSettings {
  public:
    // how we do things
    int nThreads;
    int nVariantsPerThread;
    int DEBUG;

    bool is_multi_min_allele_freq;
    vector<float> snp_multi_min_allele_freq;
    vector<float> mnp_multi_min_allele_freq;
    vector<float> indel_multi_min_allele_freq;

    ProgramControlSettings();
    void SetOpts(OptArgs &opts, Json::Value & pf_params);
    void CheckParameterLimits();
};

class ExtendParameters {
public:
  vector<string>    bams;
  string            fasta;                // -f --fasta-reference
  string            targets;              // -t --targets
  string            outputFile;
  string            blacklistFile;
  string            variantPriorsFile;
  string            postprocessed_bam;

  string            basecaller_version;
  string            tmap_version;

  // operation parameters
  bool        useDuplicateReads;        // -E --use-duplicate-reads
  int         useBestNAlleles;          // -n --use-best-n-alleles
  bool        allowIndels;              // -I --allow-indels
  bool        allowMNPs;                // -X --allow-mnps
  bool        allowComplex;             // -X --allow-complex
  int         maxComplexGap;
  bool        allowSNPs;                // -I --no-snps
  int         min_mapping_qv;           // -m --min-mapping-quality
  float       readMaxMismatchFraction;  // -z --read-max-mismatch-fraction
  int         read_snp_limit;           // -$ --read-snp-limit
  long double minAltFraction;           // -F --min-alternate-fraction
  long double minIndelAltFraction;      // Added by SU to reduce Indel Candidates for Somatic
  int         minAltCount;              // -C --min-alternate-count
  int         minAltTotal;              // -G --min-alternate-total
  int         minCoverage;              // -! --min-coverage
  int         mergeLookAhead;           // --merge-variant-lookahead
  bool        debug;                    // set if debuglevel >=1
  bool        multisample;              // multisample run


  OptArgs opts;
  ControlCallAndFilters my_controls;
  ProgramControlSettings program_flow;

  //Input files
  string outputDir;

  string sseMotifsFileName;
  bool sseMotifsProvided;

  string sampleName;
  string force_sample_name;

  string              params_meta_name;
  string              params_meta_details;

  // functions
  ExtendParameters(int argc, char** argv);

  bool ValidateAndCanonicalizePath(string &path);
  void SetupFileIO(OptArgs &opts, Json::Value& tvc_params);
  void SetFreeBayesParameters(OptArgs &opts, Json::Value& fb_params);
  void ParametersFromJSON(OptArgs &opts, Json::Value &tvc_params, Json::Value &fb_params, Json::Value &params_meta);
  void CheckParameterLimits();

};

template <class T>
bool CheckParameterLowerUpperBound(string identifier ,T &parameter, T lower_limit, T upper_limit) {
  bool is_ok = false;

  //cout << setw(35) << long_name_hyphens << " = " << setw(10) << value << " (integer, " << source << ")" << endl;
  cout << "Limit check parameter " << identifier << ": lim. "
     << lower_limit << " <= " << parameter << " <= lim. " << upper_limit << "? ";
  if (parameter < lower_limit) {
  cout << "Using " << identifier << "=" << lower_limit << " instead!";
    parameter = lower_limit;
  }
  else if (parameter > upper_limit) {
    cout << "Using " << identifier << "=" << upper_limit << " instead!";
    parameter = upper_limit;
  }
  else {
    cout << "OK!";
    is_ok = true;
  }
  cout << endl;
  return (is_ok);
}

template <class T>
bool CheckParameterLowerBound(string identifier ,T &parameter, T lower_limit) {
  bool is_ok = false;
  cout << "Limit check parameter " << identifier << ": lim. "
     << lower_limit << " <= " << parameter << "? ";
  if (parameter < lower_limit) {
  cout << "Using " << identifier << "=" << lower_limit << " instead!";
    parameter = lower_limit;
  }
    else {
    cout << "OK!";
    is_ok = true;
  }
  cout << endl;
  return (is_ok);
}

template <class T>
bool CheckParameterUpperBound(string identifier ,T &parameter, T upper_limit) {
  bool is_ok = false;
  cout << "Limit check parameter " << identifier << ": "
     << parameter << " <= lim. " << upper_limit << "? ";
  if (parameter > upper_limit) {
    cout << "Using " << identifier << "=" << upper_limit << " instead!";
    parameter = upper_limit;
  }
  else {
    cout << "OK!";
    is_ok = true;
  }
  cout << endl;
  return (is_ok);
}

bool CheckParameterStringContext(string identifier, string &parameter, const string &context, const string &default_value);

#endif // EXTENDPARAMETERS_H

