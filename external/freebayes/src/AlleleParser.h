#ifndef ALLELEPARSER_H
#define ALLELEPARSER_H

#include <string>
#include <vector>
#include <map>
#include <deque>
#include <Variant.h>
#include "ExtendParameters.h"
#include "ReferenceReader.h"
#include "BAMWalkerEngine.h"
#include "SampleManager.h"
#include "InputStructures.h"

using namespace std;

class OrderedVCFWriter;

// ====================================================================

class AlleleDetails {
public:
  AlleleDetails() : type(ALLELE_UNKNOWN), chr(0), position(0),
      ref_length(0), length(0), minimized_prefix(0), repeat_boundary(0), hp_repeat_len (0),
      initialized(false), filtered(false),
      coverage(0), coverage_fwd(0), coverage_rev(0),
      samples(1) {}

  void add_observation(const Allele& observation, int sample_index, bool is_reverse_strand, int _chr, int num_samples, int read_count) {
    int q = observation.quality_string[0] - 33;

    if (not initialized) {
      type = observation.type;
      alt_sequence.append(observation.alt_sequence, observation.alt_length);
      position = observation.position;
      chr = _chr;
      ref_length = observation.ref_length;
      initialized = true;
    }
    if (sample_index < 0) return;

    coverage += read_count;
    if ((int)samples.size() != num_samples)
      samples.resize(num_samples);
    samples[sample_index].coverage += read_count;
    samples[sample_index].alt_q += q;
    if (is_reverse_strand) {
      coverage_rev += read_count;
      samples[sample_index].coverage_rev += read_count;
      samples[sample_index].alt_q_rev += q;
    } else {
      coverage_fwd += read_count;
      samples[sample_index].coverage_fwd += read_count;
      samples[sample_index].alt_q_fwd += q;
    }
  }

  void add_reference_observation(const Allele& observation, int sample_index, bool is_reverse_strand, int chr_idx_, int read_count) {
    int q = observation.quality_string[0] - 33;

    coverage += read_count;
    samples[sample_index].coverage += read_count;
    samples[sample_index].ref_q += q;
    if (is_reverse_strand) {
      coverage_rev += read_count;
      samples[sample_index].coverage_rev += read_count;
      samples[sample_index].ref_q_rev += q;
    }
    else {
      coverage_fwd += read_count;
      samples[sample_index].coverage_fwd += read_count;
      samples[sample_index].ref_q_fwd += q;
    }
  }

  void initialize_reference(long int _position, int num_samples) {
    type = ALLELE_REFERENCE;
    initialized = true;
    position = _position;
    ref_length = 1;
    length = 0;
    coverage = 0;
    coverage_fwd = 0;
    coverage_rev = 0;
    samples.clear();
    samples.resize(num_samples);
  }


  const char *type_str() {
    if (type == ALLELE_DELETION)    return "del";
    if (type == ALLELE_INSERTION)   return "ins";
    if (type == ALLELE_COMPLEX)     return "complex";
    if (type == ALLELE_SNP)         return "snp";
    if (type == ALLELE_MNP)         return "mnp";
    if (type == ALLELE_REFERENCE)   return "ref";
    return "???";
  }

  struct SampleData {
    SampleData() :
      coverage(0), coverage_fwd(0), coverage_rev(0),
      ref_q(0), ref_q_fwd(0), ref_q_rev(0),
      alt_q(0), alt_q_fwd(0), alt_q_rev(0)
    {}
    long int coverage;      // single-sample depth at this allele's position
    long int coverage_fwd;
    long int coverage_rev;
    long int ref_q;         // The sum of base quality scores across this sample's reads supporting reference allele
    long int ref_q_fwd;
    long int ref_q_rev;
    long int alt_q;         // The sum of base quality scores across this sample's reads supporting alt allele
    long int alt_q_fwd;
    long int alt_q_rev;
  };

  AlleleType              type;                 //! type of the allele
  string                  alt_sequence;         //! allele sequence
  int                     chr;                  //! chromosome
  long int                position;             //! position 0-based against reference
  unsigned int            ref_length;           //! allele length relative to the reference
  int                     length;               //! allele length reported in LEN tag
  int                     minimized_prefix;
  long int                repeat_boundary;      //! end of homopolymer or tandem repeat for indels
  long int                hp_repeat_len;        //! length of HP for filtering
  bool                    initialized;          //! is allele type info populated?
  bool                    filtered;             //! if true, do not report this allele as candidate
  long int                coverage;             //! total allele coverage (across samples)
  long int                coverage_fwd;         //! forward strand allele coverage (across samples)
  long int                coverage_rev;         //! reverse strand allele coverage (across samples)
  vector<SampleData>      samples;              //! per-sample coverages
};

// ====================================================================

class AllelePositionCompare {
public:
  bool operator()(const Allele& a, const Allele& b) const {
    if (a.position < b.position)
      return true;
    if (a.position > b.position)
      return false;
/*
    if (alternateLength < other.alternateLength)
      return true;
    if (alternateLength > other.alternateLength)
      return false;
    if (referenceLength < other.referenceLength)
      return true;
    if (referenceLength > other.referenceLength)
      return false;

    int strdelta = strncmp(alternateSequence,other.alternateSequence,min(alternateLength,other.alternateLength));
    return strdelta < 0;
*/

    int strdelta = strncmp(a.alt_sequence,b.alt_sequence,min(a.alt_length,b.alt_length));
    if (strdelta < 0)
      return true;
    if (strdelta > 0)
      return false;
    if (a.alt_length < b.alt_length)
      return true;
    if (a.alt_length > b.alt_length)
      return false;
    return a.ref_length < b.ref_length;
  }
};


// ====================================================================

class AlleleParser {
public:

   AlleleParser(const ExtendParameters& parameters, const ReferenceReader& ref_reader,
       const SampleManager& sample_manager, OrderedVCFWriter& vcf_writer);
  ~AlleleParser();

  //! basic filters to reject a read
  bool BasicFilters(Alignment& ra) const;

  //! Populates the allele specific data in the read Alignment object
  void UnpackReadAlleles(Alignment& ra) const;

  void GenerateCandidates(deque<VariantCandidate>& variant_candidates,
      list<PositionInProgress>::iterator& position_ticket, int& haplotype_length);

private:
  void MakeAllele(deque<Allele>& alleles, AlleleType type, long int pos, int length, const char *alt_sequence, const char *quality) const;

  void PileUpAlleles(int allowed_allele_types, int haplotype_length, bool scan_haplotype, list<PositionInProgress>::iterator& position_ticket);

  void PileUpAlleles(int pos, int haplotype_length, list<PositionInProgress>::iterator& position_ticket);

  void InferAlleleTypeAndLength(AlleleDetails& allele) const;

  long ComputeRepeatBoundary(const string& seq, int chr, long position, int max_size, long &hp_repeat_len) const;

  void GenerateCandidateVariant(deque<VariantCandidate>& variant_candidates,
      list<PositionInProgress>::iterator& position_ticket, int& haplotype_length);


  // operation parameters
  bool                        only_use_input_alleles_;
  bool                        process_input_positions_only_;
  bool                        use_duplicate_reads_;      // -E --use-duplicate-reads
  int                         use_best_n_alleles_;         // -n --use-best-n-alleles
  int                         max_complex_gap_;
  unsigned int                min_mapping_qv_;                    // -m --min-mapping-quality
  float                       read_max_mismatch_fraction_;  // -z --read-max-mismatch-fraction
  int                         read_snp_limit_;            // -$ --read-snp-limit
  long double                 min_alt_fraction_;  // -F --min-alternate-fraction
  long double                 min_indel_alt_fraction_; // Added by SU to reduce Indel Candidates for Somatic
  int                         min_alt_count_;             // -C --min-alternate-count
  int                         min_alt_total_;             // -G --min-alternate-total
  int                         min_coverage_;             // -! --min-coverage
  int                         allowed_allele_types_;
  int 			      merge_lookahead_;

  // data structures
  const ReferenceReader *     ref_reader_;
  const SampleManager *       sample_manager_;
  OrderedVCFWriter *          vcf_writer_;              //! Only used as Variant factory
  int                         num_samples_;

  typedef map<Allele,AlleleDetails,AllelePositionCompare>  pileup;
  pileup                      allele_pileup_;
  AlleleDetails               ref_pileup_;
  vector<long int>           coverage_by_sample_;
  //vector<char>                black_list_strand_;
  char                        black_list_strand_; // revert to 4.2
  int                         hp_max_lenght_override_value; //! if not zero then it overrides the maxHPLenght parameter in filtering
  float                       strand_bias_override_value;   //! if below zero then it overrides the strand_bias parameter in filtering

};

#endif
