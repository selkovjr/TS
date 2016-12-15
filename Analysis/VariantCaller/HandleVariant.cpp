/* Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     HandleVariant.cpp
//! @ingroup  VariantCaller
//! @brief    HP Indel detection

#include "HandleVariant.h"
#include "SpliceVariantHypotheses.h"
#include "DecisionTreeData.h"

#define starling_streams_base int

#include "blt_common/adjust_joint_eprob.hh"
#include "blt_common/nonref_test_call.hh"
#include "blt_common/position_nonref_2allele_test.hh"
#include "blt_common/ref_context.hh"
#include "blt_common/snp_pos_info.hh"
#include "blt_util/bam_seq.hh"
#include "blt_util/depth_stream_stat_range.hh"
#include "blt_util/seq_util.hh"
#include "starling_common/gvcf_aggregator.hh"
#include "starling_common/gvcf_locus_info.hh"
#include "starling_common/pos_basecall_buffer.hh"
#include "starling_common/starling_ref_seq.hh" // get_starling_ref_seq()
#include "strelka/strelka_shared.hh" // for strelka_options
//#include "starling_common/starling_pos_processor_base.hh"

std::auto_ptr<gvcf_aggregator> _gvcfer;

/// get max-min bounds in which reads can be realigned:
static known_pos_range get_realignment_range(const pos_t pos, const stage_data& sdata) {
  const unsigned head_offset(sdata.get_stage_id_shift(STAGE::HEAD));
  const unsigned buffer_offset(sdata.get_stage_id_shift(STAGE::READ_BUFFER));
  const unsigned post_offset(sdata.get_stage_id_shift(STAGE::POST_ALIGN));
  assert(buffer_offset>head_offset);
  assert(post_offset>buffer_offset);

  const pos_t min_pos(std::max(static_cast<pos_t>(0), pos-static_cast<pos_t>(post_offset-buffer_offset)));
  const pos_t max_pos(pos+1+(buffer_offset-head_offset));
  return known_pos_range(min_pos, max_pos);
}

static void write_snp_prefix_info_file (
  const std::string& seq_name,
  const pos_t output_pos,
  const char ref,
  const unsigned n_used_calls,
  const unsigned n_unused_calls,
  std::ostream& os
) {
  os << "nonref_2allele_test\t" << seq_name << "\t"
    << output_pos << "\t"
    << n_used_calls << "\t"
    << n_unused_calls << "\t"
    << ref;
}

static void write_bsnp_diploid_allele (
  const blt_options& client_opt,
  // const blt_streams& client_io,
  const std::string& seq_name,
  const pos_t output_pos,
  const char ref,
  const unsigned n_used_calls,
  const unsigned n_unused_calls,
  const snp_pos_info& good_pi,
  const diploid_genotype& dgt,
  const unsigned hpol = 0
) {
  // std::ostream& os(*client_io.bsnp_diploid_allele_osptr());
  std::ostream& os(cout);

  os << "diploid_genotype_allele\t";
  write_snp_prefix_info_file(seq_name, output_pos, ref, n_used_calls, n_unused_calls, os);
  os << "\t";
  write_diploid_genotype_allele(client_opt, good_pi, dgt, os, hpol);
  os << "\n";
}

struct extra_position_data {
  /// stores the column of basecalls actually used for snp-calling after the
  /// mismatch density filter and other quality filters have been applied:
  snp_pos_info good_pi;

  /// stores information on approximate qscore reductions implemented to represent
  /// site-specific basecalling dependency, note this is not applied in somatic
  /// calling:
  std::vector<float> dependent_eprob;
};


// Removed read_id_counter
struct sample_info {
  sample_info (
    const starling_options& opt,
    const unsigned report_size,
    const unsigned knownref_report_size
  ):
  // indel_buff(opt.max_indel_size),
  // read_buff(ricp),
  // sample_opt(opt),
  // isync_default(indel_buff, estdepth_buff, sample_opt),
  // indel_sync_ptr(&isync_default),
  ss(report_size),
  used_ss(report_size),
  ssn(knownref_report_size),
  used_ssn(knownref_report_size)
  // wav(opt)
  {}

  // indel_synchronizer &indel_sync() { return *indel_sync_ptr; }

  // indel_buffer indel_buff;
  pos_basecall_buffer bc_buff;
  // starling_read_buffer read_buff;
  // depth_buffer estdepth_buff; // provide an early estimate of read depth before realignment.

  // starling_sample_options sample_opt;

  // indel_synchronizer isync_default;
  // indel_synchronizer* indel_sync_ptr;

  depth_stream_stat_range ss;
  depth_stream_stat_range used_ss;
  depth_stream_stat_range ssn;
  depth_stream_stat_range used_ssn;

  // regional basecall average windows:
  // win_avgs wav;

  // keep a single copy of this struct to reuse for every site to lower alloc costs:
  extra_position_data epd;
};


static void report_counts(const snp_pos_info& pi, const unsigned n_unused_calls, const pos_t output_pos, std::ostream& os) {
  unsigned base_count[N_BASE];

  for (unsigned i(0); i < N_BASE; ++i) base_count[i] = 0;

  BOOST_FOREACH(const base_call& bc, pi.calls) {
    assert(bc.base_id != BASE_ID::ANY);
    base_count[bc.base_id]++;
  }

  os << "counts\t" << output_pos << '\t';
  for (unsigned i(0); i < N_BASE; ++i) {
    os << base_count[i] << '\t';
  }
  os << n_unused_calls << '\n';
}

/*
void Evaluator::SampleLikelihood (
  PersistingThreadObjects &thread_objects,
  const InputStructures &global_context,
  const ExtendParameters &parameters,
  const ReferenceReader &ref_reader,
  int chr_idx,
  VariantCandidate &candidate_variant
) {
  bool changed_alignment;
  unsigned int  num_valid_reads = 0;
  unsigned int  num_realigned = 0;
  int  num_hyp_no_null = allele_identity_vector.size() + 1; // num alleles +1 for ref
  // generate null+ref+nr.alt hypotheses per read in the case of do_multiallele_eval
  allele_eval.alignments.resize(read_stack.size());

  cerr << "evaluating " << allele_eval.alignments.size() << " hypotheses\n";
  for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
    // --- New splicing function ---
    allele_eval.alignments[i_read].success =
      SpliceVariantHypotheses(
        *read_stack[i_read],
        *this,
        seq_context,
        thread_objects,
        allele_eval.alignments[i_read].basecall,
        allele_eval.alignments[i_read].qscore,
        allele_eval.alignments[i_read].error_prob,
        allele_eval.alignments[i_read].pos_in_read,
        changed_alignment,
        global_context,
        ref_reader,
        chr_idx
    );

    if (allele_eval.alignments[i_read].success){
      num_valid_reads++;
      if (changed_alignment)
        num_realigned++;
    }
  }
  cerr << "num_valid_reads: " << num_valid_reads << endl;
  cerr << "num_realigned: " << num_realigned << endl;

  // Check how many reads had their alignment modified
  std::ostringstream my_info;
  my_info.precision(4);
  if (doRealignment and num_valid_reads > 0) {
    cerr << "checking modified reads\n";
    float frac_realigned = (float)num_realigned / (float)num_valid_reads;
    // And re-do splicing without realignment if we exceed the threshold
    if (frac_realigned > parameters.my_controls.filter_variant.realignment_threshold) {
      my_info << "SKIPREALIGNx" << frac_realigned;
      doRealignment = false;
      for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
        allele_eval.alignments[i_read].success =
          SpliceVariantHypotheses(
            *read_stack[i_read],
            *this,
            seq_context,
            thread_objects,
            allele_eval.alignments[i_read].basecall,
            allele_eval.alignments[i_read].qscore,
            allele_eval.alignments[i_read].error_prob,
            allele_eval.alignments[i_read].pos_in_read,
            changed_alignment,
            global_context,
            ref_reader,
            chr_idx
          );
      }
    }
    else {
      my_info << "REALIGNEDx" << frac_realigned;
    }
    info_fields.push_back(my_info.str());
  }
  cerr << my_info.str() << endl;

  // Tabulate frequencies
  allele_eval.TabulateFrequencies(candidate_variant, read_stack);

  // Compute the log likelihood of joint samples
  double snv_likelihood = 0;
  for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
    double p = 0;
    const string basecall = allele_eval.alignments[i_read].basecall;
    double e = allele_eval.alignments[i_read].error_prob;
    cerr << allele_eval.alignments[i_read].sample_index << ": " << basecall << endl;
    for (const string &true_base: {"A", "C", "G", "T"}) {
      double p_b_given_a = basecall == true_base ? 1 - e : e / 3;
      cerr << "  P(" << basecall << "|" << true_base << ") = " << p_b_given_a << endl;
      cerr << "  f = " << allele_eval.freq(basecall) << endl;
      p += p_b_given_a * allele_eval.freq(basecall);
    }
    snv_likelihood += log10(p);
  }

  // Compute the log likelihood of normal sample
  double snv_likelihood_n = 0;
  for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
    string sample_name = candidate_variant.variant.sampleNames[allele_eval.alignments[i_read].sample_index];
    sample_name = sample_name.substr(0, sample_name.find("."));
    if (sample_name != "normal") continue;

    double p = 0;
    const string basecall = allele_eval.alignments[i_read].basecall;
    double e = allele_eval.alignments[i_read].error_prob;
    cerr << allele_eval.alignments[i_read].sample_index << ": " << basecall << endl;
    for (const string &true_base: {"A", "C", "G", "T"}) {
      double p_b_given_a = basecall == true_base ? 1 - e : e / 3;
      cerr << "  P(" << basecall << "|" << true_base << ") = " << p_b_given_a << endl;
      cerr << "  f = " << allele_eval.freq(basecall) << endl;
      p += p_b_given_a * allele_eval.freq(basecall);
    }
    snv_likelihood_n += log10(p);
  }

  // Compute the log likelihood of tumor sample
  double snv_likelihood_t = 0;
  for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
    string sample_name = candidate_variant.variant.sampleNames[allele_eval.alignments[i_read].sample_index];
    sample_name = sample_name.substr(0, sample_name.find("."));
    if (sample_name != "tumor") continue;

    double p = 0;
    const string basecall = allele_eval.alignments[i_read].basecall;
    double e = allele_eval.alignments[i_read].error_prob;
    cerr << allele_eval.alignments[i_read].sample_index << ": " << basecall << endl;
    for (const string &true_base: {"A", "C", "G", "T"}) {
      double p_b_given_a = basecall == true_base ? 1 - e : e / 3;
      cerr << "  P(" << basecall << "|" << true_base << ") = " << p_b_given_a << endl;
      cerr << "  f = " << allele_eval.freq(basecall) << endl;
      p += p_b_given_a * allele_eval.freq(basecall);
    }
    snv_likelihood_t += log10(p);
  }

  cerr << "total: " << allele_eval.alignments.size() << endl;
  for ( const auto &a: allele_eval.tumor_allele_count ) {
    cerr << "tumor " << a.first << " -> " << a.second << " (" << allele_eval.freq(a.first, "tumor") << ")\n";
  }
  for ( const auto &a: allele_eval.normal_allele_count ) {
    cerr << "normal " << a.first << " -> " << a.second << " (" << allele_eval.freq(a.first, "normal") << ")\n";
  }
  for (const string &base: {"A", "C", "G", "T"}) {
    cerr << base << " -> " << allele_eval.freq(base) << endl;
  }

  cerr << "joint sample likelihood: " << snv_likelihood << endl;
  cerr << "normal sample likelihood: " << snv_likelihood_n << endl;
  cerr << "tumor sample likelihood: " << snv_likelihood_t << endl;
}
*/


void Evaluator::Strelka (
  PersistingThreadObjects &thread_objects,
  const InputStructures &global_context,
  const ExtendParameters &parameters,
  const ReferenceReader &ref_reader,
  int chr_idx,
  VariantCandidate &candidate_variant
) {

  // Initialize priors and other options
  strelka_options opt;
  opt.samtools_ref_seq_file = "/data1/selkov_workdir/data/reference/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
  opt.is_samtools_ref_set = true;

  opt.bam_seq_name = candidate_variant.variant.sequenceName;

  // Strelka's realignment window for SNPs is 20 bases on either side
  opt.user_report_range.begin_pos = max(0, multiallele_window_start - 20);
  opt.user_report_range.is_begin_pos = true;
  opt.user_report_range.end_pos = multiallele_window_end + 20;
  opt.user_report_range.is_end_pos = true;

  opt.user_genome_size = parameters.my_controls.genome_size;
  opt.is_user_genome_size = true;

  opt.shared_site_error_rate = 5.0e-7;

  opt.is_counts = true;
  opt.is_compute_hapscore = true;
  opt.is_compute_VQSRmetrics = true;
  opt.is_compute_VQSRmetrics = true;
  bool _is_variant_windows(opt.variant_windows.size());
  std::set<pos_t> _variant_print_pos;

  reference_contig_segment ref;
  get_starling_ref_seq(opt, ref);
  const strelka_deriv_options dopt(opt, ref); // initializes priors
  cerr << "dopt range size: " << dopt.report_range.size() << endl;

  // pileup
  bool changed_alignment;
  unsigned int  num_valid_reads = 0;
  unsigned int  num_realigned = 0;
  int  num_hyp_no_null = allele_identity_vector.size() + 1; // num alleles +1 for ref
  // generate null+ref+nr.alt hypotheses per read in the case of do_multiallele_eval
  allele_eval.alignments.resize(read_stack.size());

  cerr << "evaluating " << allele_eval.alignments.size() << " hypotheses\n";
  for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
    // --- New splicing function ---
    allele_eval.alignments[i_read].success =
      SpliceVariantHypotheses(
        *read_stack[i_read],
        *this,
        seq_context,
        thread_objects,
        allele_eval.alignments[i_read].basecall,
        allele_eval.alignments[i_read].qscore,
        allele_eval.alignments[i_read].error_prob,
        allele_eval.alignments[i_read].pos_in_read,
        changed_alignment,
        global_context,
        ref_reader,
        chr_idx
    );

    if (allele_eval.alignments[i_read].success){
      num_valid_reads++;
      if (changed_alignment)
        num_realigned++;
    }
  }
  cerr << "num_valid_reads: " << num_valid_reads << endl;
  cerr << "num_realigned: " << num_realigned << endl;

  // Check how many reads had their alignment modified
  std::ostringstream my_info;
  my_info.precision(4);
  if (doRealignment and num_valid_reads > 0) {
    cerr << "checking modified reads\n";
    float frac_realigned = (float)num_realigned / (float)num_valid_reads;
    // And re-do splicing without realignment if we exceed the threshold
    if (frac_realigned > parameters.my_controls.filter_variant.realignment_threshold) {
      my_info << "SKIPREALIGNx" << frac_realigned;
      doRealignment = false;
      for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
        allele_eval.alignments[i_read].success =
          SpliceVariantHypotheses(
            *read_stack[i_read],
            *this,
            seq_context,
            thread_objects,
            allele_eval.alignments[i_read].basecall,
            allele_eval.alignments[i_read].qscore,
            allele_eval.alignments[i_read].error_prob,
            allele_eval.alignments[i_read].pos_in_read,
            changed_alignment,
            global_context,
            ref_reader,
            chr_idx
          );
      }
    }
    else {
      my_info << "REALIGNEDx" << frac_realigned;
    }
    info_fields.push_back(my_info.str());
  }
  cerr << my_info.str() << endl;

  // Fill in position info for strelka. The following section is roughly
  // equivalent to a portion of pileup_read_segment() for ps.type == MATCH.
  //
  int sample_no = 1;
  bool is_tier1 = true;
  long ref_pos = candidate_variant.variant.position;

  // sample_info &sif(sample(sample_no));
  const unsigned knownref_report_size(get_ref_seq_known_size(ref, dopt.report_range));
  sample_info sif(opt, dopt.report_range.size(), knownref_report_size);

  // Data structure defining parameters for a single site to be used for writing in gvcf_aggregator
  site_info _site_info = *(new site_info());  // a caching term used for gvcf

  for (unsigned int i_read = 0; i_read < allele_eval.alignments.size(); i_read++) {
    allele_eval.alignments[i_read].sample_index = read_stack[i_read]->sample_index; // this should be done while making the initial pileup
    string sample_name = candidate_variant.variant.sampleNames[allele_eval.alignments[i_read].sample_index];
    sample_name = sample_name.substr(0, sample_name.find("."));
    if (sample_name != "tumor") continue;

    Alignment al = *read_stack[i_read];
    const string basecall = allele_eval.alignments[i_read].basecall;
    // cerr << "basecall: " << basecall << endl;
    if (basecall == "N") continue;

    // double e = allele_eval.alignments[i_read].error_prob;
    int qscore = allele_eval.alignments[i_read].qscore;
    int read_pos = allele_eval.alignments[i_read].pos_in_read;

    const uint8_t call_id(bam_seq_code_to_id(get_bam_seq_code(basecall.c_str()[0])));

    // unsigned align_strand_read_pos(read_pos);
    // unsigned end_trimmed_read_len(read_end);
    // if (al.is_reverse_strand) {
    //   align_strand_read_pos = read_size - (read_pos + 1); // we have already done that with pos_in_read
    //   end_trimmed_read_len = read_size - fwd_strand_begin_skip;
    // }

    const base_call bc = base_call (
      call_id,
      qscore,
      !al.is_reverse_strand, // best_al.is_fwd_strand
      read_pos,              // align_strand_read_pos
      al.alignment.Length,   // al.end_trimmed_read_len,
      false,                 // current_call_filter,
      false,                 // is_neighbor_mismatch,
      true                  // is_tier_specific_filter
    );

    // remember bc_buff is sample-specific, so insert into the right buffer
    // we're doing tumor now (see sample_name above)
    sif.bc_buff.insert_pos_basecall(ref_pos, is_tier1, bc);


    // update mapq and rank-sum metrics
    if (opt.is_compute_VQSRmetrics) {
      int mapq = read_stack[i_read]->alignment.MapQuality;
      // insert_mapq_count(ref_pos, sample_no, mapq);
      sif.bc_buff.insert_mapq_count(ref_pos, mapq);
      // update_ranksum(ref_pos, sample_no, bc, mapq, align_strand_read_pos);
      sif.bc_buff.update_ranksums(ref.get_base(ref_pos), ref_pos, bc, mapq, read_pos);
    }

    if (opt.is_compute_hapscore) {
      string seq = string(read_stack[i_read]->read_bases);
      string qscore = string(read_stack[i_read]->read_qual);
      if (read_stack[i_read]->is_reverse_strand) {
        // Why do I need to reverse them again?
        RevComplementInPlace(seq);
        RevInPlace(qscore);
      }
      uint8_t s[seq.length() / 2 + 1];
      uint8_t qual[seq.length()];
      std::fill(s, s + seq.length() / 2 + 1, 0);
      for (unsigned i = 0; i < seq.length(); i++) {
        s[i / 2] |= (get_bam_seq_code(seq[i]) << 4*(1 - (i % 2)));
        qual[i] = qscore[i];
      }
      bam_seq bseq(s, seq.length(), 0);
      // insert_hap_cand(ref_pos, sample_no, is_tier1, bseq, qual, read_pos);
      sif.bc_buff.insert_hap_cand(ref_pos, is_tier1, bseq, qual, read_pos);
    }
  }

  // Perform the ritual from process_pos_snp_single_sample_impl()
  snp_pos_info null_pi;
  snp_pos_info* pi_ptr(sif.bc_buff.get_pos(ref_pos));
  if (NULL == pi_ptr) pi_ptr = &null_pi;
  snp_pos_info &pi(*pi_ptr);

  const unsigned n_calls(pi.calls.size());
  const unsigned n_spandel(pi.n_spandel); // number of spanning deletions
  const unsigned n_submapped(pi.n_submapped); // reads failing variant calling mapping thresholds

  const pos_t output_pos(ref_pos + 1);

  // pi.set_ref_base(ref_reader.base(chr_idx, ref_pos));
  pi.set_ref_base(ref.get_base(ref_pos));

  // for all but coverage-tests, we use a high-quality subset of the basecalls:
  //
  snp_pos_info &good_pi(sif.epd.good_pi); // stores the column of basecalls actually used for snp-calling after the
                                          // mismatch density filter and other quality filters have been applied

  // Initialize with filtered calls
  good_pi.clear();
  good_pi.set_ref_base(pi.get_ref_base());
  for (unsigned i(0); i < n_calls; ++i) {
    if (pi.calls[i].is_call_filter) continue;
    good_pi.calls.push_back(pi.calls[i]);
  }

  // Store filter counts. This structure is only used for making gvcf, but here
  // it is also used to store n_used_calls and n_unused_calls.
  _site_info.n_used_calls = (good_pi.calls.size());
  _site_info.n_unused_calls = (n_calls - _site_info.n_used_calls);

  sif.ss.update(n_calls);                      // ss:      depth_stream_stat_range
  sif.used_ss.update(_site_info.n_used_calls); // used_ss: depth_stream_stat_range
  if (pi.get_ref_base() != 'N') {
    sif.ssn.update(n_calls);                      // ssn:      depth_stream_stat_range
    sif.used_ssn.update(_site_info.n_used_calls); // used_ssn: depth_stream_stat_range
    // // wav: regional basecall average window
    // sif.wav.insert(pos, _site_info.n_used_calls, _site_info.n_unused_calls, n_spandel, n_submapped);
  }
  else {
    // sif.wav.insert_null(pos);
  }

  // note multi-sample status -- can still be called only for one sample
  // and only for sample 0. working on generalization:
  //
  // if (sample_no != 0) return; // no idea why

  if (pi.calls.empty()) return;

  for (const base_call &bc: good_pi.calls) {
    assert(bc.base_id != BASE_ID::ANY);
  }
  dependent_prob_cache dpcache;
  adjust_joint_eprob(opt, dpcache, good_pi, sif.epd.dependent_eprob, false /* _is_dependent_eprob */);

  const extended_pos_info good_epi(good_pi, sif.epd.dependent_eprob);
  cerr << "used calls: " << _site_info.n_used_calls << ", unused: " << _site_info.n_unused_calls << endl;

  // get fraction of filtered bases:
#if 0
  const double filter_fraction(static_cast<double>(n_unused_calls)/static_cast<double>(n_calls));
  const bool is_overfilter(filter_fraction > opt.max_basecall_filter_fraction);
#endif

  // delay writing any snpcalls so that anomaly tests can (optionally) be applied as filters:
  //
  nonref_test_call nrc;
  //lrt_snp_call lsc;
  _site_info.dgt.reset();
  //monoploid_genotype mgt;
  //std::auto_ptr<nploid_genotype> ngt_ptr;

  if (opt.is_counts) {
    // report_counts(good_pi, _site_info.n_unused_calls, output_pos, *_client_io.counts_osptr());
    report_counts(good_pi, _site_info.n_unused_calls, output_pos, cout);
  }

  if (opt.is_nonref_test() || opt.is_nonref_sites() || true) {
    position_nonref_2allele_test(good_pi, opt, opt.is_nonref_sites(), nrc);
#if 0
    static const bool is_mle_freq(false);

    position_nonref_test(good_pi,
        opt.nonref_variant_rate,
        opt.min_nonref_freq,
        is_mle_freq,
        nrc);
#endif

  }

#if 0
  if (opt.is_lsnp) {
    position_snp_call_lrt(opt.lsnp_alpha, good_pi, lsc);
  }
#endif
  if (true or opt.is_bsnp_diploid()) {
    dopt.pdcaller().position_snp_call_pprob_digt(opt, good_epi, _site_info.dgt, opt.is_all_sites());
    cerr << "dgt: " << _site_info.dgt << endl;
  }
#if 0
  if (opt.is_bsnp_monoploid) {
    position_snp_call_pprob_monogt(opt.bsnp_monoploid_theta,good_pi,mgt);
  }
  if (opt.is_bsnp_nploid) {
    ngt_ptr.reset(new nploid_genotype(*_ninfo));
    position_snp_call_pprob_nploid(opt.bsnp_nploid_snp_prob,good_pi,*_ninfo,*ngt_ptr);
  }
#endif

  //    const bool is_snp(nrc.is_snp || lsc.is_snp || _site_info.dgt.is_snp || mgt.is_snp || (ngt_ptr.get() && ngt_ptr->is_snp));
  const bool is_snp(nrc.is_snp || _site_info.dgt.is_snp);

  // find anomalies:
  //
#if 0
  bool is_pos_adis(false);
  bool is_pos_acov(false);

  if ((opt.is_adis_table || opt.is_adis_lrt) && is_snp) {
    if (opt.is_adis_table) {
      is_pos_adis = (is_pos_adis || position_strand_distro_anomaly(opt.adis_table_alpha,good_pi,_ws));
    }
    if (opt.is_adis_lrt) {
      is_pos_adis = (is_pos_adis || position_strand_distro_anomaly_lrt(opt.adis_lrt_alpha,good_pi));
    }
  }
  if (opt.is_acov) {
    is_pos_acov = position_strand_coverage_anomaly(opt.acov_alpha,pi);
  }
#endif

  //const bool is_anomaly(is_pos_adis || is_pos_acov);
  //const bool is_filter_snp(is_overfilter || (opt.is_filter_anom_calls && is_anomaly));

  //    const bool is_nf_snp(is_snp && (! is_filter_snp));
  if (is_snp) {
    if (opt.is_compute_hapscore) {
      _site_info.hapscore = get_hapscore(pi.hap_set);
      cerr << "hapscore: " << _site_info.hapscore << endl;
    }

    // do calculate VQSR metrics
    if (opt.is_compute_VQSRmetrics) {
      _site_info.MQ             = pi.get_rms_mq();
      _site_info.ReadPosRankSum = pi.get_read_pos_ranksum();
      _site_info.MQRankSum      = pi.get_mq_ranksum();
      _site_info.BaseQRankSum   = pi.get_baseq_ranksum();
    }

    // hpol filter
    _site_info.hpol = get_snp_hpol_size(ref_pos, ref);

    cerr << "HPLEN: " << _site_info.hpol << " MQ: " << _site_info.MQ << " ReadPosRankSum: " << _site_info.ReadPosRankSum << " MQRankSum: " << _site_info.MQRankSum << " BaseQRankSum: " << _site_info.BaseQRankSum << endl;
  }


  if (opt.is_all_sites()) {
#if 0
    const diploid_genotype* dgt_ptr(&_site_info.dgt);
    if (is_filter_snp) {
      dgt_ptr=&get_empty_dgt(pi.ref_base);
    }
#endif

    // Add site gvcf
    if (opt.is_gvcf_output()) {
      _site_info.init(ref_pos, pi.get_ref_base(), good_pi, opt.used_allele_count_min_qscore);
      _gvcfer->add_site(_site_info);
    }


    if (opt.is_bsnp_diploid_allele_file) {
      write_bsnp_diploid_allele(opt, /* _client_io,*/ opt.bam_seq_name, output_pos, pi.get_ref_base(), _site_info.n_used_calls, _site_info.n_unused_calls, good_pi, _site_info.dgt, _site_info.hpol);
    }
  }

  if (opt.is_nonref_sites()) {
    //std::ostream& bos(*_client_io.nonref_sites_osptr());
    std::ostream& bos(cout);
    write_snp_prefix_info_file(opt.bam_seq_name,output_pos, pi.get_ref_base(), _site_info.n_used_calls, _site_info.n_unused_calls, bos);
    bos << "\t";
    write_nonref_2allele_test(opt,good_pi,nrc,bos);
    bos << "\n";
  }

  // report events:
  //
  bool is_reported_event(false);

  std::ostream& report_os(std::cerr);

  if (is_snp) {
    cerr << "is_snp: " << is_snp << endl;
    if (nrc.is_snp) {
      cerr << "nrc.is_snp: " << nrc.is_snp << endl;
      // std::ostream& bos(*_client_io.nonref_test_osptr());
      std::ostream& bos(cout);
      write_snp_prefix_info_file(opt.bam_seq_name, output_pos, pi.get_ref_base(), _site_info.n_used_calls, _site_info.n_unused_calls, bos);
      bos << "\t";
      write_nonref_2allele_test(opt, good_pi, nrc, bos);
#if 0
      write_nonref_test(opt,good_pi,nrc,bos);
#endif
      bos << "\n";
    }
#if 0
    if (lsc.is_snp) {
      write_snp_prefix_info("LSNP",output_pos, pi.ref_base,_site_info.n_used_calls,_site_info.n_unused_calls,report_os);
      report_os << " " << lsc << "\n";
    }
#endif
    if (_site_info.dgt.is_snp) {
      if (opt.is_bsnp_diploid_file) {
        // std::ostream& bos(*_client_io.bsnp_diploid_osptr());
        std::ostream& bos(cout);
        write_snp_prefix_info_file(opt.bam_seq_name,output_pos, pi.get_ref_base(),_site_info.n_used_calls,_site_info.n_unused_calls,bos);
        bos << "\t";
        write_diploid_genotype_snp(opt, good_pi, _site_info.dgt, bos, _site_info.hpol);
        bos << "\n";
      }

      // this needs to be updated no matter where the snp-call is written to:
      if (_is_variant_windows) _variant_print_pos.insert(ref_pos);
    }
#if 0
    if (mgt.is_snp) {
      write_snp_prefix_info("BSNP1",output_pos, pi.ref_base,_site_info.n_used_calls,_site_info.n_unused_calls,report_os);
      report_os << " " << mgt << "\n";
    }
    if (ngt_ptr.get() && ngt_ptr->is_snp) {
      write_snp_prefix_info("BSNPN",output_pos, pi.ref_base,_site_info.n_used_calls,_site_info.n_unused_calls,report_os);
      report_os << " ";
      nploid_write(*_ninfo,*ngt_ptr,report_os);
      report_os << "\n";
    }
#endif

    is_reported_event = true;
    candidate_variant.variant.filter = "PASS";
    candidate_variant.variant.isFiltered = false;
    candidate_variant.variant.isHotSpot = false;
    candidate_variant.variant.quality = 0.0f;
    candidate_variant.variant.info["HRUN"].push_back(convertToString(_site_info.hpol));
  }

#if 0
  if (is_anomaly && (! opt.is_filter_anom_calls)) {
    if (is_pos_adis) report_os << "ANOM_DIS pos: " << output_pos << "\n";
    if (is_pos_acov) report_os << "ANOM_COV pos: " << output_pos << "\n";

    is_reported_event = true;
  }
#endif

  if (opt.is_print_all_site_evidence || (opt.is_print_evidence && is_reported_event)) {
    report_os << "EVIDENCE pos: " << output_pos << "\n"
      << "is_snp: " << is_snp << "\n"
      << pi << "\n";
  }

//  print_delayed_results(STAGE::POST_CALL, pos);
}


void SummarizeInfoFields(Evaluator &eval, vcf::Variant &candidate_variant, int _cur_allele_index, const string &sample_name) {

  float mean_ll_delta;

  eval.ScanSupportingEvidence(mean_ll_delta, _cur_allele_index);

  candidate_variant.info["MLLD"].push_back(convertToString(mean_ll_delta));

  int fwd_strand = 0;
  int rev_strand = 1;

  int var_hyp = 1;
}


// The structure InferenceResult is used for void MultiMinAlleleFreq(...) only.
// qual, gt, gq are the QUAL, GT, GQ computed for min_allele_freq, respectively
struct InferenceResult{
  float min_allele_freq;
  float qual;
  string gt;
  int gq;
  bool operator<(const InferenceResult &rhs) const { return min_allele_freq < rhs.min_allele_freq; } // use the operator "<" for "std::sort"
  bool operator==(const float &rhs) const { return min_allele_freq == rhs; } // use the operator "==" for "std::find"
};

void MultiMinAlleleFreq(Evaluator &eval, VariantCandidate &candidate_variant, int sample_index, ProgramControlSettings &program_flow, int max_detail_level){
  string sample_name = "";
  if(sample_index >= 0) {
    sample_name = candidate_variant.variant.sampleNames[sample_index];
  }

  // info tags needed for multi-min-allele-freq
  vector<string> tags_for_multi_min_allele_freq = {"MUAF", "MUQUAL", "MUGQ", "MUGT",
    "SMAF", "SMQUAL", "SMGQ", "SMGT",
    "MMAF", "MMQUAL", "MMGQ", "MMGT",
    "IMAF", "IMQUAL", "IMGQ", "IMGT",
    "HMAF", "HMQUAL", "HMGQ", "HMGT"};
  // Add the info tags for multi-min-allele-freq if we have not added yet.
  for(unsigned int i_tag = 0; i_tag < tags_for_multi_min_allele_freq.size(); ++i_tag){
    vector<string>::iterator it_format = find(candidate_variant.variant.format.begin(),
        candidate_variant.variant.format.end(),
        tags_for_multi_min_allele_freq[i_tag]);
    if(it_format == candidate_variant.variant.format.end()){
      candidate_variant.variant.format.push_back(tags_for_multi_min_allele_freq[i_tag]);
    }
  }

  // inference_results_union stores the inference results for the union of the min-allele-freq of all available variant types
  vector<InferenceResult> inference_results_union;
  bool is_snp_done = false;
  bool is_mnp_done = false;
  bool is_hotspot_done = false;
  bool is_indel_done = false;

  for(unsigned int alt_allele_index = 0; alt_allele_index < eval.allele_identity_vector.size(); ++alt_allele_index){
    // ptr_maf_vec = the pointer to the multi_min_allele_freq vector of which type of variant for this allele
    vector<float> *ptr_maf_vec = &(program_flow.snp_multi_min_allele_freq);
    string type_prefix = "S";

    // Note that no override here!
    if(eval.allele_identity_vector[alt_allele_index].status.isHotSpot){
      if(is_hotspot_done){
        continue;
      }
      ptr_maf_vec = &(program_flow.hotspot_multi_min_allele_freq);
      type_prefix = "H";
    }
    else if(eval.allele_identity_vector[alt_allele_index].ActAsSNP()){
      if(is_snp_done){
        continue;
      }
      ptr_maf_vec = &(program_flow.snp_multi_min_allele_freq);
      type_prefix = "S";
    }
    else if(eval.allele_identity_vector[alt_allele_index].ActAsMNP()){
      if(is_mnp_done){
        continue;
      }
      ptr_maf_vec = &(program_flow.mnp_multi_min_allele_freq);
      type_prefix = "M";
    }
    else if(eval.allele_identity_vector[alt_allele_index].ActAsHPIndel()){
      if(is_indel_done){
        continue;
      }
      ptr_maf_vec = &(program_flow.indel_multi_min_allele_freq);
      type_prefix = "I";
    }else{
      if(is_snp_done){
        continue;
      }
      ptr_maf_vec = &(program_flow.snp_multi_min_allele_freq);
      string type_prefix = "S";
    }

    for(unsigned int i_freq = 0; i_freq < ptr_maf_vec->size(); ++i_freq){
      float loc_min_allele_freq = ptr_maf_vec->at(i_freq);
      vector<InferenceResult>::iterator it = find(inference_results_union.begin(), inference_results_union.end(), loc_min_allele_freq);
      float loc_qual;
      int loc_gq;
      string loc_gt;

      if(it == inference_results_union.end()){ // This the first time we get loc_min_allele_freq
        int genotype_call;
        float evaluated_genotype_quality;
        // Let's do the inference for the given loc_min_allele_freq
        // eval.allele_eval.CallGermline(loc_min_allele_freq, genotype_call, evaluated_genotype_quality, loc_qual);

        vector<int> genotype_component = {eval.diploid_choice[0], eval.diploid_choice[1]}; // starts with het var

        if(true) { // genotype_call == 2){ //hom var
          genotype_component[0] = eval.diploid_choice[1];
        }
        else { // if(genotype_call == 0){ //hom ref
          genotype_component[1] = eval.diploid_choice[0];
        }

        loc_gt = convertToString(genotype_component[0]) + "/" + convertToString(genotype_component[1]);
        loc_gq = int(round(evaluated_genotype_quality)); // genotype quality is rounded as an integer.
        // append inference_results_union
        inference_results_union.push_back({loc_min_allele_freq, loc_qual, loc_gt, loc_gq});
      }
      else{ // We've seen loc_min_allele_freq before. Don't need to call CallGermline(...) again.
        loc_qual = it->qual;
        loc_gq = it->gq;
        loc_gt = it->gt;
      }

      // write the info tag for the corresponding var type
      candidate_variant.variant.samples[sample_name][type_prefix + "MUAF"].push_back(convertToString(loc_min_allele_freq));
      candidate_variant.variant.samples[sample_name][type_prefix + "MUQUAL"].push_back(convertToString(loc_qual));
      candidate_variant.variant.samples[sample_name][type_prefix + "MUGT"].push_back(loc_gt);
      candidate_variant.variant.samples[sample_name][type_prefix + "MUGQ"].push_back(convertToString(loc_gq));

      switch(type_prefix[0]){
        case 'S':
          is_snp_done = true;
          break;
        case 'M':
          is_mnp_done = true;
          break;
        case 'I':
          is_indel_done = true;
          break;
        case 'H':
          is_hotspot_done = true;
          break;
      }
    }
  }

  // sort inference_results_union according to min_allele_freq in the ascending order
  sort(inference_results_union.begin(), inference_results_union.end());
  // write the info tag for the union of min_allele_freq of the var types
  for(unsigned int i_freq = 0; i_freq < inference_results_union.size(); ++i_freq){
    candidate_variant.variant.samples[sample_name]["MUAF"].push_back(convertToString(inference_results_union[i_freq].min_allele_freq));
    candidate_variant.variant.samples[sample_name]["MUQUAL"].push_back(convertToString(inference_results_union[i_freq].qual));
    candidate_variant.variant.samples[sample_name]["MUGT"].push_back(inference_results_union[i_freq].gt);
    candidate_variant.variant.samples[sample_name]["MUGQ"].push_back(convertToString(inference_results_union[i_freq].gq));
  }
}


void GlueOutputVariant(Evaluator &eval, VariantCandidate &candidate_variant, const ExtendParameters &parameters, int _best_allele_index, int sample_index){
  string sample_name = "";
  if(sample_index >= 0) {
    sample_name = candidate_variant.variant.sampleNames[sample_index];
  }

  DecisionTreeData my_decision(*(eval.variant));
  my_decision.tune_sbias = parameters.my_controls.sbias_tune;
  my_decision.SetupFromMultiAllele(eval);

  float smallest_allele_freq = 1.0f;
  for (unsigned int _alt_allele_index = 0; _alt_allele_index < my_decision.allele_identity_vector.size(); _alt_allele_index++) {
    // for each alt allele, do my best
    // thresholds here can vary by >type< of allele
    float local_min_allele_freq = FreqThresholdByType(eval.allele_identity_vector[_alt_allele_index], parameters.my_controls,
        candidate_variant.variant_specific_params[_alt_allele_index]);

    if (local_min_allele_freq < smallest_allele_freq){
      smallest_allele_freq = local_min_allele_freq;  // choose least-restrictive amongst multialleles
    }

    /* The following piece of code seems redundant. Perhaps due to historical reasons?
       eval.ComputePosteriorGenotype(_alt_allele_index, local_min_allele_freq,
       my_decision.summary_info_vector[_alt_allele_index].genotype_call,
       my_decision.summary_info_vector[_alt_allele_index].gt_quality_score,
       my_decision.summary_info_vector[_alt_allele_index].variant_qual_score);
       */
    SummarizeInfoFields(eval, *(eval.variant), _alt_allele_index, sample_name);
  }

  my_decision.best_allele_index = _best_allele_index;
  my_decision.best_allele_set = true;

  // tell the evaluator to do a genotype
  // choose a diploid genotype
  // return it and set it so that decision tree cannot override

  //@TODO: fix this frequency to be sensible
  float local_min_allele_freq = smallest_allele_freq; // must choose a qual score relative to some frequency

  //eval.MultiAlleleGenotype(local_min_allele_freq,
  //    my_decision.eval_genotype.genotype_component,
  //    my_decision.eval_genotype.evaluated_genotype_quality,
  //    my_decision.eval_genotype.evaluated_variant_quality,
  //    parameters.my_eval_control.max_detail_level);

  my_decision.eval_genotype.genotype_already_set = true; // because we computed it here

  // Tell me what QUAL means (Is QUAL for a ref call or for a var call?) in case we won't show GT in the info tag.
  candidate_variant.variant.samples[sample_name]["QT"].push_back(convertToString(not my_decision.eval_genotype.IsReference()));

  // and I must also set for each allele so that the per-allele filter works
  for(unsigned int i_allele = 0; i_allele < my_decision.allele_identity_vector.size(); ++i_allele){
    my_decision.summary_info_vector[i_allele].variant_qual_score = my_decision.eval_genotype.evaluated_variant_quality;
    my_decision.summary_info_vector[i_allele].gt_quality_score = my_decision.eval_genotype.evaluated_genotype_quality;
  }

  // now that all the data has been gathered describing the variant, combine to produce the output
  my_decision.DecisionTreeOutputToVariant(candidate_variant, parameters, sample_index);
}

// Read and process records appropriate for this variant; positions are zero based
void Evaluator::StackUpOneVariant(const ExtendParameters &parameters, VariantCandidate &candidate_variant, const PositionInProgress& bam_position) {

  // Initialize random number generator for each stack -> ensure reproducibility
  RandSchrange RandGen(parameters.my_controls.RandSeed);

  read_stack.clear();  // reset the stack
  read_stack.reserve(2 * parameters.my_controls.downSampleCoverage);
  read_stack_n.reserve(parameters.my_controls.downSampleCoverage);
  read_stack_t.reserve(parameters.my_controls.downSampleCoverage);
  int read_counter = 0;
  int read_counter_n = 0;
  int read_counter_t = 0;

  for (Alignment* rai = bam_position.begin; rai != bam_position.end; rai = rai->next) {

    // Check global conditions to stop reading in more alignments
    if (rai->original_position > multiallele_window_start)
      break;

    if (rai->alignment.Position > multiallele_window_start)
      continue;

    if (rai->filtered)
      continue;

    if (rai->alignment.GetEndPosition() < multiallele_window_end)
      continue;

    // Get sample name with its numeric suffix stripped
    string sample_name = candidate_variant.variant.sampleNames[rai->sample_index];
    sample_name = sample_name.substr(0, sample_name.find("."));

    // Reservoir Sampling
    if (read_stack.size() < 2 * (unsigned int)parameters.my_controls.downSampleCoverage) {
      read_counter++;
      read_stack.push_back(rai);
    }
    else {
      read_counter++;
      // produces a uniformly distributed test_position between [0, read_counter - 1]
      unsigned int test_position = ((double)RandGen.Rand() / ((double)RandGen.RandMax + 1.0)) * (double)read_counter;
      if (test_position < (unsigned int)parameters.my_controls.downSampleCoverage)
        read_stack[test_position] = rai;
    }

    if (sample_name == "tumor") {
      if (read_stack_t.size() < (unsigned int)parameters.my_controls.downSampleCoverage) {
        read_counter_t++;
        read_stack_t.push_back(rai);
      }
      else {
        read_counter_t++;
        // produces a uniformly distributed test_position between [0, read_counter - 1]
        unsigned int test_position = ((double)RandGen.Rand() / ((double)RandGen.RandMax + 1.0)) * (double)read_counter;
        if (test_position < (unsigned int)parameters.my_controls.downSampleCoverage)
          read_stack_t[test_position] = rai;
      }
    }

    if (sample_name == "normal") {
      if (read_stack_n.size() < (unsigned int)parameters.my_controls.downSampleCoverage) {
        read_counter_n++;
        read_stack_n.push_back(rai);
      }
      else {
        read_counter_n++;
        // produces a uniformly distributed test_position between [0, read_counter - 1]
        unsigned int test_position = ((double)RandGen.Rand() / ((double)RandGen.RandMax + 1.0)) * (double)read_counter;
        if (test_position < (unsigned int)parameters.my_controls.downSampleCoverage)
          read_stack_n[test_position] = rai;
      }
    }
  }
}

bool ProcessOneVariant (
  PersistingThreadObjects &thread_objects,
  VariantCallerContext& vc,
  VariantCandidate &candidate_variant,
  const PositionInProgress& bam_position
) {
  int chr_idx = vc.ref_reader->chr_idx(candidate_variant.variant.sequenceName.c_str());

  Evaluator eval(candidate_variant.variant);
  eval.SetupAllAlleles(*vc.parameters, *vc.global_context, *vc.ref_reader, chr_idx);
  eval.FilterAllAlleles(vc.parameters->my_controls.filter_variant, candidate_variant.variant_specific_params); // put filtering here in case we want to skip below entries

  // We read in one stack per multi-allele variant
  eval.StackUpOneVariant(*vc.parameters, candidate_variant, bam_position);

  cerr << "read stack sizes: " << eval.read_stack_n.size() << ", " << eval.read_stack_t.size() << ", " << eval.read_stack.size() << endl;
  if (eval.read_stack_n.empty() or eval.read_stack_t.empty()) {
    string reason = "NODATA";
    string sample_name = "";
    if (eval.read_stack_n.empty()) {
      cerr << "Nonfatal: No reads found in normal sample for " << candidate_variant.variant.sequenceName << "\t" << eval.multiallele_window_start << endl;
      sample_name = "normal.1";
    }
    if (eval.read_stack_t.empty()) {
      cerr << "Nonfatal: No reads found in tumor sample for " << candidate_variant.variant.sequenceName << "\t" << eval.multiallele_window_start << endl;
      sample_name = "tumor.1";
    }
    NullFilterReason(candidate_variant.variant, sample_name);
    AddFilterReason(candidate_variant.variant, reason, sample_name);
    SetFilteredStatus(candidate_variant.variant, true);
    candidate_variant.variant.samples[sample_name]["GT"].clear();
    candidate_variant.variant.samples[sample_name]["GT"].push_back("./.");
    return false;
  }

  eval.Strelka(thread_objects, *vc.global_context, *vc.parameters, *vc.ref_reader, chr_idx, candidate_variant);
  cerr << candidate_variant.variant << endl;
  return true;

  // handle the unfortunate case in which we must try multiple alleles to be happy
  // try only ref vs alt allele here
  // leave ensemble in ref vs alt state

  // Estimate joint variant likelihood (the product of single-sample likelihoods obtained by piling reads from both samples)
  eval.SampleLikelihood(thread_objects, *vc.global_context, *vc.parameters, *vc.ref_reader, chr_idx, candidate_variant);
  cerr << candidate_variant.variant << endl;
  return true;

  // fill in quantities derived from predictions
  int num_hyp_no_null = eval.allele_identity_vector.size() + 1; // num alleles +1 for ref
  //eval.allele_eval.InitForInference(thread_objects, eval.read_stack, *vc.global_context, num_hyp_no_null, eval.allele_identity_vector);

  // output to variant
  // GlueOutputVariant(eval, candidate_variant, *vc.parameters, best_allele, sample_index);

  // output the inference results (MUQUAL, MUGT, MUGQ, etc.) if I turn on multi_min_allele_freq
  if (true or vc.parameters->program_flow.is_multi_min_allele_freq) {
    // MultiMinAlleleFreq(eval, candidate_variant, sample_index, vc.parameters->program_flow, vc.parameters->my_eval_control.max_detail_level);
  }

  cerr << "diagnostic: " << vc.parameters->program_flow.rich_json_diagnostic << endl;
  // test diagnostic output for this ensemble
  // if (vc.parameters->program_flow.rich_json_diagnostic & (!(eval.variant->isFiltered) | eval.variant->isHotSpot)) // look at everything that came through
  if (vc.parameters->program_flow.rich_json_diagnostic) // look at everything
    JustOneDiagnosis(eval, *vc.global_context, vc.parameters->program_flow.json_plot_dir, true);
  if (vc.parameters->program_flow.minimal_diagnostic & (!(eval.variant->isFiltered) | eval.variant->isHotSpot)) // look at everything that came through
    JustOneDiagnosis(eval, *vc.global_context, vc.parameters->program_flow.json_plot_dir, false);

  return true;
}





