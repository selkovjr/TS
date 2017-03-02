#define starling_streams_base int

#include "HandleVariant.h"


#include "blt_common/adjust_joint_eprob.hh"
#include "blt_common/nonref_test_call.hh"
#include "blt_common/position_nonref_2allele_test.hh"
#include "blt_common/ref_context.hh"
#include "blt_common/snp_pos_info.hh"
#include "blt_util/bam_seq.hh"
#include "blt_util/depth_stream_stat_range.hh"
#include "blt_util/seq_util.hh"
#include "starling_common/gvcf_locus_info.hh"
#include "starling_common/pos_basecall_buffer.hh"
#include "starling_common/starling_ref_seq.hh" // get_starling_ref_seq()
#include "strelka/strelka_shared.hh" // for strelka_options
//#include "starling_common/starling_pos_processor_base.hh"
#include "strelka_pos_processor.h"


#include <iostream>
using namespace std;


static void write_snp_prefix_info_file (
  const string& seq_name,
  const pos_t output_pos,
  const char ref,
  const unsigned n_used_calls,
  const unsigned n_unused_calls,
  ostream& os
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
  const string& seq_name,
  const pos_t output_pos,
  const char ref,
  const unsigned n_used_calls,
  const unsigned n_unused_calls,
  const snp_pos_info& good_pi,
  const diploid_genotype& dgt,
  const unsigned hpol = 0
) {
  // ostream& os(*client_io.bsnp_diploid_allele_osptr());
  ostream& os(cout);

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
  vector<float> dependent_eprob;
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
  // isync_default(indel_buff, eth_buff, sample_opt),
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
  // depth_buffer eth_buff; // provide an early estimate of read depth before realignment.

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


static void report_counts(const snp_pos_info& pi, const unsigned n_unused_calls, const pos_t output_pos, ostream& os) {
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
) {
  bool is_tier1 = true;
  long ref_pos = candidate_variant.variant.position;

  bool _is_variant_windows(opt.variant_windows.size());
  std::set<pos_t> _variant_print_pos;

  // sample_info &sif(sample(sample_no));
  const unsigned knownref_report_size(get_ref_seq_known_size(ref, dopt.report_range));
  sample_info sif(opt, dopt.report_range.size(), knownref_report_size);

  // Data structure defining parameters for a single site to be used for writing in gvcf_aggregator
  site_info _site_info = *(new site_info());  // a caching term used for gvcf

  // Run through the pileup at this position
  for (unsigned int i_read = 0; i_read < pileup.size(); i_read++) {
    pileup[i_read].sample_index = read_stack[i_read]->sample_index; // this should be done while making the initial pileup
    string sample_name = candidate_variant.variant.sampleNames[pileup[i_read].sample_index];
    sample_name = sample_name.substr(0, sample_name.find("."));
    cerr << "sample name: " << sample_name << endl;
    if (sample_name != "tumor") cerr << "skipping normal sample\n" << endl;
    if (sample_name != "tumor") continue;

    Alignment al = *read_stack[i_read];
    const string basecall = pileup[i_read].basecall;
    // cerr << "basecall: " << basecall << endl;
    if (basecall == "N") continue;

    // double e = pileup[i_read].error_prob;
    int qscore = pileup[i_read].qscore;
    int read_pos = pileup[i_read].pos_in_read;

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
      fill(s, s + seq.length() / 2 + 1, 0);
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
  //auto_ptr<nploid_genotype> ngt_ptr;

  if (opt.is_counts) {
    // report_counts(good_pi, _site_info.n_unused_calls, output_pos, *_client_io.counts_osptr());
    report_counts(good_pi, _site_info.n_unused_calls, output_pos, cout);
  }

  if (opt.is_nonref_test() || opt.is_nonref_sites()) {
    position_nonref_2allele_test(good_pi, opt, opt.is_nonref_sites(), nrc);
#if 0
    static const bool is_mle_freq(false);

    position_nonref_test(good_pi, opt.nonref_variant_rate, opt.min_nonref_freq, is_mle_freq, nrc);
#endif
  }

#if 0
  if (opt.is_lsnp) {
    position_snp_call_lrt(opt.lsnp_alpha, good_pi, lsc);
  }
#endif
  if (opt.is_bsnp_diploid()) { // this forces the call on a candidate
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
    //ostream& bos(*_client_io.nonref_sites_osptr());
    ostream& bos(cout);
    write_snp_prefix_info_file(opt.bam_seq_name,output_pos, pi.get_ref_base(), _site_info.n_used_calls, _site_info.n_unused_calls, bos);
    bos << "\t";
    write_nonref_2allele_test(opt,good_pi,nrc,bos);
    bos << "\n";
  }

  // report events:
  //
  bool is_reported_event(false);

  ostream& report_os(cerr);

  if (is_snp) {
    cerr << "is_snp: " << is_snp << endl;
    if (nrc.is_snp) {
      cerr << "nrc.is_snp: " << nrc.is_snp << endl;
      // ostream& bos(*_client_io.nonref_test_osptr());
      ostream& bos(cout);
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
        // ostream& bos(*_client_io.bsnp_diploid_osptr());
        ostream& bos(cout);
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
