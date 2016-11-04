/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_error.h"
#include "../../util/tmap_sam_convert.h"
#include "../../util/tmap_progress.h"
#include "../../util/tmap_sort.h"
#include "../../util/tmap_definitions.h"
#include "../../util/tmap_rand.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_sa.h"
#include "../../sw/tmap_sw.h"
#include "../../sw/tmap_fsw.h"
#include "../../sw/tmap_vsw.h"
#include "samtools/bam.h"
#include "tmap_map_opt.h"
#include "tmap_map_util.h"

#include "../../realign/realign_c_util.h"

#define tmap_map_util_reverse_query(_query, _ql, _i) \
  for(_i=0;_i<(_ql>>1);_i++) { \
    uint8_t _tmp = _query[_i]; \
    _query[_i] = _query[_ql-_i-1]; \
    _query[_ql-_i-1] = _tmp; \
  }

// sort by strand, min-seqid, min-position
#define __tmap_map_sam_sort_coord_lt(a, b) (  ((a).strand < (b).strand) \
    || ( (a).strand == (b).strand && (a).seqid < (b).seqid) \
    || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos < (b).pos ) \
    ? 1 : 0 )

#define __tmap_map_sam_sort_coord_end_lt(a, b) (  ((a).strand < (b).strand) \
    || ( (a).strand == (b).strand && (a).seqid < (b).seqid) \
    || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos + (a).target_len < (b).pos + (b).target_len) \
    ? 1 : 0 )

// sort by strand, min-seqid, min-position, max score
#define __tmap_map_sam_sort_coord_score_lt(a, b) (  ((a).strand < (b).strand) \
    || ( (a).strand == (b).strand && (a).seqid < (b).seqid) \
    || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos < (b).pos ) \
    || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos == (b).pos && (a).score > (b).score) \
    ? 1 : 0 )

// sort by max-score, min-seqid, min-position, min-strand
#define __tmap_map_sam_sort_score_coord_lt(a, b) (  ((a).score > (b).score) \
    || ((a).score == (b).score && (a).strand < (b).strand) \
    || ((a).score == (b).score && (a).strand == (b).strand && (a).seqid < (b).seqid) \
    || ((a).score == (b).score && (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos < (b).pos) \
    ? 1 : 0 )

  TMAP_SORT_INIT(tmap_map_sam_sort_coord, tmap_map_sam_t, __tmap_map_sam_sort_coord_lt)
  TMAP_SORT_INIT(tmap_map_sam_sort_coord_end, tmap_map_sam_t, __tmap_map_sam_sort_coord_end_lt)
  TMAP_SORT_INIT(tmap_map_sam_sort_coord_score, tmap_map_sam_t, __tmap_map_sam_sort_coord_score_lt)
  TMAP_SORT_INIT(tmap_map_sam_sort_score_coord, tmap_map_sam_t, __tmap_map_sam_sort_score_coord_lt)


  static int
tmap_map_get_amplicon(tmap_refseq_t *refseq,
    int32_t seqid,
    uint32_t start,
    uint32_t end,
    uint32_t *ampl_start,
    uint32_t *ampl_end,
    uint32_t strand)
{
        uint32_t *srh;
        int32_t i, j;
        uint32_t key;
        if (!refseq->bed_exist || refseq->bednum[seqid]==0) return 0;
        if (strand == 0) { srh = refseq->bedstart[seqid]; key = start;}
        else {srh = refseq->bedend[seqid]; key = end;}
        // binary search
        i = 0; j = refseq->bednum[seqid]-1;
        while (i < j) {
                uint32_t m = (i+j)/2;
                if (srh[m] == key) { i = j = m; break;}
                if (srh[m] < key) {
                        i = m+1;
                } else {
                        j = m;
                }
        }
        i--;
        if (i < 0) i = 0;
        j++;
        if (j >=refseq->bednum[seqid]) j = refseq->bednum[seqid]-1;
        uint32_t best = abs(srh[i]-key), bi= i;
        for (i++; i <= j; i++) {
                if (abs(srh[i]-key) < best) {
                        best = abs(srh[i]-key);
                        bi = i;
                }
        }
        i = bi;
        // check for overlap
        if (refseq->bedstart[seqid][i] < end && start < refseq->bedend[seqid][i]) {
                *ampl_start = refseq->bedstart[seqid][i];
                *ampl_end = refseq->bedend[seqid][i];
                return 1;
        }
        return 0;
}


static void
tmap_map_util_set_softclip(tmap_map_opt_t *opt, tmap_seq_t *seq, int32_t *softclip_start, int32_t *softclip_end)
{
  (*softclip_start) = (TMAP_MAP_OPT_SOFT_CLIP_LEFT == opt->softclip_type || TMAP_MAP_OPT_SOFT_CLIP_ALL == opt->softclip_type) ? 1 : 0;
  (*softclip_end) = (TMAP_MAP_OPT_SOFT_CLIP_RIGHT == opt->softclip_type || TMAP_MAP_OPT_SOFT_CLIP_ALL == opt->softclip_type) ? 1 : 0;
  // check if the ZB tag is present...
  if(TMAP_SEQ_TYPE_SAM == seq->type || TMAP_SEQ_TYPE_BAM == seq->type) { // SAM/BAM
      tmap_sam_t *sam = seq->data.sam;
      int32_t zb = tmap_sam_get_zb(sam);
  int32_t za = tmap_sam_get_za(sam);
      // keep 3' soft clipping on if zb tag does not exists or there are too few adapter bases.
      if(-1 != zb && opt->max_adapter_bases_for_soft_clipping < zb && -1 != za && sam->seq->l >= za-5) {
          (*softclip_end) = 0;
      }
  }
}

void
tmap_map_sam_init(tmap_map_sam_t *s)
{
  // set all to zero
  memset(s, 0, sizeof(tmap_map_sam_t));
  // nullify
  s->algo_id = TMAP_MAP_ALGO_NONE;
  s->ascore = s->score = INT32_MIN;
}


void
tmap_map_sam_malloc_aux(tmap_map_sam_t *s)
{
  switch(s->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      s->aux.map1_aux = tmap_calloc(1, sizeof(tmap_map_map1_aux_t), "s->aux.map1_aux");
      break;
    case TMAP_MAP_ALGO_MAP2:
      s->aux.map2_aux = tmap_calloc(1, sizeof(tmap_map_map2_aux_t), "s->aux.map2_aux");
      break;
    case TMAP_MAP_ALGO_MAP3:
      s->aux.map3_aux = tmap_calloc(1, sizeof(tmap_map_map3_aux_t), "s->aux.map3_aux");
      break;
    case TMAP_MAP_ALGO_MAP4:
      s->aux.map4_aux = tmap_calloc(1, sizeof(tmap_map_map4_aux_t), "s->aux.map4_aux");
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      s->aux.map_vsw_aux = tmap_calloc(1, sizeof(tmap_map_map_vsw_aux_t), "s->aux.map_vsw_aux");
      break;
    default:
      break;
  }
}

inline void
tmap_map_sam_destroy_aux(tmap_map_sam_t *s)
{
  switch(s->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      free(s->aux.map1_aux);
      s->aux.map1_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP2:
      free(s->aux.map2_aux);
      s->aux.map2_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP3:
      free(s->aux.map3_aux);
      s->aux.map3_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP4:
      free(s->aux.map4_aux);
      s->aux.map4_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      free(s->aux.map_vsw_aux);
      s->aux.map_vsw_aux = NULL;
      break;
    default:
      break;
  }
}

void
tmap_map_sam_destroy(tmap_map_sam_t *s)
{
  tmap_map_sam_destroy_aux(s);
  free(s->cigar);
  s->cigar = NULL;
  s->n_cigar = 0;
}

tmap_map_sams_t *
tmap_map_sams_init(tmap_map_sams_t *prev)
{
  tmap_map_sams_t *sams = tmap_calloc(1, sizeof(tmap_map_sams_t), "sams");
  sams->sams = NULL;
  sams->n = 0;
  if(NULL != prev) sams->max = prev->max;
  return sams;
}

void
tmap_map_sams_realloc(tmap_map_sams_t *s, int32_t n)
{
  int32_t i;
  if(n == s->n) return;
  for(i=n;i<s->n;i++) {
      tmap_map_sam_destroy(&s->sams[i]);
  }
  s->sams = tmap_realloc(s->sams, sizeof(tmap_map_sam_t) * n, "s->sams");
  for(i=s->n;i<n;i++) {
      // nullify
      tmap_map_sam_init(&s->sams[i]);
  }
  s->n = n;
}

void
tmap_map_sams_destroy(tmap_map_sams_t *s)
{
  int32_t i;
  if(NULL == s) return;
  for(i=0;i<s->n;i++) {
      tmap_map_sam_destroy(&s->sams[i]);
  }
  free(s->sams);
  free(s);
}

tmap_map_record_t*
tmap_map_record_init(int32_t num_ends)
{
  tmap_map_record_t *record = NULL;
  int32_t i;

  record = tmap_calloc(1, sizeof(tmap_map_record_t), "record");
  record->sams = tmap_calloc(num_ends, sizeof(tmap_map_sams_t*), "record->sams");
  record->n = num_ends;
  for(i=0;i<num_ends;i++) {
      record->sams[i] = tmap_map_sams_init(NULL);
  }
  return record;
}

tmap_map_record_t*
tmap_map_record_clone(tmap_map_record_t *src)
{
  int32_t i;
  tmap_map_record_t *dest = NULL;

  if(NULL == src) return NULL;

  // init
  dest = tmap_calloc(1, sizeof(tmap_map_record_t), "dest");
  dest->sams = tmap_calloc(src->n, sizeof(tmap_map_sams_t*), "dest->sams");
  dest->n = src->n;
  if(0 == src->n) return dest;

  // copy over data
  for(i=0;i<src->n;i++) {
      dest->sams[i] = tmap_map_sams_clone(src->sams[i]);
  }

  return dest;
}

void
tmap_map_record_merge(tmap_map_record_t *dest, tmap_map_record_t *src)
{
  int32_t i;
  if(NULL == src || 0 == src->n || src->n != dest->n) return;

  for(i=0;i<src->n;i++) {
      tmap_map_sams_merge(dest->sams[i], src->sams[i]);
  }
}

void
tmap_map_record_destroy(tmap_map_record_t *record)
{
  int32_t i;
  if(NULL == record) return;
  for(i=0;i<record->n;i++) {
      tmap_map_sams_destroy(record->sams[i]);
  }
  free(record->sams);
  free(record);
}

tmap_map_bam_t*
tmap_map_bam_init(int32_t n)
{
  tmap_map_bam_t *b = NULL;
  b = tmap_calloc(1, sizeof(tmap_map_bam_t), "b");
  if(0 < n) b->bams = tmap_calloc(n, sizeof(bam1_t*), "b->bams");
  b->n = n;
  return b;
}

void
tmap_map_bam_destroy(tmap_map_bam_t *b)
{
  int32_t i;
  if(NULL == b) return;
  for(i=0;i<b->n;i++) {
      bam_destroy1(b->bams[i]);
  }
  free(b->bams);
  free(b);
}

tmap_map_bams_t*
tmap_map_bams_init(int32_t n)
{
  tmap_map_bams_t *b = NULL;
  b = tmap_calloc(1, sizeof(tmap_map_bams_t), "b");
  if(0 < n) b->bams = tmap_calloc(n, sizeof(tmap_map_bam_t*), "b->bams");
  b->n = n;
  return b;
}

void
tmap_map_bams_destroy(tmap_map_bams_t *b)
{
  int32_t i;
  if(NULL == b) return;
  for(i=0;i<b->n;i++) {
      tmap_map_bam_destroy(b->bams[i]);
  }
  free(b->bams);
  free(b);
}

inline void
tmap_map_sam_copy(tmap_map_sam_t *dest, tmap_map_sam_t *src)
{
  int32_t i;
  // shallow copy
  (*dest) = (*src);
  // aux data
  tmap_map_sam_malloc_aux(dest);
  switch(src->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      (*dest->aux.map1_aux) = (*src->aux.map1_aux);
      break;
    case TMAP_MAP_ALGO_MAP2:
      (*dest->aux.map2_aux) = (*src->aux.map2_aux);
      break;
    case TMAP_MAP_ALGO_MAP3:
      (*dest->aux.map3_aux) = (*src->aux.map3_aux);
      break;
    case TMAP_MAP_ALGO_MAP4:
      (*dest->aux.map4_aux) = (*src->aux.map4_aux);
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      (*dest->aux.map_vsw_aux) = (*src->aux.map_vsw_aux);
      break;
    default:
      break;
  }
  // cigar
  if(0 < src->n_cigar && NULL != src->cigar) {
      dest->cigar = tmap_malloc(sizeof(uint32_t) * (1 + dest->n_cigar), "dest->cigar");
      for(i=0;i<dest->n_cigar;i++) {
          dest->cigar[i] = src->cigar[i];
      }
  }
}
void
tmap_map_sams_merge(tmap_map_sams_t *dest, tmap_map_sams_t *src)
{
  int32_t i, j;
  if(NULL == src || 0 == src->n) return;

  j = dest->n;
  tmap_map_sams_realloc(dest, dest->n + src->n);
  for(i=0;i<src->n;i++,j++) {
      tmap_map_sam_copy(&dest->sams[j], &src->sams[i]);
  }
}

tmap_map_sams_t *
tmap_map_sams_clone(tmap_map_sams_t *src)
{
  int32_t i;
  tmap_map_sams_t *dest = NULL;

  // init
  dest = tmap_map_sams_init(src);
  if(0 == src->n) return dest;

  // realloc
  tmap_map_sams_realloc(dest, src->n);
  // copy over data
  for(i=0;i<src->n;i++) {
      tmap_map_sam_copy(&dest->sams[i], &src->sams[i]);
  }

  return dest;
}

void
tmap_map_sam_copy_and_nullify(tmap_map_sam_t *dest, tmap_map_sam_t *src)
{
  (*dest) = (*src);
  src->n_cigar = 0;
  src->cigar = NULL;
  switch(src->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      src->aux.map1_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP2:
      src->aux.map2_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP3:
      src->aux.map3_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP4:
      src->aux.map4_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      src->aux.map_vsw_aux = NULL;
      break;
    default:
      break;
  }
}

static bam1_t*
tmap_map_sam_print(tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_map_sam_t *sam, int32_t sam_flowspace_tags, int32_t bidirectional, int32_t seq_eq,
                   int32_t nh, int32_t aln_num, int32_t end_num, int32_t mate_unmapped, tmap_map_sam_t *mate)
{
  int64_t mate_strand, mate_seqid, mate_pos, mate_tlen;
  // mate info
  mate_strand = mate_seqid = mate_pos = mate_tlen = 0;
  if(NULL != mate && 0 == mate_unmapped) {
      mate_strand = mate->strand;
      mate_seqid = mate->seqid;
      mate_pos = mate->pos;
  }
  if(NULL == sam) { // unmapped
      return tmap_sam_convert_unmapped(seq, sam_flowspace_tags, bidirectional, refseq,
                              end_num, mate_unmapped, 0,
                              mate_strand, mate_seqid, mate_pos, NULL);
  }
  else {
      // Note: samtools does not like this value
      if(INT32_MIN == sam->score_subo) {
          sam->score_subo++;
      }
      // Set the mate tlen
      if(NULL != mate && 0 == mate_unmapped) {
          // NB: assumes 5'->3' ordering of the fragments
          if(mate->pos < sam->pos) {
              mate_tlen = sam->pos + sam->target_len - mate->pos;
          }
          else {
              mate_tlen = mate->pos + mate->target_len - sam->pos;
          }
          // NB: first fragment is always positive, the rest are always negative
          if(1 == end_num) { // first end
              if(mate_tlen < 0) mate_tlen = -mate_tlen;
          }
          else { // second end
              if(0 < mate_tlen) mate_tlen = -mate_tlen;
          }
      }
      switch(sam->algo_id) {
        case TMAP_MAP_ALGO_MAP1:
          return tmap_sam_convert_mapped(seq, sam_flowspace_tags, bidirectional, seq_eq, refseq,
                                         sam->strand, sam->seqid, sam->pos, sam->target_len,  aln_num,
                                         end_num, mate_unmapped, sam->proper_pair, sam->num_stds,
                                         mate_strand, mate_seqid, mate_pos, mate_tlen,
                                         sam->mapq, sam->cigar, sam->n_cigar,
                                         sam->score, sam->ascore, sam->pscore, nh, sam->algo_id, sam->algo_stage, "");
          break;
        case TMAP_MAP_ALGO_MAP2:
          if(0 < sam->aux.map2_aux->XI) {
              return tmap_sam_convert_mapped(seq, sam_flowspace_tags, bidirectional, seq_eq, refseq,
                                             sam->strand, sam->seqid, sam->pos, sam->target_len, aln_num,
                                             end_num, mate_unmapped, sam->proper_pair, sam->num_stds,
                                             mate_strand, mate_seqid, mate_pos, mate_tlen,
                                             sam->mapq, sam->cigar, sam->n_cigar,
                                             sam->score, sam->ascore, sam->pscore, nh, sam->algo_id, sam->algo_stage,
                                             "\tXS:i:%d\tXT:i:%d\t\tXF:i:%d\tXE:i:%d\tXI:i:%d",
                                             sam->score_subo,
                                             sam->n_seeds,
                                             sam->aux.map2_aux->XF, sam->aux.map2_aux->XE,
                                             sam->aux.map2_aux->XI);
          }
          else {
              return tmap_sam_convert_mapped(seq, sam_flowspace_tags, bidirectional, seq_eq, refseq,
                                             sam->strand, sam->seqid, sam->pos, sam->target_len, aln_num,
                                             end_num, mate_unmapped, sam->proper_pair, sam->num_stds,
                                             mate_strand, mate_seqid, mate_pos, mate_tlen,
                                             sam->mapq, sam->cigar, sam->n_cigar,
                                             sam->score, sam->ascore, sam->pscore, nh, sam->algo_id, sam->algo_stage,
                                             "\tXS:i:%d\tXT:i:%d\tXF:i:%d\tXE:i:%d",
                                             sam->score_subo,
                                             sam->n_seeds,
                                             sam->aux.map2_aux->XF, sam->aux.map2_aux->XE);
          }
          break;
        case TMAP_MAP_ALGO_MAP3:
          return tmap_sam_convert_mapped(seq, sam_flowspace_tags, bidirectional, seq_eq, refseq,
                                         sam->strand, sam->seqid, sam->pos, sam->target_len, aln_num,
                                         end_num, mate_unmapped, sam->proper_pair, sam->num_stds,
                                         mate_strand, mate_seqid, mate_pos, mate_tlen,
                                         sam->mapq, sam->cigar, sam->n_cigar,
                                         sam->score, sam->ascore, sam->pscore, nh, sam->algo_id, sam->algo_stage,
                                         "\tXS:i:%d\tXT:i:%d",
                                         sam->score_subo,
                                         sam->n_seeds);
          break;
        case TMAP_MAP_ALGO_MAP4:
          return tmap_sam_convert_mapped(seq, sam_flowspace_tags, bidirectional, seq_eq, refseq,
                                         sam->strand, sam->seqid, sam->pos, sam->target_len, aln_num,
                                         end_num, mate_unmapped, sam->proper_pair, sam->num_stds,
                                         mate_strand, mate_seqid, mate_pos, mate_tlen,
                                         sam->mapq, sam->cigar, sam->n_cigar,
                                         sam->score, sam->ascore, sam->pscore, nh, sam->algo_id, sam->algo_stage,
                                         "\tXS:i:%d",
                                         sam->score_subo);
          break;
        case TMAP_MAP_ALGO_MAPVSW:
          return tmap_sam_convert_mapped(seq, sam_flowspace_tags, bidirectional, seq_eq, refseq,
                                         sam->strand, sam->seqid, sam->pos, sam->target_len, aln_num,
                                         end_num, mate_unmapped, sam->proper_pair, sam->num_stds,
                                         mate_strand, mate_seqid, mate_pos, mate_tlen,
                                         sam->mapq, sam->cigar, sam->n_cigar,
                                         sam->score, sam->ascore, sam->pscore, nh, sam->algo_id, sam->algo_stage,
                                         "\tXS:i:%d",
                                         sam->score_subo);
          break;
      }
  }
  return NULL;
}

tmap_map_bam_t*
tmap_map_sams_print(tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_map_sams_t *sams, int32_t end_num,
                    tmap_map_sams_t *mates, int32_t sam_flowspace_tags, int32_t bidirectional, int32_t seq_eq, int32_t min_al_len, double min_coverage, double min_identity, int32_t match_score, uint64_t* filtered)
{
  int32_t i;
  tmap_map_sam_t *mate = NULL;
  int32_t mate_unmapped = 0;
  tmap_map_bam_t *bams = NULL;

  if(NULL != mates) {
      if(0 < mates->n) {
          // assumes the mates are sorted by their alignment score
          mate = &mates->sams[0];
      }
      else {
          mate_unmapped = 1;
      }
  }
  // filter alignments; store filtered out ones as unmapped
  int32_t mapped_cnt = 0;
  int8_t* passed = alloca (sams->n);

  for(i=0;i<sams->n;i++)
  {
    unsigned q_len_cigar, r_len_cigar;
    seq_lens_from_bin_cigar (sams->sams [i].cigar, sams->sams [i].n_cigar, &q_len_cigar, &r_len_cigar);
    int8_t mapped = 1;
    if (mapped && r_len_cigar < min_al_len)
        mapped = 0;
    if (mapped && seq->data.sam->seq->l != 0 && ((double) r_len_cigar) / seq->data.sam->seq->l < min_coverage)
        mapped = 0;
    if (mapped && r_len_cigar != 0 && match_score != 0 && ((double) sams->sams [i].score) / (r_len_cigar*match_score) < min_identity )
        mapped = 0;
    passed [i] = mapped;
    if (mapped)
       ++ mapped_cnt;
    else
       ++ *filtered;
  }

  int32_t written_cnt = 0;
  if(mapped_cnt) {
      bams = tmap_map_bam_init(mapped_cnt);
      for(i=0;i<sams->n;i++)
      {
          if (passed [i])
             bams->bams[written_cnt++] = tmap_map_sam_print(seq, refseq, &sams->sams[i], sam_flowspace_tags, bidirectional, seq_eq, sams->max, i, end_num, mate_unmapped, mate);
      }
  }
  else {
      bams = tmap_map_bam_init(1);
      bams->bams[0] = tmap_map_sam_print(seq, refseq, NULL, sam_flowspace_tags, bidirectional, seq_eq, sams->max, 0, end_num, mate_unmapped, mate);
  }

  return bams;
}

void
tmap_map_util_keep_score(tmap_map_sams_t *sams, int32_t algo_id, int32_t score)
{
  int32_t i, j, cur_score;
  for(i=j=0;i<sams->n;i++) {
      if(TMAP_MAP_ALGO_NONE == algo_id
         || sams->sams[i].algo_id == algo_id) {
          cur_score = sams->sams[i].score;
          if(cur_score != score) { // not the best
              tmap_map_sam_destroy(&sams->sams[i]);
          }
          else {
              if(j < i) { // copy if we are not on the same index
                  tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
              }
              j++;
          }
      }
      else {
          if(j < i) { // copy if we are not on the same index
              tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
          }
          j++;
      }
  }
  // reallocate
  tmap_map_sams_realloc(sams, j);
}

void
tmap_map_sams_filter1(tmap_map_sams_t *sams, int32_t aln_output_mode, int32_t algo_id, tmap_rand_t *rand)
{
  int32_t i, j, k;
  int32_t n_best = 0;
  int32_t best_score, cur_score;
  int32_t best_subo_score;

  if(sams->n <= 1) {
      return;
  }

  for(i=j=0;i<sams->n;i++) {
      if(TMAP_MAP_ALGO_NONE == algo_id
         || sams->sams[i].algo_id == algo_id) {
          j++;
      }
  }
  if(j <= 1) {
      return;
  }

  best_score = best_subo_score = INT32_MIN;
  n_best = 0;
  for(i=0;i<sams->n;i++) {
      if(TMAP_MAP_ALGO_NONE == algo_id
         || sams->sams[i].algo_id == algo_id) {
          cur_score = sams->sams[i].score;
          if(best_score < cur_score) {
              if(0 < n_best) {
                  best_subo_score = best_score;
              }
              best_score = cur_score;
              n_best = 1;
          }
          else if(!(cur_score < best_score)) { // equal
              best_subo_score = best_score; // more than one mapping
              n_best++;
          }
          else if(best_subo_score < cur_score) {
              best_subo_score = cur_score;
          }
          // check sub-optimal
          if(TMAP_MAP_ALGO_MAP2 == sams->sams[i].algo_id
             || TMAP_MAP_ALGO_MAP3 == sams->sams[i].algo_id) {
              cur_score = sams->sams[i].score_subo;
              if(best_subo_score < cur_score) {
                  best_subo_score = cur_score;
              }
          }

      }
  }

  // adjust mapping qualities
  if(1 < n_best) {
      for(i=0;i<sams->n;i++) {
          if(TMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              sams->sams[i].mapq = 0;
          }
      }
  }
  else {
      for(i=0;i<sams->n;i++) {
          if(TMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              cur_score = sams->sams[i].score;
              if(cur_score < best_score) { // not the best
                  sams->sams[i].mapq = 0;
              }
          }
      }
  }

  // adjust suboptimal
  if(TMAP_MAP_ALGO_NONE == algo_id) {
      for(i=0;i<sams->n;i++) {
          sams->sams[i].score_subo = best_subo_score;
      }
  }

  if(TMAP_MAP_OPT_ALN_MODE_ALL == aln_output_mode) {
      return;
  }

  // copy to the front
  if(n_best < sams->n) {
      tmap_map_util_keep_score(sams, algo_id, best_score);
  }

  if(TMAP_MAP_OPT_ALN_MODE_UNIQ_BEST == aln_output_mode) {
      if(1 < n_best) { // there can only be one
          if(TMAP_MAP_ALGO_NONE == algo_id) {
              tmap_map_sams_realloc(sams, 0);
          }
          else {
              // get rid of all of them
              for(i=j=0;i<sams->n;i++) {
                  if(sams->sams[i].algo_id == algo_id) {
                      tmap_map_sam_destroy(&sams->sams[i]);
                  }
                  else {
                      if(j < i) { // copy if we are not on the same index
                          tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                      }
                      j++;
                  }
              }
              tmap_map_sams_realloc(sams, j);
          }
      }
  }
  else if(TMAP_MAP_OPT_ALN_MODE_RAND_BEST == aln_output_mode) { // get a random
      int32_t r = (int32_t)(tmap_rand_get(rand) * n_best);

      // keep the rth one
      if(TMAP_MAP_ALGO_NONE == algo_id) {
          if(0 != r) {
              tmap_map_sam_destroy(&sams->sams[0]);
              tmap_map_sam_copy_and_nullify(&sams->sams[0], &sams->sams[r]);
          }
          // reallocate
          tmap_map_sams_realloc(sams, 1);
      }
      else {
          // keep the rth one
          for(i=j=k=0;i<sams->n;i++) {
              if(sams->sams[i].algo_id == algo_id) {
                  if(k == r) { // keep
                      if(j < i) { // copy if we are not on the same index
                          tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                      }
                      j++;
                  }
                  else { // free
                      tmap_map_sam_destroy(&sams->sams[i]);
                  }
                  k++;
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                  }
                  j++;
              }
          }
          tmap_map_sams_realloc(sams, j);
      }
  }
  else if(TMAP_MAP_OPT_ALN_MODE_ALL_BEST == aln_output_mode) {
      // do nothing
  }
  else {
      tmap_bug();
  }
}

void
tmap_map_sams_filter2(tmap_map_sams_t *sams, int32_t score_thr, int32_t mapq_thr)
{
  int32_t i, j;

  // filter based on score and mapping quality
  for(i=j=0;i<sams->n;i++) {
      if(sams->sams[i].score < score_thr || sams->sams[i].mapq < mapq_thr) {
          tmap_map_sam_destroy(&sams->sams[i]);
      }
      else {
          if(j < i) { // copy if we are not on the same index
              tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
          }
          j++;
      }
  }
  tmap_map_sams_realloc(sams, j);
}

void
tmap_map_sams_filter(tmap_map_sams_t *sams, int32_t aln_output_mode, tmap_rand_t *rand)
{
  tmap_map_sams_filter1(sams, aln_output_mode, TMAP_MAP_ALGO_NONE, rand);
}

void
tmap_map_util_remove_duplicates(tmap_map_sams_t *sams, int32_t dup_window, tmap_rand_t *rand)
{
  int32_t i, next_i, j, k, end, best_score_i, best_score_n, best_score_subo;

  if(dup_window < 0 || sams->n <= 0) {
      return;
  }

  // sort
  // NB: since tmap_map_util_sw_gen_score only sets the end position of the
  // alignment, use that
  tmap_sort_introsort(tmap_map_sam_sort_coord_end, sams->n, sams->sams);

  // remove duplicates within a window
  for(i=j=0;i<sams->n;) {

      // get the change
      end = best_score_i = i;
      best_score_n = 0;
      best_score_subo = sams->sams[end].score_subo;
      while(end+1 < sams->n) {
          if(sams->sams[end].seqid == sams->sams[end+1].seqid
             && sams->sams[end].strand == sams->sams[end+1].strand
             && fabs((sams->sams[end].pos + sams->sams[end].target_len) - (sams->sams[end+1].pos + sams->sams[end+1].target_len)) <= dup_window) {
              // track the best scoring
              if(sams->sams[best_score_i].score == sams->sams[end+1].score) {
                  best_score_i = end+1;
                  best_score_n++;
              }
              else if(sams->sams[best_score_i].score < sams->sams[end+1].score) {
                  best_score_i = end+1;
                  best_score_n = 1;
              }
              if(best_score_subo < sams->sams[end+1].score_subo) {
                  best_score_subo = sams->sams[end+1].score_subo;
              }
              end++;
          }
          else {
              break;
          }
      }
      next_i = end+1;

      // randomize the best scoring
      if(1 < best_score_n) {
          k = (int32_t)(best_score_n * tmap_rand_get(rand)); // make this zero-based
          best_score_n = 0; // make this one-based
          end = i;
          while(best_score_n <= k) { // this assumes we know there are at least "best_score
              if(sams->sams[best_score_i].score == sams->sams[end].score) {
                  best_score_i = end;
                  best_score_n++;
              }
              end++;
          }
      }

      // copy over the best
      if(j != best_score_i) {
          // destroy
          tmap_map_sam_destroy(&sams->sams[j]);
          // nullify
          tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[best_score_i]);
      }

      // copy over sub-optimal score
      sams->sams[j].score_subo = best_score_subo;

      // next
      i = next_i;
      j++;
  }

  // resize
  tmap_map_sams_realloc(sams, j);
}

inline int32_t
tmap_map_util_mapq_score(int32_t seq_len, int32_t n_best, int32_t best_score, int32_t n_best_subo, int32_t best_subo_score, tmap_map_opt_t *opt)
{
  int32_t mapq;

  if(0 == n_best_subo) {
      n_best_subo = 1;
      best_subo_score = opt->score_thr;
  }
  if (opt->use_new_QV == 1) {
  int32_t x = 11*opt->score_match;
  if (best_subo_score < x) {best_subo_score = x; n_best_subo= 1;}
  double sf = (double) (best_score - best_subo_score + 1 );
  if (sf < 0) return 1.0;
  sf /= ((double)opt->score_match+opt->pen_mm);
  sf *= 7.3;
  sf -= log(n_best_subo);
  mapq = (int32_t) (sf+0.9999);
  if (mapq < 0) return 0;
  return mapq;
  }
  /*
     fprintf(stderr, "n_best=%d n_best_subo=%d\n",
     n_best, n_best_subo);
     fprintf(stderr, "best_score=%d best_subo_score=%d\n",
     best_score, best_subo_score);
     */
  // Note: this is the old calculationg, based on BWA-long
  //mapq = (int32_t)((n_best / (1.0 * n_best_subo)) * (best_score - best_subo) * (250.0 / best_score + 0.03 / opt->score_match) + .499);
  //
  double sf = 0.4; // initial scaling factor.  Note: 250 * sf is the maximum mapping quality.
  sf *= 250.0 / ((double)opt->score_match * seq_len); // scale based on the best possible alignment score
  sf *= (n_best / (1.0 * n_best_subo)); // scale based on number of sub-optimal mappings
  sf *= (double)(best_score - best_subo_score + 1 ); // scale based on distance to the sub-optimal mapping
  //sf *= (seq_len < 10) ? 1.0 : log10(seq_len); // scale based on longer reads having more information content
  mapq = (int32_t)(sf + 0.99999);
  if(mapq > 250) mapq = 250;
  if(mapq <= 0) mapq = 1;
  return mapq;
}

inline int32_t
tmap_map_util_mapq(tmap_map_sams_t *sams, int32_t seq_len, tmap_map_opt_t *opt, tmap_refseq_t *refseq)
{
  int32_t i;
  int32_t n_best = 0, n_best_subo = 0;
  int32_t best_score, cur_score, best_subo_score, best_subo_score2;
  int32_t mapq;
  int32_t algo_id = TMAP_MAP_ALGO_NONE;
  int32_t best_repetitive = 0;

  // estimate mapping quality TODO: this needs to be refined
  best_score = INT32_MIN;
  best_subo_score = best_subo_score2 = opt->score_thr;
  n_best = n_best_subo = 0;
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      if(best_score < cur_score) {
          // save sub-optimal
          best_subo_score = best_score;
          n_best_subo = n_best;
          // update
          best_score = cur_score;
          n_best = 1;
          algo_id = (algo_id == TMAP_MAP_ALGO_NONE) ? sams->sams[i].algo_id : -1;
          if(sams->sams[i].algo_id == TMAP_MAP_ALGO_MAP2 && 1 == sams->sams[i].aux.map2_aux->flag) {
              best_repetitive = 1;
          }
          else {
              best_repetitive = 0;
          }
      }
      else if(cur_score == best_score) { // qual
          if(sams->sams[i].algo_id == TMAP_MAP_ALGO_MAP2 && 1 == sams->sams[i].aux.map2_aux->flag) {
              best_repetitive = 1;
          }
          n_best++;
      }
      else {
          if(best_subo_score < cur_score) {
              best_subo_score = cur_score;
              n_best_subo = 1;
          }
          else if(best_subo_score == cur_score) {
              n_best_subo++;
          }
      }
      // get the best subo-optimal score
      cur_score = sams->sams[i].score_subo;
      if(INT32_MIN == cur_score) {
          // ignore
      }
      else if(best_subo_score < cur_score) {
          best_subo_score2 = cur_score;
      }
  }
  if(best_subo_score < best_subo_score2) {
      best_subo_score = best_subo_score2;
      if(0 == n_best_subo) n_best_subo = 1;
  }
  if(1 < n_best || best_score < best_subo_score || 0 < best_repetitive) {
      mapq = 0;
  }
  else {
      mapq = tmap_map_util_mapq_score(seq_len, n_best, best_score, n_best_subo, best_subo_score, opt);
  }
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      if(cur_score == best_score) {
          sams->sams[i].mapq = mapq;
      }
      else {
          sams->sams[i].mapq = 0;
      }
  }

  if (mapq == 0 && refseq->bed_exist) {
    int j = -1;
    // int overlap_s = 0, sec = 0;
    int sh_dis = 10000, sec_dis = 10000;
    for (i=0; i <sams->n; i++) {
      tmap_map_sam_t tmp_sam = sams->sams[i];
      uint32_t start = tmp_sam.pos+1;
      uint32_t end = start+tmp_sam.target_len-1;
      uint32_t ampl_start, ampl_end;
      if (tmp_sam.score != best_score) continue;
      if (tmap_map_get_amplicon(refseq, tmp_sam.seqid, start, end, &ampl_start, &ampl_end, tmp_sam.strand)) {
  /*
  if (abs(ampl_start-start) < 15 && abs(ampl_end-end) < 15) {
      fprintf(stderr, "%d %d %d %d %d\n", ampl_start, start, ampl_end, end, tmp_sam.seqid);
      if (j < 0) j = i;
      else {
    fprintf(stderr, "second one\n");
    return 0; // more than one hits match to an amplicon well
      }
  }
  */
  // use overlap to determine which amplicon to pick
  /*
  int ll = start;
  if (ll < ampl_start) ll = ampl_start;
  int rr = end;
  if (rr > ampl_end) rr = ampl_end;
  int ov = rr -ll;
  if (ov > overlap_s) {
      sec = overlap_s;
      overlap_s = ov;
      j = i;
  } else if (ov > sec) {
      sec = ov;
  }
  */
  // Use the distance to 5' end start
  //fprintf(stderr, "found %d %d\n", ampl_start, ampl_end);
  int dis;
  if (tmp_sam.strand == 0) dis = abs(ampl_start-start);
  else dis = abs(ampl_end-end);
  if (dis < sh_dis) {
      sec_dis = sh_dis;
      sh_dis = dis;
      j = i;
  } else if (dis < sec_dis) {
      sec_dis = dis;
  }
      }
    } // for
    if (j >= 0 && /*overlap_s - sec > 5*/ sec_dis - sh_dis > 1) {
  //fprintf(stderr, "Bump to 12 %d\n", j);
  sams->sams[j].mapq = 12;
  // make the j unique best hit
  for (i = 0; i < sams->n; i++) {
      if (i != j && sams->sams[i].score == best_score)  sams->sams[i].score = best_score-2;
  }
    }
  }
  return 0;
}

// ACGTNBDHKMRSVWYN
// This gives a mask for match iupac characters versus A/C/G/T/N
//  A  C  G  T   N  B  D  H   K  M  R  S   V  W  Y  N
static int32_t matrix_iupac_mask[89] = {
    1, 0, 0, 0,  1, 0, 1, 1,  0, 1, 1, 0,  1, 1, 0, 1, // A
    0, 1, 0, 0,  1, 1, 0, 1,  0, 1, 0, 1,  1, 0, 1, 1, // C
    0, 0, 1, 0,  1, 1, 1, 0,  1, 0, 1, 1,  1, 0, 0, 1, // G
    0, 0, 0, 1,  1, 1, 1, 1,  1, 0, 0, 0,  0, 1, 1, 1, // T
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1  // N
};

#define __map_util_gen_ap_iupac(par, opt) do { \
    int32_t i; \
    for(i=0;i<80;i++) { \
        if(0 < matrix_iupac_mask[i]) (par).matrix[i] = (opt)->score_match; \
        else (par).matrix[i] = -(opt)->pen_mm; \
    } \
    (par).gap_open = (opt)->pen_gapo; (par).gap_ext = (opt)->pen_gape; \
    (par).gap_end = (opt)->pen_gape; \
    (par).row = 16; \
    (par).band_width = (opt)->bw; \
} while(0)

static inline void
tmap_map_util_sw_get_start_and_end_pos(tmap_map_sam_t *sam, int32_t seq_len, uint8_t strand,
                                       uint32_t *start_pos, uint32_t *end_pos)
{
  if(strand == 0) {
      (*start_pos) = sam->pos + 1;
      (*end_pos) = sam->pos + sam->target_len;
  } else {
      (*start_pos) = sam->pos + 1;
      (*end_pos) = sam->pos + seq_len;
  }
}

// NB: this function unrolls banding in some cases
static int32_t
tmap_map_util_sw_gen_score_helper(tmap_refseq_t *refseq, tmap_map_sams_t *sams,
                        tmap_seq_t *seq, tmap_map_sams_t *sams_tmp,
                        int32_t *idx, int32_t start, int32_t end,
                        uint8_t strand, tmap_vsw_t *vsw,
                        int32_t seq_len, uint32_t start_pos, uint32_t end_pos,
                        int32_t *target_mem, uint8_t **target,
                        int32_t softclip_start, int32_t softclip_end,
                        int32_t prev_n_best,
                        int32_t max_seed_band, // NB: this may be modified as banding is unrolled
                        int32_t prev_score, // NB: must be greater than or equal to the scoring threshold
                        tmap_vsw_opt_t *vsw_opt,
                        tmap_rand_t *rand,
                        tmap_map_opt_t *opt)
{
  tmap_map_sam_t tmp_sam;
  uint8_t *query;
  uint32_t qlen;
  int32_t tlen, overflow = 0, is_long_hit = 0;

  // choose a random one within the window
  if(start == end) {
      tmp_sam = sams->sams[start];
  }
  else {
      int32_t r = (int32_t)(tmap_rand_get(rand) * (end - start + 1));
      r += start;
      tmp_sam = sams->sams[r];
  }

  if(0 < opt->long_hit_mult && seq_len * opt->long_hit_mult <= end_pos - start_pos + 1) {
      is_long_hit = 1;
  }

  // update the query sequence
  query = (uint8_t*)tmap_seq_get_bases(seq)->s;
  qlen = tmap_seq_get_bases_length(seq);

  // add in band width
  // one-based
  if(start_pos < opt->bw) {
      start_pos = 1;
  }
  else {
      start_pos -= opt->bw - 1;
  }
  end_pos += opt->bw - 1;
  if(refseq->annos[sams->sams[end].seqid].len < end_pos) {
      end_pos = refseq->annos[sams->sams[end].seqid].len; // one-based
  }

  // get the target sequence
  tlen = end_pos - start_pos + 1;
  if((*target_mem) < tlen) { // more memory?
      (*target_mem) = tlen;
      tmap_roundup32((*target_mem));
      (*target) = tmap_realloc((*target), sizeof(uint8_t)*(*target_mem), "target");
  }
  // NB: IUPAC codes are turned into mismatches
  if(NULL == tmap_refseq_subseq2(refseq, sams->sams[end].seqid+1, start_pos, end_pos, (*target), 1, NULL)) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }

  // reverse compliment the target
  if(1 == strand) {
      tmap_reverse_compliment_int((*target), tlen);
  }

  // Debugging
#ifdef TMAP_VSW_DEBUG
  int j;
  fprintf(stderr, "seqid:%u start_pos=%u end_pos=%u strand=%d\n", sams->sams[end].seqid+1, start_pos, end_pos, strand);
  for(j=0;j<qlen;j++) {
      fputc("ACGTN"[query[j]], stderr);
  }
  fputc('\n', stderr);
  for(j=0;j<tlen;j++) {
      fputc("ACGTN"[(*target)[j]], stderr);
  }
  fputc('\n', stderr);
#endif

  // initialize the bounds
  tmp_sam.result.query_start = tmp_sam.result.query_end = 0;
  tmp_sam.result.target_start = tmp_sam.result.target_end = 0;

  /**
   * Discussion on choosing alignments.
   *
   * In general, we would like to choose the alignment using the *most* number
   * of query bases.  The complication comes when using soft-clipping.
   *
   * When no soft-clipping occurs, then the whole query will be used.
   *
   * When there is only soft-clipping at the 3' end of the query (end of the
   * read), we can search for the alignment with the greatest query end when
   * searching in the 5'->3' direction.  Since the query start is fixed, as there
   * is no 5' start-clipping, we are then guaranteed to find the alignment with
   * the most number of query bases if we fix the query end when we align in the
   * 3'->5' direction.
   *
   * When there is only soft-clipping at the 5' end of the query (start of the
   * read), we can search for the alignment with the greatest query end when
   * searching in the 5'->3' direction, as the query end is fixed.  Then, when
   * aligning in the 3'->5' direction, we can search for the alignment with the
   * smallest query start (which in the 3'->5' direction is the one with the one
   * with the largest query end).
   *
   * We cannot guarantee anything if we are allowing soft-clipping on both the
   * 5' and 3' ends.
   */

  // NB: this aligns in the sequencing direction
  tmp_sam.score = tmap_vsw_process_fwd(vsw, query, qlen, (*target), tlen,
                                   &tmp_sam.result, &overflow, opt->score_thr, 1);

  if(1 < tmp_sam.result.n_best) {
      // What happens if soft-clipping or not soft-clipping causes two
      // alignments to be equally likely? So while we could set the scores to be
      // similar, we can also down-weight the score a little bit.
      tmp_sam.score_subo = tmp_sam.score - opt->pen_gapo - opt->pen_gape;
      //tmp_sam.score_subo = tmp_sam.score - 1;
  }
  if(1 == is_long_hit) {
      // if we have many seeds that span a large tandem repeat, then we should
      // downweight the sub-optimal score
      tmp_sam.score_subo = tmp_sam.score - opt->pen_gapo - opt->pen_gape;
  }

  /*
  fprintf(stderr, "start_pos=%d end_pos=%d start=%d end=%d score=%d n_best=%d\n",
          start_pos, end_pos, start, end, tmp_sam.score, tmp_sam.result.n_best);
          */
  if(1 == overflow) {
      return INT32_MIN;
      //tmap_error("bug encountered", Exit, OutOfRange);
  }

  if(opt->score_thr <= tmp_sam.score) { // NB: we could also specify 'prev_score <= tmp_sam.score'
      int32_t add_current = 1; // yes, by default

      // Save local variables
      // TODO: do we need to save this data in this fashion?
      int32_t t_end = end;
      int32_t t_start = start;
      int32_t t_end_pos = end_pos;
      uint32_t t_start_pos = start_pos;

      if(1 == opt->unroll_banding // unrolling is turned on
         && 0 <= max_seed_band  // do we have enough bases to unrooll
         && 1 < end - start + 1 // are there enough seeds in this group
         && ((prev_n_best < 0 && 1 < tmp_sam.result.n_best) // we need to do a first timeunroll
             || (prev_n_best == tmp_sam.result.n_best))) { // unroll was not successful, try again
          uint32_t start_pos_prev=0;
          uint32_t start_pos_cur=0, end_pos_cur=0;
          uint32_t unrolled = 0;

          while(0 == unrolled) {

              // Reset local variables
              start = t_start;
              end = t_end;

              //fprintf(stderr, "max_seed_band=%d start=%d end=%d\n", max_seed_band, start, end);
              // NB: band based on EXACT start position
              int32_t n = end + 1;
              end = start;
              while(start < n) {
                  int32_t cur_score;
                  // reset start and end position
                  if(start == end) {
                      tmap_map_util_sw_get_start_and_end_pos(&sams->sams[start], seq_len, strand, &start_pos, &end_pos);
                      start_pos_prev = start_pos;
                  }
                  if(end + 1 < n) {
                      if(sams->sams[end].strand == sams->sams[end+1].strand
                         && sams->sams[end].seqid == sams->sams[end+1].seqid) {
                          tmap_map_util_sw_get_start_and_end_pos(&sams->sams[end+1], seq_len, strand, &start_pos_cur, &end_pos_cur);


                             // fprintf(stderr, "start_pos_cur=%u start_pos_prev=%u max_seed_band=%d\n",
                             // start_pos_cur, start_pos_prev, max_seed_band);


                          // consider start positions only
                          if(start_pos_cur - start_pos_prev <= max_seed_band) { // consider start positions only
                              end++;
                              if(end_pos < end_pos_cur) {
                                  end_pos = end_pos_cur;
                              }
                              start_pos_prev = start_pos_cur;
                              continue; // there may be more to add
                          } } } // do not recurse if we did not unroll if(0 < max_seed_band && t_start == start && t_end == end) { // we did not unroll any max_seed_band = (max_seed_band >> 1); break; }
                  unrolled = 1; // we are going to unroll

                  // TODO: we could hash the alignment result, as unrolling may
                  // cause many more SWs

                  // recurse
                  cur_score = tmap_map_util_sw_gen_score_helper(refseq, sams, seq, sams_tmp, idx, start, end,
                                                                strand, vsw, seq_len, start_pos, end_pos,
                                                                target_mem, target,
                                                                softclip_start, softclip_end,
                                                                tmp_sam.result.n_best,
                                                                (max_seed_band <= 0) ? -1 : (max_seed_band >> 1),
                                                                tmp_sam.score,
                                                                vsw_opt, rand, opt);

                  if(cur_score == tmp_sam.score) add_current = 0; // do not add the current alignment, we found it during unrolling
                  // update start/end
                  end++;
                  start = end;
              }
          }
      }

      // Reset local variables
      // TODO: do we need to reset this data in this fashion?
      end = t_end;
      start = t_start;
      end_pos = t_end_pos;
      start_pos = t_start_pos;

      if(1 == add_current) {
          tmap_map_sam_t *s = NULL;

          if(sams_tmp->n <= (*idx)) {
              tmap_map_sams_realloc(sams_tmp, (*idx)+1);
              //tmap_error("bug encountered", Exit, OutOfRange);
          }

          s = &sams_tmp->sams[(*idx)];
          // shallow copy previous data
          (*s) = tmp_sam;
          //s->result = tmp_sam.result;

          // nullify the cigar
          s->n_cigar = 0;
          s->cigar = NULL;

          // adjust target length and position NB: query length is implicitly
          // stored in s->query_end (consider on the next pass)
          s->pos = start_pos - 1; // zero-based
          s->target_len = s->result.target_end + 1;

          if(1 == strand) {
              if(s->pos + tlen <= s->result.target_end) s->pos = 0;
              else s->pos += tlen - s->result.target_end - 1;
          }

          /*
          fprintf(stderr, "strand=%d pos=%d n_best=%d %d-%d %d-%d %d\n",
                  strand,
                  s->pos,
                  s->result.n_best,
                  s->query_start,
                  s->query_end,
                  s->target_start,
                  s->result.target_end,
                  s->target_len);
                  */

          // # of seeds
          s->n_seeds = (end - start + 1);
          // update aux data
          tmap_map_sam_malloc_aux(s);
          switch(s->algo_id) {
            case TMAP_MAP_ALGO_MAP1:
              (*s->aux.map1_aux) = (*tmp_sam.aux.map1_aux);
              break;
            case TMAP_MAP_ALGO_MAP2:
              (*s->aux.map2_aux) = (*tmp_sam.aux.map2_aux);
              break;
            case TMAP_MAP_ALGO_MAP3:
              (*s->aux.map3_aux) = (*tmp_sam.aux.map3_aux);
              break;
            case TMAP_MAP_ALGO_MAP4:
              (*s->aux.map4_aux) = (*tmp_sam.aux.map4_aux);
              break;
            case TMAP_MAP_ALGO_MAPVSW:
              (*s->aux.map_vsw_aux) = (*tmp_sam.aux.map_vsw_aux);
              break;
            default:
              tmap_error("bug encountered", Exit, OutOfRange);
              break;
          }
          (*idx)++;
      }
  }
  return tmp_sam.score;
}

typedef struct {
    uint32_t seqid;
    int8_t strand;
    int32_t start;
    int32_t end;
    uint32_t start_pos;
    uint32_t end_pos;
    int8_t filtered:4;
    int8_t repr_hit:4;
} tmap_map_util_gen_score_t;

tmap_map_sams_t *
tmap_map_util_sw_gen_score(tmap_refseq_t *refseq,
                           tmap_seq_t *seq,
                           tmap_map_sams_t *sams,
                           tmap_seq_t **seqs,
                           tmap_rand_t *rand,
                           tmap_map_opt_t *opt,
                           int32_t *num_after_grouping)
{
  int32_t i, j;
  int32_t start, end;
  tmap_map_sams_t *sams_tmp = NULL;
  int32_t seq_len=0, target_mem=0;
  uint8_t *target=NULL;
  int32_t best_subo_score;
  tmap_vsw_t *vsw = NULL;
  tmap_vsw_opt_t *vsw_opt = NULL;
  uint32_t start_pos, end_pos;
  int32_t softclip_start, softclip_end;
  uint32_t end_pos_prev=0;
  uint32_t start_pos_cur=0, end_pos_cur=0;
  tmap_map_util_gen_score_t *groups = NULL;
  int32_t num_groups = 0, num_groups_filtered = 0;
  double stage_seed_freqc = opt->stage_seed_freqc;
  int32_t max_group_size = 0, repr_hit, filter_ok = 0;

  if(NULL != num_after_grouping) (*num_after_grouping) = 0;

  if(0 == sams->n) {
      return sams;
  }
  // the final mappings will go here
  sams_tmp = tmap_map_sams_init(sams);
  tmap_map_sams_realloc(sams_tmp, sams->n);

  // sort by strand/chr/pos/score
  tmap_sort_introsort(tmap_map_sam_sort_coord, sams->n, sams->sams);
  tmap_map_util_set_softclip(opt, seq, &softclip_start, &softclip_end);

  // initialize opt
  vsw_opt = tmap_vsw_opt_init(opt->score_match, opt->pen_mm, opt->pen_gapo, opt->pen_gape, opt->score_thr);

  // init seqs
  seq_len = tmap_seq_get_bases_length(seqs[0]);

  // forward
  vsw = tmap_vsw_init((uint8_t*)tmap_seq_get_bases(seqs[0])->s, seq_len, softclip_start, softclip_end, opt->vsw_type, vsw_opt);

  // pre-allocate groups
  groups = tmap_calloc(sams->n, sizeof(tmap_map_util_gen_score_t), "groups");

  // determine groups
  num_groups = num_groups_filtered = 0;
  start = end = 0;
  start_pos = end_pos = 0;
  best_subo_score = INT32_MIN; // track the best sub-optimal hit
  repr_hit = 0;
  while(end < sams->n) {
      uint32_t seqid;
      uint8_t strand;

      // get the strand/start/end positions
      seqid = sams->sams[end].seqid;
      strand = sams->sams[end].strand;
      //first pass, setup start and end
      if(start == end) {
          tmap_map_util_sw_get_start_and_end_pos(&sams->sams[start], seq_len, strand, &start_pos, &end_pos);
          end_pos_prev = end_pos;
          repr_hit = sams->sams[end].repr_hit;
      }

      // sub-optimal score
      if(best_subo_score < sams->sams[end].score_subo) {
          best_subo_score = sams->sams[end].score_subo;
      }

      // check if the hits can be banded
      //fprintf(stderr, "end=%d seqid: %d start: %d end: %d\n", end, sams->sams[end].seqid, sams->sams[end].pos, (sams->sams[end].pos + sams->sams[end].target_len));
      if(end + 1 < sams->n) {
          //fprintf(stderr, "%d seed start: %d end: %d next start: %d next end: %d\n", end, sams->sams[end].pos, (sams->sams[end].pos + sams->sams[end].target_len), sams->sams[end+1].pos, (sams->sams[end+1].pos + sams->sams[end+1].target_len));
          if(sams->sams[end].strand == sams->sams[end+1].strand
             && sams->sams[end].seqid == sams->sams[end+1].seqid) {
              tmap_map_util_sw_get_start_and_end_pos(&sams->sams[end+1], seq_len, strand, &start_pos_cur, &end_pos_cur);

              if(start_pos_cur <= end_pos_prev // NB: beware of unsigned int underflow
                 || (start_pos_cur - end_pos_prev) <= opt->max_seed_band) {
                  end++;
                  if(end_pos < end_pos_cur) {
                      end_pos = end_pos_cur;
                  }
                  end_pos_prev = end_pos_cur;
                  repr_hit |= sams->sams[end].repr_hit; // representitive hit
                  continue; // there may be more to add

              }
          }
      }

      // add a group
      num_groups++;
      // update
      groups[num_groups-1].seqid = seqid;
      groups[num_groups-1].strand = strand;
      groups[num_groups-1].start = start;
      groups[num_groups-1].end = end;
      groups[num_groups-1].start_pos = start_pos;
      groups[num_groups-1].end_pos = end_pos;
      groups[num_groups-1].filtered = 0; // assume not filtered, not guilty
      groups[num_groups-1].repr_hit = repr_hit;
      if(max_group_size < end - start + 1) {
          max_group_size = end - start + 1;
      }
      // update start/end
      end++;
      start = end;
  }

  // resize
  if(num_groups < sams->n) {
      groups = tmap_realloc(groups, num_groups * sizeof(tmap_map_util_gen_score_t), "groups");
  }

  // filter groups
  // NB: if we happen to filter all the groups, iteratively relax the filter
  // settings.
  do {
      num_groups_filtered = 0;

      for(i=0;i<num_groups;i++) {
          tmap_map_util_gen_score_t *group = &groups[i];

          // NB: if match/mismatch penalties are on the opposite strands, we may
          // have wrong score
          // NOTE:  end >(sams->n * opt->seed_freqc ) comes from
          /*
           * Anatomy of a hash-based long read sequence mapping algorithm for next
           * generation DNA sequencing
           * Sanchit Misra, Bioinformatics, 2011
           * "For each read, we find the maximum of the number of q-hits in
           * all regions, say C. We keep the cutoff as a fraction f of C. Hence,
           * if a region has ≥fC q-hits, only then it is processed further. "
           */
          /*
          fprintf(stderr, "repr_hit=%d size=%d freq=%lf\n",
                  group->repr_hit,
                  (group->end - group->start + 1),
                  ( sams->n * stage_seed_freqc));
          */

          if(0 == group->repr_hit && (group->end - group->start + 1) > ( sams->n * stage_seed_freqc) ) {
              group->filtered = 0;
          }
          else {
              group->filtered = 1;
              num_groups_filtered++;
          }
      }

      /*
      fprintf(stderr, "stage_seed_freqc=%lf num_groups=%d num_groups_filtered=%d\n",
              stage_seed_freqc,
              num_groups,
              num_groups_filtered);
              */

      if(opt->stage_seed_freqc_min_groups <= num_groups
         && (num_groups - num_groups_filtered) < opt->stage_seed_freqc_min_groups) {
          stage_seed_freqc /= 2.0;
          // reset
          for(i=0;i<num_groups;i++) {
              tmap_map_util_gen_score_t *group = &groups[i];
              group->filtered = 0;
          }
          num_groups_filtered = 0;
          if(stage_seed_freqc < 1.0 / sams->n) { // the filter will accept all hits
              break;
          }
      }
      else { // the filter was succesfull
          filter_ok = 1;
          break;
      }
  } while(1);

  if(num_groups < opt->stage_seed_freqc_min_groups) { // if too few groups, always keep
      // reset
      for(i=0;i<num_groups;i++) {
          tmap_map_util_gen_score_t *group = &groups[i];
          group->filtered = 0;
      }
      num_groups_filtered = 0;
  }
  else if(0 == filter_ok) { // we could not find an effective filter
      int32_t cur_max;
      // Try filtering based on retaining the groups with most number of seeds.
      cur_max = (max_group_size * 2) - 1; // NB: this will be reduced
      // Assume all are filtered
      num_groups_filtered = num_groups;
      while(num_groups - num_groups_filtered < opt->stage_seed_freqc_min_groups) { // too many groups were filtered
          num_groups_filtered = 0;
          // halve the minimum group size
          cur_max /= 2;
          // go through the groups
          for(i=0;i<num_groups;i++) {
              tmap_map_util_gen_score_t *group = &groups[i];
              if(group->end - group->start + 1 < cur_max) { // filter if there were too few seeds
                  group->filtered = 1;
                  num_groups_filtered++;
              }
              else {
                  group->filtered = 0;
              }
          }
      }
  }

  /*
  for(i=j=0;i<num_groups;i++) { // go through each group
      tmap_map_util_gen_score_t *group = &groups[i];
      if(1 == group->filtered && 0 == group->repr_hit) continue;
      j++;
  }
  fprintf(stderr, "num_groups=%d num_groups_filtered=%d j=%d\n", num_groups, num_groups_filtered, j);
  */

  /*
  fprintf(stderr, "final sams->n=%d num_groups=%d num_groups_filtered=%d\n",
          sams->n,
          num_groups,
          num_groups_filtered);
          */

  // process unfiltered...
  for(i=j=0;i<num_groups;i++) { // go through each group
      tmap_map_util_gen_score_t *group = &groups[i];
      /*
      fprintf(stderr, "start=%d end=%d num=%d seqid=%u strand=%d start_pos=%u end_pos=%u repr_hit=%u filtered=%d\n",
              group->start,
              group->end,
              group->end - group->start + 1,
              group->seqid,
              group->strand,
              group->start_pos,
              group->end_pos,
              group->repr_hit,
              group->filtered);
              */
      if(1 == group->filtered && 0 == group->repr_hit) continue;
      /*
      fprintf(stderr, "start=%d end=%d num=%d seqid=%u strand=%d start_pos=%u end_pos=%u filtered=%d\n",
              group->start,
              group->end,
              group->end - group->start + 1,
              group->seqid,
              group->strand,
              group->start_pos,
              group->end_pos,
              group->filtered);
              */

      // generate the score
      tmap_map_util_sw_gen_score_helper(refseq, sams, seqs[0], sams_tmp, &j, group->start, group->end,
                                        group->strand, vsw, seq_len, group->start_pos, group->end_pos,
                                        &target_mem, &target,
                                        softclip_start, softclip_end,
                                        -1, // this is our first call
                                        opt->max_seed_band, // NB: this may be modified as banding is unrolled
                                        opt->score_thr-1,
                                        vsw_opt, rand, opt);
      // save the number of groups
      if(NULL != num_after_grouping) (*num_after_grouping)++;
  }
  //fprintf(stderr, "unfiltered # = %d\n", j);

  // only if we applied the freqc filter
  if(0 < stage_seed_freqc && 0 < opt->stage_seed_freqc_rand_repr) {
      // check if we filtered too many groups, and so if we should keep
      // representative hits.
      /*
      fprintf(stderr, "Should we get representatives %lf < %lf\n",
              opt->stage_seed_freqc_group_frac,
              (double)num_groups_filtered / num_groups);
              */
      if(opt->stage_seed_freqc_group_frac <= (double)num_groups_filtered / num_groups) {
          int32_t best_score = INT32_MIN, best_score_filt = INT32_MIN;
          double pr = 0.0;
          int32_t k, l, c, n;

          // get the best score from a group that passed the filter
          for(k=0;k<j;k++) { // NB: keep k for later
              if(best_score < sams_tmp->sams[k].score) {
                  best_score = sams_tmp->sams[k].score;
              }
          }

          // the # added
          n = 0;

          // pick the one with the most # of seeds
          c = l = -1;
          for(i=0;i<num_groups;i++) {
              tmap_map_util_gen_score_t *group = &groups[i];
              if(0 == group->filtered) continue;
              if(c < group->end - group->start + 1) {
                  c = group->end - group->start + 1;
                  l = i;
              }
          }
          if(0 <= l) {
              tmap_map_util_gen_score_t *group = &groups[l];
              // generate the score
              tmap_map_util_sw_gen_score_helper(refseq, sams, seqs[0], sams_tmp, &j, group->start, group->end,
                                                group->strand, vsw, seq_len, group->start_pos, group->end_pos,
                                                &target_mem, &target,
                                                softclip_start, softclip_end,
                                                -1, // this is our first call
                                                opt->max_seed_band, // NB: this may be modified as banding is unrolled
                                                opt->score_thr-1,
                                                vsw_opt, rand, opt);
              group->filtered = 0; // no longer filtered
          }

          // now, choose uniformly
          pr = num_groups / (1+opt->stage_seed_freqc_rand_repr);
          for(i=c=n=0;i<num_groups;i++,c++) {
              tmap_map_util_gen_score_t *group = &groups[i];
              if(0 == group->filtered) continue;
              if(pr < c) {
                  if(n == opt->stage_seed_freqc_rand_repr) break;
                  // generate the score
                  tmap_map_util_sw_gen_score_helper(refseq, sams, seqs[0], sams_tmp, &j, group->start, group->end,
                                                    group->strand, vsw, seq_len, group->start_pos, group->end_pos,
                                                    &target_mem, &target,
                                                    softclip_start, softclip_end,
                                                    -1, // this is our first call
                                                    opt->max_seed_band, // NB: this may be modified as banding is unrolled
                                                    opt->score_thr-1,
                                                    vsw_opt, rand, opt);
                  group->filtered = 0; // no longer filtered
                  n++;
                  // reset
                  c = 0;
              }
          }

          // get the best from the filtered
          for(;k<j;k++) { // NB: k should not have been modified
              if(best_score_filt < sams_tmp->sams[k].score) {
                  best_score_filt = sams_tmp->sams[k].score;
              }
          }

          /*
          fprintf(stderr, "best_score=%d best_score_filt=%d\n",
                  best_score,
                  best_score_filt);
                  */
          // check if we should redo ALL the filtered groups
          if(best_score < best_score_filt) {
              // NB: do not redo others...
              for(i=0;i<num_groups;i++) {
                  tmap_map_util_gen_score_t *group = &groups[i];
                  if(0 == group->filtered) continue;
                  // generate the score
                  tmap_map_util_sw_gen_score_helper(refseq, sams, seqs[0], sams_tmp, &j, group->start, group->end,
                                                    group->strand, vsw, seq_len, group->start_pos, group->end_pos,
                                                    &target_mem, &target,
                                                    softclip_start, softclip_end,
                                                    -1, // this is our first call
                                                    opt->max_seed_band, // NB: this may be modified as banding is unrolled
                                                    opt->score_thr-1,
                                                    vsw_opt, rand, opt);
              }
          }
      }
  }

  // realloc
  tmap_map_sams_realloc(sams_tmp, j);

  // sub-optimal score
  for(i=0;i<sams_tmp->n;i++) {
      if(best_subo_score < sams_tmp->sams[i].score_subo) {
          best_subo_score = sams_tmp->sams[i].score_subo;
      }
  }
  // update the sub-optimal
  for(i=0;i<sams_tmp->n;i++) {
      sams_tmp->sams[i].score_subo = best_subo_score;
  }

  // free memory
  tmap_map_sams_destroy(sams);
  free(target);
  tmap_vsw_opt_destroy(vsw_opt);
  tmap_vsw_destroy(vsw);
  free(groups);

  return sams_tmp;
}

static void
tmap_map_util_keytrim(uint8_t *query, int32_t qlen,
                      uint8_t *target, int32_t tlen,
                      int8_t strand, uint8_t key_base,
                      tmap_map_sam_t *s)
{
  int32_t j, k, l;
  int32_t op, op_len;
  //NB: we may only need to look at the first cigar
  op = op_len = 0;
  if(0 == strand) { // forward
      for(j = k = l = 0; j < s->n_cigar; j++) {
          // get the cigar
          op = TMAP_SW_CIGAR_OP(s->cigar[j]);
          op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[j]);
          if(op == BAM_CDEL) break; // looking for mismatch/insertion
          if(op == BAM_CSOFT_CLIP) break; // already start trimming

          while(0 < op_len) {
              if(query[k] != key_base) break;
              if(op == BAM_CINS) {
              }
              else if(op == BAM_CMATCH && target[l] != query[k]) {
                  l++;
              }
              else {
                  break;
              }
              op_len--;
              k++; // since we can only have match/mismatch/insertion
          }
          if(0 == k) {
              // no trimming
              break;
          }
          else if(0 < op_len) {
              if(j == 0) {
                  // reallocate
                  s->n_cigar++;
                  s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*s->n_cigar, "s->cigar");
                  for(k=s->n_cigar-1;0<k;k--) { // shift up
                      s->cigar[k] = s->cigar[k-1];
                  }
                  s->cigar[0] = 0;
                  j++; // reflect the shift in j
              }
              // NB: 0 < j
              // add to the leading soft-clip
              TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[j]) - op_len + TMAP_SW_CIGAR_LENGTH(s->cigar[0]));
              // reduce the cigar length
              TMAP_SW_CIGAR_STORE(s->cigar[j], op, op_len);
              break;
          }
          else { // NB: the full op was removed
              if(0 == j) {
                  TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[0]));
              }
              else {
                  // add to the leading soft-clip
                  TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[j]) + TMAP_SW_CIGAR_LENGTH(s->cigar[0]));
                  // shift down and overwrite the current cigar
                  for(k=j;k<s->n_cigar-1;k++) {
                      s->cigar[k] = s->cigar[k+1];
                  }
                  s->n_cigar--;
                  j--; // since j is incremented later
              }
          }
      }
      // update the position based on the number of reference bases we
      // skipped
      s->pos += l;
  }
  else { // reverse
      for(j = s->n_cigar-1, k = qlen-1, l = tlen-1; 0 <= j; j--) {
          // get the cigar
          op = TMAP_SW_CIGAR_OP(s->cigar[j]);
          op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[j]);
          if(op == BAM_CDEL) break; // looking for mismatch/insertion
          if(op == BAM_CSOFT_CLIP) break; // already start trimming
          while(0 < op_len) {
              if(query[k] != (3 - key_base)) break;
              if(op == BAM_CINS) {
              }
              else if(op == BAM_CMATCH && target[l] != query[k]) {
                  l--;
              }
              else {
                  break;
              }
              op_len--;
              k--; // since we can only have match/mismatch/insertion
          }
          if(qlen-1 == k) {
              // no trimming
              break;
          }
          else if(0 < op_len) {
              if(j == s->n_cigar-1) {
                  // reallocate
                  s->n_cigar++;
                  s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*s->n_cigar, "s->cigar");
                  s->cigar[s->n_cigar-1] = 0;
              }
              // add to the ending soft-clip
              TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar-1], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[j]) - op_len + TMAP_SW_CIGAR_LENGTH(s->cigar[s->n_cigar-1]));
              // reduce the cigar length
              TMAP_SW_CIGAR_STORE(s->cigar[j], op, op_len);
              break;
          }
          else { // NB: the full op was removed
              if(j == s->n_cigar-1) {
                  TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar-1], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[s->n_cigar-1]));
              }
              else {
                  // add to the ending soft-clip
                  TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar-1], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[s->n_cigar-1]));
                  // shift down and overwrite the current cigar
                  for(k=j;k<s->n_cigar-1;k++) {
                      s->cigar[k] = s->cigar[k+1];
                  }
                  s->n_cigar--;
                  j++; // since j is decremented later
              }
          }
      }
  }
}


static void
tmap_map_util_merge_adjacent_cigar_operations(tmap_map_sam_t *s)
{
  int32_t i, j;
  // merge adjacent cigar operations that have the same value
  for(i=s->n_cigar-2;0<=i;i--) {
      if(TMAP_SW_CIGAR_OP(s->cigar[i]) == TMAP_SW_CIGAR_OP(s->cigar[i+1])) {
          TMAP_SW_CIGAR_STORE(s->cigar[i], TMAP_SW_CIGAR_OP(s->cigar[i]), TMAP_SW_CIGAR_LENGTH(s->cigar[i]) + TMAP_SW_CIGAR_LENGTH(s->cigar[i+1]));
          // shift down
          for(j=i+1;j<s->n_cigar-1;j++) {
              s->cigar[j] = s->cigar[j+1];
          }
          s->n_cigar--;
      }
  }
}

//ZZ: finding #S#M#I or #S#M#D alignment for the two sequences allowing at most M mismatch.
static int
tmap_map_util_one_gap(uint8_t *query, int32_t qlen, uint8_t *target, int32_t tlen, int32_t min_size, int max_q, int max_t, int max_MM, int max_indel,
      int8_t strand, int *n_mismatch, int *indel_size, int *is_ins, int *softclip)
{
    int i, j, q, qq, t, tt, inc;
    int qini, tini;
    if (strand == 0) {
  qini = qlen-1; tini = tlen-1;
  inc = -1;
    } else {
  qini = tini = 0; inc = 1;
    }
    if (qlen < min_size || tlen < min_size) return 0;
    if (qlen-min_size < max_q) max_q= qlen-min_size;
    if (tlen-min_size < max_t) max_t = tlen-min_size;
    if (qlen*3 < max_indel) max_indel = qlen*3;
    if (/*max_MM > 1 &&*/ qlen-max_q+1 < 8*max_MM) max_MM = (qlen-max_q+1)/8; // 7-14 bp allow 1 MM, 15 allow 2 etc.
    if (qlen-tlen> max_indel) return 0;
    if (tlen-max_t-qlen > max_indel) return 0;
    //debug
    /*
    char code[] = "ACGT";
    fprintf(stderr, "tlen=%d qlen=%d\n", tlen, qlen);
    for (i = 0; i < tlen; i++) {
  fprintf(stderr, "%c", code[target[i]]);
    }
    fprintf(stderr, "\n");
    for (i = 0; i < qlen; i++) {
  fprintf(stderr,"%c", code[query[i]]);
    }
    fprintf(stderr, "\n");
    */
    for (i = 0, t = tini; i < max_t; i++, t+= inc) {
  int last_M = -1;
  int nM = 0;
  for (j = i, q = 0, qq = qini, tt = t; j < tlen && q < qlen; j++, q++, qq += inc, tt+= inc) {
      if (query[qq] == target[tt]) continue;
      if (q < max_q) last_M = q;
      else {
    nM++;
    if (nM > max_MM) break;
      }
  }
  if (j == tlen || q == qlen) {
      // find one
      *softclip = last_M+1;
      *n_mismatch = nM;
      if (j == tlen) {
    *indel_size = qlen-q;
    *is_ins = 1;
      } else {
    *indel_size = tlen-j;
    *is_ins = 0;
      }
      if (*indel_size > max_indel) return 0;
      if (*indel_size == 0) return 0;
      return 1;
  }
    }
    return 0;
}

static void
tmap_map_util_end_repair(tmap_seq_t *seq, uint8_t *query, int32_t qlen,
                         uint8_t *target_prev, int32_t tlen,
                         int8_t strand,
                         tmap_sw_path_t *path,
                         tmap_refseq_t *refseq,
                         tmap_map_sam_t *s,
                         tmap_map_opt_t *opt, int repair5p)
{
  int32_t i, cigar_i;
  int32_t op, op_len, cur_len;
  int32_t cur_op, cur_op_len, cur_cigar_i, cur_cur_len, target_adj;
  int32_t softclip_start, softclip_end;
  int32_t old_score, new_score;
  int32_t start_pos, end_pos;
  // scoring matrix
  int32_t matrix[25], matrix_iupac[80];
  tmap_sw_param_t par, par_iupac;
  int32_t path_len = 0;
  uint32_t *cigar = NULL;
  int32_t n_cigar = 0;
  int32_t found = 0;

  // TODO: use outside memory
  uint8_t *target = NULL;
  int32_t target_mem = 0;
  int32_t conv = 0;

  // PARAMETERS
  int32_t max_match_len = 5; // the maximum # of query bases to consider for the sub-alignment
  int32_t num_prev_bases = 2; // the maximum number of bases prior to the alignment to include in the target

  // check if we allow soft-clipping on the 5' end
  tmap_map_util_set_softclip(opt, seq, &softclip_start, &softclip_end);
  if (opt->end_repair > 2) {
    if (repair5p == 0 && 1 == softclip_end) return; // allow softclip, then no need to repair.
    if (repair5p == 1 && 1 == softclip_start) return;
    if (repair5p) strand = (!strand);
    old_score = 0;
    found = 0;
    int worst = 0;
    int ind = 0, nMM = 0, /*nc = 0,*/ nqb = 0, ntb = 0, nCC, target_red = 0;
    int tb, qb, inc;
    int i_c;
    int i, nn_cigar=0, total_scl = 0;
    int max_gap = 0, num_gap = 0; // ZZ, for TS-10540
    if (0 != strand) { tb = qb = 0; inc = 1;}
    else {tb = tlen-1; qb = qlen-1; inc = -1;}
    for (i_c = 0; i_c < s->n_cigar; i_c++) {
  cigar_i = (0 != strand) ? i_c : s->n_cigar-1-i_c;
  op = TMAP_SW_CIGAR_OP(s->cigar[cigar_i]);
      op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[cigar_i]);
  if (op == BAM_CMATCH) {
      for(i=0;i<op_len;i++, qb+=inc, tb += inc, nqb++, ntb++) {
    if(query[qb] == target_prev[tb]) old_score += opt->score_match;
            else old_score -= opt->pen_mm, found++;
    if (old_score < worst) { worst = old_score; ind = i+1; nMM = found; target_red =ntb+1; nn_cigar = i_c; total_scl = nqb+1;}
      }
  } else {
      if (op == BAM_CDEL) {
    ntb += op_len;
    tb += inc*op_len;
            } else if (op == BAM_CINS) {
    nqb += op_len;
    qb  += inc*op_len;
      } else break;
      old_score -= opt->pen_gapo+opt->pen_gape*op_len;
      found += op_len;
      if (old_score < worst) { worst = old_score; ind = 0; nMM = found; target_red =ntb; nn_cigar = i_c; total_scl = nqb;
    num_gap++; if (op_len > max_gap) max_gap = op_len; // ZZ, for TS-10540
      }
  }
  if (old_score > 6*opt->score_match) break; // no need go further.
    }
    if (worst >= 0) return;
    nCC = target_red+total_scl+1;
    if (num_gap < 3 && max_gap >= 3) {
  int del = max_gap/15; //increase the error by 1 for each 15 bp of the gap.
  int adj = max_gap-1-del;
  /* int adj = max_gap-1;*/
  nMM -= adj; nCC-= adj;
    }
    nCC /= 2;
    int maxMM = (nCC-2)*opt->end_repair/100+2;
    //if (repair5p) maxMM = (nCC-2)*50/100+2; // temp, try symmetric
    if (maxMM < 2) maxMM = 2; // no need? always at least 2
    if (repair5p && total_scl <= 5) maxMM = nMM+1; // for 5', do not repair short stretch
    if (nCC < maxMM) maxMM = nCC;
    if (nMM >= maxMM) {
  s->score -= worst;
  cigar_i = (0 != strand) ? nn_cigar : s->n_cigar-1-nn_cigar;
  op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[cigar_i]);
  if (ind == op_len || ind == 0) {
      if (nn_cigar == s->n_cigar-1) return; // all are removed, something not quite right. Give up
      if (nn_cigar == s->n_cigar-2) {
    if (TMAP_SW_CIGAR_OP(s->cigar[cigar_i+inc]) == BAM_CSOFT_CLIP) return; // something wrong, will become XSYS??
      }
      TMAP_SW_CIGAR_STORE(s->cigar[cigar_i], BAM_CSOFT_CLIP, total_scl);
      s->n_cigar -= nn_cigar;
  } else {
      TMAP_SW_CIGAR_STORE(s->cigar[cigar_i], BAM_CMATCH, op_len-ind);
      if (nn_cigar == 0) {
    s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
    if (0 != strand) {
                    for (i = s->n_cigar-1; i >= 0; i--) s->cigar[i+1] = s->cigar[i];
                    TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, total_scl);
                } else {
                    TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar],  BAM_CSOFT_CLIP, total_scl);
                }
    s->n_cigar++;
      } else {
    cigar_i -= inc; // move back one
    TMAP_SW_CIGAR_STORE(s->cigar[cigar_i], BAM_CSOFT_CLIP, total_scl);
    s->n_cigar -= nn_cigar-1;
      }
  }
  if (0 != strand && cigar_i > 0) {
      for (i = 0; i < s->n_cigar; i++) {
    s->cigar[i] = s->cigar[i+cigar_i];
      }
  }
  if (0 != strand) s->pos += target_red;
  s->target_len -= target_red;
  // do one large gap label:oneGapZZ
  uint32_t ampl_start, ampl_end;
  if (refseq->bed_exist && total_scl > 5 &&
    tmap_map_get_amplicon(refseq, s->seqid, s->pos, s->pos+s->target_len-1, &ampl_start, &ampl_end, strand)) {
    //fprintf(stderr, "%d %d \n", ampl_start, ampl_end);
    int need_do = 1;
    int half_buffer = /*6*/ opt->amplicon_overrun; // outside freedom.
    int buffer = half_buffer+3;
    uint32_t start, end;
    if (strand == 0) {
        start = s->pos + s->target_len + 1;
        end = ampl_end+half_buffer;
        if (end > refseq->annos[s->seqid].len) end = refseq->annos[s->seqid].len;
        query += qlen-total_scl;
    } else {
        if (ampl_start <= half_buffer) start = 1;
        else start = ampl_start-half_buffer;
        end = s->pos;
    }
    tlen = end-start+1;
    //fprintf(stderr, "%d %d \n", start, end);
    if (tlen > 5) {
         if(target_mem < tlen) { // more memory?
            target_mem = tlen;
            tmap_roundup32(target_mem);
            target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
           }
         if(NULL == tmap_refseq_subseq2(refseq, s->seqid+1, start, end, target, 1, NULL)) {
                  tmap_bug();
                 }
    } else need_do = 0;
    int softclip, target_adj, n_mismatch, indel_size, is_ins;
    if (need_do && tmap_map_util_one_gap(query, total_scl, target, tlen,  opt->min_anchor_large_indel_rescue /*anchor size*/, 2/*int max_q*/, buffer /*int max_t*/, 3 /*maxNM*/, opt->max_one_large_indel_rescue, /*int max_indel*/
                        strand, &n_mismatch, &indel_size, &is_ins, &softclip)) {
        //adjust cigar
        int match_size;
        int target_adj; // how many additional target bases aligned.
        if (is_ins) match_size = total_scl-softclip-indel_size;
        else match_size = total_scl-softclip;
        target_adj = match_size;
        if (!is_ins) target_adj += indel_size;
        int s_adj = (match_size-n_mismatch)*opt->score_match-n_mismatch*opt->pen_mm-opt->pen_gapo-opt->pen_gape;
        s->score += s_adj;
        int inc_cigar = 1;
        if (softclip> 0) inc_cigar++;
        s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(inc_cigar + s->n_cigar), "s->cigar");
        if (strand == 0) {
          TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar-1], is_ins? BAM_CINS: BAM_CDEL , indel_size);
      TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar], BAM_CMATCH, match_size);
      s->n_cigar++;
      if (softclip> 0) {
          TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar],  BAM_CSOFT_CLIP, softclip);
          s->n_cigar++;
      }
        } else {
      for (i = s->n_cigar-1; i >= 0; i--) s->cigar[i+inc_cigar] = s->cigar[i];
      i = 0;
      if (softclip> 0) TMAP_SW_CIGAR_STORE(s->cigar[i++],  BAM_CSOFT_CLIP, softclip);
      TMAP_SW_CIGAR_STORE(s->cigar[i++], BAM_CMATCH, match_size);
      TMAP_SW_CIGAR_STORE(s->cigar[i], is_ins? BAM_CINS: BAM_CDEL , indel_size);
      s->pos -= target_adj;
      s->n_cigar+= inc_cigar;
        }
        s->target_len += target_adj;
        free(target);
    }
  }
    }

    /*
    cigar_i = (0 != strand) ? 0 : (s->n_cigar - 1);
    op = TMAP_SW_CIGAR_OP(s->cigar[cigar_i]);
    op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[cigar_i]);
    if(op != BAM_CMATCH) return; // right now just check terminal matches, will extend to other.
    old_score = 0;
    found = 0; // did we find mismatches?
    int worst = 0;
    cur_len = op_len;
    int ind = 0, nMM = 0;
    if(0 != strand) { // reverse
      for(i=0;i<cur_len;i++) { // go through the match/mismatches
          if(query[i] == target_prev[i]) old_score += opt->score_match;
          else old_score -= opt->pen_mm, found++;
    if (old_score < worst) { worst = old_score; ind = i+1; nMM = found;}
      }
    }
    else { // forward
      for(i=0;i<cur_len;i++) { // go through the match/mismatches
          if(query[qlen-i-1] == target_prev[tlen-i-1]) old_score += opt->score_match;
          else old_score -= opt->pen_mm, found++;
    if (old_score < worst) { worst = old_score; ind = i+1; nMM = found;}
      }
    }
    if (worst >= 0) return; // local alignment is ok
    int maxMM = 2;
    if (ind > 5) maxMM = (ind-2)*opt->end_repair/100+2;
    if (ind < maxMM) maxMM = ind;
    if (nMM >= maxMM) { // trim
  s->score -= worst;
  if (ind == op_len) { // remove the whole cigar
      if (s->n_cigar == 1) return;
      TMAP_SW_CIGAR_STORE(s->cigar[cigar_i], BAM_CSOFT_CLIP, op_len);
  } else {
      TMAP_SW_CIGAR_STORE(s->cigar[cigar_i], BAM_CMATCH, op_len-ind);
      s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
      if (0 != strand) {
    int i;
    for (i = s->n_cigar-1; i >= 0; i--) s->cigar[i+1] = s->cigar[i];
    TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, ind);
      } else {
    TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar],  BAM_CSOFT_CLIP, ind);
      }
      s->n_cigar++;
  }
  if (0 != strand) s->pos += ind;
  s->target_len -= ind;
  //fprintf(stderr, "ZZ tlen%d %d %d %d worst %d MM=%d %d\n", s->target_len, ind, s->score, s->pos, worst, maxMM, nMM);
    }
    */
    return;
  }
  // do not perform if so
  if (repair5p) return;
  if(1 == softclip_start) return;

  // check the first cigar op
  cigar_i = (0 == strand) ? 0 : (s->n_cigar - 1);
  op = TMAP_SW_CIGAR_OP(s->cigar[cigar_i]);
  op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[cigar_i]);
  if(op != BAM_CMATCH) return; // ignore, since we care only about turning mismatches into indels

  // get the amount of bases to adjust
  cur_len = (max_match_len < op_len) ? max_match_len : op_len;
  // NB: must include full HPs
  while(cur_len < qlen && query[cur_len-1] == query[cur_len]) { // TODO: what about reference HPs?
      cur_len++;
  }
  if(cur_len < qlen) {
      cur_len++; // NB: include one more base to get left-justification correct
  }

  // compute the old score of this sub-alignment
  old_score = 0;
  // NB: assumes the op is a match/mismatch
  found = 0; // did we find mismatches?
  if(0 == strand) { // forward
      for(i=0;i<cur_len;i++) { // go through the match/mismatches
          if(query[i] == target_prev[i]) old_score += opt->score_match;
          else old_score -= opt->pen_mm, found = 1;
      }
  }
  else { // reverse
      for(i=0;i<cur_len;i++) { // go through the match/mismatches
          if(query[qlen-i-1] == target_prev[tlen-i-1]) old_score += opt->score_match;
          else old_score -= opt->pen_mm, found = 1;
      }
  }
  if(0 == found) return; // no mismatches, so what are you worrying about?

  // see how many bases we should include from the reference
  target_adj = 0;
  cur_cur_len = cur_len;
  cur_cigar_i = (0 == strand) ? 0 : (s->n_cigar - 1);
  cur_op = TMAP_SW_CIGAR_OP(s->cigar[cur_cigar_i]);
  cur_op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[cur_cigar_i]);
  while(0 < cur_cur_len) { // while read bases
      if(0 == cur_op_len) { // update cigar
          if(0 == strand) {
              cur_cigar_i++;
          }
          else {
              cur_cigar_i--;
          }
          cur_op = TMAP_SW_CIGAR_OP(s->cigar[cur_cigar_i]);
          cur_op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[cur_cigar_i]);
      }

      switch(cur_op) {
        case BAM_CMATCH:
          cur_cur_len--;
          target_adj++;
          break;
        case BAM_CINS:
          cur_cur_len--;
          break;
        case BAM_CDEL:
          target_adj++;
          break;
        default:
          tmap_bug();
          break;
      }
      cur_op_len--;
  }

  // get more target upstream
  if(0 == strand) { // forward
      start_pos = s->pos + 1; // one-based
      if(start_pos < num_prev_bases) start_pos = 1;
      else start_pos -= num_prev_bases;
      end_pos = s->pos + target_adj; // one-based
  }
  else {
      start_pos = s->pos + s->target_len; // 1-based
      if(start_pos <= target_adj) start_pos = 1;
      else start_pos -= target_adj - 1;
      end_pos = s->pos + s->target_len + num_prev_bases;
      if(refseq->annos[s->seqid].len < end_pos) end_pos = refseq->annos[s->seqid].len; // bound
  }
  tlen = end_pos - start_pos + 1;
  if(target_mem < tlen) { // more memory?
      target_mem = tlen;
      tmap_roundup32(target_mem);
      target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
  }
  // NB: IUPAC codes are turned into mismatches
  if(NULL == tmap_refseq_subseq2(refseq, s->seqid+1, start_pos, end_pos, target, 0, &conv)) {
      tmap_bug();
  }

  // set the scoring parameters
  if(0 < conv) {
      par_iupac.matrix = matrix_iupac;
      for(i=0;i<80;i++) {
          if(0 < matrix_iupac_mask[i]) (par_iupac).matrix[i] = (opt)->score_match;
          else (par_iupac).matrix[i] = -cur_len * (opt)->pen_mm;
      }
      (par_iupac).row = 16;
      (par_iupac).band_width = (opt)->bw;

      (par_iupac).gap_open = 0;
      (par_iupac).gap_ext = (opt)->pen_mm;
      (par_iupac).gap_end = (opt)->pen_mm;

  }
  else {
      par.matrix = matrix;
      __map_util_gen_ap(par, opt);
      for(i=0;i<25;i++) {
          (par).matrix[i] = -cur_len * opt->pen_mm; // TODO: is this enough?
      }
      for(i=0;i<4;i++) {
          (par).matrix[i*5+i] = (opt)->score_match;
      }
      (par).row = 5;
      (par).band_width = (opt)->bw;
      (par).gap_open = 0;
      (par).gap_ext = (opt)->pen_mm;
      (par).gap_end = (opt)->pen_mm;
   }

  // adjust query for the reverse strand
  if(1 == strand) { // reverse
      query = query + qlen - cur_len; // NB: do not use this otherwise
  }

  /*
  for(i=0;i<cur_len;i++) {
      fputc("ACGTN"[query[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<tlen;i++) {
      fputc("ACGTN"[target[i]], stderr);
  }
  fputc('\n', stderr);
  */

  // map the target against the the query
  if(0 == strand) {
      new_score = tmap_sw_clipping_core2(target, tlen, query, cur_len,
                                         (0 < conv) ? &par_iupac : &par,
                                         1, 0, // do not skip the end of the target
                                         0, 0, // no soft-clipping
                                         path, &path_len, 0);
  }
  else {
      new_score = tmap_sw_clipping_core2(target, tlen, query, cur_len,
                                         (0 < conv) ? &par_iupac : &par,
                                         0, 1, // do not skip the start of the target
                                         0, 0, // no soft-clipping
                                         path, &path_len, 0);
  }
  //fprintf(stderr, "new_score=%d old_score=%d\n", new_score, old_score);

  found = 0; // reset
  if(old_score < new_score || (2 == opt->end_repair && old_score == new_score)) { // update
      // get the cigar
      cigar = tmap_sw_path2cigar(path, path_len, &n_cigar);
      if(0 == n_cigar) {
          tmap_bug();
      }
      // TODO: what if it starts with an indel

      /*
      if(0 == strand) { // reverse the cigar
          for(i=0;i<n_cigar>>1;i++) {
              int32_t j = cigar[i];
              cigar[i] = cigar[n_cigar-i-1];
              cigar[n_cigar-i-1] = j;
          }
      }
      */

      // check if the cigars match
      if(1 == n_cigar && TMAP_SW_CIGAR_OP(cigar[0]) == BAM_CMATCH) found = 0;
      else found = 1; // do not match, hmmn
  }

  // should we update?
  if(1 == found) {
      /*
      fprintf(stderr, "strand=%d\n", strand);
      fprintf(stderr, "cur_len=%d\n", cur_len);
      fprintf(stderr, "OLD:\n");
      for(i=0;i<s->n_cigar;i++) { // old
          fprintf(stderr, "i=%d %d%c\n", i,
                  bam_cigar_oplen(s->cigar[i]),
                  bam_cigar_opchr(s->cigar[i]));
      }
      fprintf(stderr, "ADDITION:\n");
      for(i=0;i<n_cigar;i++) { // addition
          fprintf(stderr, "i=%d %d%c\n", i,
                  bam_cigar_oplen(cigar[i]),
                  bam_cigar_opchr(cigar[i]));
      }
      */
      // update the previous cigar operators
      int32_t sum = 0, n = 0;
      int32_t num_ref_added = 0, num_ref_removed = 0;
      i = (0 == strand) ? 0 : s->n_cigar-1;
      while((0 == strand && i < s->n_cigar) || (1 == strand && 0 <= i)) { // gets the number to delete, and modifies the last one if necessary in-place
          // query
          switch(TMAP_SW_CIGAR_OP(s->cigar[i])) {
            case BAM_CMATCH:
            case BAM_CINS:
              sum += bam_cigar_oplen(s->cigar[i]);
            default:
              break;
          }
          // target
          switch(TMAP_SW_CIGAR_OP(s->cigar[i])) {
            case BAM_CMATCH:
            case BAM_CDEL:
              if(cur_len <= sum) num_ref_removed += (cur_len - (sum - bam_cigar_oplen(s->cigar[i])));
              else num_ref_removed += bam_cigar_oplen(s->cigar[i]);
            default:
              break;
          }

          //fprintf(stderr, "i=%d cur_len=%d sum=%d\n", i, cur_len, sum);
          if(cur_len <= sum) {
              if(cur_len < sum) { // adjust current cigar
                  TMAP_SW_CIGAR_STORE(s->cigar[i], TMAP_SW_CIGAR_OP(s->cigar[i]), sum - cur_len); // store the difference
                  i = (0 == strand) ? (i-1) : (i+1);
              }
              break;
          }
          i = (0 == strand) ? (i+1) : (i-1);
      }
      n = (0 == strand) ? (i+1) : (s->n_cigar - i); // the number of cigar operations to delete
      //fprintf(stderr, "to delete: %d\n", n);
      if(0 < n) { // delete the cigars
          if(0 == strand) { // shift down, delete the first cigar operator
              for(i=0;i<s->n_cigar-n;i++) {
                  s->cigar[i] = s->cigar[i+n];
              }
          }
          s->n_cigar -= n;
          s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*s->n_cigar, "s->cigar"); // reduce the size
      }
      /*
      fprintf(stderr, "CUT DOWN:\n");
      for(i=0;i<s->n_cigar;i++) { // old
          fprintf(stderr, "i=%d %d%c\n", i,
                  bam_cigar_oplen(s->cigar[i]),
                  bam_cigar_opchr(s->cigar[i]));
      }
      */

      // get more cigar
      s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(n_cigar+s->n_cigar), "s->cigar");
      if(0 == strand) { // forward
          // shift up
          for(i=s->n_cigar-1;0<=i;i--) {
              s->cigar[i+n_cigar] = s->cigar[i];
          }
      } // otherwise append to the end

      // add new cigar and update the position
      for(i=0;i<n_cigar;i++) {
          if(0 == strand) { // forward
              s->cigar[i] = cigar[i];
          }
          else { // reverse
              s->cigar[i + s->n_cigar] = cigar[i];
          }
          switch(TMAP_SW_CIGAR_OP(cigar[i])) {
            case BAM_CMATCH:
            case BAM_CDEL:
              num_ref_added += TMAP_SW_CIGAR_LENGTH(cigar[i]);
            default:
              break;
          }
          /*
          switch(TMAP_SW_CIGAR_OP(cigar[i])) {
            case BAM_CDEL:
              s->target_len -= TMAP_SW_CIGAR_LENGTH(cigar[i]);
              if(0 == strand) s->pos -= TMAP_SW_CIGAR_LENGTH(cigar[i]);
              break;
            case BAM_CINS:
              s->target_len += TMAP_SW_CIGAR_LENGTH(cigar[i]);
              if(0 == strand) s->pos += TMAP_SW_CIGAR_LENGTH(cigar[i]);
              break;
            default:
              break;
          }
          */
      }
      s->n_cigar += n_cigar;
      s->target_len += (num_ref_added - num_ref_removed);
      if(0 == strand) s->pos += (num_ref_removed - num_ref_added);

      // merge adjacent cigar operations that have the same value
      tmap_map_util_merge_adjacent_cigar_operations(s);
      /*
      for(i=0;i<s->n_cigar;i++) {
          fprintf(stderr, "i=%d %d%c\n", i,
                  bam_cigar_oplen(s->cigar[i]),
                  bam_cigar_opchr(s->cigar[i]));
      }
      */

      // update the score
      // TODO: is this correct?
      s->score += new_score;
      s->score -= old_score;

      /*
      fprintf(stderr, "NEW:\n");
      for(i=0;i<s->n_cigar;i++) { // old
          fprintf(stderr, "i=%d %d%c\n", i,
                  bam_cigar_oplen(s->cigar[i]),
                  bam_cigar_opchr(s->cigar[i]));
      }
      */
      // Check that the cigar is valid
      for(i=cur_len=0;i<s->n_cigar;i++) {
          switch(TMAP_SW_CIGAR_OP(s->cigar[i])) { // NB: do not include soft-clip bases
            case BAM_CMATCH:
            case BAM_CINS:
              cur_len += TMAP_SW_CIGAR_LENGTH(s->cigar[i]);
            default:
              break;
          }
      }
      //fprintf(stderr, "cur_len=%d qlen=%d\n", cur_len, qlen);
      if(cur_len != qlen) tmap_bug();
  }

  // free
  free(cigar);
  free(target);
}

void
tmap_map_rmquerybases_cigar(tmap_map_sam_t *s, uint32_t left_rm_bases, uint32_t right_rm_bases)
{


}

tmap_map_sams_t *
tmap_map_util_sw_gen_cigar(tmap_refseq_t *refseq,
                 tmap_map_sams_t *sams,
                 tmap_seq_t *seq,
                 tmap_seq_t **seqs,
                 tmap_map_opt_t *opt)
{
  int32_t i, j, matrix[25], matrix_iupac[80];
  // int32_t start, end;
  int32_t end;
  tmap_map_sams_t *sams_tmp = NULL;
  tmap_sw_param_t par, par_iupac;
  tmap_sw_path_t *path = NULL;
  int32_t path_len, path_mem=0;
  int32_t seq_len=0, tlen, target_mem=0;
  uint8_t *target=NULL;
  tmap_vsw_t *vsw = NULL;
  tmap_vsw_opt_t *vsw_opt = NULL;
  uint32_t start_pos, end_pos;
  int32_t overflow, softclip_start, softclip_end;
  uint8_t key_base = 0;
  int32_t iupac_init = 0;

  if(0 == sams->n) {
      return sams;
  }

  // the final mappings will go here
  sams_tmp = tmap_map_sams_init(sams);
  tmap_map_sams_realloc(sams_tmp, sams->n);

  // scoring matrix
  par.matrix = matrix;
  __map_util_gen_ap(par, opt);

  // sort by strand/chr/pos/score
  tmap_sort_introsort(tmap_map_sam_sort_coord, sams->n, sams->sams);
  tmap_map_util_set_softclip(opt, seq, &softclip_start, &softclip_end);
  if (0 && softclip_end == 0 && refseq->bed_exist) {
    // ZZ: only disable 3' end softclip if the read aligned all the way close to the end of amplicon
    end = 0;
    while(end < sams->n) {
      uint8_t strand; // *query=NULL, *query_rc=NULL, *tmp_target=NULL;
      // int32_t qlen;
      tmap_map_sam_t tmp_sam;
      // int32_t query_start, query_end;
      // int32_t conv = 0;

      // do not band when generating the cigar
      tmp_sam = sams->sams[end];

      // get the strand/start/end positions
      strand = tmp_sam.strand;
      start_pos = tmp_sam.pos + 1;
      end_pos = start_pos + tmp_sam.target_len - 1;
      uint32_t ampl_start, ampl_end;
      if (tmap_map_get_amplicon(refseq, tmp_sam.seqid, start_pos, end_pos, &ampl_start, &ampl_end, strand)) {
  if (strand == 0 && abs(ampl_end-end_pos) < 15) break;
  if (strand != 0 && abs(ampl_start-start_pos) < 15) break;
      }
      end++;
    }
    if (end >= sams->n) softclip_end = 1;
  }

  // initialize opt
  vsw_opt = tmap_vsw_opt_init(opt->score_match, opt->pen_mm, opt->pen_gapo, opt->pen_gape, opt->score_thr);

  // init seqs
  seq_len = tmap_seq_get_bases_length(seqs[0]);

  // reverse compliment query
  vsw = tmap_vsw_init((uint8_t*)tmap_seq_get_bases(seqs[1])->s, seq_len, softclip_end, softclip_start, opt->vsw_type, vsw_opt);

  if(1 == opt->softclip_key) {
      // get first base
      if(NULL == seq->ks) {
          key_base = INT8_MAX;
      }
      else {
          key_base = tmap_nt_char_to_int[(uint8_t)seq->ks[strlen(seq->ks)-1]]; // NB: last base of the key
      }
  }

  // i = start = end = 0;
  i = end = 0;
  start_pos = end_pos = 0;
  while(end < sams->n) {
      uint8_t strand, *query=NULL, *query_rc=NULL, *tmp_target=NULL;
      int32_t qlen;
      tmap_map_sam_t tmp_sam;
      tmap_map_sam_t *s = NULL;
      int32_t query_start, query_end;
      int32_t conv = 0;

      // do not band when generating the cigar
      tmp_sam = sams->sams[end];

      // get the strand/start/end positions
      strand = tmp_sam.strand;
      start_pos = tmp_sam.pos + 1;
      end_pos = start_pos + tmp_sam.target_len - 1;
      uint32_t left_added_bases = 0, right_added_bases = 0; // ZZ to be removed

      /**
       * Step 1: find the start of the alignment
       */

      // retrieve the reverse complimented query sequence
      query_rc = (uint8_t*)tmap_seq_get_bases(seqs[1])->s;
      query_rc += seq_len - tmp_sam.result.query_end - 1; // offset query
      qlen = tmp_sam.result.query_end + 1; // adjust based on the query end

      // get the target sequence
      tlen = tmp_sam.result.target_end + 1; // adjust based on the target end
      if(target_mem < tlen) { // more memory?
          target_mem = tlen;
          tmap_roundup32(target_mem);
          target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }
      // NB: IUPAC codes are turned into mismatches
      if(NULL == tmap_refseq_subseq2(refseq, sams->sams[end].seqid+1, start_pos, end_pos, target, 1, &conv)) {
          tmap_bug();
      }

      // reverse compliment the target
      if(0 == strand) {
          tmap_reverse_compliment_int(target, tlen);
      }

      // NB: if match/mismatch penalties are on the opposite strands, we may
      // have wrong scores
      // NB: this aligns in the opposite direction than sequencing (3'->5')
      tmp_sam.score = tmap_vsw_process_rev(vsw, query_rc, qlen, target, tlen,
                                       &tmp_sam.result, &overflow, opt->score_thr,
                                       (1 ==  softclip_end && 1 == softclip_start) ? 0 : 1); // NB: to guarantee correct soft-clipping if both ends are clipped

      if(1 == overflow) {
          tmap_bug();
      }

      // reverse back compliment the target
      if(0 == strand) {
          tmap_reverse_compliment_int(target, tlen);
      }

      if(tmp_sam.score < opt->score_thr) { // this could happen if VSW fails
          end++;
          continue;
      }

      // ZZ, figure out amplicon condition to add primers or similar
      uint32_t ampl_start = 0, ampl_end = 0;
      uint32_t start_pos_a = start_pos, end_pos_a = end_pos;
      if (strand == 0) {
  start_pos_a += tmp_sam.result.target_start;
      } else {
  end_pos_a = start_pos+ (tmp_sam.result.target_end-tmp_sam.result.target_start);
      }
      if (0 && softclip_end == 0 && tmap_map_get_amplicon(refseq, tmp_sam.seqid, start_pos_a, end_pos_a, &ampl_start, &ampl_end, strand)) {
  uint32_t ds = abs(start_pos_a-ampl_start);
  uint32_t de = abs(end_pos_a-ampl_end);
  if (ds < 10 && de < 10 && (de > 2 || ds > 2)) {
  //if ((strand == 0 && ds < 10 && start_pos_a >= ampl_start+2) || (strand != 0 && de < 10 &&  end_pos_a <= ampl_end-2)) {
  //if ((strand == 0 && ds < 10) || (strand != 0 && de < 10)) {
        //int left = 1,  right = 1;
        //if (strand == 0) left = 0; else right = 0;
        //if (ds < 10 && de < 10 /*&& start_pos_a >= ampl_start+2 && end_pos_a <= ampl_end-2*/) left = right = 0;
        // target mem
        int pad = 20;
        tlen = ampl_end-ampl_start+1+pad*2;
              if(target_mem < tlen) { // more memory?
                  target_mem = tlen;
                  tmap_roundup32(target_mem);
                  target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
              }
        query = (uint8_t*)tmap_seq_get_bases(seqs[strand])->s;
              // Get the new target
              // NB: IUPAC codes are turned into mismatches
              s = &sams_tmp->sams[i];
             (*s) = tmp_sam;

              if(NULL == tmap_refseq_subseq2(refseq, s->seqid+1, ampl_start-pad, ampl_end+pad, target, 0, &conv)) {
                  tmap_bug();
              }
        // finding the digestion sites.
        int x;
        for (x = tlen-pad; x < tlen; x++) {
    if (target[x] == 0) break;
        }
        x--;
        ampl_end = x+(ampl_start-pad);

        for (x = pad-2; x >=0; x--) {
    if (target[x] == 3) break;
        }
        x++;
        tmp_target = target;
        target += x;
        ampl_start = ampl_start-pad+x;
        tlen = ampl_end-ampl_start+1;
              if(0 < conv && 0 == iupac_init) {
                  par_iupac.matrix = matrix_iupac;
                  __map_util_gen_ap_iupac(par_iupac, opt);
                  iupac_init = 1;
              }

        // path memory
        if(path_mem <= tlen + seq_len) { // lengthen the path
            path_mem = tlen + seq_len;
            tmap_roundup32(path_mem);
            path = tmap_realloc(path, sizeof(tmap_sw_path_t)*path_mem, "path");
        }

              s->score = 20*opt->score_match+ tmap_sw_clipping_core2(target, tlen, query, seq_len,
                                            (0 < conv) ? &par_iupac : &par,
                                            0, 0, //whole amplicon must be aligned ZZ
                                            1, 1, /*0, 0,*/ // whole global, may use a better algorithm ZZ
                                            path, &path_len, 0);
        //fprintf(stderr, "ZZ= score=%d %d seqid%d %d %d %d %d %d %d strand=%d %d %d\n", s->score, tmp_sam.score, s->seqid+1, ampl_start, ampl_end, start_pos_a, end_pos_a, tmp_sam.result.target_start, tmp_sam.result.target_end, strand, start_pos, end_pos);
        //if (s->score > tmp_sam.score) {
    int tail_clip =  seq_len-path[0].j;
    int begin_clip = path[path_len-1].j-1;
    if (TMAP_SW_FROM_D == path[path_len-1].ctype) begin_clip++;
          s->cigar = tmap_sw_path2cigar(path, path_len, &s->n_cigar);
    // need to add soft clip, possibly on both sides
          if(0 == s->n_cigar) {
            tmap_bug();
                }
          s->pos = ampl_start-1;
          s->target_len = tlen;
    int f = (begin_clip >0) ? 1:0;
    int e = (tail_clip > 0) ? 1:0;
    if (f+e > 0) {
      s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(f+e+s->n_cigar), "s->cigar");
      if (f > 0) {
          int i;
          for (i = s->n_cigar-1; i>=0; i--) {
        s->cigar[i+1] = s->cigar[i];
          }
          TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, begin_clip);
      }
      if (e > 0) TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar+f+e-1], BAM_CSOFT_CLIP, tail_clip);
      s->n_cigar += e+f;
    }
    // check query length and cigar
    int tttlen = 0;
    for (j = 0; j < s->n_cigar; j++) {
      switch(TMAP_SW_CIGAR_OP(s->cigar[j])) {
                        case BAM_CMATCH:
                        case BAM_CINS:
      case BAM_CSOFT_CLIP:
                          tttlen += TMAP_SW_CIGAR_LENGTH(s->cigar[j]);
        break;
      }
    }
    if (tttlen != seq_len) {
      fprintf(stderr, "%d seqlen=%d tttlen=%d\t tail_skip=%d  begin_clip=%d n_cigar=%d\n", path[path_len-1].ctype, seq_len, tttlen, tail_clip, begin_clip, s->n_cigar);
    }

        // May need adjust start??
          target = tmp_target;
        end++; i++;
        // start = end;
          // update aux data
          tmap_map_sam_malloc_aux(s);
          switch(s->algo_id) {
            case TMAP_MAP_ALGO_MAP1:
                (*s->aux.map1_aux) = (*tmp_sam.aux.map1_aux);
                break;
            case TMAP_MAP_ALGO_MAP2:
                (*s->aux.map2_aux) = (*tmp_sam.aux.map2_aux);
                break;
            case TMAP_MAP_ALGO_MAP3:
                (*s->aux.map3_aux) = (*tmp_sam.aux.map3_aux);

            case TMAP_MAP_ALGO_MAP4:
                (*s->aux.map4_aux) = (*tmp_sam.aux.map4_aux);
                break;
            case TMAP_MAP_ALGO_MAPVSW:
                (*s->aux.map_vsw_aux) = (*tmp_sam.aux.map_vsw_aux);
                break;
            default:
                tmap_bug();
                break;
          }
    continue;
       //}
  }
      }
      /*
      fprintf(stderr, "qlen=%d tlen=%d\n", qlen, tlen);
      fprintf(stderr, "tmp_sam.result.query_start=%d\n", tmp_sam.result.query_start);
      fprintf(stderr, "tmp_sam.result.query_end=%d\n", tmp_sam.result.query_end);
      fprintf(stderr, "tmp_sam.result.target_start=%d\n", tmp_sam.result.target_start);
      fprintf(stderr, "tmp_sam.result.target_end=%d\n", tmp_sam.result.target_end);
      fprintf(stderr, "tmp_sam.result.n_best=%d\n", tmp_sam.result.n_best);
      fprintf(stderr, "tmp_sam.pos=%u\n", tmp_sam.pos);
      */

      // NB: target_start and target_end are inverted for the reverse strand but we
      // only care about their difference.

      // adjust the query and target based off of the start of the alignment
      query = (uint8_t*)tmap_seq_get_bases(seqs[strand])->s; // forward genomic strand
      tmp_target = target;
      if(0 == strand) {
          query_start = tmp_sam.result.query_start;
          query_end = tmp_sam.result.query_end;
          query += tmp_sam.result.query_start; // offset query
          target += tmp_sam.result.target_start;
          tmp_sam.pos += tmp_sam.result.target_start; // keep it zero based
      }
      else {
          // <-------< QE <-------< QS <-------<
          query_start = seq_len - tmp_sam.result.query_end - 1;
          query_end = seq_len - tmp_sam.result.query_start - 1;
          query += seq_len - tmp_sam.result.query_end - 1; // offset query
          // does not affect target offset
          // does not affect tmp_sam.pos; this was updated earlier
      }
      qlen = tmp_sam.result.query_end - tmp_sam.result.query_start + 1; // update query length
      tlen = tmp_sam.result.target_end - tmp_sam.result.target_start + 1;

      /**
       * Step 2: generate the cigar
       */
      if(0 < conv) { // NB: there were IUPAC bases
          // init iupac parameters, if necessary
          if(0 == iupac_init) {
              par_iupac.matrix = matrix_iupac;
              __map_util_gen_ap_iupac(par_iupac, opt);
              iupac_init = 1;
          }
          // update start/end
          start_pos = tmp_sam.pos + 1; // one-based
          end_pos = start_pos + tlen - 1;
          target = tmp_target; // reset target in memory
          if(target_mem < tlen) { // more memory?
              target_mem = tlen;
              tmap_roundup32(target_mem);
              target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
          }
          // Get the new target
          // NB: IUPAC codes are turned into mismatches
          if(NULL == tmap_refseq_subseq2(refseq, sams->sams[end].seqid+1, start_pos, end_pos, target, 0, NULL)) {
              tmap_bug();
          }
          tmp_target = target; // store for later
      }

      // path memory
      if(path_mem <= tlen + seq_len) { // lengthen the path
          path_mem = tlen + seq_len;
          tmap_roundup32(path_mem);
          path = tmap_realloc(path, sizeof(tmap_sw_path_t)*path_mem, "path");
      }

      /*
      // Debugging
      for(j=0;j<qlen;j++) {
          fputc("ACGTN"[query[j]], stderr);
      }
      fputc('\n', stderr);
      for(j=0;j<tlen;j++) {
          fputc(tmap_iupac_int_to_char[target[j]], stderr);
      }
      fputc('\n', stderr);
      */

      // path memory
      if(path_mem <= tlen + qlen) { // lengthen the path
          path_mem = tlen + qlen;
          tmap_roundup32(path_mem);
          path = tmap_realloc(path, sizeof(tmap_sw_path_t)*path_mem, "path");
      }

      s = &sams_tmp->sams[i];

      // shallow copy previous data
      (*s) = tmp_sam;

      // Smith Waterman with banding
      // NB: we store the score from the banded version, which does not allow
      // ins then del, or del then ins.  The vectorized version does.
      // NB: left genomic indel justification is facilitated by always using the
      // forward strand target/query combination.
      // NB: iupac bases may also increase the score
      if(0 < conv) { // NB: there were IUPAC bases
          s->score = tmap_sw_global_banded_core(target, tlen, query, qlen, &par_iupac,
                                                tmp_sam.result.score_fwd, path, &path_len, 0);
      }
      else {
          s->score = tmap_sw_global_banded_core(target, tlen, query, qlen, &par,
                                                tmp_sam.result.score_fwd, path, &path_len, 0);
      }

      // ZZ: There can be case when the global score is the same as the vectorized version
      // whileas the alignments are different. We are not addressing all of them, but the
      // cases where the resulting alignment starts or ends with deletion shall be avoided.
      if( (s->score != tmp_sam.score && 1 < tmp_sam.result.n_best)
    || (path_len > 0 && (TMAP_SW_FROM_D == path[0].ctype || TMAP_SW_FROM_D == path[path_len-1].ctype)))
      {
          /*
          fprintf(stderr, "s->score=%d tmp_sam.score=%d\n", s->score, tmp_sam.score);
          fprintf(stderr, "tmp_sam.result.n_best=%d\n", tmp_sam.result.n_best);
          */

    // TS-9081 (DVK): catch the state of the input parameters
#ifdef CHECK_SWCORE_PARAMS
    // qlen and query
    int32_t saved_query_len = qlen;
    uint8_t* saved_query = alloca (qlen+1);
    memcpy (saved_query, query, qlen); saved_query [qlen] = 0;
    // tlen and target
    int32_t saved_target_len = tlen;
    unit8_t* saved_target = alloca (tlen+1);
    memcpy (saved_target, target, tlen); saved_target [tlen] = 0;
    // SW parameters
    tmap_sw_param_t saved_params;
    memcpy (&saved_params, (0 < conv)? &par_iupac : &par, sizeof (tmap_sw_param_t));
    saved_params.matrix = alloca (saved_params.row * saved_params.row * sizeof (int32_t));
    memcpy (saved_params.matrix, (0 < conv)? par_iupac.matrix : par.matrix, saved_params.row * saved_params.row * sizeof (int32_t));
#endif

          // explicitly fit the query into the target
          if(0 < conv) { // NB: there were IUPAC bases
              s->score = tmap_sw_fitting_core(target, tlen, query, qlen, &par_iupac, path, &path_len, 0);
          }
          else {
              s->score = tmap_sw_fitting_core(target, tlen, query, qlen, &par, path, &path_len, 0);
          }

#ifdef CHECK_SWCORE_PARAMS
    tmap_sw_param_t* pp = (0 < conv)? &par_iupac : &par;
    if (qlen != saved_query_len ||
        tlen != saved_target_len ||
        memcmp (saved_query_len, query, qlen) ||
        memcmp (saved_target_len, target, tlen) ||
        pp->gap_open != saved_params.gap_open ||
        pp->gap_ext != saved_params.gap_ext ||
        pp->gap_end != saved_params.gap_end ||
        pp->row != saved_params.row ||
        pp->band_width != saved_params.gap_width ||
        memcmp (pp->matrix, saved_params.matrix, saved_params.row * saved_params.row * sizeof (int32_t)))
    {
      fprintf (stderr, "\nUpper stack data changed during the function call:\n");
      fprintf (stderr, "Read: %s (%ld bp)\n", qry_name, qry_seq->l);
      fprintf (stderr, "Full Query (%4ld bp): %s\n", qry_seq->l, qry_seq->s);
      fprintf (stderr, "Query      (%4d bp): %s\n", qlen, query);
      fprintf (stderr, "Subject    (%4d bp): ", tlen);
      for (iii = 0; iii != tlen; ++iii)
    fputc ("ACGTN" [(target [iii] < 4) ? target [iii] : 4], stderr);
      fputc ('\n', stderr);
      fprintf (stderr, "score is %d; tmp_sam.score is %d\n", s->score, tmp_sam.score);
      tmap_bug ();
    }

#endif

    int32_t leading, trailing;
    while(1) {
      // check to see if there are leading/trailing deletions, then redo
      leading = trailing = 0;
      // leading
      for(j=path_len-1;0<=j;j--) {
        if(TMAP_SW_FROM_D == path[j].ctype) leading++;
        else break;
      }
      // trailing
      for(j=0;j<path_len;j++) {
        if(TMAP_SW_FROM_D == path[j].ctype) trailing++;
        else break;
      }
      // do we need to continue?
      if(0 == leading && 0 == trailing) break; // no need


 { // TS-9081 (DVK): reporting the erratic (presumably impossible) tail-deletion condition

    unsigned iii;
    const char* qry_name= tmap_seq_get_name (seq)->s;
    tmap_string_t* qry_seq = tmap_seq_get_bases (seq);

    fprintf (stderr, "\nLead/trail delete on read: %s (%ld bp), leading = %d, trailing = %d\n", qry_name, qry_seq->l, leading, trailing);
    fprintf (stderr, "Full Query (%4ld bp): %s\n", qry_seq->l, qry_seq->s);
    fprintf (stderr, "Query      (%4d bp): %s\n", qlen, query);
    // target in alignment
    fprintf (stderr, "Subject    (%4d bp): ", tlen);
    for (iii = 0; iii != tlen; ++iii)
        fputc ("ACGTN" [(target [iii] < 4) ? target [iii] : 4], stderr);
    fputc ('\n', stderr);
    fprintf (stderr, "score is %d; tmp_sam.score is %d\n", s->score, tmp_sam.score);
    fprintf (stderr, "Alignment path (%d positions):\n", path_len);
    for (iii = 0; iii != path_len; ++iii)
      fprintf (stderr, "%3d", path [iii].i);
    fputc ('\n', stderr);
    for (iii = 0; iii != path_len; ++iii)
      fprintf (stderr, "%3c", "MIDSN" [(path [iii].ctype < 4) ? path [iii].ctype : 4]);
    fputc ('\n', stderr);
    for (iii = 0; iii != path_len; ++iii)
      fprintf (stderr, "%3d", path [iii].j);
    fputc ('\n', stderr);
    fprintf (stderr, "Reference Id: %d, strand: %d, position: %d", sams->sams[end].seqid, sams->sams[end].strand, sams->sams[end].pos);

    fprintf (stderr, "n_best: %d\n", tmp_sam.result.n_best);
    fprintf (stderr, "Converted IUPAC bases in reference: %d\n", conv);

    tmap_sw_param_t *pp = (0 == conv)? &par : &par_iupac;
    fprintf (stderr, "SW parameters:\n");
    fprintf (stderr, "    band width = %d\n", pp->band_width); /*!< for Smith-Waterman banding */
    fprintf (stderr, "    gap open = %d\n",  pp->gap_open); /*!< gap open penalty (positive) */
    fprintf (stderr, "    gap_ext = %d\n", pp->gap_ext); /*!< gap extension penalty (positive) */
    fprintf (stderr, "    gap_end = %d\n", pp->gap_end); /*!< gap end penalty (positive */
    fprintf (stderr, "    alphabet size = %d\n", pp->row);
    unsigned rr, cc;
    for (rr = 0; rr != pp->row; ++rr)
    {
      fprintf (stderr, "    %3u ", rr);
      for (cc = 0; cc != pp->row; ++cc)
      {
          fprintf (stderr, "%5d", pp->matrix [pp->row*rr + cc]);
      }
      fprintf (stderr, "\n");
    }
 }

              tmap_bug(); // NB: this should not happen with tmap_sw_fitting_core
          }
          // TODO: should we skip otherwise? No, since VSW can differfrom GSW
      }

      s->pos = s->pos + (path[path_len-1].i-1); // zero-based
      if(path[path_len-1].ctype == TMAP_SW_FROM_I) { // don't quite understand this ZZ
          s->pos++;
      }
      s->cigar = tmap_sw_path2cigar(path, path_len, &s->n_cigar);
      if(0 == s->n_cigar) {
          tmap_bug();
      }

      s->target_len = 0;
      for(j=0;j<s->n_cigar;j++) {
          switch(TMAP_SW_CIGAR_OP(s->cigar[j])) {
            case BAM_CMATCH:
            case BAM_CDEL:
              s->target_len += TMAP_SW_CIGAR_LENGTH(s->cigar[j]);
              break;
            default:
              break;
          }
      }

      /**
       * Try to detect large indels present as soft-clipped prefixes of the
       * alignment.
       */
  uint32_t ampl_exist = 0;
      while(1) { // NB: so we can break out at any time
          if(0 < opt->pen_gapl && 0 < query_start) { // start of the alignment
        // ZZ:We may check amplicon here as well. May try this only when it is near start of an amplicon.
        // Here the reference mapping position is known by s->pos to s->pos+s->target_len
        // We can check the amplicon with these.
        //   A. When bed file is not present, current logic holds.
        //   B.i. If bed file is given, when an amplicon is present, try use all the bases to the begin/end of amplicon, very small or no long gap penalty.
        //   B.ii.If bed file is given and the hit position is not in an amplicon, normal logic apply.

              int32_t del_score = 0;
              int32_t ins_score = 0;
              uint32_t op, op_len;
              int32_t new_score;
              uint32_t *cigar = NULL;
              int32_t n_cigar;
              int32_t pos_adj = 0;

              query = (uint8_t*)tmap_seq_get_bases(seqs[strand])->s; // forward genomic strand

              // get the target sequence before the start of the alignment
              tlen = query_start + opt->gapl_len;
              if(s->pos <= tlen) start_pos = 1;
              else start_pos = s->pos - tlen + 1;
              end_pos = s->pos;
              if(end_pos < start_pos) break; // exit the loop
              tlen = end_pos - start_pos + 1;
        int32_t long_gap = opt->pen_gapl;
    //ZZ: check amplicon
                if (tmap_map_get_amplicon(refseq, s->seqid, s->pos, s->pos+s->target_len, &ampl_start, &ampl_end, strand)) {
                        ampl_exist = 1;
      if (ampl_start >= s->pos) break;
      if (start_pos > ampl_start) {
          start_pos = ampl_start;
          tlen = end_pos - start_pos + 1;
      }
      //Alternative, align only to the begin of the amplicon.
      /*
      start_pos = ampl_start;
      tlen = end_pos - start_pos + 1;
      */
      long_gap = 0;
                }

              // target mem
              target = tmp_target; // reset target in memory
              if(target_mem < tlen) { // more memory?
                  target_mem = tlen;
                  tmap_roundup32(target_mem);
                  target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
                  tmp_target = target; // store for later
              }
              // Get the new target
              // NB: IUPAC codes are turned into mismatches
              if(NULL == tmap_refseq_subseq2(refseq, s->seqid+1, start_pos, end_pos, target, 0, &conv)) {
                  tmap_bug();
              }
              if(0 < conv && 0 == iupac_init) {
                  par_iupac.matrix = matrix_iupac;
                  __map_util_gen_ap_iupac(par_iupac, opt);
                  iupac_init = 1;
              }

              // try to add a larger deletion
              del_score = tmap_sw_clipping_core2(target, tlen, query, query_start,
                                                (0 < conv) ? &par_iupac : &par,
                                                1, 1, // allow the start and end of the target to be skipped
                                                1, 0, // deletions
                                                NULL, 0, 0);

              // try to add a larger insertion
              // TODO: how do we enforce that the target starts immediately where
              // we want it to?
              ins_score = tmap_sw_clipping_core2(target, tlen, query, query_start,
                                                (0 < conv) ? &par_iupac : &par,
                                                1, 0, // only allow the start of the target to be skipped
                                                1, 1, // insertions,
                                                NULL, 0, 0);
              // reset start/end position
              start_pos = s->pos + 1;
              end_pos = s->pos + s->target_len;

              if(del_score <= 0 && ins_score <= 0) {
                  // do nothing
                  new_score = -1;
              }
              else {
                  if(del_score < ins_score) {
                      // get the path
                      tmap_sw_clipping_core2(target, tlen, query, query_start,
                                            (0 < conv) ? &par_iupac : &par,
                                            1, 0, // only allow the start of the target to be skipped
                                            1, 1, // insertions,
                                            path, &path_len, 0);
                      op = BAM_CINS;
                      op_len = query_start - path[0].j;
                      new_score = ins_score;
                      pos_adj = 0;
                  }
                  else {
                      // get the path
                      tmap_sw_clipping_core2(target, tlen, query, query_start,
                                            (0 < conv) ? &par_iupac : &par,
                                            1, 1, // allow the start and end of the target to be skipped
                                            1, 0, // deletions
                                            path, &path_len, 0);
                      op = BAM_CDEL;
                      op_len = tlen - path[0].i;
                      new_score = del_score;
                      pos_adj = op_len;
                  }
                  if(path[path_len-1].ctype == TMAP_SW_FROM_I) pos_adj++; // TODO: is this correct?
              }
              if(0 < new_score && 0 <= new_score - long_gap) {
                  int32_t n_cigar_op = (0 < op_len) ? 1 : 0;
                  // get the cigar
                  cigar = tmap_sw_path2cigar(path, path_len, &n_cigar);
                  if(0 == n_cigar) {
                      tmap_bug();
                  }
                  // re-allocate the cigar
                  s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(n_cigar_op+n_cigar+s->n_cigar), "s->cigar");
                  // shift up
                  for(j=s->n_cigar-1;0<=j;j--) {
                      s->cigar[j+n_cigar_op+n_cigar] = s->cigar[j];
                  }
                  // add the operation
                  if(1 == n_cigar_op) { // 0 < op_len
                      TMAP_SW_CIGAR_STORE(s->cigar[n_cigar], op, op_len);
                  }
                  // add the rest of the cigar
                  for(j=0;j<n_cigar;j++) {
                      s->cigar[j] = cigar[j];
                  }
                  s->n_cigar += n_cigar + n_cigar_op;
                  for(j=0;j<n_cigar;j++) {
                      switch(TMAP_SW_CIGAR_OP(cigar[j])) {
                        case BAM_CMATCH:
                        case BAM_CDEL:
                          s->target_len += TMAP_SW_CIGAR_LENGTH(cigar[j]);
                          s->pos -= TMAP_SW_CIGAR_LENGTH(cigar[j]); // NB: must adjust the position at the start of the alignment
                          break;
                        default:
                          break;
                      }
                  }
                  free(cigar); cigar = NULL;
                  // update the query end
                  query_start = path[path_len-1].j - 1; // query_start is one-based
                  if(query_start < 0) tmap_bug();
                  s->score += new_score - long_gap;
                  s->pos -= pos_adj; // adjust the position further if there was a deletion
                  // merge adjacent cigar operations
                  tmap_map_util_merge_adjacent_cigar_operations(s);
              }

              // reset start/end position
              start_pos = s->pos + 1;
              end_pos = s->pos + s->target_len;
              if(end_pos < start_pos) tmap_bug();

              // retrieve the target
              //TODO: is this always necessary (keytrim)?
              // target mem
              tlen = end_pos - start_pos + 1;
              if(tlen <= 0) tmap_bug();
              target = tmp_target; // reset target in memory
              if(target_mem < tlen) { // more memory?
                  target_mem = tlen;
                  tmap_roundup32(target_mem);
                  target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
                  tmp_target = target; // store for later
              }
              // Get the new target
              // NB: IUPAC codes are turned into mismatches
              if(NULL == tmap_refseq_subseq2(refseq, s->seqid+1, start_pos, end_pos, target, 0, &conv)) {
                  tmap_bug();
              }
              if(0 < conv && 0 == iupac_init) {
                  par_iupac.matrix = matrix_iupac;
                  __map_util_gen_ap_iupac(par_iupac, opt);
                  iupac_init = 1;
              }

              // reset query
              query = (uint8_t*)tmap_seq_get_bases(seqs[strand])->s; // forward genomic strand
              query += query_start;
          }
          break;
      }

      /**
       * Try to detect large indels present as soft-clipped suffixes of the
       * alignment.
       */
      while(1) { // NB: so we can break out at any time
          if(0 < opt->pen_gapl && query_end < seq_len-1) {
              int32_t del_score = 0;
              int32_t ins_score = 0;
              uint32_t op, op_len;
              int32_t new_score = 0;
              uint32_t *cigar = NULL;
              int32_t n_cigar;

              // get the target sequence after the end of the alignment.
              tlen = seq_len - 1 - query_end + opt->gapl_len;
              start_pos = s->pos + s->target_len + 1;
              end_pos = start_pos + tlen - 1;
              if(refseq->annos[s->seqid].len < end_pos) { // bound
                  end_pos = refseq->annos[s->seqid].len; // one-based
              }
              if(end_pos < start_pos) break; // exit the loop
              tlen = end_pos - start_pos + 1;
        int32_t long_gap = opt->pen_gapl;

              //ZZ: check amplicon
                if (ampl_exist || tmap_map_get_amplicon(refseq, s->seqid, s->pos, s->pos+s->target_len, &ampl_start, &ampl_end, strand)) {
                        ampl_exist = 1;
                        if (ampl_end <= s->pos+s->target_len + 1) break;
                        if (end_pos < ampl_end) {
                            end_pos = ampl_end;
                            tlen = end_pos - start_pos + 1;
                        }
                        //Alternative, align only to the end of the amplicon.
                        /*
                        end_pos = ampl_end;
                        tlen = end_pos - start_pos + 1;
                        */
                        long_gap = 0;
                }

              // target mem
              target = tmp_target; // reset target in memory
              if(target_mem < tlen) { // more memory?
                  target_mem = tlen;
                  tmap_roundup32(target_mem);
                  target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
                  tmp_target = target; // store for later
              }
              // Get the new target
              // NB: IUPAC codes are turned into mismatches
              if(NULL == tmap_refseq_subseq2(refseq, s->seqid+1, start_pos, end_pos, target, 0, &conv)) {
                  tmap_bug();
              }
              if(0 < conv && 0 == iupac_init) {
                  par_iupac.matrix = matrix_iupac;
                  __map_util_gen_ap_iupac(par_iupac, opt);
                  iupac_init = 1;
              }

              // try to add a larger deletion
              del_score = tmap_sw_clipping_core2(target, tlen, query + query_end + 1, seq_len - query_end - 1,
                                                (0 < conv) ? &par_iupac : &par,
                                                1, 1, // allow the start and end of the target to be skipped
                                                0, 1, // deletions
                                                NULL, 0, 0);

              // try to add a larger insertion
              ins_score = tmap_sw_clipping_core2(target, tlen, query + query_end + 1, seq_len - query_end - 1,
                                                (0 < conv) ? &par_iupac : &par,
                                                0, 1, // only allow the end of the target to be skipped
                                                1, 1, // insertions,
                                                NULL, 0, 0);

              if(del_score <= 0 && ins_score <= 0) {
                  // do nothing
                  new_score = -1;
              }
              else {
                  if(del_score < ins_score) {
                      // get the path
                      tmap_sw_clipping_core2(target, tlen, query + query_end + 1, seq_len - query_end - 1,
                                            (0 < conv) ? &par_iupac : &par,
                                            0, 1, // only allow the end of the target to be skipped
                                            1, 1, // insertions,
                                            path, &path_len, 0);
                      op = BAM_CINS;
                      op_len = path[path_len-1].j - 1;
                      new_score = ins_score;
                  }
                  else {
                      // get the path
                      tmap_sw_clipping_core2(target, tlen, query + query_end + 1, seq_len - query_end - 1,
                                            (0 < conv) ? &par_iupac : &par,
                                            1, 1, // allow the start and end of the target to be skipped
                                            0, 1, // deletions
                                            path, &path_len, 0);
                      op = BAM_CDEL;
                      op_len = path[path_len-1].i - 1;
                      new_score = del_score;
                  }
              }
              if(0 < new_score && 0 <= new_score - long_gap) {
                  int32_t n_cigar_op = (0 < op_len) ? 1 : 0;
                  // get the cigar
                  cigar = tmap_sw_path2cigar(path, path_len, &n_cigar);
                  if(0 == n_cigar) {
                      tmap_bug();
                  }
                  // re-allocate the cigar
                  s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(n_cigar_op+n_cigar+s->n_cigar), "s->cigar");
                  // add the operation
                  if(1 == n_cigar_op) { // 0 < op_len
                      TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar], op, op_len);
                  }
                  // add the rest of the cigar
                  for(j=0;j<n_cigar;j++) {
                      s->cigar[j+s->n_cigar+n_cigar_op] = cigar[j];
                  }
                  s->n_cigar += n_cigar + n_cigar_op;
                  for(j=0;j<n_cigar;j++) {
                      switch(TMAP_SW_CIGAR_OP(cigar[j])) {
                        case BAM_CMATCH:
                        case BAM_CDEL:
                          s->target_len += TMAP_SW_CIGAR_LENGTH(cigar[j]);
                          break;
                        default:
                          break;
                      }
                  }
                  free(cigar); cigar = NULL;
                  // update the query end
                  query_end += path[0].j; // query_end is zero-based
                  if(seq_len <= query_end) tmap_bug();
                  s->score += new_score - long_gap;
                  // merge adjacent cigar operations
                  tmap_map_util_merge_adjacent_cigar_operations(s);
              }

              // reset start/end position
              start_pos = s->pos + 1;
              end_pos = s->pos + s->target_len;

              // retrieve the target
              //TODO: is this always necessary (keytrim)?
              // target mem
              tlen = end_pos - start_pos + 1;
              target = tmp_target; // reset target in memory
              if(target_mem < tlen) { // more memory?
                  target_mem = tlen;
                  tmap_roundup32(target_mem);
                  target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
                  tmp_target = target; // store for later
              }
              // Get the new target
              // NB: IUPAC codes are turned into mismatches
              if(NULL == tmap_refseq_subseq2(refseq, s->seqid+1, start_pos, end_pos, target, 0, &conv)) {
                  tmap_bug();
              }
              if(0 < conv && 0 == iupac_init) {
                  par_iupac.matrix = matrix_iupac;
                  __map_util_gen_ap_iupac(par_iupac, opt);
                  iupac_init = 1;
              }
          }
          break;
      }

      // add soft clipping
      if(0 < query_start) {
          // soft clip the front of the read
          s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1+s->n_cigar), "s->cigar");
          for(j=s->n_cigar-1;0<=j;j--) { // shift up
              s->cigar[j+1] = s->cigar[j];
          }
          TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, query_start);
          s->n_cigar++;
      }
      if(query_end < seq_len-1) {
          // soft clip the end of the read
          s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1+s->n_cigar), "s->cigar");
          TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar], BAM_CSOFT_CLIP, seq_len - query_end - 1);
          s->n_cigar++;
      }

      // update aux data
      tmap_map_sam_malloc_aux(s);
      switch(s->algo_id) {
        case TMAP_MAP_ALGO_MAP1:
          (*s->aux.map1_aux) = (*tmp_sam.aux.map1_aux);
          break;
        case TMAP_MAP_ALGO_MAP2:
          (*s->aux.map2_aux) = (*tmp_sam.aux.map2_aux);
          break;
        case TMAP_MAP_ALGO_MAP3:
          (*s->aux.map3_aux) = (*tmp_sam.aux.map3_aux);
          break;
        case TMAP_MAP_ALGO_MAP4:
          (*s->aux.map4_aux) = (*tmp_sam.aux.map4_aux);
          break;
        case TMAP_MAP_ALGO_MAPVSW:
          (*s->aux.map_vsw_aux) = (*tmp_sam.aux.map_vsw_aux);
          break;
        default:
          tmap_bug();
          break;
      }

      // key trim the data
      if(1 == opt->softclip_key && INT8_MAX != key_base) {
          // get the new target
          tmap_map_util_keytrim(query, qlen, target, tlen, strand, key_base, s);
      }

      // end repair
      if(0 != opt->end_repair) {
          tmap_map_util_end_repair(seq, query, qlen, target, tlen, strand, path, refseq, s, opt, 0);
    tmap_map_util_end_repair(seq, query, qlen, target, tlen, strand, path, refseq, s, opt, 1); // do 5' end repair??
      }
      // TODO: ZZ: If the primer sequences are added, we remove the cigar that corresponds to the sequence, at the same time calculate the score
      //           of this part of alignment, this will be substracted from the alignment score there.
      // to be removed if the other path works.
      if (left_added_bases + right_added_bases > 0) {
  tmap_map_rmquerybases_cigar(s, left_added_bases, right_added_bases);
      }

      i++;

      // update start/end
      end++;
      // start = end;
      target = tmp_target;
  }

  // realloc
  tmap_map_sams_realloc(sams_tmp, i);

  // free memory
  tmap_map_sams_destroy(sams);
  free(path);
  free(target);
  tmap_vsw_destroy(vsw);
  tmap_vsw_opt_destroy(vsw_opt);

  // sort by max score, then min coordinate
  if(1 < sams_tmp->n) {
      tmap_sort_introsort(tmap_map_sam_sort_score_coord, sams_tmp->n, sams_tmp->sams);
  }

  return sams_tmp;
}

// TODO: make sure the "longest" read alignment is found
int32_t
tmap_map_util_fsw(tmap_seq_t *seq, tmap_map_sams_t *sams, tmap_refseq_t *refseq,
                  int32_t bw, int32_t softclip_type, int32_t score_thr,
                  int32_t score_match, int32_t pen_mm, int32_t pen_gapo,
                  int32_t pen_gape, int32_t fscore, int32_t use_flowgram)
{
  int32_t i, j, k, l;
  uint8_t *target = NULL;
  int32_t target_mem = 0, target_len = 0;
  int32_t was_int = 1;

  tmap_fsw_flowseq_t *fseq = NULL;
  tmap_fsw_path_t *path = NULL;
  int32_t path_mem = 0, path_len = 0;
  tmap_fsw_param_t param;
  int32_t matrix[25];
  uint8_t *flow_order = NULL;
  int32_t flow_order_len = 0;
  uint8_t *key_seq = NULL;
  int32_t key_seq_len = 0;

  if(0 == sams->n || NULL == seq->fo || NULL == seq->ks) return 0;

  // flow order
  flow_order_len = strlen(seq->fo);
  flow_order = tmap_malloc(sizeof(uint8_t)*(flow_order_len+1), "flow_order");
  memcpy(flow_order, seq->fo, flow_order_len);
  tmap_to_int((char*)flow_order, flow_order_len);

  // key sequence
  key_seq_len = strlen(seq->ks);
  key_seq = tmap_malloc(sizeof(uint8_t)*(key_seq_len+1), "key_seq");
  memcpy(key_seq, seq->ks, key_seq_len);
  tmap_to_int((char*)key_seq, key_seq_len);

  // generate the alignment parameters
  param.matrix = matrix;
  param.band_width = 0;
  param.offset = TMAP_MAP_OPT_FSW_OFFSET; // this sets the hp difference
  __tmap_fsw_gen_ap1(param, score_match, pen_mm, pen_gapo, pen_gape, fscore);

  was_int = tmap_seq_is_int(seq);
  if(0 == tmap_seq_is_int(seq)) {
      tmap_seq_to_int(seq);
  }

  // get flow sequence
  fseq = tmap_fsw_flowseq_from_seq(NULL, seq, flow_order, flow_order_len, key_seq, key_seq_len, use_flowgram);

  if(NULL == fseq) {
      free(flow_order);
      free(key_seq);
      return 0;
  }

  // go through each hit
  for(i=0;i<sams->n;i++) {
      tmap_map_sam_t *s = &sams->sams[i];
      uint32_t ref_start, ref_end;

      // get the reference end position
      // NB: soft-clipping at the start may cause ref_start to be moved
      param.band_width = 0;
      ref_start = ref_end = s->pos + 1;
      for(j=0;j<s->n_cigar;j++) {
          int32_t op, op_len;

          op = TMAP_SW_CIGAR_OP(s->cigar[j]);
          op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[j]);

          switch(op) {
            case BAM_CMATCH:
              ref_end += op_len;
              break;
            case BAM_CDEL:
              if(param.band_width < op_len) param.band_width += op_len;
              ref_end += op_len;
              break;
            case BAM_CINS:
              if(param.band_width < op_len) param.band_width += op_len;
              break;
            case BAM_CSOFT_CLIP:
              if(0 == j) {
                  if(ref_start <= op_len) {
                      ref_start = 1;
                  }
                  else {
                      ref_start = ref_start - op_len;
                  }
              }
              else ref_end += op_len;
              break;
            default:
              // ignore
              break;
          }
      }
      ref_end--; // NB: since we want this to be one-based

      // check bounds
      if(ref_start < 1) ref_start = 1;
      if(refseq->annos[s->seqid].len < ref_end) {
          ref_end = refseq->annos[s->seqid].len;
      }
      else if(ref_end < 1) ref_end = 1;

      // get the target sequence
      target_len = ref_end - ref_start + 1;
      if(target_mem < target_len) {
          target_mem = target_len;
          tmap_roundup32(target_mem);
          target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }
      target_len = tmap_refseq_subseq(refseq, ref_start + refseq->annos[s->seqid].offset, target_len, target);
      /*
      // NB: IUPAC codes are turned into mismatches
      if(NULL == tmap_refseq_subseq2(refseq, sams->sams[end].seqid+1, start_pos, end_pos, target, 1, NULL)) {
          tmap_bug();
      }
      */

      if(1 == s->strand) { // reverse compliment
          tmap_reverse_compliment_int(target, target_len);
      }

      // add to the band width
      param.band_width += 2 * bw;

      // make sure we have enough memory for the path
      while(path_mem <= target_len + fseq->num_flows) { // lengthen the path
          path_mem = target_len + fseq->num_flows + 1;
          tmap_roundup32(path_mem);
          path = tmap_realloc(path, sizeof(tmap_fsw_path_t)*path_mem, "path");
      }

      /*
      fprintf(stderr, "strand=%d\n", s->strand);
      fprintf(stderr, "ref_start=%d ref_end=%d\n", ref_start, ref_end);
      fprintf(stderr, "base_calls:\n");
      for(j=0;j<fseq->num_flows;j++) {
          for(k=0;k<fseq->base_calls[j];k++) {
              fputc("ACGTN"[fseq->flow_order[j % fseq->flow_order_len]], stderr);
          }
      }
      fputc('\n', stderr);
      fprintf(stderr, "target:\n");
      for(j=0;j<target_len;j++) {
          fputc("ACGTN"[target[j]], stderr);
      }
      fputc('\n', stderr);
      for(j=0;j<fseq->flow_order_len;j++) {
          fputc("ACGTN"[fseq->flow_order[j]], stderr);
      }
      fputc('\n', stderr);
      */

      // re-align
      s->ascore = s->score;
      path_len = path_mem;
      //fprintf(stderr, "old score=%d\n", s->score);
      switch(softclip_type) {
        case TMAP_MAP_OPT_SOFT_CLIP_ALL:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq, &param,
                                            1, 1, s->strand, path, &path_len);
          break;
        case TMAP_MAP_OPT_SOFT_CLIP_LEFT:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq, &param,
                                            1, 0, s->strand, path, &path_len);
          break;
        case TMAP_MAP_OPT_SOFT_CLIP_RIGHT:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq, &param,
                                            0, 1, s->strand, path, &path_len);
          break;
        case TMAP_MAP_OPT_SOFT_CLIP_NONE:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq, &param,
                                            0, 0, s->strand, path, &path_len);
          break;
        default:
          tmap_error("soft clipping type was not recognized", Exit, OutOfRange);
          break;
      }
      s->score_subo = INT32_MIN;
      //fprintf(stderr, "new score=%d path_len=%d\n", s->score, path_len);

      if(0 < path_len) { // update

          /*
          for(j=0;j<path_len;j++) {
              fprintf(stderr, "j=%d path[j].i=%d path[j].j=%d path[j].type=%d\n", j, path[j].i, path[j].j, path[j].ctype);
          }
          */

          // score
          s->score = (int32_t)((s->score + 99.99)/100.0);

          // position
          s->pos = (ref_start-1); // zero-based
          // NB: must be careful of leading insertions and strandedness
          if(0 == s->strand) {
              if(0 <= path[path_len-1].j) {
                  s->pos += (path[path_len-1].j);
              }
          }
          else {
              if(path[0].j < target_len) {
                  s->pos += target_len - path[0].j - 1;
              }
          }
          if(refseq->len < s->pos) {
              tmap_bug();
          }


          // new cigar
          free(s->cigar);
          s->cigar = tmap_fsw_path2cigar(path, path_len, &s->n_cigar, 1);

          // reverse the cigar
          if(1 == s->strand) {
              for(i=0;i<s->n_cigar>>1;i++) {
                  uint32_t t = s->cigar[i];
                  s->cigar[i] = s->cigar[s->n_cigar-i-1];
                  s->cigar[s->n_cigar-i-1] = t;
              }
          }

          int32_t skipped_start, skipped_end;
          skipped_start = skipped_end = 0;

          // soft-clipping
          if(0 < path[path_len-1].i) { // skipped beginning flows
              // get the number of bases to clip
              for(j=0;j<path[path_len-1].i;j++) {
                  skipped_start += fseq->base_calls[j];
              }
          }
          if(path[0].i+1 < fseq->num_flows) { // skipped ending flows
              // get the number of bases to clip
              for(j=path[0].i+1;j<fseq->num_flows;j++) {
                  skipped_end += fseq->base_calls[j];
              }
          }
          if(1 == s->strand) { // swap
              k = skipped_start;
              skipped_start = skipped_end;
              skipped_end = k;
          }
          if(0 < skipped_start) { // start soft clip
              s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
              for(l=s->n_cigar-1;0<=l;l--) {
                  s->cigar[l+1] = s->cigar[l];
              }
              TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, skipped_start);
              s->n_cigar++;
          }
          if(0 < skipped_end) { // end soft clip
              s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
              s->cigar[s->n_cigar] = (skipped_end << 4) | 4;
              s->n_cigar++;
          }
      }
  }
  // free
  free(target);
  free(path);
  free(flow_order);
  free(key_seq);

  if(0 == was_int) {
      tmap_seq_to_char(seq);
  }
  return 1;
}
