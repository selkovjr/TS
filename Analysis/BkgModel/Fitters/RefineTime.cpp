/* Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved */
#include "RefineTime.h"
#include "DNTPRiseModel.h"
#include "BkgModSingleFlowFitTMidNuc.h"
#include "TraceCorrector.h"

RefineTime::RefineTime (SignalProcessingMasterFitter &_bkg) :
    bkg (_bkg)
{

}


void RefineTime::RefinePerFlowTimeEstimate (float *t_mid_nuc_shift_per_flow)
{
  if (!bkg.global_defaults.signal_process_control.per_flow_t_mid_nuc_tracking)
  {
    memset (t_mid_nuc_shift_per_flow,0,sizeof (float[NUMFB]));
  }
  else
  {
    if (bkg.region_data->my_trace.AlreadyAdjusted())
    {
      // estimate t0 per flow from the average 1-mer trace if possible
      FitAverage1MerPerFlow (t_mid_nuc_shift_per_flow,false);

      // really? does this help in the slightest?
      RezeroUsingLocalShift(t_mid_nuc_shift_per_flow);

      // really?
      // re-estimate t0 per flow from the average 1-mer trace if possible
      FitAverage1MerPerFlow (t_mid_nuc_shift_per_flow,true);
    }
    else
    {
      // do version that doesn't require pre-corrected background
      // @TODO: note that t_mid_nuc_shift has not been set to anything but zero yet?
      //@TODO: SIDE EFFECTS of changing the traces
      RezeroUsingLocalShift(t_mid_nuc_shift_per_flow);

      FitAverage1MerAllFlows(t_mid_nuc_shift_per_flow,true);
    }
  }
}

void RefineTime::RezeroUsingLocalShift(float *t_mid_nuc_shift_per_flow)
{
      // re-zero the traces now that we have a better idea where t0 is in each flow
    for (int fnum=0;fnum < NUMFB;fnum++)
    {
      bkg.region_data->RezeroTraces (bkg.region_data->time_c.time_start,GetTypicalMidNucTime (&bkg.region_data->my_regions.rp.nuc_shape) +t_mid_nuc_shift_per_flow[fnum],
                                      bkg.region_data->my_regions.rp.nuc_shape.sigma,MAGIC_OFFSET_FOR_EMPTY_TRACE,fnum);
    }
}


#define MIN_1MER_TMIDNUC_EST_THRESHOLD (50)
float RefineTime::FitAverage1Mer (int fnum,bool debug_output)
{
  int len = bkg.region_data->time_c.npts();
  float avg_1mer[len];
  bead_params avg_bead;

  int trc_cnt = FindAvg1MerFromSingleFlowAdjustedData(avg_1mer,len, avg_bead,fnum);

  return(FitSingleFlowTimeShiftFromOneMer(avg_1mer,len, trc_cnt,&avg_bead, fnum, debug_output));
}

int RefineTime::FindAvg1MerFromSingleFlowAdjustedData(float *avg_1mer, int len, bead_params &avg_bead, int fnum)
{
    int trc_cnt = 0;
  float ampl;

  memset (avg_1mer,0,sizeof (float[len]));

  params_SetBeadZeroValue(&avg_bead);

  for (int ibd=0;ibd < bkg.region_data->my_beads.numLBeads;ibd++)
  {
    bead_params *beadp = &bkg.region_data->my_beads.params_nn[ibd];
    ampl = beadp->Ampl[fnum];

    if ( ((ampl > 0.5f) && (ampl < 1.5f)) & FitBeadLogic(beadp)) // must be willing to fit this bead to include it
    {
      bkg.region_data->my_trace.AccumulateSignal (avg_1mer,ibd,fnum,bkg.region_data->time_c.npts());
      params_AccumulateBeadValue(&avg_bead,beadp);

      trc_cnt++;
    }
  }
  if (trc_cnt>0)
  {
    MultiplyVectorByScalar (avg_1mer,1.0f/ ( (float) trc_cnt),bkg.region_data->time_c.npts());
    params_ScaleBeadValue(&avg_bead, 1.0f/float(trc_cnt));
  }
  
  return(trc_cnt);
}


float RefineTime::FitSingleFlowTimeShiftFromOneMer(float *avg_1mer, int len, int trc_cnt, bead_params *avg_bead, int fnum, bool debug_output)
{
  float ret = 0.0f;
  float new_kmult=1.0f;
  // make sure we at least have enough to form a good average before proceeding
  if (trc_cnt > MIN_1MER_TMIDNUC_EST_THRESHOLD)
  {
    float evect[bkg.region_data->time_c.npts()];
    bkg.region_data->emphasis_data.CustomEmphasis (evect,1.0f);
    int NucID = bkg.region_data->my_flow.flow_ndx_map[fnum];
    BkgModSingleFlowFitTMidNucParams param_min,param_max;

    param_min.Ampl = 0.5f;
    param_max.Ampl = 1.5f;
    param_min.delta_t_mid_nuc = -3.0f;
    param_max.delta_t_mid_nuc = 3.0f;
    param_min.kmult = 0.2f;
    param_max.kmult = 2.0f;

    BkgModSingleFlowFitTMidNuc t_mid_nuc_fit (bkg.region_data->time_c.npts(),
        &bkg.region_data->time_c.frameNumber[0],
        &bkg.region_data->time_c.deltaFrame[0],
        &bkg.region_data->time_c.deltaFrameSeconds[0],
        bkg.math_poiss);

    // TODO: get average dmult and kmult for this call
    t_mid_nuc_fit.SetWellRegionParams (avg_bead->Copies,avg_bead->R,avg_bead->gain,1.0f,1.0f,&bkg.region_data->my_regions.rp,fnum,NucID,bkg.region_data->my_flow.buff_flow[fnum]);

    t_mid_nuc_fit.SetParamMin (param_min);
    t_mid_nuc_fit.SetParamMax (param_max);

    t_mid_nuc_fit.Fit (100,avg_1mer);

//   printf("Flow %d, found t_mid_nuc shift of %8.5f\n",my_flow.buff_flow[fnum],t0fit.params.delta_t_mid_nuc);
    ret = t_mid_nuc_fit.params.delta_t_mid_nuc;
    new_kmult = t_mid_nuc_fit.params.kmult;
  }
  DebugOneMerOutput(avg_1mer, len, trc_cnt, ret,new_kmult,debug_output);

  return(ret);
}

void RefineTime::DebugOneMerOutput(float *avg_1mer, int len, int trc_cnt, float delta_mid_nuc, float new_kmult, bool debug_output)
{
    if ( (debug_output) && bkg.my_debug.region_1mer_trace_file)
  {
   if (trc_cnt < 10)
     memset (avg_1mer,0,sizeof (float[len]));
    char sep_char = '\t';
    // I'm going to cheat and stuff the t0 offset into the first value of the average 1-mer trace
    fprintf (bkg.my_debug.region_1mer_trace_file,"%8.5f%c",delta_mid_nuc,sep_char);
    fprintf (bkg.my_debug.region_1mer_trace_file,"%8.5f%c",new_kmult,sep_char);
    for (int i=2;i < bkg.region_data->time_c.npts();i++)
    {
      if (i== (bkg.region_data->time_c.npts()-1))
        sep_char='\n';
      fprintf (bkg.my_debug.region_1mer_trace_file,"%8.5f%c",avg_1mer[i],sep_char);
    }
  }
}

void RefineTime::FitAverage1MerPerFlow (float *t_mid_nuc_shift_per_flow,bool debug_output)
{
  for (int fnum=0;fnum < NUMFB;fnum++)
  {
    t_mid_nuc_shift_per_flow[fnum] = FitAverage1Mer (fnum,debug_output);
  }
}

// do this using unadjusted data
void RefineTime::FitAverage1MerAllFlows(float *t_mid_nuc_shift_per_flow, bool debug_output)
{
  float block_bkg_corrected_avg_signal[bkg.region_data->my_scratch.bead_flow_t];
  float *avg_1_mer = NULL;
  int cur_count[NUMFB];
  bead_params avg_bead_by_flow[NUMFB];
  int local_offset = bkg.region_data->time_c.npts();

  ConstructMultiFlowOneMers(block_bkg_corrected_avg_signal,cur_count,avg_bead_by_flow);
  for (int fnum=0; fnum<NUMFB; fnum++)
  {
    avg_1_mer = &block_bkg_corrected_avg_signal[fnum*local_offset];
    t_mid_nuc_shift_per_flow[fnum] = FitSingleFlowTimeShiftFromOneMer(avg_1_mer,local_offset, cur_count[fnum],&avg_bead_by_flow[fnum], fnum, debug_output);
  }
}

void RefineTime::ConstructMultiFlowOneMers(float *block_bkg_corrected_avg_signal, int *cur_count, bead_params *avg_bead_by_flow)
{
//  int cur_count[NUMFB]; // entries per flow of "1-mer" signals
//  bead_params avg_bead_by_flow[NUMFB]; //goofy, why would these be appreciably different?

  // assume uncorrected data in buffers
  // assume we're building traces that are background corrected
  int buffer_size = bkg.region_data->my_scratch.bead_flow_t;
  int local_offset = bkg.region_data->time_c.npts();
  float cur_bkg_corrected_signal[buffer_size];

  memset(block_bkg_corrected_avg_signal,0, sizeof(float[buffer_size]));

  for (int fnum=0; fnum<NUMFB; fnum++)
    cur_count[fnum] = 0;
  // typical beads had >better< be exactly what we expect, but just in case
  for (int fnum=0; fnum<NUMFB; fnum++)
    params_SetBeadZeroValue(&avg_bead_by_flow[fnum]);

 int beadmax = bkg.region_data->GetNumLiveBeads();
 for (int ibd=0; ibd<beadmax; ibd++) 
  {
    bead_params *p = &bkg.region_data->my_beads.params_nn[ibd];
    if (FitBeadLogic(p))  // might be a sampling of some type
    {
      //@TODO: should be able to average un-adjusted data, then adjust the average
      bkg.trace_bkg_adj->ReturnBackgroundCorrectedSignal(cur_bkg_corrected_signal,ibd);
      for (int fnum=0; fnum<NUMFB; fnum++)
      {
        if ((p->Ampl[fnum] > 0.5f) && (p->Ampl[fnum] < 1.5f))
        {
          AccumulateVector(&block_bkg_corrected_avg_signal[fnum*local_offset],&cur_bkg_corrected_signal[fnum*local_offset],local_offset);
          cur_count[fnum] +=1;
          params_AccumulateBeadValue(&avg_bead_by_flow[fnum], p);
        }
      }
    }
  }
  // iterated over everything
  for (int fnum=0; fnum<NUMFB; fnum++)
  {
    if (cur_count[fnum]>0)
    {
      MultiplyVectorByScalar(&block_bkg_corrected_avg_signal[fnum*local_offset],1.0f/(cur_count[fnum]),local_offset);
      params_ScaleBeadValue(&avg_bead_by_flow[fnum],1.0f/(cur_count[fnum]));
    }
  }
  // okay, set up with average bead, set up with block signal and counts
}

