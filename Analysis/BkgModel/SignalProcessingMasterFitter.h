/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef SIGNALPROCESSINGMASTERFITTER_H
#define SIGNALPROCESSINGMASTERFITTER_H

#include <stdio.h>
#include <set>
#include <vector>

#include "Image.h"
#include "BkgFitMatrixPacker.h"
#include "Region.h"
#include "FlowBuffer.h"
//#include "ExpBkgWGPFit.h"
#include "DiffEqModel.h"
#include "MultiFlowModel.h"
#include "Utils.h"
#include "BkgMagicDefines.h"
#include "BeadParams.h"
#include "BeadTracker.h"
#include "RegionParams.h"
#include "RegionTracker.h"
#include "EmphasisVector.h"
#include "GlobalDefaultsForBkgModel.h"
#include "CrossTalkSpec.h"
#include "TimeCompression.h"
#include "BkgFitOptim.h"
#include "FitControl.h"


#include "BkgTrace.h"
#include "EmptyTraceTracker.h"
#include "LevMarState.h"
#include "BkgSearchAmplitude.h"
#include "DataCube.h"
#include "XtalkCurry.h"
#include "SynchDat.h"
#include "TraceChunk.h"

#include "DebugWriter.h"
#include "GlobalWriter.h"

#include "RegionalizedData.h"

//#include "BkgDataPointers.h"

//#define FIT_DOUBLE_TAP_DATA

// Initialize lev mar sparse matrices
void InitializeLevMarSparseMatrices (int *my_nuc_block);

// SHOULD be called after all SignalProcessingMasterFitter objects are deleted to free up memory
void CleanupLevMarSparseMatrices (void);


// forward declaration to make function calls happy
class BkgModelCuda;
class MultiFlowLevMar;
class Axion;
class RefineFit;
class SpatialCorrelator;
class RefineTime;
class TraceCorrector;

/*
// declare save_construct_data
class SignalProcessingMasterFitter;
namespace boost {namespace serialization {
    template<class Archive>
      inline void save_construct_data(Archive& ar, const SignalProcessingMasterFitter * o, const unsigned int version);
    template<class Archive>
      inline void load_construct_data(Archive& ar, SignalProcessingMasterFitter * o, const unsigned int version);
  }}
*/

class SignalProcessingMasterFitter
{
  public:

    // Forward declaration of the CUDA model which needs private access to the CPU model
    // too many friends!  Need to clean up the relationships between all these objects
    // i.e. >data< vs >control flow< vs >fitting<
    friend class BkgModelCuda;
    friend class MultiFlowLevMar;
    friend class Axion;
    friend class RefineFit;
    friend class SpatialCorrelator;
    friend class RefineTime;
    friend class TraceCorrector;

    friend class debug_collection;

    // constructor used by Analysis pipeline

    SignalProcessingMasterFitter (RegionalizedData *local_patch, GlobalDefaultsForBkgModel &_global_defaults, char *_results_folder, Mask *_mask, PinnedInFlow *_pinnedInFlow, RawWells *_rawWells, Region *_region, std::set<int>& sample,
              std::vector<float>& sep_t0_est,bool debug_trace_enable,
              int _rows, int _cols, int _frames, int _uncompFrames, int *_timestamps,  EmptyTraceTracker *emptyTraceTracker,
              float sigma_guess=2.5,float t0_guess=35,
	      SequenceItem* seqList=NULL,int numSeqListItems=2,
              bool restart=false);

    // constructor used for testing outside of Analysis pipeline (doesn't require mask, region, or RawWells obects)
    SignalProcessingMasterFitter (GlobalDefaultsForBkgModel &_global_defaults, int numLBeads, int numFrames,
              float sigma_guess=2.5,float t0_guess=35);
    void SetUpFitObjects();

    void SetComputeControlFlags (bool enable_xtalk_correction=true);

    void   UpdateBeadBufferingFromExternalEstimates (std::vector<float> *tauB, std::vector<float> *tauE);

    virtual ~SignalProcessingMasterFitter();

    // image process entry point used by Analysis pipeline
    // trigger computation
    //bool TestAndExecuteBlock (int flow, bool last);
    //void FitUpstreamModel (int flow, bool last);
    void MultiFlowRegionalFitting ( int flow, bool last );
    void FitAllBeadsForInitialFlowBlock();
    void RemainingFitStepsForInitialFlowBlock();
    void FitEmbarassinglyParallelRefineFit();
    void PreWellCorrectionFactors();
    void ExportAllAndReset(int flow, bool last);
    bool TestAndTriggerComputation(bool last);
    //void ExecuteBlock (int flow, bool last);
   
    // break apart image processing and computation
    bool  ProcessImage (Image *img, int flow);
    bool ProcessImage (SynchDat &data, int flow);

    // image process entry point for testing outside of Analysis pipeline
    // (doesn't require Image object, takes well data and background data separately)
    // bool  ProcessImage (short *img, short *bkg, int flow, bool last);

    void SetupTimeAndBuffers (float sigma_guess,
                              float t_mid_nuc_guess,
                              float t0_offset);

    void SetPointers (BkgDataPointers *ptrs)
    {
      global_state.mPtrs = ptrs;
    }

    void SetImageParams (int _rows, int _cols, int _frames, int _uncompFrames, int *_timestamps)
    {
      region_data->my_trace.SetImageParams (_rows,_cols,_frames,_uncompFrames,_timestamps);
    }


    void GetRegParams (struct reg_params *pOut)
    {
      // only for backward compatibility of Obsolete/bkgFit.cpp
      if (pOut != NULL)
      {
        memcpy (pOut,& (region_data->my_regions.rp),sizeof (struct reg_params));
      }
    }

    void GetRegParams (struct reg_params &pOut)
    {
      //ION_ASSERT(sizeof(pOut) == sizeof(reg_params), "sizeof(pOut) != sizeof(reg_params)");
      if (&pOut != NULL)
      {
        // avoid the variable sizeof(reg_params) problem, if there's one here, due to the align(16) statement in the struc
        if (sizeof (pOut) == sizeof (reg_params))
          memcpy (&pOut,& (region_data->my_regions.rp),sizeof (struct reg_params));
        else   // if size diff due to align(16) in reg_params
        {
          std::cout << "Warning in SignalProcessingMasterFitter::GetRegParams(): sizeof(pOut)=" << sizeof (pOut) << " != sizeof(reg_params)=" << sizeof (reg_params) << std::endl;
          pOut = region_data->my_regions.rp;
        }
      }
    }

    Region *GetRegion (void)
    {
      return region_data->region;
    }

    // evaluates the model function for a set of well and region parameters using the background
    // data already stored in the SignalProcessingMasterFitter object from previous calls to ProcessImage
    int GetModelEvaluation (int iWell,struct bead_params *p,struct reg_params *rp,
                            float **fg,float **bg,float **feval,float **isig,float **pf);

    void SetXtalkName (char *my_name)
    {
      char xtalk_name[512];
      strcpy (xtalk_name,my_name);
      xtalk_spec.ReadCrossTalkFromFile (xtalk_name); //happens after initialization
    };


    void DumpExemplarBead (FILE *my_fp, bool debug_only)
    {
      if (region_data->region!=NULL)
        region_data->my_beads.DumpBeads (my_fp,debug_only, region_data->region->col, region_data->region->row);
    }
    void DumpDarkMatterTitle (FILE *my_fp)
    {
      region_data->my_regions.missing_mass.DumpDarkMatterTitle (my_fp);
    }
    void DumpDarkMatter (FILE *my_fp)
    {
      if (region_data->region!=NULL)
        region_data->my_regions.missing_mass.DumpDarkMatter (my_fp,region_data->region->col,region_data->region->row,region_data->my_regions.rp.darkness[0]);
    }

    void DumpExemplarBeadDcOffset (FILE *my_fp, bool debug_only)
    {
      if (region_data->region!=NULL)
        region_data->my_trace.DumpBeadDcOffset (my_fp, debug_only, region_data->my_beads.DEBUG_BEAD, region_data->region->col,region_data->region->row,region_data->my_beads);
    }
    void DumpTimeAndEmphasisByRegion (FILE *my_fp)
    {
      region_data->DumpTimeAndEmphasisByRegion (my_fp);
    };
    void DumpTimeAndEmphasisByRegionH5 (int r)
    {
      global_state.DumpTimeAndEmphasisByRegionH5 (r,region_data->time_c,region_data->emphasis_data);
    };
    void DumpBkgModelBeadFblkInfo (int r);

    
    void SetPoissonCache (PoissonCDFApproxMemo *poiss)
    {
      math_poiss = poiss;
    }

    // functions to access the private members of SignalProcessingMasterFitter
    int get_trace_imgFrames()
    {
      return region_data->my_trace.imgFrames;
    }

    int get_time_c_npts()
    {
      return region_data->time_c.npts();
    }
    float get_t_mid_nuc_start()
    {
      return region_data->t_mid_nuc_start;
    }
    float get_sigma_start()
    {
      return region_data->sigma_start;
    }
    // why are we not accessing the region directly here?
    int get_region_col()
    {
      return region_data->region->col;
    }
    int get_region_row()
    {
      return region_data->region->row;
    }

    GlobalDefaultsForBkgModel getGlobalDefaultsForBkgModel()
    {
      return global_defaults;
    }

// making this public for temporary simplicity
    // Data and parameters here --------------
    RegionalizedData *region_data;


  private:

    bool LoadOneFlow (Image *img, int flow);
    bool LoadOneFlow (SynchDat &data, int flow);

    void AllocateRegionDataIfNeeded(Image *img);
    void AllocateRegionDataIfNeeded(SynchDat &data);

    bool TriggerBlock (bool last);
    bool NeverProcessRegion();
    bool IsFirstBlock(int flow);
    
    void DoPreComputationFiltering();    
    void PostModelProtonCorrection();
    void CompensateAmplitudeForEmptyWellNormalization();

    void NothingInit();

    void  BkgModelInit (bool debug_trace_enable,float sigma_guess,
                        float t0_guess,std::vector<float>& sep_t0_est,std::set<int>& sample,SequenceItem* _seqList,int _numSeqListItems, bool restart);


    void NothingFitters();
    void SetFittersIfNeeded();
    void DestroyFitObjects();

    void InitXtalk();

    // when we're ready to do something to a block of Flows
    // void FitModelForBlockOfFlows (int flow, bool last);

    // export data of various types per flow
    void ExportStatusToMask();
    void ExportDataToWells();
    void ExportDataToDataCubes (bool last);

    void UpdateClonalFilterData (int flow);
    void ResetForNextBlockOfData();

    // emphasis vector stuff
    /*  void    GenerateEmphasis(float amult,float width, float ampl);*/
    void    ReadEmphasisVectorFromFile (void);
    void    SetTimeAndEmphasis (float t0_offset);

    /* Older function set */
    void CPUxEmbarassinglyParallelRefineFit();
    
    /* Relevant functions to integrate GPU multi flow fitting into signal processing pipeline*/
    void RefineAmplitudeEstimates (double &elapsed_time, Timer &fit_timer);
    void ApproximateDarkMatter(bool isSampled);
    void FitAmplitudeAndDarkMatter (double &elapsed_time, Timer &fit_timer);
    void PostKeyFit (double &elapsed_time, Timer &fit_timer);
    void PostKeyFitAllWells(double &elapsed_time, Timer &fit_timer);
    void FitWellParametersConditionalOnRegion (double &elapsed_time, Timer &fit_timer);
    void BootUpModel (double &elapsed_time,Timer &fit_timer);
    void PickRepresentativeHighQualityWells();
    void GuessCrudeAmplitude (double &elapsed_time, Timer &fit_timer);
    void FitTimeVaryingRegion (double &elapsed_time, Timer &fit_timer);
    void RegionalFittingForInitialFlowBlock();
    void RegionalFittingForLaterFlowBlock();

    // debugging functions

    void    DumpRegionParameters();


    // end debugging functions


    // things that fit the model after here ----------------

    // pointer to the actual common defaults entity at top level
    GlobalDefaultsForBkgModel &global_defaults;

    // talking to external world
    extern_links global_state;

    debug_collection my_debug;
    
    // cache math that all bkgmodel objects need to execute
    PoissonCDFApproxMemo *math_poiss;

// local region cross-talk parameters - may vary across the chip by region
    CrossTalkSpecification xtalk_spec;
    XtalkCurry xtalk_execute;

    // corrector for proton
    SpatialCorrelator *correct_spatial;
    RefineTime *refine_time_fit;
    TraceCorrector *trace_bkg_adj;

    // optimizers for generating fits to data
    SearchAmplitude my_search; // crude/quick guess at amplitude
    MultiFlowLevMar *lev_mar_fit; // allocated size of beads

    // this controls refining the fit for a collection of beads
    RefineFit *refine_fit;

    // specialized residual fitter
    Axion *axion_fit;
    
 private:
    SignalProcessingMasterFitter();
    /*
    SignalProcessingMasterFitter(GlobalDefaultsForBkgModel& obj)
      : global_defaults (obj) {
        math_poiss = NULL;
	lev_mar_fit = NULL;
	refine_fit = NULL;
	axion_fit = NULL;
	correct_spatial = NULL;
	refine_time_fit = NULL;
	region_data = NULL;
	trace_bkg_adj = NULL;
    }

    // Serialization section
    friend class boost::serialization::access;

    template<typename Archive>
      friend void boost::serialization::save_construct_data(Archive& ar, const SignalProcessingMasterFitter * o, const unsigned int version);
    template<class Archive>
      friend void boost::serialization::load_construct_data(Archive& ar, SignalProcessingMasterFitter * o, const unsigned int version);

    template<typename Archive>
      void serialize(Archive& ar, const unsigned int version) const {
      fprintf(stdout, "Serialize SignalProcessingMasterFitter\n");
        ar & 
	  const_cast<extern_links &>(global_state) &
	  // my_debug &
	  // math_poiss &
	  const_cast<CrossTalkSpecification &>(xtalk_spec);
	// xtalk_execute &
	fprintf(stdout, "Serialization, SignalProcessingMasterFitter: need to rebuild my_debug (debug_collection), math_poiss (PoisonCDFApproxMemo), xtalk_execute (XtalkCurry)\n");
    }
    */
};

/*
namespace boost { namespace serialization {
    template<typename Archive>
      inline void save_construct_data(Archive& ar, const SignalProcessingMasterFitter * o, const unsigned int version){
      fprintf(stdout, "Serialization: save_construct_data SignalProcessingMasterFitter\n");
      GlobalDefaultsForBkgModel * gb = & o->global_defaults;
      ar << gb;
    }

    template<typename Archive>
      inline void load_construct_data(Archive& ar, SignalProcessingMasterFitter * o, const unsigned int version){
      // will this give us scoping issues???
      GlobalDefaultsForBkgModel *global_defaults_ref;
      ar >> global_defaults_ref;
      ::new(o) SignalProcessingMasterFitter(*global_defaults_ref);
    }
  }
}
*/

#endif // SIGNALPROCESSINGMASTERFITTER_H
