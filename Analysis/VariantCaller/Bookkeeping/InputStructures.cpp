/* Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     InputStructures.cpp
//! @ingroup  VariantCaller
//! @brief    HP Indel detection

#include "InputStructures.h"
#include "ExtendedReadInfo.h"
#include "json/json.h"

InputStructures::InputStructures()
{
  DEBUG = 0;
#ifdef __SSE3__
  use_SSE_basecaller  = true;
#else
  use_SSE_basecaller  = false;
#endif
  resolve_clipped_bases = false;
}

// ------------------------------------------------------------------------------------

void InputStructures::Initialize(ExtendParameters &parameters, const ReferenceReader& ref_reader, const SamHeader &bam_header)
{
  DEBUG                 = parameters.program_flow.DEBUG;

  use_SSE_basecaller    = parameters.program_flow.use_SSE_basecaller;
  resolve_clipped_bases = parameters.program_flow.resolve_clipped_bases;

  if (parameters.sseMotifsProvided) {
    cout << "Loading systematic error contexts." << endl;
    read_error_motifs(parameters.sseMotifsFileName);
    cout << "Loaded." << endl;
  }
}



// =====================================================================================
// We create one basecaller object per unique flow order
// This prevents us from having to rebuild and initialize the whole object for every read

PersistingThreadObjects::PersistingThreadObjects(const InputStructures &global_context)
    : use_SSE_basecaller(global_context.use_SSE_basecaller), realigner(50, 1)
{
#ifdef __SSE3__
    if (use_SSE_basecaller) {
	  for (unsigned int iFO=0; iFO < global_context.flow_order_vector.size(); iFO++){
        TreephaserSSE*    treephaser_sse = new TreephaserSSE(global_context.flow_order_vector.at(iFO), DPTreephaser::kWindowSizeDefault_);
        treephaserSSE_vector.push_back(treephaser_sse);
      }
    }
    else {
#endif
      for (unsigned int iFO=0; iFO < global_context.flow_order_vector.size(); iFO++){
        DPTreephaser      dpTreephaser(global_context.flow_order_vector.at(iFO));
        dpTreephaser_vector.push_back(dpTreephaser);
      }
#ifdef __SSE3__
    }
#endif
};

// ------------------------------------------------------------------------------------



