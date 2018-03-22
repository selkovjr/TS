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
}

// ------------------------------------------------------------------------------------

void InputStructures::Initialize(ExtendParameters &parameters, const ReferenceReader& ref_reader, const SamHeader &bam_header)
{
  DEBUG                 = parameters.program_flow.DEBUG;
}



// // =====================================================================================
// // We create one basecaller object per unique flow order
// // This prevents us from having to rebuild and initialize the whole object for every read
// PersistingThreadObjects::PersistingThreadObjects(const InputStructures &global_context)
//     : use_SSE_basecaller(global_context.use_SSE_basecaller), realigner(50, 1)
PersistingThreadObjects::PersistingThreadObjects(const InputStructures &global_context) : realigner(50, 1) {};


