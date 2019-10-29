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
  DEBUG = parameters.program_flow.DEBUG;
}

