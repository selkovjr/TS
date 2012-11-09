/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include "ImageSpecClass.h"
#include "SynchDatSerialize.h"
#include "H5File.h"

ImageSpecClass::ImageSpecClass()
: rows(0)
, cols(0)
, scale_of_chip(0)
, uncompFrames(0)
, timestamps(0)
, vfr_enabled(true)
, n_timestamps(0)
{
}

int ImageSpecClass::LeadTimeForChipSize()
{
  int totalMem = totalMemOnTorrentServer();
  if ( totalMem > ( 24*1024*1024 ) )
  {
    if ( ( rows * cols ) > 10000000 )
      return ( 8 ); // 318
    if ( ( rows * cols ) > 2000000 )
      return ( 10 ); // 316
    else
      return ( 20 );
  }
  else
  {
    if ( ( rows * cols ) > 10000000 )
      return ( 2 ); // 318
    if ( ( rows * cols ) > 2000000 )
      return ( 4 ); // 316
    else
      return ( 20 );
  }
}

void ImageSpecClass::TimeStampsFromImage ( Image &img, ImageControlOpts &img_control )
{
  vfr_enabled= img.VFREnabled();

  // this globally limits the Loadraw method and other methods to process only this many frames (for speed)
  img_control.totalFrames = img.GetMaxFrames();
  if ( img_control.maxFrames != 0 )
  {
    // command-line override of the frames to analyze
    img_control.maxFrames = ( img_control.maxFrames > img_control.totalFrames ? img_control.totalFrames
                              : img_control.maxFrames );
    img.SetMaxFrames ( img_control.maxFrames );
  }
  else
  {
    img_control.maxFrames = img_control.totalFrames; // set to total frames in image.
  }
  uncompFrames = img.GetUnCompFrames();
  timestamps = new int[img_control.maxFrames];
  memcpy ( timestamps,img.GetImage()->timestamps,sizeof ( int ) *img_control.maxFrames );
  n_timestamps = img_control.maxFrames;
  // Deallocate image memory
}

void ImageSpecClass::DimensionsFromImage(Image &img, SpatialContext &loc_context)
{
  // grab rows & cols here - as good a place as any
  // @TODO: this is no longer obvious(!)
  rows = img.GetImage()->rows;
  cols = img.GetImage()->cols;
  scale_of_chip = rows*cols;

  loc_context.chip_offset_x = img.GetImage()->chip_offset_x;
  loc_context.chip_offset_y = img.GetImage()->chip_offset_y;
  loc_context.rows = rows;
  loc_context.cols = cols;
}

/********************************************************************
             // Open an image file to get some parameters for the dataset
             // use the first nuke flow file instead of beadfind file.
 ********************************************************************/
//@TODO: fix side effects on "command line options"

void ImageSpecClass::DeriveSpecsFromDat ( SystemContext &sys_context, ImageControlOpts &img_control, SpatialContext &loc_context )
{
  Image img;

  ReadFirstImage(img, sys_context, img_control, loc_context);

  TimeStampsFromImage(img,img_control);
  DimensionsFromImage(img,loc_context);

  img.Close();
}

void ImageSpecClass::ReadFirstImage(Image &img, SystemContext &sys_context, ImageControlOpts &img_control, SpatialContext &loc_context )
{
  img.SetImgLoadImmediate ( false );
  img.SetNumAcqFiles ( 1 );
  img.SetIgnoreChecksumErrors ( img_control.ignoreChecksumErrors );
  bool datIsSdat = false;
  char *firstAcqFile = ( char * ) malloc ( strlen ( sys_context.dat_source_directory ) + strlen (
                         img_control.acqPrefix ) + 10 );
  sprintf ( firstAcqFile, "%s/%s%04d.dat", sys_context.dat_source_directory, img_control.acqPrefix, 0 );
  /* Check to see if our dat file is really an sdat file in disguise, if so updat our suffix and doing sdats. */
  if (H5File::IsH5File(firstAcqFile)) {
    datIsSdat = true;
    img_control.sdatSuffix = "dat";
    img_control.doSdat = true;
  }
  char *firstSdatAcqFile = ( char * ) malloc ( strlen ( sys_context.dat_source_directory ) + strlen (
                         img_control.acqPrefix ) + 11 );
  sprintf ( firstSdatAcqFile, "%s/%s%04d.%s", sys_context.dat_source_directory, img_control.acqPrefix, 0, img_control.sdatSuffix.c_str() );

  if (isFile(firstAcqFile) && !datIsSdat) {
    if ( !img.LoadRaw ( firstAcqFile, 0, true, false ) )
      {
        exit ( EXIT_FAILURE );
      }
  }
  else if(isFile(firstSdatAcqFile)) {
    TraceChunkSerializer reader;
    SynchDat sdat;
    reader.Read(firstSdatAcqFile, sdat);
    img.InitFromSdat(&sdat);
    if (loc_context.regionXSize == 0){ // not explicitly set, set now
      loc_context.regionXSize = sdat.GetColStep();
      loc_context.regionYSize = sdat.GetRowStep();
    }
    if ((sdat.GetColStep() != (size_t)loc_context.regionXSize) ||
	(sdat.GetRowStep() != (size_t)loc_context.regionYSize)) {
      fprintf(stdout, "Warning: Analysis region sizes width=%d, height=%d do not match sdat region sizes width=%d, height=%d\n", loc_context.regionXSize, loc_context.regionYSize, (int)sdat.GetColStep(), (int)sdat.GetRowStep());
    }
  }
  else {
    exit ( EXIT_FAILURE );
  }
  img.SetOffsetFromChipOrigin ( firstAcqFile );
  free ( firstAcqFile );
  free ( firstSdatAcqFile );

  if ( !loc_context.IsSetRegionXYSize() ){ // not yet set, has to be set now
    loc_context.SetRegionXYSize(50, 50);
  }

  img.SetDir ( sys_context.results_folder );
  img.SetFlowOffset ( img_control.flowTimeOffset );
}

ImageSpecClass::~ImageSpecClass()
{
  if ( timestamps!=NULL ) delete[] timestamps;
}

