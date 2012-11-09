/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef T0CALC_H
#define T0CALC_H

#include <iostream>
#include <set>
#include "GridMesh.h"
#include "FindSlopeChange.h"
#include "T0Model.h"
#include "PJobQueue.h"
#include "Mask.h"
#include "Traces.h"

/** Some ideas about where we think t0 should be for a region. */
class T0Prior {

public:
  T0Prior() {
    mTimeStart = 0;
    mTimeEnd = std::numeric_limits<size_t>::max();
    mT0Prior = 0;
    mT0Weight = 0;
  }

  size_t mTimeStart, mTimeEnd; ///< Range start and end to investigate.
  float mT0Prior;
  float mT0Weight;
};

/** 
    Class to hold information and algorithms for
    calculating time zero (when nuc hits well) for a region. 
*/
class T0Calc {

 public:

  T0Calc() {
    mMask = NULL;
    mRow = mCol = mFrame = mRowStep = mColStep = 0;
    mWindowSize = 6; 
    mMaxFirstHingeSlope = 0;
    mMinFirstHingeSlope = 0;
    mMaxSecondHingeSlope = 0;
    mMinSecondHingeSlope = 0;
    mAvgNth = 1;
  }

  void SetTimeStamps(int *timestamps, int size) {
    assert(size = mFrame);
    mTimeStamps.resize(size);
    copy(&timestamps[0], &timestamps[0] + size, mTimeStamps.begin());
  }

  /* Don't accumulate every well for average just every Nth. */
  void SetStepSize(int size) { mAvgNth = size; }
  void SetWindowSize(int size) { mWindowSize = size; }
  void SetMaxFirstHingeSlope(float maxSlope) { mMaxFirstHingeSlope = maxSlope; }
  void SetMinFirstHingeSlope(float minSlope) { mMinFirstHingeSlope = minSlope; }
  void SetMaxSecondHingeSlope(float maxSlope) { mMaxSecondHingeSlope = maxSlope; }
  void SetMinSecondHingeSlope(float minSlope) { mMinSecondHingeSlope = minSlope; }

  void SetMask(Mask *mask) { mMask = mask; }

  /** Setup our t0 calculator. */
  void Init(size_t nRow, size_t nCol, size_t nFrame,
	    size_t nRowStep, size_t nColStep, size_t nThreads) {
    mRow = nRow;
    mCol = nCol;
    mFrame = nFrame;
    mRowStep = nRowStep;
    mColStep = nColStep;
    /* Init our grids. */
    mRegionSum.Init(mRow, mCol, mRowStep, mColStep);
    mT0.Init(mRow, mCol, mRowStep, mColStep);
    mSlope.Init(mRow, mCol, mRowStep, mColStep);
    mT0Prior.Init(mRow, mCol, mRowStep, mColStep);

    /* Set some default values. */
    size_t numBin = mRegionSum.GetNumBin();
    for (size_t bIx = 0; bIx < numBin; bIx++) {
      float &t0 = mT0.GetItem(bIx);
      t0 = 0;
      std::pair<size_t, std::vector<float> > &bin = mRegionSum.GetItem(bIx);
      bin.first = 0;
      bin.second.resize(mFrame);
      std::fill(bin.second.begin(), bin.second.end(), 0.0);
    }
  }

  bool isOk(size_t wIx) {
    if(mMask == NULL || !(mMask->Match(wIx, MaskPinned) || mMask->Match(wIx, MaskExclude) || mMask->Match(wIx, MaskIgnore))) {
      return true;
    }
    return false;
  }

  /** Algorithm to fit t0 for this trace using a two piece linear model. */
  static void CalcT0(T0Finder &finder, std::vector<float> &trace, 
                     std::vector<float> &timestamps,
                     T0Prior &prior, float &t0) {
    /* Convert the t0 prior in time into frames. */
    std::vector<float>::iterator tstart;
    int frameEnd = trace.size();
    tstart = std::lower_bound(timestamps.begin(), timestamps.end(), prior.mTimeEnd);
    int frame = tstart - timestamps.begin();
    if (tstart != timestamps.end() &&  frame < frameEnd) {
      frameEnd = frame;
    }
    int frameStart = 0;
    tstart = std::lower_bound(timestamps.begin(), timestamps.end(), prior.mTimeStart);
    frame = tstart - timestamps.begin();
    if (tstart != timestamps.end() && (frame) -1 > frameStart) {
      frameStart = frame - 1;
    }
    finder.SetSearchRange(frameStart, frameEnd);
    bool ok = finder.FindT0Time(&trace[0], &timestamps[0], trace.size());
    if (ok) {
      t0 = finder.GetT0Est();
      t0 = (prior.mT0Weight * prior.mT0Prior) + t0;
      t0 = t0 / (prior.mT0Weight + 1);
    }
    else {
      t0 = -1;
    }
  }

  void CalcIndividualT0(std::vector<float> &t0, int useMeshNeighbors) {
    std::vector<double> dist(7);
    std::vector<float *> values;

    t0.resize(mRow*mCol);
    fill(t0.begin(), t0.end(), -1);
    for (size_t rowIx = 0; rowIx < mRow; rowIx++) {
      for (size_t colIx = 0; colIx < mCol; colIx++) {
        int idx = rowIx * mCol + colIx;
        mT0.GetClosestNeighbors(rowIx, colIx, useMeshNeighbors, dist, values);
        double distWeight = 0;
        double startX = 0;
        for (size_t i = 0; i < values.size(); i++) {
          if (*(values[i]) > 0) {
            double w = Traces::WeightDist(dist[i]); //1/sqrt(dist[i] + 1);
            distWeight += w;
            startX += w * (*(values[i]));
          }
        }
        if (distWeight > 0 && startX >= 0) {
          t0[idx]  = startX / distWeight;
        }
      }
    }
  }

  /** Calculate the t0 from accumulated traces. */
  void CalcT0FromSum() {
    size_t numBin = mRegionSum.GetNumBin();
    T0Finder mFinder;
    mFinder.SetWindowSize(mWindowSize);
    mFinder.SetFirstSlopeRange(mMinFirstHingeSlope, mMaxFirstHingeSlope);
    mFinder.SetSecondSlopeRange(mMinSecondHingeSlope, mMaxSecondHingeSlope);
    for (size_t bIx = 0; bIx < numBin; bIx++) {
      std::pair<size_t, std::vector<float> > &bin = mRegionSum.GetItem(bIx);
      if (bin.first < 20) {
        continue;
      }
      for (size_t fIx = 0; fIx < bin.second.size(); fIx++) {
        bin.second[fIx] = bin.second[fIx] / bin.first;
      }
      T0Prior &prior = mT0Prior.GetItem(bIx);
      float &t0 = mT0.GetItem(bIx);
      T0Calc::CalcT0(mFinder, bin.second, mTimeStamps, prior, t0);
    }
  }
  
  /** Sum up the traces for the regions. */
  template <typename shortvec> void CalcAllSumTrace(shortvec &data) {
    size_t numBin = mRegionSum.GetNumBin();
    int rowStart, rowEnd, colStart, colEnd;
    for (size_t bIx = 0; bIx < numBin; bIx++) {
      mRegionSum.GetBinCoords(bIx, rowStart, rowEnd, colStart, colEnd);
      std::pair<size_t, std::vector<float> > &bin = mRegionSum.GetItem(bIx);
      CalcSumTrace(data, rowStart, rowEnd, colStart, colEnd, bin);
    }
  }

  /** Sum up the traces for this region. */
  template <typename shortvec> void CalcSumTrace(shortvec &data, 
		    int rowStart, int rowEnd,
		    int colStart, int colEnd,
                    std::pair<size_t, std::vector<float> > &bin) {
    size_t frameStep = mRow * mCol;
    int notOk = 0;
    for (int row = rowStart; row < rowEnd; row++) {
      for (int col = colStart; col < colEnd; col++) {
        int wIx = row * mCol + col;
        if (isOk(wIx)) {
          bin.first++;
          for (size_t frameIx = 0; frameIx < mFrame; frameIx++) {
            bin.second[frameIx] += (float)(data[frameStep * frameIx + wIx]) - data[wIx];
          }
        }
        else {
          notOk++;
        }
      }
    }
    /* if (notOk > 0 || bin.first > 0) { */
    /*   int ok = bin.first; */
    /*   cout << "Ok: " << bin.first << "not ok: " << notOk << endl; */
    /*   if (ok == notOk) { */
    /*     cout << "weird." << endl; */
    /*   } */
      
    /* } */
    /* if(bin.first > 0) { */
    /*   cout << rowStart << "\t" << rowEnd << "\t" << colStart << "\t" << colEnd; */
    /*   for (size_t f =0; f < mFrame; f++) { */
    /*     cout << "\t" << bin.second[f]; */
    /*   } */
    /*   cout << endl; */
    /*   bin.first = -1 * bin.first; */
    /*   bin.first = -1 * bin.first; */
    /* } */
  }
/*     int frameIx = 0; */
/*     /\*  */
/*        Loop through first frame to count wells in each region. This */
/*        loop is optimized to loop through the memory in order of the */
/*        frames in memory while minimizing the number of lookups into */
/*        the regionSum grid (once every mColStep) */
/*     *\/ */
/*     for (int rowIx = rowStart; rowIx < rowEnd; rowIx++) { */
/*       for (int colIx = colStart; colIx < (int)colEnd; colIx+= mColStep) { */
/*         int wIx = rowIx * mCol + colIx; */
/* 	int end = min(mCol, mColStep + colIx);  */
/* 	std::pair<size_t, std::vector<float> > &bin = mRegionSum.GetItemByRowCol(rowIx, colIx); */
/* 	for (int cIx = colIx; cIx < end; cIx++) { */
/* k	  if (isOk(wIx)) { */
/* 	    bin.first++; */
/* 	    bin.second[frameIx] += (float)data[wIx] - data[wIx]; */
/*           } */
/* 	  wIx++; */
/* 	} */
/*       } */
/*     } */

/*     /\* Loop through rest of the frames. *\/ */
/*     for (int frameIx = 1; frameIx < (int)mFrame; frameIx++) { */
/*       for (int rowIx = rowStart; rowIx < rowEnd; rowIx++) { */
/*         for (int colIx = colStart; colIx < (int)colEnd; colIx+= mColStep) { */
/*           int fwIx = rowIx * mCol + colIx; */
/*           if (fwIx % mAvgNth == 0) { */
/*             int wIx = fwIx + frameIx * mRow * mCol; */
/*             int end = min(mCol, mColStep + colIx);  */
/*             std::pair<size_t, std::vector<float> > &bin = mRegionSum.GetItemByRowCol(rowIx, colIx); */
/*             for (int cIx = colIx; cIx < end; cIx++) { */
/*               if (isOk(fwIx)) { */
/*                 bin.second[frameIx] += (float) data[wIx] - data[fwIx]; */
/*               } */
/*               fwIx++; */
/*               wIx++; */
/*             } */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*  } */

  /** Write out our fits and average trace as a text file if requested. */
  void WriteResults(std::ostream &out) {
    int rowStart, rowEnd, colStart, colEnd;
    for (int bIx = 0; bIx < (int)mRegionSum.GetNumBin(); bIx++) {
      float t0 = mT0.GetItem(bIx);
      float t0Frame = GetInterpolatedFrameForTime(t0);
      mRegionSum.GetBinCoords(bIx, rowStart, rowEnd, colStart, colEnd);
      std::pair<size_t, std::vector<float> > &vec = mRegionSum.GetItem(bIx);
      std::vector<float> &v = vec.second;
      int first = t0Frame + 1;
      int second = t0Frame + 10;
      first = min(first, (int) v.size());
      second = min(second, (int) v.size());
      float slope = GetSlope(bIx);// (v[second] - v[first])/(mTimeStamps[second] - mTimeStamps[first]);
      float fitSlope = 0; //mSlope.GetItem(bIx);
      out << rowStart << '\t' << rowEnd << '\t' << colStart << '\t' << colEnd << '\t' << vec.first << '\t' << t0 << '\t' << t0Frame << '\t' << fitSlope << '\t' << slope;
      for (size_t i = 0; i < vec.second.size(); i++) {
        out << '\t' << vec.second.at(i);
      }
      out << std::endl;
    }
  }
  
  float GetT0(int bin) { return mT0.GetItem(bin); }

  void SetT0(int bin, float val) { mT0.GetItem(bin) = val; }

  float GetSlope(int bin) { return mSlope.GetItem(bin); }

  float GetT0(int row, int col) { return mT0.GetItemByRowCol(row, col); }

  void CalculateSlopePostT0(int neighbors) {
    for (size_t binIx = 0; binIx < mSlope.GetNumBin(); binIx++) {
      float &s = mSlope.GetItem(binIx);
      s = GetSlopePostT0(binIx);
    }
    GridMesh<float> smoothed = mSlope;
    std::vector<double> dist(7);
    std::vector<float *> values;
    int rowStart, rowEnd, colStart, colEnd;
    for (size_t binIx = 0; binIx < mSlope.GetNumBin(); binIx++) {
      float &smooth = smoothed.GetItem(binIx);
      float nosmooth = mSlope.GetItem(binIx);
      mSlope.GetBinCoords(binIx, rowStart, rowEnd, colStart, colEnd);
      int rowIx = (rowEnd + rowStart)/2;
      int colIx = (colEnd + colStart)/2;
      mSlope.GetClosestNeighbors(rowIx, colIx, neighbors, dist, values);
      double distWeight = 0;
      double startX = 0;
      for (size_t i = 0; i < values.size(); i++) {
        if (*(values[i]) > 0) {
          double w = Traces::WeightDist(dist[i]); //1/sqrt(dist[i] + 1);
          distWeight += w;
          startX += w * (*(values[i]));
        }
      }
      if (distWeight > 0 && startX >= 0) {
        smooth  = startX / distWeight;
      }
      else {
        smooth = nosmooth;
      }
    }
    mSlope = smoothed;
  }

  float GetSlopePostT0(size_t binIx) {
    std::pair<size_t, std::vector<float> > &vec = mRegionSum.GetItem(binIx);
    std::vector<float> &v = vec.second;
    float t0 = GetT0(binIx);
    int t0Frame = GetFrameForTime(t0);
    int first = t0Frame + 1;
    int second = t0Frame + 10;
    first = min(first, (int) v.size());
    second = min(second, (int) v.size());
    float slope = (v[second] - v[first])/(mTimeStamps[second] - mTimeStamps[first]);
    return slope;
  }

  int GetFrameForTime(float time) {
    if (time <= 0) {
      return -1;
    }
    // @todo - use lower_bound()
    for (size_t i = 0; i < mTimeStamps.size() - 1; i++) {
      if (time < mTimeStamps[i+1]) {
        return i;
      }
    }
    assert(0);
    return -1;
  }

  float GetInterpolatedFrameForTime(float time) {
    if (time <= 0) {
      return -1.0f;
    }
    // @todo - use lower_bound()
    for (size_t i = 0; i < mTimeStamps.size() - 1; i++) {
      if (time < mTimeStamps[i+1]) {
        return i + ((time - mTimeStamps[i])/(mTimeStamps[i+1]-mTimeStamps[i]));
      }
    }
    assert(0);
    return -1.0f;
  }

  size_t GetNumRow() { return mRow; }
  size_t GetNumCol() { return mCol; }

  size_t GetNumRegions() { return mT0.GetNumBin(); }
  
  void GetRegionCoords(int regionIdx,int &rowStart, int &rowEnd, int &colStart, int &colEnd) {
    return mT0.GetBinCoords(regionIdx, rowStart, rowEnd, colStart, colEnd);
  }

  void FillInT0Prior(GridMesh<T0Prior> &priorMesh, int preMillis, int postMillis, float weight) {
    priorMesh.Init(mRow, mCol, mRowStep, mColStep);
    assert(priorMesh.GetNumBin() == mT0.GetNumBin());
    for (size_t bIx = 0; bIx < mT0.GetNumBin(); bIx++) {
      T0Prior &prior = priorMesh.GetItem(bIx);
      float t0 = GetT0(bIx);
      prior.mT0Weight = 0.0f;
      prior.mTimeStart = 0;
      prior.mTimeEnd = std::numeric_limits<size_t>::max();
      prior.mT0Prior = -1.0f;
      if (t0 > 0) {
        prior.mTimeStart = max(t0 - preMillis, (float) prior.mTimeStart);
        prior.mTimeEnd = min(t0 + preMillis, (float) prior.mTimeEnd);
        prior.mT0Prior = t0;
        prior.mT0Weight = weight;
      }
    }
  }

  void SetT0Prior(GridMesh<T0Prior> &prior) { mT0Prior = prior; }

  template <typename ShortVec> void ComputeT0(ShortVec& img)
  {
    size_t numBin = mT0.GetNumBin();

    int rowStart, rowEnd, colStart, colEnd;
    size_t frameStep = mRow * mCol;
    //std::vector<bool> badWells(frameStep, false); //frameStep is also the number of wells

    for (size_t bin = 0; bin < numBin; ++bin)
    {
        size_t numWells;
        std::vector<float> averageTrace(mFrame);
        std::set<size_t> badWells;
        //int numBadWells = 0;
        long firstFrameTotal = 0;

        mT0.GetBinCoords(bin, rowStart, rowEnd, colStart, colEnd);
        int regionSize = ((rowEnd-rowStart)*(colEnd-colStart));

        //first loop is a special case (have to normalize averageTrace vector and tally good wells
        for (int col = colStart; col < colEnd; ++col)
        {
            for (int row = rowStart; row < rowEnd; ++row)
            {
                size_t well = (row * mCol) + col;

                if (!isOk(well))
                {
                    badWells.insert(well); //if the well is bad, remember to ignore it in the future
                    //badWells[well] = true;
                    //++numBadWells;
                }
                else
                {
                    firstFrameTotal += img[well];
                }
            }
        }
        numWells = regionSize - badWells.size(); //the total number of good wells is the number in the region minus the bad ones
        //numWells = regionSize - numBadWells;
        averageTrace[0] = 0; //subtracting the first frame from itself always results in zero

        //std::cout << "Bin " << bin << " has " << numWells << " live wells of " << regionSize << " total: " << firstFrameTotal << std::endl;
        if (numWells < 10) continue; //dead region...let's not waste our time


        //compute an average trace from all the wells
        std::set<size_t>::iterator badWellsEnd = badWells.end();
        //std::cout << "Trace: ";
        for (size_t frame = 1; frame < mFrame; ++frame)
        {
            long frameTotal = 0;
            for (int col = colStart; col < colEnd; ++col)
            {
                for (int row = rowStart; row < rowEnd; ++row)
                {
                    size_t well = (row* mCol) + col;
                    if (badWells.find(well) == badWellsEnd)
                    //if (!badWells[well])
                    {
                        frameTotal += img[(frameStep * frame) + well];
                    }
                }
            }

            averageTrace[frame] = (float)(frameTotal - firstFrameTotal) / (float)numWells;
            //std::cout << averageTrace[frame] << ":" << frameTotal << " ";
        }
        //std::cout << std::endl;

        //use the average trace to get a T0 and slope estimate
        T0Prior& prior = mT0Prior.GetItem(bin);
        float& t0 = mT0.GetItem(bin);
        T0Finder finder;
        finder.SetWindowSize(mWindowSize);
        finder.SetFirstSlopeRange(mMinFirstHingeSlope, mMaxFirstHingeSlope);
        finder.SetSecondSlopeRange(mMinSecondHingeSlope, mMaxSecondHingeSlope);
        T0Calc::CalcT0(finder, averageTrace, mTimeStamps, prior, t0);
    }
  }

  //compute slopes to use for smoothing step
  template <typename ShortVec> void ComputeSlopes(ShortVec& img)
  {
    size_t numBin = mSlope.GetNumBin();

    int rowStart, rowEnd, colStart, colEnd;
    size_t frameStep = mRow * mCol;
    //std::vector<bool> badWells(frameStep, false); //frameStep is also the number of wells

    for (size_t bin = 0; bin < numBin; ++bin)
    {
        size_t numWells;
        std::vector<float> averageTrace(mFrame);
        std::set<size_t> badWells;
        //size_t numBadWells = 0;
        long firstFrameTotal = 0;

        mSlope.GetBinCoords(bin, rowStart, rowEnd, colStart, colEnd);
        int regionSize = ((rowEnd-rowStart)*(colEnd-colStart));

        //first loop is a special case (have to normalize averageTrace vector and tally good wells
        for (int col = colStart; col < colEnd; ++col)
        {
            for (int row = rowStart; row < rowEnd; ++row)
            {
                size_t well = (row * mCol) + col;

                if (!isOk(well))
                {
                    badWells.insert(well); //if the well is bad, remember to ignore it in the future
                    //badWells[well] = true;
                    //++numBadWells;
                }
                else
                {
                    firstFrameTotal += img[well];
                }
            }
        }
        numWells = regionSize - badWells.size(); //the total number of good wells is the number in the region minus the bad ones
        //numWells = regionSize - numBadWells;
        averageTrace[0] = 0; //subtracting the first frame from itself always results in zero

        if (numWells < 10) continue; //dead region...let's not waste our time


        //compute an average trace from all the wells
        std::set<size_t>::iterator badWellsEnd = badWells.end();
        for (size_t frame = 1; frame < mFrame; ++frame)
        {
            long frameTotal = 0;
            for (int col = colStart; col < colEnd; ++col)
            {
                for (int row = rowStart; row < rowEnd; ++row)
                {
                    size_t well = (row * mCol) + col;
                    if (badWells.find(well) == badWellsEnd)
                    //if (!badWells[well])
                    {
                        frameTotal += img[(frameStep * frame) + well];
                    }
                }
            }

            averageTrace[frame] = (float)(frameTotal - firstFrameTotal) / (float)numWells;
        }

        //use the average trace to get a T0 and slope estimate
        T0Prior& prior = mT0Prior.GetItem(bin);
        float& t0 = mT0.GetItem(bin);
        T0Finder finder;
        finder.SetWindowSize(mWindowSize);
        finder.SetFirstSlopeRange(mMinFirstHingeSlope, mMaxFirstHingeSlope);
        finder.SetSecondSlopeRange(mMinSecondHingeSlope, mMaxSecondHingeSlope);
        T0Calc::CalcT0(finder, averageTrace, mTimeStamps, prior, t0);

        //GetSlopePostT0
        size_t t0Frame = GetFrameForTime(t0);
        size_t first = t0Frame+1;
        size_t second = t0Frame + 10;
        first = min(first, mFrame);
        second = min(second, mFrame);
        mSlope.GetItem(bin) = (averageTrace[second] - averageTrace[first])/(mTimeStamps[second] - mTimeStamps[first]);
    }
  }

  void EstimateSigma(int neighbors, GridMesh<SigmaEst>& sigmaTMid, SigmaTMidNucEstimation& sigmaEst)
  {
      std::vector<double> dist(7);
      std::vector<float*> values;
      int badCount = 0;
      for (size_t bin = 0; bin < mSlope.GetNumBin(); ++bin)
      {
          int rowStart, rowEnd, colStart, colEnd;
          float smooth = mSlope.GetItem(bin);

          mSlope.GetBinCoords(bin, rowStart, rowEnd, colStart, colEnd);

          int midRow = (rowStart + rowEnd) / 2;
          int midCol = (colStart + colEnd) / 2;

          mSlope.GetClosestNeighbors(midRow, midCol, neighbors, dist, values);

          double distWeight = 0;
          double startX = 0;
          for (size_t i = 0; i < values.size(); ++i)
          {
              double w = Traces::WeightDist(dist[i]);
              distWeight += w;
              startX += w * (*(values[i]));
          }

          if (distWeight > 0 && startX >= 0)
          {
              smooth = startX / distWeight;
          }

          //armed with a smoothed slope, we can now estimate sigma
          float t0Est = GetT0(bin);
          SigmaEst& est = sigmaTMid.GetItem(bin);
          est.mSigma = 0;
          est.mTMidNuc = 0;
          if (t0Est > 0 && smooth > 0)
          {
              sigmaEst.Predict(t0Est, smooth, est.mSigma, est.mTMidNuc);
          }
          else
          {
	    badCount++;
	    //              std::cout << "Couldn't estimate sigma for bin: " << bin << std::endl;
          }
      }
      std::cout << "Couldn't estimate sigma for: " << badCount << " of " << mSlope.GetNumBin() << " bins." << std::endl;
  }

 private:
  /** Pointer to the mask characterizing wells. */
  Mask *mMask;
  /** Dimensions of chip. */ 
  size_t mRow, mCol, mFrame;
  /** Timestamps for each frame. */
  std::vector<float> mTimeStamps; 
  /** Region height and width */
  size_t mRowStep, mColStep;
  /** Actual t0 as calculated. */
  GridMesh<float> mT0;
  GridMesh<float> mSlope;
  /** Prior for t0 */
  GridMesh<T0Prior>  mT0Prior;
  /** Pair of number of wells seen and vector average of the region. */
  GridMesh<std::pair<size_t, std::vector<float> > > mRegionSum;
  /** Window size for edges of the hinges. */
  int mWindowSize; 
  /** What step size for wells to average. */
  int mAvgNth;
  /* Constraints for slopes on first and second portion of hinge. */
  float mMaxFirstHingeSlope;
  float mMinFirstHingeSlope;
  float mMaxSecondHingeSlope;
  float mMinSecondHingeSlope;

};


#endif // T0CALC_H

