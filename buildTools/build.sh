#!/bin/bash
# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved

# Version of build is now
# determined from "version" file
# in each module folder by the
# cmake build system

BUILD_ROOT=`pwd`

# The only remaining module in this repo is Analysis
MODULE=Analysis
if [ ! -d "$MODULE" ]; then
  echo "Must run $0 from the root folder which includes ./$MODULE:"
  exit -1;
fi

ERR=0
ERRMSG=""

echo "=================================================="
echo " Building module $MODULE"
echo "=================================================="
mkdir -p build/$MODULE
(
  LOCALERR=0
  find build/$MODULE -name \*.deb | xargs rm -f
  cd build/$MODULE
  cmake $@ -G 'Unix Makefiles' $BUILD_ROOT/$MODULE -DION_AVX:BOOL=FALSE
  if [ "$?" != 0 ]; then LOCALERR=1; fi
    make -j13
  if [ "$?" != 0 ]; then LOCALERR=1; fi
    # make test
  if [ "$?" != 0 ]; then LOCALERR=1; fi

  make package

  if [ "$?" != 0 ]; then LOCALERR=1; fi
  find . -name _CPack_Packages | xargs rm -rf
  # do not delete; only used for official builds
  #    if [ -x ../../$MODULE/srcmkr.sh ]; then
  #      ../../$MODULE/srcmkr.sh
  #    fi
  if [ "$LOCALERR" != 0 ]; then
    false
  else
    true
  fi
)
if [ "$?" != 0 ]; then
  ERR=$(($ERR + 1))
  ERRMSG="${ERRMSG}Build of module $MODULE failed.\n"
fi
echo "=================================================="
echo

if [ $ERR != 0 ]; then
  echo -e $ERRMSG
  echo "FAILURES: $ERR modules failed to build."
  exit $ERR
else
  echo "SUCCESS: All modules built."
fi
