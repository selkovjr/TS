# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved
cmake_minimum_required (VERSION 2.6)
project (ion-analysis)
set(PROJECT_DESCRIPTION "Ion Torrent Data Processing Pipeline")
enable_testing()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckCCompilerFlag)
check_c_compiler_flag( -malign-double HAS_ALIGN_DOUBLE )

set(CMAKE_CXX_FLAGS -g3)
add_compile_options(-fprofile-arcs -ftest-coverage -Wfatal-errors)
set(CMAKE_CXX_STANDARD_LIBRARIES -lgcov)

include(../buildTools/cmake/CMakeLists.version.txt)
include(../buildTools/cmake/CMakeLists.compiler.txt)

FIND_PACKAGE(PythonInterp)
# Determine python site-packages or dist-packages directory
execute_process (
        COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process (
        COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_lib; print get_python_lib(prefix='/usr/local')"
        OUTPUT_VARIABLE PYTHON_LOCAL_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process (
        COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_version; print get_python_version()"
        OUTPUT_VARIABLE PYTHON_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
)


# Where to install - override with: -DION_INSTALL_PREFIX
set(ION_INSTALL_PREFIX "/opt/ion" CACHE PATH "Ion Install Prefix")

set(CPACK_PACKAGING_INSTALL_PREFIX ${ION_INSTALL_PREFIX})
set(CPACK_SET_DESTDIR "ON")

configure_file (
    "${PROJECT_SOURCE_DIR}/debian/postinst.in"
    "${PROJECT_BINARY_DIR}/debian/postinst" @ONLY
)

include(../buildTools/cmake/CMakeLists.dependencies.txt)

add_custom_target(compileall ALL
    COMMAND ${PYTHON_EXECUTABLE} -m compileall -q -f "${PROJECT_SOURCE_DIR}/pynvml/"
)

include_directories("${PROJECT_SOURCE_DIR}/Util")

# Ion Analysis Library
add_library(ion-analysis

    ${PROJECT_BINARY_DIR}/IonVersion.cpp

    Util/IonErr.cpp
    Util/Utils.cpp
    Util/HandleExpLog.cpp
    Util/OptArgs.cpp
    Util/OptionArgsBase.cpp
    Util/OptionArgs.cpp
    Util/LoggerBase.cpp
    Util/Logger.cpp
    Util/WorkerInfoQueue.cpp
    Util/RingBuffer.cpp
    Util/SeqUtils.cpp

    Region.cpp
    LinuxCompat.cpp

    SamUtils/types/Cigar.cpp
    SamUtils/types/Qual.cpp
    SamUtils/types/Sequence.cpp
    SamUtils/types/BAMRead.cpp
    SamUtils/types/MD.cpp
    SamUtils/BAMReader.cpp
    SamUtils/BAMUtils.cpp
    SamUtils/alignStats.cpp

    ${ION_TS_EXTERNAL}/jsoncpp-src-amalgated0.6.0-rc1/jsoncpp.cpp
)
add_dependencies(ion-analysis IONVERSION samtools bamtools)

target_link_libraries(ion-analysis ${EXTRA_LIBS} z)

set_target_properties(ion-analysis PROPERTIES
    VERSION   "${ION_VERSION_MAJOR}.${ION_VERSION_MINOR}.${ION_VERSION_RELEASE}.${ION_VERSION_BUILDNUM}"
)
install(TARGETS ion-analysis DESTINATION lib)


## BAM walker derived from Torrent Variant Caller, tvc
set(ION_VCFLIB_DIR    ${ION_TS_EXTERNAL}/vcflib)
set(ION_FREEBAYES_DIR ${ION_TS_EXTERNAL}/freebayes)

include_directories(${ION_VCFLIB_DIR})
include_directories(${ION_FREEBAYES_DIR}/src)
include_directories("${PROJECT_SOURCE_DIR}/VariantCaller")
include_directories("${PROJECT_SOURCE_DIR}/VariantCaller/Bookkeeping")
include_directories("${PROJECT_SOURCE_DIR}/VariantCaller/Reads")
include_directories("${PROJECT_SOURCE_DIR}/realignment")

set(tvcSRCS
  VariantCaller/VariantCaller.cpp
  VariantCaller/BAMWalkerEngine.cpp
  VariantCaller/OrderedBAMWriter.cpp
  VariantCaller/SampleManager.cpp
  VariantCaller/TargetsManager.cpp
  VariantCaller/MetricsManager.cpp

  VariantCaller/Bookkeeping/MiscUtil.cpp
  VariantCaller/Bookkeeping/ExtendParameters.cpp
  VariantCaller/Bookkeeping/InputStructures.cpp
  VariantCaller/Bookkeeping/VcfFormat.cpp

  VariantCaller/Reads/ExtendedReadInfo.cpp

  # TODO: Actually build vcflib as a static library and link to variant caller.
  # TODO2: Resolve bgzf.c collisions between vcflib and bamtools
  ${ION_VCFLIB_DIR}/Variant.cpp
  ${ION_VCFLIB_DIR}/split.cpp
  ${ION_VCFLIB_DIR}/tabixpp/tabix.cpp
  ${ION_VCFLIB_DIR}/tabixpp/index.c
  ${ION_VCFLIB_DIR}/tabixpp/bgzf.c
  ${ION_VCFLIB_DIR}/smithwaterman/LeftAlign.cpp
  ${ION_VCFLIB_DIR}/smithwaterman/Repeats.cpp
  ${ION_VCFLIB_DIR}/smithwaterman/IndelAllele.cpp
  ${ION_VCFLIB_DIR}/smithwaterman/SmithWatermanGotoh.cpp
  ${ION_FREEBAYES_DIR}/src/AlleleParser.cpp

  Util/OptArgs.cpp
  realignment/Realigner.cpp
  ${ION_TS_EXTERNAL}/jsoncpp-src-amalgated0.6.0-rc1/jsoncpp.cpp
  ${PROJECT_BINARY_DIR}/IonVersion.cpp
)


add_executable(bamwalker ${tvcSRCS})
target_link_libraries(bamwalker ${ION_BAMTOOLS_LIBS} ${EXTRA_LIBS} z pthread)
add_dependencies(bamwalker IONVERSION bamtools htslib_proj)
install(TARGETS bamwalker DESTINATION bin)


# Setup for testing and gtest in particular, only needed once in file before tests.
enable_testing()
find_package(GTest QUIET)
if(GTEST_FOUND)

# Copy over test data...
   configure_file (
     "utest/reference.txt"
     "${CMAKE_CURRENT_BINARY_DIR}/reference.txt"
     COPYONLY
    )
   configure_file (
     "utest/traces.txt"
     "${CMAKE_CURRENT_BINARY_DIR}/traces.txt"
      COPYONLY
   )
   configure_file (
     "utest/reference.txt"
     "${CMAKE_CURRENT_BINARY_DIR}/reference.txt"
     COPYONLY
   )

        configure_file (
          "utest/comp_clean_test_data.txt"
          "${CMAKE_CURRENT_BINARY_DIR}/comp_clean_test_data.txt"
          COPYONLY
        )

        configure_file (
          "utest/comp_clean_test_ev_3.txt"
          "${CMAKE_CURRENT_BINARY_DIR}/comp_clean_test_ev_3.txt"
          COPYONLY
        )

        add_executable(IonErr_Test utest/IonErr_Test.cpp)
        target_link_libraries(IonErr_Test ion-analysis ${GTEST_BOTH_LIBRARIES} pthread)
        add_test(IonErrTest IonErr_Test --gtest_output=xml:./)

        add_executable(OptArgs_Test utest/OptArgs_Test.cpp)
        target_link_libraries(OptArgs_Test ion-analysis ${GTEST_BOTH_LIBRARIES} pthread)
        add_test(OptArgsTest OptArgs_Test --gtest_output=xml:./)

        add_executable(OptionArgs_Test utest/OptionArgs_Test.cpp)
        target_link_libraries(OptionArgs_Test ion-analysis ${GTEST_BOTH_LIBRARIES} pthread)
        add_test(OptionArgsTest OptionArgs_Test --gtest_output=xml:./)

endif()



include(../buildTools/cmake/CMakeLists.cpack.txt)
#% nice to automate this
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libatlas3gf-base,
    libblas3gf,
    libc6,
    libgcc1,
    libstdc++6,
")

set(CPACK_PACKAGE_DESCRIPTION "ion-analysis provides the command line workflow and supporting tools to process raw DAT files into bam and  fastq files.")
set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA
    ${CMAKE_CURRENT_BINARY_DIR}/debian/postinst
    ${CMAKE_CURRENT_BINARY_DIR}/debian/prerm
)
include(CPack)
