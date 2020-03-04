# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved
cmake_minimum_required (VERSION 2.6)
project (ion-analysis)
set(PROJECT_DESCRIPTION "Tempus variant candidate generator (extracted from Torrent Variant Caller)")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckCCompilerFlag)
check_c_compiler_flag( -malign-double HAS_ALIGN_DOUBLE )

set(CMAKE_CXX_FLAGS -g3)

### To generate coverage report:
###
###   rm -rf build; ./buildTools/build.sh -DCMAKE_BUILD_TYPE=Debug
###   cd build/Analysis
###   make code_coverage
###
### The results will be written to ./code_coverage
###
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
#include(CodeCoverage)
#APPEND_COVERAGE_COMPILER_FLAGS()
#setup_target_for_coverage_lcov(NAME code_coverage EXECUTABLE tvc-candidate-generator --output-dir . --debug 2 --reference /data3/selkov_workdir/data/reference/human_g1k_v37.fasta --target-file ../../xt_probe_gene_targets.bed --parameters-file "/data3/selkov_workdir/src/TS/plugin/variantCaller/pluginMedia/configs/somatic_low_stringency.json" --use-best-n-alleles 100 --min-var-freq 0.001 --snp-min-allele-freq 0.001 --gen-min-alt-allele-freq 0.001 --gen-min-indel-alt-allele-freq 0.001 --allow-snps on --allow-mnps off --input-bam ../../test.bam --output-vcf test.vcf DEPENDENCY tvc-candidate-generator)

include(../buildTools/cmake/CMakeLists.version.txt)
include(../buildTools/cmake/CMakeLists.compiler.txt)


# Where to install - override with: -DION_INSTALL_PREFIX
set(ION_INSTALL_PREFIX "/opt/ion" CACHE PATH "Ion Install Prefix")

set(CPACK_PACKAGING_INSTALL_PREFIX ${ION_INSTALL_PREFIX})
set(CPACK_SET_DESTDIR "ON")

configure_file (
  "${PROJECT_SOURCE_DIR}/debian/postinst.in"
  "${PROJECT_BINARY_DIR}/debian/postinst" @ONLY
)

include(../buildTools/cmake/CMakeLists.dependencies.txt)

include_directories("${PROJECT_SOURCE_DIR}/Util")

# Ion Analysis Library
add_library(ion-analysis
  ${PROJECT_BINARY_DIR}/TVCVersion.cpp

  Util/IonErr.cpp
  Util/OptArgs.cpp
  Util/OptionArgsBase.cpp
  Util/OptionArgs.cpp
  Util/LoggerBase.cpp
  Util/Logger.cpp

  LinuxCompat.cpp

  ${ION_TS_EXTERNAL}/jsoncpp-src-amalgated0.6.0-rc1/jsoncpp.cpp
)
add_dependencies(ion-analysis TVCVERSION samtools bamtools)

target_link_libraries(ion-analysis z)

set_target_properties(ion-analysis PROPERTIES
  VERSION   "${TVC_VERSION_MAJOR}.${TVC_VERSION_MINOR}.${TVC_VERSION_RELEASE}.${TVC_VERSION_BUILDNUM}"
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

set(SOURCESS
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
  ${ION_TS_EXTERNAL}/jsoncpp-src-amalgated0.6.0-rc1/jsoncpp.cpp
  ${PROJECT_BINARY_DIR}/TVCVersion.cpp
)


add_executable(tvc-candidate-generator ${SOURCESS})
target_link_libraries(tvc-candidate-generator ${ION_BAMTOOLS_LIBS} z pthread)
add_dependencies(tvc-candidate-generator TVCVERSION bamtools htslib_proj)
install(TARGETS tvc-candidate-generator DESTINATION bin)


include(../buildTools/cmake/CMakeLists.cpack.txt)
#% nice to automate this
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libatlas3gf-base,
  libblas3gf,
  libc6,
  libgcc1,
  libstdc++6,
")

set(CPACK_PACKAGE_DESCRIPTION "ion-analysis provides the command line workflow and supporting tools to process raw DAT files into bam and fastq files.")
set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA
  ${CMAKE_CURRENT_BINARY_DIR}/debian/postinst
  ${CMAKE_CURRENT_BINARY_DIR}/debian/prerm
)
include(CPack)
