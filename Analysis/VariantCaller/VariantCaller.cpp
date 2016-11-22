/* Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     VariantCaller.cpp
//! @ingroup  VariantCaller
//! @brief    HP Indel detection


#include <string>
#include <vector>
#include <stdio.h>
#include <pthread.h>
#include <armadillo>

#include "ExtendParameters.h"

#include "InputStructures.h"
#include "HandleVariant.h"
#include "ReferenceReader.h"
#include "OrderedVCFWriter.h"
#include "BAMWalkerEngine.h"
#include "SampleManager.h"
#include "ExtendedReadInfo.h"
#include "TargetsManager.h"
#include "HotspotReader.h"
#include "MetricsManager.h"
#include "DecisionTreeData.h"

#include "IonVersion.h"

#include <boost/math/distributions/poisson.hpp>
#include "tvcutils/viterbi.h"
#include "tvcutils/unify_vcf.h"

#include "IndelAssembly/IndelAssembly.h"

using namespace std;


void TheSilenceOfTheArmadillos(ofstream &null_ostream)
{
  // Disable armadillo warning messages.
  arma::set_stream_err1(null_ostream);
  arma::set_stream_err2(null_ostream);
}

void * VariantCallerWorker(void *input);

int main(int argc, char* argv[])
{

  printf("tvcgen %s-%s (%s) - Torrent Variant Caller (hacked)\n\n",
         IonVersion::GetVersion().c_str(), IonVersion::GetRelease().c_str(), IonVersion::GetGitHash().c_str());

  // stolen from "Analysis" to silence error messages from Armadillo library
  ofstream null_ostream("/dev/null"); // must stay live for entire scope, or crash when writing
  TheSilenceOfTheArmadillos(null_ostream);

  time_t start_time = time(NULL);


  ExtendParameters parameters(argc, argv);


  mkdir(parameters.outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  if (parameters.program_flow.rich_json_diagnostic || parameters.program_flow.minimal_diagnostic) {
    // make output directory "side effect bad"
    parameters.program_flow.json_plot_dir = parameters.outputDir + "/json_diagnostic/";
    mkdir(parameters.program_flow.json_plot_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }

  ReferenceReader ref_reader;
  ref_reader.Initialize(parameters.fasta);

  TargetsManager targets_manager;
  targets_manager.Initialize(ref_reader, parameters.targets);

  BAMWalkerEngine bam_walker;
  bam_walker.Initialize(ref_reader, targets_manager, parameters.bams, parameters.postprocessed_bam);
  bam_walker.GetProgramVersions(parameters.basecaller_version, parameters.tmap_version);

  SampleManager sample_manager;
  string sample_name = parameters.sampleName;
  sample_manager.Initialize(bam_walker.GetBamHeader(), parameters.sampleName, parameters.force_sample_name);
  if ((sample_name == "") and (sample_manager.num_samples_ > 1)) {
    parameters.multisample = true;
    cout << "Multisample run detected." << endl;
    cerr << "FATAL ERROR: Multisample runs are not currently supported." << endl;
    exit(-1);
  }

  InputStructures global_context;
  global_context.Initialize(parameters, ref_reader, bam_walker.GetBamHeader());

  OrderedVCFWriter vcf_writer;
  vcf_writer.Initialize(parameters.outputDir + "/" + parameters.outputFile, parameters, ref_reader, sample_manager);

  OrderedBAMWriter bam_writer;

  HotspotReader hotspot_reader;
  hotspot_reader.Initialize(ref_reader, parameters.variantPriorsFile);
  if (!parameters.blacklistFile.empty()) {
    if (parameters.variantPriorsFile.empty()) {hotspot_reader.Initialize(ref_reader);}
    hotspot_reader.MakeHintQueue(parameters.blacklistFile);
  }
  string parameters_file = parameters.opts.GetFirstString('-', "parameters-file", "");

  IndelAssemblyArgs parsed_opts;
  parsed_opts.setReference(parameters.fasta);
  parsed_opts.setBams(parameters.bams);
  parsed_opts.setTargetFile(parameters.targets);
  parsed_opts.setOutputVcf(parameters.outputDir + "/indel_assembly.vcf");
  parsed_opts.setParametersFile(parameters_file);
  parsed_opts.setSampleName(sample_name);
  // Print the indel_assembly parameters if do_indel_assembly = true
  if(parameters.program_flow.do_indel_assembly)
    parsed_opts.processParameters(parameters.opts);
  else
    cout<<"TVC: Indel assembly off."<< endl;

  IndelAssembly indel_assembly(&parsed_opts, &ref_reader, &sample_manager, &targets_manager);

  // set up producer of variants
  AlleleParser candidate_generator(parameters, ref_reader, sample_manager, vcf_writer, hotspot_reader);

  MetricsManager metrics_manager;

  VariantCallerContext vc;
  vc.ref_reader = &ref_reader;
  vc.targets_manager = &targets_manager;
  vc.bam_walker = &bam_walker;
  vc.parameters = &parameters;
  vc.global_context = &global_context;
  vc.candidate_generator = &candidate_generator;
  vc.vcf_writer = &vcf_writer;
  vc.bam_writer = &bam_writer;
  vc.metrics_manager = &metrics_manager;
  vc.sample_manager  = &sample_manager;
  vc.indel_assembly = &indel_assembly;

  pthread_mutex_init(&vc.candidate_generation_mutex, NULL);
  pthread_mutex_init(&vc.read_loading_mutex, NULL);
  pthread_mutex_init(&vc.bam_walker_mutex, NULL);
  pthread_mutex_init(&vc.read_removal_mutex, NULL);
  pthread_cond_init(&vc.memory_contention_cond, NULL);
  pthread_cond_init(&vc.alignment_tail_cond, NULL);

  vc.bam_walker->openDepth(parameters.outputDir + "/depth.txt");

  vc.candidate_counter = 0;
  vc.dot_time = time(NULL) + 30;
  //vc.candidate_dot = 0;

  pthread_t worker_id[parameters.program_flow.nThreads];
  for (int worker = 0; worker < parameters.program_flow.nThreads; worker++)
    if (pthread_create(&worker_id[worker], NULL, VariantCallerWorker, &vc)) {
      printf("*Error* - problem starting thread\n");
      exit(-1);
    }

  for (int worker = 0; worker < parameters.program_flow.nThreads; worker++)
    pthread_join(worker_id[worker], NULL);

  pthread_mutex_destroy(&vc.candidate_generation_mutex);
  pthread_mutex_destroy(&vc.read_loading_mutex);
  pthread_mutex_destroy(&vc.bam_walker_mutex);
  pthread_mutex_destroy(&vc.read_removal_mutex);
  pthread_cond_destroy(&vc.memory_contention_cond);
  pthread_cond_destroy(&vc.alignment_tail_cond);

  vector<MergedTarget>::iterator depth_target = targets_manager.merged.begin();
  Alignment* save_list = vc.bam_writer->process_new_entries(vc.bam_walker->alignments_first_);
  vc.bam_walker->SaveAlignments(save_list, vc, depth_target);
  vc.bam_walker->FinishReadRemovalTask(save_list, -1);
  save_list = vc.bam_writer->flush();
  //vc.bam_walker->SaveAlignments(save_list, vc, depth_target);
  vc.bam_walker->FinishReadRemovalTask(save_list, -1);

  vc.bam_walker->closeDepth(vc.ref_reader);

  vcf_writer.Close();
  bam_walker.Close();
  metrics_manager.FinalizeAndSave(parameters.outputDir + "/tvc_metrics.json");

  cerr << endl;
  cout << endl;
  cout << "[tvc] Processing time: " << (time(NULL)-start_time) << " seconds." << endl;


  if(vc.parameters->program_flow.do_indel_assembly){
    indel_assembly.onTraversalDone();
  }
  else{
    indel_assembly.out.close();
  }

  string novel_vcf = parameters.outputDir + "/" + parameters.outputFile;
  string assembly_vcf = parameters.outputDir + "/indel_assembly.vcf";
  string hotspots_file = parameters.variantPriorsFile;
  string output_vcf = parameters.outputDir + "/TSVC_variants.vcf";
  string tvc_metrics = parameters.outputDir + "/tvc_metrics.json";
  string input_depth = parameters.outputDir + "/depth.txt";
  string output_genome_vcf = parameters.outputDir + "/TSVC_variants.genome.vcf";

  if (!parameters.postprocessed_bam.empty()) {
      int return_code = system(string("samtools index " + parameters.postprocessed_bam).c_str());
  }

  { // block serves as isolation of merging and building tabix index
    // Prepare merger object
    VcfOrderedMerger merger(novel_vcf, assembly_vcf, hotspots_file, output_vcf, tvc_metrics, input_depth, output_genome_vcf, ref_reader, targets_manager, 10, max(0, parameters.minCoverage), true);

    merger.perform();
  }

  build_index(parameters.outputDir + "/TSVC_variants.vcf");
  build_index(parameters.outputDir + "/TSVC_variants.genome.vcf");

  cerr << endl;
  cout << endl;
  cout << "[tvc] Normal termination. Processing time: " << (time(NULL)-start_time) << " seconds." << endl;

  return 0;
}



void * VariantCallerWorker(void *input)
{
  BamAlignment alignment;
  VariantCallerContext& vc = *static_cast<VariantCallerContext*>(input);

  vector<MergedTarget>::iterator indel_target = vc.targets_manager->merged.begin();
  vector<MergedTarget>::iterator depth_target = vc.targets_manager->merged.begin();

  deque<VariantCandidate> variant_candidates;
  const static int kReadBatchSize = 40;
  Alignment * new_read[kReadBatchSize];
  bool success[kReadBatchSize];
  list<PositionInProgress>::iterator position_ticket;
  PersistingThreadObjects  thread_objects(*vc.global_context);
  bool more_positions = true;

  pthread_mutex_lock(&vc.bam_walker_mutex);
  MetricsAccumulator& metrics_accumulator = vc.metrics_manager->NewAccumulator();
  pthread_mutex_unlock(&vc.bam_walker_mutex);

  while (true /*more_positions*/) {

    // Opportunistic read removal

    if (vc.bam_walker->EligibleForReadRemoval()) {
      if (pthread_mutex_trylock(&vc.read_removal_mutex) == 0) {

        Alignment *removal_list = NULL;
        pthread_mutex_lock(&vc.bam_walker_mutex);
        vc.bam_walker->RequestReadRemovalTask(removal_list);
        pthread_mutex_unlock(&vc.bam_walker_mutex);

        // In rare cases, the Eligible check passes, but another thread gets to
        // remove the reads first. When this thread has the lock, it finds
        // there are no reads to remove. That cause unexpected behavior in
        // SaveAlignment(): when it sees NULL, it saves all reads and removes
        // ZM tags.  To prevent that, we need to check for an empty read set.
        if (removal_list) {
          Alignment* save_list = vc.bam_writer->process_new_entries(removal_list);
          vc.bam_walker->SaveAlignments(save_list, vc, depth_target);
          pthread_mutex_lock(&vc.bam_walker_mutex);
          vc.bam_walker->FinishReadRemovalTask(save_list);
          pthread_mutex_unlock(&vc.bam_walker_mutex);
        }

        pthread_mutex_unlock(&vc.read_removal_mutex);

        pthread_cond_broadcast(&vc.memory_contention_cond);
      }
    }

    // If too many reads are in memory and at least one candidate evaluator is in progress, pause this thread.
    // Wake up when the oldest candidate evaluator task is complete.
    pthread_mutex_lock(&vc.bam_walker_mutex);
    if (vc.bam_walker->MemoryContention()) {
      pthread_cond_wait (&vc.memory_contention_cond, &vc.bam_walker_mutex);
      pthread_mutex_unlock(&vc.bam_walker_mutex);
      continue;
    }
    pthread_mutex_unlock(&vc.bam_walker_mutex);


    //
    // Task dispatch: Decide whether to load more reads or to generate more variants
    //
    bool ready_for_next_position = false;
    if (vc.bam_walker->EligibleForGreedyRead()) {
      // Greedy reading allowed: if candidate generation is in progress, just grab a new read
      if (pthread_mutex_trylock(&vc.candidate_generation_mutex) == 0) {
        pthread_mutex_lock(&vc.bam_walker_mutex);
        ready_for_next_position = vc.bam_walker->ReadyForNextPosition();
        if (not ready_for_next_position) {
          pthread_mutex_unlock(&vc.bam_walker_mutex);
          pthread_mutex_unlock(&vc.candidate_generation_mutex);
        }
      }

    }
    else {
      // Greedy reading disallowed: if candidate generation is in progress,
      // wait for it to finish before deciding what to do.
      pthread_mutex_lock(&vc.candidate_generation_mutex);
      pthread_mutex_lock(&vc.bam_walker_mutex);
      ready_for_next_position = vc.bam_walker->ReadyForNextPosition();
      if (not ready_for_next_position) {
        pthread_mutex_unlock(&vc.bam_walker_mutex);
        pthread_mutex_unlock(&vc.candidate_generation_mutex);
      }
    }

    //
    // Dispatch outcome: Load and process more reads
    //
    if (not ready_for_next_position) {
      pthread_mutex_lock(&vc.read_loading_mutex);
      if (not vc.bam_walker->HasMoreAlignments()) {
        pthread_mutex_unlock(&vc.read_loading_mutex);
        pthread_cond_broadcast(&vc.memory_contention_cond);
        break;
      }

      pthread_mutex_lock(&vc.bam_walker_mutex);
      for (int i = 0; i < kReadBatchSize; ++i) {
        vc.bam_walker->RequestReadProcessingTask(new_read[i]);
        success[i] = false;
      }
      pthread_mutex_unlock(&vc.bam_walker_mutex);

      for (int i = 0; i < kReadBatchSize; ++i) {
        success[i] = vc.bam_walker->GetNextAlignmentCore(new_read[i], vc, indel_target);
        if (not success[i])
          break;
      }
      pthread_mutex_unlock(&vc.read_loading_mutex);


      for (int i = 0; i < kReadBatchSize and success[i]; ++i) {
        // 1) Fill in read body information and apply the initial set of read filters
        if (not vc.candidate_generator->BasicFilters(*new_read[i])) {
          // cerr << "filtered by BasicFilters(): " << new_read[i]->alignment.Name << endl;
          continue;
        }

        // // 2) Alignment information-altering methods (can also filter a read)
        // do_something(new_read[i]);
        // if (new_read[i]->filtered) {
        //   cerr << "filtered by do_something(): " << new_read[i]->alignment.Name << endl;
        //   continue;
        // }

        // 3) Parsing alignment: Read filtering & populating allele specific data types in an Alignment object
        // cerr << "UnpackReadAlleles( " << new_read[i]->alignment.Name << " )\n" << flush;
        vc.candidate_generator->UnpackReadAlleles(*new_read[i]);

        // 4) Unpacking read metadata for the evaluator
        UnpackOnLoad(new_read[i], *vc.global_context);
      }

      pthread_mutex_lock(&vc.bam_walker_mutex);
      for (int i = 0; i < kReadBatchSize; ++i) {
        vc.bam_walker->FinishReadProcessingTask(new_read[i], success[i]);
      }
      pthread_mutex_unlock(&vc.bam_walker_mutex);

      continue;
    }

    //
    // Dispatch outcome: Generate candidates at next position
    //
    if (not vc.bam_walker->HasMoreAlignments() and not more_positions) {
      pthread_mutex_unlock(&vc.bam_walker_mutex);
      pthread_mutex_unlock(&vc.candidate_generation_mutex);
      pthread_cond_broadcast(&vc.memory_contention_cond);
      break;
    }

    vc.bam_walker->BeginPositionProcessingTask(position_ticket);
    pthread_mutex_unlock(&vc.bam_walker_mutex);

    int sample_num = 1;

    int haplotype_length = 1;
    vc.candidate_generator->GenerateCandidates(variant_candidates, position_ticket, haplotype_length);

    pthread_mutex_lock(&vc.bam_walker_mutex);
    int next_hotspot_chr = -1;
    long next_hotspot_position = -1;
    if (vc.candidate_generator->GetNextHotspotLocation(next_hotspot_chr, next_hotspot_position))
      more_positions = vc.bam_walker->AdvancePosition(haplotype_length, next_hotspot_chr, next_hotspot_position);
    else
      more_positions = vc.bam_walker->AdvancePosition(haplotype_length);
    pthread_mutex_unlock(&vc.bam_walker_mutex);

    if (not variant_candidates.empty()) {
      // cerr << " Got candidates!" << endl;

      int vcf_writer_slot = vc.vcf_writer->ReserveSlot();
      vc.candidate_counter += variant_candidates.size();
      //while (vc.candidate_counter > vc.candidate_dot) {
      //  cerr << ".";
      //  /*
      //  pthread_mutex_lock(&vc.bam_walker_mutex);
      //  vc.bam_walker->PrintStatus();
      //  pthread_mutex_unlock(&vc.bam_walker_mutex);
      //  */
      //  vc.candidate_dot += 50;
      //}
      if (time(NULL) > vc.dot_time) {
        cerr << '.';
        vc.dot_time += 30;
      }

      pthread_mutex_unlock(&vc.candidate_generation_mutex);
      // Separate the queuing of variants from actual work of calling variants
      for (deque<VariantCandidate>::iterator v = variant_candidates.begin(); v != variant_candidates.end(); ++v) {
        if (vc.parameters->multisample) {
          bool pass = false; // if pass == false there are no reads for the candidate
          bool filter = true;
          for (int sample_index = 0; (sample_index < vc.sample_manager->num_samples_); ++sample_index) {
            pass = true;
            if (!v->variant.isFiltered) {filter = false;}
            v->variant.isFiltered = false;
          }
          if (filter) {
            cerr << "------------------- NOCALL: filtered candidate ----------------\n";
            v->variant.filter = "NOCALL";
            v->variant.isFiltered = true;
          }
          else {
            v->variant.filter = "PASS";
            v->variant.isFiltered = false;
          }
          if (!pass) {
            for (int sample_index = 0; (sample_index < vc.sample_manager->num_samples_); ++sample_index) {
              AutoFailTheCandidate(v->variant, vc.parameters->my_controls.use_position_bias, v->variant.sampleNames[sample_index]);
            }
          }
        } /* multisample */
        else {
          if (!ProcessOneVariant(thread_objects, vc, *v, *position_ticket, 0)) {
            AutoFailTheCandidate(v->variant, vc.parameters->my_controls.use_position_bias, v->variant.sampleNames[0]);
          }
        }
      }

      // cerr << "writing candidates\n";
      vc.vcf_writer->WriteSlot(vcf_writer_slot, variant_candidates);

      variant_candidates.clear();
    }
    else {
      pthread_mutex_unlock(&vc.candidate_generation_mutex);
    }

    metrics_accumulator.CollectMetrics(position_ticket, haplotype_length, vc.ref_reader);

    pthread_mutex_lock(&vc.bam_walker_mutex);
    bool signal_contention_removal = vc.bam_walker->IsEarlierstPositionProcessingTask(position_ticket);
    vc.bam_walker->FinishPositionProcessingTask(position_ticket);
    pthread_mutex_unlock(&vc.bam_walker_mutex);

    if (signal_contention_removal)
      pthread_cond_broadcast(&vc.memory_contention_cond);

  }
  return NULL;
}





