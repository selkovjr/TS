/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

#include "DiagnosticJson.h"
//----------------------output some diagnostic information below---------------

void DiagnosticWriteJson(const Json::Value & json, const std::string& filename_json) {
  std::ofstream outJsonFile(filename_json.c_str(), std::ios::out);
  if (outJsonFile.good())
    outJsonFile << json.toStyledString();
  else
    std::cerr << "[tvc] diagnostic unable to write JSON file " << filename_json << std::endl;
  outJsonFile.close();
}

void DiagnosticJsonReadStack(Json::Value &json, const vector<const Alignment *>& read_stack, const InputStructures &global_context) {

  for (unsigned int i_read = 0; i_read < read_stack.size(); i_read++) {
    json["MapQuality"][i_read] = read_stack[i_read]->alignment.MapQuality;
  }
}

void DiagnosticJsonCrossHypotheses(Json::Value &json, const CrossHypotheses &my_cross) {
  // output relevant data from the individual read hypothesis tester
  json["strand"] = my_cross.strand_key;
  json["success"] = my_cross.success ? 1 : 0;

  json["heavy"] = my_cross.heavy_tailed;

  // sequence, responsibility (after clustering), etc
  for (unsigned int i_hyp = 0; i_hyp < my_cross.responsibility.size(); i_hyp++) {
    json["instancebystate"][i_hyp] = my_cross.instance_of_read_by_state[i_hyp];
    json["responsibility"][i_hyp] = my_cross.responsibility[i_hyp];
    json["loglikelihood"][i_hyp] = my_cross.log_likelihood[i_hyp];
    json["scaledlikelihood"][i_hyp] = my_cross.scaled_likelihood[i_hyp];
  }
}

void DiagnosticJsonCrossStack(Json::Value &json, const HypothesisStack &hypothesis_stack) {
  for (unsigned int i_read = 0; i_read < hypothesis_stack.total_theory.my_hypotheses.size(); i_read++) {
    DiagnosticJsonCrossHypotheses(json["Cross"][i_read], hypothesis_stack.total_theory.my_hypotheses[i_read]);
  }
}


void DiagnosticJsonHistory(Json::Value &json, const HypothesisStack &hypothesis_stack){
  for (unsigned int i_start=0; i_start<hypothesis_stack.ll_record.size(); i_start++){
    json["LLrecord"][i_start] = hypothesis_stack.ll_record[i_start];
  }
}

void TinyDiagnosticOutput(const vector<const Alignment *>& read_stack, const HypothesisStack &hypothesis_stack,
    const string& variant_contig, int variant_position, const string& ref_allele, const string& var_allele,
    const InputStructures &global_context, const string &out_dir){
  string outFile;
  Json::Value diagnostic_json;

  outFile = out_dir + variant_contig + "."
            + convertToString(variant_position) + "."
            + ref_allele + "." + var_allele + ".tiny.json";
  // just a little bit of data
  DiagnosticJsonReadStack(diagnostic_json["ReadStack"], read_stack, global_context);
  // TinyDiagnosticJsonCrossStack(diagnostic_json["CrossHypotheses"], hypothesis_stack);
  // write it out
  DiagnosticWriteJson(diagnostic_json, outFile);
}

void RichDiagnosticOutput(const vector<const Alignment *>& read_stack, const HypothesisStack &hypothesis_stack,
    const string& variant_contig, int variant_position, const string& ref_allele, const string& var_allele,
    const InputStructures &global_context, const string  &out_dir) {
  string outFile;
  Json::Value diagnostic_json;

  outFile = out_dir + variant_contig + "."
            + convertToString(variant_position) + "."
            + ref_allele + "." + var_allele + ".diagnostic.json";

  cerr << "RichDiagnosticOutput() -> " << outFile << endl;

  diagnostic_json["MagicNumber"] = 12;

  DiagnosticJsonReadStack(diagnostic_json["ReadStack"], read_stack, global_context);
  DiagnosticJsonCrossStack(diagnostic_json["CrossHypotheses"], hypothesis_stack);
  DiagnosticJsonHistory(diagnostic_json["History"],hypothesis_stack);

  DiagnosticWriteJson(diagnostic_json, outFile);
}

void JustOneDiagnosis(const Evaluator &eval, const InputStructures &global_context,
    const string &out_dir, bool rich_diag)
{
  //diagnose one particular variant
  // check against a list?
  // only do this if using a small VCF for input of variants

  cerr << "JustOneDiagnosis() rich = " << rich_diag << endl;

  // build a unique identifier to write out diagnostics
  int variant_position = eval.variant->position;
  string ref_allele = eval.variant->ref;
  string var_allele = eval.variant->alt[0];
  for (unsigned int i_allele=1; i_allele<eval.variant->alt.size(); i_allele++) {
    var_allele += ',';
    var_allele += eval.variant->alt[i_allele];
  }
  string variant_contig =  eval.variant->sequenceName;

  if (rich_diag)
    RichDiagnosticOutput(eval.read_stack, eval.allele_eval,
        variant_contig, variant_position, ref_allele, var_allele, global_context, out_dir);
  else
    TinyDiagnosticOutput(eval.read_stack, eval.allele_eval,
        variant_contig, variant_position, ref_allele, var_allele, global_context, out_dir);
}
