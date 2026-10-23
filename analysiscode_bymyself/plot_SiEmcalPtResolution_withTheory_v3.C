#include "/sphenix/tg/tg01/commissioning/INTT/repositories/macros/common/sPhenixStyle.C"
#include "/sphenix/u/tomoki/work_tomoki/Signature_pp.hh"

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
constexpr double kPtCoefficient = 0.21;
constexpr double kPtExponent = 0.986;

struct Cuts
{
  double max_abs_dz = 10.0;
  double max_abs_dphi = 0.06;
  double max_track_chi2 = 5000.0;
  double min_calo_chi2 = -100.0;
  double max_calo_chi2 = 48.0;
};

struct SampleConfig
{
  std::string tag;
  std::string title;
  int pid_selection = 0;
  std::string center_pattern;
  std::string first_pattern;
  std::string position_pattern;
};

struct FitValue
{
  bool valid = false;
  double mean = 0.0;
  double sigma = 0.0;
  double sigma_error = 0.0;
  double half68 = 0.0;
  double entries = 0.0;
};

struct Histograms
{
  std::vector<std::unique_ptr<TH1D>> silicon;
  std::vector<std::unique_ptr<TH1D>> center;
  std::vector<std::unique_ptr<TH1D>> first;
  std::vector<std::unique_ptr<TH1D>> cluster_phi;
};

struct PanelResult
{
  SampleConfig config;
  Histograms histograms;
  std::unique_ptr<TGraphErrors> silicon;
  std::unique_ptr<TGraphErrors> center;
  std::unique_ptr<TGraphErrors> first;
  std::unique_ptr<TGraphErrors> position_prediction;
};

int AddFiles(TChain& chain, const std::string& pattern,
             int first_process, int num_files, const std::string& label)
{
  int added = 0;
  for (int offset = 0; offset < num_files; ++offset)
  {
    const std::string path = Form(pattern.c_str(), first_process + offset);
    if (gSystem->AccessPathName(path.c_str())) continue;
    if (chain.Add(path.c_str()) > 0) ++added;
  }
  std::cout << label << ": added " << added << " / " << num_files
            << " files, entries=" << chain.GetEntries() << std::endl;
  return added;
}

bool HasBranches(TChain& chain, const std::vector<std::string>& names,
                 const std::string& label)
{
  for (const std::string& name : names)
  {
    if (!chain.GetBranch(name.c_str()))
    {
      std::cerr << label << ": missing branch " << name << std::endl;
      return false;
    }
  }
  return true;
}

size_t FindPtBin(double pt, const std::vector<double>& edges)
{
  const auto upper = std::upper_bound(edges.begin(), edges.end(), pt);
  if (upper == edges.begin() || upper == edges.end()) return edges.size();
  return static_cast<size_t>(upper - edges.begin() - 1);
}

bool PassCuts(float match_dphi, float match_dz,
              float track_chi2, float calo_chi2, const Cuts& cuts)
{
  return std::isfinite(match_dphi) && std::isfinite(match_dz) &&
         std::isfinite(track_chi2) && std::isfinite(calo_chi2) &&
         std::abs(match_dphi) < cuts.max_abs_dphi &&
         std::abs(match_dz) < cuts.max_abs_dz &&
         track_chi2 < cuts.max_track_chi2 &&
         calo_chi2 > cuts.min_calo_chi2 &&
         calo_chi2 < cuts.max_calo_chi2;
}

std::vector<std::unique_ptr<TH1D>> BookHistograms(
    const std::string& prefix, const std::vector<double>& edges,
    int bins, double minimum, double maximum)
{
  std::vector<std::unique_ptr<TH1D>> result;
  result.reserve(edges.size() - 1);
  for (size_t bin = 0; bin + 1 < edges.size(); ++bin)
  {
    auto histogram = std::make_unique<TH1D>(
        Form("h_%s_ptbin%zu", prefix.c_str(), bin),
        Form("%.2f #leq p_{T}^{truth} < %.2f GeV/c",
             edges[bin], edges[bin + 1]),
        bins, minimum, maximum);
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    result.emplace_back(std::move(histogram));
  }
  return result;
}

FitValue FitCore(TH1D& histogram, int minimum_entries = 20)
{
  FitValue value;
  value.entries = histogram.GetEntries();
  if (value.entries < minimum_entries) return value;

  const double probabilities[3] = {0.16, 0.50, 0.84};
  double quantiles[3] = {0.0, 0.0, 0.0};
  histogram.GetQuantiles(3, quantiles, probabilities);
  value.half68 = 0.5*(quantiles[2] - quantiles[0]);
  if (!std::isfinite(value.half68) || value.half68 <= 0.0) return value;

  const double low = std::max(histogram.GetXaxis()->GetXmin(),
                              quantiles[1] - 2.5*value.half68);
  const double high = std::min(histogram.GetXaxis()->GetXmax(),
                               quantiles[1] + 2.5*value.half68);
  if (!(high > low)) return value;

  TF1 first(Form("f_first_%s", histogram.GetName()), "gaus", low, high);
  first.SetParameters(histogram.GetMaximum(), quantiles[1], value.half68);
  TFitResultPtr first_result = histogram.Fit(&first, "Q0RS");
  if (static_cast<int>(first_result) != 0 ||
      !std::isfinite(first.GetParameter(2)) || first.GetParameter(2) <= 0.0)
    return value;

  const double first_mean = first.GetParameter(1);
  const double first_sigma = std::abs(first.GetParameter(2));
  TF1 second(Form("f_core_%s", histogram.GetName()), "gaus",
             first_mean - 2.0*first_sigma,
             first_mean + 2.0*first_sigma);
  second.SetParameters(first.GetParameter(0), first_mean, first_sigma);
  TFitResultPtr second_result = histogram.Fit(&second, "Q0RS");
  if (static_cast<int>(second_result) != 0) return value;

  value.mean = second.GetParameter(1);
  value.sigma = std::abs(second.GetParameter(2));
  value.sigma_error = second.GetParError(2);
  value.valid = std::isfinite(value.mean) && std::isfinite(value.sigma) &&
                std::isfinite(value.sigma_error) && value.sigma > 0.0 &&
                value.sigma_error >= 0.0;
  return value;
}

void FillPtResiduals(TChain& chain, int pid_selection,
                     std::vector<std::unique_ptr<TH1D>>* silicon,
                     std::vector<std::unique_ptr<TH1D>>& emcal,
                     const std::vector<double>& edges, const Cuts& cuts,
                     const std::string& label)
{
  const std::vector<std::string> required = {
      "truth_pid", "truth_pt", "track_pt", "calo_pt",
      "match_dphi", "match_dz", "track_chi2", "calo_chi2"};
  if (chain.GetEntries() <= 0 || !HasBranches(chain, required, label)) return;

  int truth_pid = 0;
  float truth_pt = 0.F, track_pt = 0.F, calo_pt = 0.F;
  float match_dphi = 0.F, match_dz = 0.F;
  float track_chi2 = 0.F, calo_chi2 = 0.F;
  chain.SetBranchStatus("*", 0);
  for (const std::string& branch : required)
    chain.SetBranchStatus(branch.c_str(), 1);
  chain.SetBranchAddress("truth_pid", &truth_pid);
  chain.SetBranchAddress("truth_pt", &truth_pt);
  chain.SetBranchAddress("track_pt", &track_pt);
  chain.SetBranchAddress("calo_pt", &calo_pt);
  chain.SetBranchAddress("match_dphi", &match_dphi);
  chain.SetBranchAddress("match_dz", &match_dz);
  chain.SetBranchAddress("track_chi2", &track_chi2);
  chain.SetBranchAddress("calo_chi2", &calo_chi2);

  Long64_t accepted = 0;
  const Long64_t entries = chain.GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    chain.GetEntry(entry);
    if (pid_selection != 0 && truth_pid != pid_selection) continue;
    if (!std::isfinite(truth_pt) || !std::isfinite(track_pt) ||
        !std::isfinite(calo_pt) || truth_pt <= 0.F ||
        track_pt <= 0.F || calo_pt <= 0.F) continue;
    if (!PassCuts(match_dphi, match_dz, track_chi2, calo_chi2, cuts)) continue;
    const size_t bin = FindPtBin(truth_pt, edges);
    if (bin >= emcal.size()) continue;
    if (silicon)
      silicon->at(bin)->Fill((track_pt - truth_pt)/truth_pt);
    emcal[bin]->Fill((calo_pt - truth_pt)/truth_pt);
    ++accepted;
  }
  std::cout << label << ": accepted " << accepted << " / " << entries
            << " entries" << std::endl;
  chain.ResetBranchAddresses();
}

void FillPhiResiduals(TChain& chain, int pid_selection,
                      std::vector<std::unique_ptr<TH1D>>& residuals,
                      const std::vector<double>& edges, const Cuts& cuts,
                      const std::string& label)
{
  const std::vector<std::string> required = {
      "truth_pid", "truth_pt", "truth_cemc_valid", "emcal_truth_dphi",
      "match_dphi", "match_dz", "track_chi2", "calo_chi2"};
  if (chain.GetEntries() <= 0 || !HasBranches(chain, required, label)) return;

  int truth_pid = 0, truth_cemc_valid = 0;
  float truth_pt = 0.F, emcal_truth_dphi = 0.F;
  float match_dphi = 0.F, match_dz = 0.F;
  float track_chi2 = 0.F, calo_chi2 = 0.F;
  chain.SetBranchStatus("*", 0);
  for (const std::string& branch : required)
    chain.SetBranchStatus(branch.c_str(), 1);
  chain.SetBranchAddress("truth_pid", &truth_pid);
  chain.SetBranchAddress("truth_pt", &truth_pt);
  chain.SetBranchAddress("truth_cemc_valid", &truth_cemc_valid);
  chain.SetBranchAddress("emcal_truth_dphi", &emcal_truth_dphi);
  chain.SetBranchAddress("match_dphi", &match_dphi);
  chain.SetBranchAddress("match_dz", &match_dz);
  chain.SetBranchAddress("track_chi2", &track_chi2);
  chain.SetBranchAddress("calo_chi2", &calo_chi2);

  Long64_t accepted = 0;
  const Long64_t entries = chain.GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    chain.GetEntry(entry);
    if (pid_selection != 0 && truth_pid != pid_selection) continue;
    if (!truth_cemc_valid || !std::isfinite(truth_pt) ||
        !std::isfinite(emcal_truth_dphi) || truth_pt <= 0.F) continue;
    if (!PassCuts(match_dphi, match_dz, track_chi2, calo_chi2, cuts)) continue;
    const size_t bin = FindPtBin(truth_pt, edges);
    if (bin >= residuals.size()) continue;
    residuals[bin]->Fill(emcal_truth_dphi);
    ++accepted;
  }
  std::cout << label << ": accepted " << accepted << " / " << entries
            << " truth-entry residuals" << std::endl;
  chain.ResetBranchAddresses();
}

FitValue FitPtCore(TH1D& histogram, double fit_min, double fit_max,
                   int minimum_entries = 20)
{
  FitValue value;
  value.entries = histogram.GetEntries();
  if (value.entries < minimum_entries) return value;
  const int first_bin = histogram.GetXaxis()->FindFixBin(fit_min);
  const int last_bin = histogram.GetXaxis()->FindFixBin(fit_max);
  if (histogram.Integral(first_bin, last_bin) < minimum_entries) return value;
  int peak_bin = first_bin;
  for (int bin = first_bin + 1; bin <= last_bin; ++bin)
    if (histogram.GetBinContent(bin) > histogram.GetBinContent(peak_bin))
      peak_bin = bin;
  const double peak = histogram.GetXaxis()->GetBinCenter(peak_bin);
  TF1 gaussian(Form("f_pt_%s", histogram.GetName()), "gaus", fit_min, fit_max);
  gaussian.SetParameters(histogram.GetBinContent(peak_bin), peak, 0.05);
  gaussian.SetParLimits(
      0, 0.0, std::max(1.0, 10.0*histogram.GetBinContent(peak_bin)));
  gaussian.SetParLimits(1, fit_min, fit_max);
  gaussian.SetParLimits(2, 0.002, 0.5);
  const TFitResultPtr fit_result = histogram.Fit(&gaussian, "RQ0S");
  const int status = static_cast<int>(fit_result);
  const int covariance = fit_result.Get() ? fit_result->CovMatrixStatus() : 0;
  value.mean = gaussian.GetParameter(1);
  value.sigma = std::abs(gaussian.GetParameter(2));
  value.sigma_error = gaussian.GetParError(2);
  value.valid = status == 0 && covariance >= 2 &&
      std::isfinite(value.mean) && value.mean >= fit_min &&
      value.mean <= fit_max && std::isfinite(value.sigma) &&
      value.sigma > 0.002 && value.sigma < 0.5 &&
      std::isfinite(value.sigma_error) && value.sigma_error > 0.0;
  return value;
}

std::unique_ptr<TGraphErrors> MakeResolutionGraph(
    std::vector<std::unique_ptr<TH1D>>& histograms,
    const std::vector<double>& edges, const std::string& name,
    int color, int marker, double fit_min, double fit_max)
{
  auto graph = std::make_unique<TGraphErrors>();
  graph->SetName(name.c_str());
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(1.05);
  graph->SetLineWidth(2);
  for (size_t bin = 0; bin < histograms.size(); ++bin)
  {
    const FitValue fit = FitPtCore(*histograms[bin], fit_min, fit_max);
    if (!fit.valid) continue;
    const int point = graph->GetN();
    graph->SetPoint(point, 0.5*(edges[bin] + edges[bin + 1]), fit.sigma);
    graph->SetPointError(point, 0.0,
                         fit.sigma_error);
    std::cout << name << " pT=[" << edges[bin] << "," << edges[bin + 1]
              << ") N=" << fit.entries << " mean=" << fit.mean
              << " sigma=" << fit.sigma << " RMS68=" << fit.half68
              << std::endl;
  }
  return graph;
}

double BendAngleForPt(double pt)
{
  if (!(pt > 0.0)) return 0.0;
  return std::pow(kPtCoefficient/pt, 1.0/kPtExponent);
}

double PropagatePhiToRelativePt(double pt, double sigma_phi)
{
  const double bend = BendAngleForPt(pt);
  if (!(bend > 0.0) || !(sigma_phi >= 0.0)) return 0.0;
  return kPtExponent*sigma_phi/bend;
}

std::unique_ptr<TGraphErrors> MakePositionPrediction(
    std::vector<std::unique_ptr<TH1D>>& phi_histograms,
    const std::vector<double>& edges, const std::string& name,
    int color, int marker, int line_style)
{
  auto graph = std::make_unique<TGraphErrors>();
  graph->SetName(name.c_str());
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineStyle(line_style);
  graph->SetLineWidth(3);
  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(0.9);
  for (size_t bin = 0; bin < phi_histograms.size(); ++bin)
  {
    const FitValue fit = FitCore(*phi_histograms[bin]);
    if (!fit.valid) continue;
    const double pt = 0.5*(edges[bin] + edges[bin + 1]);
    const double scale = kPtExponent/BendAngleForPt(pt);
    const int point = graph->GetN();
    graph->SetPoint(point, pt, scale*fit.sigma);
    graph->SetPointError(point, 0.0,
                         scale*fit.sigma_error);
  }
  return graph;
}


PanelResult BuildPanel(const SampleConfig& config,
                       int first_process, int num_files,
                       const std::vector<double>& edges, const Cuts& cuts)
{
  PanelResult result;
  result.config = config;
  result.histograms.silicon = BookHistograms(
      config.tag + "_silicon", edges, 400, -2.0, 2.0);
  result.histograms.center = BookHistograms(
      config.tag + "_center", edges, 400, -2.0, 2.0);
  result.histograms.first = BookHistograms(
      config.tag + "_first", edges, 400, -2.0, 2.0);
  result.histograms.cluster_phi = BookHistograms(
      config.tag + "_cluster_phi", edges, 300, -0.15, 0.15);

  TChain center("SiEmcalTopoTree");
  TChain first("SiEmcalTopoTree");
  TChain position("SiEmcalTopoTree");
  AddFiles(center, config.center_pattern, first_process, num_files,
           config.title + " center");
  AddFiles(first, config.first_pattern, first_process, num_files,
           config.title + " first-position");
  AddFiles(position, config.position_pattern, first_process, num_files,
           config.title + " truth-position");

  FillPtResiduals(center, config.pid_selection, &result.histograms.silicon,
                  result.histograms.center, edges, cuts,
                  config.title + " center");
  FillPtResiduals(first, config.pid_selection, nullptr,
                  result.histograms.first, edges, cuts,
                  config.title + " first-position");
  FillPhiResiduals(position, config.pid_selection,
                   result.histograms.cluster_phi, edges, cuts,
                   config.title + " position");

  result.silicon = MakeResolutionGraph(
      result.histograms.silicon, edges, "g_" + config.tag + "_silicon",
      kBlack, 20, -0.20, 0.20);
  result.center = MakeResolutionGraph(
      result.histograms.center, edges, "g_" + config.tag + "_center",
      kRed + 1, 21, -0.08, 0.08);
  result.first = MakeResolutionGraph(
      result.histograms.first, edges, "g_" + config.tag + "_first",
      kBlue + 1, 22, -0.08, 0.08);
  result.position_prediction = MakePositionPrediction(
      result.histograms.cluster_phi, edges,
      "g_" + config.tag + "_center_position_prediction",
      kRed + 1, 24, 2);
  return result;
}

void DrawPanel(PanelResult& result, const std::vector<double>& edges,
               double y_max, bool compact)
{
  gPad->SetLeftMargin(compact ? 0.14 : 0.13);
  gPad->SetRightMargin(0.035);
  gPad->SetBottomMargin(compact ? 0.13 : 0.12);
  gPad->SetTopMargin(compact ? 0.13 : 0.12);
  gPad->SetTicks(1, 1);

  TH1* frame = gPad->DrawFrame(edges.front(), 0.0, edges.back(), y_max);
  frame->GetXaxis()->SetTitle("p_{T}^{truth} [GeV/c]");
  frame->GetYaxis()->SetTitle(
      "#sigma_{core}((p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth})");
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleOffset(compact ? 1.30 : 1.15);
  if (compact)
  {
    frame->GetXaxis()->SetTitleSize(0.048);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);
  }

  result.position_prediction->Draw("LP SAME");
  result.silicon->Draw("PE1 SAME");
  result.center->Draw("PE1 SAME");
  result.first->Draw("PE1 SAME");

  const double text_size = compact ? 0.026 : 0.030;
  auto* legend = new TLegend(
      compact ? 0.39 : 0.50, compact ? 0.55 : 0.60, 0.91, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(text_size);
  legend->SetHeader(result.config.title.c_str(), "C");
  legend->AddEntry(result.silicon.get(), "Silicon", "pe");
  legend->AddEntry(result.center.get(), "Silicon + EMCal center", "pe");
  legend->AddEntry(result.first.get(), "Silicon + EMCal first position", "pe");
  legend->AddEntry(result.position_prediction.get(),
                  "From measured center #sigma_{#phi}", "lp");
  legend->Draw();


  signature::text_size = compact ? 0.032 : 0.038;
  signature::pos_x_left = 0.14;
  signature::pos_x_date = 0.80;
  signature::pos_y_top = 0.965;
  signature::WritesPhenix(false);
  signature::WriteDate();
}

void WritePanel(TFile& output, PanelResult& result)
{
  output.cd();
  for (auto* collection : {
           &result.histograms.silicon, &result.histograms.center,
           &result.histograms.first,
           &result.histograms.cluster_phi})
    for (auto& histogram : *collection) histogram->Write();
  result.silicon->Write();
  result.center->Write();
  result.first->Write();
  result.position_prediction->Write();
}
}  // namespace

void plot_SiEmcalPtResolution_withTheory_v3(
    int num_files = 10000,
    int first_process = 0,
    const std::string& output_directory =
        "/sphenix/u/tomoki/work_tomoki/SiCalo/analysiscode_bymyself/"
        "pdf/emcal_pt_resolution_theory_v3",
    double y_max = 0.35,
    double pt_min = 0.0,
    double pt_max = 5.0,
    int num_pt_bins = 50,
    double max_abs_dz = 10.0,
    double max_abs_dphi = 0.06,
    double max_track_chi2 = 5000.0,
    double min_calo_chi2 = -100.0,
    double max_calo_chi2 = 48.0)
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  TH1::AddDirectory(false);
  gSystem->mkdir(output_directory.c_str(), true);

  if (num_files <= 0 || num_pt_bins <= 0 || !(pt_max > pt_min))
  {
    std::cerr << "Invalid file count, pT range, or pT-bin count." << std::endl;
    return;
  }
  std::vector<double> pt_edges(static_cast<size_t>(num_pt_bins) + 1);
  for (int bin = 0; bin <= num_pt_bins; ++bin)
    pt_edges[static_cast<size_t>(bin)] =
        pt_min + (pt_max - pt_min)*bin/num_pt_bins;
  const Cuts cuts{max_abs_dz, max_abs_dphi, max_track_chi2,
                  min_calo_chi2, max_calo_chi2};

  const std::string base =
      "/sphenix/u/tomoki/work_tomoki/SiCalo/output/";
  const std::vector<SampleConfig> configurations = {
      {"single_e", "single gun e^{-}", 11,
       base + "ana_eminus_sicalo_v2/SiEmcalTopoTree_e-_%06d_v4.root",
       base + "ana_eminus_sicalo_v3/SiEmcalTopoTree_e-_%06d_v3.root",
       base + "ana_eminus_sicalo_position_v4/SiEmcalTopoTree_e-_%06d_v4.root"},
      {"single_pi", "single gun #pi^{-}", -211,
       base + "ana_0803piminus_sicalo_v2/SiEmcalTopoTree_pi-_%06d_v4.root",
       base + "ana_0803piminus_sicalo_v3/SiEmcalTopoTree_pi-_%06d_v3.root",
       base + "ana_0803piminus_sicalo_position_v4/SiEmcalTopoTree_pi-_%06d_v4.root"},
      {"pp_e", "p+p e^{-}", 11,
       base + "ana_pythia_sicalo_v2/SiEmcalTopoTree_0000000029-%06d_v4.root",
       base + "ana_pythia_sicalo_v3/SiEmcalTopoTree_0000000029-%06d_v3.root",
       base + "ana_pythia_sicalo_position_v5/SiEmcalTopoTree_0000000029-%06d_v5.root"},
      {"pp_pi", "p+p #pi^{-}", -211,
       base + "ana_pythia_sicalo_v2/SiEmcalTopoTree_0000000029-%06d_v4.root",
       base + "ana_pythia_sicalo_v3/SiEmcalTopoTree_0000000029-%06d_v3.root",
       base + "ana_pythia_sicalo_position_v5/SiEmcalTopoTree_0000000029-%06d_v5.root"}};

  std::vector<PanelResult> results;
  results.reserve(configurations.size());
  for (const SampleConfig& config : configurations)
    results.emplace_back(BuildPanel(
        config, first_process, num_files, pt_edges, cuts));

  TFile output((output_directory +
                "/SiEmcalPtResolution_withTheory_v3.root").c_str(),
               "RECREATE");
  for (PanelResult& result : results) WritePanel(output, result);
  output.Close();

  TCanvas combined("c_combined", "pT resolution with position prediction",
                   1700, 1300);
  combined.Divide(2, 2, 0.002, 0.002);
  for (size_t index = 0; index < results.size(); ++index)
  {
    combined.cd(static_cast<int>(index + 1));
    DrawPanel(results[index], pt_edges, y_max, true);
  }
  combined.SaveAs((output_directory +
                   "/SiEmcalPtResolution_withTheory_4panel_v3.pdf").c_str());
  combined.SaveAs((output_directory +
                   "/SiEmcalPtResolution_withTheory_4panel_v3.png").c_str());

  for (PanelResult& result : results)
  {
    TCanvas canvas(Form("c_%s", result.config.tag.c_str()),
                   result.config.title.c_str(), 1100, 850);
    DrawPanel(result, pt_edges, y_max, false);
    canvas.SaveAs((output_directory + "/SiEmcalPtResolution_withTheory_" +
                   result.config.tag + "_v3.pdf").c_str());
    canvas.SaveAs((output_directory + "/SiEmcalPtResolution_withTheory_" +
                   result.config.tag + "_v3.png").c_str());
  }

  std::cout << "Saved four-panel and individual plots in "
            << output_directory << std::endl;
  std::cout << "Prediction curve uses the measured Gaussian-core EMCal-center phi width "
            << "relative to the truth CEMC crossing and "
            << "pT=" << kPtCoefficient << "*|dphi|^(-" << kPtExponent
            << "). No fixed 2 mm assumption is used." << std::endl;
}
