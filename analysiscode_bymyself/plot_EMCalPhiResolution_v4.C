#include "/sphenix/tg/tg01/commissioning/INTT/repositories/macros/common/sPhenixStyle.C"

#include <TCanvas.h>
#include <TChain.h>
#include <TDatime.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
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
constexpr double kEMCalClusterRadius = 104.0;  // cm

struct Sample
{
  std::string name;
  std::string legend;
  int abs_pid;
  int color;
  int marker;
  std::vector<std::shared_ptr<TH1D>> residuals;
};

int AddFiles(TChain& chain, const std::string& pattern,
             int first_process, int num_files, const std::string& name)
{
  int added = 0;
  for (int offset = 0; offset < num_files; ++offset)
  {
    const std::string path = Form(
        pattern.c_str(), first_process + offset);
    if (gSystem->AccessPathName(path.c_str())) continue;
    if (chain.Add(path.c_str()) > 0) ++added;
  }
  std::cout << name << ": added " << added << " / " << num_files
            << " files" << std::endl;
  return added;
}

bool HasRequiredBranches(TChain& chain, const std::string& name)
{
  for (const char* branch : {
           "truth_pid", "truth_pt", "truth_cemc_valid",
           "truth_cemc_r", "emcal_truth_dphi", "match_dphi",
           "match_dz", "track_chi2", "calo_chi2"})
  {
    if (!chain.GetBranch(branch))
    {
      std::cerr << name << ": missing branch " << branch
                << ". Regenerate with SiliconSeedsAna + "
                   "make_SiEmcalTopoTree_v4.C." << std::endl;
      return false;
    }
  }
  return true;
}

void BookResiduals(Sample& sample, const std::vector<double>& pt_edges)
{
  sample.residuals.clear();
  for (size_t ibin = 0; ibin + 1 < pt_edges.size(); ++ibin)
  {
    const std::string title = Form(
        "%s, %.2f #leq p_{T}^{truth} < %.2f GeV/c;"
        "#phi_{EMCal cluster}-#phi_{truth CEMC entry} [rad];Counts",
        sample.legend.c_str(), pt_edges[ibin], pt_edges[ibin + 1]);
    sample.residuals.emplace_back(std::make_shared<TH1D>(
        Form("h_dphi_%s_ptbin%zu", sample.name.c_str(), ibin),
        title.c_str(), 240, -0.12, 0.12));
    sample.residuals.back()->Sumw2();
  }
}

void FillResiduals(TChain& chain, std::vector<Sample*>& samples,
                   const std::vector<double>& pt_edges,
                   bool apply_quality_cuts,
                   double max_abs_dz, double max_abs_dphi,
                   double max_track_chi2,
                   double min_calo_chi2, double max_calo_chi2)
{
  if (chain.GetEntries() <= 0) return;

  int truth_pid = 0;
  int truth_cemc_valid = 0;
  float truth_pt = 0.F;
  float truth_cemc_r = 0.F;
  float residual = 0.F;
  float match_dphi = 0.F;
  float match_dz = 0.F;
  float track_chi2 = 0.F;
  float calo_chi2 = 0.F;

  chain.SetBranchStatus("*", 0);
  for (const char* branch : {
           "truth_pid", "truth_pt", "truth_cemc_valid",
           "truth_cemc_r", "emcal_truth_dphi", "match_dphi",
           "match_dz", "track_chi2", "calo_chi2"})
    chain.SetBranchStatus(branch, 1);

  chain.SetBranchAddress("truth_pid", &truth_pid);
  chain.SetBranchAddress("truth_pt", &truth_pt);
  chain.SetBranchAddress("truth_cemc_valid", &truth_cemc_valid);
  chain.SetBranchAddress("truth_cemc_r", &truth_cemc_r);
  chain.SetBranchAddress("emcal_truth_dphi", &residual);
  chain.SetBranchAddress("match_dphi", &match_dphi);
  chain.SetBranchAddress("match_dz", &match_dz);
  chain.SetBranchAddress("track_chi2", &track_chi2);
  chain.SetBranchAddress("calo_chi2", &calo_chi2);

  Long64_t valid = 0;
  Long64_t cut_rejected = 0;
  const Long64_t entries = chain.GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    chain.GetEntry(entry);
    if (!truth_cemc_valid || !std::isfinite(truth_pt) ||
        !std::isfinite(truth_cemc_r) || !std::isfinite(residual) ||
        truth_pt <= 0.F || truth_cemc_r <= 0.F) continue;

    if (apply_quality_cuts &&
        (!std::isfinite(match_dphi) || !std::isfinite(match_dz) ||
         !std::isfinite(track_chi2) || !std::isfinite(calo_chi2) ||
         std::abs(match_dz) >= max_abs_dz ||
         std::abs(match_dphi) >= max_abs_dphi ||
         track_chi2 >= max_track_chi2 ||
         calo_chi2 <= min_calo_chi2 ||
         calo_chi2 >= max_calo_chi2))
    {
      ++cut_rejected;
      continue;
    }

    const auto upper = std::upper_bound(
        pt_edges.begin(), pt_edges.end(), truth_pt);
    if (upper == pt_edges.begin() || upper == pt_edges.end()) continue;
    const size_t ptbin =
        static_cast<size_t>(upper - pt_edges.begin() - 1);
    for (Sample* sample : samples)
    {
      if (sample && std::abs(truth_pid) == sample->abs_pid)
        sample->residuals[ptbin]->Fill(residual);
    }
    ++valid;
  }
  std::cout << chain.GetName() << ": entries=" << entries
            << ", valid truth CEMC crossings=" << valid
            << ", quality-cut rejected=" << cut_rejected << std::endl;
  chain.ResetBranchAddresses();
}

bool FitCore(TH1D& histogram, double& sigma, double& sigma_error,
             double& mean, double& half68)
{
  if (histogram.GetEntries() < 30) return false;

  const double probabilities[3] = {0.16, 0.50, 0.84};
  double quantiles[3] = {0., 0., 0.};
  histogram.GetQuantiles(3, quantiles, probabilities);
  half68 = 0.5*(quantiles[2] - quantiles[0]);
  if (!std::isfinite(half68) || half68 <= 0.) return false;

  const double low = std::max(histogram.GetXaxis()->GetXmin(),
                              quantiles[1] - 2.5*half68);
  const double high = std::min(histogram.GetXaxis()->GetXmax(),
                               quantiles[1] + 2.5*half68);
  TF1 first(Form("f_first_%s", histogram.GetName()), "gaus", low, high);
  first.SetParameters(histogram.GetMaximum(), quantiles[1], half68);
  TFitResultPtr first_result = histogram.Fit(&first, "Q0RS");
  if (static_cast<int>(first_result) != 0 ||
      !std::isfinite(first.GetParameter(2)) ||
      first.GetParameter(2) <= 0.) return false;

  const double first_mean = first.GetParameter(1);
  const double first_sigma = std::abs(first.GetParameter(2));
  TF1 second(Form("f_core_%s", histogram.GetName()), "gaus",
             first_mean - 2.0*first_sigma,
             first_mean + 2.0*first_sigma);
  second.SetParameters(first.GetParameter(0), first_mean, first_sigma);
  TFitResultPtr second_result = histogram.Fit(&second, "Q0RS");
  if (static_cast<int>(second_result) != 0) return false;

  mean = second.GetParameter(1);
  sigma = std::abs(second.GetParameter(2));
  sigma_error = second.GetParError(2);
  return std::isfinite(mean) && std::isfinite(sigma) &&
         std::isfinite(sigma_error) && sigma > 0.;
}

std::unique_ptr<TGraphErrors> MakeGraph(
    Sample& sample, const std::vector<double>& pt_edges,
    double scale, const char* suffix, TFile& output)
{
  auto graph = std::make_unique<TGraphErrors>();
  graph->SetName(Form("g_sigma_%s_%s", sample.name.c_str(), suffix));
  graph->SetTitle(sample.legend.c_str());
  graph->SetMarkerColor(sample.color);
  graph->SetLineColor(sample.color);
  graph->SetMarkerStyle(sample.marker);
  graph->SetMarkerSize(1.15);
  graph->SetLineWidth(2);

  for (size_t ibin = 0; ibin < sample.residuals.size(); ++ibin)
  {
    TH1D& histogram = *sample.residuals[ibin];
    output.cd();
    histogram.Write();
    double sigma = 0., error = 0., mean = 0., half68 = 0.;
    if (!FitCore(histogram, sigma, error, mean, half68)) continue;
    const int point = graph->GetN();
    const double center = 0.5*(pt_edges[ibin] + pt_edges[ibin + 1]);
    const double xerror = 0.5*(pt_edges[ibin + 1] - pt_edges[ibin]);
    graph->SetPoint(point, center, scale*sigma);
    graph->SetPointError(point, xerror, scale*error);
    std::cout << sample.name << " pt=[" << pt_edges[ibin] << ","
              << pt_edges[ibin + 1] << ") N=" << histogram.GetEntries()
              << " mean=" << mean << " sigma=" << sigma
              << " half68=" << half68 << std::endl;
  }
  output.cd();
  graph->Write();
  return graph;
}

void DrawHeader()
{
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(0.045);
  text.SetTextAlign(13);
  text.DrawLatex(0.13, 0.965, "#it{#bf{sPHENIX}} Simulation");
  TDatime date;
  text.SetTextAlign(33);
  text.SetTextSize(0.032);
  text.DrawLatex(0.95, 0.965,
                 Form("%d/%d/%d", date.GetMonth(), date.GetDay(),
                      date.GetYear()));
}

void DrawResolution(
    const std::vector<std::unique_ptr<TGraphErrors>>& graphs,
    const std::vector<Sample>& samples,
    double pt_min, double pt_max, double y_max,
    const char* y_title, const std::string& output_path)
{
  TCanvas canvas("c_resolution", "", 1200, 900);
  canvas.SetTopMargin(0.075);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.04);
  canvas.SetBottomMargin(0.12);
  TH1D frame("frame", "", 100, pt_min, pt_max);
  frame.SetMinimum(0.);
  frame.SetMaximum(y_max);
  frame.GetXaxis()->SetTitle("p_{T}^{truth} [GeV/c]");
  frame.GetYaxis()->SetTitle(y_title);
  frame.GetXaxis()->SetTitleOffset(1.05);
  frame.GetYaxis()->SetTitleOffset(1.10);
  frame.Draw("AXIS");

  TLegend legend(0.55, 0.58, 0.94, 0.91);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.030);
  legend.SetNColumns(2);
  for (size_t index = 0; index < graphs.size(); ++index)
  {
    graphs[index]->Draw("P SAME");
    legend.AddEntry(graphs[index].get(), samples[index].legend.c_str(), "pe");
  }
  legend.Draw();
  DrawHeader();
  canvas.SaveAs(output_path.c_str());
}

}  // namespace

void plot_EMCalPhiResolution_v4(
    int num_files = 10000,
    int first_process = 0,
    bool apply_quality_cuts = true,
    const std::string& eminus_pattern =
        "/sphenix/u/tomoki/work_tomoki/SiCalo/output/"
        "ana_eminus_sicalo_position_v4/SiEmcalTopoTree_e-_%06d_v4.root",
    const std::string& piminus_pattern =
        "/sphenix/u/tomoki/work_tomoki/SiCalo/output/"
        "ana_0803piminus_sicalo_position_v4/SiEmcalTopoTree_pi-_%06d_v4.root",
    const std::string& pp_pattern =
        "/sphenix/u/tomoki/work_tomoki/SiCalo/output/"
        "ana_pythia_sicalo_position_v5/SiEmcalTopoTree_0000000029-%06d_v5.root",
    const std::string& output_directory =
        "/sphenix/u/tomoki/work_tomoki/SiCalo/"
        "analysiscode_bymyself/pdf/emcal_position_resolution_v4",
    double max_abs_dz = 10.0,
    double max_abs_dphi = 0.06,
    double max_track_chi2 = 5000.0,
    double min_calo_chi2 = -100.0,
    double max_calo_chi2 = 48.0)
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gSystem->mkdir(output_directory.c_str(), true);

  const std::vector<double> pt_edges = {
      0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0};
  std::vector<Sample> samples = {
      {"sg_eminus", "single gun e^{-}", 11, kRed + 1, 20, {}},
      {"sg_piminus", "single gun #pi^{-}", 211, kBlue + 1, 21, {}},
      {"pp_e", "p+p e^{#pm}", 11, kMagenta + 2, 22, {}},
      {"pp_pi", "p+p #pi^{#pm}", 211, kCyan + 2, 23, {}},
      {"pp_k", "p+p K^{#pm}", 321, kOrange + 7, 33, {}},
      {"pp_p", "p+p p/#bar{p}", 2212, kGreen + 2, 34, {}},
      {"pp_mu", "p+p #mu^{#pm}", 13, kGray + 2, 29, {}}};
  for (Sample& sample : samples) BookResiduals(sample, pt_edges);

  TChain eminus("SiEmcalTopoTree");
  TChain piminus("SiEmcalTopoTree");
  TChain pp("SiEmcalTopoTree");
  AddFiles(eminus, eminus_pattern, first_process, num_files, "single e-");
  AddFiles(piminus, piminus_pattern, first_process, num_files, "single pi-");
  AddFiles(pp, pp_pattern, first_process, num_files, "p+p");

  if ((eminus.GetEntries() > 0 && !HasRequiredBranches(eminus, "single e-")) ||
      (piminus.GetEntries() > 0 && !HasRequiredBranches(piminus, "single pi-")) ||
      (pp.GetEntries() > 0 && !HasRequiredBranches(pp, "p+p")))
    return;

  std::vector<Sample*> e_targets = {&samples[0]};
  std::vector<Sample*> pi_targets = {&samples[1]};
  std::vector<Sample*> pp_targets = {
      &samples[2], &samples[3], &samples[4], &samples[5], &samples[6]};
  FillResiduals(eminus, e_targets, pt_edges, apply_quality_cuts,
                max_abs_dz, max_abs_dphi, max_track_chi2,
                min_calo_chi2, max_calo_chi2);
  FillResiduals(piminus, pi_targets, pt_edges, apply_quality_cuts,
                max_abs_dz, max_abs_dphi, max_track_chi2,
                min_calo_chi2, max_calo_chi2);
  FillResiduals(pp, pp_targets, pt_edges, apply_quality_cuts,
                max_abs_dz, max_abs_dphi, max_track_chi2,
                min_calo_chi2, max_calo_chi2);

  TFile output(
      (output_directory + "/EMCal_phi_resolution_v4.root").c_str(),
      "RECREATE");
  std::vector<std::unique_ptr<TGraphErrors>> rad_graphs;
  std::vector<std::unique_ptr<TGraphErrors>> cm_graphs;
  for (Sample& sample : samples)
  {
    rad_graphs.emplace_back(
        MakeGraph(sample, pt_edges, 1.0, "rad", output));
    cm_graphs.emplace_back(
        MakeGraph(sample, pt_edges, kEMCalClusterRadius, "cm", output));
  }

  DrawResolution(
      rad_graphs, samples, pt_edges.front(), pt_edges.back(), 0.05,
      "#sigma_{core}(#phi_{cluster}-#phi_{truth entry}) [rad]",
      output_directory + "/EMCal_phi_resolution_vs_truthPt_v4.pdf");
  DrawResolution(
      cm_graphs, samples, pt_edges.front(), pt_edges.back(), 5.2,
      "R_{cluster} #sigma_{core}(#Delta#phi) [cm]",
      output_directory + "/EMCal_phi_resolution_cm_vs_truthPt_v4.pdf");
  output.Close();

  std::cout << "Saved results in " << output_directory << std::endl;
  std::cout << "Note: apply_quality_cuts=" << apply_quality_cuts
            << ". match_dphi is a reco matching cut and can bias the tails; "
               "rerun with false for an untruncated diagnostic." << std::endl;
}

