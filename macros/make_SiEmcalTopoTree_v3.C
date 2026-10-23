#include <TFile.h>
#include <TTree.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
constexpr float kMissing = -9999.F;
constexpr float kCaloOffsetX = -0.7F, kCaloOffsetY = -2.3F, kCaloPhiOffset = 0.03F;

// x_proj_emc and x_proj_oemc are the track states at the front (93.5 cm)
// and back (114 cm) of the EMCal.  The EMCal cluster position represents
// the shower near R=104 cm, so interpolate the track state to that radius.
constexpr double kEMCalClusterRadius = 104.0;
constexpr double kEMCalRadiusTolerance = 1.0;
constexpr double kEMCalMaxAbsEta = 1.1232;
constexpr double kEMCalMaxAbsZ = 142.971;
constexpr double kMaxMatchAbsPhi = 0.06;
constexpr double kMaxMatchAbsZ = 10.0;
constexpr int kRejectedOHCalTowerCount = 5;
constexpr int kOHCal = 2;

float FlatTreeWrapPhi(float dphi)
{
  // ROOT emits an error for Phi_mpi_pi(NaN).  Do not let invalid detector
  // values enter the matching calculation.
  if (!std::isfinite(dphi)) return kMissing;
  return TVector2::Phi_mpi_pi(dphi);
}

float FiniteOrMissing(float value)
{
  return std::isfinite(value) ? value : kMissing;
}

// Interpolate existing front/back EMCal track states at R=104 cm.
bool InterpolateProjectionAtRadius(
    const float x_front, const float y_front, const float z_front,
    const float x_back, const float y_back, const float z_back,
    const double target_r,
    float& x_target, float& y_target, float& z_target)
{
  if (!std::isfinite(x_front) || !std::isfinite(y_front) ||
      !std::isfinite(z_front) || !std::isfinite(x_back) ||
      !std::isfinite(y_back) || !std::isfinite(z_back)) return false;

  const double r_front = std::hypot(x_front, y_front);
  const double r_back = std::hypot(x_back, y_back);
  const double dr = r_back - r_front;
  if (!std::isfinite(r_front) || !std::isfinite(r_back) ||
      std::abs(dr) < 1.e-6 ||
      target_r < std::min(r_front, r_back) ||
      target_r > std::max(r_front, r_back)) return false;

  const double fraction = (target_r - r_front) / dr;
  const double phi_front = std::atan2(y_front, x_front);
  const double dphi = TVector2::Phi_mpi_pi(
      std::atan2(y_back, x_back) - phi_front);
  const double phi_target = phi_front + fraction*dphi;
  x_target = static_cast<float>(target_r*std::cos(phi_target));
  y_target = static_cast<float>(target_r*std::sin(phi_target));
  z_target = static_cast<float>(z_front + fraction*(z_back - z_front));
  return std::isfinite(x_target) && std::isfinite(y_target) &&
         std::isfinite(z_target);
}

void GetEmcalClusterXY(const float phi, const bool is_mc,
                       float& x, float& y)
{
  // Simulation uses ideal geometry; offsets are data alignment corrections.
  if (is_mc)
  {
    x = kEMCalClusterRadius*std::cos(phi);
    y = kEMCalClusterRadius*std::sin(phi);
  }
  else
  {
    x = kEMCalClusterRadius*std::cos(phi + kCaloPhiOffset) + kCaloOffsetX;
    y = kEMCalClusterRadius*std::sin(phi + kCaloPhiOffset) + kCaloOffsetY;
  }
}

float CalculateCaloPt(const int charge, const float dphi)
{
  const float x = -charge*dphi;
  const float pt_inv = -0.0002030F + 4.97834F*x + 0.38221F*x*x - 1.27968F*x*x*x;
  return (std::isfinite(pt_inv) && pt_inv > 0.F) ? 1.F/pt_inv : 0.F;
}

// Number of EMCal towers common to one calo cluster and one topo cluster.
int CountSharedEmcalTowers(
    const std::vector<float>& calo_x,
    const std::vector<float>& calo_y,
    const std::vector<float>& calo_z,
    const std::vector<float>& calo_e,
    const std::vector<float>& topo_x,
    const std::vector<float>& topo_y,
    const std::vector<float>& topo_z,
    const std::vector<float>& topo_e,
    const std::vector<int>& topo_caloid,
    const double tol_x, const double tol_y,
    const double tol_z, const double tol_e)
{
  int nshared = 0;
  const size_t ncalo = std::min(std::min(calo_x.size(), calo_y.size()),
                                std::min(calo_z.size(), calo_e.size()));
  const size_t ntopo = std::min(
      std::min(topo_x.size(), topo_y.size()),
      std::min(topo_z.size(), topo_e.size()));

  for (size_t ic = 0; ic < ncalo; ++ic)
  {
    for (size_t it = 0; it < ntopo && it < topo_caloid.size(); ++it)
    {
      if (topo_caloid[it] != 1) continue;  // EMCal only

      if (std::fabs(calo_x[ic] - topo_x[it]) < tol_x &&
          std::fabs(calo_y[ic] - topo_y[it]) < tol_y &&
          std::fabs(calo_z[ic] - topo_z[it]) < tol_z &&
          std::fabs(calo_e[ic] - topo_e[it]) < tol_e)
      {
        ++nshared;
        break;  // one calo tower is counted at most once
      }
    }
  }
  return nshared;
}

// One-to-one calo--topo assignment, starting with the largest shared-tower count.
std::vector<int> MatchCaloToTopo(
    const std::vector<std::vector<float>>* calo_tower_x,
    const std::vector<std::vector<float>>* calo_tower_y,
    const std::vector<std::vector<float>>* calo_tower_z,
    const std::vector<std::vector<float>>* calo_tower_e,
    const std::vector<std::vector<float>>* topo_tower_x,
    const std::vector<std::vector<float>>* topo_tower_y,
    const std::vector<std::vector<float>>* topo_tower_z,
    const std::vector<std::vector<float>>* topo_tower_e,
    const std::vector<std::vector<int>>* topo_tower_caloid,
    const double tol_x, const double tol_y,
    const double tol_z, const double tol_e, const int min_shared_towers)
{
  if (!calo_tower_x || !calo_tower_y || !calo_tower_z || !calo_tower_e ||
      !topo_tower_x || !topo_tower_y || !topo_tower_z || !topo_tower_e ||
      !topo_tower_caloid)
    return {};

  const size_t ncalo = calo_tower_x->size();
  const size_t ntopo = topo_tower_x->size();
  std::vector<int> result(ncalo, -1);
  std::vector<std::vector<int>> shared(ncalo, std::vector<int>(ntopo, 0));

  for (size_t ic = 0; ic < ncalo; ++ic)
  {
    if (ic >= calo_tower_y->size() || ic >= calo_tower_z->size() ||
        ic >= calo_tower_e->size()) continue;
    for (size_t it = 0; it < ntopo; ++it)
    {
      if (it >= topo_tower_y->size() || it >= topo_tower_z->size() ||
          it >= topo_tower_e->size() || it >= topo_tower_caloid->size()) continue;
      shared[ic][it] = CountSharedEmcalTowers(
          calo_tower_x->at(ic), calo_tower_y->at(ic),
          calo_tower_z->at(ic), calo_tower_e->at(ic),
          topo_tower_x->at(it), topo_tower_y->at(it),
          topo_tower_z->at(it), topo_tower_e->at(it), topo_tower_caloid->at(it),
          tol_x, tol_y, tol_z, tol_e);
    }
  }

  std::vector<bool> calo_used(ncalo, false), topo_used(ntopo, false);
  while (true)
  {
    int best_count = 0;
    int best_calo = -1, best_topo = -1;
    for (size_t ic = 0; ic < ncalo; ++ic)
    {
      if (calo_used[ic]) continue;
      for (size_t it = 0; it < ntopo; ++it)
      {
        if (!topo_used[it] && shared[ic][it] > best_count)
        {
          best_count = shared[ic][it];
          best_calo = static_cast<int>(ic);
          best_topo = static_cast<int>(it);
        }
      }
    }
    if (best_count < min_shared_towers) break;
    result[best_calo] = best_topo;
    calo_used[best_calo] = true;
    topo_used[best_topo] = true;
  }
  return result;
}
}  // namespace

// One output entry = one Si-track -> EMCal -> TopoClus match.
// No pair reconstruction and no QA histograms are created.
void make_SiEmcalTopoTree_v3(
    int nevt = 0,
    const std::string& filename="../output/ana_0803piminus/topo_ana_marged.root",
    const std::string& outname="../output/ana_0803piminus/SiEmcalTopoTree_pi-_vscodexv2.root",
    //const std::string& filename="../output/ana_pythia2/ana_0000000029_merged.root",
    //const std::string& outname="../output/ana_pythia2/ana_0000000029_merged_calomatchtowerv4.root",
    bool isMC = true,
    int min_shared_towers = 1,
    double tol_x = 1.e-3, double tol_y = 1.e-3,
    double tol_z = 1.e-3, double tol_e = 1.e-3)
{
  TFile* input = TFile::Open(filename.c_str(), "READ");
  if (!input || input->IsZombie())
  {
    std::cerr << "Cannot open input file: " << filename << std::endl;
    return;
  }

  auto* trackTree = dynamic_cast<TTree*>(input->Get("trackTree"));
  auto* caloTree = dynamic_cast<TTree*>(input->Get("caloTree"));
  auto* siClusTree = dynamic_cast<TTree*>(input->Get("SiClusTree"));
  auto* topoTree = dynamic_cast<TTree*>(input->Get("TopoClusTree"));
  auto* truthTree = isMC ? dynamic_cast<TTree*>(input->Get("truthTree")) : nullptr;
  if (!trackTree || !caloTree || !siClusTree || !topoTree || (isMC && !truthTree))
  {
    std::cerr << "Required tree is missing (trackTree, caloTree, SiClusTree, TopoClusTree"
              << (isMC ? ", truthTree" : "") << ")." << std::endl;
    input->Close();
    return;
  }

  int evt = -1, calo_evt = -1, topo_evt = -1;
  std::vector<float> *track_pt = nullptr, *track_eta = nullptr;
  std::vector<float> *track_pt_proj = nullptr, *track_chi2 = nullptr;
  std::vector<float> *sicalo_pt = nullptr, *sicalo_dphi = nullptr, *sicalo_dz = nullptr;
  std::vector<int> *track_id = nullptr, *track_charge = nullptr;
  std::vector<float> *track_x_emc = nullptr, *track_y_emc = nullptr, *track_z_emc = nullptr;
  std::vector<float> *track_x_oemc = nullptr, *track_y_oemc = nullptr, *track_z_oemc = nullptr;
  std::vector<int>* track_id_truth = nullptr;
  trackTree->SetBranchAddress("evt", &evt);
  trackTree->SetBranchAddress("track_id", &track_id);
  trackTree->SetBranchAddress("charge", &track_charge);
  trackTree->SetBranchAddress("pt0", &track_pt);
  trackTree->SetBranchAddress("eta0", &track_eta);
  trackTree->SetBranchAddress("pt_proj_emc", &track_pt_proj);
  trackTree->SetBranchAddress("chi2ndf", &track_chi2);
  trackTree->SetBranchAddress("sicalo_pt", &sicalo_pt);
  trackTree->SetBranchAddress("sicalo_dphi", &sicalo_dphi);
  trackTree->SetBranchAddress("sicalo_dz", &sicalo_dz);
  trackTree->SetBranchAddress("x_proj_emc", &track_x_emc);
  trackTree->SetBranchAddress("y_proj_emc", &track_y_emc);
  trackTree->SetBranchAddress("z_proj_emc", &track_z_emc);
  trackTree->SetBranchAddress("x_proj_oemc", &track_x_oemc);
  trackTree->SetBranchAddress("y_proj_oemc", &track_y_oemc);
  trackTree->SetBranchAddress("z_proj_oemc", &track_z_oemc);
  if (isMC) trackTree->SetBranchAddress("track_id_truth", &track_id_truth);

  int siclus_evt = -1;
  std::vector<int> *siclus_trackid = nullptr, *siclus_layer = nullptr;
  std::vector<float> *siclus_x = nullptr, *siclus_y = nullptr;
  siClusTree->SetBranchAddress("evt", &siclus_evt);
  siClusTree->SetBranchAddress("Siclus_trackid", &siclus_trackid);
  siClusTree->SetBranchAddress("Siclus_layer", &siclus_layer);
  siClusTree->SetBranchAddress("Siclus_x", &siclus_x);
  siClusTree->SetBranchAddress("Siclus_y", &siclus_y);

  std::vector<float> *calo_phi = nullptr, *calo_eta = nullptr, *calo_r = nullptr;
  std::vector<float> *calo_z = nullptr, *calo_energy = nullptr, *calo_chi2 = nullptr, *calo_prob = nullptr;
  std::vector<std::vector<float>> *calo_tower_e = nullptr, *calo_tower_x = nullptr;
  std::vector<std::vector<float>> *calo_tower_y = nullptr, *calo_tower_z = nullptr;
  caloTree->SetBranchAddress("calo_evt", &calo_evt);
  caloTree->SetBranchAddress("phi", &calo_phi);
  caloTree->SetBranchAddress("eta", &calo_eta);
  caloTree->SetBranchAddress("r", &calo_r);
  caloTree->SetBranchAddress("z", &calo_z);
  caloTree->SetBranchAddress("energy", &calo_energy);
  caloTree->SetBranchAddress("chi2", &calo_chi2);
  caloTree->SetBranchAddress("prob", &calo_prob);
  caloTree->SetBranchAddress("Calo_tower_e", &calo_tower_e);
  caloTree->SetBranchAddress("Calo_tower_x", &calo_tower_x);
  caloTree->SetBranchAddress("Calo_tower_y", &calo_tower_y);
  caloTree->SetBranchAddress("Calo_tower_z", &calo_tower_z);

  // Every TopoClusTree variable used in the original macro is retained below.
  std::vector<float> *clus_e = nullptr, *clus_phi = nullptr, *clus_eta = nullptr, *clus_r = nullptr;
  std::vector<float> *clus_emcal_e = nullptr, *clus_ihcal_e = nullptr, *clus_ohcal_e = nullptr;
  std::vector<std::vector<int>>* tower_caloid = nullptr;
  std::vector<std::vector<float>> *tower_e = nullptr, *tower_x = nullptr, *tower_y = nullptr;
  std::vector<std::vector<float>> *tower_z = nullptr, *tower_phi = nullptr;
  topoTree->SetBranchAddress("evt", &topo_evt);
  topoTree->SetBranchAddress("clus_e", &clus_e);
  topoTree->SetBranchAddress("clus_phi", &clus_phi);
  topoTree->SetBranchAddress("clus_eta", &clus_eta);
  topoTree->SetBranchAddress("clus_r", &clus_r);
  topoTree->SetBranchAddress("clus_emcal_e", &clus_emcal_e);
  topoTree->SetBranchAddress("clus_ihcal_e", &clus_ihcal_e);
  topoTree->SetBranchAddress("clus_ohcal_e", &clus_ohcal_e);
  // tower branches are vector<vector<T>>.  Check first so a schema mismatch
  // is reported cleanly rather than emitting ROOT's "unknown branch" errors.
  if (topoTree->GetBranch("tower_caloid")) topoTree->SetBranchAddress("tower_caloid", &tower_caloid);
  else std::cerr << "TopoClusTree has no branch named tower_caloid" << std::endl;
  if (topoTree->GetBranch("tower_e")) topoTree->SetBranchAddress("tower_e", &tower_e);
  else std::cerr << "TopoClusTree has no branch named tower_e" << std::endl;
  if (topoTree->GetBranch("tower_x")) topoTree->SetBranchAddress("tower_x", &tower_x);
  else std::cerr << "TopoClusTree has no branch named tower_x" << std::endl;
  if (topoTree->GetBranch("tower_y")) topoTree->SetBranchAddress("tower_y", &tower_y);
  else std::cerr << "TopoClusTree has no branch named tower_y" << std::endl;
  if (topoTree->GetBranch("tower_z")) topoTree->SetBranchAddress("tower_z", &tower_z);
  else std::cerr << "TopoClusTree has no branch named tower_z" << std::endl;
  if (topoTree->GetBranch("tower_phi")) topoTree->SetBranchAddress("tower_phi", &tower_phi);
  else std::cerr << "TopoClusTree has no branch named tower_phi" << std::endl;

  std::vector<int> *truth_id = nullptr, *truth_pid = nullptr;
  std::vector<float>* truth_pt = nullptr;
  if (truthTree)
  {
    truthTree->SetBranchAddress("truth_id", &truth_id);
    truthTree->SetBranchAddress("truth_pid", &truth_pid);
    truthTree->SetBranchAddress("truth_pt", &truth_pt);
  }

  TFile output(outname.c_str(), "RECREATE");
  TTree outTree("SiEmcalTopoTree", "Si track matched to EMCal and TopoClus");
  int out_evt, track_index, calo_index, topo_index, out_truth_pid, n_shared_emcal_towers;
  int out_track_charge;
  float out_track_pt, out_track_eta, out_track_chi2, out_track_phi_proj, out_track_z_proj;
  float out_track_r_proj, out_calo_pt_radius;
  float out_truth_pt;
  float out_calo_phi, out_calo_eta, out_calo_r, out_calo_z, out_calo_energy, out_calo_chi2, out_calo_prob;
  float out_match_dphi, out_match_dz, out_match_dr;
  float out_phi_intt, out_phi_calo, out_calo_dphi, out_calo_pt;
  // Preserve the original TopoClusTree cluster branches with their input types.
  std::vector<float> out_clus_e, out_clus_phi, out_clus_eta, out_clus_r;
  std::vector<float> out_clus_emcal_e, out_clus_ihcal_e, out_clus_ohcal_e;
  float out_dphi_calo_topo, out_deta_calo_topo, out_de_calo_topo;
  // Keep the TopoClusTree tower branches without flattening them: every
  // branch remains vector<vector<T>>, including towers not selected by topo_index.
  std::vector<std::vector<int>> out_tower_caloid;
  std::vector<std::vector<float>> out_tower_e, out_tower_x, out_tower_y, out_tower_z, out_tower_phi;

  outTree.Branch("evt", &out_evt);
  outTree.Branch("track_index", &track_index);
  outTree.Branch("calo_index", &calo_index);
  outTree.Branch("topo_index", &topo_index);
  outTree.Branch("track_pt", &out_track_pt);
  outTree.Branch("track_charge", &out_track_charge);
  outTree.Branch("track_eta", &out_track_eta);
  // Additive diagnostic; all vscodexv2 branches are retained unchanged.
  outTree.Branch("track_chi2", &out_track_chi2);
  outTree.Branch("track_phi_proj_emc", &out_track_phi_proj);
  outTree.Branch("track_z_proj_emc", &out_track_z_proj);
  outTree.Branch("track_r_proj_emc", &out_track_r_proj);
  outTree.Branch("truth_pid", &out_truth_pid);
  outTree.Branch("truth_pt", &out_truth_pt);
  outTree.Branch("calo_phi", &out_calo_phi);
  outTree.Branch("calo_eta", &out_calo_eta);
  outTree.Branch("calo_r", &out_calo_r);
  outTree.Branch("calo_z", &out_calo_z);
  outTree.Branch("calo_energy", &out_calo_energy);
  outTree.Branch("calo_chi2", &out_calo_chi2);
  outTree.Branch("calo_prob", &out_calo_prob);
  outTree.Branch("match_dphi", &out_match_dphi);
  outTree.Branch("match_dz", &out_match_dz);
  outTree.Branch("match_dr", &out_match_dr);
  outTree.Branch("phi_intt", &out_phi_intt);
  outTree.Branch("phi_calo", &out_phi_calo);
  outTree.Branch("calo_dphi", &out_calo_dphi);
  outTree.Branch("calo_pt", &out_calo_pt);
  outTree.Branch("calo_pt_radius", &out_calo_pt_radius);
  outTree.Branch("n_shared_emcal_towers", &n_shared_emcal_towers);
  outTree.Branch("clus_e", &out_clus_e);
  outTree.Branch("clus_phi", &out_clus_phi);
  outTree.Branch("clus_eta", &out_clus_eta);
  outTree.Branch("clus_r", &out_clus_r);
  outTree.Branch("clus_emcal_e", &out_clus_emcal_e);
  outTree.Branch("clus_ihcal_e", &out_clus_ihcal_e);
  outTree.Branch("clus_ohcal_e", &out_clus_ohcal_e);
  outTree.Branch("dphi_calo_topo", &out_dphi_calo_topo);
  outTree.Branch("deta_calo_topo", &out_deta_calo_topo);
  outTree.Branch("de_calo_topo", &out_de_calo_topo);
  outTree.Branch("tower_caloid", &out_tower_caloid);
  outTree.Branch("tower_e", &out_tower_e);
  outTree.Branch("tower_x", &out_tower_x);
  outTree.Branch("tower_y", &out_tower_y);
  outTree.Branch("tower_z", &out_tower_z);
  outTree.Branch("tower_phi", &out_tower_phi);

  const Long64_t nentries = std::min(std::min(trackTree->GetEntries(), siClusTree->GetEntries()),
      std::min(caloTree->GetEntries(), topoTree->GetEntries()));

  // TopoClusTree::evt is one-based in some productions.  Determine the
  // constant numbering offset from the first aligned physical event.
  int topo_evt_offset = 0;
  if (nentries > 0)
  {
    trackTree->GetEntry(0);
    siClusTree->GetEntry(0);
    caloTree->GetEntry(0);
    topoTree->GetEntry(0);
    if (evt == siclus_evt && evt == calo_evt)
      topo_evt_offset = topo_evt - evt;
  }
  std::cout << "TopoClusTree evt numbering offset: " << topo_evt_offset
            << std::endl;

  Long64_t rejected_event_mismatch = 0;
  Long64_t rejected_projection = 0;
  Long64_t rejected_match_window = 0;
  Long64_t rejected_topo_energy = 0;
  Long64_t rejected_large_ohcal = 0;
  for (Long64_t ievt = 0; ievt < nentries && (nevt <= 0 || ievt < nevt); ++ievt)
  {
    trackTree->GetEntry(ievt);
    siClusTree->GetEntry(ievt);
    caloTree->GetEntry(ievt);
    topoTree->GetEntry(ievt);
    if (truthTree) truthTree->GetEntry(ievt);
    if (evt != calo_evt || evt != siclus_evt ||
        evt != topo_evt - topo_evt_offset)
    {
      std::cerr << "Skipping entry " << ievt << ": physical evt mismatch "
                << "track=" << evt << ", calo=" << calo_evt
                << ", SiClus=" << siclus_evt << ", Topo="
                << topo_evt - topo_evt_offset << std::endl;
      ++rejected_event_mismatch;
      continue;
    }
    if (!track_id || !track_charge || !track_pt || !track_pt_proj || !track_chi2 || !track_eta ||
        !sicalo_pt || !sicalo_dphi || !sicalo_dz ||
        !track_x_emc || !track_y_emc || !track_z_emc ||
        !track_x_oemc || !track_y_oemc || !track_z_oemc ||
        !siclus_trackid || !siclus_layer || !siclus_x || !siclus_y ||
        !calo_phi || !calo_eta || !calo_r || !calo_z || !calo_energy || !calo_chi2 || !calo_prob ||
        !clus_e || !clus_phi || !clus_eta || !clus_r || !clus_emcal_e || !clus_ihcal_e || !clus_ohcal_e ||
        !tower_caloid || !tower_e || !tower_x || !tower_y || !tower_z || !tower_phi)
      continue;

    std::unordered_map<int, int> truth_index;
    if (truth_id)
      for (size_t i = 0; i < truth_id->size(); ++i)
        truth_index[truth_id->at(i)] = static_cast<int>(i);

    const std::vector<int> calo_to_topo = MatchCaloToTopo(
        calo_tower_x, calo_tower_y, calo_tower_z, calo_tower_e,
        tower_x, tower_y, tower_z, tower_e, tower_caloid,
        tol_x, tol_y, tol_z, tol_e, min_shared_towers);

    const size_t ntracks = std::min(
        std::min(std::min(track_pt_proj->size(), track_chi2->size()),
                 std::min(track_pt->size(), track_eta->size())),
        std::min(std::min(
                     std::min(track_x_emc->size(), std::min(track_y_emc->size(), track_z_emc->size())),
                     std::min(track_x_oemc->size(), std::min(track_y_oemc->size(), track_z_oemc->size()))),
                 std::min(std::min(track_id->size(), track_charge->size()),
                          std::min(sicalo_pt->size(), std::min(sicalo_dphi->size(), sicalo_dz->size())))));
    const size_t ncalo = std::min(calo_phi->size(), std::min(calo_eta->size(),
        std::min(calo_r->size(), std::min(calo_z->size(), std::min(calo_energy->size(),
        std::min(calo_chi2->size(), calo_prob->size()))))));

    for (size_t itrk = 0; itrk < ntracks; ++itrk)
    {
      // PHActsTrackProjection created this state directly at R=104 cm.
      // Do not interpolate between the standard EMCal front/back surfaces.
      const float projected_x = track_x_emc->at(itrk);
      const float projected_y = track_y_emc->at(itrk);
      const float projected_z = track_z_emc->at(itrk);
      if (!std::isfinite(projected_x) || !std::isfinite(projected_y) ||
          !std::isfinite(projected_z) || !std::isfinite(sicalo_pt->at(itrk)) ||
          !std::isfinite(sicalo_dphi->at(itrk)) || !std::isfinite(sicalo_dz->at(itrk)))
      {
        ++rejected_projection;
        continue;
      }
      const double projected_r = std::hypot(projected_x, projected_y);
      const double projected_eta = projected_r > 0.0 ?
          std::asinh(projected_z / projected_r) :
          std::numeric_limits<double>::infinity();
      if (std::abs(projected_r - kEMCalClusterRadius) > kEMCalRadiusTolerance ||
          std::abs(projected_eta) >= kEMCalMaxAbsEta ||
          std::abs(projected_z) >= kEMCalMaxAbsZ)
      {
        ++rejected_projection;
        continue;
      }
      const float phi_proj = std::atan2(projected_y, projected_x);
      const float z_proj = projected_z;
      size_t best_calo = ncalo;
      float best_dphi = kMissing, best_dz = kMissing, best_dr = std::numeric_limits<float>::max();

      for (size_t icalo = 0; icalo < ncalo; ++icalo)
      {
        if (!std::isfinite(calo_phi->at(icalo)) || !std::isfinite(calo_z->at(icalo)) ||
            !std::isfinite(calo_energy->at(icalo)) || calo_energy->at(icalo) < 0.2F) continue;
        const float dphi = FlatTreeWrapPhi(calo_phi->at(icalo) - phi_proj);
        const float dz = calo_z->at(icalo) - z_proj;
        const float dr = std::hypot(
            dphi * static_cast<float>(kEMCalClusterRadius) / 3.F, dz);
        if (std::isfinite(dr) && dr < best_dr)
        {
          best_calo = icalo;
          best_dphi = dphi;
          best_dz = dz;
          best_dr = dr;
        }
      }
      // Keep every association made by SiliconCaloMatching. Analysis cuts on
      // match_dphi and match_dz are deliberately deferred to plotting.
      if (best_calo == ncalo || best_calo >= calo_to_topo.size())
      {
        ++rejected_match_window;
        continue;
      }
      const int matched_topo = calo_to_topo[best_calo];
      if (matched_topo < 0)
      {
        ++rejected_match_window;
        continue;
      }
      const size_t itopo = static_cast<size_t>(matched_topo);
      if (itopo >= clus_e->size() || itopo >= clus_phi->size() || itopo >= clus_eta->size() ||
          itopo >= clus_r->size() || itopo >= clus_emcal_e->size() || itopo >= clus_ihcal_e->size() ||
          itopo >= clus_ohcal_e->size() || itopo >= tower_caloid->size() || itopo >= tower_e->size() ||
          itopo >= tower_x->size() || itopo >= tower_y->size() || itopo >= tower_z->size() || itopo >= tower_phi->size())
        continue;
      if (!(clus_e->at(itopo) > 0.0F))
      {
        ++rejected_topo_energy;
        continue;
      }
      const int ohcal_tower_count = static_cast<int>(std::count(
          tower_caloid->at(itopo).begin(), tower_caloid->at(itopo).end(),
          kOHCal));
      if (ohcal_tower_count >= kRejectedOHCalTowerCount)
      {
        ++rejected_large_ohcal;
        continue;
      }

      out_evt = evt;
      track_index = static_cast<int>(itrk);
      calo_index = static_cast<int>(best_calo);
      topo_index = matched_topo;
      // pt0 is updated by SiliconCaloMatching. The already-created R=104
      // projection state retains the silicon-only pT for comparison.
      out_track_pt = track_pt_proj->at(itrk);
      out_track_chi2 = track_chi2->at(itrk);
      out_track_charge = track_charge->at(itrk);
      out_track_eta = track_eta->at(itrk);
      out_track_phi_proj = phi_proj;
      out_track_z_proj = z_proj;
      out_track_r_proj = static_cast<float>(projected_r);
      out_truth_pid = static_cast<int>(kMissing);
      out_truth_pt = kMissing;
      if (track_id_truth && itrk < track_id_truth->size())
      {
        const auto found = truth_index.find(track_id_truth->at(itrk));
        if (found != truth_index.end())
        {
          const std::size_t truth_position =
              static_cast<std::size_t>(found->second);
          if (truth_pid && truth_position < truth_pid->size())
            out_truth_pid = truth_pid->at(truth_position);
          if (truth_pt && truth_position < truth_pt->size())
            out_truth_pt = truth_pt->at(truth_position);
        }
      }
      out_calo_phi = calo_phi->at(best_calo);
      out_calo_eta = calo_eta->at(best_calo);
      out_calo_r = calo_r->at(best_calo);
      out_calo_z = calo_z->at(best_calo);
      out_calo_energy = calo_energy->at(best_calo);
      out_calo_chi2 = calo_chi2->at(best_calo);
      out_calo_prob = calo_prob->at(best_calo);
      out_match_dphi = sicalo_dphi->at(itrk);
      out_match_dz = sicalo_dz->at(itrk);
      out_match_dr = std::hypot(
          out_match_dphi * static_cast<float>(kEMCalClusterRadius) / 3.F,
          out_match_dz);
      out_phi_intt = kMissing;
      float inner_x = 0.F, inner_y = 0.F, outer_x = 0.F, outer_y = 0.F;
      bool have_inner = false, have_outer = false;
      const size_t nsiclus = std::min(std::min(siclus_trackid->size(), siclus_layer->size()),
                                     std::min(siclus_x->size(), siclus_y->size()));
      for (size_t iclus = 0; iclus < nsiclus; ++iclus)
      {
        if (siclus_trackid->at(iclus) != track_id->at(itrk)) continue;
        if (siclus_layer->at(iclus) < 5)
        { inner_x = siclus_x->at(iclus); inner_y = siclus_y->at(iclus); have_inner = true; }
        else
        { outer_x = siclus_x->at(iclus); outer_y = siclus_y->at(iclus); have_outer = true; }
      }
      if (have_inner && have_outer)
        out_phi_intt = std::atan2(outer_y - inner_y, outer_x - inner_x);
      // pT reconstruction uses the matched caloTree EMCal cluster, never the
      // TopoCluster centroid. Evaluate it at the shower radius R=104 cm;
      // retain calo_r above only as an input diagnostic.
      float emcal_cluster_x = kMissing, emcal_cluster_y = kMissing;
      GetEmcalClusterXY(out_calo_phi, isMC, emcal_cluster_x, emcal_cluster_y);
      out_phi_calo = have_outer ?
          std::atan2(emcal_cluster_y - outer_y, emcal_cluster_x - outer_x) :
          kMissing;
      out_calo_dphi = (have_inner && have_outer) ? FlatTreeWrapPhi(out_phi_calo - out_phi_intt) : kMissing;
      // This is the pT calculated by SiliconCaloMatching in the DST job.
      out_calo_pt = FiniteOrMissing(sicalo_pt->at(itrk));
      out_calo_pt_radius = static_cast<float>(kEMCalClusterRadius);
      out_clus_e = *clus_e;
      out_clus_phi = *clus_phi;
      out_clus_eta = *clus_eta;
      out_clus_r = *clus_r;
      out_clus_emcal_e = *clus_emcal_e;
      out_clus_ihcal_e = *clus_ihcal_e;
      out_clus_ohcal_e = *clus_ohcal_e;
      out_dphi_calo_topo = FlatTreeWrapPhi(out_calo_phi - clus_phi->at(itopo));
      out_deta_calo_topo = FiniteOrMissing(out_calo_eta - clus_eta->at(itopo));
      out_de_calo_topo = FiniteOrMissing(out_calo_energy - clus_e->at(itopo));
      n_shared_emcal_towers = CountSharedEmcalTowers(
          calo_tower_x->at(best_calo), calo_tower_y->at(best_calo),
          calo_tower_z->at(best_calo), calo_tower_e->at(best_calo),
          tower_x->at(itopo), tower_y->at(itopo), tower_z->at(itopo), tower_e->at(itopo),
          tower_caloid->at(itopo), tol_x, tol_y, tol_z, tol_e);
      out_tower_caloid = *tower_caloid;
      out_tower_e = *tower_e;
      out_tower_x = *tower_x;
      out_tower_y = *tower_y;
      out_tower_z = *tower_z;
      out_tower_phi = *tower_phi;
      outTree.Fill();
    }
  }

  output.cd();
  outTree.Write();
  output.Close();
  input->Close();
  std::cout << "Selection summary:"
            << " event_mismatch=" << rejected_event_mismatch
            << ", projection=" << rejected_projection
            << ", match_window_or_no_topo=" << rejected_match_window
            << ", nonpositive_topo_energy=" << rejected_topo_energy
            << ", oHCal_towers_ge_" << kRejectedOHCalTowerCount << "="
            << rejected_large_ohcal << std::endl;
  std::cout << "Saved " << outTree.GetEntries() << " entries to " << outname << std::endl;
}