
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

namespace SiEmcalTopoV4
{
constexpr float kMissing = -9999.F;
constexpr int kEmcalId = 1;
constexpr float kMinEnergy = 1.e-8F;

float SafeRatio(const float numerator, const float denominator)
{
  if (!std::isfinite(numerator) || !std::isfinite(denominator) || denominator <= kMinEnergy)
    return kMissing;
  return numerator / denominator;
}

float EtaFromXYZ(const float x, const float y, const float z)
{
  const float r = std::hypot(x, y);
  return (std::isfinite(r) && r > 0.F && std::isfinite(z)) ? std::asinh(z / r) : kMissing;
}

int DeltaIPhi(const int lhs, const int rhs)
{
  // sPHENIX CEMC has 256 phi bins.  The second expression handles the
  // boundary between bins 0 and 255.
  const int direct = std::abs(lhs - rhs);
  return std::min(direct, 256 - std::min(direct, 256));
}

bool AddV4Branches(const std::string& input_name, const std::string& output_name)
{
  TFile input(input_name.c_str(), "READ");
  auto* track_tree = dynamic_cast<TTree*>(input.Get("trackTree"));
  auto* topo_tree = dynamic_cast<TTree*>(input.Get("TopoClusTree"));
  if (!track_tree || !topo_tree)
  {
    std::cerr << "add_SiEmcalTopoTreeFeatures_v4: trackTree or TopoClusTree is missing" << std::endl;
    return false;
  }

  const char* required_track[] = {
      "evt", "pt0", "eta0", "phi0", "z0",
      "px_proj_emc", "py_proj_emc", "pz_proj_emc"};
  const char* required_topo[] = {
      "evt", "clus_e", "clus_x", "clus_y", "clus_z", "clus_eta",
      "clus_phi", "clus_prob", "clus_emcal_e", "clus_ihcal_e",
      "clus_ohcal_e", "tower_caloid",
      "tower_ieta", "tower_iphi", "tower_e", "tower_x", "tower_y",
      "tower_z", "tower_phi"};
  for (const char* name : required_track)
    if (!track_tree->GetBranch(name))
    {
      std::cerr << "add_SiEmcalTopoTreeFeatures_v4: trackTree is missing " << name << std::endl;
      return false;
    }
  for (const char* name : required_topo)
    if (!topo_tree->GetBranch(name))
    {
      std::cerr << "add_SiEmcalTopoTreeFeatures_v4: TopoClusTree is missing " << name << std::endl;
      return false;
    }

  int track_evt = -1, topo_evt = -1;
  std::vector<float> *input_track_pt = nullptr, *input_track_eta = nullptr;
  std::vector<float> *input_track_phi = nullptr, *input_track_z = nullptr;
  std::vector<float> *px_emc = nullptr, *py_emc = nullptr, *pz_emc = nullptr;
  std::vector<float> *clus_e = nullptr, *clus_x = nullptr, *clus_y = nullptr, *clus_z = nullptr;
  std::vector<float> *clus_eta_source = nullptr, *clus_phi_source = nullptr;
  std::vector<float> *clus_prob = nullptr, *clus_emcal_e_source = nullptr;
  std::vector<float> *clus_ihcal_e_source = nullptr, *clus_ohcal_e_source = nullptr;
  std::vector<std::vector<int>> *tower_caloid = nullptr, *tower_ieta = nullptr, *tower_iphi = nullptr;
  std::vector<std::vector<float>> *tower_e = nullptr, *tower_x = nullptr, *tower_y = nullptr;
  std::vector<std::vector<float>> *tower_z = nullptr, *tower_phi = nullptr;

  track_tree->SetBranchAddress("evt", &track_evt);
  track_tree->SetBranchAddress("pt0", &input_track_pt);
  track_tree->SetBranchAddress("eta0", &input_track_eta);
  track_tree->SetBranchAddress("phi0", &input_track_phi);
  track_tree->SetBranchAddress("z0", &input_track_z);
  track_tree->SetBranchAddress("px_proj_emc", &px_emc);
  track_tree->SetBranchAddress("py_proj_emc", &py_emc);
  track_tree->SetBranchAddress("pz_proj_emc", &pz_emc);
  topo_tree->SetBranchAddress("evt", &topo_evt);
  topo_tree->SetBranchAddress("clus_e", &clus_e);
  topo_tree->SetBranchAddress("clus_x", &clus_x);
  topo_tree->SetBranchAddress("clus_y", &clus_y);
  topo_tree->SetBranchAddress("clus_z", &clus_z);
  topo_tree->SetBranchAddress("clus_eta", &clus_eta_source);
  topo_tree->SetBranchAddress("clus_phi", &clus_phi_source);
  topo_tree->SetBranchAddress("clus_prob", &clus_prob);
  topo_tree->SetBranchAddress("clus_emcal_e", &clus_emcal_e_source);
  topo_tree->SetBranchAddress("clus_ihcal_e", &clus_ihcal_e_source);
  topo_tree->SetBranchAddress("clus_ohcal_e", &clus_ohcal_e_source);
  topo_tree->SetBranchAddress("tower_caloid", &tower_caloid);
  topo_tree->SetBranchAddress("tower_ieta", &tower_ieta);
  topo_tree->SetBranchAddress("tower_iphi", &tower_iphi);
  topo_tree->SetBranchAddress("tower_e", &tower_e);
  topo_tree->SetBranchAddress("tower_x", &tower_x);
  topo_tree->SetBranchAddress("tower_y", &tower_y);
  topo_tree->SetBranchAddress("tower_z", &tower_z);
  topo_tree->SetBranchAddress("tower_phi", &tower_phi);

  TFile output(output_name.c_str(), "UPDATE");
  auto* tree = dynamic_cast<TTree*>(output.Get("SiEmcalTopoTree"));
  if (!tree)
  {
    std::cerr << "add_SiEmcalTopoTreeFeatures_v4: output SiEmcalTopoTree is missing" << std::endl;
    return false;
  }

  int evt = -1, track_index = -1, topo_index = -1;
  float track_pt = kMissing, track_eta = kMissing;
  std::vector<float> *clus_emcal_e = nullptr, *clus_ihcal_e = nullptr, *clus_ohcal_e = nullptr;
  std::vector<float> *clus_eta = nullptr;
  tree->SetBranchAddress("evt", &evt);
  tree->SetBranchAddress("track_index", &track_index);
  tree->SetBranchAddress("topo_index", &topo_index);
  tree->SetBranchAddress("track_pt", &track_pt);
  tree->SetBranchAddress("track_eta", &track_eta);
  tree->SetBranchAddress("clus_eta", &clus_eta);
  tree->SetBranchAddress("clus_emcal_e", &clus_emcal_e);
  tree->SetBranchAddress("clus_ihcal_e", &clus_ihcal_e);
  tree->SetBranchAddress("clus_ohcal_e", &clus_ohcal_e);

  float track_p, emcal_e_over_p, ihcal_over_emcal, ohcal_over_emcal;
  float incidence_phi, incidence_z;
  float emcal_shower_sigma_eta, emcal_shower_sigma_phi;
  float emcal_max_tower_fraction, emcal_core_energy_fraction;
  float cluster_energy, cluster_et, cluster_eta_scalar, cluster_phi;
  float emcal_energy, ihcal_energy, ohcal_energy, cluster_prob;
  float track_phi0, track_z0;
  int cluster_ntower;
  int n_emcal_towers;

  TBranch* b_track_p = tree->Branch("track_p", &track_p);
  TBranch* b_eop = tree->Branch("emcal_e_over_p", &emcal_e_over_p);
  TBranch* b_ihem = tree->Branch("ihcal_over_emcal", &ihcal_over_emcal);
  TBranch* b_ohem = tree->Branch("ohcal_over_emcal", &ohcal_over_emcal);
  TBranch* b_inphi = tree->Branch("incidence_phi", &incidence_phi);
  TBranch* b_inz = tree->Branch("incidence_z", &incidence_z);
  TBranch* b_seta = tree->Branch("emcal_shower_sigma_eta", &emcal_shower_sigma_eta);
  TBranch* b_sphi = tree->Branch("emcal_shower_sigma_phi", &emcal_shower_sigma_phi);
  TBranch* b_max = tree->Branch("emcal_max_tower_fraction", &emcal_max_tower_fraction);
  TBranch* b_core = tree->Branch("emcal_core_energy_fraction", &emcal_core_energy_fraction);
  TBranch* b_ntower = tree->Branch("n_emcal_towers", &n_emcal_towers);
  TBranch* b_cluster_energy = tree->Branch("cluster_energy", &cluster_energy);
  TBranch* b_cluster_et = tree->Branch("cluster_et", &cluster_et);
  TBranch* b_cluster_eta = tree->Branch("cluster_eta", &cluster_eta_scalar);
  TBranch* b_cluster_phi = tree->Branch("cluster_phi", &cluster_phi);
  TBranch* b_cluster_ntower = tree->Branch("cluster_ntower", &cluster_ntower);
  TBranch* b_emcal_energy = tree->Branch("emcal_energy", &emcal_energy);
  TBranch* b_ihcal_energy = tree->Branch("ihcal_energy", &ihcal_energy);
  TBranch* b_ohcal_energy = tree->Branch("ohcal_energy", &ohcal_energy);
  TBranch* b_cluster_prob = tree->Branch("cluster_prob", &cluster_prob);
  TBranch* b_track_phi0 = tree->Branch("track_phi0", &track_phi0);
  TBranch* b_track_z0 = tree->Branch("track_z0", &track_z0);

  Long64_t source_entry = -1;
  std::set<int> used_track_indices;
  for (Long64_t row = 0; row < tree->GetEntries(); ++row)
  {
    tree->GetEntry(row);
    track_p = emcal_e_over_p = ihcal_over_emcal = ohcal_over_emcal = kMissing;
    incidence_phi = incidence_z = kMissing;
    emcal_shower_sigma_eta = emcal_shower_sigma_phi = kMissing;
    emcal_max_tower_fraction = emcal_core_energy_fraction = kMissing;
    cluster_energy = cluster_et = cluster_eta_scalar = cluster_phi = kMissing;
    emcal_energy = ihcal_energy = ohcal_energy = cluster_prob = kMissing;
    track_phi0 = track_z0 = kMissing;
    cluster_ntower = 0;
    n_emcal_towers = 0;

    if (std::isfinite(track_pt) && std::isfinite(track_eta) && track_pt > 0.F)
      track_p = track_pt * std::cosh(track_eta);

    const size_t itopo = topo_index >= 0 ? static_cast<size_t>(topo_index) : std::numeric_limits<size_t>::max();
    if (clus_emcal_e && clus_ihcal_e && clus_ohcal_e &&
        itopo < clus_emcal_e->size() && itopo < clus_ihcal_e->size() && itopo < clus_ohcal_e->size())
    {
      const float emcal_e = clus_emcal_e->at(itopo);
      emcal_e_over_p = SafeRatio(emcal_e, track_p);
      ihcal_over_emcal = SafeRatio(clus_ihcal_e->at(itopo), emcal_e);
      ohcal_over_emcal = SafeRatio(clus_ohcal_e->at(itopo), emcal_e);
    }

    // hadd inputs can repeat evt numbers. Follow output order and track
    // kinematics to recover the exact source entry.
    Long64_t matched_source_entry = -1;
    for (Long64_t candidate = std::max<Long64_t>(0, source_entry);
         candidate < track_tree->GetEntries(); ++candidate)
    {
      if (candidate != source_entry) used_track_indices.clear();
      track_tree->GetEntry(candidate);
      topo_tree->GetEntry(candidate);
      if (track_evt != evt || track_index < 0 ||
          used_track_indices.count(track_index) || !input_track_eta ||
          static_cast<size_t>(track_index) >= input_track_eta->size())
        continue;
      // pt0 is updated by SiliconCaloMatching, whereas output track_pt is
      // silicon-only. Match using event/index and the unchanged eta.
      if (std::fabs(input_track_eta->at(track_index) - track_eta) > 1.e-5F)
        continue;
      matched_source_entry = candidate;
      source_entry = candidate;
      used_track_indices.insert(track_index);
      break;
    }

    if (matched_source_entry >= 0)
    {
      const bool valid_track_index = track_index >= 0 &&
          input_track_phi && input_track_z &&
          static_cast<size_t>(track_index) < input_track_phi->size() &&
          static_cast<size_t>(track_index) < input_track_z->size();
      if (valid_track_index)
      {
        track_phi0 = input_track_phi->at(track_index);
        track_z0 = input_track_z->at(track_index);
      }

      if (clus_e && clus_eta_source && clus_phi_source && clus_prob &&
          clus_emcal_e_source && clus_ihcal_e_source && clus_ohcal_e_source &&
          tower_e && itopo < clus_e->size() &&
          itopo < clus_eta_source->size() &&
          itopo < clus_phi_source->size() && itopo < clus_prob->size() &&
          itopo < clus_emcal_e_source->size() &&
          itopo < clus_ihcal_e_source->size() &&
          itopo < clus_ohcal_e_source->size() && itopo < tower_e->size())
      {
        cluster_energy = clus_e->at(itopo);
        cluster_eta_scalar = clus_eta_source->at(itopo);
        cluster_phi = clus_phi_source->at(itopo);
        cluster_prob = clus_prob->at(itopo);
        emcal_energy = clus_emcal_e_source->at(itopo);
        ihcal_energy = clus_ihcal_e_source->at(itopo);
        ohcal_energy = clus_ohcal_e_source->at(itopo);
        cluster_ntower = static_cast<int>(tower_e->at(itopo).size());
        if (std::isfinite(cluster_energy) && cluster_energy > 0.F &&
            std::isfinite(cluster_eta_scalar))
          cluster_et = cluster_energy / std::cosh(cluster_eta_scalar);
      }

      if (track_evt == evt && track_index >= 0 &&
          px_emc && py_emc && pz_emc && static_cast<size_t>(track_index) < px_emc->size() &&
          static_cast<size_t>(track_index) < py_emc->size() && static_cast<size_t>(track_index) < pz_emc->size() &&
          clus_x && clus_y && clus_z && itopo < clus_x->size() && itopo < clus_y->size() && itopo < clus_z->size())
      {
        const float px = px_emc->at(track_index), py = py_emc->at(track_index), pz = pz_emc->at(track_index);
        const float pnorm = std::sqrt(px*px + py*py + pz*pz);
        const float x = clus_x->at(itopo), y = clus_y->at(itopo);
        const float r = std::hypot(x, y);
        if (std::isfinite(pnorm) && pnorm > 0.F && std::isfinite(r) && r > 0.F)
        {
          const float pn = (px*x + py*y) / (pnorm*r);
          const float pphi = (-px*y + py*x) / (pnorm*r);
          const float pz_local = pz / pnorm;
          incidence_phi = std::atan2(pphi, pn);
          incidence_z = std::atan2(pz_local, pn);
        }
      }

      if (track_evt == evt && tower_caloid && tower_ieta && tower_iphi &&
          tower_e && tower_x && tower_y && tower_z && tower_phi && clus_eta && clus_x && clus_y &&
          itopo < clus_x->size() && itopo < clus_y->size() && itopo < tower_caloid->size() && itopo < tower_ieta->size() && itopo < tower_iphi->size() &&
          itopo < tower_e->size() && itopo < tower_x->size() && itopo < tower_y->size() &&
          itopo < tower_z->size() && itopo < tower_phi->size() && itopo < clus_eta->size())
      {
        const auto& ids = tower_caloid->at(itopo);
        const auto& ietas = tower_ieta->at(itopo);
        const auto& iphis = tower_iphi->at(itopo);
        const auto& energies = tower_e->at(itopo);
        const auto& xs = tower_x->at(itopo);
        const auto& ys = tower_y->at(itopo);
        const auto& zs = tower_z->at(itopo);
        const auto& phis = tower_phi->at(itopo);
        const size_t ntower = std::min({ids.size(), ietas.size(), iphis.size(), energies.size(),
                                        xs.size(), ys.size(), zs.size(), phis.size()});
        float sum_e = 0.F, sum_deta2 = 0.F, sum_dphi2 = 0.F, max_e = -1.F;
        int max_ieta = 0, max_iphi = 0;
        for (size_t i = 0; i < ntower; ++i)
        {
          if (ids[i] != kEmcalId || !std::isfinite(energies[i]) || energies[i] <= 0.F) continue;
          const float eta = EtaFromXYZ(xs[i], ys[i], zs[i]);
          const float cluster_phi = std::atan2(clus_y->at(itopo), clus_x->at(itopo));
          const float relative_phi = TVector2::Phi_mpi_pi(phis[i] - cluster_phi);
          const float deta = eta - clus_eta->at(itopo);
          if (!std::isfinite(eta) || !std::isfinite(relative_phi) || !std::isfinite(deta)) continue;
          ++n_emcal_towers;
          sum_e += energies[i];
          sum_deta2 += energies[i] * deta * deta;
          sum_dphi2 += energies[i] * relative_phi * relative_phi;
          if (energies[i] > max_e)
          {
            max_e = energies[i];
            max_ieta = ietas[i];
            max_iphi = iphis[i];
          }
        }
        if (sum_e > kMinEnergy && max_e >= 0.F)
        {
          float core_e = 0.F;
          for (size_t i = 0; i < ntower; ++i)
            if (ids[i] == kEmcalId && std::isfinite(energies[i]) && energies[i] > 0.F &&
                std::abs(ietas[i] - max_ieta) <= 1 && DeltaIPhi(iphis[i], max_iphi) <= 1)
              core_e += energies[i];
          emcal_shower_sigma_eta = std::sqrt(std::max(0.F, sum_deta2 / sum_e));
          emcal_shower_sigma_phi = std::sqrt(std::max(0.F, sum_dphi2 / sum_e));
          emcal_max_tower_fraction = max_e / sum_e;
          emcal_core_energy_fraction = core_e / sum_e;
        }
      }
    }

    b_track_p->Fill(); b_eop->Fill(); b_ihem->Fill(); b_ohem->Fill();
    b_inphi->Fill(); b_inz->Fill(); b_seta->Fill(); b_sphi->Fill();
    b_max->Fill(); b_core->Fill(); b_ntower->Fill();
    b_cluster_energy->Fill(); b_cluster_et->Fill();
    b_cluster_eta->Fill(); b_cluster_phi->Fill(); b_cluster_ntower->Fill();
    b_emcal_energy->Fill(); b_ihcal_energy->Fill(); b_ohcal_energy->Fill();
    b_cluster_prob->Fill(); b_track_phi0->Fill(); b_track_z0->Fill();
  }

  tree->Write("", TObject::kOverwrite);
  output.Close();
  input.Close();
  return true;
}
}  // namespace SiEmcalTopoV4

