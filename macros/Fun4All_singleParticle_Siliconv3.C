#ifndef MACRO_FUN4ALL_SINGLEPARTICLE_SILICONV3_C
#define MACRO_FUN4ALL_SINGLEPARTICLE_SILICONV3_C

// pi- single-gun control sample using a charge-compatible constituent tower
// on the CEMC front face as the sixth point of the pT reconstruction.

#include <GlobalVariables.C>

#include <G4Setup_sPHENIX.C>
#include <G4_Global.C>
#include <G4_Input.C>
#include <G4_Mbd.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>

#include <caloreco/RawClusterBuilderTemplate.h>

#include <Trkr_Clustering.C>
#include <Trkr_Reco.C>
#include <Trkr_TruthTables.C>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <siliconseedsana/ChargeRecalc.h>
#include <siliconseedsana/SiliconCaloMatching.h>
#include <siliconseedsana/SiliconSeedsAna.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackreco/PHActsTrackProjection.h>

#include <TROOT.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libsiliconseedsana.so)
R__LOAD_LIBRARY(libtrack_reco.so)

int Fun4All_singleParticle_Siliconv3(
    const int processID = 0,
    const int nEvents = 10,
    const std::string outputBase =
        "/sphenix/u/tomoki/work_tomoki/SiCalo/output",
    const std::string particleName = "e-")
{
  if (processID < 0 || nEvents <= 0)
  {
    std::cerr << "processID must be >= 0 and nEvents must be > 0."
              << std::endl;
    return 1;
  }

  // RawTowerGeom reports CEMC constituent-tower centers on this front face.
  constexpr float emcalFirstTouchRadius = 93.5F;
  constexpr float emcalLowEnergyCut = 0.2F;
  const bool useTopologicalCluster = false;

  // Change only particleName above (or pass it as the fourth argument).
  // These are Geant4 particle names; the tag is used only in path names.
  std::string particleTag = particleName;
  if (particleName == "pi-") particleTag = "piminus";
  else if (particleName == "pi+") particleTag = "piplus";
  else if (particleName == "e-") particleTag = "eminus";
  else if (particleName == "e+") particleTag = "eplus";
  else if (particleName == "mu-") particleTag = "muminus";
  else if (particleName == "mu+") particleTag = "muplus";
  else if (particleName == "kaon-") particleTag = "kminus";
  else if (particleName == "kaon+") particleTag = "kplus";
  else if (particleName == "anti_proton") particleTag = "antiproton";
  const std::string directoryTag = (particleName == "pi-") ? "0803piminus" : particleTag;

  const std::string process = std::to_string(processID);
  const std::string process6 = Form("%06d", processID);
  const std::string dstDir = outputBase + "/DST_" + directoryTag + "_sicalo_v3";
  const std::string anaDir = outputBase + "/ana_" + directoryTag + "_sicalo_v3";
  gSystem->mkdir(dstDir.c_str(), true);
  gSystem->mkdir(anaDir.c_str(), true);

  const std::string dstOutput =
      dstDir + "/DST_SiliconOnly_single_" + particleTag + "_sicalo_v3_" +
      process + ".root";
  const std::string analyzerOutput =
      anaDir + "/topo_ana_" + process + "_v3.root";
  const std::string flatOutput =
      anaDir + "/SiEmcalTopoTree_" + particleName + "_" + process6 + "_v3.root";

  std::cout << "ana release: use ana.532" << std::endl;
  std::cout << "particle: " << particleName
            << ", pT=0.3--3.0 GeV/c, |eta|<1" << std::endl;
  std::cout << "pT sixth point: charge-compatible first-touch CEMC tower"
            << " at R=" << emcalFirstTouchRadius << " cm" << std::endl;
  std::cout << "events: " << nEvents << ", process: " << processID
            << std::endl;
  std::cout << "analyzer output: " << analyzerOutput << std::endl;
  std::cout << "flat-tree output: " << flatOutput << std::endl;

  auto* se = Fun4AllServer::instance();
  se->Verbosity(0);
  PHRandomSeed::Verbosity(0);

  auto* rc = recoConsts::instance();
  Input::VERBOSITY = 0;
  Input::SIMPLE = true;
  InputInit();

  INPUTGENERATOR::SimpleEventGenerator[0]->add_particles(particleName, 1);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(
      PHG4SimpleEventGenerator::Gaus,
      PHG4SimpleEventGenerator::Gaus,
      PHG4SimpleEventGenerator::Gaus);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(
      0.0, 0.0, 0.0);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-1.0, 1.0);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(0.3, 3.0);
  InputRegister();

  if (!Input::READHITS)
  {
    rc->set_IntFlag("RUNNUMBER", 1);
    se->registerSubsystem(new SyncReco());
    se->registerSubsystem(new HeadReco());
  }
  se->registerSubsystem(new FlagHandler());

  Enable::MBDFAKE = true;
  Enable::PIPE = true;
  Enable::PIPE_ABSORBER = true;
  Enable::MVTX = true;
  Enable::INTT = true;
  Enable::TPC = true;
  Enable::MICROMEGAS = true;
  Enable::MAGNET = true;
  Enable::MAGNET_ABSORBER = true;
  Enable::PLUGDOOR_ABSORBER = true;

  Enable::CEMC = true;
  Enable::CEMC_ABSORBER = true;
  Enable::CEMC_CELL = true;
  Enable::CEMC_TOWER = true;
  Enable::CEMC_CLUSTER = true;

  Enable::HCALIN = true;
  Enable::HCALIN_ABSORBER = true;
  Enable::HCALIN_CELL = true;
  Enable::HCALIN_TOWER = true;
  Enable::HCALOUT = true;
  Enable::HCALOUT_ABSORBER = true;
  Enable::HCALOUT_CELL = true;
  Enable::HCALOUT_TOWER = true;

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", CDB::global_tag);
  rc->set_uint64Flag("TIMESTAMP", CDB::timestamp);

  MagnetInit();
  G4Init();
  G4Setup();

  Mbd_Reco();
  Mvtx_Cells();
  Intt_Cells();
  ACTSGEOM::ActsGeomInit();
  Mvtx_Clustering();
  Intt_Clustering();

  CEMC_Cells();
  CEMC_Towers();
  CEMC_Clusters();
  auto* towerInfoCluster =
      new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplateTowerInfo");
  towerInfoCluster->setInputTowerNodeName("TOWERINFO_CALIB_CEMC");
  towerInfoCluster->Detector("CEMC");
  towerInfoCluster->set_threshold_energy(0.070);
  std::string emcalProfile = getenv("CALIBRATIONROOT");
  emcalProfile += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  towerInfoCluster->LoadProfile(emcalProfile);
  towerInfoCluster->set_UseTowerInfo(1);
  towerInfoCluster->set_UseAltZVertex(1);
  se->registerSubsystem(towerInfoCluster);
  HCALInner_Cells();
  HCALOuter_Cells();
  HCALInner_Towers();
  HCALOuter_Towers();
  TopoClusterReco();

  auto* siliconSeeding = new PHActsSiliconSeeding();
  siliconSeeding->Verbosity(0);
  siliconSeeding->setinttZSearchWindow(1.0);
  siliconSeeding->seedAnalysis(false);
  se->registerSubsystem(siliconSeeding);

  auto* merger = new PHSiliconSeedMerger();
  merger->Verbosity(0);
  se->registerSubsystem(merger);

  auto* converter = new TrackSeedTrackMapConverter();
  converter->setTrackSeedName("SiliconTrackSeedContainer");
  converter->setFieldMap(G4MAGNET::magfield_tracking);
  converter->Verbosity(0);
  se->registerSubsystem(converter);

  auto* vertexFinder = new PHSimpleVertexFinder();
  vertexFinder->Verbosity(0);
  vertexFinder->setDcaCut(0.1);
  vertexFinder->setTrackPtCut(0.0);
  vertexFinder->setBeamLineCut(1);
  vertexFinder->setTrackQualityCut(20);
  vertexFinder->setNmvtxRequired(3);
  vertexFinder->setOutlierPairCut(0.1);
  se->registerSubsystem(vertexFinder);
  Global_Reco();

  // The reconstructed charge must be fixed before selecting the EMCal edge.
  auto* chargeRecalc = new ChargeRecalc("ChargeRecalc");
  chargeRecalc->setVerbosity(0);
  se->registerSubsystem(chargeRecalc);

  auto* projection = new PHActsTrackProjection("CaloProjection");
  projection->setLayerRadius(SvtxTrack::CEMC, emcalFirstTouchRadius);
  se->registerSubsystem(projection);

  auto* matching = new SiliconCaloMatching("SiliconCaloMatching");
  matching->setMC(true);
  matching->setEMCalClusterContainerName("CLUSTERINFO_CEMC");
  matching->setEmcalLowEcut(emcalLowEnergyCut);
  matching->setEmcalRadius(emcalFirstTouchRadius);
  matching->setEmcalInOutFrac(1.0F);
  matching->setUseFirstTouchTower(true);
  matching->setUpdate_pTinSvtxTrack(true);
  se->registerSubsystem(matching);

  build_truthreco_tables();

  auto* analyzer = new SiliconSeedsAna("SiliconSeedsAna");
  analyzer->setMC(true);
  analyzer->setVtxSkip(true);
  analyzer->setTopoCluster(useTopologicalCluster);
  analyzer->setEMCalClusterContainerName("CLUSTERINFO_CEMC");
  analyzer->setEMcalRadius(emcalFirstTouchRadius);
  analyzer->setStartEventNumber(nEvents * processID);
  analyzer->setOutputFileName(analyzerOutput);
  se->registerSubsystem(analyzer);

  auto* dst = new Fun4AllDstOutputManager("DSTOUT", dstOutput);
  se->registerOutputManager(dst);

  se->run(nEvents);
  se->End();
  delete se;

  // Keep the flat-tree schema identical to the center-cluster production.
  // calo_pt is the value already written by SiliconCaloMatching; the
  // converter also records the chosen first-touch tower diagnostics.
  const std::string convertCommand = Form(
      ".x /sphenix/u/tomoki/work_tomoki/SiCalo/macros/"
      "make_SiEmcalTopoTree_firsttouch.C(0,\"%s\",\"%s\",true)",
      analyzerOutput.c_str(), flatOutput.c_str());
  const Long_t convertStatus = gROOT->ProcessLine(convertCommand.c_str());
  if (convertStatus != 0 || gSystem->AccessPathName(flatOutput.c_str()))
  {
    std::cerr << "Flat-tree conversion failed: " << flatOutput << std::endl;
    return 2;
  }

  std::cout << "All done" << std::endl;
  std::cout << "Use plot input pattern: " << anaDir
            << "/SiEmcalTopoTree_" << particleName << "_%06d_v3.root"
            << std::endl;
  return 0;
}

#endif
