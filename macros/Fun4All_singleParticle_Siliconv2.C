#ifndef MACRO_FUN4ALL_SINGLEPARTICLE_SILICONV2_C
#define MACRO_FUN4ALL_SINGLEPARTICLE_SILICONV2_C

// pi- single-gun control sample for the Silicon+EMCal pT study.
//
// This macro deliberately performs generation, detector simulation, silicon
// tracking, R=104 cm projection, EMCal matching and pT update in one job.  It
// therefore avoids the ROOT schema/release dependence of rereading an old DST.

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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libsiliconseedsana.so)
R__LOAD_LIBRARY(libtrack_reco.so)

int Fun4All_singleParticle_Siliconv2(
    const int processID = 0,
    const int nEvents = 10,
    const std::string outputBase =
        "/sphenix/u/tomoki/work_tomoki/SiCalo/output",
    const std::string particleName = "anti_proton")
{
  if (processID < 0 || nEvents <= 0)
  {
    std::cerr << "processID must be >= 0 and nEvents must be > 0." << std::endl;
    return 1;
  }

  constexpr float emcalRadius = 104.0F;
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
  else if (particleName == "proton") particleTag = "proton";
  else if (particleName == "anti_proton") particleTag = "antiproton";
  const std::string directoryTag = (particleName == "pi-") ? "0803piminus" : particleTag;

  const std::string process = std::to_string(processID);
  const std::string process6 = Form("%06d", processID);
  const std::string dstDir = outputBase + "/DST_" + directoryTag + "_sicalo_v2";
  const std::string anaDir = outputBase + "/ana_" + directoryTag + "_sicalo_v2";
  gSystem->mkdir(dstDir.c_str(), true);
  gSystem->mkdir(anaDir.c_str(), true);

  const std::string dstOutput =
      dstDir + "/DST_SiliconOnly_single_" + particleTag + "_sicalo_v2_" +
      process + ".root";
  const std::string analyzerOutput =
      anaDir + "/topo_ana_" + process + "_v2.root";
  const std::string flatOutput =
      anaDir + "/SiEmcalTopoTree_" + particleName + "_" + process6 + "_v4.root";

  std::cout << "ana release: use ana.532" << std::endl;
  std::cout << "particle: " << particleName
            << ", pT=0.3--3.0 GeV/c, |eta|<1" << std::endl;
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

  // Build both the standard CEMC clusters and the reference TowerInfo-based
  // CLUSTERINFO_CEMC collection used by SiliconCaloMatching.
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

  // Keep the reference ordering from Fun4All_DST_SiliconSeedAnaKFP.C:
  // charge -> projection -> matching/update -> truth table -> analyzer.
  auto* chargeRecalc = new ChargeRecalc("ChargeRecalc");
  chargeRecalc->setVerbosity(0);
  se->registerSubsystem(chargeRecalc);

  auto* projection = new PHActsTrackProjection("CaloProjection");
  projection->setLayerRadius(SvtxTrack::CEMC, emcalRadius);
  se->registerSubsystem(projection);

  auto* matching = new SiliconCaloMatching("SiliconCaloMatching");
  matching->setEMCalClusterContainerName("CLUSTERINFO_CEMC");
  matching->setEmcalLowEcut(emcalLowEnergyCut);
  matching->setEmcalRadius(emcalRadius);
  matching->setEmcalInOutFrac(1.0F);
  matching->setUpdate_pTinSvtxTrack(true);
  se->registerSubsystem(matching);

  build_truthreco_tables();

  auto* analyzer = new SiliconSeedsAna("SiliconSeedsAna");
  analyzer->setMC(true);
  analyzer->setVtxSkip(true);
  analyzer->setTopoCluster(useTopologicalCluster);
  analyzer->setEMCalClusterContainerName("CLUSTERINFO_CEMC");
  analyzer->setEMcalRadius(emcalRadius);
  analyzer->setStartEventNumber(nEvents * processID);
  analyzer->setOutputFileName(analyzerOutput);
  se->registerSubsystem(analyzer);

  auto* dst = new Fun4AllDstOutputManager("DSTOUT", dstOutput);
  se->registerOutputManager(dst);

  se->run(nEvents);
  se->End();
  delete se;

  // Build the EMCal-center tree, then append the v4 electron-ID features.
  const std::string convertCommand = Form(
      ".x /sphenix/u/tomoki/work_tomoki/SiCalo/macros/"
      "make_SiEmcalTopoTree_v3.C(0,\"%s\",\"%s\",true)",
      analyzerOutput.c_str(), flatOutput.c_str());
  const Long_t convertStatus = gROOT->ProcessLine(convertCommand.c_str());
  if (convertStatus != 0 || gSystem->AccessPathName(flatOutput.c_str()))
  {
    std::cerr << "Flat-tree conversion failed: " << flatOutput << std::endl;
    return 2;
  }

  gROOT->ProcessLine(
      ".L /sphenix/u/tomoki/work_tomoki/SiCalo/macros/"
      "add_SiEmcalTopoTreeFeatures_v4.C");
  const std::string featureCommand = Form(
      "SiEmcalTopoV4::AddV4Branches(\"%s\",\"%s\")",
      analyzerOutput.c_str(), flatOutput.c_str());
  if (gROOT->ProcessLine(featureCommand.c_str()) == 0)
  {
    std::cerr << "Adding v4 electron-ID features failed: "
              << flatOutput << std::endl;
    return 3;
  }

  std::cout << "All done" << std::endl;
  std::cout << "Use plot input pattern: " << anaDir
            << "/SiEmcalTopoTree_" << particleName << "_%06d_v4.root"
            << std::endl;
  return 0;
}

#endif
