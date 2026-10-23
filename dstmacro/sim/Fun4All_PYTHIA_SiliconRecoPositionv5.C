#ifndef MACRO_FUN4ALL_PYTHIA_SILICONRECOPOSITIONV5_C
#define MACRO_FUN4ALL_PYTHIA_SILICONRECOPOSITIONV5_C

#include <fun4all/Fun4AllUtils.h>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <G4_Mbd.C>
#include <GlobalVariables.C>
#include <QA.C>
#include <G4_TopoClusterReco.C>
#include <G4_RunSettings.C>
#include <Trkr_Clustering.C>
// #include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>
#include <Trkr_TruthTables.C>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>


#include <phool/recoConsts.h>
#include <cdbobjects/CDBTTree.h>
#include <tpccalib/PHTpcResiduals.h>
#include <trackreco/PHActsTrackProjection.h>
#include <trackbase_historic/SvtxTrack.h>
#include <TROOT.h>
#include <kfparticle_sphenix/KFParticle_sPHENIX.h>

#include <tpcqa/TpcRawHitQA.h>
#include <trackingqa/InttClusterQA.h>
#include <trackingqa/MicromegasClusterQA.h>
#include <trackingqa/MvtxClusterQA.h>
#include <trackingqa/TpcClusterQA.h>
#include <trackingqa/TpcSeedsQA.h>
#include <trackingqa/SiliconSeedsQA.h>
#include <trackingqa/VertexQA.h>
#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>
// #include <trackingdiagnostics/KshortReconstruction.h>
#include <simqa_modules/QAG4SimulationTracking.h>


// #include <caloreco/CaloGeomMapping.h>
// #include <caloreco/CaloGeomMappingv2.h>
// #include <caloreco/RawClusterBuilderTemplate.h>
// #include <caloreco/RawClusterBuilderTopo.h>

#include <siliconseedsana/SiliconSeedsAna.h>
#include <siliconseedsana/SiliconCaloMatching.h>
#include <siliconseedsana/ChargeRecalc.h>

#include <stdio.h>
#include <sys/stat.h> // stat, mkdir
#include <limits.h>   // PATH_MAX

#pragma GCC diagnostic push

#pragma GCC diagnostic ignored "-Wundefined-internal"

// #include <kfparticle_sphenix/KshortReconstruction_local.h>

#pragma GCC diagnostic pop

R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtrack_reco.so)
// R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libtpcqa.so)
R__LOAD_LIBRARY(libsiliconseedsana.so)
// For HepMC Hijing
// try inputFile = /sphenix/sim/sim01/sphnxpro/sHijing_HepMC/sHijing_0-12fm.dat


const std::string basedir  = "/sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/";
const int         runnum   = 29;

const std::string outbasedir = "/sphenix/u/tomoki/work_tomoki/SiCalo/output/";
//const std::string outbasedir = "./datafile/";
//
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/calocluster/run0029/detroit/DST_CALO_CLUSTER_pythia8_Detroit-0000000029-000000.root
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/global/run0029/detroit/DST_GLOBAL_pythia8_Detroit-0000000029-000000.root
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/jets/run0029/detroit/DST_TRUTH_JET_pythia8_Detroit-0000000029-000000.root
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/mbdepd/run0029/detroit/DST_MBD_EPD_pythia8_Detroit-0000000029-000000.root
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/trackseeds/run0029/detroit/DST_TRACKSEEDS_pythia8_Detroit-0000000029-000000.root
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/tracks/run0029/detroit/DST_TRACKS_pythia8_Detroit-0000000029-000000.root
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/trkrcluster/run0029/detroit/DST_TRKR_CLUSTER_pythia8_Detroit-0000000029-000000.root
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/trkrhit/run0029/detroit/DST_TRKR_HIT_pythia8_Detroit-0000000029-000000.root
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/trkrhit/run0029/detroit/DST_TRUTH_pythia8_Detroit-0000000029-000000.root
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/nopileup/truthreco/run0029/detroit/DST_TRUTH_RECO_pythia8_Detroit-0000000029-000000.root
///sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/g4hits/run0029/detroit/G4Hits_pythia8_Detroit-0000000029-000000.root
//

int Fun4All_PYTHIA_SiliconRecoPositionv5(
  const int processID = 0,
  const int nEvents   = 10
)
{
  bool useTopologicalCluster = false;

  //--std::ostringstream pid;
  //--pid << std::setw(6) << std::setfill('0') << processID;
  //--std::string pid_str = pid.str();

  // input file
  const std::string trkrclsdir = basedir+"nopileup/trkrcluster/"+Form("run%04d", runnum)+"/detroit/";
  const std::string trkseedsdir= basedir+"nopileup/trackseeds/" +Form("run%04d", runnum)+"/detroit/";
  const std::string trkrdir    = basedir+"nopileup/tracks/"     +Form("run%04d", runnum)+"/detroit/";
  const std::string calodir    = basedir+"nopileup/calocluster/"+Form("run%04d", runnum)+"/detroit/";
  const std::string globaldir  = basedir+"nopileup/global/"     +Form("run%04d", runnum)+"/detroit/";
  const std::string trkrhitdir = basedir+"nopileup/trkrhit/"    +Form("run%04d", runnum)+"/detroit/";
  const std::string hittruthdir= trkrhitdir;
  const std::string truthdir   = basedir+"g4hits/"              +Form("run%04d", runnum)+"/detroit/";

  std::string s_run_pid     = Form("%010d-%06d", runnum, processID);
  std::string trkrclsfilename = "DST_TRKR_CLUSTER_pythia8_Detroit-" + s_run_pid + ".root";
  std::string trkrhitfilename = "DST_TRKR_HIT_pythia8_Detroit-"     + s_run_pid + ".root";
  std::string hittruthfilename= "DST_TRUTH_pythia8_Detroit-"        + s_run_pid + ".root";
  std::string trkseedsfilename= "DST_TRACKSEEDS_pythia8_Detroit-"   + s_run_pid + ".root";
  std::string trksfilename    = "DST_TRACKS_pythia8_Detroit-"       + s_run_pid + ".root";
  std::string calofilename    = "DST_CALO_CLUSTER_pythia8_Detroit-" + s_run_pid + ".root";
  std::string globalfilename  = "DST_GLOBAL_pythia8_Detroit-"       + s_run_pid + ".root";
  std::string truthfilename   = "G4Hits_pythia8_Detroit-"           + s_run_pid + ".root";

  std::string inputtrkrclsFile     = trkrclsdir  + trkrclsfilename;
  std::string inputtrkrhitFile     = trkrhitdir  + trkrhitfilename;
  std::string inputhittruthFile    = hittruthdir + hittruthfilename;
  std::string inputtrkseedsFile    = trkseedsdir + trkseedsfilename;
  std::string inputtrksFile        = trkrdir     + trksfilename;
  std::string inputcaloclusterFile = calodir     + calofilename;
  std::string inputglobalFile      = globaldir   + globalfilename;
  std::string inputtruthFile       = truthdir    + truthfilename;

  std::cout<<"Nevents : "<<nEvents<<std::endl;
  std::cout<<inputtrkrclsFile    <<std::endl;
  std::cout<<inputtrkrhitFile    <<std::endl;
  std::cout<<inputhittruthFile   <<std::endl;
  std::cout<<inputtrkseedsFile   <<std::endl;
  //std::cout<<inputtrksFile       <<std::endl;
  std::cout<<inputcaloclusterFile<<std::endl;
  //std::cout<<inputglobalFile     <<std::endl;
  std::cout<<inputtruthFile      <<std::endl;


  // output dir
  //std::string outDir  = outbasedir + "/dst/condor/";
  //std::string outDir2 = outbasedir + "/ana/condor/";
  //--std::string outDir  = outbasedir + "/dst_pythia1/";
  //--std::string outDir2 = outbasedir + "/ana_pythia1/";
  //std::string outDir  = outbasedir + "/dst_pythia2/"; // with topo
  //std::string outDir2 = outbasedir + "/ana_pythia2/";
  // Do not overwrite the existing production. Both products stay in output/.
  std::string outDir  = outbasedir + "/dst_pythia_sicalo_position_v5/";
  std::string outDir2 = outbasedir + "/ana_pythia_sicalo_position_v5/";
  gSystem->mkdir(outDir.c_str(), true);
  gSystem->mkdir(outDir2.c_str(), true);

  std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(inputtruthFile);
  int runnumber = runseg.first;
  int segment = runseg.second;

  cout<<" run : "<<runnum<<" "<<runnumber<<" "<<segment<<endl;


  RunSettings(runnum);

  
  // Enable::MVTX_APPLYMISALIGNMENT = true;
  // ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;
  // ACTSGEOM::mvtxMisalignment = 100;
  // ACTSGEOM::inttMisalignment = 100.;
  // ACTSGEOM::tpotMisalignment = 100.;
  auto se = Fun4AllServer::instance();
  se->Verbosity(0);

  auto rc = recoConsts::instance();

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  rc->set_uint64Flag("TIMESTAMP", runnum);
  rc->set_IntFlag("RUNNUMBER", runnum);
  //--rc->set_uint64Flag("TIMESTAMP", 6);
  //--rc->set_IntFlag("RUNNUMBER", 1);

  SyncReco *sync = new SyncReco();
  se->registerSubsystem(sync);

  HeadReco *head = new HeadReco();
  se->registerSubsystem(head);
  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  // geometry init done in TrackingInit()
  // std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  // Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  // ingeo->AddFile(geofile);
  // se->registerInputManager(ingeo);

  G4MAGNET::magfield_rescale = 1;
  TrackingInit();

  auto hitsin_truth   = new Fun4AllDstInputManager("DSTin_truth");
  hitsin_truth->fileopen(inputtruthFile);
  se->registerInputManager(hitsin_truth);

  auto hitsin_trkrhit = new Fun4AllDstInputManager("DSTin_trkrhit");
  hitsin_trkrhit->fileopen(inputtrkrhitFile);
  se->registerInputManager(hitsin_trkrhit);

  auto hitsin_hittruth = new Fun4AllDstInputManager("DSTin_hittruth");
  hitsin_hittruth->fileopen(inputhittruthFile);
  se->registerInputManager(hitsin_hittruth);

  auto hitsin_trkrcls = new Fun4AllDstInputManager("DSTin_trkrclus");
  hitsin_trkrcls->fileopen(inputtrkrclsFile);
  se->registerInputManager(hitsin_trkrcls);

  auto hitsin_trkseeds= new Fun4AllDstInputManager("DSTin_seeds");
  hitsin_trkseeds->fileopen(inputtrkseedsFile);
  se->registerInputManager(hitsin_trkseeds);

  auto hitsin_cluster = new Fun4AllDstInputManager("DSTin_cluster");
  // The waveform nodes are unused here and have a newer schema than ana.532.
  // Disable only those branches before fileopen to avoid huge warning logs.
  hitsin_cluster->BranchSelect("DST#CEMC#WAVEFORM_CEMC", 0);
  hitsin_cluster->BranchSelect("DST#HCALIN#WAVEFORM_HCALIN", 0);
  hitsin_cluster->BranchSelect("DST#HCALOUT#WAVEFORM_HCALOUT", 0);
  hitsin_cluster->fileopen(inputcaloclusterFile);
  se->registerInputManager(hitsin_cluster);


  //auto hitsin_trks    = new Fun4AllDstInputManager("DSTin_trks");
  //hitsin_trks->fileopen(inputtrksFile);
  //se->registerInputManager(hitsin_trks);

  //auto hitsin_global  = new Fun4AllDstInputManager("DSTin_global");
  //hitsin_global->fileopen(inputglobalFile);
  //se->registerInputManager(hitsin_global);

  // RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");

  // ClusterBuilder->Detector("CEMC");
  // ClusterBuilder->set_threshold_energy(0.030); // This threshold should be the same as in CEMCprof_Thresh**.root file below
  // std::string emc_prof = getenv("CALIBRATIONROOT");
  // emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  // ClusterBuilder->LoadProfile(emc_prof);
  // ClusterBuilder->set_UseTowerInfo(1); // just use towerinfo
  // //    ClusterBuilder->set_UseTowerInfo(1); // to use towerinfo objects rather than old RawTower
  // se->registerSubsystem(ClusterBuilder);
  TopoClusterReco();

  //--auto silicon_Seeding = new PHActsSiliconSeeding;
  //--// silicon_Seeding->set_track_map_name("SvtxTrackSeedContainer");
  //--silicon_Seeding->Verbosity(0);
  //--// silicon_Seeding->searchInIntt();
  //--// silicon_Seeding->setinttRPhiSearchWindow(0.4);
  //--// silicon_Seeding->setinttZSearchWindow(1.6);
  //--silicon_Seeding->setinttRPhiSearchWindow(0.2);
  //--silicon_Seeding->setinttZSearchWindow(1.0);
  //--silicon_Seeding->seedAnalysis(false);
  //--se->registerSubsystem(silicon_Seeding);

  //--auto merger = new PHSiliconSeedMerger;
  //--// merger->trackMapName("SvtxTrackSeedContainer");
  //--merger->Verbosity(0);
  //--se->registerSubsystem(merger);

  // Default set to full SvtxTrackSeeds. Can be set to
  // SiliconTrackSeedContainer or TpcTrackSeedContainer

  auto converter = new TrackSeedTrackMapConverter("SiliconSeedConverter");
  converter->setTrackSeedName("SiliconTrackSeedContainer");
  // converter->setTrackMapName("SiliconSvtxTrackMap");
  //
  converter->setFieldMap(G4MAGNET::magfield_tracking);
  converter->Verbosity(0);
  se->registerSubsystem(converter);
 


  auto finder = new PHSimpleVertexFinder;
  finder->Verbosity(0);
  finder->setDcaCut(0.075);
  finder->setTrackPtCut(0.1);
  finder->setBeamLineCut(1);
  finder->setTrackQualityCut(300);
  finder->setNmvtxRequired(3);
  finder->setOutlierPairCut(0.075);

  // finder->setTrackMapName("SiliconSvtxTrackMap");
  // finder->setVertexMapName("SiliconSvtxVertexMap");
  se->registerSubsystem(finder);

  // charge recal before projection
  // --- charge recalculation (before projection) ---
  auto chrecalc = new ChargeRecalc("ChargeRecalc");
  chrecalc->setVerbosity(0);
  se->registerSubsystem(chrecalc);
  

  // register GlobalVertexReco
  Global_Reco(); // Definition is below..

  // track projection CALO
  auto projection = new PHActsTrackProjection("CaloProjection");
  constexpr float emcalProjectionRadius = 104.0F;
  projection->setLayerRadius(SvtxTrack::CEMC, emcalProjectionRadius);
  se->registerSubsystem(projection);


  // Match the reference chain: build the high-level truth/reco association
  // before SiliconCaloMatching and every evaluator/analyzer.
  build_truthreco_tables();
  bool runTruth = false;

  //std::string s_pid = Form("%d", processID);
  std::string outputName  = outDir  + "/DST_SiliconOnly_PHYTIA_Info_" + s_run_pid + ".root";
  std::string outputName2 = outDir2 + "/ana_"                         + s_run_pid + ".root";
  std::string outputName3 = outDir2 + "/SiEmcalTopoTree_"             + s_run_pid + "_v5.root";

  std::cout<<outputName <<std::endl;
  std::cout<<outputName2<<std::endl;
  std::cout<<outputName3<<std::endl;


  // SecondaryVertexTagger* myTagger = new SecondaryVertexTagger();
  // myTagger->Verbosity(2);
  // if (runTruth)
  //{
  //   myTagger->truthMode();
  // }
  // else
  //{
  //   myTagger->truthRecoMatch();
  // }
  // myTagger->setOutputName(outputName.c_str());
  // myTagger->saveOutput();
  // se->registerSubsystem(myTagger);
  // outputName = processID + "test.root";

  //--Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", outputName);
  //--se->registerOutputManager(out);

  //--QAG4SimulationTracking *qa = new QAG4SimulationTracking();
  //--qa->Verbosity(0);
  //--  se->registerSubsystem(qa);
  //-- outputName2 = processID + "ana.root";
  //--
 
  // Si-Calo Matching
  SiliconCaloMatching* sicalomatch = new SiliconCaloMatching();
  // Keep the reference default crossing selection. SiliconSeedsAna below
  // remains in MC mode so truth quantities are written to the output tree.
  sicalomatch->setEMCalClusterContainerName("CLUSTERINFO_CEMC");
  sicalomatch->setEmcalLowEcut(0.2);
  sicalomatch->setEmcalRadius(emcalProjectionRadius);
  sicalomatch->setEmcalInOutFrac(1.0F); // reference: inner R=104 cm state only
  sicalomatch->setUpdate_pTinSvtxTrack(true);
  se->registerSubsystem(sicalomatch);
 
  int startnumber = nEvents * processID;
  auto siana = new SiliconSeedsAna("SiliconSeedsAna");
  siana->setMC(true);
  siana->setTrackMapName("SvtxTrackMap");
  siana->setVertexMapName("SvtxVertexMap");
  siana->setTopoCluster(useTopologicalCluster);
  // siana->setVtxSkip(true);
  siana->setOutputFileName(outputName2);
  siana->setStartEventNumber(startnumber);
  siana->setEMCalClusterContainerName("CLUSTERINFO_CEMC");
  siana->setEMcalRadius(emcalProjectionRadius);
  se->registerSubsystem(siana);


  //  se->skip(900);

  se->run(nEvents);

  //std::string qaFile = outDir + "/qa/singleElectron_QA_VTXfixed-" + s_run_pid + ".root";
  //QAHistManagerDef::saveQARootFile(qaFile);

  //-----
  // Exit
  //-----

  //  CDBInterface::instance()->Print(); // print used DB files
  se->End();
  std::cout << "All done" << std::endl;
  delete se;

  // Matching diagnostics are saved without the requested analysis cuts.
  // Loading the converter here (instead of including it above) preserves
  // ROOTs automatic invocation of this Fun4All macro by filename.
  const std::string convertCommand = Form(
      ".x /sphenix/u/tomoki/work_tomoki/SiCalo/macros/make_SiEmcalTopoTree_v4.C(0,\"%s\",\"%s\",true)",
      outputName2.c_str(), outputName3.c_str());
  const Long_t convertStatus = gROOT->ProcessLine(convertCommand.c_str());
  if (convertStatus != 0 || gSystem->AccessPathName(outputName3.c_str()))
  {
    std::cerr << "Base SiEmcalTopoTree conversion failed" << std::endl;
    gSystem->Exit(2);
    return 2;
  }

  gROOT->ProcessLine(
      ".L /sphenix/u/tomoki/work_tomoki/SiCalo/macros/"
      "add_SiEmcalTopoTreeFeatures_v4.C");
  const std::string featureCommand = Form(
      "SiEmcalTopoV4::AddV4Branches(\"%s\",\"%s\")",
      outputName2.c_str(), outputName3.c_str());
  if (gROOT->ProcessLine(featureCommand.c_str()) == 0)
  {
    std::cerr << "Adding v4 SiEmcalTopoTree features failed" << std::endl;
    gSystem->Exit(3);
    return 3;
  }

  gSystem->Exit(0);
  return 0;
}
#endif
