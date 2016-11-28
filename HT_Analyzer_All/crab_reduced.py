from WMCore.Configuration import Configuration
config = Configuration()

subScript = "HT_Analyzer_All_JFFCorr2_PbPb_Reduced.C"
#subScript = "HT_Analyzer_All_JFFCorr2.C"

config.section_("General")
config.General.requestName = 'PbPb_Pythia6Hydjet_10EvtMixed_histoFiles_Reduced_RecoJetRecoTrack'
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'runScript_PbPb.sh'
#config.JobType.scriptExe = 'runScript.sh'
config.JobType.scriptArgs = ['script='+subScript]
config.JobType.inputFiles = ['FrameworkJobReport.xml','VzCentReweights.root','VertexCentReweightingFits.root','TrkCorr_Jun7_Iterative_PbPb_etaLT2p4.tar.gz','hist_class_def_HT.h',subScript]
#config.JobType.inputFiles = ['FrameworkJobReport.xml','hist_class_def_HT.h',subScript,'TrkCorr_July22_Iterative_pp_eta2p4.tar.gz','PbPb_JetSpectra_JFFCorr.root','pp_JetSpectra_JFFCorr.root']
config.JobType.outputFiles = ['HydJet15_PbPb.root']
config.JobType.maxJobRuntimeMin = 1800

config.section_("Data")
#config.Data.userInputFiles = open('PbPb_MC_Pythia6Hydjet_5TeV.txt').readlines()
config.Data.userInputFiles = open('PbPb_MC_Pythia6Hydjet_5TeV.txt').readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'PbPb_MC_Histograms'
config.Data.outLFNDirBase = '/store/user/kjung/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Purdue']
config.Site.storageSite = 'T2_US_Purdue'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

