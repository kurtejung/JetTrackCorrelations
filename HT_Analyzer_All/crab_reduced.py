from WMCore.Configuration import Configuration
config = Configuration()

subScript = "HT_Analyzer_All_JFFCorr2_PbPb_Reduced.C"
#subScript = "HT_Analyzer_All_JFFCorr2.C"

config.section_("General")
config.General.requestName = 'pp_5TeVMC_RecoReco_ak4PFJets_taggedB_corrHistos_xiaoLooseppTrkCorrs'
config.General.workArea = config.General.requestName
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
#config.JobType.scriptExe = 'runScript_PbPb.sh'
config.JobType.scriptExe = 'runScript.sh'
config.JobType.scriptArgs = ['script='+subScript]
config.JobType.inputFiles = ['FrameworkJobReport.xml','VzCentReweights.root','VertexCentReweightingFits.root','xiaoTrkCorr.tar','TrkCorr_Jun7_Iterative_PbPb_etaLT2p4.tar.gz','TrkCorr_July22_Iterative_pp_eta2p4.tar.gz','hist_class_def_HT.h',subScript]
config.JobType.outputFiles = ['Pythia_pp.root']
#config.JobType.outputFiles = ['HydJet15_PbPb.root']
#config.JobType.outputFiles = ['Data2015_pp_PbPb.root']
config.JobType.maxJobRuntimeMin = 1800

config.section_("Data")
#config.Data.userInputFiles = open('PbPbData_nCSreskim_July2017.txt').readlines()
#config.Data.userInputFiles = open('pp_Data_5TeV_finalJFFs.txt').readlines()
config.Data.userInputFiles = open('ppMC_Pythia6_MC_finalJFFlist.txt').readlines()
#config.Data.userInputFiles = open('PbPbMC_Pythia6HydjetCymbal_list.txt').readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'pp_MC_Histograms'
config.Data.outLFNDirBase = '/store/user/kjung/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Purdue']
config.Site.storageSite = 'T2_US_Purdue'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

