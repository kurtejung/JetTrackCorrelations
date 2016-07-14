from WMCore.Configuration import Configuration
config = Configuration()

subScript = "make_ntuples2.C"

config.section_("General")
config.General.requestName = 'PbPb_5TeV_skimJetTrack'
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'runScript.sh'
config.JobType.scriptArgs = ['script='+subScript]
config.JobType.inputFiles = ['FrameworkJobReport.xml','class_def/GenParticles.h','class_def/run2/HiTree.h','class_def/run2/HLT.h','class_def/run2/JetAna.h','class_def/run2/pfcand.h','class_def/run2/Skim.h','class_def/run2/Tracks.h','fragmentation_JEC.h',subScript]
config.JobType.outputFiles = ['Data2015.root']

config.section_("Data")
config.Data.userInputFiles = open('PbPbForest_HiRun2015_MITForest.txt').readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'PbPb_5TeV_JetTrackCorrCrabSkim'
config.Data.outLFNDirBase = '/store/user/kjung/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_MIT']
config.Site.storageSite = 'T2_US_Purdue'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

