from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
#config.General.transferLogs = True
config.General.requestName = 'MINIAOD_2022M7_TnP_v3'

config.section_('JobType')
config.JobType.psetName = '../run_muonAnalyzer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['output.root']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')
config.Data.inputDataset = '/ParkingDoubleMuonLowMass7/Run2022F-22Sep2023-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 7
config.Data.splitting = 'LumiBased'
#config.Data.outLFNDirBase = '/store/user/yik/mymultilepNtpl/muonia/2017Bv1Nov3017ReReco'
#config.Data.outLFNDirBase = '/store/group/lpcbphy/hwen/Jpsidata/Charmonium/2018/2018Aver2v1feb2020ReReco'
config.Data.outLFNDirBase = '/store/group/lpcbphy/zhipeng/miniAOD_TnP/MINI2022_v3/2022F/Mass7'
config.Data.lumiMask = './Cert_Collisions2022_355100_362760_Muon.json'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_299614-299617_13TeV_EOY2017ReReco_Collisions17_50ns_JSON_MuonPhys.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_MuonPhys.txt'
config.Data.outputDatasetTag = 'RUN3_22FM7_TnP_v3'

config.section_('User')

#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T2_CH_CERN','T2_US_*','T2_IT_Bari']
#config.Site.blacklist = ['T2_IT_Legnaro']

config.section_('Site')
# config.Site.storageSite = 'cmseos.fnal.gov'
config.Site.storageSite = 'T3_US_FNALLPC'
# config.Site.storageSite = 'T3_CH_CERNBOX' 
#config.Site.whitelist = ['T1_UK_RAL']
