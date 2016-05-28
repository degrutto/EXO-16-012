#!/usr/bin/env python

filename = "./GrowTree_noReg.C"

from ROOT import gROOT

flag = gROOT.ProcessLine('.L %s++O' % filename)
import os
print '-'
print

   
processes = { 



'ZH_HToBB_ZToNuNu_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'ZnnH125' ,
'WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' :  'WlnH125', 
     'ggZH_HToBB_ZToNuNu_M125_13TeV_amcatnlo_pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'ggZH125',
	      'WW_TuneCUETP8M1_13TeV-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'WW',
	      'WZ_TuneCUETP8M1_13TeV-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'WZ',
              'ZZ_TuneCUETP8M1_13TeV-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v3' : 'ZZ',
'ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'ZqqZnunu',
'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'T_s_comb_lep',
'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'Tbar_t_lep',
'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'T_t_lep',
'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'Tbar_tW',
'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'T_tW',



              #                    'T_t', 'Tbar_t', 'T_s', 'Tbar_s', 'T_tW', 'Tbar_tW',



              'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'Tbar_tW',
	      'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'T_tW',                      
              #  'WJetsIncl', 
              'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'WJetsHT100', 
              'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'WJetsHT200', 
              'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v3' : 'WJetsHT400',  
              'WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1': 'WJetsHT600', 
              
              # 'ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400',  'ZJetsHT600', 
              'ZJetsToNuNu_HT-100To200_13TeV-madgraph__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1':  'ZJetsHT100',
              'ZJetsToNuNu_HT-200To400_13TeV-madgraph__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1':  'ZJetsHT200',
              'ZJetsToNuNu_HT-400To600_13TeV-madgraph__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1':  'ZJetsHT400',
              'ZJetsToNuNu_HT-600ToInf_13TeV-madgraph__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1':  'ZJetsHT600',
	      'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2': 'TTMadDiLep' ,
	      'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-ext1-v1': 'TTMadDiLep' ,

	      'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1': 'TTMadSLepT' ,
	      'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_ext1-v1': 'TTMadSLepT' ,

	      'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2': 'TTMadSLepTbar' ,
'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-ext1-v1': 'TTMadSLepTbar' ,
              'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2' : 'QCDHT100',
              'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2': 'QCDHT200', 
              'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2' : 'QCDHT300', 
              'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'QCDHT500', 
              'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'QCDHT700',
              'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'QCDHT1000',
              'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'QCDHT1500',
              'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1' : 'QCDHT2000',
          'TT_TuneCUETP8M1_13TeV-powheg-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2' : 'TTPow'             
                # fix QCD HT1000
              }

outfiles = []
for p in processes.keys():
        print 'root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\"\)' % (filename,p,processes[p])
#	os.system('root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\"\)' % (filename,p,processes[p]))
        
