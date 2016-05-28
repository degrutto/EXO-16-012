#!/usr/bin/env python

filename = "./GrowTree_noReg.C"

from ROOT import gROOT

flag = gROOT.ProcessLine('.L %s++O' % filename)
import os
print '-'
print

   
processes = { 



	'ZH_HToBB_ZToNuNu_M125_13TeV_amcatnloFXFX_madspin_pythia8' : 'ZnnH125' ,
	'WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8' :  'WlnH125', 
	'ggZH_HToBB_ZToNuNu_M125_13TeV_amcatnlo_pythia8' : 'ggZH125',
	'WW_TuneCUETP8M1_13TeV-pythia8' : 'WW',
	'WZ_TuneCUETP8M1_13TeV-pythia8' : 'WZ',
	'ZZ_TuneCUETP8M1_13TeV-pythia8' : 'ZZ',

	'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1' : 'T_s_comb_lep',
#	'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1' : 'Tbar_t_lep',
#ZvvHighPt_V21_ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root
	'ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1' : 'T_t_comb_lep',
#ZvvHighPt_V21_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root
	'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1' : 'Tbar_tW',
	'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1' : 'T_tW',
	


              #                    'T_t', 'Tbar_t', 'T_s', 'Tbar_s', 'T_tW', 'Tbar_tW',




#        'ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1' : 'T_t_comb_lep'
#	'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1' : 'Tbar_tW',
#	'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1' : 'T_tW',                      
	#  'WJetsIncl', 
	'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : 'WJetsHT100', 
	'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : 'WJetsHT200', 
	'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : 'WJetsHT400',  
	'WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8': 'WJetsHT600', 
	
	# 'ZJetsHT100', 'ZJetsHT200', 'ZJetsHT400',  'ZJetsHT600', 
	'ZJetsToNuNu_HT-100To200_13TeV-madgraph':  'ZJetsHT100',
	'ZJetsToNuNu_HT-200To400_13TeV-madgraph':  'ZJetsHT200',
	'ZJetsToNuNu_HT-400To600_13TeV-madgraph':  'ZJetsHT400',
	'ZJetsToNuNu_HT-600ToInf_13TeV-madgraph':  'ZJetsHT600',

#	'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8': 'TTMadDiLep' ,
#	'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8': 'TTMadDiLep' ,
	
#	'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8': 'TTMadSLepT' ,
#	'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1': 'TTMadSLepT' ,
	
#	'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM': 'TTMadSLepTbar' ,
#	'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8': 'TTMadSLepTbar' ,
	'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : 'QCDHT100',
	'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8': 'QCDHT200', 
	'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : 'QCDHT300', 
	'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : 'QCDHT500', 
	'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : 'QCDHT700',
	'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : 'QCDHT1000',
	'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : 'QCDHT1500',
	'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : 'QCDHT2000', 
#	'TT_TuneCUETP8M1_13TeV-powheg-pythia8' : 'TTPow',
        'TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8': 'TTMad'

#        'MonoHToBBarMZp-600GeV_MA0-300GeV' : 'monoH',             
#        'MonoHToBBarMZp-800GeV_MA0-300GeV' : 'monoH',             
#        'MonoHToBBarMZp-1000GeV_MA0-300GeV' : 'monoH',             
#        'MonoHToBBarMZp-1200GeV_MA0-300GeV' : 'monoH',             


#	'ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8' : 'ZqqZnunu',

                # fix QCD HT1000
	}

outfiles = []
for p in processes.keys():
        print 'root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\"\)' % (filename,p,processes[p])
#	os.system('root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\"\)' % (filename,p,processes[p]))
        
