#! /usr/bin/env python
import sys,os,shutil,ROOT

# TO BE CUSTOMIZED
#eosmount ~/eos

#prefix = '/afs/cern.ch/user/d/degrutto/eos/cms/'
prefix = ''
#outdir = '/store/group/cmst3/user/degrutto/VHbbV12'
inputdir = '/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_v14/'
#inputdir = '/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Spring15_PU20bx25/skimV12_v2/'
#/store/group/phys_higgs/hbb/ntuples/V12/'
outputdir = '/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_v14/skimmedhaddedv2/'
#/store/group/cmst3/user/degrutto/VHbbV12/hadded/'


filelist = 'listOfZnnFiles.txt'
datasets = os.popen('more '+filelist).read()
datasets = datasets.split('\n')
datasets = [w.replace('\r','') for w in datasets]
datasets = filter(None, datasets) # fastest
# print datasets

count=1

for dataset in datasets:
#    files = os.popen('ls '+ prefix + inputdir +"/"+ dataset).read()
#    files = files.split('\n') 
#    print dataset, files
#    for f in files:
    os.system('hadd  -f '+ prefix+outputdir+'/skim_'+dataset+'.root'+ ' '+prefix+inputdir+dataset+'/skim_'+dataset+'tree*.root')

#    os.system('hadd  -f '+ prefix+outputdir+dataset+"/'/skim_'"'.root'+ ' '+prefix+inputdir+'/skim_'+dataset+'tree*.root')
    

#
#    files = [w.replace('\r','') for w in files]
#    files = filter(None, files) # fastest
#    files =  [elem for elem in files if "/store" in elem]
#    # print dataset
#    # print files
#    for file in files:
#      subdir = file.split('/')[6].replace("VHBB_HEPPY_V12_", "")
#      outfile = file.split('/')[6]+'_'+file.split('/')[9]
#      outfile_reduced = file.split('/')[9]
##      print file.split('/')[9]
##      check = float(os.popen('cmsLs '+outdir+'|wc -l').read())
##      check = float(os.popen('lcg-ls -b -D srmv2 '+outdir+'|wc -l').read())
#      check = float(os.popen('lcg-ls -b -D srmv2 '+'srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=/eos/cms'+outdir+'/'+subdir+'/'+outfile_reduced+' |wc -l').read())
#      print outfile,check,check != 1 
#      if check != 1:
#        os.system('cmsMkdir '+outdir+subdir)  # < -- needed the first time
#        print 'running',os.popen('ps aux | grep $USER | grep lcg-cp |wc -l').read()
#        while float(os.popen('ps aux | grep $USER | grep lcg-cp |wc -l').read())>4: 
#          print os.popen('ps aux | grep $USER | grep lcg-cp |wc -l').read(),' running, waiting 30 seconds'; 
#          os.system('sleep 30')
#        # if(count%6==0): print 'waiting 60 sec'; os.system('sleep 60')
#        print 'copy',outfile
#        # os.system('lcg-cp root://cms-xrd-global.cern.ch/'+file+' '+outdir+'/'+outfile+'&')
#        os.system('lcg-cp --verbose --connect-timeout=300 -b -D srmv2 "srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/'+file+'"'+ ' "srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=/eos/cm#s'+outdir+'/'+subdir+'/'+outfile_reduced+'"'+"&")
##        os.system('echo lcg-cp -b -D srmv2 srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/'+file+' '+outdir+'/'+outfile+'&')
#        count = count+1
#      else: print 'skipping' 
#      # sys.exit()
#    # infile = ROOT.TFile.Open(dataset)
#    # if infile.IsZombie(): 
#        # print dataset,'Is Zombie!'
#        # zombies = zombies+1

print 'finished'
# print 'checked',len(datasets),'datasets,',zombies,'zombies'
