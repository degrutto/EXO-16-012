#! /usr/bin/env python
import sys,os,shutil,ROOT

# TO BE CUSTOMIZED
#eosmount ~/eos
#outdir = '/store/group/cmst3/user/degrutto/VHbbV12'
#inputdir = '/store/group/cmst3/user/degrutto/VHbbV12/hadded/'
inputdir = '/afs/cern.ch/user/d/degrutto/eos/cms/store/group/cmst3/user/degrutto/ZnnHbbV13/skimmedMet150/'
outputdir = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_v13/"
filelist = 'dataset_V13_skimmed.txt'
#listOfZnnFiles.txt'
datasets = os.popen('more '+filelist).read()
datasets = datasets.split('\n')
datasets = [w.replace('\r','') for w in datasets]
datasets = filter(None, datasets) # fastest
# print datasets

count=1

#os.system("eosmount ~/eos")

for dataset in datasets:
  
  
#    print dataset
    inputname = inputdir + "skim_" + dataset + "*tree*.root"          
    outname = outputdir + "skim_" + dataset + ".root"
#    outname = outname.replace("skimmedMet150", "skimmedMet150hadded")
    print inputname, outname
    os.system("hadd -f %s %s"%(outname, inputname))
 

#    cmd= "bsub -q 8nh -G CMS_CERN01_YODA -J job%s <batchJob_%s"%(outname.replace('.root','').rstrip(),outname.replace('.root',''))
#    open("batchJob_final","a").write(cmd)
 
    
 


#   'bsub -q 8nh  -o std_output.txt -G CMS_CERN01_YODA -J testMDGrelval 

    #    last = outname.split()[-1]
#    last = f.rstrip().split("/")[-1]
#    print 'last', last
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\",\\\"%s\\\"\)" % (filename,f.rstrip().split('/')[-1] ))
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\",\\\"%s\\\"\)" % (filename, outname.rstrip()))
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\"\)" % (filename) )
              #   root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\",%i,%i\)
#                print "root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\",%i,%i\)" % (filename,p,method,begin,end)
              
