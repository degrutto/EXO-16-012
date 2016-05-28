#! /usr/bin/env python
import sys,os,shutil,ROOT

# TO BE CUSTOMIZED
#eosmount ~/eos

prefix = ''
#outdir = '/store/group/cmst3/user/degrutto/VHbbV12'
#inputdir = '/store/group/cmst3/user/degrutto/VHbbV12/hadded/'
inputdir = '/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_v13/step3/'






filelist = 'listStep3.txt'
#listOfZnnFiles.txt'
datasets = os.popen('more '+filelist).read()
datasets = datasets.split('\n')
datasets = [w.replace('\r','') for w in datasets]
datasets = filter(None, datasets) # fastest
# print datasets

count=1

#os.system("eosmount ~/eos")

for dataset in datasets:
    print dataset
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\",\\\"%s\\\"\)" % (prefix+inputdir+dataset+".root",dataset+".root"))                                                                               
    os.system("root -b -l -q SkimStep3.C++\(\\\"%s\\\",\\\"%s\\\"\)" % (prefix+inputdir+dataset, dataset  ))                                                                                                     


    
 




    #    last = outname.split()[-1]
#    last = f.rstrip().split("/")[-1]
#    print 'last', last
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\",\\\"%s\\\"\)" % (filename,f.rstrip().split('/')[-1] ))
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\",\\\"%s\\\"\)" % (filename, outname.rstrip()))
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\"\)" % (filename) )
              #   root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\",%i,%i\)
#                print "root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\",%i,%i\)" % (filename,p,method,begin,end)
              
