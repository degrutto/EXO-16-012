#! /usr/bin/env python
import sys,os,shutil,ROOT

# TO BE CUSTOMIZED
#eosmount ~/eos

prefix = ''
#outdir = '/store/group/cmst3/user/degrutto/VHbbV12'
#inputdir = '/store/group/cmst3/user/degrutto/VHbbV12/hadded/'
#inputdir = '/afs/cern.ch/user/d/degrutto/eos/cms/store/group/cmst3/user/degrutto/ZnnHbbV14/'
inputdir = ''






filelist = 'listMC.txt'
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
    outputname= dataset.replace('/', '')
    os.system("root -b -l -q SkimCopyTree.C++\(\\\"%s\\\",\\\"%s\\\"\)" % (prefix+inputdir+dataset, outputname.replace('/', '')  ))                                                                                                     


    
 




    #    last = outname.split()[-1]
#    last = f.rstrip().split("/")[-1]
#    print 'last', last
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\",\\\"%s\\\"\)" % (filename,f.rstrip().split('/')[-1] ))
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\",\\\"%s\\\"\)" % (filename, outname.rstrip()))
#    os.system("root -b -l -q Skim.C++\(\\\"%s\\\"\)" % (filename) )
              #   root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\",%i,%i\)
#                print "root -b -l -q %s+\(\\\"%s\\\",\\\"%s\\\",%i,%i\)" % (filename,p,method,begin,end)
              
