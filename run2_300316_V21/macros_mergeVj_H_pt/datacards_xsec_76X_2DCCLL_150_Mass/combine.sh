cd ~/scratch3/ATLASCMSVhbbComboFix/CMSSW_7_1_15/src/
cmsenv
cd -
#rm zp600_mass_SpluB_2DCC.txt
#combineCards.py SR_1=zp600met150250_mass_1J.txt  SR_2=zp600met250350_mass_1J.txt  SR_3=zp600met3509999_mass_1J.txt   TTCR=zp600met1509999_mass_TT.txt WbCR=zp600met1509999_mass_Wbb.txt  ZbCR=zp600met1509999_mass_SB.txt >> zp600_mass_SpluB_2DCC.txt 
#more rateParam.txt >> zp600_mass_SpluB_2DCC.txt
#combine -M MaxLikelihoodFit zp600_mass_SpluB_2DCC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0  
#combine -M Asymptotic -t -1 zp600_mass_SpluB_2DCC.txt



#rm zp800_mass_SpluB_2DCC.txt
#combineCards.py SR_1=zp800met150200_mass_1J.txt  SR_2=zp800met200350_mass_1J.txt  SR_3=zp800met3509999_mass_1J.txt   TTCR=zp800met1509999_mass_TT.txt WbCR=zp800met1509999_mass_Wbb.txt >> zp800_mass_SpluB_2DCC.txt 
#more rateParam.txt >> zp800_mass_SpluB_2DCC.txt
#combine -M MaxLikelihoodFit zp800_mass_SpluB_2DCC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0  
#combine -M Asymptotic -t -1 zp800_mass_SpluB_2DCC.txt


rm zp1000_mass_SpluB_2DCC.txt
combineCards.py SR_1=zp1000met150200_mass_1J.txt  SR_2=zp1000met200350_mass_1J.txt  SR_3=zp1000met3509999_mass_1J.txt SR_1_1b=zp1000met150200_mass_1b.txt  SR_2_1b=zp1000met200350_mass_1b.txt  SR_3_1b=zp1000met3509999_mass_1b.txt    TTCR=zp1000met1509999_mass_TT.txt   WbCR=zp1000met1509999_mass_Wbb.txt >> zp1000_mass_SpluB_2DCC.txt
more rateParam.txt >> zp1000_mass_SpluB_2DCC.txt
combine -M MaxLikelihoodFit zp1000_mass_SpluB_2DCC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0  
combine -M Asymptotic -t -1 zp1000_mass_SpluB_2DCC.txt




