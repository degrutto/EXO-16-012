cd ~/scratch3/ATLASCMSVhbbComboFix/CMSSW_7_1_15/src/
cmsenv
cd -





rm zp1000_met_SpluB_CC.txt
combineCards.py SR_1=zp600CCt150250_mass_1J.txt  SR_2=zp600CCt250350_mass_1J.txt  SR_3=zp600CCt3509999_mass_1J.txt SR_1_1b=zp600CCt150250_mass_1b.txt  SR_2_1b=zp600CCt250350_mass_1b.txt  SR_3_1b=zp600CCt3509999_mass_1b.txt  SB_1=zp600CCt150250_mass_SB.txt SB_2=zp600CCt250350_mass_SB.txt SB_3=zp600CCt3509999_mass_SB.txt    TTCR=zp600met1509999_mass_TT.txt   WbCR=zp600met1509999_mass_Wbb.txt >> zp1000_met_SpluB_CC.txt
more rateParam.txt >> zp1000_met_SpluB_CC.txt
combine -M MaxLikelihoodFit zp600_met_SpluB_CC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0  
combine -M Asymptotic -t -1 zp1000_met_SpluB_CC.txt
cp higgsCombineTest.Asymptotic.mH120.root zp1000_met_SpluB_CC_Asymptotic.root







#  LocalWords:  txt
