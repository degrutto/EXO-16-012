cd ~/scratch3/ATLASCMSVhbbComboFix/CMSSW_7_1_15/src/
cmsenv
cd -





rm Scaledzp1000_met_SpluB_CC.txt
combineCards.py SR_1=Scaledzp600CCt150250_mass_1J.txt  SR_2=Scaledzp600CCt250350_mass_1J.txt  SR_3=Scaledzp600CCt3509999_mass_1J.txt SR_1_1b=Scaledzp600CCt150250_mass_1b.txt  SR_2_1b=Scaledzp600CCt250350_mass_1b.txt  SR_3_1b=Scaledzp600CCt3509999_mass_1b.txt  SB_1=Scaledzp600CCt150250_mass_SB.txt SB_2=Scaledzp600CCt250350_mass_SB.txt SB_3=Scaledzp600CCt3509999_mass_SB.txt    TTCR=Scaledzp600met1509999_mass_TT.txt   WbCR=Scaledzp600met1509999_mass_Wbb.txt >> Scaledzp1000_met_SpluB_CC.txt
more rateParam.txt >> Scaledzp1000_met_SpluB_CC.txt
combine -M MaxLikelihoodFit Scaledzp600_met_SpluB_CC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0  
combine -M Asymptotic -t -1 Scaledzp1000_met_SpluB_CC.txt
cp higgsCombineTest.Asymptotic.mH120.root Scaledzp1000_met_SpluB_CC_Asymptotic.root







#  LocalWords:  txt
