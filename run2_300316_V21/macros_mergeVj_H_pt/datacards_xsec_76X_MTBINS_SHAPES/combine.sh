cd ~/scratch3/ATLASCMSVhbbComboFix/CMSSW_7_1_15/src/
cmsenv
cd -





#rm zp600_mt_SpluB_shape.txt
#combineCards.py SR_1=zp600met150250_mt_1J.txt  SR_2=zp600met250350_mt_1J.txt  SR_3=zp600met3509999_mt_1J.txt SR_1_1b=zp600met150250_mt_1b.txt  SR_2_1b=zp600met250350_mt_1b.txt  SR_3_1b=zp600met3509999_mt_1b.txt  SB_1=zp600met150250_mt_SB.txt SB_2=zp600met250350_mt_SB.txt SB_3=zp600met3509999_mt_SB.txt    TTCR_1=zp600met150250_mt_TT.txt TTCR_2=zp600met250350_mt_TT.txt TTCR_3=zp600met3509999_mt_TT.txt   WbCR_1=zp600met150250_mt_Wbb.txt  WbCR_2=zp600met250350_mt_Wbb.txt  WbCR_3=zp600met3509999_mt_Wbb.txt  >> zp600_mt_SpluB_shape.txt 
#more rateParam.txt >> zp600_mt_SpluB_shape.txt 
#combine -M MaxLikelihoodFit zp600_mt_SpluB_shape.txt  --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0  
#combine -M Asymptotic -t -1 zp600_mt_SpluB_shape.txt 
#cp higgsCombineTest.Asymptotic.mH120.root zp600_mt_SpluB_shape_Asymptotic.root



#rm zp800_mt_SpluB_shape.txt
#combineCards.py SR_1=zp800met150250_mt_1J.txt  SR_2=zp800met250350_mt_1J.txt  SR_3=zp800met3509999_mt_1J.txt SR_1_1b=zp800met150250_mt_1b.txt  SR_2_1b=zp800met250350_mt_1b.txt  SR_3_1b=zp800met3509999_mt_1b.txt  SB_1=zp800met150250_mt_SB.txt SB_2=zp800met250350_mt_SB.txt SB_3=zp800met3509999_mt_SB.txt    TTCR_1=zp800met150250_mt_TT.txt TTCR_2=zp800met250350_mt_TT.txt TTCR_3=zp800met3509999_mt_TT.txt   WbCR_1=zp800met150250_mt_Wbb.txt  WbCR_2=zp800met250350_mt_Wbb.txt  WbCR_3=zp800met3509999_mt_Wbb.txt  >> zp800_mt_SpluB_shape.txt 
#more rateParam.txt >> zp800_mt_SpluB_shape.txt 
#combine -M MaxLikelihoodFit zp800_mt_SpluB_shape.txt  --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0  
#combine -M Asymptotic -t -1 zp800_mt_SpluB_shape.txt 
#cp higgsCombineTest.Asymptotic.mH120.root zp800_mt_SpluB_shape_Asymptotic.root


rm zp1000_mt_SpluB_shape.txt
combineCards.py SR_1=zp1000met150250_mt_1J.txt  SR_2=zp1000met250350_mt_1J.txt  SR_3=zp1000met3509999_mt_1J.txt SR_1_1b=zp1000met150250_mt_1b.txt  SR_2_1b=zp1000met250350_mt_1b.txt  SR_3_1b=zp1000met3509999_mt_1b.txt  SB_1=zp1000met150250_mt_SB.txt SB_2=zp1000met250350_mt_SB.txt SB_3=zp1000met3509999_mt_SB.txt    TTCR_1=zp1000met150250_mt_TT.txt TTCR_2=zp1000met250350_mt_TT.txt TTCR_3=zp1000met3509999_mt_TT.txt   WbCR_1=zp1000met150250_mt_Wbb.txt  WbCR_2=zp1000met250350_mt_Wbb.txt  WbCR_3=zp1000met3509999_mt_Wbb.txt  >> zp1000_mt_SpluB_shape.txt 
more rateParam.txt >> zp1000_mt_SpluB_shape.txt 
combine -M MaxLikelihoodFit zp1000_mt_SpluB_shape.txt  --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0  
combine -M Asymptotic -t -1 zp1000_mt_SpluB_shape.txt 
cp higgsCombineTest.Asymptotic.mH120.root zp1000_mt_SpluB_shape_Asymptotic.root




#  LocalWords:  txt
