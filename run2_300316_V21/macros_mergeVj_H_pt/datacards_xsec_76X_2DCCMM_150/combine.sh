combineCards.py SR_1=zp600met150250_met_1J.txt  SR_2=zp600met250350_met_1J.txt  SR_3=zp600met3509999_met_1J.txt   TTCR=zp600met1509999_met_TT.txt WbCR=zp600met1509999_met_Wbb.txt  ZbCR=zp600met1509999_met_SB.txt >> zp600_met_SpluB_2DCC.txt 
more rateParam.txt >> zp600_met_SpluB_2DCC.txt
combine -M MaxLikelihoodFit zp600_met_SpluB_2DCC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0  
combine -M Asymptotic -t -1 zp600_met_SpluB_2DCC.txt





#combineCards.py SR_1=zp1000met150250_met_1J.txt   TTCR_1=zp1000met150250_met_TT.txt WbCR_1=zp1000met150250_met_Wbb.txt  ZbCR_1=zp1000met150250_met_SB.txt SR_2=zp1000met250350_met_1J.txt   TTCR_2=zp1000met250350_met_TT.txt WbCR_2=zp1000met250350_met_Wbb.txt  ZbCR_2=zp1000met250350_met_SB.txt  SR_3=zp1000met3509999_met_1J.txt   TTCR_3=zp1000met3509999_met_TT.txt WbCR_3=zp1000met3509999_met_Wbb.txt  ZbCR_3=zp1000met3509999_met_SB.txt >> zp1000_met_SpluB_2DCC.txt 
#more rateParam.txt >> zp1000_met_SpluB_2DCC.txt
#combine -M MaxLikelihoodFit zp1000_met_SpluB_2DCC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0
#combine -M Asymptotic -t -1 zp1000_met_SpluB_2DCC.txt


#combineCards.py SR_1=zp800met150250_met_1J.txt   TTCR_1=zp800met150250_met_TT.txt WbCR_1=zp800met150250_met_Wbb.txt  ZbCR_1=zp800met150250_met_SB.txt SR_2=zp800met250350_met_1J.txt   TTCR_2=zp800met250350_met_TT.txt WbCR_2=zp800met250350_met_Wbb.txt  ZbCR_2=zp800met250350_met_SB.txt  SR_3=zp800met3509999_met_1J.txt   TTCR_3=zp800met3509999_met_TT.txt WbCR_3=zp800met3509999_met_Wbb.txt  ZbCR_3=zp800met3509999_met_SB.txt >> ! zp800_met_SpluB_2DCC.txt 
#more rateParam.txt >> zp800_met_SpluB_2DCC.txt
#combine -M MaxLikelihoodFit zp800_met_SpluB_2DCC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0
#combine -M Asymptotic -t -1 zp800_met_SpluB_2DCC.txt
