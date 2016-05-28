combineCards.py SR_1=zp600met200300_mt_1J.txt   TTCR_1=zp600met200300_mt_TT.txt WbCR_1=zp600met200300_mt_Wbb.txt  ZbCR_1=zp600met200300_mt_SB.txt SR_2=zp600met300400_mt_1J.txt   TTCR_2=zp600met300400_mt_TT.txt WbCR_2=zp600met300400_mt_Wbb.txt  ZbCR_2=zp600met300400_mt_SB.txt  SR_3=zp600met4009999_mt_1J.txt   TTCR_3=zp600met4009999_mt_TT.txt WbCR_3=zp600met4009999_mt_Wbb.txt  ZbCR_3=zp600met4009999_mt_SB.txt >> ! zp600_mt_SpluB_2DCC.txt 
more rateParam.txt >> zp600_mt_SpluB_2DCC.txt
combine -M MaxLikelihoodFit zp600_mt_SpluB_2DCC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0
combine -M Asymptotic -t -1 zp600_mt_SpluB_2DCC.txt





#combineCards.py SR_1=zp1000met200300_mt_1J.txt   TTCR_1=zp1000met200300_mt_TT.txt WbCR_1=zp1000met200300_mt_Wbb.txt  ZbCR_1=zp1000met200300_mt_SB.txt SR_2=zp1000met300400_mt_1J.txt   TTCR_2=zp1000met300400_mt_TT.txt WbCR_3=zp1000met300400_mt_Wbb.txt  ZbCR_2=zp1000met300400_mt_SB.txt  SR_3=zp1000met4009999_mt_1J.txt   TTCR_3=zp1000met4009999_mt_TT.txt WbCR_3=zp1000met4009999_mt_Wbb.txt  ZbCR_3=zp1000met4009999_mt_SB.txt >> zp1000_mt_SpluB_2DCC.txt 
#more rateParam.txt >> zp1000_mt_SpluB_2DCC.txt
#combine -M MaxLikelihoodFit zp1000_mt_SpluB_2DCC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0
#combine -M Asymptotic -t -1 zp1000_mt_SpluB_2DCC.txt


#combineCards.py SR_1=zp800met200300_mt_1J.txt   TTCR_1=zp800met200300_mt_TT.txt WbCR_1=zp800met200300_mt_Wbb.txt  ZbCR_1=zp800met200300_mt_SB.txt SR_2=zp800met300400_mt_1J.txt   TTCR_2=zp800met300400_mt_TT.txt WbCR_3=zp800met300400_mt_Wbb.txt  ZbCR_2=zp800met300400_mt_SB.txt  SR_3=zp800met4009999_mt_1J.txt   TTCR_3=zp800met4009999_mt_TT.txt WbCR_3=zp800met4009999_mt_Wbb.txt  ZbCR_3=zp800met4009999_mt_SB.txt >> ! zp800_mt_SpluB_2DCC.txt 
#more rateParam.txt >> zp800_mt_SpluB_2DCC.txt
#combine -M MaxLikelihoodFit zp800_mt_SpluB_2DCC.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0
#combine -M Asymptotic -t -1 zp800_mt_SpluB_2DCC.txt
