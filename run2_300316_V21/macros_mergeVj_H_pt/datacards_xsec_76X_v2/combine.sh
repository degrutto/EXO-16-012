combineCards.py SR=zp800_mt_1J.txt   TTCR=zp800_mt_TT.txt WbCR=zp800_mt_Wbb.txt  ZbCR=zp800_mt_SB.txt  > zp800_mt_SpluB.txt
more rateParam.txt >> zp800_mt_SpluB.txt
combine -M MaxLikelihoodFit zp800_mt_SpluB.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0
combine -M Asymptotic -t -1 zp800_mt_SpluB.txt




combineCards.py SR=zp600_mt_1J.txt   TTCR=zp600_mt_TT.txt WbCR=zp600_mt_Wbb.txt  ZbCR=zp600_mt_SB.txt  > zp600_mt_SpluB.txt
more rateParam.txt >> zp600_mt_SpluB.txt
combine -M MaxLikelihoodFit zp600_mt_SpluB.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0
combine -M Asymptotic -t -1 zp600_mt_SpluB.txt




combineCards.py SR=zp1000_mt_1J.txt   TTCR=zp1000_mt_TT.txt WbCR=zp1000_mt_Wbb.txt  ZbCR=zp1000_mt_SB.txt  > zp1000_mt_SpluB.txt
more rateParam.txt >> zp1000_mt_SpluB.txt
combine -M MaxLikelihoodFit zp1000_mt_SpluB.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0
combine -M Asymptotic -t -1 zp1000_mt_SpluB.txt



