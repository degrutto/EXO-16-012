combineCards.py SR=zp800_mt_SR.txt   TTCR=zp800_mt_TT.txt WbCR=zp800_mt_Wbb.txt  ZbCR=zp800_mt_Zbb.txt  > zp800_mt_SpluB.txt
more rateParam.txt >> zp800_mt_SpluB.txt
combine -M MaxLikelihoodFit zp800_mt_SpluB.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL --expectSignal=0
combine -M Asymptotic -t -1 zp800_mt_SpluB.txt

