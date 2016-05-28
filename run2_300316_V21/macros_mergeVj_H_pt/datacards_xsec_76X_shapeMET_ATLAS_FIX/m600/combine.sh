cd ~/scratch3/ATLASCMSVhbbComboFix/CMSSW_7_1_15/src/
cmsenv
cd -



# zp600metshape1509999_Wbb.txt
# zp600metshape1509999_TT.txt
# zp600metshape1509999_SB.txt
# zp600metshape1509999_1J.txt
# zp600metshape1509999_1b.txt


sed 's/plotsmonoH\_SF\_shapesMET\_ATLAS\///g' zp600metshape1509999_1J.txt > zp600metshape1509999_1J_new.txt
sed 's/plotsmonoH\_SF\_shapesMET\_ATLAS\///g' zp600metshape1509999_Wbb.txt > zp600metshape1509999_Wbb_new.txt 
sed 's/plotsmonoH\_SF\_shapesMET\_ATLAS\///g' zp600metshape1509999_TT.txt > zp600metshape1509999_TT_new.txt 
sed 's/plotsmonoH\_SF\_shapesMET\_ATLAS\///g' zp600metshape1509999_SB.txt > zp600metshape1509999_SB_new.txt
sed 's/plotsmonoH\_SF\_shapesMET\_ATLAS\///g' zp600metshape1509999_1b.txt > zp600metshape1509999_1b_new.txt


rm zp600_metshapes.txt
combineCards.py SR=zp600metshape1509999_1J_new.txt  SR_1b=zp600metshape1509999_1b_new.txt   TTCR=zp600metshape1509999_TT_new.txt   WbCR=zp600metshape1509999_Wbb_new.txt SB=zp600metshape1509999_SB_new.txt>> zp600_metshapes.txt
more rateParam.txt >> zp600_metshapes.txt
combine -M MaxLikelihoodFit zp600_metshapes.txt --saveShapes --saveWithUncertainties -v 3 --forceRecreateNLL -t -1  
combine -M Asymptotic -t -1 zp600_metshapes.txt




