import os

def Process(file, value ):
       newcard = "Scaled%s" % file
       os.system('more %s | grep "rate " > a' % file)
       os.system("more a |  awk '{print $2}'  > aMonoH")
 
       fMonoH = open("aMonoH", "r")

       sMonoH = fMonoH.read().rstrip('\n')

       MonoH = float(sMonoH)

       bMonoH = (value)*MonoH

       os.system("sed 's/%s/%f/g' a > b" % (sMonoH,bMonoH))

       fRO = open("a", "r")
       sRO = fRO.read().rstrip('\n')
       fRF = open("b", "r")
       sRF = fRF.read().rstrip('\n')
       os.system("sed 's/%s/%s/g' %s > %s" % (sRO,sRF,file,newcard))


for item in [
'zp600met1509999_mass_Zlight.txt',
'zp600met1509999_mass_Wlight.txt',
'zp600met1509999_mass_Wbb.txt',
'zp600met1509999_mass_TT.txt',
'zp600met1509999_mass_SB.txt',
'zp600met1509999_mass_1b.txt',
'zp600met1509999_mass_1J.txt',
'zp600met1509999_mass_0J.txt',
'zp600CCt3509999_mass_Zlight.txt',
'zp600CCt3509999_mass_Wlight.txt',
'zp600CCt3509999_mass_Wbb.txt',
'zp600CCt3509999_mass_TT.txt',
'zp600CCt3509999_mass_SB.txt',
'zp600CCt3509999_mass_1b.txt',
'zp600CCt3509999_mass_1J.txt',
'zp600CCt3509999_mass_0J.txt',
'zp600CCt250350_mass_Zlight.txt',
'zp600CCt250350_mass_Wlight.txt',
'zp600CCt250350_mass_Wbb.txt',
'zp600CCt250350_mass_TT.txt',
'zp600CCt250350_mass_SB.txt',
'zp600CCt250350_mass_1b.txt',
'zp600CCt250350_mass_1J.txt',
'zp600CCt250350_mass_0J.txt',
'zp600CCt150250_mass_Zlight.txt',
'zp600CCt150250_mass_Wlight.txt',
'zp600CCt150250_mass_Wbb.txt',
'zp600CCt150250_mass_TT.txt',
'zp600CCt150250_mass_SB.txt',
'zp600CCt150250_mass_1b.txt',
'zp600CCt150250_mass_1J.txt',
'zp600CCt150250_mass_0J.txt',]:



# scale for new CMS/ATLAS numbers                                                                                                                                                           
#     hboostedZH->Scale( 42.386 * 0.577 / 0.026 );                                                                                                                                         
#      hboostedZH->Scale( 372.2 * 0.577 / 0.026 );                                                                                                                                          

##      hboostedZH_2->Scale( 45.097 * 0.577 / 0.0288 );                                                                                                                                      
##hboostedZH_2->Scale( 230.67 * 0.577 / 0.0288 );                                                                                                                                            


#hboostedZH_2->Scale( 35.444 * 0.577 / 0.02337 );                                                                                                                                     
#hboostedZH_2->Scale( 119.34 * 0.577 / 0.02337 );                                                                                                                                           




  Process(item, 0.3722 * 0.577 / 0.026)
