import sys
import os


NewAltCodes = {"UsingMueller.html" :
               {b"_images/UsingMueller_20_0.svg" :
                b"Noisy curve starting at 1 and decaying to a level of between about 0.1 and 0.2.",
                b"_images/UsingMueller_26_0.svg" :
                b"Four curves starting, one starting at 1 and decaying to about 0.35 and three"
                b"starting at 0 and rising to levels around about 0.25."
               },
               "UsingFresnel.html" :
               {b"_images/UsingFresnel_4_0.svg" :
                b"Optical constants for several materials as functions of wavelength",
                b"_images/UsingFresnel_12_0.svg" :
                b"Reflectance of silicon as a function of angle",
                b"_images/UsingFresnel_14_0.svg" :
                b"Reflectance of silicon as a function of wavelength",
                b"_images/UsingFresnel_16_0.svg" :
                b"Ellipsometric parameters for silicon as a function of wavelength",
                b"_images/UsingFresnel_18_0.svg" :
                b"Ellipsometric parameters for films on silicon as a function of wavelength",
                b"_images/UsingFresnel_20_0.svg" :
                b"Reflectance and transmittance for films on silicon as a function of angle",
                b"_images/UsingFresnel_22_0.svg" :
                b"Reflectance for different Bragg reflectors as a function of wavelength",
                b"_images/UsingFresnel_24_0.svg" :
                b"Reflectance for a notch filter as a function of wavelength"
               },
               "UsingBRDF.html" :
               {b"_images/UsingBRDF_21_0.svg" :
                b"BRDF as a function of angle with a sharp peak at 60 degrees",
                b"_images/UsingBRDF_23_0.svg" :
                b"Three BRDFs as a function of angle with sharp peaks at 0, 30, and 60 degrees",
                b"_images/UsingBRDF_25_0.svg" :
                b"Five BRDFs as a function of angle with sharp peaks at 60 degrees"
               },
               "UsingLocal.html" :
               {b"_images/UsingLocal_8_0.svg" :
                b"Differential scattering cross section as a function of angle with two broad peaks",
                b"_images/UsingLocal_9_0.svg" :
                b"Four differential scattering cross sections as a function of angle"
               },
               "UsingIntegrate.html" :
               {b"_images/UsingIntegrate_7_0.svg" :
                b"The unit circle with unevenly spaced points inside",
                b"_images/UsingIntegrate_7_1.svg" :
                b"The unit circle with unevenly spaced points inside",
                b"_images/UsingIntegrate_7_2.svg" :
                b"The unit circle with evenly spaced points inside",
                b"_images/UsingIntegrate_9_0.svg" :
                b"The unit circle with points in an ellipse inside of it",
                b"_images/UsingIntegrate_11_0.svg" :
                b"The unit circle with points in an ellipse inside of it",
                b"_images/UsingIntegrate_13_0.svg" :
                b"The unit circle with points in a tilted ellipse inside of it",
                b"_images/UsingIntegrate_15_0.svg" :
                b"The unit circle with points in a rectangle inside of it",
                b"_images/UsingIntegrate_17_0.svg" :
                b"The unit circle with points in a tilted rectangle inside of it",
                b"_images/UsingIntegrate_19_0.svg" :
                b"The unit circle with points in a bow tie inside of it",
                b"_images/UsingIntegrate_21_0.svg" :
                b"The unit circle with points filling it, except in two small circles",
                b"_images/UsingIntegrate_23_0.svg" :
                b"The unit circle with points filling it, except in two small notches and an annulus",
                b"_images/UsingIntegrate_25_0.svg" :
                b"The unit circle with points filling it, except along a horizontal strip",
                b"_images/UsingIntegrate_27_0.svg" :
                b"The unit circle with points filling it, except for a central circle",
                b"_images/UsingIntegrate_29_0.svg" :
                b"A straight line",
                b"_images/UsingIntegrate_31_0.svg" :
                b"The unit circle with points filling it, except for a central circle",
                b"_images/UsingIntegrate_33_0.svg" :
                b"Two curves",
                b"_images/UsingIntegrate_35_0.svg" :
                b"The unit circle with a complicated color mapping on it",
                b"_images/UsingIntegrate_37_0.svg" :
                b"Single curve representing the cross section versus diameter."},
               "UsingScatterer.html" :
               {b"_images/UsingScatterer_8_0.svg" :
                b"Polar graph with three curves",
                b"_images/UsingScatterer_10_0.svg" :
                b"A peak around 0.520 micrometers",
                b"_images/UsingScatterer_13_0.svg" :
                b"Three graphs, each with seven curves, showing scattering, extinction, and absorption."
               },
                "UsingRCW.html" :
               {b"_images/UsingRCW_6_0.svg" :
                b"Three complicated curves labeled N, C, and S as a function of wavelength.",
                b"_images/UsingRCW_8_0.svg" :
                b"Four by four array of panes, each with one curve.",
                b"_images/UsingRCW_10_0.svg" :
                b"Cross section of a trench profile.",
                b"_images/UsingRCW_12_0.svg" :
                b"The unit circle with circles of different sizes lying along a straight line.",
                b"_images/UsingRCW_14_0.svg" :
                b"Two curves",
                b"_images/UsingRCW_14_1.svg" :
                b"Cross section of a blazed grating.",
                b"_images/UsingRCW_16_0.svg" :
                b"Two curves",
                b"_images/UsingRCW_16_1.svg" :
                b"Cross section of a sinusoidal grating."
               },
                "UsingCrossRCW.html" :
               {b"_images/UsingCrossRCW_6_0.svg" :
                b"Three complicated curves labeled N, C, and S as a function of wavelength.",
                b"_images/UsingCrossRCW_8_0.svg" :
                b"Four by four array of panes, each with one curve.",
                b"_images/UsingCrossRCW_10_0.svg" :
                b"The unit circle with circles of different sizes on hexagonal grid.",
               }
              }                

os.system("make clean")
os.system("make html")

directory = 'build/html/'
files = [directory+f for f in os.listdir(directory) if ".html" in f]

prefix = b"Graph showing results of preceding Python code: "

def printRed(*skk) : print("\033[91m",*skk,"\033[00m")

def doReplacement(contents,filename):
    """
    Look for pattern '<img alt="_images/* "' and replaces 
    it with 
    <img alt="Graph showing results of preceding Python code"
    """
    n2=0
    while True:
        n1 = contents.find(b'<img alt="_images/',n2)
        if n1 == -1: break
        n2 = contents.find(b'" ',n1)
        original = contents[n1:n2+2]
        tag = contents[n1+10:n2]
        try:
            print(filename[len(directory):])
            newAltCode = NewAltCodes[filename[len(directory):]][tag]
            contents = contents.replace(contents[n1:n2+2], b'<img alt="' + prefix + newAltCode + b'" ')
            print(original,"replaced at ",n1,":",newAltCode)
        except:
            printRed(original,"not replaced at ",n1," in ",filename)
            
    return contents

for f in files:
    contents = open(f,"rb").read() 
    contents = doReplacement(contents,f)
    open(f,"wb").write(contents)


#os.system(r"robocopy . \\elwood\68_pml\685\internal\Public\Germer\pySCATMECH * /e")


SphinxConfigFile = "source/conf.py"
InstallSetupFile = "../setup.py"

print("Checking version consistency:")
contents = open(SphinxConfigFile,"r").read()
n1 = contents.find("release")
n2 = contents.find("\n",n1)
print(SphinxConfigFile,'says',contents[n1:n2])

contents = open(InstallSetupFile,"r").read()
n1 = contents.find("version")
n2 = contents.find(",",n1)
print(InstallSetupFile,'says',contents[n1:n2])

reply = input("Copy docs? (y/n)")
if reply=='y':
    os.system(r'robocopy build/html "c:/Programming/ExternalGit/pySCATMECH (docs)" * /e')
    
