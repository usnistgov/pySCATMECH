#
# Test script for pySCATMECH. This will create lots of graphs and output. Should not crash!
#

print("""
This script exercises several features in each module of pySCATMECH.
All temporary files should be deleted at the end. 
Numerous windows will open up. Close each to continue.
Press return now to start.""")
input()

def WORKING():
    print("Working...",end="",flush=True)

def DONE():
    print("Done")


print("==========================")
print("Testing pySCATMECH.mueller")
print("==========================")

from pySCATMECH.mueller import *
import random 
import matplotlib.pyplot as plt

print("Creating some unit MuellerMatrix")
m1 = MuellerMatrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
m2 = MuellerMatrix([1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1])
m3 = MuellerMatrix(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1)
m4 = JonesMueller([[1,0],[0,1]])
m5 = MuellerUnit()
print(m5)

print("Mueller matrix maps Stokes vectors to Stokes vectors:",m1.valid())
print("Mueller matrix is a convex sum of non-depolarizing Mueller matrices:",m1.physically_valid())

# Returns a random number with a normal distribution...
ran = lambda : random.gauss(0,1)

print("Creating 20 totally random Mueller matrices.")
m1 = [MuellerMatrix([[ran(), ran(), ran(), ran()],
                    [ran(), ran(), ran(), ran()],
                    [ran(), ran(), ran(), ran()],
                    [ran(), ran(), ran(), ran()]]) for i in range(20)]

print("Number valid() is True:", sum([m.valid() for m in m1]))
print("Number physically_valid() is True:", sum([m.physically_valid() for m in m1]))

# Returns a random complex number with a normal distribution...
ranc = lambda : random.gauss(0,1)+random.gauss(0,1)*1j

print("Creating 20 random Jones-Mueller matrices")
m1 = [JonesMueller([[ranc(), ranc()], [ranc(), ranc()]]) for i in range(20)]

# See if they are valid...
print("Number value() is True:", sum([m.valid() for m in m1]))
print("Number physically valid() is True:", sum([m.physically_valid() for m in m1]))

def RandomMuellerMatrix(n):
    """
    Creates a random Mueller matrix by averaging `n` Jones-Mueller matrices.
    """
    ranc = lambda : (random.gauss(0,1)+random.gauss(0,1)*1j)/2
    return sum([JonesMueller([[ranc(),ranc()],[ranc(),ranc()]]) for i in range(n)])/n

print("A random Mueller matrix generated from sum of four random nondepolarizing Mueller matrices:")
m = RandomMuellerMatrix(4)
print(m)
print("m is a valid Stokes to Stokes mapper:",m.valid())
print("m is the sum of Jones-Mueller matrices:",m.physically_valid())

print("Unpolarized Stokes vector:", StokesVector((1,0,0,0)) )

print("Unpolarized Stokes vector:", StokesVector(1,0,0,0) )

print("s-polarized unit polarization from a Jones vector:", JonesStokes((1.,0.)) )

print("p-polarized unit Stokes vector:", Polarization('p') )

print("s-polarized unit Stokes vector:", Polarization(state = 's') )

print("Left-handed circularly polarized unit Stokes vector:", Polarization("L") )

print("Right-handed circularly polarized Stokes vector with intensity 2:", Polarization("R", I=2) )

print("Unpolarized Stokes vector:", Polarization("U") )

print("Unpolarized Stokes vector:", Polarization('s', DOP=0) )

print("Linear polarization at 45 degrees:", Polarization(angle = 45*deg) )

print("Elliptically polarized radiation:", Polarization(angle = 34*deg, ellipticity = 20*deg) )

print("Partially polarized radiation:", Polarization("S", DOP = 0.8) )

print("Partially polarized elliptical polarized radiation with intensity 3:\n", 
      Polarization(angle = 34*deg, ellipticity = 20*deg, DOP = 0.7, I = 3) )

def RandomStokesVector(n):
    ranc = lambda : (random.gauss(0,1) + random.gauss(0,1)*1j)/2
    return sum([JonesStokes([ranc(),ranc()]) for i in range(n)])/n

print("Creating random Stokes vectors with the sum of n random polarized Stokes vectors")
plt.figure()
plt.plot(range(1,100),[RandomStokesVector(i).DOP() for i in range(1,100)])
plt.title("Degree of polarization vs. number of Jones-Stokes vectors")
plt.xlabel("n")
plt.ylabel("DOP")
plt.show()


jmax = 100
imax = 100
print("Creating",imax*jmax,"random Mueller matrices")

WORKING()

mCD00sums = []
for i in range(1,imax):
    mCD00sum = np.array([0,0,0,0])
    for j in range(jmax):
        M = RandomMuellerMatrix(i)
        mCD = Cloude_Decomposition(M)
        mCD00sum = mCD00sum + np.array([m[0,0] for m in mCD])
    mCD00sums.append(mCD00sum/jmax)
DONE()
plt.figure()
plt.plot(mCD00sums)
plt.xlabel("Number of Jones-Mueller matrices averaged")
plt.ylabel("Cloude decomposition")
plt.show()

def printMatrix(name,M):
    print("\n",name,"=\n")
    print(M)

print("Creating a random Mueller matrix.")
M = RandomMuellerMatrix(4)

MCD = M.Cloude_Decomposition()
printMatrix("M",M)

print("\nCloude decomposition")
print("--------------------\n")

for m in MCD:
    print(m,"\n")

depol, ret, diatten = Lu_Chipman_Decomposition(M)
print("\nLu-Chipman Decomposition")
print("------------------------\n")

printMatrix("Depolarizer",depol)
printMatrix("Retarder",ret)
printMatrix("Diattenuator",diatten)
printMatrix("The product Depolarizer . Retarder .  Diattenuator",depol @ ret @ diatten)

diatten, ret, depol = Reverse_Lu_Chipman_Decomposition(M)
print("\nReverse Lu-Chipman Decomposition")
print("--------------------------------\n")
printMatrix("Diattenuator",diatten)
printMatrix("Retarder",ret)
printMatrix("Depolarizer",depol)
printMatrix("The product Diattenuator . Retarder .  Depolarizer", diatten @ ret @ depol)

diatten2, ret2, depol, ret1, diatten1 = Symmetric_Decomposition(M)
print("\nSymmetric decomposition")
print("-----------------------\n")
printMatrix("diatten2",diatten2)
printMatrix("ret2",ret2)
printMatrix("depol",depol)
printMatrix("ret1",ret1)
printMatrix("diatten1",diatten1)
printMatrix("The product Diattenuator2 . Retarder2 .  Depolarizer . Retarder1 . Diattenuator1",
            diatten2 @ ret2 @ depol @ ret1 @ diatten1)

L = MuellerLog(M)
print("\nLogarithmic decomposition")
print("-------------------------\n")
printMatrix("L",L)
printMatrix("MuellerExp(L)",MuellerExp(L))

print("\nSome other characteristic parameters")
print("------------------------------------\n")
print("Tmax = ",M.Tmax())
print("Tmin = ",M.Tmin())
print("diattenuation = ",M.diattenuation())
print("linear diattenuation = ",M.linear_diattenuation())
print("polarization dependent loss = ",M.polarization_dependent_loss())
print("polarizance = ",M.polarizance())
print("depolarization index = ",M.depolarization_index())
print("extinction_ratio = ",M.extinction_ratio())

parametersM = CharacterizedMueller(M)
print(parametersM)

print("==========================")
print("Testing pySCATMECH.fresnel")
print("==========================")

from pySCATMECH.fresnel import *

# We will be plotting some stuff, so we need matplotlib.pyplot
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

wavelengths = np.arange(0.2,2.5,0.01)

# The following looks for a file named "silicon"
# Assume that "silicon.txt" exists on your computer...
silicon = OpticalFunction("silicon.txt") 

# A fixed real index:
glass = OpticalFunction(1.5) 

# A metal, i.e., a complex index of refraction:
aMetal = OpticalFunction(1.4+2.1j)

# A material with a Cauchy dispersion, given as a lambda function
# When providing a function, one must provide an iterable list of wavelengths because 
# the function is tabulated.
SiO2 = OpticalFunction(lambda L: 1.4580 + 0.00354/L**2, wavelengths)

# A couple materials with Sellmeier functions...
MgF = OpticalFunction(lambda L: math.sqrt(1+0.48755108*L**2/(L**2-0.04338408**2)
                                           +0.39875031*L**2/(L**2-0.09461442**2)
                                           +2.3120353*L**2/(L**2-23.793604**2)),
                                             wavelengths)

ZnO = OpticalFunction(lambda L: math.sqrt(1+1.347091*L**2/(L**2-0.062543**2)
                                           +2.117788*L**2/(L**2-0.166739**2)
                                           +9.452943*L**2/(L**2-24.320570**2)),
                                             wavelengths)

# It is useful to have one for vacuum:
vacuum = OpticalFunction(1)

wavelengths = np.arange(0.25,1.5,0.001)

# Plot the curves...
plt.figure()
plt.plot(wavelengths, [silicon(w).real for w in wavelengths], label="Si n")
plt.plot(wavelengths, [silicon(w).imag for w in wavelengths], label="Si k")
plt.plot(wavelengths, [glass(w) for w in wavelengths], label="glass")
plt.plot(wavelengths, [SiO2(w) for w in wavelengths], label="SiO2")
plt.plot(wavelengths, [MgF(w) for w in wavelengths], label="MgF")
plt.plot(wavelengths, [ZnO(w) for w in wavelengths], label="ZnO")
plt.legend()
plt.show()

# Create a film of ZnO (defined above) of thickness 0.5 micrometers...
film1 = Film(ZnO, thickness = 0.5)

# Create a quarter-wave thickness of ZnO for wavelength 0.5 micrometers...
film2 = Film(ZnO, waves = 1/4, wavelength = 0.5)

# create a quarter-wave thickness of ZnO for wavelength 0.5 micrometers for 45 degrees incidence...
film3 = Film(ZnO, waves = 1/4, wavelength = 0.5, angle = 45*deg)

print(film1)
print(film2)
print(film3)

H = Film(ZnO, waves = 1/4, wavelength = 0.5)
L = Film(MgF, waves = 1/4, wavelength = 0.5)

notch_filter = 5*[H,L] + [2*H] + 5*[L,H]

stack = FilmStack()

wavelength = 0.500 

# Return a Jones matrix
r = stack.reflectionCoefficient(60*deg, wavelength, vacuum, silicon)
print("Jones matrix = \n",r,"\n")

# Return a Mueller matrix
R = stack.R(60*deg, wavelength, vacuum, silicon)
print("Mueller matrix = \n",R,"\n")

# Return the s and p reflectance
Rs = stack.Rs(60*deg, wavelength, vacuum, silicon) 
Rp = stack.Rp(60*deg, wavelength, vacuum, silicon) 
print("Rs = ",Rs,", Rp = ",Rp)

plt.figure()
thetas=[t for t in np.arange(0,90)]
Rs = [stack.Rs(t*deg, wavelength, vacuum, silicon) for t in thetas]
Rp = [stack.Rp(t*deg, wavelength, vacuum, silicon) for t in thetas]
plt.plot(thetas, Rs, label="$R_\mathrm{s}$")
plt.plot(thetas, Rp, label="$R_\mathrm{p}$")
plt.legend()
plt.title("Reflectance of silicon at " + str(wavelength) + " µm")
plt.xlabel(r"$\theta$ (degrees)")
plt.ylabel("Reflectance")
plt.show()

plt.figure()
wavelengths = [L for L in np.arange(0.1, 1.5, 0.01)]
Rs = [stack.Rs(70*deg, wavelength, vacuum, silicon) for wavelength in wavelengths]
Rp = [stack.Rp(70*deg, wavelength, vacuum, silicon) for wavelength in wavelengths]
plt.plot(wavelengths, Rs, label="$R_\mathrm{s}(\lambda)$")
plt.plot(wavelengths, Rp, label="$R_\mathrm{p}(\lambda)$")
plt.legend()
plt.title("Reflectance of silicon at 70° incidence")
plt.xlabel("$\lambda$ (µm)")
plt.ylabel("Reflectance")
plt.show()

plt.figure()
wavelengths = [L for L in np.arange(0.1, 1.5, 0.01)]
M = [stack.R(70*deg, wavelength, vacuum, silicon) for wavelength in wavelengths]
delta = [(m @ Polarization(angle=45*deg)).delta() for m in M]
psi = [(m @  Polarization(angle=45*deg)).psi() for m in M]
plt.plot(wavelengths, delta, label="$\Delta(\lambda)$")
plt.plot(wavelengths, psi, label="$\Psi(\lambda)$")
plt.legend()
plt.title("Ellipsometric parameters")
plt.xlabel("$\lambda$ (µm)")
plt.ylabel('$\Psi$, $\Delta$ (rad)')
plt.show()

stack = FilmStack()
film = Film(SiO2, thickness = 0.01)
stack.grow(film)

plt.figure()
wavelengths = [L for L in np.arange(0.1, 1.5, 0.01)]

for t in [0., 0.01, 0.02, 0.03]:
    stack.clean()
    stack.grow(Film(SiO2, thickness = t))
    
    M = [stack.R(70*deg, wavelength, vacuum, silicon) for wavelength in wavelengths]
    delta = [(m @ StokesVector(1,0,1,0)).delta() for m in M]
    psi = [(m @ StokesVector(1,0,1,0)).psi() for m in M]
    plt.plot(wavelengths, delta)
    plt.plot(wavelengths, psi)

plt.title("Ellipsometric parameters")
plt.xlabel("$\lambda$ (µm)")
plt.ylabel('$\Psi$, $\Delta$ (rad)')
plt.show()

stack = FilmStack()

substrate = SiO2
wavelength = 0.500
stack.grow(Film(MgF, waves=1/4, wavelength=wavelength))

plt.figure()
ax = plt.subplot()
thetas=np.array([t for t in np.arange(0,90)])
Rs = [stack.Rs(theta*deg, wavelength, vacuum, substrate) for theta in thetas]
Rp = [stack.Rp(theta*deg, wavelength, vacuum, substrate) for theta in thetas]
Ts = [stack.Ts(theta*deg, wavelength, vacuum, substrate) for theta in thetas]
Tp = [stack.Tp(theta*deg, wavelength, vacuum, substrate) for theta in thetas]
plt.plot(thetas, Rs, label = "Rs")
plt.plot(thetas, Rp, label = "Rp")
plt.plot(thetas, Ts, label = "Ts")
plt.plot(thetas, Tp, label = "Tp")
plt.legend()
plt.xlim([-5, 95])
plt.xticks(np.linspace(0, 90, 10))
plt.xlabel("Angle")
plt.ylabel("Reflectance, Transmittance")
plt.title("Reflectance and Transmittance of a $\lambda/4$ MgF film on SiO$_2$")
class AddDegrees(ticker.Formatter):
    def __call__(self,x,pos):
        return "$%g \degree$"% x
ax.xaxis.set_major_formatter(AddDegrees())
class AddPercent(ticker.Formatter):
    def __call__(self,x,pos):
        return ("%g" % (x*100))+" %"
ax.yaxis.set_major_formatter(AddPercent())
plt.show()

L = Film(MgF, waves=1/4, wavelength=0.500)
H = Film(ZnO, waves=1/4, wavelength=0.500)

substrate = glass

def pltFilm(stack,label=""):
    stack = FilmStack(stack)
    wavelengths = np.arange(0.2, 1.0, 0.001)
    Rs = [stack.Rs(0*deg, L, vacuum, substrate) for L in wavelengths]
    plt.plot(wavelengths, Rs, label=label)

WORKING()
plt.figure()
plt.title("Buildup of a high reflector")
pltFilm(2*[L,H], label="2 LH")
pltFilm(3*[L,H], label="3 LH")
pltFilm(4*[L,H], label="4 LH")
pltFilm(10*[L,H], label="10 LH")
pltFilm(20*[L,H], label="20 LH")
plt.legend()
plt.xlabel("Wavelength (µm)")
plt.ylabel("Reflectance")
DONE()
plt.show()

WORKING()
plt.figure()
plt.title("Notch filter")
pltFilm(5*[H,L] + [2*H] + 5*[L,H],"5(H,L)+2H+5(L,H)")
plt.legend()
plt.xlabel("Wavelength (µm)")
plt.ylabel("Reflectance")
DONE()
plt.show()

print("========================")
print("Testing pySCATMECH.model")
print("========================")

from pySCATMECH.model import *

modelDictionary = getModelDictionary('Roughness_BRDF_Model')
printDictOfDict(modelDictionary)

print("Creating an ABC_PSD_Function")
psdFunction = Model("ABC_PSD_Function")
print("printParameters() returns:")
print(psdFunction.printParameters())

print("\ngetParameters() returns:")
print(psdFunction.getParameters())

print("\ngetParameterDictionary() returns:")
print(printDictOfDict(psdFunction.getParameterDictionary()))

psdFunction.setParameter('A',0.2)
psdFunction.getParameters()

parameters ={None : 'Bobbert_Vlieger_BRDF_Model',
              'lambda': '0.633',
              'substrate': '(1.56,0.2)',
              'type': '0',
              'density': '1',
              'sphere': '(1.59,0)',
              'radius': '1',
              'spherecoat': 'No_StackModel',
              'stack': 'SingleFilm_StackModel',
              'delta': '0',
              'quiet': 1} 
model = Model(parameters)
print(model)

printDictOfDict(model.getParameterDictionary())

model.getParameter('radius')

model.getModelName()



print("A dialog box should open. Choose a model, by double clicking on it.")
model = DialogGetModel('BRDF_Model')
print("""
A dialog box should open showing that model's parameters. Change parameters
as you like and then close window.""")
m = DialogGetModelParameters(Model(model))
print(model)

print("=======================")
print("Testing pySCATMECH.brdf")
print("=======================")

from pySCATMECH.brdf import *
import matplotlib.pyplot as plt

model = BRDF_Model('Correlated_Roughness_Stack_BRDF_Model')
print(model)

psd = {None : "ABC_PSD_Function",
       "A" : 0.01,
       "B" : 360,
       "C" : 2.4}

stack = {None : 'SingleFilm_StackModel',
         'material' : 1.59,
         'thickness': 0.05
        }

parameters = {'lambda' : 0.532,
              'substrate' : 4.05+0.05j,
              'type' : 0,
              'stack' : stack,
              'psd' : psd}

model = BRDF_Model('Correlated_Roughness_Stack_BRDF_Model',parameters)
print(model)

parameters['lambda'] = 0.600
psd['C'] = 2.3
stack['thickness'] = 0.1
model.setParameters(parameters)
print(model)

# Needed to use OpticalFunction, Film, and FilmStack...
from pySCATMECH.fresnel import *

# Define two materials...
SiO2 = OpticalFunction(lambda L: 1.4580 + 0.00354/L**2,np.arange(0.2,1.5,0.1))
Si = OpticalFunction('silicon')

# Define the power spectral density function...
psd = Model("ABC_PSD_Function", A=0.01, B=360, C=2.4)

# Define the stack...
stack =  FilmStack([Film(SiO2,thickness =0.05)])

model.setParameters(wavelength = 0.600, substrate=Si, stack = stack, psd = psd)
print(model)

model.getParameters()

print(model.MuellerBRDF(60*deg, 30*deg))

mBRDF = model.MuellerBRDF(60*deg, 30*deg)
inc = Polarization('p')
sens = Sensitivity('u')
print(sens @ mBRDF @ inc)

print(model.BRDF(60*deg, 30*deg, inc=Polarization('p'), sens=Sensitivity('u')))

print(model.BRDF(60*deg, 30*deg))

thetaslist = np.arange(-89, 90, 1)
BRDFlist = [model.BRDF(60*deg, thetas*deg, inc=[1,-1,0,0], sens=[1,0,0,0]) for thetas in thetaslist]

plt.figure()
plt.yscale('log')
plt.plot(thetaslist,BRDFlist)
plt.ylabel('BRDF (sr$^{-1}$)')
plt.xlabel(r'$\theta_\mathrm{r}$ (degrees)')
plt.xticks(np.linspace(-90,90,7))
plt.show()

thetaslist = np.arange(-89,90,1) 
plt.figure()
plt.yscale('log')
for thetai in [0,30,60]:
    BRDFlist = [model.BRDF(thetai*deg,thetas*deg,inc=[1,-1,0,0],sens=[1,0,0,0]) for thetas in thetaslist]
    plt.plot(thetaslist,BRDFlist,label = r"$\theta_i=%d^\circ$" % thetai)

plt.xticks(np.linspace(-90,90,7))
plt.legend()
plt.show()

thetaslist = np.linspace(-89,89,100) 

plt.figure()
for t in [0.1,0.2,0.3,0.4,0.5]:
    stack = FilmStack([Film(SiO2,thickness = t)])
    model.setParameters(stack=stack)
    BRDFlist = [model.BRDF(60*deg,thetas*deg,inc=[1,-1,0,0],sens=[1,0,0,0]) for thetas in thetaslist]
    plt.plot(thetaslist,BRDFlist,label = r"$D = %g$" % t)
    
plt.yscale('log')
plt.xticks(np.linspace(-90,90,7))
plt.legend()
plt.show()

print("========================")
print("Testing pySCATMECH.local")
print("========================")

import matplotlib.pyplot as plt
from pySCATMECH.local import *

spherecoat = {None : 'No_StackModel'}

stack = {None : 'No_StackModel'}

parameters = {'lambda' : 0.532,
              'substrate' : 4.06+0.05j,
              'type' : 0,
              'sphere' : 1.59,
              'radius' : 0.05,
              'spherecoat' : spherecoat,
              'stack' : stack,
              'delta' : 0,
              'lmax' : 0,
              'order' : -1,
              'Norm_Inc_Approx' : 0,
              'improve' : 3,
              'quiet' : 1}

model = Local_BRDF_Model("Bobbert_Vlieger_BRDF_Model",parameters)
print(model)

# Needed to use OpticalFunction, Film, and FilmStack...
from pySCATMECH.fresnel import *

# Define two materials...
PS = OpticalFunction(lambda L: (1+1.4435*L**2/(L**2-0.020216))**(1/2), np.arange(0.4356,1.052,0.01))
SiO2 = OpticalFunction(lambda L:(1+0.6961663*L**2/(L**2-0.0684043**2)+
                                   0.4079426*L**2/(L**2-0.1162414**2)+
                                   0.8974794*L**2/(L**2-9.896161**2))**(1/2), 
                                    np.exp(np.linspace(math.log(0.210),math.log(6.7),50)))

Si = OpticalFunction('silicon')

# Define the stack...
spherecoat =  FilmStack()
stack = FilmStack([Film(SiO2,thickness =0.0017)])

parameters = {'lambda' : 0.305,
              'substrate' : Si,
              'type' : 0,
              'sphere' : PS,
              'radius' : 0.05,
              'spherecoat' : spherecoat,
              'stack' : stack,
              'delta' : 0,
              'lmax' : 0,
              'order' : -1,
              'Norm_Inc_Approx' : 0,
              'improve' : 3}

model.setParameters(parameters)
print(model)

mDSC = model.MuellerDSC(60*deg, 30*deg)
print(mDSC)

thetaslist = np.linspace(-89, 89, 100)
DSClist = [model.DSC(60*deg, thetas*deg, inc=Polarization('p'), sens=Sensitivity('u')) for thetas in thetaslist]

plt.plot(thetaslist,DSClist)
plt.xlabel(r"$\theta_s$ [degrees]")
plt.ylabel(r"DSC [$\mathrm{\mu m}^2/\mathrm{sr}$]")
plt.xlim((-90,90))
plt.xticks(np.linspace(-90,90,7))
plt.show()

thetaslist = np.linspace(-89, 89, 100)

plt.figure()
for D in [0.1, 0.2, 0.5, 1.0]:
    parameters["radius"] = D/2
    model.setParameters(parameters)
    DSClist = [model.DSC(60*deg, thetas*deg, inc=Polarization('p'), sens=Sensitivity('u')) for thetas in thetaslist]
    plt.plot(thetaslist, DSClist, label = "D = %g µm" % D)

plt.xlabel(r"$\theta_s$ [degrees]")
plt.ylabel(r"DSC [$\mathrm{\mu m}^2/\mathrm{sr}$]")
plt.xlim((-90,90))
plt.yscale('log')
plt.xticks(np.linspace(-90,90,7))
plt.legend()
plt.show()


print("============================")
print("Testing pySCATMECH.integrate")
print("============================")

import matplotlib.pyplot as plt
from pySCATMECH.integrate import *

hemi = Hemisphere()

integrator1 = Integrator(2*deg, hemi, type=1)
integrator1.PlotSamplingPoints()

integrator2 = Integrator(2*deg, hemi, type=2)
integrator2.PlotSamplingPoints()

integrator3 = Integrator(2*deg, hemi, type=3)
integrator3.PlotSamplingPoints()

cone = CircularCone(theta=45*deg, phi=90*deg, alpha=20*deg)
    
integrator = Integrator(1*deg, cone)
integrator.PlotSamplingPoints()

cone = EllipticalCone(theta=0*deg, phi=0*deg, alpha=(20*deg, 45*deg))
    
integrator = Integrator(1*deg, cone)
integrator.PlotSamplingPoints()

cone = EllipticalCone(theta=0*deg, phi=0*deg, alpha=(20*deg, 45*deg), gamma=45*deg)
    
integrator = Integrator(1*deg, cone)
integrator.PlotSamplingPoints()


cone = RectangularCone(theta=0*deg, phi=0*deg, alpha=(20*deg, 45*deg))
    
integrator = Integrator(1*deg, cone)
integrator.PlotSamplingPoints()

cone = RectangularCone(theta=0*deg, phi=0*deg, alpha=(20*deg, 45*deg), gamma=45*deg)
    
integrator = Integrator(1*deg, cone)
integrator.PlotSamplingPoints()

cone = ProjectedPolygon([(0.5,0.5), 
                         (0.5,-0.5), 
                         (0,-0.2), 
                         (-0.5,-0.5), 
                         (-0.5,0.5), 
                         (0,0.2)])
    
integrator = Integrator(1*deg, cone)
integrator.PlotSamplingPoints()

hole1 = CircularCone(theta=70*deg, phi=0, alpha=5*deg)
hole2 = CircularCone(theta=70*deg, phi=180*deg, alpha=5*deg)

detector =  Hemisphere() & ~(hole1 | hole2)

integrator = Integrator(1*deg, detector)
integrator.PlotSamplingPoints()

circle = CircularCone(theta=0, phi=0, alpha=70*deg)
hole1 = CircularCone(theta=70*deg, phi=0, alpha=5*deg)
hole2 = CircularCone(theta=70*deg, phi=180*deg, alpha=5*deg)
hole3 = CircularCone(theta=0, phi=0, alpha=30*deg)
center = CircularCone(theta=0, phi=0, alpha=20*deg)

detector = circle & ~(hole1 | hole2 | hole3) | center

integrator = Integrator(1*deg, detector)
integrator.PlotSamplingPoints()

slot = RectangularCone(theta=0, phi=0, alpha=(90*deg, 6*deg))

detector3 = ~slot

integrator = Integrator(1*deg, detector3)
integrator.PlotSamplingPoints()

# The model parameters as a dictionary...
parameters = {'lambda' : 0.532,
              'substrate' : 4.05+0.05j,
              'type' : 0,
              'psd' : {None : 'Gaussian_PSD_Function',
                       'sigma' : 0.05,
                       'length' : 1
                       }
             }

# Define the model and set its parameters...
model = BRDF_Model("Microroughness_BRDF_Model", **parameters)

# Define the incident angle...
thetaIncident = 6*deg

# The detector is the hemisphere with a hole cut out of it to let in the 
# incident light and to let out the specular reflection...
detector = Hemisphere() & ~CircularCone(theta=0*deg, alpha=10*deg)

# Get the integration points and show them...
integrator = Integrator(1*deg, detector)
integrator.PlotSamplingPoints()

# Calculate and print the reflectance ...
ref = integrator.Reflectance(model, thetaIncident) 
print("Integrated reflectance = %g" % ref)

# Get the RMS roughnesses in logarithmic spacing...
sigmas = np.exp(np.linspace(np.log(0.0001), np.log(.01), 20)) 

# Initialize a list of reflectances...
refs = []

# s-polarized incident radiation...
spol = StokesVector(1, 1, 0, 0)

WORKING()
for sigma in sigmas:
    # Set the RMS roughness...
    model.setParameter("psd.sigma", sigma)
    
    # Calculate the reflectance...
    ref = integrator.Reflectance(model, thetaIncident, incpol=spol) 
    
    # Add it to the list..
    refs.append(ref)
DONE()

# Plot it...
plt.figure()
plt.loglog(sigmas,refs)
plt.xlabel("$\sigma$")
plt.ylabel("Reflectance")
plt.title("Reflectance versus RMS Roughness")
plt.show()

BVparameters = {'lambda' : 0.532,
              'substrate' : "silicon",
              'type' : 0, 
              'sphere' : 1.59,
              'radius' : 0.05,
              'spherecoat' : "No_StackModel",
              'stack' : {None : 'SingleFilm_StackModel',
                        'material' : '(1.5,0)',
                        'thickness' : '0.0016'},                         
                'delta' : 0,
                'quiet' : 1} # use defaults for operational parameters

BVmodel = Local_BRDF_Model("Bobbert_Vlieger_BRDF_Model",**BVparameters)

thetaIncident = 70*deg

circle = CircularCone(theta=0, phi=0, alpha=70*deg)
hole1 = CircularCone(theta=70*deg, phi=0, alpha=5*deg)
hole2 = CircularCone(theta=70*deg, phi=180*deg, alpha=5*deg)
hole3 = CircularCone(theta=0*deg, phi=0, alpha=30*deg)

detector = circle & ~(hole1 | hole2 | hole3)

integrator = Integrator(2*deg, detector)
integrator.PlotSamplingPoints()

cs = integrator.CrossSection(BVmodel,thetaIncident,incpol=[1,-1,0,0])

print("Cross section = %g µm^2" % cs)

diameters = np.exp(np.linspace(np.log(0.05),np.log(2),20)) #logarithmic spacing

refp = []
refs = []
spol = Polarization('s')
ppol = Polarization('p')
WORKING()
for d in diameters:
    BVmodel.setParameters(radius=d/2)
    refs.append( integrator.CrossSection(BVmodel, thetaIncident, incpol=spol) )
    refp.append( integrator.CrossSection(BVmodel, thetaIncident, incpol=ppol) )
DONE()
plt.figure()
plt.loglog(diameters,refs,label="s")
plt.loglog(diameters,refp,label="p")
plt.xlabel("Diameter / $\mathrm{\mu m}$")
plt.ylabel("Cross section / $\mathrm{\mu m}^2$")
plt.title("Cross Section versus Diameter")
plt.legend()
plt.show()

hemi = Hemisphere()
integrator = Integrator(1*deg, hemi)
dsc = integrator.DSC(BVmodel, thetaIncident, incpol=ppol)
integrator.PlotSamplingPoints(color=np.log(dsc), expand=3)

detector = CircularCone(theta=45*deg, phi=90*deg, alpha=30*deg, sensitivity=Sensitivity('s'))
integrator = Integrator(1*deg, detector)

diameters = np.exp(np.linspace(np.log(0.05), np.log(2),20)) #logarithmic spacing

WORKING()
ref = []
for d in diameters:
    BVmodel.setParameters(radius=d/2)
    ref.append( integrator.CrossSection(BVmodel, thetaIncident, incpol=Polarization('p')) )
DONE()

plt.figure()
plt.loglog(diameters, ref)
plt.xlabel("Diameter / $\mathrm{\mu m}$")
plt.ylabel("Cross section / $\mathrm{\mu m}^2$")
plt.title("Cross Section versus Diameter")
plt.show()

print("============================")
print("Testing pySCATMECH.scatterer")
print("============================")

import matplotlib.pyplot as plt
import numpy as np
from pySCATMECH.scatterer import *

parameters = {
            'lambda' : 0.532,
            'medium' : 1,
            'radius' : 0.1,
            'sphere' : 1.59}

model = Free_Space_Scatterer("MieScatterer",parameters)

print(model)

theta = 45*deg

vin = [0,0,1] 
vout = [np.sin(theta),0,np.cos(theta)]

m = model.DifferentialScatteringCrossSection(vin,vout)
print(m)

def dsc(theta):
    vin = [0,0,1]
    vout = [np.sin(theta),0,np.cos(theta)]
    return model.DifferentialScatteringCrossSection(vin,vout)
    
angles = np.linspace(0,2*pi,180)
results = [dsc(angle) for angle in angles]

spol = Polarization('S') # perpendicular to the scattering plane
ppol = Polarization('P') # parallel to the scattering plane
upol = Polarization('U') # unpolarized 
usens = Sensitivity('U') # unpolarized sensitivity

plt.figure()
plt.polar(angles, [usens @ r @ upol for r in results], label = "unpol")
plt.polar(angles, [usens @ r @ spol for r in results], label = "perp")
plt.polar(angles, [usens @ r @ ppol for r in results], label = "par")
plt.legend()
plt.show()

parameters = {
            'lambda' : 0.532,
            'medium' : 1,
            'radius' : 0.03,
            'sphere' : 'gold'}

model = Free_Space_Scatterer("MieScatterer", parameters)

def integrate_scatter(model, wavelength, steps):
    "Integrate `model` at `wavelength`, performing the integral in `steps` steps."
    integ = 0
    dtheta = pi / (steps - 1)
    parameters['lambda'] = wavelength
    model.setParameters(parameters)
    vin = [0, 0, 1]
    angles = np.arange(0, pi, dtheta)
    # Integrate angles from 0 to pi...
    for theta in angles:
        vout = [np.sin(theta), 0, np.cos(theta)]
        m = model.DifferentialScatteringCrossSection(vin, vout)
        integ += m[0][0] * np.abs(np.sin(theta))
    integ *= 2 * pi * dtheta
    return integ

wavelengths = np.linspace(0.380,0.780,200)

plt.figure()

# Diameter 0.04 µm...
D = 0.04
parameters['radius'] = D/2
model.setParameters(parameters)

crosssection = [integrate_scatter(model,L,50) for L in wavelengths]

plt.plot(wavelengths, crosssection, label = "D = %g µm"%D)

plt.legend(frameon=False)
plt.show()

parameters = {
            'lambda' : 0.532,
            'medium' : 'water', 
            'radius' : 0.03,
            'sphere' : 'glass',
            'coating' : 'gold',
            'thickness' : 0.005}

model = Free_Space_Scatterer("CoatedMieScatterer", parameters)

def extinction(model, wavelength):
    "Calculate the extinction coefficient at a given `wavelength` for a given `model`."
    parameters['lambda'] = wavelength
    model.setParameters(parameters)
    vin = [0, 0, 1]
    return model.Extinction(vin)[0][0]

def integrated_scatter(model,wavelength,steps):
    "Integrate the scatter over polar angle"
    integ = 0
    dtheta = pi/(steps - 1)
    parameters['lambda'] = wavelength
    model.setParameters(parameters)
    vin = [0, 0, 1]
    for theta in np.linspace(0, pi, steps):
        vout = [np.sin(theta), 0, np.cos(theta)]
        m = model.DifferentialScatteringCrossSection(vin, vout)
        integ += m[0][0]*np.abs(np.sin(theta))
    integ *= 2*pi*dtheta
    return integ

wavelengths = np.linspace(0.500, 1.300, 200)

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

WORKING()

for t in [0.002,0.003,0.004,0.005,0.006,0.007,0.008]:
    parameters['thickness'] = t
    sca = np.array([integrated_scatter(model,L,50) for L in wavelengths])
    ext = np.array([extinction(model,L) for L in wavelengths])
    abs = ext-sca
    ax1.plot(wavelengths, sca, label = "t = %g"%t)
    ax2.plot(wavelengths, ext, label = "t = %g"%t)
    ax3.plot(wavelengths, abs, label = "t = %g"%t)
DONE()

ax1.legend()
ax1.set_title("Scattering")
ax2.set_title("Extinction")
ax3.set_title("Absorption")
plt.show()

parameters = {
            'lambda' : 0.589,
            'medium' : 1,
            'radius' : 1,
            'materialright' : 1.56+ 6.9E-7,
            'materialleft' : 1.56}

model = Free_Space_Scatterer("Optically_Active_Sphere_Scatterer",parameters)

def extinction(model,wavelength):
    model.setParameters(wavelength=wavelength)
    vin = [0,0,1]
    return model.Extinction(vin)

extmatrix = model.Extinction([0,0,1])

print(extmatrix)

shape = {None : 'Ellipsoid_Axisymmetric_Shape',
        'npoints' : '100',
        'vertical' : '0.010',
        'horizontal' : '0.003'
        }

tmatrix = {'lambda' : 0.8,
        'medium' : 1.4,
        'Shape' : shape,
        'particle' : 'gold',
        'lmax' : '30',
        'mmax' : '0'}

model = Free_Space_Scatterer("TMatrix_Axisymmetric_Scatterer",tmatrix)

extCoeff100 = model.Extinction([1,0,0])
extCoeff001 = model.Extinction([0,0,1])

print("Mueller matrix transmittance perpendicular to ellipsoid axis:")
print(MuellerExp(-extCoeff100 * 10**5),'\n')
print("Mueller matrix transmittance parallel to ellipsoid axis:")
print(MuellerExp(-extCoeff001 * 10**5),'\n')

print("======================")
print("Testing pySCATMECH.rcw")
print("======================")


import matplotlib.pyplot as plt
import numpy as np
from pySCATMECH.rcw import *


gratingParameters = {
    None : 'Single_Line_Grating', # The Grating model can be assigned with the `None` key
    'period' : 0.2,               # The period  
    'medium_i' : 1,               # The incident medium0
    'medium_t' : 'silicon',       # The transmitting (substrate) material
    'material' : 'silicon',       # The grating material
    'space' : 1,                  # The space material
    'height' : 0.2,               # The height
    'topwidth' : 0.1,             # Top width
    'bottomwidth' : 0.15,         # Bottom width
    'offset' : 0,                 # Offset of top, relative to bottom
    'nlevels' : 10}               # Number of layers in the staircase approximation

rcwParameters = {
    'order' : 25,                 # number of Floquet orders considered
    'type' : 0,                   # 0 indicates reflection
    'lambda' : 0.532,             # wavelength
    'thetai' : 70,                # incident angle in degrees 
    'rotation' : 0,               # rotation angle in degrees
    'grating' : gratingParameters,
    'quiet' : 1
    }

model = RCW_Model(rcwParameters)

print(model)

# 250 wavelengths from 0.2 to 1.5...
wavelengths = np.linspace(0.2, 1.5, 250)


def rcw(wavelength):
    rcwParameters['lambda'] = wavelength
    model.setParameters(rcwParameters)
    return model.DiffractionEfficiency(0)

WORKING()
muellerMatrices = [rcw(wavelength) for wavelength in wavelengths]
DONE()

N = [m[0,1] for m in muellerMatrices]
C = [m[2,2] for m in muellerMatrices]
S = [m[2,3] for m in muellerMatrices]

plt.plot(wavelengths*1000, N, label="N")
plt.plot(wavelengths*1000, C, label="C")
plt.plot(wavelengths*1000, S, label="S")
plt.xlabel('$\lambda$ / nm')
plt.legend()
plt.show()

wavelengths = np.linspace(0.2, 1.5, 75)

rcwParameters['rotation'] = 45
model.setParameters(rcwParameters)


def rcw(wavelength):
    rcwParameters['lambda'] = wavelength
    model.setParameters(rcwParameters)
    return model.DiffractionEfficiency(0)

WORKING()
muellerMatrices = [rcw(wavelength) for wavelength in wavelengths]
DONE()

plt.figure(figsize = [13,13])
for element in range(16):
    plt.subplot(4,4,element+1)
    if element==0:
        plt.ylim(0,1)
    else:
        plt.ylim(-1,1)
    plt.plot(wavelengths,[m[element//4,element%4] for m in muellerMatrices],'b-')

plt.show()

def showGrating(model,x,z):
    xx, zz = np.meshgrid(x, z, sparse=True)

    def eps(x,z):
        return model.getEpsilon(x,z).real

    epsilon = np.array([[eps(x0,z0) for x0 in x] for z0 in z])

    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111)
    h = ax.contourf(x,z,epsilon,256,cmap="brg")
    plt.colorbar(h,ax=ax)

    plt.show()
    
showGrating(model,np.linspace(-0.2, 0.2, 200),np.linspace(-0.3,.1,200))

gratingParameters = {
    None : 'Single_Line_Grating',
    'period' : 10,
    'medium_i' : 1,
    'medium_t' : 'silicon',  
    'material' : 'silicon',
    'space' : 1,
    'height' : '0.5',
    'topwidth' : '0.4',
    'bottomwidth' : '0.15',
    'offset' : '0',
    'nlevels' : '10'}

rcwParameters = {
    None : 'RCW_Model',
    'order' : '25',
    'type' : '0',
    'lambda' : '0.532',
    'thetai' : 45,
    'rotation' : 45,
    'grating' : gratingParameters,
    'quiet' : 1
    }

model.setParameters(rcwParameters)

def plotDiffractionPattern(model):
    WORKING()
    order = int(model.getParameter("order"))
    x, y, I = [], [], []
    for i in range(-order,order):
        direction = model.Direction(i)
        if direction[2] != 0 :
            x.append(direction[0])
            y.append(direction[1])
            I.append(model.DiffractionEfficiency(i)[0,0])
       
    thetas = np.linspace(0, 2*pi, 100)
    circlex = np.cos(thetas)
    circley = np.sin(thetas)
    
    plt.figure(figsize = (10,10))
    ax = plt.subplot(111)
    ax.scatter(x, y, s=np.array(I)*100) 
    ax.plot(circlex, circley, c = 'r')
    DONE()
    plt.show()

plotDiffractionPattern(model)

# Set the lines per millimeter:
lpm = 600
# Set the nominal wavelength in µm:
blazeLambda = 1.00

##########################
period = 1000/lpm
sinBlaze = blazeLambda/2/period
height = period*sinBlaze
##########################

# Description of the grating...
gratingParameters = {None: 'Triangular_Grating', 
                     'period': 1000/lpm, 
                     'medium_i': 1, 
                     'medium_t': 'aluminum', 
                     'material': 'aluminum',
                     'amplitude' : height,
                     'aspect' : 1,
                     'nlevels' : 20}

# The RCW parameters
rcwParameters = {
    None : 'RCW_Model',
    'order' : 25,
    'type' : 0,
    'lambda' : 0.500,
    'thetai' : 0,
    'rotation' : 0,
    'grating' : gratingParameters,
    'quiet' : 1
    }
model = RCW_Model(rcwParameters)

def efficiency(wavelength):
    """Set the wavelength and calculate the efficiency"""
    thetai=math.asin(wavelength/2/period)
    model.setParameters(wavelength=wavelength,thetai=thetai)
    return model.DiffractionEfficiency(1)

wavelengths = np.arange(0.2, 1.5, 0.02)

WORKING()
eff = [efficiency(L) for L in wavelengths]
DONE()

effs = [Sensitivity('u') @ e @ Polarization('s') for e in eff]
effp = [Sensitivity('u') @ e @ Polarization('p') for e in eff]

plt.figure()
plt.plot(wavelengths, effs, label = "s-pol")
plt.plot(wavelengths, effp, label = "p-pol")
plt.legend()
plt.show()

showGrating(model,np.linspace(-6, 6, 200),np.linspace(-6,6,200))

# Set the lines per millimeter:
lpm = 1000
# Sinusoid amplitude:
amplitude = .5
# Incident angle:
thetai = 45*deg


period = 1000/lpm
from pySCATMECH.fresnel import OpticalFunction

FusedSilica = OpticalFunction(lambda L: 1.4580 + 0.00354/L**2, np.arange(0.2,1.5,0.05))

# Description of the grating...
gratingParameters = {None: 'Sinusoidal_Relief_Grating', 
                     'period': period, 
                     'medium_i': 1, 
                     'medium_t': FusedSilica, 
                     'material': FusedSilica,
                     'amplitude' : amplitude,
                     'base' : 0,
                     'option' : 0,
                     'nlevels' : 20}

# The RCW parameters
rcwParameters = {
    None : 'RCW_Model',
    'order' : 10,
    'type' : 1,
    'lambda' : 0.500,
    'thetai' : 0,
    'rotation' : 0,
    'grating' : gratingParameters,
    'quiet' : 1
    }
model = RCW_Model(rcwParameters)

def efficiency(wavelength):
    """Set the wavelength and calculate the efficiency"""
    rcwParameters['lambda'] = wavelength
    model.setParameters(rcwParameters)
    return model.DiffractionEfficiency(1)

wavelengths = np.arange(0.2, 1.5, 0.02)

WORKING()
eff = [efficiency(L) for L in wavelengths]
DONE()

effs = [Sensitivity('u') @ e @ Polarization('s') for e in eff]
effp = [Sensitivity('u') @ e @ Polarization('p') for e in eff]

plt.figure()
plt.plot(wavelengths, effs, label = "s-pol")
plt.plot(wavelengths, effp, label = "p-pol")
plt.legend()
plt.show()

showGrating(model,np.linspace(-6, 6, 200),np.linspace(-6,6,200))

print("===========================")
print("Testing pySCATMECH.crossrcw")
print("===========================")


import matplotlib.pyplot as plt
import numpy as np
from pySCATMECH.crossrcw import *

period = 0.5
order = 5
theta = 70

gratingParameters = {
    None : 'Cylinder_CrossGrating',
              'medium_i' : 1,
              'medium_t' : "fusedsilica",
              'zeta' : 0,
              'd1' : period,
              'd2' : period,
              'grid1' : 256,
              'grid2' : 256,
              'rtop' : 0.1,
              'rbottom' : 0.1,
              'thickness' : 0.1,
              'nlevels' : 1,
              'inside' : 1,
              'outside' : "fusedsilica"}

rcwParameters = {
    'thetai' : theta,
    'rotation' : 0,
    'lambda' : 0.532,
    'type' : '0',
    'order1' : order,
    'order2' : order,
    'grating' : gratingParameters
    }

model = CrossRCW_Model(rcwParameters)
print(model)

wavelengths = np.linspace(0.2, 1.5, 100)

def rcw(wavelength):
    model.setParameter('lambda',wavelength)
    return model.DiffractionEfficiency(0,0)

WORKING()
muellerMatrices = [rcw(wavelength) for wavelength in wavelengths]
DONE()

N = [m[0,1] for m in muellerMatrices]
C = [m[2,2] for m in muellerMatrices]
S = [m[2,3] for m in muellerMatrices]

plt.figure()
plt.plot(wavelengths*1000, N, label="N")
plt.plot(wavelengths*1000, C, label="C")
plt.plot(wavelengths*1000, S, label="S")
plt.xlabel('$\lambda$/nm')
plt.legend()
plt.show()

um2eV = 1.239841775 
energies = np.linspace(um2eV/0.2,um2eV/1.5,10)

rcwParameters['rotation'] = 22.5

def rcwA(energy):
    rcwParameters["lambda"] = um2eV/energy
    model.setParameters(rcwParameters)
    m = model.DiffractionEfficiency(0,0)
    m, m[0,0] = m/m[0,0],m[0,0] # Normalize the matrix, keeping the m[0,0] element
    return m

WORKING()
muellerMatrices = [rcwA(E) for E in energies]
DONE()

plt.figure(figsize = [13,13])
for element in range(16):
    plt.subplot(4,4,element+1)
    plt.ylim(-1,1)
    plt.plot(energies,[muellerMatrices[i][element//4,element%4] for i in range(len(energies))],'b:')

plt.show()

period = 2
order = 10
theta = 25

gratingParameters = {
    None : 'Cylinder_CrossGrating',
              'medium_i' : 1, 
              'medium_t' : 1, # Free standing
              'zeta' : 30,    # Hexagonal pattern
              'd1' : period,
              'd2' : period,
              'grid1' : 256,
              'grid2' : 256,
              'rtop' : period/4,
              'rbottom' : period/4,
              'thickness' : 100, # Very deep
              'nlevels' : 1,
              'inside' : 1,     # Holes in
              'outside' : 1.5   # glass
            }

rcwParameters = {
    'thetai' : theta,
    'rotation' : 0,
    'lambda' : 0.532,
    'type' :   1,  # in transmission
    'order1' : order,
    'order2' : order,
    'grating' : gratingParameters
    }

model.setParameters(rcwParameters)

def plotDiffractionPattern(model):
    WORKING()
    order1 = int(model.getParameter("order1"))
    order2 = int(model.getParameter("order2"))
    x, y, I = [], [], []
    for i in range(-order1,order1):
        for j in range(-order2,order2):
            # Direction() function gets the direction 
            direction = model.Direction(i,j) 
            if direction[2] != 0 :
                x.append(direction[0])
                y.append(direction[1])
                I.append(model.DiffractionEfficiency(i,j)[0,0])
       
    thetas = np.linspace(0, 2*pi, 100)
    circlex = np.cos(thetas)
    circley = np.sin(thetas)

    plt.figure(figsize = (10,10))
    ax = plt.subplot(111)
    ax.scatter(x, y, s=np.array(I)*2000) # You may need to adjust this size scaling
    ax.plot(circlex, circley, c = 'r')
    DONE()
    plt.show()

plotDiffractionPattern(model)

import os
# The following file got created by TMatrix_Axisymmetric_Scatterer and should be removed...
os.remove('shape.dat')

print("========")
print("Finished")
print("========")
