import numpy as np
from pySCATMECH.model import *
from pySCATMECH.mueller import *
import SCATPY 
import os

class OpticalFunction():
    """
    Class for handling optical functions.

    Parameters
    ----------
    name: string, complex, float, or function returning complex or float
        The complex index of refraction. If a string, it can be a 
        "(n,k)" pair or the name of a file containing a table of 
        wavelengths, n, and k.  

    """
    def __init__(self, name, iter=None):
        # if a function, store the function
        if callable(name):
            self.function = name
            # create a temporary file
            from random import randint
            import string
            # Make a 12 character random name starting with a letter
            chars = string.ascii_letters
            self.filename = chars[randint(0,len(chars)-1)]
            # Numbers can follow...
            chars = chars + string.digits
            for i in range(11):
                self.filename += chars[randint(0,len(chars)-1)]
            self.filename += ".tmp"
            try:
                with open(self.filename,"w") as file:
                    for L in iter:
                        nk = self.function(L)
                        file.write("%.16g\t%.16g\t%.16g\n" %
                                   (L,nk.real,nk.imag))
            except Exception as e:
                # if there was an exception, delete the file
                os.remove(self.filename)
                del self.filename
                raise
            self.name = self.filename
        else: 
            self.name = name
            
    def __del__(self):
        if hasattr(self,'filename'):
            try: 
                os.remove(self.filename)
            except:
                print("Problem removing",self.filename)
               

    def __call__(self, wavelength):
        """
        Evaluate the optical constants n and k at wavelength.
        
        Parameters
        ----------

        wavelength: float 
                    The wavelength
        Returns
        -------
                    The complex optical constant
        """
        if hasattr(self,"function"):
            return self.function(wavelength)
        if type(self.name) is str:
            return SCATPY.GetOpticalConstant(wavelength,self.name)
        return self.name
            
    def nkstr(self,wavelength=1):
        """
        Evaluate the optical constants n and k at wavelength and returns a 
        string that can be used by SCATMECH::dielectric_function. 

        Parameters
        ----------
        
        wavelength: float
                    The wavelength
        Returns
        -------
            : str
            The complex optical constant as a string "(n,k)" or the name of
            a file.
        """        
        if isinstance(self.name,str):
            return self.name
        nk = self(wavelength)
        return str((nk.real,nk.imag)).replace(" ","")
    
    def __str__(self):
        return str(self.name)
    
    def writeToFile(self,filename,wavelengths):
        """
        Write the optical function to a file that can be used by SCATMECH 
        routines as a SCATMECH::dielectric_function.  

        Parameters
        ----------

        filename: str
                  The name of the file
        wavelengths: iterable of float
                  The wavelengths that will be written.

        Returns
        -------
            str
            The filename, if successful, 
            `(n,k)` at the first wavelength if not successful
        """
        try:
            with open(filename,'w') as file:
                for wavelength in wavelengths:
                    nk = complex(self(wavelength))
                    file.write(str(wavelength) + "\t" +
                               str(nk.real)+ "\t" + str(nk.imag) + "\n")
            return filename
        except Exception as e:
            print(e)
            return self.nkstr(wavelengths[0])
        
class Film():
    """
    Class for handling an optical film. 

    Film stores an OpticalFunction and a thickness. One feature is that 
    multiplication by a scalar returns a Film with its thickness multiplied 
    by that scalar.

    Parameters
    ----------
    material: OpticalFunction
              The optical function for the material

    thickness: float, optional
               The thickness of the film

    waves : float, optional
            The number of waves for the layer

    wavelength : float, optional
                 The wavelength the layer should be optimized for, if  
                 `waves` is set 

    angle : float, optional
            Vacuum angle in radians the layer should be optimized for, if 
            `waves` is set (default is 0)
    """
    def __init__(self, material, thickness=None, waves=None, wavelength=None, angle=0 ):

        if type(material) is OpticalFunction:
            self.material = material
        else:
            self.material = OpticalFunction(material)
            
        if thickness is not None:
            self.thickness = thickness
            return
        if waves is not None:
            n = self.material(wavelength)
            angle_inside = cmath.asin(cmath.sin(angle)/n)
            t = wavelength/n.real*waves/np.cos(angle_inside)
            self.thickness = t.real    

    def __mul__(self,a):
        """Returns a Film with a thickness multiplied by `a`"""
        return Film(self.material,thickness = self.thickness*a)
    
    def __rmul__(self,a):
        """Returns a Film with a thickness multiplied by `a`"""
        return Film(self.material,thickness = self.thickness*a)

    def __str__(self):
        return str((self.material,self.thickness))
   
    
class FilmStack(Model):
    """A class for handling stacks of films"""
 
    def __init__(self,films=[]):
        """
        Parameters
        ----------

        films : list of Film
                A list of films. films[0] is closest to the substrate.

        """
        self.newer = SCATPY.Get_StackModel
        self.freer = SCATPY.Free_StackModel
        self.handle = SCATPY.Get_StackModel("Stack_StackModel")
        try: 
            self.films = []
            for f in films:
                self.films.append(f)
        except TypeError:
            self.films = [films]
        
    def grow(self,film):
        """
        Grow a film on the stack.

        Parameters
        ----------

        film : Film
               The film to be grown
        """
        self.films.append(film)
        
    def clean(self):
        """
        Clear the stack of all films.
        """
        self.films = []
        
    def getStackCommand(self,wavelength):
        """
        Get the list of films suitable for a SCATMECH::Stack_StackModel

        
        Parameters
        ----------

        wavelength : float
                     The wavelength to evaluate the complex index at.

        Returns
        -------
        
        command : str
                  List of films suitable for a SCATMECH::Stack_StackModel
        """
        self.command = ""
        for f in self.films:
            materialstring = str((f.material.nkstr(wavelength)).replace(" ",""))
            thicknessstring = str(f.thickness)
            self.command = (self.command + materialstring +
                            " " + thicknessstring + " ")
        return self.command

    def getModelDict(self):
        """
        Return a Model parameter dictionary. Used primarily internally so that
        a FilmStack can be passed as a parameter to a Model for a stack.
        """
        dictionary = {None: "Stack_StackModel"}
        dictionary["stack"] = self.getStackCommand(1)
        return dictionary
    
    def __str__(self):
        dictionary = self.getModelDict()
        return str(dictionary)
        
    def reflectionCoefficient(self, theta, wavelength, no, nt,type="12"):
        """
        Return the reflection coefficent.

        Parameters
        ----------

        theta : float
                Incident angle in radians

        wavelength : float
                     Wavelength in vacuum

        no : OpticalFunction
             Optical constants for the incident medium

        nt : OpticalFunction
             Optical constants for the transmitted medium

        type : str, optional
               If `12`, then the incident medium is the medium above the 
               substrate and theta is the incident angle, evaluated as if 
               the radiation were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `21`, then the incident medium is the substrate and theta 
               is the incident angle evaluated is evaluated as if the radiation
               were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `12i`, then the incident medium is the medium above the 
               substrate and theta is the incident angle in that medium.
               (The component of the wavevector parallel to the surface is 
               2*pi*sin(theta)/no/lambda)
               If `21i`, then the incident medium is the substrate and theta 
               is the incident angle in the substrate. (The component of the 
               wavevector parallel to the surface is 2*pi*sin(theta)/n0/lambda)
               (Default is `12`)

        Returns
        -------

        coeff : list of list of complex
                The Jones matrix reflection coefficient, which is the
                linear relationship between the electric field amplitudes.
                coeff[0][0] is the s-polarized reflection coefficient.
                coeff[1][1] is the p-polarized reflection coefficient.
                coeff[1][0] = coeff[0][1] = 0
        """
        self.setParameter("stack",self.getStackCommand(wavelength))
        return SCATPY.ReflectionCoefficient(
                             self.handle, theta, wavelength, 
                             no.nkstr(), 
                             nt.nkstr(),
                             type)
    
    def transmissionCoefficient(self, theta, wavelength, no, nt,type="12"):
        """
        Return the transmission coefficient.

        Parameters
        ----------

        theta : float
                Incident angle in radians

        wavelength : float
                     Wavelength in vacuum

        no : OpticalFunction
             Optical constants for the incident medium

        nt : OpticalFunction
             Optical constants for the transmitted medium

        type : str, optional
               If `12`, then the incident medium is the medium above the 
               substrate and theta is the incident angle, evaluated as if 
               the radiation were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `21`, then the incident medium is the substrate and theta 
               is the incident angle evaluated is evaluated as if the radiation
               were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `12i`, then the incident medium is the medium above the 
               substrate and theta is the incident angle in that medium.
               (The component of the wavevector parallel to the surface is 
               2*pi*sin(theta)/no/lambda)
               If `21i`, then the incident medium is the substrate and theta 
               is the incident angle in the substrate. (The component of the 
               wavevector parallel to the surface is 2*pi*sin(theta)/n0/lambda)
               (Default is `12`)

        Returns
        -------

        coeff : list of list of complex
                The Jones matrix transmission coefficient, which is the
                linear relationship between the electric field amplitudes.
                coeff[0][0] is the s-polarized transmission coefficient.
                coeff[1][1] is the p-polarized transmission coefficient.
                coeff[1][0] = coeff[0][1] = 0
        """
        self.setParameter("stack",self.getStackCommand(wavelength))
        return SCATPY.TransmissionCoefficient(
                             self.handle, theta, wavelength, 
                             no.nkstr(), 
                             nt.nkstr(),
                             type)

    def R(self, theta, wavelength, no, nt, type="12"):
        """
        Return the reflectance.

        Parameters
        ----------

        theta : float
                Incident angle in radians

        wavelength : float
                     Wavelength in vacuum

        no : OpticalFunction
             Optical constants for the incident medium

        nt : OpticalFunction
             Optical constants for the transmitted medium

        type : str, optional
               If `12`, then the incident medium is the medium above the 
               substrate and theta is the incident angle, evaluated as if 
               the radiation were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `21`, then the incident medium is the substrate and theta 
               is the incident angle evaluated is evaluated as if the radiation
               were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `12i`, then the incident medium is the medium above the 
               substrate and theta is the incident angle in that medium.
               (The component of the wavevector parallel to the surface is 
               2*pi*sin(theta)/no/lambda)
               If `21i`, then the incident medium is the substrate and theta 
               is the incident angle in the substrate. (The component of the 
               wavevector parallel to the surface is 2*pi*sin(theta)/n0/lambda)
               (Default is `12`)

        Returns
        -------

        refl : MuellerMatrix
               The Mueller matrix reflectance, which is the
               linear relationship between radiant fluxes.
        """
        return JonesMueller(self.reflectionCoefficient(theta,
                                                       wavelength,
                                                       no, nt,type))

    def T(self, theta, wavelength, no, nt, type="12"):
        """
        Return the transmittance.

        Parameters
        ----------

        theta : float
                Incident angle in radians

        wavelength : float
                     Wavelength in vacuum

        no : OpticalFunction
             Optical constants for the incident medium

        nt : OpticalFunction
             Optical constants for the transmitted medium

        type : str, optional
               If `12`, then the incident medium is the medium above the 
               substrate and theta is the incident angle, evaluated as if 
               the radiation were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `21`, then the incident medium is the substrate and theta 
               is the incident angle evaluated is evaluated as if the radiation
               were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `12i`, then the incident medium is the medium above the 
               substrate and theta is the incident angle in that medium.
               (The component of the wavevector parallel to the surface is 
               2*pi*sin(theta)/no/lambda)
               If `21i`, then the incident medium is the substrate and theta 
               is the incident angle in the substrate. (The component of the 
               wavevector parallel to the surface is 2*pi*sin(theta)/n0/lambda)
               (Default is `12`)

        Returns
        -------

        refl : MuellerMatrix
               The Mueller matrix transmittance, which is the
               linear relationship between radiant fluxes.
        """
        No = no(wavelength)
        Nt = nt(wavelength)
        thetat = np.arcsin(complex(np.sin(theta)*No/Nt))
        factor = abs(Nt*np.cos(thetat)/np.abs(np.cos(theta))/No)
        return JonesMueller(self.transmissionCoefficient(
                            theta, wavelength, no, nt,type))*factor
    
    def Rs(self,theta,wavelength,no,nt,type="12"):
        """
        Return the reflectance for s-polarized radiation.

        Parameters
        ----------

        theta : float
                Incident angle in radians

        wavelength : float
                     Wavelength in vacuum

        no : OpticalFunction
             Optical constants for the incident medium

        nt : OpticalFunction
             Optical constants for the transmitted medium

        type : str, optional
               If `12`, then the incident medium is the medium above the 
               substrate and theta is the incident angle, evaluated as if 
               the radiation were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `21`, then the incident medium is the substrate and theta 
               is the incident angle evaluated is evaluated as if the radiation
               were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `12i`, then the incident medium is the medium above the 
               substrate and theta is the incident angle in that medium.
               (The component of the wavevector parallel to the surface is 
               2*pi*sin(theta)/no/lambda)
               If `21i`, then the incident medium is the substrate and theta 
               is the incident angle in the substrate. (The component of the 
               wavevector parallel to the surface is 2*pi*sin(theta)/n0/lambda)
               (Default is `12`)

        Returns
        -------

        refl : float
               The reflectance for s-polarization, which is the
               linear relationship between radiant fluxes.
        """
        return (self.R(theta,wavelength,no,nt,type)*StokesVector([1,1,0,0]))[0]

    def Rp(self,theta,wavelength,no,nt,type="12"):
        """
        Return the reflectance for p-polarized radiation.

        Parameters
        ----------

        theta : float
                Incident angle in radians

        wavelength : float
                     Wavelength in vacuum

        no : OpticalFunction
             Optical constants for the incident medium

        nt : OpticalFunction
             Optical constants for the transmitted medium

        type : str, optional
               If `12`, then the incident medium is the medium above the 
               substrate and theta is the incident angle, evaluated as if 
               the radiation were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `21`, then the incident medium is the substrate and theta 
               is the incident angle evaluated is evaluated as if the radiation
               were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `12i`, then the incident medium is the medium above the 
               substrate and theta is the incident angle in that medium.
               (The component of the wavevector parallel to the surface is 
               2*pi*sin(theta)/no/lambda)
               If `21i`, then the incident medium is the substrate and theta 
               is the incident angle in the substrate. (The component of the 
               wavevector parallel to the surface is 2*pi*sin(theta)/n0/lambda)
               (Default is `12`)

        Returns
        -------

        refl : float
               The reflectance for p-polarization, which is the
               linear relationship between radiant fluxes.
        """
        return (self.R(theta,wavelength,no,nt,type)*StokesVector([1,-1,0,0]))[0]

    def Ts(self,theta,wavelength,no,nt,type="12"):
        """
        Return the transmittance for s-polarized radiation.

        Parameters
        ----------

        theta : float
                Incident angle in radians

        wavelength : float
                     Wavelength in vacuum

        no : OpticalFunction
             Optical constants for the incident medium

        nt : OpticalFunction
             Optical constants for the transmitted medium

        type : str, optional
               If `12`, then the incident medium is the medium above the 
               substrate and theta is the incident angle, evaluated as if 
               the radiation were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `21`, then the incident medium is the substrate and theta 
               is the incident angle evaluated is evaluated as if the radiation
               were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `12i`, then the incident medium is the medium above the 
               substrate and theta is the incident angle in that medium.
               (The component of the wavevector parallel to the surface is 
               2*pi*sin(theta)/no/lambda)
               If `21i`, then the incident medium is the substrate and theta 
               is the incident angle in the substrate. (The component of the 
               wavevector parallel to the surface is 2*pi*sin(theta)/n0/lambda)
               (Default is `12`)

        Returns
        -------

        refl : float
               The transmittance for s-polarization, which is the
               linear relationship between radiant fluxes.
        """
        return (self.T(theta,wavelength,no,nt,type)*StokesVector([1,1,0,0]))[0]

    def Tp(self,theta,wavelength,no,nt,type="12"):
        """
        Return the transmittance for p-polarized radiation.

        Parameters
        ----------

        theta : float
                Incident angle in radians

        wavelength : float
                     Wavelength in vacuum

        no : OpticalFunction
             Optical constants for the incident medium

        nt : OpticalFunction
             Optical constants for the transmitted medium

        type : str, optional
               If `12`, then the incident medium is the medium above the 
               substrate and theta is the incident angle, evaluated as if 
               the radiation were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `21`, then the incident medium is the substrate and theta 
               is the incident angle evaluated is evaluated as if the radiation
               were in vacuum. (The component of the wavevector
               parallel to the surface is 2*pi*sin(theta)/lambda)
               If `12i`, then the incident medium is the medium above the 
               substrate and theta is the incident angle in that medium.
               (The component of the wavevector parallel to the surface is 
               2*pi*sin(theta)/no/lambda)
               If `21i`, then the incident medium is the substrate and theta 
               is the incident angle in the substrate. (The component of the 
               wavevector parallel to the surface is 2*pi*sin(theta)/n0/lambda)
               (Default is '12')

        Returns
        -------

        refl : float
               The transmittance for p-polarization, which is the
               linear relationship between radiant fluxes.
        """
        return (self.T(theta,wavelength,no,nt,type)*StokesVector([1,-1,0,0]))[0]


