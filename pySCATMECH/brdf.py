import SCATPY
from pySCATMECH.model import *
from pySCATMECH.mueller import *

class BRDF_Model(Model):
    """
    A class for handling BRDF Models. 

    Arguments *args and **kwargs can contain parameters passed to the model. 

    One form is to pass a dict of {"parameter" : value, ... }.

    The keyword argument 'wavelength' is automatically changed to 'lambda',
    since `lambda` is a Python keyword.  However, as a str, it should be 
    'lambda'.

    Parameters
    ----------
    name : str
        The name of a SCATMECH model inheriting BRDF_Model

    """
    
    def __init__(self,*args, **kwargs):

        self.newer = SCATPY.Get_BRDF_Model
        self.freer = SCATPY.Free_BRDF_Model
        if len(args)>0 and type(args[0])==str:
                self.handle = self.newer(args[0])
                args = args[1::]
        else:
            self.handle = self.newer("Microroughness_BRDF_Model")
        self.setParameters(*args, **kwargs)


    def __del__(self):
        SCATPY.Free_BRDF_Model(self.handle)

    def MuellerBRDF(self,
                    thetai=0, thetas=0, phis=0, rotation=0, coords="psps"):
        """
        Evaluate the Mueller matrix BRDF in a given geometry, and coordinate 
        system.

        Parameters
        ----------
        thetai : float
            Polar angle of incidence in radians.
        thetas : float
            Polar angle of viewing in radians.
        phis: float, optional
            Azimuthal angle of viewing in radians.
            (Note: Specular occurs when thetai = thetas, phis = 0)
            (Default is 0)
        rotation : float, optional
            Sample rotation angle in radians. (Default is 0)
        coords : string, optional
            The coordinate system for the Mueller matrix, 
            either "psps", "xyxy", or "plane." (Default is "psps")

        Returns
        -------
        BRDF : MuellerMatrix 
            The Mueller matrix BRDF
        """
        return MuellerMatrix(SCATPY.BRDF(self.handle,
                                         thetai, thetas, phis,
                                         rotation, coords))

    def BRDF(self,
             thetai=0, thetas=0, phis =0, rotation=0,
             coords="psps",
             inc=StokesVector([1,0,0,0]),
             sens=StokesVector([1,0,0,0])):
        """
        Return the scalar BRDF in a given geometry.

        Parameters
        ----------
        thetai : float
            Polar angle of incidence in radians.
        thetas : float
            Polar angle of viewing in radians. 
        phis : float, optional
            Azimuthal angle of viewing in radians.
            (Note: Specular occurs when thetai = thetas, phis = 0)
            (default is 0)
        rotation : float, optional
            Sample rotation angle in radians. (Default is 0)
        coords : string, optional
            The coordinate system for the Mueller matrix, 
            either "psps", "xyxy", or "plane". (Default is `psps`)
        inc : StokesVector or 4-list of float or int, optional
            Incident polarization as Stokes vector.
            (Default is unpolarized)
        sens : StokesVector or 4-list of float or int, optional
            Polarization sensitivity of the viewer as a Stokes vector. 
            (Default is unpolarized)

        Returns
        -------
        BRDF : float
            The Mueller matrix BRDF for the given incident and viewing 
            directions and polarizations.
        """
        return float((StokesVector(sens) @
                      MuellerMatrix(SCATPY.BRDF(self.handle, thetai, thetas,
                                                phis, rotation, coords))) @
                     StokesVector(inc))

