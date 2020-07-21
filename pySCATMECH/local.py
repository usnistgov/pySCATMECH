import SCATPY
from pySCATMECH.model import *
from pySCATMECH.mueller import *

class Local_BRDF_Model(Model):
    """
    A class for handling local BRDF Models (ones that return differential 
    scattering cross section of isolated particles
    """
    def __init__(self, *args, **kwargs):
        """
        Create and allocate a SCATMECH Local_BRDF_Model 
        
        Parameters
        ----------
        name : string 
            The name of a SCATMECH model inheriting Local_BRDF_Model

        parameters : dict of string : string
            Initializing parameters for the model.  The dict items can be
            of the form 'param' : 'value', 'param.subparam' : 'value', 
            or 'param' : dict, where the first dict item is 
            None : 'submodelname' and 
            subsequent items are parameters of submodelname.
        
        Returns
        -------
        r : Local_BRDF_Model
            The object containing a handle to the requested Local_BRDF_Model
        """
        
        self.newer = SCATPY.Get_Local_BRDF_Model
        self.freer = SCATPY.Free_Local_BRDF_Model
        if len(args)>0 and type(args[0])==str:
                self.handle = self.newer(args[0])
                args = args[1::]
        else:
            self.handle = self.newer("Bobbert_Vlieger_BRDF_Model")
        self.setParameters(*args, **kwargs)

    def __del__(self):
        SCATPY.Free_Local_BRDF_Model(self.handle)

    def MuellerDSC(self, thetai=0, thetas=0, phis=0, rotation=0, coords="psps"):
        """
        Evaluate the Mueller matrix differential scattering cross section 
        in a given geometry, and coordinate system

        Parameters
        ----------
        thetai : float
            Polar angle of incidence in radians

        thetas : float
            Polar angle of viewing in radians 

        phis: float
            Azimuthal angle of viewing in radians 
            (Note: Specular occurs when thetai = thetas, phis = 0)

        rotation : float
            Sample rotation angle in radians

        coords : string
            The coordinate system for the Mueller matrix, 
            either "psps", "xyxy", or "plane"

        Returns
        -------
        r : MuellerMatrix 
            The Mueller matrix differential scattering cross section
        """
        return MuellerMatrix(
            SCATPY.LocalDSC(self.handle, thetai, thetas, phis,
                            rotation, coords))

    def DSC(self, thetai=0, thetas=0, phis=0, rotation=0, coords="psps",
            inc=StokesVector([1,0,0,0]),sens=StokesVector([1,0,0,0])):
        """
        Return the differential scattering cross section in a given geometry.

        Parameters
        ----------
        thetai : float
            Polar angle of incidence in radians

        thetas : float
            Polar angle of viewing in radians 

        phis : float
            Azimuthal angle of viewing in radians 
            (Note: Specular occurs when thetai = thetas, phis = 0)

        rotation : float
            Sample rotation angle in radians

        coords : string
            The coordinate system for the Mueller matrix, 
            either "psps", "xyxy", or "plane"

        inc : StokesVector
            Incident polarization as Stokes vector

        sens : StokesVector
            Polarization sensitivity of the viewer as a Stokes vector 

        Returns
        -------
        r : float
            The Mueller matrix differential scattering cross section for the 
            given incident and viewing directions and polarizations
        """
        return float((StokesVector(sens) @
                      MuellerMatrix(
                          SCATPY.LocalDSC(self.handle,
                                          thetai, thetas, phis, rotation,
                                          coords))) @ StokesVector(inc))

