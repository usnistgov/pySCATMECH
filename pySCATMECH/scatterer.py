import SCATPY
from pySCATMECH.model import *
from pySCATMECH.mueller import *

class Free_Space_Scatterer(Model):
    """A class for handling free space scatterer models"""

    def __init__(self, *args, **kwargs):
        """
        Creates and allocates a SCATMECH Free_Space_Scatterer 
        
        Parameters
        ----------
        name : str
            The name of a SCATMECH model inheriting Free_Space_Scatterer 

        args and kwargs : dict of str : str
            Initializing parameters for the model.  The dict items can be
            of the form 'param' : 'value', 'param.subparam' : 'value', 
            or 'param' : dict, where the first dict item is 
            None : 'submodelname' and 
            subsequent items are parameters of submodelname.
        
        Returns
        -------
        r : FSScatterer
            The object containing a handle to the requested Free_Space_Scatterer
        """

        self.newer = SCATPY.Get_Free_Space_Scatterer
        self.freer = SCATPY.Free_Free_Space_Scatterer
        if len(args)>0 and type(args[0])==str:
                self.handle = self.newer(args[0])
                args = args[1::]
        else:
            self.handle = self.newer("MieScatterer")
        self.setParameters(*args, **kwargs)

    def __del__(self):
        SCATPY.Free_Free_Space_Scatterer(self.handle)

    def DifferentialScatteringCrossSection(self, vi, vo):
        """
        Returns the Mueller matrix differential scattering cross section for 
        a given geometry

        Arguments
        ---------
        vi : 3-vector 
             The incident direction
        
        vo : 3-vector 
             The scattering direction

        Returns
        -------
        DSC : MuellerMatrix
              Mueller matrix differential scattering cross section
        """
        return JonesMueller(SCATPY.FSSjones(self.handle,
                                        vi[0], vi[1], vi[2],
                                        vo[0], vo[1], vo[2]))

    def ScatteringMatrix(self, vi, vo):
        """
        Returns the Mueller matrix differential scattering cross section for 
        a given geometry

        Arguments
        ---------
        vi : 3-vector 
             The incident direction
        
        vo : 3-vector 
             The scattering direction

        Returns
        -------
        s : numpy.array
              Scattering matrix as a Jones matrix
        """
        return np.array(SCATPY.FSSjones(self.handle,
                               vi[0], vi[1], vi[2],
                               vo[0], vo[1], vo[2]))

    def Extinction(self, v):
        """
        Returns the Mueller matrix extinction cross section for 
        a given geometry

        Arguments
        ---------
        v : 3-vector 
             The incident direction
        
        Returns
        -------
        e : MuellerMatrix
            Mueller matrix extinction cross section
        """
        return MuellerMatrix(SCATPY.FSSext(self.handle,v[0],v[1],v[2]))

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
        vin = [-np.sin(thetai),0,-np.cos(thetai)]
        vout = [np.cos(phis)*np.sin(thetas), np.sin(phis)*np.sin(thetas), np.cos(thetas)]
        return float((StokesVector(sens) @ 
                      self.DifferentialScatteringCrossSection(vin,vout)) @ 
                            StokesVector(inc).rotate(phis))