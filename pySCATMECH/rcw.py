import SCATPY
from pySCATMECH.model import *
from pySCATMECH.mueller import *

class RCW_Model(Model):
    """
    Class handling rigorous coupled wave (RCW) analysis for 1D gratings.
    """
    def __init__(self, *args, **kwargs):
        """
        args and kwargs are forwarded to Model.setParameters() and contain
        model parameters
        """
        self.newer = SCATPY.Get_RCW_Model
        self.freer = SCATPY.Free_RCW_Model
        self.handle = SCATPY.Get_RCW_Model()
        self.setParameters(*args, **kwargs)
        
    def __del__(self):
        SCATPY.Free_RCW_Model(self.handle)

    def DiffractionEfficiency(self, i):
        """
        Returns the Mueller matrix diffraction efficiency for 
        the i-th diffraction order.

        Parameters
        ----------

        i : int
            Diffraction order

        Returns
        -------

        efficiency : MuellerMatrix
                     The Mueller matrix diffraction efficiency for order i

        """
        return MuellerMatrix(SCATPY.RCWDiffractionEfficiency(self.handle,i))

    def DiffractionAmplitude(self, i):
        """
        Returns the Jones matrix diffraction amplitude for 
        the i-th diffraction order.

        Parameters
        ----------

        i : int
            Diffraction order

        Returns
        -------

        amplitude : JonesMatrix
                     The Jones matrix diffraction amplitude for order i

        """
        return np.array(SCATPY.RCWDiffractionAmplitude(self.handle,i))

    def Direction(self, i):
        """
        Returns a 3D unit vector [x,y,z] in the direction of 
        propagation of the i-th diffraction order
        Parameters
        ----------

        i : int
            Diffraction order

        Returns
        -------
        direction : list of float
                    Directional cosines of diffracted order

        """
        return SCATPY.RCWDirection(self.handle, i)

    def getEpsilon(self, x, z, direction='x'):
        """
        Returns the grating dielectric constant at a specific location
        
        Parameters
        ----------

        x, z : float
               Location coordinates 

        direction : str, optional
                    If the grating is anisotropic, then direction can be set
                    to choose a particular field direction. Can be 'x', 'y', or 
                    'z'. (Default is 'x').

        Returns
        -------

        epsilon : complex
                  The dielectric constant evuated at the current wavelength 
                  and grating.
        """
        if direction=='x': d=0
        elif direction=='y': d=1
        elif direction=='z': d=2
        else: d=0

        return SCATPY.GetGratingEpsilon(self.handle,x,z,d)


    def getGratingDefinition(self):
        """
        Returns the grating definition.


        Returns
        -------

        Dictionary containing the following keys:

        "position"    : list of lists of x-positions of the transitions for 
                        each of the layers in the grating.
        "materialx"   : list of lists of the dielectric constants (xx element)
                        in each of the layers of the grating.
        "materialy"   : list of lists of the dielectric constants (yy element) 
                        in each of the layers of the grating.
        "materialz"   : list of lists of the dielectric constants (zz element) 
                        in each of the layers of the grating.
        "materialmux" : list of lists of the magnetic susceptibility (xx 
                        element) in each of the layers of the grating.
        "materialmuy" : list of lists of the magnetic susceptibility (yy 
                        element) in each of the layers of the grating.
        "materialmuz" : list of lists of the magnetic susceptibility (zz 
                        element) in each of the layers of the grating.
        "thickness"   : list of thicknesses of each layer of the grating.
        "heights"     : list of the heights of each interface of the grating.

        If the grating is a Generic_Grating, then an additional key "segments" 
        contains a list of line segments for each of the boundaries, each 
        element a dict containing keys "x1", "y1", "x2", "y2", "mat1", and
        "mat2".

        """

        definition =  SCATPY.GetGratingDefinition(self.handle)

        definition['heights'] = [0]
        accum = 0
        for thickness in definition['thickness']:
            accum = accum - thickness
            definition['heights'].append(accum)

        return definition
            
    
