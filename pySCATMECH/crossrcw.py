import SCATPY
from pySCATMECH.model import *
from pySCATMECH.mueller import *


class CrossRCW_Model(Model):
    """
    Class handling rigorous coupled wave (RCW) analysis for cross gratings.
    """
    def __init__(self, *args, **kwargs):
        """
        Arguments and keyword arguments are passed to Model.setParameters
        and contain parameters and values
        """
        self.handle = SCATPY.Get_CrossRCW_Model()
        self.setParameters(*args, **kwargs)
        self.newer =  SCATPY.Get_CrossRCW_Model
        self.freer =  SCATPY.Free_CrossRCW_Model
        
    def __del__(self):
        SCATPY.Free_CrossRCW_Model(self.handle)

    def DiffractionEfficiency(self, i, j):
        """
        Returns the Mueller matrix diffraction efficiency for the (i,j)-th 
        diffraction order

        Parameters
        ----------

        i, j : int
               Diffraction order

        Returns
        -------

        efficiency : MuellerMatrix
                     Diffraction efficiency, as a Mueller matrix, at order (i,j)

        """
        return MuellerMatrix(
            SCATPY.CrossRCWDiffractionEfficiency(self.handle, i, j))

    def Direction(self,i,j):
        """
        Returns a 3D unit vector [x,y,z] in the direction of propagation of 
        the (i,j)-th diffraction order.

        Parameters
        ----------

        i, j : int
               Diffraction order
        
        Returns
        -------

        direction : list of float
                    Directional cosines of diffracted direction

        """
        return SCATPY.CrossRCWDirection(self.handle, i, j)

