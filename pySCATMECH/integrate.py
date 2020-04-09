from pySCATMECH.brdf import *
from pySCATMECH.local import *
import matplotlib.pyplot as plt
import numpy as np

class SolidAngle:
    """
    Abstract class for a solid angle on a hemisphere.

    Inherited classes are expected to define inside(self,v) where
    v contains directional cosines and returns a tuple of 
    bool (`True` if (x,y) is inside solid angle, `False` otherwise) and
    the Stokes vector sensitivity.

    The parent class provides __and__, __or__, and __invert__
    functionality, so that solid angles can be defined by combining
    other solid angle classes.
    """
    def __and__(self, other):
        """
        Combines two SolidAngle with an `and` operation.
        Returns a new SolidAngle.
        """
        class newSolidAngle(SolidAngle):
            def __init__(s, a, b):
                s.a = a
                s.b = b

            def inside(s, v):
                ai = s.a.inside(v)
                bi = s.b.inside(v)
                if ai[0] and not bi[0]:
                    return False, ai[1]
                elif not ai[0] and bi[0]:
                    return False, bi[1]
                elif not ai[0] and not bi[0]:
                    return False, (ai[1] + bi[1])/2
                else:
                    return True, (ai[1] + bi[1])/2
        return newSolidAngle(self, other)
    
    def __or__(self, other):
        """
        Combines two SolidAngle with an `or` operation.
        Returns a new SolidAngle.
        """
        class newSolidAngle(SolidAngle):
            def __init__(s, a, b):
                s.a = a
                s.b = b
            def inside(s, v):
                ai = s.a.inside(v)
                bi = s.b.inside(v)
                if ai[0] and not bi[0]:
                    return True, ai[1]
                elif not ai[0] and bi[0]:
                    return True, bi[1]
                elif not ai[0] and not bi[0]:
                    return False, (ai[1] + bi[1])/2
                else: 
                    return True, (ai[1] + bi[1])/2
        return newSolidAngle(self, other)

    def __invert__(self):
        """
        Inverts a SolidAngle.
        Returns a new SolidAngle.
        """
        class newSolidAngle(SolidAngle):
            def __init__(s, a):
                s.a = a
            def inside(s, v):
                ai = s.a.inside(v)
                return not ai[0], ai[1] 
        return newSolidAngle(self)

class Hemisphere(SolidAngle):
    """
    A solid angle defined by the entire hemisphere.
    """
    def __init__(self, sensitivity = StokesVector(1,0,0,0), polphi = None):
        """
        Parameters
        ----------
        
        sensitivity : StokesVector, optional
                      The Stokes vector sensitivity, where signal is assumed
                      to be the inner product of this value with the incident
                      Stokes vector. The Stokes vector is represented in an s-p
                      basis set.
                      (default: No polarization sensitivity)

        polphi : None or float, optional
                 The azimuthal angle in radians where the polarization
                 sensitivity is defined. If None (default), the polarization 
                 sensitivity is defined in an x-y basis set.
        """
        if polphi is None:
            self.sensitivity = sensitivity
        else:
            self.sensitivity = MuellerRotator(polphi) @ sensitivity
        
    def inside(self, v):
        """
        Parameters
        ----------

        v : tuple of float
            Direction represented as directional cosines

        Returns
        -------
        is_inside : bool
                    True if v is inside solid angle, False otherwise

        sensitivity : StokesVector
                      The Stokes vector sensitivity
        """
        v2 = v[0]**2 + v[1]**2
        if v2>=1: return False, self.sensitivity
        return True, self.sensitivity
    
class ProjectedPolygon(SolidAngle):
    """
    A solid angle defined by an array of (x,y) points in projected cosine space.
    """
    def __init__(self, boundary,
                 sensitivity = StokesVector(1,0,0,0),
                 polphi = None):
        """
        Parameters
        ----------
        
        sensitivity : StokesVector, optional
                      The Stokes vector sensitivity, where signal is assumed
                      to be the inner product of this value with the incident
                      Stokes vector. 
                      (default: No polarization sensitivity)

        polphi : None, float, or str, optional
                 If None or "sp" (default), the polarization sensitivity is 
                 defined in an s-p basis at the centroid of the projected 
                 solid angle.
                 If float, the azimuthal angle in radians where the polarization
                 sensitivity is defined in the s-p basis.
                 If "xy", the sensitivity is defined in an x-y basis.
        """
        self.boundary = boundary
        self.boundary.append(boundary[0])

        # Find centroid:
        area = 0
        centroidx = 0
        centroidy = 0
        for a,b in zip(boundary,boundary[1:]+boundary[:1]):
            temp = a[0]*b[1]-a[1]*b[0]
            area += temp
            centroidx += (a[0]+b[0])*temp
            centroidy += (a[1]+b[1])*temp
        centroidx /= 6*area
        centroidy /= 6*area

        centphi = math.atan2(centroidy,centroidx)
        
        if polphi is None or polphi == 'sp':
            self.sensitivity = MuellerRotator(-centphi) @ sensitivity
        elif polphi=='xy':
            self.sensitivity = sensitivity
        else:
            self.sensitivity = MuellerRotator(-polphi) @ sensitivity

            
    def inside(self, v):
        """
        Parameters
        ----------

        v : tuple of float
            Direction represented as directional cosines

        Returns
        -------
        is_inside : bool
                    True if v is inside solid angle, False otherwise

        sensitivity : StokesVector
                      The Stokes vector sensitivity
        """
        v2 = v[0]**2 + v[1]**2        
        if v2>1:
            return False, self.sensitivity
        count = 0
        b0 = self.boundary[0]
        x, y = v[0], v[1]
        for b1 in self.boundary[1:]:
            x0, x1, y0, y1 = b0[0], b1[0], b0[1], b1[1]
            b0 = b1
            if x0 < x and x1 <= x: pass
            elif x0 > x and x1 >= x: pass
            elif y0 < y and y1 <= y: pass
            elif y0 >= y and y1 > y:
                if (x0 >= x and x1 < x) or (x0 <= x and x1 > x):
                    count = count + 1
            elif x1 == x0:
                if (y0 >= x and y1 < y) or ( y0 <= x and y1 > y):
                    count = count + 1
            else:
                yc = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
                if yc > y:
                    count = count + 1
        return count % 2 == 1, self.sensitivity
        
class EllipticalCone(SolidAngle):
    """
    A solid angle defined as an elliptical circular cone.
    """
    def __init__(self, theta = 0, phi = 0, alpha = 0, gamma = 0,
                 sensitivity = StokesVector(1,0,0,0), polphi = None):
        """
        A solid angle defined as an elliptical circular cone.

        Parameters
        ----------
        
        theta : float, optional
                Polar angle of center of cone [rad]
                (Default is 0)

        phi : float, optional
              Azimuth angle of center of cone [rad]
              (Default is 0)

        alpha : float or tuple of float
                If a float, half angle of cone [rad]
                If a tuple of float, (alphax, alphay), where
                alphax is the half angle of cone along the x-z plane before 
                rotation [rad], and alphay is the half angle of cone along the 
                y-z plane before rotation [rad]

        gamma : float, optional
                Rotation about the z axis [rad]
                (Default is 0)

        sensitivity : StokesVector, optional
                      The Stokes vector sensitivity, where signal is assumed
                      to be the inner product of this value with the incident
                      Stokes vector. The Stokes vector is represented in an s-p
                      basis set.
                      (default: No polarization sensitivity)

        polphi : None, float, or str, optional
                 If None or "sp" (default), the polarization sensitivity is 
                 defined in an s-p basis at the centroid of the projected 
                 solid angle.
                 If float, the azimuthal angle in radians where the polarization
                 sensitivity is defined in the s-p basis.
                 If "xy", the sensitivity is defined in an x-y basis.

        """
        if hasattr(alpha, '__len__'):
            alphax = alpha[0]
            alphay = alpha[1]
        else:
            alphax = alphay = alpha

        if polphi is None or polphi == 'sp':
            self.sensitivity = MuellerRotator(-phi) @ sensitivity
        elif polphi == 'xy':
            self.sensitivity = sensitivity
        else:
            self.sensitivity = MuellerRotator(-polphi) @ sensitivity

        sin = np.sin
        cos = np.cos
        
        self.theta = theta
        self.phi = phi
        self.alphax = alphax
        self.alphay = alphay
        self.gamma = gamma
        v = np.array([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])
        s = perpto(v, (1, 0, 0))
        p = perpto(v, s)
        self.rotation_matrix = np.array([p, s, v])
        self.sinalphax = sin(alphax)
        self.sinalphay = sin(alphay)
        v = np.array([cos(gamma), sin(-gamma), 0])
        s = perpto(v, (0, 0, 1))
        p = perpto(v, s)
        self.rotation_matrix = np.array([v, s, p]) @ self.rotation_matrix
        
    def inside(self, v):
        """
        Parameters
        ----------

        v : tuple of float
            Direction represented as directional cosines

        Returns
        -------
        is_inside : bool
                    True if v is inside solid angle, False otherwise

        sensitivity : StokesVector
                      The Stokes vector sensitivity
        """
        v2 = v[0]**2 + v[1]**2
        if v2 > 1: return False, self.sensitivity
        v = np.array([v[0], v[1], np.sqrt(1 - v2)])
        w = self.rotation_matrix @ v
        if w[0]**2 / self.sinalphax**2 + w[1]**2 / self.sinalphay**2 < 1:
            return True, self.sensitivity
        return False, self.sensitivity

# A pseudonym 
CircularCone = EllipticalCone

class RectangularCone(SolidAngle):
    """
    A solid angle defined by a rectangular collection cone.
    """
    def __init__(self, theta = 0, phi = 0, alpha = 0, gamma = 0,
                 sensitivity = StokesVector(1,0,0,0), polphi = None):
        """
        A solid angle defined by a rectangular collection cone.

        Parameters
        ----------
        
        theta : float, optional
                Polar angle of center of cone [rad]
                (Default is 0)

        phi : float, optional
              Azimuth angle of center of cone [rad]
              (Default is 0)

        alpha : float or tuple of float
                If a float, half angle of cone [rad]
                If a tuple of float, (alphax, alphay), where
                alphax is the half angle of cone along the x-z plane before 
                rotation [rad], and alphay is the half angle of cone along the 
                y-z plane before rotation [rad]

        gamma : float, optional
                Rotation about the z axis [rad]
                (Default is 0)
        
        sensitivity : StokesVector, optional
                      The Stokes vector sensitivity, where signal is assumed
                      to be the inner product of this value with the incident
                      Stokes vector. The Stokes vector is represented in an s-p
                      basis set.
                      (default: No polarization sensitivity)

        polphi : None, float, or str, optional
                 If None or "sp" (default), the polarization sensitivity is 
                 defined in an s-p basis at the centroid of the projected 
                 solid angle.
                 If float, the azimuthal angle in radians where the polarization
                 sensitivity is defined in the s-p basis.
                 If "xy", the sensitivity is defined in an x-y basis.
        """
        if hasattr(alpha, '__len__'):
            alphax = alpha[0]
            alphay = alpha[1]
        else:
            alphax = alphay = alpha

        if polphi is None or polphi == 'sp':
            self.sensitivity = MuellerRotator(-phi) @ sensitivity
        elif polphi == 'xy':
            self.sensitivity = sensitivity
        else:
            self.sensitivity = MuellerRotator(-polphi) @ sensitivity

        sin = np.sin
        cos = np.cos

        self.theta = theta
        self.phi = phi
        self.alphax = alphax
        self.alphay = alphay
        v = np.array([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])
        s = perpto(v, (1, 0, 0))
        p = perpto(v, s)
        self.rotation_matrix = np.array([p, s, v])
        self.sinalphax = sin(alphax)
        self.sinalphay = sin(alphay)
        v = np.array([cos(gamma), -sin(gamma), 0])
        s = perpto(v, (0, 0, 1))
        p = perpto(v, s)
        self.rotation_matrix = np.array([v, s, p]) @ self.rotation_matrix
    
    def inside(self, v):
        """
        Parameters
        ----------

        v : tuple of float
            Direction represented as directional cosines

        Returns
        -------
        is_inside : bool
                    True if v is inside solid angle, False otherwise

        sensitivity : StokesVector
                      The Stokes vector sensitivity
        """
        v2 = v[0]**2 + v[1]**2
        if v2 > 1: return False, self.sensitivity
        
        v = np.array([v[0], v[1], np.sqrt(1-v2)])
        w = self.rotation_matrix @ v
        if abs(w[0]) < self.sinalphax and abs(w[1]) < self.sinalphay:
            return True, self.sensitivity
        return False, self.sensitivity

# A pseudonym 
SquareCone = RectangularCone

class Integrator():
    """Class designed to integrate scattering models.

    Class designed to integrate a BRDF_Model (to obtain reflectance) or
    a Local_BRDF_Model (to obtain scattering cross section). 

 
    Attributes
    ----------
    x, y, z : ndarray(float)
           The x and y projected cosines

    w    : ndarray(float)
           The weights

    sens : list of StokesVectors
           The Stokes sensitivities

    theta : ndarray(float)
            The polar angles in radians.

    phi : ndarray(float)
          The azimuthal angles in radians.
 
    """    
    def __init__(self, gridsize, detector, type=1):
        """
        Parameters
        ----------
    
        gridsize : float 
               The size of the integration step in radians. In the direction 
               of the surface normal, the step size is asin(gridsize). 
        
        detector : SolidAngle
               A description of the detection system.

        type : int
           The method used for sampling type
           1 for evenly spaced integration in angle space, 
             [thetas*cos(phis),thetas*sin(phis)]
           2 for evenly spaced integration in solid angle, 
             [sqrt(1-cos(thetas))*cos(phis),sqrt(1-cos(thetas))*sin(phis)]
           3 for evenly spaced integration in projected solid angle, 
             [sin(thetas)*cos(phis),sin(thetas)*sin(phis)]
        """
        if (type == 1):
            def sinc(t):
                if t!=0: return math.sin(t)/t
                return 1
    
            def trans(x, y):
                theta = math.sqrt(x**2 + y**2) * pi / 2
                if theta > pi/2: return (1, 1, 0)
                sintheta = math.sin(theta)
                phi = math.atan2(y, x)
                return (sintheta * math.cos(phi),
                        sintheta * math.sin(phi),
                        math.cos(theta) * sinc(theta) * (gridsize*pi/2)**2)
        
            temp = [trans(x, y)
                    for x in np.arange(-1, 1, gridsize)
                    for y in np.arange(-1, 1, gridsize)
                    if x**2 + y**2 < 1]

        if (type == 2):
        
            def trans(x, y):
                sqrt2 = math.sqrt(2)
                costheta = 1 - x**2 - y**2
                sintheta = math.sqrt(1 - costheta**2)
                phi = math.atan2(y, x)
                return (sintheta * math.cos(phi), 
                        sintheta * math.sin(phi),
                        2 * costheta * gridsize**2)
        
            temp = [trans(x, y)
                    for x in np.arange(-1, 1, gridsize)
                    for y in np.arange(-1, 1, gridsize)
                    if x**2 + y**2 < 1]

        if (type==3):

            temp = [(x, y, gridsize**2)
                    for x in np.arange(-1, 1, gridsize)
                    for y in np.arange(-1, 1, gridsize)
                    if x**2 + y**2 < 1]

        dtemp = [(x, y, w, detector.inside((x, y))) for x, y, w in temp]
        dtemp = [(x, y, w, s[1]) for x, y, w, s in dtemp if s[0]]
                 
        self.x = np.array([t[0] for t in dtemp])
        self.y = np.array([t[1] for t in dtemp])
        self.z = np.sqrt(1 - self.x**2 - self.y**2)
        self.theta = np.arccos(self.z)
        self.phi = np.arctan2(self.y,self.x)
        self.w = np.array([t[2] for t in dtemp])
        self.sens = [t[3] for t in dtemp]
    
    def PlotSamplingPoints(self, figsize=5,
             expand=1, color='b', circlecolor='r'):
        """
        Graphs the integration points on the projected hemisphere.
    
        Parameters
        ----------
    
        figsize : float, optional
              Size of figure, which will be square (default is 5)
              
        expand : float, optional
             Scaling of sampling point size (default is 1) 
    
        color : str, optional
            Color to use for the sampling points (default is 'b')
            
        circlecolor : str, optional
                  Color to use for the horizon (default is 'r') 
          
        Requires: matplotlib 
        """
        import matplotlib.pyplot as plt
        
        # Create a curve to show the horizon...
        thetas = np.linspace(0, 2*pi, 100)
        circlex = np.cos(thetas)
        circley = np.sin(thetas)
    
        plt.figure(figsize = (figsize, figsize))
        ax = plt.subplot(111)
        ax.scatter(self.x, self.y,
                   s = self.w * 300 * figsize**2 * expand, c=color) 
        ax.plot(circlex, circley, c = circlecolor)
        plt.show()
    
    def Reflectance(self, model, thetai, 
                    incpol = StokesVector(1,0,0,0)):
        """
        Integrates a BRDF_Model and returns reflectance.
        
        Parameters
        ----------

        model :  BRDF_Model 
                 The model to be integrated.

        thetai : float
                 The incident angle in radians

        incpol : StokesVector
                 The incident Stokes vector

        Returns
        -------

        reflectance : float
                      The integrated reflectance
        """
        weights = self.w
        return self.BRDF(model, thetai, incpol) @ weights

    def BRDF(self, model, thetai, incpol = StokesVector(1,0,0,0)):
        """
        Returns the BRDF at each integration point.

        Parameters
        ----------

        model : BRDF_Model
                The BRDF model to be integrated

        thetai : float
                The incident angle in radians

        incpol : StokesVector
                 The incident Stokes vector

        Returns
        -------

        BRDF : np.array 
               An array of BRDF
        """
        return np.array([model.BRDF(thetai, theta, phi,
                                    coords="xyxy", inc=incpol,
                                    sens=sens)
                         for theta,phi,sens in zip(self.theta,
                                                   self.phi,
                                                   self.sens)])

    def CrossSection(self, model, thetai,
                     incpol = StokesVector(1,0,0,0)):
        """
        Integrates a Local_BRDF_Model to obtain differential scattering 
        cross section.

        Parameters
        ----------

        model:  Local_BRDF_Model 
                The model to be integrated.

        thetai: float
                The incident angle in radians

        incpol : StokesVector
                 Incident Stokes vector

        Returns
        -------

        crosssection : float
                       The integrated cross section
        """
        weights = self.w / self.z
        return self.DSC(model, thetai, incpol) @ weights

    def DSC(self, model, thetai, incpol = StokesVector(1,0,0,0)):
        """
        Returns the differential scattering cross section at each 
        integration point.

        Parameters
        ----------

        model : Local_BRDF_Model
                The BRDF model to be integrated

        thetai: float
                The incident angle in radians

        incpol : StokesVector
                 Incident Stokes vector

        Returns
        -------

        DSC : np.array 
              An array of differential scattering cross-section.
        """
        return np.array([model.DSC(thetai, theta, phi,
                                   coords="xyxy",
                                   inc=incpol, sens=sens)
                         for theta, phi, sens in zip(self.theta,
                                                     self.phi,
                                                     self.sens)])
