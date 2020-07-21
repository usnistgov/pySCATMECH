import numpy as np
import math
import cmath

def unit(a):
    """
    Return the unit matrix associated with a vector.
    
    Parameters
    ----------
    a : np.array
        Unnormalized vector

    Returns
    -------
    unit : np.array  
           Vector such that a.dot(a) == 1.

    """
    return a/np.linalg.norm(a)

def perpto(a,b):
    """
    Return a unit vector perpendicular to two vectors.

    Parameters
    ----------
    a, b : np.array
           Two 3-vectors

    Returns
    -------
    c : np.array
        A unit vector perpendicular to both a and b.

    """
    c = np.cross(a,b)
    if np.linalg.norm(a)!=0:
        return unit(np.cross(a,b))
    else:
        c = np.cross(a,np.array([1,0,0]))
        if np.linalg.norm(a)!=0:
            return unit(np.cross(a,b))
        else:
            c = np.cross(a,np.array([0,0,1]))
            if np.linalg.norm(a)!=0:
                return unit(np.cross(a,b))
            else:
                return [1,0,0]
            
pi = math.pi # A useful symbol
deg = pi/180 # Another useful symbol
            
class MuellerMatrix(np.ndarray):
    """
    A class for handling Mueller matrices 

    Parameters
    ----------
    x : 4x4 list or array of float or int, optional
        Default is a zero Mueller matrix

    """
    def __new__(cls,*x):
        if (len(x)==0):
            m = np.array(((0.,0.,0.,0.),
                          (0.,0.,0.,0.),
                          (0.,0.,0.,0.),
                          (0.,0.,0.,0.))).view(MuellerMatrix)
            return m
        if (len(x)==16):
            m = np.asarray(x,dtype=float).view(MuellerMatrix)
            m = np.resize(m,(4,4))
            return m
        if (len(x)==4):
            for xx in x:
                if (len(xx)!=4):
                    raise Exception("Mueller matrix must be 4x4 real matrix")
            x = [x[0],x[1],x[2],x[3]]
            m = np.asarray(x,dtype=float).view(MuellerMatrix)
            return m
        if (len(x)==1):
            m = np.asarray(x,dtype=float).view(MuellerMatrix)
            m = m.reshape((4,4))
            return m
        raise Exception("Mueller matrix must be a 4x4 real matrix")

 
    def __mul__(a,b):
        """Multiply two Mueller matrices, or a Mueller matrix with 
           a Stokes vector, or multiply by a scalar 
        """
        if isinstance(b,MuellerMatrix):
            return a.dot(b)
        if isinstance(b,StokesVector):
            return StokesVector(a.dot(b))
        if isinstance(b,(float,int)):
            return MuellerMatrix(a.dot(b))
        raise Exception("Mueller matrices should only be multiplied "
                "by Mueller matrices, Stokes vectors, floats, or ints")

    __matmul__ = __mul__
    
    def __truediv__(a,b):
        "Divide a Mueller matrix by a scalar"
        if isinstance(b,int) or isinstance(b,float):
            return super().__truediv__(b)
        raise Exception("Can only divide a Mueller matrix by a float or int")

    def __pow__(a,b):
        """
        Mueller matrix to the power b.

        Parameters
        ----------

        b : float
            The power

        Returns
        -------

            : MuellerMatrix
              Matrix power operation

        """
        Q, W = np.linalg.eig(a)
        Q = Q**b
        return W.dot(np.diag(Q)).dot(np.linalg.inv(W)).real
        
    def Tmax(self):
        """
        Return the maximum transmittance.
        
        Returns
        -------

        Tmax : float
               The maximum transmittance of the MuellerMatrix
        """
        return self[0,0]+math.sqrt(self[0,1]**2+self[0,2]**2+self[0,3]**2)

    def Tmin(self):
        """
        Return the minimum transmittance.
        
        Returns
        -------

        Tmax : float
               The minimum transmittance of the MuellerMatrix
        """
        return self[0,0]-math.sqrt(self[0,1]**2+self[0,2]**2+self[0,3]**2)

    def diattenuation(self):
        """
        Return the diattenuation.
        
        Returns
        -------

        diattenuation : float
                        The diattenuation
        """
        return math.sqrt(self[0,1]**2+self[0,2]**2+self[0,3]**2)/self[0,0]

    def linear_diattenuation(self):
        """
        Return the linear diattenuation.
        
        Returns
        -------

        diattenuation : float
                        The linear diattenuation
        """
        return math.sqrt(self[0,1]**2+self[0,2]**2)/self[0,0]

    def polarization_dependent_loss(self):
        """
        Return the polarization dependent loss.

        Returns
        -------

        pdl : float
              The polarization dependent loss
        """
        return 10.*math.log(self.Tmax()/self.Tmin())/math.log(10.);

    def polarizance(self):
        """
        Return the polarizance.

        Returns
        -------

        polarizance : float
                      The polarizance
        """
        return math.sqrt(self[1,0]**2+self[2,0]**2+self[3,0]**2)/self[0,0]

    def depolarization_index(self):
        """
        Return the depolarization index.

        Returns
        -------

        PI : float
             The depolarization index
        """
        mmt = (self.T).dot(self)
        return math.sqrt((mmt[0,0]+mmt[1,1]+mmt[2,2]+mmt[3,3]-
                          self[0,0]**2)/(3*self[0,0]**2))

    def extinction_ratio(self):
        """
        Return the extinction ratio.

        Returns
        -------
        
        ER : float
             The extinction ratio
        """
        return self.Tmax()/self.Tmin()

    def rotate(self,angle):
        """
        Return a Mueller matrix rotated by an angle.

        The angle is measured clockwise looking into the beam.

        Parameters
        ----------

        angle : float
                The angle in radians

        Returns
        -------

        M : MuellerMatrix
            Mueller matrix rotated by angle
        """
        sina = math.sin(angle)
        cosa = math.cos(angle)
        R = JonesMueller([[cosa,sina],[-sina,cosa]])
        return R @ self @ R.T

    def parity(self):
        """
        Mueller matrix with a parity conversion.

        Returns
        -------

        M : MuellerMatrix
            Mueller matrix for a parity conversion
        """
        mm = self.copy()
        mm[0,2] = -mm[0,2]
        mm[0,3] = -mm[0,3]
        mm[1,2] = -mm[1,2]
        mm[1,3] = -mm[1,3]
        mm[2,0] = -mm[2,0]
        mm[2,1] = -mm[2,1]
        mm[3,0] = -mm[3,0]
        mm[3,1] = -mm[3,1]
        return mm

    def transpose(self):
        """
        Mueller matrix transpose.

        Returns
        -------

        M : MuellerMatrix
            The transposed Mueller matrix
        """
        return self.T

    def inverse(self):
        """
        Mueller matrix inverse.

        Returns
        -------

        M : MuellerMatrix
            The inverted Mueller matrix
        """
        return np.linalg.inv(self)

    def valid(self):
        """
        Test if a matrix is a Stokes-Stokes mapping matrix.

        Uses the condition described in C.R. Givens and A.B. Kostinski, 
        'A simple necessary and sufficient condition on physically realizable 
        Mueller matrices,' Journal of Modern Optics 40, 471-481 (1993).

        Returns 
        -------
            bool
            True if the matrix will always map a valid Stokes vector to a 
            valid Stokes vector; otherwise, False.
        """
        M = self + 1E-10*self[0,0]*MuellerDiagonal([1,0,0,0])
        G = MinkowskiG
        MM = G @ M.T @ G @ M
        Q, W = np.linalg.eig(MM)
        W = np.array(W.T)
        
        imax=0
        max = -1E308
        ww = W[0]
        a = (Q @ np.conj(Q)) *1E-14
        for q,w in zip(Q,W):
            if abs(q.imag)> a:
                return False
            if (abs(q)>max):
                max = abs(q)
                ww = w.copy()

        ww = ww/ww[0]
        for i in range(1,4):
            if abs(ww[i].imag) > a:
                return False
        if 1 < ww[1].real**2+ww[2].real**2+ww[3].real**2:
            return False
        return True

    def physically_valid(self):
        """ 
        Test if a matrix is physically valid.

        The function physically_valid() calculates the Mueller matrix coherency
        matrix, which is the coherency matrix for the four Jones-Mueller 
        matrices generated by the Pauli matrices. This coherency matrix must be
        a valid coherency matrix (positive semi-definite) if the matrix is the
        convex sum of Jones-Mueller matrices.  This test is a more stringent 
        test than valid() for the validity of the Mueller matrix.

        See B.N. Simon, et al, "A complete characterization of pre-Mueller and 
        Mueller matrices in polarization optics," 
        J. Opt. Soc. Am. A 27(2), 188-199 (2010).

        Returns 
        -------
            bool
            True if the matrix is physically valid, that is, if it is a convex
            sum of Jones-Mueller matrices.
        """
        h = MuellerToHermitian(self)
        Q, W = np.linalg.eig(h)
        a = (Q @ np.conj(Q)) *1E-8
        for q in Q:
            if q.real < -a:
                return False
            if abs(q.imag) > a:
                return False
        return True
    
    def Closest_NonDepolarizing(self):
        """
        Return the largest matrix in the Cloude decomposition, which may be
        interpreted as the closest non-depolarizing Mueller matrix.

        Returns
        -------

        M : MuellerMatrix
            The largest matrix in the Cloude decomposition.
        """
        M = self.Cloude_Decomposition()
        return M[0-rank]

    def normalized(self):
        """
        The normalized Mueller matrix.

        Returns
        -------

        M : MuellerMatrix
            The Mueller matrix normalized so that its [0,0] element is 1.
        """
        return self/self[0,0]

    def entropy(self):
        """
        The polarization entropy.

        Returns
        -------

        entropy : float
                  Returns the polarization entropy defined in S.R. Cloude and 
                  E. Pottier, "Concept of polarization entropy in optical
                  scattering," Opt. Eng. 34(6) 1599-1610 (1995).
        """
        h = MuellerToHermitian(self)
        Q, W =  np.linalg.eig(h)
        result = 0.
        sumQ = np.sum(Q).real
        for q in Q:
            Pi = abs(q)/sumQ
            if Pi>0:
                result += Pi*np.log(Pi)
        return -result/np.log(4.)


class StokesVector(np.ndarray):
    """
    A class for handling Stokes vectors

    Parameters
    ----------
        
    x : 4-element list or array of float or int, optional
        Default is zero Stokes vector
        
    """
    def __new__(cls,*x):
        if (len(x)==0):
            return np.asarray((0,0,0,0)).view(StokesVector)

        if (len(x)==1):
            s = np.asarray(x).view(StokesVector)
            s = s.reshape((4,))
            return s
        
        if (len(x)==4):
            s = np.asarray(x).view(StokesVector)
            return s
        
        raise Exception("Stokes vector must be 4 elements and real")

    def __mul__(a,b):
        """
        Multiply Stokes vector by a Mueller matrix, another Stokes vector 
        (inner product), or a scalar
        """
        if isinstance(b,MuellerMatrix):
            return a.dot(b)
        if isinstance(b,StokesVector):
            return float(a.dot(b))
        if isinstance(b,(float,int)):
            return StokesVector(super()*b)
    
        raise Exception("Stokes vector can only be multiplied by Mueller "
                        "matrix, a Stokes vector, or a scalar")

    __matmul__ = __mul__

    def __truediv__(a,b):
        "Divide Stokes vector by a scalar"
        if isinstance(b,(int,float)):
            return super().__truediv__(b)
        raise Exception("Stokes vector can only be divided by int or float")

    def I(self):
        return self[0]
    def Q(self):
        return self[1]
    def U(self):
        return self[2]
    def V(self):
        return self[3]
    def q(self):
        return self[1]/self[0]
    def u(self):
        return self[2]/self[0]
    def v(self):
        return self[3]/self[0]
    
    def rotate(self,angle):
        """
        Return a Stokes vector rotated by an angle.

        Parameters
        ----------

        angle : float
                Angle to rotate the Stokes vector in radians

        Returns
        -------
        rot : StokesVector
              Stokes vector rotated by angle
        """
        return JonesMueller([[math.cos(angle),math.sin(angle)],
                             [-math.sin(angle),math.cos(angle)]]) @ self

    def eta(self):
        """
        Return the principal angle of the polarization ellipse

        Returns
        -------

        eta : float
              The principal angle of the polarization ellipse
        """
        return math.atan2(self[2],self[1])/2


    def intensity(self):
        """
        Return the intensity.

        Returns
        -------
        
        I : float
            The intensity
        """
        return self[0]

    def DOLP(self):
        """
        Return the degree of linear polarization.

        Returns
        -------

        DOLP : float
               The degree of linear polarization
        """
        return np.sqrt(self[1]**2+self[2]**2)/self[0]

    def DOP(self):
        """
        Return the degree of polarization.

        Returns
        -------

        DOP : float
              The degree of polarization
        """
        return np.sqrt(self[1]**2+self[2]**2+self[3]**2)/self[0]
    
    def DOCP(self):
        """
        Return the degree of circular polarization.

        Returns
        -------

        DOCP : float
               The degree of circular polarization
        """
        return self[3]/self[0]

    def ellipticity(self):
        """
        Return the ellipticity of the field (ratio of minor to major axes).

        Returns
        -------
        
        e : float
            The ellipticity of the field (ratio of the minor to major axes 
            of the polarization ellipse.
        """
        return self[3]/(self[0]+np.sqrt(self[1]**2+self[3]**2))

    def delta(self):
        """
        Return the phase difference between the two components.

        Returns
        -------
        
        delta : float
                The phase difference between the two components (s and p) 
                in radians. Uses the polarized part of the Stokes vector.
        """
        j = self.JonesVector()
        return cmath.phase(j[1])-cmath.phase(j[0])

    def psi(self):
        """
        Return arc tangent of the amplitude ratio between the two components.

        Returns
        -------
        
        psi : float
              The arc tangent in radians of the ratio between the two 
              components of the polarized part of the light. The returned 
              value is zero for p-polarized light.  
        """
        j = self.JonesVector()
        return math.atan(np.abs(j[1])/np.abs(j[0]));

    def eccentricity(self):
        """
        Return the eccentricity.

        Returns
        -------
       
        eccentricity : float
                       The eccentricity
        """
        return np.sqrt(1-self.e()**2)

    def valid(self):
        """
        Checks if the Stokes vector is valid.

        Returns
        -------

        valid : bool
                True if the Stokes vector is valid.
        """
        return self[0]>=math.sqrt(self[1]**2+self[2]**2+self[3]**2)

    def pol_part(self):
        """
        Return the polarized part of the Stokes vector.

        Returns
        -------
        
        s : StokesVector
            The polarized part of the Stokes vector.
        """
        s = self.copy()
        s[0] = self[0]*self.DOP()
        return s;

    def unpol_part(self):
        """
        Return the unpolarized part of the Stokes vector.

        Returns
        -------
        
        s : StokesVector
            The unpolarized part of the Stokes vector.
        """
        return StokesVector(self[0]*(1-self.DOP()),0,0,0)
        return s;
    
    def JonesVector(self):
        """
        Return a Jones vector associated with the polarized part of the 
        Stokes vector.

        Returns
        -------

        j : 2 element np.array
            A Jones vector that, if converted to a Stokes vector, yields
            the polarized part of the Stokes vector.
        """
        xx = self.pol_part()
        arg = 0 if xx[3]==0. and xx[2]==0. else math.atan2(xx[3],xx[2])
        mag = np.sqrt((xx[0]-xx[1])/2.)
        return np.array([np.sqrt((xx[0]+xx[1])/2.),
                         mag*np.cos(arg) + 1j*mag*np.sin(arg)])

    def normalized(self):
        """
        Return normalized Stokes vector.

        Returns
        -------

        s : StokesVector
            The Stokes vector normalized so that its 0 element is 1.
        """
        return self/self[0]

    def __str__(self):
        return ("StokesVector(" +
                str(self[0]) + "," +
                str(self[1]) + "," +
                str(self[2]) + "," +
                str(self[3]) + ")" )

    def __repr__(self):
        return ("StokesVector(" +
                str(self[0]) + "," +
                str(self[1]) + "," +
                str(self[2]) + "," +
                str(self[3]) + ")" )
    
def JonesMueller(jones): 
    """Convert a Jones matrix to a Mueller matrix
    
    Parameters
    ----------

    jones : 2x2 list or array of int, float, or complex
            A Jones matrix

    Returns
    -------
   
    mueller : MuellerMatrix
              The Mueller matrix associated with jones.

    """
    A = np.array(((1,0,0,1),(1,0,0,-1),(0,1,1,0),(0,1j,-1j,0)))
    Ainv = np.array(((0.5,0.5,0,0),
                     (0,0,0.5,-0.5j),
                     (0,0,0.5,0.5j),
                     (0.5,-0.5,0,0)))
    return MuellerMatrix((A @ np.kron(jones,np.conj(jones)) @ Ainv).real)

def JonesStokes(jones):
    """Convert a Jones vector to a Stokes vector
    
    Parameters
    ----------

    jones : 2 element list or array of int, float, or complex

    """
    x = np.array(jones)
    if x.shape==(2,):
        cx = x.conjugate()
        return StokesVector([(x[0]*cx[0]).real+(x[1]*cx[1]).real,
                             (x[0]*cx[0]).real-(x[1]*cx[1]).real,
                             2*(cx[0]*x[1]).real,
                             2*(cx[0]*x[1]).imag])
    raise Exception("Jones vector must be a 2-element vector")
    
def MuellerToHermitian(m):
    """
    Convert a Mueller matrix into its equivalent Hermitian coherency matrix.
    
    Parameters
    ----------
    m : MuellerMatrix
        The Mueller matrix to be converted.

    Returns
    -------

    h : 4x4 complext np.array
        Hermitian coherence matrix 

    """
    h = np.ndarray((4,4),complex)

    h[0,0]=(m[0,0]+m[1,1]+m[0,1]+m[1,0])/2
    h[0,1]=complex(m[0,2]+m[1,2],m[0,3]+m[1,3])/2
    h[0,2]=complex(m[2,0]+m[2,1],-m[3,0]-m[3,1])/2
    h[0,3]=complex(m[2,2]+m[3,3],m[2,3]-m[3,2])/2
    h[1,1]=(m[0,0]-m[1,1]-m[0,1]+m[1,0])/2
    h[1,2]=complex(m[2,2]-m[3,3],-m[2,3]-m[3,2])/2
    h[1,3]=complex(m[2,0]-m[2,1],-m[3,0]+m[3,1])/2
    h[2,2]=(m[0,0]-m[1,1]+m[0,1]-m[1,0])/2
    h[2,3]=complex(m[0,2]-m[1,2],m[0,3]-m[1,3])/2
    h[3,3]=(m[0,0]+m[1,1]-m[0,1]-m[1,0])/2
    h[1,0]=np.conj(h[0,1])
    h[2,0]=np.conj(h[0,2])
    h[3,0]=np.conj(h[0,3])
    h[2,1]=np.conj(h[1,2])
    h[3,1]=np.conj(h[1,3])
    h[3,2]=np.conj(h[2,3])

    return h

def HermitianToMueller(h):
    """
    Convert a Hermitian coherency matrix to its equivalent Mueller matrix.

    Parameters
    ----------
    h : ndarray((4,4),complex)
        A Hermitian coherency matrix

    Returns
    -------
    m : MuellerMatrix
        The Mueller matrix form of the coherency matrix.

    """
    mm = MuellerMatrix()
    mm[0,0] = (h[0,0]+h[1,1]+h[2,2]+h[3,3]).real/2
    mm[0,1] = (h[0,0]-h[1,1]+h[2,2]-h[3,3]).real/2
    mm[0,2] = (h[0,1]+h[1,0]+h[2,3]+h[3,2]).real/2
    mm[0,3] = (h[0,1]-h[1,0]+h[2,3]-h[3,2]).imag/2
    mm[1,0] = (h[0,0]+h[1,1]-h[2,2]-h[3,3]).real/2
    mm[1,1] = (h[0,0]-h[1,1]-h[2,2]+h[3,3]).real/2
    mm[1,2] = (h[0,1]+h[1,0]-h[2,3]-h[3,2]).real/2
    mm[1,3] = (h[0,1]-h[1,0]-h[2,3]+h[3,2]).imag/2
    mm[2,0] = (h[0,2]+h[2,0]+h[1,3]+h[3,1]).real/2
    mm[2,1] = (h[0,2]+h[2,0]-h[1,3]-h[3,1]).real/2
    mm[2,2] = (h[0,3]+h[3,0]+h[1,2]+h[2,1]).real/2
    mm[2,3] = (h[0,3]-h[3,0]-h[1,2]+h[2,1]).imag/2
    mm[3,0] = (h[2,0]-h[0,2]-h[1,3]+h[3,1]).imag/2
    mm[3,1] = (h[2,0]-h[0,2]+h[1,3]-h[3,1]).imag/2
    mm[3,2] = (h[3,0]-h[0,3]+h[2,1]-h[1,2]).imag/2
    mm[3,3] = (h[0,3]+h[3,0]-h[1,2]-h[2,1]).real/2

    return mm


def Cloude_Decomposition(M):
    """
    Return a list of four non-depolarizing Mueller matrices that sum to the 
    Mueller matrix.

    Parameters
    ----------
        M : MuellerMatrix
            Mueller matrix to be decomposed
    
    Returns
    -------
        M1, M2, M3, M4 : MuellerMatrix 
                         The decomposed Jones-Mueller matrices sorted so that
                         M1[0,0] >= M2[0,0] >= M3[0,0] >= M4[0,0]

    """
    h = MuellerToHermitian(M)
    Q, W = np.linalg.eig(h)
    Winv = np.linalg.inv(W)
    
    MM = [HermitianToMueller(Q[i] * np.outer( W[:,i] , Winv[i,:] ))
         for i in range(4)]
    
    for i in range(0,3):
        for j in range(3,i,-1):
            if MM[j-1][0,0]<MM[j][0,0]:
                MM[j-1], MM[j] = MM[j], MM[j-1]
    return MM

MuellerMatrix.Cloude_Decomposition = Cloude_Decomposition

def Lu_Chipman_Decomposition(M):
    """
    Return a depolarizer, a retarder, and a diattenuator whose ordered product is the
    Mueller matrix. 

    See: S.-Y. Lu and R.A. Chipman, "Interpretation of Mueller matrices based on 
    polar decomposition," J. Opt. Soc. Am. A **13**, 1106-1113 (1996).

    Parameters
    ----------
       M : MuellerMatrix
           The matrix to be decomposed

    Returns
    -------
       depolarizer : MuellerMatrix
       retarder : MuellerMatrix
       diattenuator : MuellerMatrix
            M = depolarizer @ retarder @ diattenuator
    """
    diattenuator = MuellerMatrix()
    retarder = MuellerMatrix()
    depolarizer = MuellerMatrix()
    
    # The diattenuator takes the net transmittance of the matrix.
    # See L&C, Eq. (18)
    Tu = M[0,0]

    # If the transmittance is zero, then the decomposition is trivial.
    if (Tu==0):
        return (MuellerMatrix(np.diag([1.,1.,1.,1.])),
                MuellerMatrix(np.diag([1.,1.,1.,1.])),
                MuellerMatrix())

    # From L&C, Eq. (36), Darrow and Parrow are the diattenuation and
    # polarizance, respectively.
    Darrow = M[0,1:4]/Tu 
    Parrow = M[1:4,0]/Tu 
    D = math.sqrt(Darrow[0]**2+Darrow[1]**2+Darrow[2]**2)

    # Thus, the diattenuator is given by L&C, Eq. (18),
    if D>1E-15:
        diattenuator = MuellerMatrix()
        diattenuator[0,0] = 1
        for i in range(3):
            diattenuator[i+1,0] = diattenuator[0,i+1] = Darrow[i]
        dd = math.sqrt(1-D**2)
        for i in range(3):
            for j in range(3):
                 diattenuator[i+1,j+1] = dd if i==j else 0
                 diattenuator[i+1,j+1] += (1-dd)*Darrow[i]*Darrow[j]/D**2

        del dd                    
    else:
        diattenuator = MuellerUnit()

    diattenuator = diattenuator*Tu

    norm = np.linalg.norm
    cross = np.cross
    sqrt = np.sqrt
                                    
    inverse = np.linalg.inv

    # If the matrix is non-depolarizing, then use
    if 1-M.depolarization_index() < 1E-8:
        depolarizer=MuellerUnit()
        # If the diattenuation is unity (or very close!)...
        if 1-D < 1E-8:
            # Get the retardance vector direction...
            Rarrow = cross(Parrow,Darrow)
            normRarrow = norm(Rarrow)
            Rhat = unit(Rarrow) if normRarrow!=0 else perpto(Parrow,[0,0,1])   
            cosR = Parrow*Darrow/norm(Parrow)/norm(Darrow)
            sinR = sqrt(1-cosR**2) if abs(cosR)<1 else 0

            # The retarder is determined from L&C, Eqs. (14) and (15)...
            retarder[0,0]=1
            retarder[0,1]=retarder[0,2]=retarder[0,3]=0
            retarder[1,0]=retarder[2,0]=retarder[3,0]=0
            for i in range(1,4):
                for j in range(1,4):
                    retarder[i,j] = ((cosR if i==j else 0)
                                     + Rhat[i-1]*Rhat[j-1]*(1-cosR))
                    for k in range(1,4):
                        retarder[i,j] += LeviCivita([i,j,k])*Rhat[k-1]*sinR

        else:
            # If the matrix does not have unit diattenutation, then
            # just use L&C, Eq. (35)...
            retarder = M*np.linalg.inv(diattenuator)


        return depolarizer, retarder, diattenuator

    # L&C, Eq. (47)...
    Mprime = M.dot(np.linalg.inv(diattenuator))

    # L&C, Eq. (48)...
    mprime = np.diag([0.,0.,0.])
    for i in range(3):
        for j in range(3):
            mprime[i,j] = Mprime[i+1,j+1]

    # Also from L&C, Eq. (48)...
    PDelta = Mprime[1:4,0]

    # Get the determinant of mprime...
    mprimedeterminant = np.linalg.det(mprime)

    # If mprime is not singular...
    if mprimedeterminant!=0: 
        # To evaluate L&C, Eq. (52), we need the eigenvalues of
        # m'(m')^T ...
        mmT = mprime.dot(mprime.T)

        # Calculate eigenvalues...
        Q, W = np.linalg.eig(mmT)

        # These are the square roots of the eigenvalues
        lambda1 = sqrt(Q[0].real)
        lambda2 = sqrt(Q[1].real)
        lambda3 = sqrt(Q[2].real)

        ident = np.diag([1.,1.,1.])

        # Finally, this is L&C, Eq. (52)...
        mDelta = (
            np.linalg.inv(mmT+(lambda1*lambda2+lambda2
                               *lambda3+lambda3*lambda1)*ident).dot(
            ((lambda1+lambda2+lambda3)*mmT+lambda1*lambda2*lambda3*ident))
            )

        # Use the right sign, using the determinant...
        if mprimedeterminant<0:
            mDelta = -mDelta

        # Thus, the depolarizer is given by L&C, Eq. (48)...
        depolarizer[0,0] = 1
        depolarizer[0,1] = depolarizer[0,2] = depolarizer[0,3] = 0
        depolarizer[1,0] = PDelta[0]
        depolarizer[2,0] = PDelta[1]
        depolarizer[3,0] = PDelta[2]
        for i in range(1,4):
            for j in range(1,4):
                depolarizer[i,j] = mDelta[i-1,j-1]

        # And the retarder is determined by L&C, Eq. (53)...
        retarder = np.linalg.inv(depolarizer)*Mprime

    else:
        # If the matrix is singular, then we have to use the
        # formalism described in L&C, Appendix B...

        # First, perform the singular value decomposition
        # of mprime...
        U, S, V = np.linalg.svd(mprime)

        # The singular values are...
        lambda1 = S[0]
        lambda2 = S[1]
        lambda3 = S[2]


        mR = np.diag([0.,0.,0.])
        # How one determines the retardance depends upon how many of the
        # singular values are zero...
        if abs(lambda2/lambda1)>1E-15:
            # One of them is zero...
            for i in range(3):
                for j in range(3):
                    # L&C, Eq. (B3)...
                    mR[i,j] = (V[i,0]*U[j,0]+V[i,1]*U[j,1]+V[i,2]*U[j,2])
        elif lambda1!=0:
            # Two of them are zero...
            v = V[:,0]
            u = U[:,0]

            pv = perpto(v,[1,0,0])
            ppv = perpto(pv,v)
            pu = perpto(u,[1,0,0])
            ppu = perpto(pu,u)
            # Essentially, L&C, Eq. (B9)...
            mR = np.outer(v,u)+np.outer(pv,pu)+np.outer(ppv,ppu)
        else:
            # The case of all of them zero is trivial...
            mR = np.diag([1.,1.,1.])

        # The full retarder Mueller matrix is...
        for i in range(4):
            for j in range(4):
                if i==0 and j==0: retarder[i,j]=1
                elif i==0 or j==0: retarder[i,j]=0
                else: retarder[i,j] = mR[i-1,j-1]

        # The depolarizer is determined from...
        # (The matrix mR is always invertible.)
        mDelta = mprime.dot(np.linalg.inv(mR))

        depolarizer[0,0] = 1
        depolarizer[0,1] = depolarizer[0,2] = depolarizer[0,3] = 0
        depolarizer[1,0] = PDelta[0];
        depolarizer[2,0] = PDelta[1];
        depolarizer[3,0] = PDelta[2];
        for i in range(1,4):
            for j in range(1,4):
                depolarizer[i,j] = mDelta[i-1,j-1];


    return depolarizer, retarder, diattenuator

MuellerMatrix.Lu_Chipman_Decomposition = Lu_Chipman_Decomposition

class CharacterizedMueller():
    """
    Characterize the Mueller matrix according to 
    
    N. Ghosh, M.F.G. Wood, and I.A. Vitkin,
    'Mueller matrix decomposition for extraction of individual polarization 
    parameters from complex turbid media exhibiting multiple scattering, optical
    activity, and linear birefringence,'' J. Biomedical Opt. 13, 044036 (2008).
    
    A few other parameters are also calculated.
    
    Parameters
    ----------
    M : MuellerMatrix
        The Mueller matrix to be characterized
    
    Attributes
    ----------
    Mdepol : MuellerMatrix
             Depolarizer in Lu-Chipman decomposition
    
    Mret : MuellerMatrix
           Retarder in Lu-Chipman decomposition
    
    Mdiatten : MuellerMatrix
               Diattenuator in Lu-Chipman decomposition
    
    DiattenuationVector : numpy.array
                          The 3-element diattenuation vector
    
    Diattenuation : float
                    Diattenuation
    
    CircularDiattenuation : float
                            Circular diattenuation
    
    LinearDiattenuation : float
                          Linear diattenuation
    
    DiattenuationAngle : float
                         Angle of diattenuation in radians
    
    Polarizance : numpy.array
                  The 3-element polarizance vector
    
    DepolarizationCoefficient : float
                                Depolarization coefficient 
                                (1 = nondepolarizing, 0 = depolarizing)
    
    LinearRetardance : float
                       Linear retardance in radians
                       
    OpticalRotation : float
                      Optical rotation in radians
                      
    Retardance : float
                 Total retardance in radians
    
    RetardanceAngle : float
                      Angle of retardance in radians
                          
    """
    def __init__(self,M):
        Mdepol,Mret,Mdiatten = Lu_Chipman_Decomposition(M)
        
        self.Mdepol = Mdepol
        self.Mret = Mret
        self.Mdiatten = Mdiatten
        
        self.DiattenuationVector = np.array(Mdiatten[0,1:4])
        self.Diattenuation = math.sqrt(sum(self.DiattenuationVector**2))
        self.CircularDiattenuation = self.DiattenuationVector[2] 
        self.LinearDiattenuation = math.sqrt(sum(np.array(Mdiatten[0,1:3])**2))
        self.DiattenuationAngle = math.atan2(Mdiatten[2,0],Mdiatten[1,0])/2
        
        self.PolarizanceVector = np.array(Mdepol[1:4,0])
        self.Polarizance = math.sqrt(sum(self.PolarizanceVector**2))
        self.DepolarizationCoefficient = abs(np.trace(Mdepol)-1)/3

        self.LinearRetardance = math.acos(math.sqrt((Mret[1,1]+Mret[2,2])**2+
                                                    (Mret[2,1]-Mret[1,2])**2)-1)
        self.OpticalRotation = math.atan((Mret[2,1]-Mret[1,2])/
                                         (Mret[1,1]+Mret[2,2]))/2
        self.Retardance = math.acos(np.trace(Mret)/2-1)
        psi = self.OpticalRotation
        Mcircret = MuellerMatrix([[1,0,0,0],
                                  [0,math.cos(2*psi),math.sin(2*psi),0],
                                  [0,-math.sin(2*psi),math.cos(2*psi),0],
                                  [0,0,0,1]])
        Mlinret = Mret @ Mcircret.inverse()
        self.RetardanceAngle = math.atan2((Mlinret[3,1]-Mlinret[1,3]),
                                     (Mlinret[2,3]-Mlinret[3,2]))/2

    def __str__(self):
        s = ("Mdepol = \n" + str(self.Mdepol)
             + "\nMret = \n" + str(self.Mret)
             + "\nMdiatten = \n" + str(self.Mdiatten)
             + "\nDiattenuationVector = " + str(self.DiattenuationVector)
             + "\nDiattenuation = " + str(self.Diattenuation) 
             + "\nCircularDiattenuation = " + str(self.CircularDiattenuation)
             + "\nLinearDiattenuation = " + str(self.LinearDiattenuation)
             + "\nDiattenuationAngle = " + str(self.DiattenuationAngle)
             + " rad (" + str(self.DiattenuationAngle/deg) + " deg)"
             + "\nPolarizanceVector = " + str(self.PolarizanceVector)
             + "\nPolarizance = " + str(self.Polarizance)
             + "\nDepolarizationCoefficient = "
             + str(self.DepolarizationCoefficient)
             + "\nLinearRetardance = " + str(self.LinearRetardance)
             + " rad (" + str(self.LinearRetardance/deg) + " deg)"
             + "\nOpticalRotation = " + str(self.OpticalRotation)
             + " rad (" + str(self.OpticalRotation/deg) + " deg)"
             + "\nRetardance = " + str(self.Retardance)
             + " rad (" + str(self.Retardance/deg) + "deg)"
             + "\nRetardanceAngle = " + str(self.RetardanceAngle)
             + " rad (" + str(self.RetardanceAngle/deg) + " deg)")
        return s
    

def Reverse_Lu_Chipman_Decomposition(M):
    """
    Return a diattenuator, a retarder, and a depolarizer whose ordered product is the
    Mueller matrix. 

    See: S.-Y. Lu and R.A. Chipman, "Interpretation of Mueller matrices based on 
    polar decomposition," J. Opt. Soc. Am. A **13**, 1106-1113 (1996).

    Parameters
    ----------
       M : MuellerMatrix
           The matrix to be decomposed

    Returns
    -------
       diattenuator : MuellerMatrix
       retarder : MuellerMatrix
       depolarizer : MuellerMatrix
             M = diattenuator @ retarder @ depolarizer

    """
    depolarizer, retarder, diattenuator = Lu_Chipman_Decomposition(M.T)
    return diattenuator.T, retarder.T, depolarizer.T

MuellerMatrix.Reverse_Lu_Chipman_Decomposition = (
                                    Reverse_Lu_Chipman_Decomposition)

def Symmetric_Decomposition(M):
    """
    Perform the symmetric decomposition described by Ossikovski.

    See: R. Ossikovski, J. Opt. Soc. Am. A **26**, 1109-1118 (2009).

    Parameters
    ----------
       M : MuellerMatrix

    Returns
    -------
       diatten2 : MuellerMatrix
       ret2 : MuellerMatrix
       depol : MuellerMatrix
       ret1 : MuellerMatrix
       diatten1 : MuellerMatrix
             M = diatten2 @ ret2 @ depol @ ret1 @ diatten1
    """
    def diattenFromVec(D):
        """
        Construct the diattenuator from the diattenuation veotor.
        """
        MD = np.zeros((4,4))
        MD[0,0] = 1
        for i in range(1,4):
            MD[0,i] = D[i]
            MD[i,0] = D[i]
        Dhat = np.array([D[1],D[2],D[3]])
        Dsqr = Dhat @ Dhat
        Dhat = Dhat/np.sqrt(Dsqr)
        if Dsqr<1:
            mD = (np.sqrt(1-Dsqr)*np.diag([1,1,1])+
                  (1-np.sqrt(1-Dsqr))*np.outer(Dhat,Dhat))
        else:
            mD = np.outer(Dhat,Dhat)

        for i in range(3):
            for j in range(3):
                MD[i+1,j+1] = mD[i,j]
        return MD

    def retFromSubMatrix(m):
        """
        Constructs the retarder from the 3x3 submatrix
        """
        MR = np.zeros((4,4))
        MR[0,0]=1
        for i in range(3):
            for j in range(3):
                MR[i+1,j+1] = m[i,j]
        return MR

    def FindValid(Q,W):
        """
        Given the eigenvalues Q and eigenvector matrix W, returns
        the eigenvector associated with the valid Stokes vector
        """
        W = W.T

        imax=0
        max = -1E308
        ww = W[0]
        a = (Q @ np.conj(Q)) *1E-14
        for q,w in zip(Q,W):
            if abs(q.imag)> a:
                raise ValueError("Not a valid Mueller matrix")
            if (abs(q)>max):
                max = abs(q)
                ww = w.copy()

        if ww[0].real**2 < ww[1].real**2+ww[2].real**2+ww[3].real**2:
            raise ValueError("Not a valid Mueller matrix")

        return ww

    # Step 1 
    M = np.array(M) + np.diag([1,0,0,0])*M[0,0]*1E-10
    G = np.array(MinkowskiG)
    MTGMG = M.T @ G @ M @ G
    Q, W = np.linalg.eig(MTGMG)

    S1 = FindValid(Q,W)
    S1 = S1/S1[0]
    S2 = M @ G @ S1
    S2 = S2/S2[0]

    # Step 2
    diatten1 = diattenFromVec(S1)
    diatten2 = diattenFromVec(S2)

    # Step 3
    Mprime = np.linalg.inv(diatten2) @ M @ np.linalg.inv(diatten1)

    # Step 4
    mprime = Mprime[1:4,1:4]
    mR2, mdepol, mR1T = np.linalg.svd(mprime)

    ret1 = retFromSubMatrix(mR1T)
    ret2 = retFromSubMatrix(mR2)

    depol = np.diag([Mprime[0,0],mdepol[0],mdepol[1],mdepol[2]])

    return (MuellerMatrix(diatten2),
            MuellerMatrix(ret2),
            MuellerMatrix(depol),
            MuellerMatrix(ret1),
            MuellerMatrix(diatten1))

MuellerMatrix.Symmetric_Decomposition = Symmetric_Decomposition

def MuellerLog(M):
    """
    Return matrix logarithm (logarithmic decomposition).

    Returns
    -------
    4x4 numpy.array 
        The matrix logarithm of M 
    """
    Q, W = np.linalg.eig(np.array(M))
    QQ = []
    for q in Q:
        if q<=0:
            QQ.append(np.log(q*(1+1E-16j)))
        else:
            QQ.append(np.log(q))  
    L = W @ np.diag(QQ) @ np.linalg.inv(W)
    isreal = True
    for el in L.flatten():
        if abs(el.imag)>1E-14:
            isreal = False
    if isreal: return L.real
    else: return L

MuellerMatrix.MuellerLog = MuellerLog

def MuellerExp(L):
    """
    Return matrix exponential.

    Returns
    -------

    L : 4x4 np.array 
        The matrix exponent of M
    """
    Q, W = np.linalg.eig(np.array(L))
    Q = np.exp(Q)
    M = W @ np.diag(Q) @ np.linalg.inv(W)
    isreal = True
    for m in M.flatten():
        if abs(m.imag)>1E-14:
            isreal = False
    if isreal: return M.real
    else: return M

MuellerMatrix.MuellerExp = MuellerExp

MinkowskiG = MuellerMatrix(np.diag([1,-1,-1,-1]))

def LeviCivita(ijk):
    """
    Return the Levi-Civita symbol
    Argument: 3-element list of int
    """
    result = 1
    if max(ijk)-min(ijk)!=len(ijk)-1: return 0
    for i in range(0,len(ijk)-1):
        for j in range(len(ijk)-1,i,-1):
            if ijk[j-1]==ijk[j]: return 0
            if ijk[j-1]>ijk[j]:
                ijk[j-1], ijk[j] = ijk[j], ijk[j-1]
                result = -result
    return result;

PauliMatrices = [np.array([[1,0],
                           [0,1]]),
                 np.array([[1,0],
                           [0,-1]]),
                 np.array([[0,1],
                           [1,0]]),
                 np.array([[0,-1j],
                           [1j,0]])]

def JonesZero():
    """
    Return a zero Jones matrix
    """
    return np.array([[0,0],[0,0]])

def JonesUnit():
    """
    Return a unit Jones matrix
    """
    return np.array([[1,0],[0,1]])

def JonesRotator(angle):
    """
    Return a Jones matrix that rotates the polarization by angle"

    Parameters
    ----------

    angle : float
            Angle that matrix will rotate in radians

    Returns
    -------

    R : 2x2 np.array
        The rotation matrix 
    """
    cosa = np.cos(angle)
    sina = np.sin(angle)
    return np.array([[cosa,sina],[-sina,cosa]])

def MuellerRotator(angle):
    """
    Return a Mueller matrix that rotates the polarization by angle"

    Parameters
    ----------

    angle : float
            Angle that matrix will rotate in radians

    Returns
    -------

    R : MuellerMatrix
        The rotation matrix 
    """
    cos2a = np.cos(2*angle)
    sin2a = np.sin(2*angle)
    return MuellerMatrix(((1,0,0,0),
                          (0,cos2a,sin2a,0),
                          (0,-sin2a,cos2a,0),
                          (0,0,0,1)))

def JonesBasis(angle=0, DOCP=None, jones=None):
    """
    Return the basis vectors for a Jones Matrix using a Jones vector or a 
    degree of circular polarization and an angle. The returned vectors are 
    orthogonal to one another and of unit intensity. This function
    is primarily used by JonesRetarder and JonesDiattenuator
    
    Parameters
    ----------
    
    angle : float
            Angle (in radians) of axis (default is 0)
    
    DOCP : float
           Degree of circular polarization of the axis.
           Cannot define both DOCP and jones. 
    
    jones : 2-element list of complex
            Jones vector associated with axis
            Cannot define both DOCP and jones. (Default is [1,0])
    
    Returns
    -------
    
    jones1, jones2 : 2-element lists of complex
                     Two Jones vectors that are mutually orthogonal 
                     and of unit intensity
    """
    if jones == None:
        try:
            jones1 = np.array([1,(1/DOCP-np.sqrt(1/DOCP**2-1))*1j])
        except:
            jones1 = np.array([1,0])
    else:
        if DOCP == None:
            jones1 = jones/np.sqrt(jones @ np.conj(jones))
        else:
            raise Error("DOCP and jones both defined")
    
    jones2 = np.array([np.conj(jones1[1]),-np.conj(jones1[0])])
    jones1 = jones1/np.sqrt(jones1 @ np.conj(jones1))
    jones2 = jones2/np.sqrt(jones2 @ np.conj(jones2))

    R = JonesRotator(angle)
    return (jones1 @ R,jones2 @ R)

def JonesRetarder(phase=None, angle=0, DOCP=None, jones=None):
    """
    Return the Jones matrix for a homogeneous retarder with a given phase 
    and axes. 
    
    Parameters
    ----------
    
    angle : float
            Angle (in radians) of axis (default is 0)
    
    DOCP : float
           Degree of circular polarization of the axis.
           Cannot define both DOCP and jones. 
    
    jones : 2-element list of complex
            Jones vector associated with axis
            Cannot define both DOCP and jones. (Default is [1,0])
    
    Returns
    -------
    
    jones : 2x2 numpy array of complex
            Jones matrix retarder
    
    """
    j1, j2 = JonesBasis(angle, DOCP, jones)
    return np.outer(j1,np.conj(j1)) + np.outer(j2,np.conj(j2))*np.exp(1j*phase)

        
def JonesDiattenuator(diatten=None, angle=0, DOCP=None, jones=None):
    """
    Return the Jones matrix for a homogeneous diattenuator with a given 
    phase and axes. 
    
    Parameters
    ----------
    
    angle : float
            Angle (in radians) of axis (default is 0)
    
    DOCP : float
           Degree of circular polarization of the axis.
           Cannot define both DOCP and jones. 
    
    jones : 2-element list of complex
            Jones vector associated with axis
            Cannot define both DOCP and jones. (Default is [1,0])
    
    Returns
    -------
    
    jones : 2x2 numpy array of complex
            Jones matrix diattenuator
    
    """
    j1, j2 = JonesBasis(angle, DOCP, jones)
    return np.outer(j1,np.conj(j1)) + np.outer(j2,np.conj(j2))*np.sqrt(diatten)

def MuellerRetarder(phase=0, angle=0, DOCP=None, jones=None):
    """
    Return the Mueller matrix for a homogeneous retarder with a given 
    phase and axis. 
    
    Parameters
    ----------
    
    angle : float
            Angle (in radians) of axis (default is 0)
    
    DOCP : float
           Degree of circular polarization of the axis.
           Cannot define both DOCP and jones. 
    
    jones : 2-element list of complex
            Jones vector associated with axis
            Cannot define both DOCP and jones. (Default is [1,0])
    
    Returns
    -------
    
    jones : 2x2 numpy array of complex
            Jones matrix retarder
    
    """
    return JonesMueller(JonesRetarder(phase,angle,DOCP,jones))

def MuellerDiattenuator(diatten=0, angle=0, DOCP=None, jones=None):
    """
    Return the Mueller matrix for a homogeneous diattenuator with a 
    given phase and axis. 
    
    Parameters
    ----------
    
    angle : float
            Angle (in radians) of axis (default is 0)
    
    DOCP : float
           Degree of circular polarization of the axis.
           Cannot define both DOCP and jones. 
    
    jones : 2-element list of complex
            Jones vector associated with axis
            Cannot define both DOCP and jones. (Default is [1,0])
    
    Returns
    -------
    
    MuellerMatrix diattenuator
    
    """    
    return JonesMueller(JonesDiattenuator(diatten,angle,DOCP,jones))

JonesPolarizer = JonesDiattenuator
MuellerPolarizer = MuellerDiattenuator

def MuellerZero():
    """
    Return a zero Mueller matrix

    Returns
    -------

    m : MuellerMatrix
        Zero matrix
    """
    return MuellerMatrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])

def MuellerUnit(attenuation=1):
    """
    Return a unit Mueller matrix

    Arguments
    --------
    attenuation : float, optional
                  Default is 1
    """
    a = attenuation
    return MuellerMatrix([[a,0,0,0],[0,a,0,0],[0,0,a,0],[0,0,0,a]])

def MuellerDiagonal(a):
    """
    Return a diagonal Mueller matrix

    Arguments
    ---------
    a : 4 element list or array of float 
    """
    return MuellerMatrix([[a[0],0,0,0],[0,a[1],0,0],[0,0,a[2],0],[0,0,0,a[3]]])

def StokesZero():
    """
    Return zero Stokes vector

    """
    return StokesVector([0,0,0,0])
    
def Polarization(*args, **kwargs):
    """
    Return a StokesVector with a given name
   
    Parameters
    ----------

    args : str or float
           Can be any of:
               `P` or `p` for p-polarization
               `S` or `s` for s-polarization
               `L` or `l` for LCP-polarization
               `R` or `r` for RCP-polarization
               `U` or `u` for unpolarized
               any float for linear polarization at angle `name` in radians
    
    state : str
            Any of the aboves strings that can be used in args

    I : float, optional
        The intensity. (Default is 1)
    
    DOP : float, optional
          The degree of polarization. (Default is 1)

    angle : float, optional
            The principal angle of the polarization. (Default is set by
            args or state.)

    ellipticity : float
                  The ellipticity of the polarization

    Returns
    -------

    s : StokesVector
        The Stokes vector with the given parameters
    """
    state = None
    sensitivity = False
    if 'sensitivity' in kwargs:
        sensitivity = kwargs.pop('sensitivity')

    try:
        if len(args)!=0:
            s = StokesVector(*args)
            if sensitivity:
                return s.pol_part()/2+s.unpol_part()
            else:
                return s
    except Exception:
        pass

    if len(args)==1:
        state = args[0]
    
    I = 1
    p = 1
    chi = 0
    
    if 'state' in kwargs:
        state = kwargs.pop('state')
    if 'I' in kwargs:
        I = kwargs.pop('I')
    if 'DOP' in kwargs:
        p = kwargs.pop('DOP')
    if 'angle' in kwargs:
        psi = kwargs.pop('angle')
    if 'ellipticity' in kwargs:
        chi = kwargs.pop('ellipticity')
    if len(kwargs) != 0:
        raise Exception("Unrecognized keyword")

    if sensitivity:
        s0 = (1-p/2)
        s1 = p/2
        p = 1
    else:
        s0 = 1
        s1 = 1
        
    if state == "P" or state == 'p':
        return StokesVector(I*s0,-I*p*s1,0,0)
    if state == "S" or state == 's':
        return StokesVector(I*s0,I*p*s1,0,0)
    if state == "U" or state == 'u':
        return StokesVector(I,0,0,0)
    if state == "L" or state == 'l':
        return StokesVector(I*s0,0,0,I*p*s1)
    if state == "R" or state == 'r':
        return StokesVector(I*s0,0,0,-I*p*s1)
    
    return StokesVector(I*s0,
                        I*p*s1*math.cos(2*psi)*math.cos(2*chi),
                        I*p*s1*math.sin(2*psi)*math.cos(2*chi),
                        I*p*s1*math.sin(2*chi))


def Sensitivity(*args, **kwargs):
    """
    Return a StokesVector with a given name that 
    corresponds to a sensitivity.  
   
    Parameters
    ----------

    name : str or float
           Can be any of:
               `P` or `p` for p-polarization
               `S` or `s` for s-polarization
               `L` or `l` for LCP-polarization
               `R` or `r` for RCP-polarization
               `U` or `u` for unpolarized
               any float for linear polarization at angle `name`
    """
    if 'sensitivity' not in kwargs:
        kwargs['sensitivity'] = True

    return Polarization(*args, **kwargs)

    
