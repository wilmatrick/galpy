###############################################################################
#   ChandrasekharDynamicalFriction: Class that implements the Chandrasekhar 
#                                   dynamical friction
###############################################################################
import numpy as nu
from scipy import special
from Potential import Potential
_INVSQRTTWO= 1./nu.sqrt(2.)
_INVSQRTPI= 1./nu.sqrt(nu.pi)
def isothermalrhor(r):
    return 1./r**2.
def isothermalsigmar(r):
    return _INVSQRTTWO
class ChandrasekharDynamicalFriction(Potential):
    """Class that implements the Chandrasekhar dynamical friction"""
    def __init__(self,amp=1.,lnLambda=2.,ms=.1,rhor=isothermalrhor,
                 sigmar=isothermalsigmar):
        """
        NAME:
           __init__
        PURPOSE:
           initialize a Chandrasekhar Dynamical Friction force
        INPUT:
           amp - amplitude to be applied to the potential (default: 1)
           ms - satellite mass
           lnLambda - Coulomb integral
           rhor - function that gives the density as a function of r
           sigmar - function that gives the velocity dispersion as a function 
                    of r
        OUTPUT:
           (none)
        HISTORY:
           2011-12-26 - Started - Bovy (NYU)
        """
        Potential.__init__(self,amp=amp)
        self._lnLambda= lnLambda
        self._ms= ms
        self._rhor= rhor
        self._sigmar= sigmar
        return None

    def _evaluate(self,R,z,phi=0.,t=0.,dR=0,dphi=0,v=None):
        """
        NAME:
           _evaluate
        PURPOSE:
           evaluate the potential at R,z
        INPUT:
           R - Galactocentric cylindrical radius
           z - vertical height
           phi - azimuth
           t - time
           dR=, dphi=
           v= current velocity as [vR,vphi,vz]
        OUTPUT:
           Phi(R,z)
        HISTORY:
           2010-04-02 - Started - Bovy (NYU)
           2010-04-30 - Adapted for R,z - Bovy (NYU)
        """
        if dR == 0 and dphi == 0:
            raise NotImplementedError("'_evaluate' not implemented for ChandrasekharDynamicalFriction")
        elif dR == 1 and dphi == 1:
            return self._Rforce(R,z,phi=phi,t=t,v=v)
        elif dR == 0 and dphi == 1:
            return self._phiforce(R,z,phi=phi,t=t,v=v)
        else:
            raise NotImplementedError("'_evaluate' not implemented for ChandrasekharDynamicalFriction")
            pass

    def _Rforce(self,R,z,phi=0.,t=0.,v=None):
        """
        NAME:
           _Rforce
        PURPOSE:
           evaluate the radial force for this potential
        INPUT:
           R - Galactocentric cylindrical radius
           z - vertical height
           phi - azimuth
           t - time
           v= current velocity
        OUTPUT:
           the radial force
        HISTORY:
        """
        vT= v[1]*R
        vs= nu.sqrt(v[0]**2.+vT+v[2]**2.)
        r= nu.sqrt(R**2.+z**2.)
        X= vs*_INVSQRTTWO/self._sigmar(r)
        Xfactor= special.erf(X)-2.*X*_INVSQRTPI*nu.exp(-X**2.)
        force= -self._ms*self._rhor(r)/vs**3.*Xfactor
        return force*v[0]

    def _phiforce(self,R,z,phi=0.,t=0.,v=None):
        """
        NAME:
           _phiforce
        PURPOSE:
           evaluate the azimuthal force for this potential
        INPUT:
           R - Galactocentric cylindrical radius
           z - vertical height
           phi - azimuth
           t - time
           v= current velocity
        OUTPUT:
           the azimuthal force
        HISTORY:
        """
        vT= v[1]*R
        vs= nu.sqrt(v[0]**2.+vT+v[2]**2.)
        r= nu.sqrt(R**2.+z**2.)
        X= vs*_INVSQRTTWO/self._sigmar(r)
        Xfactor= special.erf(X)-2.*X*_INVSQRTPI*nu.exp(-X**2.)
        force= -self._ms*self._rhor(r)/vs**3.*Xfactor
        return force*v[1]

    def _zforce(self,R,z,phi=0.,t=0.,v=None):
        """
        NAME:
           _zforce
        PURPOSE:
           evaluate the vertical force for this potential
        INPUT:
           R - Galactocentric cylindrical radius
           z - vertical height
           phi - azimuth
           t - time
           v= current velocity
        OUTPUT:
           the vertical force
        HISTORY:
        """
        vT= v[1]*R
        vs= nu.sqrt(v[0]**2.+vT+v[2]**2.)
        r= nu.sqrt(R**2.+z**2.)
        X= vs*_INVSQRTTWO/self._sigmar(r)
        Xfactor= special.erf(X)-2.*X*_INVSQRTPI*nu.exp(-X**2.)
        force= -self._ms*self._rhor(r)/vs**3.*Xfactor
        return force*v[2]

