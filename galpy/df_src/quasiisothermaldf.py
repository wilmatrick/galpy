#A 'Binney' quasi-isothermal DF
import math
import warnings
import numpy
from scipy import optimize, interpolate, integrate
from galpy import potential
from galpy import actionAngle
from galpy.actionAngle import actionAngleIsochrone
from galpy.potential import IsochronePotential
from galpy.orbit import Orbit
from galpy.util import galpyWarning
_NSIGMA=4
_DEFAULTNGL=10
_DEFAULTNGL2=20
class quasiisothermaldf:
    """Class that represents a 'Binney' quasi-isothermal DF"""
    def __init__(self,hr,sr,sz,hsr,hsz,pot=None,aA=None,
                 cutcounter=False,
                 _precomputerg=True,_precomputergrmax=None,
                 _precomputergnLz=51,
                 ro=1.,lo=10./220./8.):
        """
        NAME:

           __init__

        PURPOSE:

           Initialize a quasi-isothermal DF

        INPUT:

           hr - radial scale length

           sr - radial velocity dispersion at the solar radius

           sz - vertical velocity dispersion at the solar radius

           hsr - radial-velocity-dispersion scale length

           hsz - vertial-velocity-dispersion scale length

           pot= Potential instance or list thereof

           aA= actionAngle instance used to convert (x,v) to actions

           cutcounter= if True, set counter-rotating stars' DF to zero

           ro= reference radius for surface mass and sigmas

           lo= reference angular momentum below where there are significant numbers of retrograde stars

        OTHER INPUTS:

           _precomputerg= if True (default), pre-compute the rL(L)

           _precomputergrmax= if set, this is the maximum R for which to pre-compute rg (default: 5*hr)

           _precomputergnLz if set, number of Lz to pre-compute rg for (default: 51)

        OUTPUT:

           object

        HISTORY:

           2012-07-25 - Started - Bovy (IAS@MPIA)

        """
        self._hr= hr
        self._sr= sr
        self._sz= sz
        self._hsr= hsr
        self._hsz= hsz
        self._ro= ro
        self._lo= lo
        self._lnsr= math.log(self._sr)
        self._lnsz= math.log(self._sz)
        if pot is None:
            raise IOError("pot= must be set")
        self._pot= pot
        if aA is None:
            raise IOError("aA= must be set")
        self._aA= aA
        if not self._aA._pot == self._pot:
            if not isinstance(self._aA,actionAngleIsochrone):
                raise IOError("Potential in aA does not appear to be the same as given potential pot")
            elif isinstance(self._pot,IsochronePotential) and \
                    not self._aA.b == self._pot.b and \
                    not self._aA.amp == self._pot._amp:
                raise IOError("Potential in aA does not appear to be the same as given potential pot")
        self._cutcounter= cutcounter
        if _precomputerg:
            if _precomputergrmax is None:
                _precomputergrmax= 5*self._hr
            self._precomputergrmax= _precomputergrmax
            self._precomputergnLz= _precomputergnLz
            self._precomputergLzmin= 0.01
            self._precomputergLzmax= self._precomputergrmax\
                *potential.vcirc(self._pot,self._precomputergrmax)
            self._precomputergLzgrid= numpy.linspace(self._precomputergLzmin,self._precomputergLzmax,self._precomputergnLz)
            self._rls= numpy.array([potential.rl(self._pot,l) for l in self._precomputergLzgrid])
            #Spline interpolate
            self._rgInterp= interpolate.InterpolatedUnivariateSpline(self._precomputergLzgrid,self._rls,k=3)
        else:
            self._precomputergrmax= 0.
            self._rgInterp= None
            self._rls= None
            self._precomputergnr= None
            self._precomputergLzgrid= None
            self._precomputergLzmin= \
                numpy.finfo(numpy.dtype(numpy.float64)).max
            self._precomputergLzmax= \
                numpy.finfo(numpy.dtype(numpy.float64)).min
        self._precomputerg= _precomputerg
        self._glxdef, self._glwdef= \
            numpy.polynomial.legendre.leggauss(_DEFAULTNGL)
        self._glxdef2, self._glwdef2= \
            numpy.polynomial.legendre.leggauss(_DEFAULTNGL2)
        self._glxdef12, self._glwdef12= \
            numpy.polynomial.legendre.leggauss(_DEFAULTNGL/2)
        return None

    def __call__(self,*args,**kwargs):
        """
        NAME:
           __call__
        PURPOSE:
           return the DF
        INPUT:
           Either:
              a)(jr,lz,jz) tuple
                 where:
                    jr - radial action
                    lz - z-component of angular momentum
                    jz - vertical action
              b) R,vR,vT,z,vz
              c) Orbit instance: initial condition used if that's it, orbit(t)
                 if there is a time given as well

           log= if True, return the natural log

           +scipy.integrate.quadrature kwargs

           func= function of (jr,lz,jz) to multiply f with (useful for moments)

        OUTPUT:
           value of DF
        HISTORY:
           2012-07-25 - Written - Bovy (IAS@MPIA)
        NOTE:
           For Miyamoto-Nagai/adiabatic approximation this seems to take 
           about 30 ms / evaluation in the extended Solar neighborhood
           For a MWPotential/adiabatic approximation this takes about 
           50 ms / evaluation in the extended Solar neighborhood

           For adiabatic-approximation grid this seems to take 
           about 0.67 to 0.75 ms / evaluation in the extended Solar 
           neighborhood (includes some out of the grid)

           up to 200x faster when called with vector R,vR,vT,z,vz
        """
        #First parse log
        if kwargs.has_key('log'):
            log= kwargs['log']
            kwargs.pop('log')
        else:
            log= False
        if kwargs.has_key('_return_actions'):
            _return_actions= kwargs['_return_actions']
            kwargs.pop('_return_actions')
        else:
            _return_actions= False
        if kwargs.has_key('_return_freqs'):
            _return_freqs= kwargs['_return_freqs']
            kwargs.pop('_return_freqs')
        else:
            _return_freqs= False
        if kwargs.has_key('rg'):
            thisrg= kwargs['rg']
            kwargs.pop('rg')
            kappa= kwargs['kappa']
            kwargs.pop('kappa')
            nu= kwargs['nu']
            kwargs.pop('nu')
            Omega= kwargs['Omega']
            kwargs.pop('Omega')
        else:
            thisrg= None
            kappa= None
            nu= None
            Omega= None
        #First parse args
        if len(args) == 1 and not isinstance(args[0],Orbit): #(jr,lz,jz)
            jr,lz,jz= args[0]
        else:
            #Use self._aA to calculate the actions
            try:
                jr,lz,jz= self._aA(*args,**kwargs)
            except actionAngle.UnboundError:
                if log: return -numpy.finfo(numpy.dtype(numpy.float64)).max
                else: return 0.
            #if isinstance(jr,(list,numpy.ndarray)) and len(jr) > 1: jr= jr[0]
            #if isinstance(jz,(list,numpy.ndarray)) and len(jz) > 1: jz= jz[0]
        if not isinstance(lz,numpy.ndarray) and self._cutcounter and lz < 0.:
            if log: return -numpy.finfo(numpy.dtype(numpy.float64)).max
            else: return 0.
        #First calculate rg
        if thisrg is None:
            thisrg= self.rg(lz)
            #Then calculate the epicycle and vertical frequencies
            kappa, nu= self._calc_epifreq(thisrg), self._calc_verticalfreq(thisrg)
            Omega= numpy.fabs(lz)/thisrg/thisrg
        #calculate surface-densities and sigmas
        lnsurfmass= (self._ro-thisrg)/self._hr
        lnsr= self._lnsr+(self._ro-thisrg)/self._hsr
        lnsz= self._lnsz+(self._ro-thisrg)/self._hsz
        #Calculate func
        if kwargs.has_key('func'):
            if log:
                funcTerm= numpy.log(kwargs['func'](jr,lz,jz))
            else:
                funcFactor= kwargs['func'](jr,lz,jz)
        #Calculate fsr
        else:
            if log:
                funcTerm= 0.
            else:
                funcFactor= 1.            
        if log:
            lnfsr= numpy.log(Omega)+lnsurfmass-2.*lnsr-math.log(math.pi)\
                -numpy.log(kappa)\
                +numpy.log(1.+numpy.tanh(lz/self._lo))\
                -kappa*jr*numpy.exp(-2.*lnsr)
            lnfsz= numpy.log(nu)-math.log(2.*math.pi)\
                -2.*lnsz-nu*jz*numpy.exp(-2.*lnsz)
            out= lnfsr+lnfsz+funcTerm
            if isinstance(lz,numpy.ndarray):
                out[numpy.isnan(out)]= -numpy.finfo(numpy.dtype(numpy.float64)).max
                if self._cutcounter: out[(lz < 0.)]= -numpy.finfo(numpy.dtype(numpy.float64)).max
            elif numpy.isnan(out): out= -numpy.finfo(numpy.dtype(numpy.float64)).max
        else:
            srm2= numpy.exp(-2.*lnsr)
            fsr= Omega*numpy.exp(lnsurfmass)*srm2/math.pi/kappa\
                *(1.+numpy.tanh(lz/self._lo))\
                *numpy.exp(-kappa*jr*srm2)
            szm2= numpy.exp(-2.*lnsz)
            fsz= nu/2./math.pi*szm2*numpy.exp(-nu*jz*szm2)
            out= fsr*fsz*funcFactor
            if isinstance(lz,numpy.ndarray):
                out[numpy.isnan(out)]= 0.
                if self._cutcounter: out[(lz < 0.)]= 0.
            elif numpy.isnan(out): out= 0.
        if _return_actions and _return_freqs:
            return (out,jr,lz,jz,thisrg,kappa,nu,Omega)
        elif _return_actions:
            return (out,jr,lz,jz)
        elif _return_freqs:
            return (out,thisrg,kappa,nu,Omega)
        else:
            return out

    def estimate_hr(self,R,z=0.,dR=10.**-8.,**kwargs):
        """
        NAME:
           estimate_hr
        PURPOSE:
           estimate the exponential scale length at R
        INPUT:
           R - Galactocentric radius

           z= height (default: 0 pc)

           dR- range in R to use

           density kwargs
        OUTPUT:
           estimated hR
        HISTORY:
           2012-09-11 - Written - Bovy (IAS)
           2013-01-28 - Re-written - Bovy
        """
        Rs= [R-dR/2.,R+dR/2.]
        if z is None:
            sf= numpy.array([self.surfacemass_z(r,**kwargs) for r in Rs])
        else:
            sf= numpy.array([self.density(r,z,**kwargs) for r in Rs])
        lsf= numpy.log(sf)
        return -dR/(lsf[1]-lsf[0])

    def estimate_hz(self,R,z,dz=10.**-8.,**kwargs):
        """
        NAME:
           estimate_hz
        PURPOSE:
           estimate the exponential scale height at R
        INPUT:
           R - Galactocentric radius

           dz - z range to use

           density kwargs
        OUTPUT:
           estimated hz
        HISTORY:
           2012-08-30 - Written - Bovy (IAS)

           2013-01-28 - Re-written - Bovy
        """
        if z == 0.:
            zs= [z,z+dz]
        else:
            zs= [z-dz/2.,z+dz/2.]
        sf= numpy.array([self.density(R,zz,**kwargs) for zz in zs])
        lsf= numpy.log(sf)
        return -dz/(lsf[1]-lsf[0])

    def estimate_hsr(self,R,z=0.,dR=10.**-8.,**kwargs):
        """
        NAME:
           estimate_hsr
        PURPOSE:
           estimate the exponential scale length of the radial dispersion at R
        INPUT:
           R - Galactocentric radius

           z= height (default: 0 pc)

           dR- range in R to use

           density kwargs
        OUTPUT:
           estimated hsR
        HISTORY:
           2013-03-08 - Written - Bovy (IAS)
        """
        Rs= [R-dR/2.,R+dR/2.]
        sf= numpy.array([self.sigmaR2(r,z,**kwargs) for r in Rs])
        lsf= numpy.log(sf)/2.
        return -dR/(lsf[1]-lsf[0])

    def estimate_hsz(self,R,z=0.,dR=10.**-8.,**kwargs):
        """
        NAME:
           estimate_hsz
        PURPOSE:
           estimate the exponential scale length of the vertical dispersion at R
        INPUT:
           R - Galactocentric radius

           z= height (default: 0 pc)

           dR- range in R to use

           density kwargs
        OUTPUT:
           estimated hsz
        HISTORY:
           2013-03-08 - Written - Bovy (IAS)
        """
        Rs= [R-dR/2.,R+dR/2.]
        sf= numpy.array([self.sigmaz2(r,z,**kwargs) for r in Rs])
        lsf= numpy.log(sf)/2.
        return -dR/(lsf[1]-lsf[0])

    def surfacemass_z(self,R,nz=7,zmax=1.,fixed_quad=True,fixed_order=8,
                      **kwargs):
        """
        NAME:
           surfacemass_z
        PURPOSE:
           calculate the vertically-integrated surface density
        INPUT:
           R - Galactocentric radius

           fixed_quad= if True (default), use Gauss-Legendre integration

           fixed_order= (20), order of GL integration to use

           nz= number of zs to use to estimate

           zmax=m minimum z to use

           density kwargs
        OUTPUT:
           \Sigma(R)
        HISTORY:
           2012-08-30 - Written - Bovy (IAS)
        """
        if fixed_quad:
            return 2.*integrate.fixed_quad(lambda x: self.density(R*numpy.ones(fixed_order),x),
                                           0.,.5,n=fixed_order)[0]
        zs= numpy.linspace(0.,zmax,nz)
        sf= numpy.array([self.density(R,z,**kwargs) for z in zs])
        lsf= numpy.log(sf)
        #Interpolate
        lsfInterp= interpolate.UnivariateSpline(zs,
                                                lsf,
                                                k=3)
        #Integrate
        return 2.*integrate.quad((lambda x: numpy.exp(lsfInterp(x))),
                                 0.,1.)[0]

    def vmomentdensity(self,R,z,n,m,o,nsigma=None,mc=False,nmc=10000,
                       _returnmc=False,_vrs=None,_vts=None,_vzs=None,
                       _rawgausssamples=False,
                       gl=False,ngl=_DEFAULTNGL,_returngl=False,_glqeval=None,
                       _return_actions=False,_jr=None,_lz=None,_jz=None,
                       _return_freqs=False,
                       _rg=None,_kappa=None,_nu=None,_Omega=None,
                       _sigmaR1=None,_sigmaz1=None,
                       **kwargs):
        """
        NAME:
           vmomentdensity
        PURPOSE:
           calculate the an arbitrary moment of the velocity distribution 
           at R times the density
        INPUT:
           R - radius at which to calculate the moment(/ro)

           n - vR^n

           m - vT^m

           o - vz^o

        OPTIONAL INPUT:
           nsigma - number of sigma to integrate the velocities over (when doing explicit numerical integral)

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

           gl= use Gauss-Legendre

           _returngl= if True, return the evaluated DF

           _return_actions= if True, return the evaluated actions (does not work with _returngl currently)

           _return_freqs= if True, return the evaluated frequencies and rg (does not work with _returngl currently)

        OUTPUT:
           <vR^n vT^m  x density> at R,z
        HISTORY:
           2012-08-06 - Written - Bovy (IAS@MPIA)
        """
        if isinstance(R,numpy.ndarray):
            return numpy.array([self.vmomentdensity(r,zz,n,m,o,nsigma=nsigma,
                                                    mc=mc,nmc=nmc,
                                                    gl=gl,ngl=ngl,**kwargs) for r,zz in zip(R,z)])
        if isinstance(self._aA,(actionAngle.actionAngleAdiabatic,
                                actionAngle.actionAngleAdiabaticGrid)):
            if n % 2 == 1. or o % 2 == 1.:
                return 0. #we know this must be the case
        if nsigma == None:
            nsigma= _NSIGMA
        if _sigmaR1 is None:
            sigmaR1= self._sr*numpy.exp((self._ro-R)/self._hsr)
        else:
            sigmaR1= _sigmaR1
        if _sigmaz1 is None:
            sigmaz1= self._sz*numpy.exp((self._ro-R)/self._hsz)
        else:
            sigmaz1= _sigmaz1
        thisvc= potential.vcirc(self._pot,R)
        #Use the asymmetric drift equation to estimate va
        gamma= numpy.sqrt(0.5)
        va= sigmaR1**2./2./thisvc\
            *(gamma**2.-1. #Assume close to flat rotation curve, sigphi2/sigR2 =~ 0.5
               +R*(1./self._hr+2./self._hsr))
        if math.fabs(va) > sigmaR1: va = 0.#To avoid craziness near the center
        if gl:
            if ngl % 2 == 1:
                raise ValueError("ngl must be even")
            if not _glqeval is None and ngl != _glqeval.shape[0]:
                _glqeval= None
            #Use Gauss-Legendre integration for all
            if ngl == _DEFAULTNGL:
                glx, glw= self._glxdef, self._glwdef
                glx12, glw12= self._glxdef12, self._glwdef12
            elif ngl == _DEFAULTNGL2:
                glx, glw= self._glxdef2, self._glwdef2
                glx12, glw12= self._glxdef, self._glwdef
            else:
                glx, glw= numpy.polynomial.legendre.leggauss(ngl)
                glx12, glw12= numpy.polynomial.legendre.leggauss(ngl/2)
            #Evaluate everywhere
            if isinstance(self._aA,(actionAngle.actionAngleAdiabatic,
                                    actionAngle.actionAngleAdiabaticGrid)):
                vRgl= 4.*sigmaR1/2.*(glx+1.)
                vzgl= 4.*sigmaz1/2.*(glx+1.)
                vRglw= glw
                vzglw= glw
            else:
                vRgl= 4.*sigmaR1/2.*(glx12+1.)
                #vRgl= 1.5/2.*(glx12+1.)
                vRgl= list(vRgl)
                vRgl.extend(-4.*sigmaR1/2.*(glx12+1.))
                #vRgl.extend(-1.5/2.*(glx12+1.))
                vRgl= numpy.array(vRgl)
                vzgl= 4.*sigmaz1/2.*(glx12+1.)
                #vzgl= 1.5/2.*(glx12+1.)
                vzgl= list(vzgl)
                vzgl.extend(-4.*sigmaz1/2.*(glx12+1.))
                #vzgl.extend(-1.5/2.*(glx12+1.))
                vzgl= numpy.array(vzgl)
                vRglw= glw12
                vRglw= list(vRglw)
                vRglw.extend(glw12)
                vRglw= numpy.array(vRglw)
                vzglw= glw12
                vzglw= list(vzglw)
                vzglw.extend(glw12)
                vzglw= numpy.array(vzglw)
            vTgl= 1.5/2.*(glx+1.)
            #Tile everything
            vTgl= numpy.tile(vTgl,(ngl,ngl,1)).T
            vRgl= numpy.tile(numpy.reshape(vRgl,(1,ngl)).T,(ngl,1,ngl))
            vzgl= numpy.tile(vzgl,(ngl,ngl,1))
            vTglw= numpy.tile(glw,(ngl,ngl,1)).T #also tile weights
            vRglw= numpy.tile(numpy.reshape(vRglw,(1,ngl)).T,(ngl,1,ngl))
            vzglw= numpy.tile(vzglw,(ngl,ngl,1))
            #evaluate
            if _glqeval is None and _jr is None:
                logqeval, jr, lz, jz, rg, kappa, nu, Omega= self(R+numpy.zeros(ngl*ngl*ngl),
                                           vRgl.flatten(),
                                           vTgl.flatten(),
                                           z+numpy.zeros(ngl*ngl*ngl),
                                           vzgl.flatten(),
                                           log=True,
                                           _return_actions=True,
                                                                 _return_freqs=True)
                logqeval= numpy.reshape(logqeval,(ngl,ngl,ngl))
            elif not _jr is None and _rg is None:
                logqeval, jr, lz, jz, rg, kappa, nu, Omega= self((_jr,_lz,_jz),
                                                                 log=True,
                                                                 _return_actions=True,
                                                                 _return_freqs=True)
                logqeval= numpy.reshape(logqeval,(ngl,ngl,ngl))
            elif not _jr is None and not _rg is None:
                logqeval, jr, lz, jz, rg, kappa, nu, Omega= self((_jr,_lz,_jz),
                                           rg=_rg,kappa=_kappa,nu=_nu,
                                           Omega=_Omega,
                                           log=True,
                                           _return_actions=True,
                                                                 _return_freqs=True)
                logqeval= numpy.reshape(logqeval,(ngl,ngl,ngl))
            else:
                logqeval= _glqeval
            if _returngl:
                return (numpy.sum(numpy.exp(logqeval)*vRgl**n*vTgl**m*vzgl**o
                                  *vTglw*vRglw*vzglw)*sigmaR1*sigmaz1*3.,
                        logqeval)
            elif _return_actions and _return_freqs:
                return (numpy.sum(numpy.exp(logqeval)*vRgl**n*vTgl**m*vzgl**o
                                  *vTglw*vRglw*vzglw)*sigmaR1*sigmaz1*3.,
                        jr,lz,jz,
                        rg,kappa,nu,Omega)
            elif _return_actions:
                return (numpy.sum(numpy.exp(logqeval)*vRgl**n*vTgl**m*vzgl**o
                                  *vTglw*vRglw*vzglw)*sigmaR1*sigmaz1*3.,
                        jr,lz,jz)
            else:
                return numpy.sum(numpy.exp(logqeval)*vRgl**n*vTgl**m*vzgl**o
                                 *vTglw*vRglw*vzglw*sigmaR1*sigmaz1*3.)
        elif mc:
            mvT= (thisvc-va)/gamma/sigmaR1
            if _vrs is None:
                vrs= numpy.random.normal(size=nmc)
            else:
                vrs= _vrs
            if _vts is None:
                vts= numpy.random.normal(size=nmc)+mvT
            else:
                if _rawgausssamples:
                    vts= _vts+mvT
                else:
                    vts= _vts
            if _vzs is None:
                vzs= numpy.random.normal(size=nmc)
            else:
                vzs= _vzs
            Is= _vmomentsurfaceMCIntegrand(vzs,vrs,vts,numpy.ones(nmc)*R,
                                           numpy.ones(nmc)*z,
                                           self,sigmaR1,gamma,sigmaz1,mvT,
                                           n,m,o)
            if _returnmc:
                if _rawgausssamples:
                    return (numpy.mean(Is)*sigmaR1**(2.+n+m)*gamma**(1.+m)*sigmaz1**(1.+o),
                        vrs,vts-mvT,vzs)
                else:
                    return (numpy.mean(Is)*sigmaR1**(2.+n+m)*gamma**(1.+m)*sigmaz1**(1.+o),
                            vrs,vts,vzs)
            else:
                return numpy.mean(Is)*sigmaR1**(2.+n+m)*gamma**(1.+m)*sigmaz1**(1.+o)
        else: #pragma: no cover because this is too slow; a warning is shown
            warnings.warn("Calculations using direct numerical integration using tplquad is not recommended and extremely slow; it has also not been carefully tested",galpyWarning)
            return integrate.tplquad(_vmomentsurfaceIntegrand,
                                     1./gamma*(thisvc-va)/sigmaR1-nsigma,
                                     1./gamma*(thisvc-va)/sigmaR1+nsigma,
                                     lambda x: 0., lambda x: nsigma,
                                     lambda x,y: 0., lambda x,y: nsigma,
                                     (R,z,self,sigmaR1,gamma,sigmaz1,n,m,o),
                                     **kwargs)[0]*sigmaR1**(2.+n+m)*gamma**(1.+m)*sigmaz1**(1.+o)
        
    def jmomentdensity(self,R,z,n,m,o,nsigma=None,mc=True,nmc=10000,
                       _returnmc=False,_vrs=None,_vts=None,_vzs=None,
                       **kwargs):
        """
        NAME:
           jmomentdensity
        PURPOSE:
           calculate the an arbitrary moment of an action
           of the velocity distribution 
           at R times the surfacmass
        INPUT:
           R - radius at which to calculate the moment(/ro)

           n - jr^n

           m - lz^m

           o - jz^o

        OPTIONAL INPUT:
           nsigma - number of sigma to integrate the velocities over (when doing explicit numerical integral)

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

        OUTPUT:
           <jr^n lz^m jz^o  x density> at R
        HISTORY:
           2012-08-09 - Written - Bovy (IAS@MPIA)
        """
        if nsigma == None:
            nsigma= _NSIGMA
        sigmaR1= self._sr*numpy.exp((self._ro-R)/self._hsr)
        sigmaz1= self._sz*numpy.exp((self._ro-R)/self._hsz)
        thisvc= potential.vcirc(self._pot,R)
        #Use the asymmetric drift equation to estimate va
        gamma= numpy.sqrt(0.5)
        va= sigmaR1**2./2./thisvc\
            *(gamma**2.-1. #Assume close to flat rotation curve, sigphi2/sigR2 =~ 0.5
               +R*(1./self._hr+2./self._hsr))
        if math.fabs(va) > sigmaR1: va = 0.#To avoid craziness near the center
        if mc:
            mvT= (thisvc-va)/gamma/sigmaR1
            if _vrs is None:
                vrs= numpy.random.normal(size=nmc)
            else:
                vrs= _vrs
            if _vts is None:
                vts= numpy.random.normal(size=nmc)+mvT
            else:
                vts= _vts
            if _vzs is None:
                vzs= numpy.random.normal(size=nmc)
            else:
                vzs= _vzs
            Is= _jmomentsurfaceMCIntegrand(vzs,vrs,vts,numpy.ones(nmc)*R,numpy.ones(nmc)*z,self,sigmaR1,gamma,sigmaz1,mvT,n,m,o)
            if _returnmc:
                return (numpy.mean(Is)*sigmaR1**2.*gamma*sigmaz1,
                        vrs,vts,vzs)
            else:
                return numpy.mean(Is)*sigmaR1**2.*gamma*sigmaz1
        else: #pragma: no cover because this is too slow; a warning is shown
            warnings.warn("Calculations using direct numerical integration using tplquad is not recommended and extremely slow; it has also not been carefully tested",galpyWarning)
            return integrate.tplquad(_jmomentsurfaceIntegrand,
                                     1./gamma*(thisvc-va)/sigmaR1-nsigma,
                                     1./gamma*(thisvc-va)/sigmaR1+nsigma,
                                     lambda x: 0., lambda x: nsigma,
                                     lambda x,y: 0., lambda x,y: nsigma,
                                     (R,z,self,sigmaR1,gamma,sigmaz1,n,m,o),
                                     **kwargs)[0]*sigmaR1**2.*gamma*sigmaz1
        
    def density(self,R,z,nsigma=None,mc=False,nmc=10000,
                gl=True,ngl=_DEFAULTNGL,**kwargs):
        """
        NAME:
           density
        PURPOSE:
           calculate the density at R,z by marginalizing over velocity
        INPUT:

           R - radius at which to calculate the density

           z - height at which to calculate the density

        OPTIONAL INPUT:
           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

           gl= if True, calculate using Gauss-Legendre integration

           ngl= if gl, use ngl-th order Gauss-Legendre integration for each dimension
        OUTPUT:
           density at (R,z)
        HISTORY:
           2012-07-26 - Written - Bovy (IAS@MPIA)
        """
        return self.vmomentdensity(R,z,0.,0.,0.,
                                   nsigma=nsigma,mc=mc,nmc=nmc,
                                   gl=gl,ngl=ngl,
                                   **kwargs)
    
    def sigmaR2(self,R,z,nsigma=None,mc=False,nmc=10000,
                gl=True,ngl=_DEFAULTNGL,**kwargs):
        """
        NAME:
           sigmaR2
        PURPOSE:
           calculate sigma_R^2 by marginalizing over velocity
        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this
        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

           gl= if True, calculate using Gauss-Legendre integration

           ngl= if gl, use ngl-th order Gauss-Legendre integration for each dimension

        OUTPUT:
           sigma_R^2
        HISTORY:
           2012-07-30 - Written - Bovy (IAS@MPIA)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            return self.vmomentdensity(R,z,2.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                                             **kwargs)/surfmass
        elif gl:
            surfmass, glqeval= self.vmomentdensity(R,z,0.,0.,0.,
                                                       gl=gl,ngl=ngl,
                                                       _returngl=True,
                                                       **kwargs)
            return self.vmomentdensity(R,z,2.,0.,0.,
                                           ngl=ngl,gl=gl,
                                           _glqeval=glqeval,
                                           **kwargs)/surfmass
        else: #pragma: no cover because this is too slow; a warning is shown
            return (self.vmomentdensity(R,z,2.,0.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/
                    self.vmomentdensity(R,z,0.,0.,0.,
                                            nsigma=nsigma,mc=mc,nmc=nmc,
                                            **kwargs))
        
    def sigmaRz(self,R,z,nsigma=None,mc=False,nmc=10000,
                gl=True,ngl=_DEFAULTNGL,**kwargs):
        """
        NAME:
           sigmaRz
        PURPOSE:
           calculate sigma_RZ^2 by marginalizing over velocity
        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this

        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

           gl= if True, calculate using Gauss-Legendre integration

           ngl= if gl, use ngl-th order Gauss-Legendre integration for each dimension

        OUTPUT:
           sigma_Rz^2
        HISTORY:
           2012-07-30 - Written - Bovy (IAS@MPIA)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            return self.vmomentdensity(R,z,1.,0.,1.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                                             **kwargs)/surfmass
        elif gl:
            surfmass, glqeval= self.vmomentdensity(R,z,0.,0.,0.,
                                                       gl=gl,ngl=ngl,
                                                       _returngl=True,
                                                       **kwargs)
            return self.vmomentdensity(R,z,1.,0.,1.,
                                           ngl=ngl,gl=gl,
                                           _glqeval=glqeval,
                                           **kwargs)/surfmass
        else: #pragma: no cover because this is too slow; a warning is shown
            return (self.vmomentdensity(R,z,1.,0.,1.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/
                    self.vmomentdensity(R,z,0.,0.,0.,
                                            nsigma=nsigma,mc=mc,nmc=nmc,
                                            **kwargs))
        
    def tilt(self,R,z,nsigma=None,mc=False,nmc=10000,
             gl=True,ngl=_DEFAULTNGL,**kwargs):
        """
        NAME:
           tilt
        PURPOSE:
           calculate the tilt of the velocity ellipsoid by marginalizing over velocity
        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this
        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

           gl= if True, calculate using Gauss-Legendre integration

           ngl= if gl, use ngl-th order Gauss-Legendre integration for each dimension

        OUTPUT:
           tilt in degree
        HISTORY:
           2012-12-23 - Written - Bovy (IAS)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            tsigmar2= self.vmomentdensity(R,z,2.,0.,0.,
                                              nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                              _vrs=vrs,_vts=vts,_vzs=vzs,
                                              **kwargs)/surfmass
            tsigmaz2= self.vmomentdensity(R,z,0.,0.,2.,
                                              nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                              _vrs=vrs,_vts=vts,_vzs=vzs,
                                              **kwargs)/surfmass
            tsigmarz= self.vmomentdensity(R,z,1.,0.,1.,
                                              nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                              _vrs=vrs,_vts=vts,_vzs=vzs,
                                              **kwargs)/surfmass
            return 0.5*numpy.arctan(2.*tsigmarz/(tsigmar2-tsigmaz2))/numpy.pi*180.
        elif gl:
            surfmass, glqeval= self.vmomentdensity(R,z,0.,0.,0.,
                                                       gl=gl,ngl=ngl,
                                                       _returngl=True,
                                                       **kwargs)
            tsigmar2= self.vmomentdensity(R,z,2.,0.,0.,
                                             ngl=ngl,gl=gl,
                                             _glqeval=glqeval,
                                             **kwargs)/surfmass
            tsigmaz2= self.vmomentdensity(R,z,0.,0.,2.,
                                              ngl=ngl,gl=gl,
                                              _glqeval=glqeval,
                                              **kwargs)/surfmass
            tsigmarz= self.vmomentdensity(R,z,1.,0.,1.,
                                              ngl=ngl,gl=gl,
                                              _glqeval=glqeval,
                                              **kwargs)/surfmass
            return 0.5*numpy.arctan(2.*tsigmarz/(tsigmar2-tsigmaz2))/numpy.pi*180.
        else:
            raise NotImplementedError("Use either mc=True or gl=True")
        
    def sigmaz2(self,R,z,nsigma=None,mc=False,nmc=10000,
                gl=True,ngl=_DEFAULTNGL,**kwargs):
        """
        NAME:
           sigmaz2
        PURPOSE:
           calculate sigma_z^2 by marginalizing over velocity
        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this

        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

           gl= if True, calculate using Gauss-Legendre integration

           ngl= if gl, use ngl-th order Gauss-Legendre integration for each dimension

        OUTPUT:
           sigma_z^2
        HISTORY:
           2012-07-30 - Written - Bovy (IAS@MPIA)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            return self.vmomentdensity(R,z,0.,0.,2.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                                             **kwargs)/surfmass
        elif gl:
            surfmass, glqeval= self.vmomentdensity(R,z,0.,0.,0.,
                                                       gl=gl,ngl=ngl,
                                                       _returngl=True,
                                                       **kwargs)
            return self.vmomentdensity(R,z,0.,0.,2.,
                                           ngl=ngl,gl=gl,
                                           _glqeval=glqeval,
                                           **kwargs)/surfmass
        else: #pragma: no cover because this is too slow; a warning is shown
            return (self.vmomentdensity(R,z,0.,0.,2.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/
                    self.vmomentdensity(R,z,0.,0.,0.,
                                            nsigma=nsigma,mc=mc,nmc=nmc,
                                            **kwargs))
        
    def meanvT(self,R,z,nsigma=None,mc=False,nmc=10000,
                gl=True,ngl=_DEFAULTNGL,**kwargs):
        """
        NAME:

           meanvT

        PURPOSE:

           calculate the mean rotational velocity by marginalizing over velocity 

        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this

        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

           gl= if True, calculate using Gauss-Legendre integration

           ngl= if gl, use ngl-th order Gauss-Legendre integration for each dimension

        OUTPUT:
           meanvT
        HISTORY:
           2012-07-30 - Written - Bovy (IAS@MPIA)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            return self.vmomentdensity(R,z,0.,1.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                                             **kwargs)/surfmass
        elif gl:
            surfmass, glqeval= self.vmomentdensity(R,z,0.,0.,0.,
                                                       gl=gl,ngl=ngl,
                                                       _returngl=True,
                                                       **kwargs)
            return self.vmomentdensity(R,z,0.,1.,0.,
                                           ngl=ngl,gl=gl,
                                           _glqeval=glqeval,
                                           **kwargs)/surfmass
        else: #pragma: no cover because this is too slow; a warning is shown
            return (self.vmomentdensity(R,z,0.,1.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/
                    self.vmomentdensity(R,z,0.,0.,0.,
                                            nsigma=nsigma,mc=mc,nmc=nmc,
                                            **kwargs))
        
    def meanvR(self,R,z,nsigma=None,mc=False,nmc=10000,
               gl=True,ngl=_DEFAULTNGL,**kwargs):
        """
        NAME:
           meanvR
        PURPOSE:
           calculate the mean radial velocity by marginalizing over velocity
        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this

        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

           gl= if True, calculate using Gauss-Legendre integration

           ngl= if gl, use ngl-th order Gauss-Legendre integration for each dimension

        OUTPUT:
           meanvR
        HISTORY:
           2012-12-23 - Written - Bovy (IAS)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            return self.vmomentdensity(R,z,1.,0.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                                             **kwargs)/surfmass
        elif gl:
            surfmass, glqeval= self.vmomentdensity(R,z,0.,0.,0.,
                                                       gl=gl,ngl=ngl,
                                                       _returngl=True,
                                                       **kwargs)
            return self.vmomentdensity(R,z,1.,0.,0.,
                                           ngl=ngl,gl=gl,
                                           _glqeval=glqeval,
                                           **kwargs)/surfmass
        else: #pragma: no cover because this is too slow; a warning is shown
            return (self.vmomentdensity(R,z,1.,0.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/
                    self.vmomentdensity(R,z,0.,0.,0.,
                                            nsigma=nsigma,mc=mc,nmc=nmc,
                                            **kwargs))
        
    def meanvz(self,R,z,nsigma=None,mc=False,nmc=10000,
               gl=True,ngl=_DEFAULTNGL,**kwargs):
        """
        NAME:
           meanvz
        PURPOSE:
           calculate the mean vertical velocity by marginalizing over velocity
        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this

        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

           gl= if True, calculate using Gauss-Legendre integration

           ngl= if gl, use ngl-th order Gauss-Legendre integration for each dimension

        OUTPUT:
           meanvz
        HISTORY:
           2012-12-23 - Written - Bovy (IAS)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            return self.vmomentdensity(R,z,0.,0.,1.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                                             **kwargs)/surfmass
        elif gl:
            surfmass, glqeval= self.vmomentdensity(R,z,0.,0.,0.,
                                                       gl=gl,ngl=ngl,
                                                       _returngl=True,
                                                       **kwargs)
            return self.vmomentdensity(R,z,0.,0.,1.,
                                           ngl=ngl,gl=gl,
                                           _glqeval=glqeval,
                                           **kwargs)/surfmass
        else: #pragma: no cover because this is too slow; a warning is shown
            return (self.vmomentdensity(R,z,0.,0.,1.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/
                    self.vmomentdensity(R,z,0.,0.,0.,
                                            nsigma=nsigma,mc=mc,nmc=nmc,
                                            **kwargs))
        
    def sigmaT2(self,R,z,nsigma=None,mc=False,nmc=10000,
                gl=True,ngl=_DEFAULTNGL,**kwargs):
        """
        NAME:
           sigmaT2
        PURPOSE:
           calculate sigma_T^2 by marginalizing over velocity
        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this

        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

           gl= if True, calculate using Gauss-Legendre integration

           ngl= if gl, use ngl-th order Gauss-Legendre integration for each dimension

        OUTPUT:
           sigma_T^2
        HISTORY:
           2012-07-30 - Written - Bovy (IAS@MPIA)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            mvt= self.vmomentdensity(R,z,0.,1.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                                             **kwargs)/surfmass
            return self.vmomentdensity(R,z,0.,2.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                           **kwargs)/surfmass\
                                           -mvt**2.
        elif gl:
            surfmass, glqeval= self.vmomentdensity(R,z,0.,0.,0.,
                                                       gl=gl,ngl=ngl,
                                                       _returngl=True,
                                                       **kwargs)
            mvt= self.vmomentdensity(R,z,0.,1.,0.,
                                         ngl=ngl,gl=gl,
                                         _glqeval=glqeval,
                                         **kwargs)/surfmass
            return self.vmomentdensity(R,z,0.,2.,0.,
                                           ngl=ngl,gl=gl,
                                           _glqeval=glqeval,
                                           **kwargs)/surfmass-mvt**2.

        else: #pragma: no cover because this is too slow; a warning is shown
            surfmass= self.vmomentdensity(R,z,0.,0.,0.,
                                              nsigma=nsigma,mc=mc,nmc=nmc,
                                              **kwargs)
            return (self.vmomentdensity(R,z,0.,2.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/surfmass\
                        -(self.vmomentdensity(R,z,0.,2.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/surfmass)**2.)

    def meanjr(self,R,z,nsigma=None,mc=True,nmc=10000,**kwargs):
        """
        NAME:
           meanjr
        PURPOSE:
           calculate the mean radial action by marginalizing over velocity
        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this

        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

        OUTPUT:
           meanjr
        HISTORY:
           2012-08-09 - Written - Bovy (IAS@MPIA)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            return self.jmomentdensity(R,z,1.,0.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                                             **kwargs)/surfmass
        else: #pragma: no cover because this is too slow; a warning is shown
            return (self.jmomentdensity(R,z,1.,0.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/
                    self.vmomentdensity(R,z,0.,0.,0.,
                                            nsigma=nsigma,mc=mc,nmc=nmc,
                                            **kwargs))
        
    def meanlz(self,R,z,nsigma=None,mc=True,nmc=10000,**kwargs):
        """
        NAME:
           meanlz
        PURPOSE:
           calculate the mean angular momemtum by marginalizing over velocity
        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this

        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

        OUTPUT:
           meanlz
        HISTORY:
           2012-08-09 - Written - Bovy (IAS@MPIA)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            return self.jmomentdensity(R,z,0.,1.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                                             **kwargs)/surfmass
        else: #pragma: no cover because this is too slow; a warning is shown
            return (self.jmomentdensity(R,z,0.,1.,0.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/
                    self.vmomentdensity(R,z,0.,0.,0.,
                                            nsigma=nsigma,mc=mc,nmc=nmc,
                                            **kwargs))
        
    def meanjz(self,R,z,nsigma=None,mc=True,nmc=10000,**kwargs):
        """
        NAME:
           meanjz
        PURPOSE:
           calculate the mean vertical action by marginalizing over velocity
        INPUT:

           R - radius at which to calculate this

           z - height at which to calculate this

        OPTIONAL INPUT:

           nsigma - number of sigma to integrate the velocities over

           scipy.integrate.tplquad kwargs epsabs and epsrel

           mc= if True, calculate using Monte Carlo integration

           nmc= if mc, use nmc samples

        OUTPUT:
           meanjz
        HISTORY:
           2012-08-09 - Written - Bovy (IAS@MPIA)
        """
        if mc:
            surfmass, vrs, vts, vzs= self.vmomentdensity(R,z,0.,0.,0.,
                                                             nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=True,
                                                             **kwargs)
            return self.jmomentdensity(R,z,0.,0.,1.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,_returnmc=False,
                                           _vrs=vrs,_vts=vts,_vzs=vzs,
                                                             **kwargs)/surfmass
        else: #pragma: no cover because this is too slow; a warning is shown
            return (self.jmomentdensity(R,z,0.,0.,1.,
                                           nsigma=nsigma,mc=mc,nmc=nmc,
                                           **kwargs)/
                    self.vmomentdensity(R,z,0.,0.,0.,
                                            nsigma=nsigma,mc=mc,nmc=nmc,
                                            **kwargs))
        
    def sampleV(self,R,z,n=1):
        """
        NAME:
           sampleV
        PURPOSE:
           sample a radial, azimuthal, and vertical velocity at R,z
        INPUT:

           R - Galactocentric distance

           z - height

           n= number of distances to sample

        OUTPUT:
           list of samples
        HISTORY:
           2012-12-17 - Written - Bovy (IAS)
        """
        #Determine the maximum of the velocity distribution
        maxVR= 0.
        maxVz= 0.
        maxVT= optimize.fmin_powell((lambda x: -self(R,0.,x,z,0.,log=True)),
                                    1.)
        logmaxVD= self(R,maxVR,maxVT,z,maxVz,log=True)
        #Now rejection-sample
        vRs= []
        vTs= []
        vzs= []
        while len(vRs) < n:
            nmore= n-len(vRs)+1
            #sample
            propvR= numpy.random.normal(size=nmore)*2.*self._sr
            propvT= numpy.random.normal(size=nmore)*2.*self._sr+maxVT
            propvz= numpy.random.normal(size=nmore)*2.*self._sz
            VDatprop= self(R+numpy.zeros(nmore),
                           propvR,propvT,z+numpy.zeros(nmore),
                           propvz,log=True)-logmaxVD
            VDatprop-= -0.5*(propvR**2./4./self._sr**2.+propvz**2./4./self._sz**2.\
                                 +(propvT-maxVT)**2./4./self._sr**2.)
            VDatprop= numpy.reshape(VDatprop,(nmore))
            indx= (VDatprop > numpy.log(numpy.random.random(size=nmore))) #accept
            vRs.extend(list(propvR[indx]))
            vTs.extend(list(propvT[indx]))
            vzs.extend(list(propvz[indx]))
        out= numpy.empty((n,3))
        out[:,0]= vRs[0:n]
        out[:,1]= vTs[0:n]
        out[:,2]= vzs[0:n]
        return out

    def pvR(self,vR,R,z,gl=True,ngl=_DEFAULTNGL2):
        """
        NAME:
           pvR
        PURPOSE:
           calculate the marginalized vR probability at this location (NOT normalized by the density)
        INPUT:

           vR - radial velocity (/vo)

           R - radius (/ro)

           z - height (/ro)

           gl - use Gauss-Legendre integration (True, currently the only option)

           ngl - order of Gauss-Legendre integration

        OUTPUT:
           p(vR,R,z)
        HISTORY:
           2012-12-22 - Written - Bovy (IAS)
        """
        sigmaz1= self._sz*numpy.exp((self._ro-R)/self._hsz)
        if gl:
            if ngl % 2 == 1:
                raise ValueError("ngl must be even")
            #Use Gauss-Legendre integration for all
            if ngl == _DEFAULTNGL:
                glx, glw= self._glxdef, self._glwdef
                glx12, glw12= self._glxdef12, self._glwdef12
            elif ngl == _DEFAULTNGL2:
                glx, glw= self._glxdef2, self._glwdef2
                glx12, glw12= self._glxdef, self._glwdef
            else:
                glx, glw= numpy.polynomial.legendre.leggauss(ngl)
                glx12, glw12= numpy.polynomial.legendre.leggauss(ngl/2)
            #Evaluate everywhere
            if isinstance(self._aA,(actionAngle.actionAngleAdiabatic,
                                    actionAngle.actionAngleAdiabaticGrid)):
                vzgl= 4.*sigmaz1/2.*(glx+1.)
                vzglw= glw
            else:
                vzgl= 4.*sigmaz1/2.*(glx12+1.)
                vzgl= list(vzgl)
                vzgl.extend(-4.*sigmaz1/2.*(glx12+1.))
                vzgl= numpy.array(vzgl)
                vzglw= glw12
                vzglw= list(vzglw)
                vzglw.extend(glw12)
                vzglw= numpy.array(vzglw)
            vTgl= 1.5/2.*(glx+1.)
            #Tile everything
            vTgl= numpy.tile(vTgl,(ngl,1)).T
            vzgl= numpy.tile(vzgl,(ngl,1))
            vTglw= numpy.tile(glw,(ngl,1)).T #also tile weights
            vzglw= numpy.tile(vzglw,(ngl,1))
            #evaluate
            logqeval= numpy.reshape(self(R+numpy.zeros(ngl*ngl),
                                         vR+numpy.zeros(ngl*ngl),
                                         vTgl.flatten(),
                                         z+numpy.zeros(ngl*ngl),
                                         vzgl.flatten(),
                                         log=True),
                                    (ngl,ngl))
            return numpy.sum(numpy.exp(logqeval)*vTglw*vzglw*sigmaz1)*1.5

    def pvT(self,vT,R,z,gl=True,ngl=_DEFAULTNGL2):
        """
        NAME:
           pvT
        PURPOSE:
           calculate the marginalized vT probability at this location (NOT normalized by the density)
        INPUT:

           vT - tangential velocity (/vo)

           R - radius (/ro)

           z - height (/ro)

           gl - use Gauss-Legendre integration (True, currently the only option)

           ngl - order of Gauss-Legendre integration

        OUTPUT:
           p(vT,R,z)
        HISTORY:
           2012-12-22 - Written - Bovy (IAS)
        """
        sigmaR1= self._sr*numpy.exp((self._ro-R)/self._hsr)
        sigmaz1= self._sz*numpy.exp((self._ro-R)/self._hsz)
        if gl:
            if ngl % 2 == 1:
                raise ValueError("ngl must be even")
            #Use Gauss-Legendre integration for all
            if ngl == _DEFAULTNGL:
                glx, glw= self._glxdef, self._glwdef
                glx12, glw12= self._glxdef12, self._glwdef12
            elif ngl == _DEFAULTNGL2:
                glx, glw= self._glxdef2, self._glwdef2
                glx12, glw12= self._glxdef, self._glwdef
            else:
                glx, glw= numpy.polynomial.legendre.leggauss(ngl)
                glx12, glw12= numpy.polynomial.legendre.leggauss(ngl/2)
            #Evaluate everywhere
            if isinstance(self._aA,(actionAngle.actionAngleAdiabatic,
                                    actionAngle.actionAngleAdiabaticGrid)):
                vRgl= 4.*sigmaR1/2.*(glx+1.)
                vzgl= 4.*sigmaz1/2.*(glx+1.)
                vRglw= glw
                vzglw= glw
            else:
                vRgl= 4.*sigmaR1/2.*(glx12+1.)
                vRgl= list(vRgl)
                vRgl.extend(-4.*sigmaR1/2.*(glx12+1.))
                vRgl= numpy.array(vRgl)
                vzgl= 4.*sigmaz1/2.*(glx12+1.)
                vzgl= list(vzgl)
                vzgl.extend(-4.*sigmaz1/2.*(glx12+1.))
                vzgl= numpy.array(vzgl)
                vRglw= glw12
                vRglw= list(vRglw)
                vRglw.extend(glw12)
                vRglw= numpy.array(vRglw)
                vzglw= glw12
                vzglw= list(vzglw)
                vzglw.extend(glw12)
                vzglw= numpy.array(vzglw)
            #Tile everything
            vRgl= numpy.tile(vRgl,(ngl,1)).T
            vzgl= numpy.tile(vzgl,(ngl,1))
            vRglw= numpy.tile(vRglw,(ngl,1)).T #also tile weights
            vzglw= numpy.tile(vzglw,(ngl,1))
            #evaluate
            logqeval= numpy.reshape(self(R+numpy.zeros(ngl*ngl),
                                         vRgl.flatten(),
                                         vT+numpy.zeros(ngl*ngl),
                                         z+numpy.zeros(ngl*ngl),
                                         vzgl.flatten(),
                                         log=True),
                                    (ngl,ngl))
            return numpy.sum(numpy.exp(logqeval)*vRglw*vzglw*sigmaR1*sigmaz1)

    def pvz(self,vz,R,z,gl=True,ngl=_DEFAULTNGL2,
            _return_actions=False,_jr=None,_lz=None,_jz=None,
            _return_freqs=False,
            _rg=None,_kappa=None,_nu=None,_Omega=None,
            _sigmaR1=None):
        """
        NAME:
           pvz
        PURPOSE:
           calculate the marginalized vz probability at this location (NOT normalized by the density)
        INPUT:
           vz - vertical velocity (/vo)

           R - radius (/ro)

           z - height (/ro)

           gl - use Gauss-Legendre integration (True, currently the only option)

           ngl - order of Gauss-Legendre integration

        OUTPUT:
           p(vz,R,z)
        HISTORY:
           2012-12-22 - Written - Bovy (IAS)
        """
        if _sigmaR1 is None:
            sigmaR1= self._sr*numpy.exp((self._ro-R)/self._hsr)
        else:
            sigmaR1= _sigmaR1
        if gl:
            if ngl % 2 == 1:
                raise ValueError("ngl must be even")
            #Use Gauss-Legendre integration for all
            if ngl == _DEFAULTNGL:
                glx, glw= self._glxdef, self._glwdef
                glx12, glw12= self._glxdef12, self._glwdef12
            elif ngl == _DEFAULTNGL2:
                glx, glw= self._glxdef2, self._glwdef2
                glx12, glw12= self._glxdef, self._glwdef
            else:
                glx, glw= numpy.polynomial.legendre.leggauss(ngl)
                glx12, glw12= numpy.polynomial.legendre.leggauss(ngl/2)
            #Evaluate everywhere
            if isinstance(self._aA,(actionAngle.actionAngleAdiabatic,
                                    actionAngle.actionAngleAdiabaticGrid)):
                vRgl= (glx+1.)
                vRglw= glw
            else:
                vRgl= (glx12+1.)
                vRgl= list(vRgl)
                vRgl.extend(-(glx12+1.))
                vRgl= numpy.array(vRgl)
                vRglw= glw12
                vRglw= list(vRglw)
                vRglw.extend(glw12)
                vRglw= numpy.array(vRglw)
            vTgl= 1.5/2.*(glx+1.)
            #Tile everything
            vTgl= numpy.tile(vTgl,(ngl,1)).T
            vRgl= numpy.tile(vRgl,(ngl,1))
            vTglw= numpy.tile(glw,(ngl,1)).T #also tile weights
            vRglw= numpy.tile(vRglw,(ngl,1))
            #If inputs are arrays, tile
            if isinstance(R,numpy.ndarray):
                nR= len(R)
                R= numpy.tile(R,(ngl,ngl,1)).T.flatten()
                z= numpy.tile(z,(ngl,ngl,1)).T.flatten()
                vz= numpy.tile(vz,(ngl,ngl,1)).T.flatten()
                vTgl= numpy.tile(vTgl,(nR,1,1)).flatten()
                vRgl= numpy.tile(vRgl,(nR,1,1)).flatten()
                vTglw= numpy.tile(vTglw,(nR,1,1))
                vRglw= numpy.tile(vRglw,(nR,1,1))
                scalarOut= False
            else:
                R= R+numpy.zeros(ngl*ngl)
                z= z+numpy.zeros(ngl*ngl)
                vz= vz+numpy.zeros(ngl*ngl)
                nR= 1
                scalarOut= True
                vRgl= vRgl.flatten()
            vRgl*= numpy.tile(4.*sigmaR1/2.,(ngl,ngl,1)).T.flatten()
            #evaluate
            if _jr is None and _rg is None:
                logqeval, jr, lz, jz, rg, kappa, nu, Omega= self(R,
                                                                 vRgl.flatten(),
                                                                 vTgl.flatten(),
                                                                 z,
                                                                 vz,
                                                                 log=True,
                                                                 _return_actions=True,
                                                                 _return_freqs=True)
                logqeval= numpy.reshape(logqeval,(nR,ngl*ngl))
            elif not _jr is None and not _rg is None:
                logqeval, jr, lz, jz, rg, kappa, nu, Omega= self((_jr,_lz,_jz),
                                           rg=_rg,kappa=_kappa,nu=_nu,
                                           Omega=_Omega,
                                           log=True,
                                           _return_actions=True,
                                           _return_freqs=True)
                logqeval= numpy.reshape(logqeval,(nR,ngl*ngl))
            elif not _jr is None and _rg is None:
                logqeval, jr, lz, jz, rg, kappa, nu, Omega= self((_jr,_lz,_jz),
                                           log=True,
                                           _return_actions=True,
                                           _return_freqs=True)
                logqeval= numpy.reshape(logqeval,(nR,ngl*ngl))
            elif _jr is None and not _rg is None:
                logqeval, jr, lz, jz, rg, kappa, nu, Omega= self(R,
                                             vRgl.flatten(),
                                             vTgl.flatten(),
                                             z,
                                             vz,
                                             rg=_rg,kappa=_kappa,nu=_nu,
                                             Omega=_Omega,
                                             log=True,
                                             _return_actions=True,
                                             _return_freqs=True)
                logqeval= numpy.reshape(logqeval,(nR,ngl*ngl))
            vRglw= numpy.reshape(vRglw,(nR,ngl*ngl))
            vTglw= numpy.reshape(vTglw,(nR,ngl*ngl))
            if scalarOut:
                result= numpy.sum(numpy.exp(logqeval)*vTglw*vRglw,axis=1)[0]*sigmaR1*1.5
            else:
                result= numpy.sum(numpy.exp(logqeval)*vTglw*vRglw,axis=1)*sigmaR1*1.5
            if _return_actions and _return_freqs:
                return (result,
                        jr,lz,jz,
                        rg, kappa, nu, Omega)          
            elif _return_freqs:
                return (result,
                        rg, kappa, nu, Omega)
            elif _return_actions:
                return (result,
                        jr,lz,jz)
            else:
                return result

    def pvRvT(self,vR,vT,R,z,gl=True,ngl=_DEFAULTNGL2):
        """
        NAME:
           pvRvT
        PURPOSE:
           calculate the marginalized (vR,vT) probability at this location (NOT normalized by the density)
        INPUT:

           vR - radial velocity (/vo)

           vT - tangential velocity (/vo)

           R - radius (/ro)

           z - height (/ro)

           gl - use Gauss-Legendre integration (True, currently the only option)

           ngl - order of Gauss-Legendre integration

        OUTPUT:
           p(vR,vT,R,z)
        HISTORY:
           2013-01-02 - Written - Bovy (IAS)
        """
        sigmaz1= self._sz*numpy.exp((self._ro-R)/self._hsz)
        if gl:
            if ngl % 2 == 1:
                raise ValueError("ngl must be even")
            #Use Gauss-Legendre integration for all
            if ngl == _DEFAULTNGL:
                glx, glw= self._glxdef, self._glwdef
                glx12, glw12= self._glxdef12, self._glwdef12
            elif ngl == _DEFAULTNGL2:
                glx, glw= self._glxdef2, self._glwdef2
                glx12, glw12= self._glxdef, self._glwdef
            else:
                glx, glw= numpy.polynomial.legendre.leggauss(ngl)
                glx12, glw12= numpy.polynomial.legendre.leggauss(ngl/2)
            #Evaluate everywhere
            if isinstance(self._aA,(actionAngle.actionAngleAdiabatic,
                                    actionAngle.actionAngleAdiabaticGrid)):
                vzgl= 4.*sigmaz1/2.*(glx+1.)
                vzglw= glw
            else:
                vzgl= 4.*sigmaz1/2.*(glx12+1.)
                vzgl= list(vzgl)
                vzgl.extend(-4.*sigmaz1/2.*(glx12+1.))
                vzgl= numpy.array(vzgl)
                vzglw= glw12
                vzglw= list(vzglw)
                vzglw.extend(glw12)
                vzglw= numpy.array(vzglw)
            #evaluate
            logqeval= self(R+numpy.zeros(ngl),
                           vR+numpy.zeros(ngl),
                           vT+numpy.zeros(ngl),
                           z+numpy.zeros(ngl),
                           vzgl,
                           log=True)
            return numpy.sum(numpy.exp(logqeval)*vzglw*sigmaz1)
        
    def pvTvz(self,vT,vz,R,z,gl=True,ngl=_DEFAULTNGL2):
        """
        NAME:
           pvTvz
        PURPOSE:
           calculate the marginalized (vT,vz) probability at this location (NOT normalized by the density)
        INPUT:

           vT - tangential velocity (/vo)

           vz - vertical velocity (/vo)

           R - radius (/ro)

           z - height (/ro)

           gl - use Gauss-Legendre integration (True, currently the only option)

           ngl - order of Gauss-Legendre integration

        OUTPUT:
           p(vT,vz,R,z)
        HISTORY:
           2012-12-22 - Written - Bovy (IAS)
        """
        sigmaR1= self._sr*numpy.exp((self._ro-R)/self._hsr)
        if gl:
            if ngl % 2 == 1:
                raise ValueError("ngl must be even")
            #Use Gauss-Legendre integration for all
            if ngl == _DEFAULTNGL:
                glx, glw= self._glxdef, self._glwdef
                glx12, glw12= self._glxdef12, self._glwdef12
            elif ngl == _DEFAULTNGL2:
                glx, glw= self._glxdef2, self._glwdef2
                glx12, glw12= self._glxdef, self._glwdef
            else:
                glx, glw= numpy.polynomial.legendre.leggauss(ngl)
                glx12, glw12= numpy.polynomial.legendre.leggauss(ngl/2)
            #Evaluate everywhere
            if isinstance(self._aA,(actionAngle.actionAngleAdiabatic,
                                    actionAngle.actionAngleAdiabaticGrid)):
                vRgl= 4.*sigmaR1/2.*(glx+1.)
                vRglw= glw
            else:
                vRgl= 4.*sigmaR1/2.*(glx12+1.)
                vRgl= list(vRgl)
                vRgl.extend(-4.*sigmaR1/2.*(glx12+1.))
                vRgl= numpy.array(vRgl)
                vRglw= glw12
                vRglw= list(vRglw)
                vRglw.extend(glw12)
                vRglw= numpy.array(vRglw)
            #evaluate
            logqeval= self(R+numpy.zeros(ngl),
                           vRgl,
                           vT+numpy.zeros(ngl),
                           z+numpy.zeros(ngl),
                           vz+numpy.zeros(ngl),
                           log=True)
            return numpy.sum(numpy.exp(logqeval)*vRglw*sigmaR1)

    def pvRvz(self,vR,vz,R,z,gl=True,ngl=_DEFAULTNGL2):
        """
        NAME:
           pvR
        PURPOSE:
           calculate the marginalized (vR,vz) probability at this location (NOT normalized by the density)
        INPUT:

           vR - radial velocity (/vo)

           vz - vertical velocity (/vo)

           R - radius (/ro)

           z - height (/ro)

           gl - use Gauss-Legendre integration (True, currently the only option)

           ngl - order of Gauss-Legendre integration

        OUTPUT:
           p(vR,vz,R,z)
        HISTORY:
           2013-01-02 - Written - Bovy (IAS)
        """
        if gl:
            if ngl % 2 == 1:
                raise ValueError("ngl must be even")
            #Use Gauss-Legendre integration for all
            if ngl == _DEFAULTNGL:
                glx, glw= self._glxdef, self._glwdef
                glx12, glw12= self._glxdef12, self._glwdef12
            elif ngl == _DEFAULTNGL2:
                glx, glw= self._glxdef2, self._glwdef2
                glx12, glw12= self._glxdef, self._glwdef
            else:
                glx, glw= numpy.polynomial.legendre.leggauss(ngl)
                glx12, glw12= numpy.polynomial.legendre.leggauss(ngl/2)
            #Evaluate everywhere
            vTgl= 1.5/2.*(glx+1.)
            vTglw= glw
            #If inputs are arrays, tile
            if isinstance(R,numpy.ndarray):
                nR= len(R)
                R= numpy.tile(R,(ngl,1)).T.flatten()
                z= numpy.tile(z,(ngl,1)).T.flatten()
                vR= numpy.tile(vR,(ngl,1)).T.flatten()
                vz= numpy.tile(vz,(ngl,1)).T.flatten()
                vTgl= numpy.tile(vTgl,(nR,1)).flatten()
                vTglw= numpy.tile(vTglw,(nR,1))
                scalarOut= False
            else:
                R= R+numpy.zeros(ngl)
                vR= vR+numpy.zeros(ngl)
                z= z+numpy.zeros(ngl)
                vz= vz+numpy.zeros(ngl)
                nR= 1
                scalarOut= True
            #evaluate
            logqeval= numpy.reshape(self(R,
                                         vR,
                                         vTgl,
                                         z,
                                         vz,
                                         log=True),
                                    (nR,ngl))
            out= numpy.sum(numpy.exp(logqeval)*vTglw,axis=1)
            if scalarOut: return out[0]
            else: return out

    def _calc_epifreq(self,r):
        """
        NAME:
           _calc_epifreq
        PURPOSE:
           calculate the epicycle frequency at r
        INPUT:
           r - radius
        OUTPUT:
           kappa
        HISTORY:
           2012-07-25 - Written - Bovy (IAS@MPIA)
        NOTE:
           takes about 0.1 ms for a Miyamoto-Nagai potential
        """
        return potential.epifreq(self._pot,r)

    def _calc_verticalfreq(self,r):
        """
        NAME:
           _calc_verticalfreq
        PURPOSE:
           calculate the vertical frequency at r
        INPUT:
           r - radius
        OUTPUT:
           nu
        HISTORY:
           2012-07-25 - Written - Bovy (IAS@MPIA)
        NOTE:
           takes about 0.05 ms for a Miyamoto-Nagai potential
        """
        return potential.verticalfreq(self._pot,r)

    def rg(self,lz):
        """
        NAME:
           rg
        PURPOSE:
           calculate the radius of a circular orbit of Lz
        INPUT:
           lz - Angular momentum
        OUTPUT:
           radius
        HISTORY:
           2012-07-25 - Written - Bovy (IAS@MPIA)
        NOTE:
           seems to take about ~0.5 ms for a Miyamoto-Nagai potential; 
           ~0.75 ms for a MWPotential
           about the same with or without interpolation of the rotation curve

           Not sure what to do about negative lz...
        """
        if isinstance(lz,numpy.ndarray):
            indx= (lz > self._precomputergLzmax)*(lz < self._precomputergLzmin)
            indxc= True-indx
            out= numpy.empty(lz.shape)
            out[indxc]= self._rgInterp(lz[indxc])
            out[indx]= numpy.array([potential.rl(self._pot,lz[indx][ii]) for ii in range(numpy.sum(indx))])
            return out
        else:
            if lz > self._precomputergLzmax or lz < self._precomputergLzmin:
                return potential.rl(self._pot,lz)
            return numpy.atleast_1d(self._rgInterp(lz))

def _vmomentsurfaceIntegrand(vz,vR,vT,R,z,df,sigmaR1,gamma,sigmaz1,n,m,o): #pragma: no cover because this is too slow; a warning is shown
    """Internal function that is the integrand for the vmomentsurface mass integration"""
    return vR**n*vT**m*vz**o*df(R,vR*sigmaR1,vT*sigmaR1*gamma,z,vz*sigmaz1)

def _vmomentsurfaceMCIntegrand(vz,vR,vT,R,z,df,sigmaR1,gamma,sigmaz1,mvT,n,m,o):
    """Internal function that is the integrand for the vmomentsurface mass integration"""
    return vR**n*vT**m*vz**o*df(R,vR*sigmaR1,vT*sigmaR1*gamma,z,vz*sigmaz1)*numpy.exp(vR**2./2.+(vT-mvT)**2./2.+vz**2./2.)

def _jmomentsurfaceIntegrand(vz,vR,vT,R,z,df,sigmaR1,gamma,sigmaz1,n,m,o): #pragma: no cover because this is too slow; a warning is shown
    """Internal function that is the integrand for the vmomentsurface mass integration"""
    return df(R,vR*sigmaR1,vT*sigmaR1*gamma,z,vz*sigmaz1,
              func= (lambda x,y,z: x**n*y**m*z**o))

def _jmomentsurfaceMCIntegrand(vz,vR,vT,R,z,df,sigmaR1,gamma,sigmaz1,mvT,n,m,o):
    """Internal function that is the integrand for the vmomentsurface mass integration"""
    return df(R,vR*sigmaR1,vT*sigmaR1*gamma,z,vz*sigmaz1,
              func=(lambda x,y,z: x**n*y**m*z**o))\
              *numpy.exp(vR**2./2.+(vT-mvT)**2./2.+vz**2./2.)

