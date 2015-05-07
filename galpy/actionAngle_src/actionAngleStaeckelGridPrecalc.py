#_____import packages_____
import time
from galpy.actionAngle import actionAngleStaeckel,actionAngleStaeckelGrid
from galpy.actionAngle_src.actionAngleStaeckel_c import _ext_loaded as ext_loaded

class actionAngleStaeckelGridPrecalc(actionAngleStaeckelGrid):
    """Action-angle formalism for axisymmetric potentials using 
       Binney (2012)'s Staeckel approximation, grid-based interpolation. 
       This class, derived from actionAngleStaeckelGrid, can initialize 
       a actionAngleStaeckelGrid Object using pre-calculated grids."""
    def __init__(self,pot=None,delta=None,
                      Lzs=None,ERL=None,ERa=None,
                      mu0=None,mjr=None,mjz=None,
                      numcores=1,**kwargs):

        #____keywords needed for super class_____
        if pot is None:
            raise IOError("Must specify pot= for actionAngleStaeckelGridPrecalc")
        self._pot= pot
        if delta is None:
            raise IOError("Must specify delta= for actionAngleStaeckelGridPrecalc")
        self._delta = delta
        self._nLz   = len(Lzs)                     #precalculated!!!
        self._nE    = len(mu0)/self._nLz           #precalculated!!!
        self._npsi  = len(mjr)/self._nLz/self._nE  #precalculated!!!

        #____c???_____
        if ext_loaded and 'c' in kwargs and kwargs['c']:
            self._c= True
        else:
            self._c= False
        

        #____initialize super-class_____
        start = time.time()
        actionAngleStaeckelGrid.__init__(self,
                    pot=pot,delta=self._delta,
                    Rmax=None,  #not needed
                    nE=self._nE,npsi=self._npsi,nLz=self._nLz,numcores=numcores,**kwargs)
        print time.time()-start

        #_____Grid____
        self._Lzmin= 0.01
        self._Lzs = Lzs #precalculated!!!
        self._Lzmax= self._Lzs[-1]

        #_____E_c(R=RL), energy of circular orbit_____
        self._ERL       = ERL #precalculated!!!
        self._ERLmax    = numpy.amax(self._ERL)+1.
        self._ERLInterp= interpolate.InterpolatedUnivariateSpline(self._Lzs,
                                                                  numpy.log(-(self._ERL-self._ERLmax)),k=3)

        self._ERa       = ERa #precalculated!!!
        self._ERamax    = numpy.amax(self._ERa)+1.
        self._ERaInterp = interpolate.InterpolatedUnivariateSpline(self._Lzs,
                                                                  numpy.log(-(self._ERa-self._ERamax)),k=3)
        y = numpy.linspace(0.,1.,self._nE)


        #_____setup u0_____
        u0 = numpy.reshape(mu0,(self._nLz,self._nE)) #mu0 precalcualted!!!

        #_____Setup and interpolate actions_____
        jr= numpy.reshape(mjr,(self._nLz,self._nE,self._npsi))  #mjr precalcualted!!!
        jz= numpy.reshape(mjz,(self._nLz,self._nE,self._npsi))  #mjz precalcualted!!!
        jrLzE= numpy.zeros((self._nLz))
        jzLzE= numpy.zeros((self._nLz))
        for ii in range(self._nLz):
            jrLzE[ii]= numpy.nanmax(jr[ii,(jr[ii,:,:] != 9999.99)])#:,:])
            jzLzE[ii]= numpy.nanmax(jz[ii,(jz[ii,:,:] != 9999.99)])#:,:])
        jrLzE[(jrLzE == 0.)]= numpy.nanmin(jrLzE[(jrLzE > 0.)])
        jzLzE[(jzLzE == 0.)]= numpy.nanmin(jzLzE[(jzLzE > 0.)])
        for ii in range(self._nLz):
            jr[ii,:,:]/= jrLzE[ii]
            jz[ii,:,:]/= jzLzE[ii]
        #Deal w/ 9999.99
        jr[(jr > 1.)]= 1.
        jz[(jz > 1.)]= 1.
        #Deal w/ NaN
        jr[numpy.isnan(jr)]= 0.
        jz[numpy.isnan(jz)]= 0.
        #First interpolate the maxima
        self._jr= jr
        self._jz= jz
        self._jrLzInterp= interpolate.InterpolatedUnivariateSpline(self._Lzs,
                                                                   numpy.log(jrLzE+10.**-5.),k=3)
        self._jzLzInterp= interpolate.InterpolatedUnivariateSpline(self._Lzs,
                                                                   numpy.log(jzLzE+10.**-5.),k=3)

        #_____Interpolate u0_____
        self._logu0Interp= interpolate.RectBivariateSpline(self._Lzs,
                                                           y,
                                                           numpy.log(u0),
                                                           kx=3,ky=3,s=0.)

        #_____spline filter jr and jz_____
        #such that they can be used with ndimage.map_coordinates
        self._jrFiltered= ndimage.spline_filter(numpy.log(self._jr+10.**-10.),order=3)
        self._jzFiltered= ndimage.spline_filter(numpy.log(self._jz+10.**-10.),order=3)

        #_____Set up the actionAngleStaeckel object for special cases______
        self._aA= actionAngleStaeckel.actionAngleStaeckel(pot=self._pot,delta=self._delta,c=self._c)


