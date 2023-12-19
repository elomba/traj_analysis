from __future__ import print_function
import math

import numpy as np

class LAMMPSPairPotential(object):
    def __init__(self):
        self.pmap=dict()
        self.units='lj'
    def map_coeff(self,name,ltype):
        self.pmap[ltype]=name
    def check_units(self,units):
        if (units != self.units):
           raise Exception("Conflicting units: %s vs. %s" % (self.units,units))

class Harmonic(LAMMPSPairPotential):
    def __init__(self):
        super(Harmonic,self).__init__()
        self.units = 'real'
        # set coeffs: K, r0
        self.coeff = {'A'  : {'A'  : (0.2,9.0),
                              'B'  : (math.sqrt(0.2*0.4),9.0)},
                      'B'  : {'A'  : (math.sqrt(0.2*0.4),9.0),
                              'B'  : (0.4,9.0)}}

    def compute_force(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r = math.sqrt(rsq)
        delta = coeff[1]-r
        if (r <= coeff[1]):
          return 2.0*delta*coeff[0]/r
        else:
          return 0.0

    def compute_energy(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r = math.sqrt(rsq)
        delta = coeff[1]-r
        if (r <= coeff[1]):
          return delta*delta*coeff[0]
        else:
          return 0.0

class LJCutReal(LAMMPSPairPotential):
    def __init__(self):
        super(LJCutReal,self).__init__()
        # set coeffs: 48*eps*sig**12, 24*eps*sig**6,
        #              4*eps*sig**12,  4*eps*sig**6
        self.units = 'real'
        eps = 0.838
        sig = 4.0
        self.coeff = {'1'  : {'1'  : (48.0*eps*sig**12,24.0*eps*sig**6,
                                         4.0*eps*sig**12, 4.0*eps*sig**6)}, '2'  : {'1'  : (48.0*eps*sig**12,24.0*eps*sig**6,
                                         4.0*eps*sig**12, 4.0*eps*sig**6)}}

    def compute_force(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj1 = coeff[0]
        lj2 = coeff[1]
        return (r6inv * (lj1*r6inv - lj2))*r2inv

    def compute_energy(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj3 = coeff[2]
        lj4 = coeff[3]
        return (r6inv * (lj3*r6inv - lj4))


class LJCutSPCE(LAMMPSPairPotential):
    def __init__(self):
        super(LJCutSPCE,self).__init__()
        self.units='real'
        # SPCE oxygen LJ parameters in real units
        eps=0.15535
        sig=3.166
        self.coeff = {'OW'  : {'OW'  : (48.0*eps*sig**12,24.0*eps*sig**6,
                                         4.0*eps*sig**12, 4.0*eps*sig**6),
                               'HW'  : (0.0,0.0, 0.0,0.0)},
                      'HW'  : {'OW'  : (0.0,0.0, 0.0,0.0),
                               'HW'  : (0.0,0.0, 0.0,0.0)}}

    def compute_force(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj1 = coeff[0]
        lj2 = coeff[1]
        return (r6inv * (lj1*r6inv - lj2))*r2inv

    def compute_energy(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj3 = coeff[2]
        lj4 = coeff[3]
        return (r6inv * (lj3*r6inv - lj4))

class SALR(LAMMPSPairPotential):
    def __init__(self):
        super(SALR,self).__init__()
        # set coeffs: 48*eps*sig**12, 24*eps*sig**6,
        #              4*eps*sig**12,  4*eps*sig**6
        self.units = 'real'
        eps = 0.238
        sig = 4.0
        sig2 = sig*sig
        sig10=sig2*sig2*sig2*sig2*sig2
        k1=1.0
        a1=0.1/sig
        k2=2.0
        a2=0.25/sig
        rc = 100.0
        r2inv  = 1/rc**2
        r10inv  = r2inv*r2inv*r2inv*r2inv*r2inv
        urc = sig10*r10inv+k1*np.exp(-a1*rc)-k2*np.exp(-a2*rc)
        self.coeff = {'1'  : {'1'  : (eps,sig10,k1, k2, a1, a2, urc),
                              '2'  : (eps,sig10,k1, k2, a1, a2, urc)},
                      '2'  : {'1'  : (eps,sig10,k1, k2, a1, a2, urc),
                              '2'  : (eps,sig10,k1, k2, a1, a2, urc)}}

    def compute_force(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r=np.sqrt(rsq)
        r10inv  = r2inv*r2inv*r2inv*r2inv*r2inv
        eps = coeff[0]
        sig10 = coeff[1]
        k1 = coeff[2]
        k2 = coeff[3]
        a1 = coeff[4]
        a2 = coeff[5]
        duslr = 10*sig10*r10inv*r2inv+(k1*a1*np.exp(-a1*r)-k2*a2*np.exp(-a2*r))/r
        return (eps*duslr)
        # r6inv  = r2inv*r2inv*r2inv

        # return (r6inv * (lj1*r6inv - lj2))*r2inv

    def compute_energy(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r=np.sqrt(rsq)
        r10inv  = r2inv*r2inv*r2inv*r2inv*r2inv
        eps = coeff[0]
        sig10 = coeff[1]
        k1 = coeff[2]
        k2 = coeff[3]
        a1 = coeff[4]
        a2 = coeff[5]
        urc = coeff[6]
        uslr = sig10*r10inv+k1*np.exp(-a1*r)-k2*np.exp(-a2*r)-urc
        return (eps*uslr)
