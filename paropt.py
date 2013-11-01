import math
import numpy, sympy
from numpy import *
from functools import reduce
from pprint import pprint
import itertools
from sympy import symbols, diff, Symbol

#DEBUG = 20
DEBUG = 0
PROFILE = True
Ht = 0.01
import time

def constant(function):
    if CACHING:
        def wrapper(self, *args, **kwargs):
            c=self.__cache__
            if function in c:
                return c[function]
            else:
                rv = function(self, *args, **kwargs)
                c[function] = rv
                return rv
        return wrapper
    else:
        return function

def compile_method(name, f):
    s="""
def %s(self, t, x, u):
    return %s
    """ % (name, f)

    f=compile(s, "<diffunctions>", 'exec')
    f=exec(f)

def _vcomp(exp):
    sexp="array(%s)" % str(exp)
    return compile(sexp, "<-%s->" % sexp, "eval")

def _comp(exp):
    sexp=str(exp)
    return compile(sexp, "<-%s->" % sexp, "eval")

def _eval(f, t, xc, uc, V):
    g={'t':t, 'array':array}
    for i, x in enumerate(V.x[0]):
        g.setdefault(str(x),xc[0][i])
    for i, u in enumerate(V.u[0]):
        g.setdefault(str(u),uc[0][i])

    rc = eval(f, g)
    return rc

def _vdiff(F, V, num=1):
    answer=[]
    DF=[]

    for f in F[0]:
        dfs=[]
        for v in V[0]:
            df=diff(f, v, num)
            dfs.append(df)
        DF.append(dfs)
    return DF

def _diff(f, V, num=1):
    answer=[]
    DF=[]

    dfs=[]
    for v in V[0]:
        df=diff(f, v, num)
        dfs.append(df)
    DF.append(dfs)
    return DF

def _teval(d, T,X,U, V): # d is a code of function to be evaluated at all T time instance.
    rc=[]
    for t in T:
        rc.append(_eval(d, t, X[t], U[t], V))
    return rc


class Helper():
    pass

class ParOptModel(object):
    """Parametric Optimised model class.
    See docs for individual methods on usage ways.
    """

    def __init__(self, X0, N, M):
        self.N=N     # Dimention of x
        self.M=M     # Dimention of u
        self.X0=X0
        self.v=Helper()
        self.v.t=Symbol('t')
        a=[]
        for i in range(N):
            a.append([Symbol('x%s' % i)])
        self.v.x=a
        a=[]
        for i in range(M):
            a.append([Symbol('u%s' % i)])
        self.v.u=a

        t, x, u = self.v.t, self.v.x, self.v.u

        self._f=self.f(t, x, u)
        self._f0=self.f0(t, x, u)
        self._F=self.F(x)

        # -----------------
        self._dfdx=_vcomp(_vdiff(self._f, x))
        self._df0dx=_vcomp(_diff(self._f0, x))
        self._dFdx=_vcomp(_diff(self._F, x))

        self._dfdu=_vcomp(_vdiff(self._f, u))
        self._df0du=_vcomp(_diff(self._f0, u))

        if DEBUG>=5:
            print (self._f, self._f0, self._dfdx, self._df0dx, sep=", ")

        #self.dfdu=compile_method("dfdu", self._dfdu)
        if DEBUG>=5:
            print (self._F, self._dFdx, sep=", ")
            print (self._f, self._f0, self._dfdu, self._df0du, sep=", ")
        self.__cache__={}

    def I(self, X, U):
        def _add(acc, t):
            return acc + self.f0(t, X[t], U[t])
        return reduce(_add, range(len(X)-1), self.F(X[-1]))

    def F(self, x):
        """x - ending state of a trajectory
        """
        raise RuntimeError("should be implemented by subclass")

    def f(self, t, x0, U, dt=1):
        raise RuntimeError("should be implemented by subclass")

    def f0(self, t, x, u, dt=1):
        raise RuntimeError("should be implemented by subclass")

    def dFdx(self, xe):
        return eval(self._dFdx, {'x':xe, 'array':array})

    def dfdx(self, T, X, U, dt=1):
        return _teval(self._dfdx,  T,X,U, self.v)

    def df0dx(self, T, X, U, dt=1):
        return _teval(self._df0dx, T,X,U, self.v)

    def dfdu(self, T, X, U, dt=1):
        return _teval(self._dfdu, T,X,U, self.v)

    def df0du(self, T, X, U, dt=1):
        return _teval(self._df0du, T,X,U, self.v)

    def Psi(self, t, X, U, alpha):

        psie = -self.dFdx(X[-1])

        psi=[psie]

        _df0dx=self.df0dx(t[:-1], X[:-1], U) # last element is useless
        _dfdx =self.dfdx (t[:-1], X[:-1], U) # last element is useless


        tt=t[:-1]

        # -----
        # rt.reverse()
        j=len(tt)-1
        p=psie

        while j>=1:
            i=tt[j]
            pp=p
            _dfdx_t = transpose(_dfdx[i])
            _df0dx_t = transpose(_df0dx[i])
            pn = dot(_dfdx_t, pp) - alpha*_df0dx_t
            psi.append(pn)
            p=pn
            j-=1

        # ----

        psi=array(psi)
        #lp=len(psi)
        #for i in range(len(psi) / 2 ):
        #    p=psi[i]
        #    psi[i]=psi[lp-i-1]
        #    psi[lp-i-1]=p
        # so Psi[0] is Psi[t+1] at the start t
        return psi[::-1]

    def dHdu(self, t, X, U, Psi):
        _df0du=self.df0du(t[:-1], X[:-1], U)
        _dfdu =self.dfdu (t[:-1], X[:-1], U)
        p=Psi
        if DEBUG>10:
            print ("dfdu:", len(_dfdu))
            print ("psi:", len(p))
            print ("df0du", len(_df0du))
        _s=[dot(transpose(Psi[i]),_dfdu[i]) for i in range(len(X)-1)] # Psi[i] is shifted left to 1 step.
        return array(_s) + _df0du # FIXME Check dot operation.

    def start_control(self):
        raise RuntimeError("should be implemented by subclass")

    #------------ These functions must be defined for Second Order Improvement Process -----

    def dPsidalpha(self, t, X, U):
        psie = -self.dFdx(X[-1])
        psi=[psie]
        _df0dx=self.df0dx(t, X, U) # last element is useless
        _dfdx =self.dfdx (t, X, U) # last element is useless

        tt=t[:-1]

        # -----
        # rt.reverse()
        j=len(tt)-1
        p=psie

        while j>=0:
            i=tt[j]
            pp=p
            _dfdx_t = transpose(_dfdx[i])
            _df0dx_t = transpose(_df0dx[i])
            pn = dot(_dfdx_t, pp) - _df0dx_t
            psi.append(pn)
            p=pn
            j-=1
        # ----


        psi=array(psi)

        return psi[::-1] # Really it is dPsidalpha

    def dSigmadalpha(self, t, X, U):
        sige  = -self.dFdxdx(X[-1])
        sig=[sige]

        _df0dx=self.df0dx(t, X, U) # last element is useless
        _dfdx =self.dfdx (t, X, U) # last element is useless

        tt=t[:-1]

        # -----
        # rt.reverse()
        j=len(tt)-1
        s=sige

        while j>=0:
            i=tt[j]
            sp=s
            _dfdx_t = transpose(_dfdx[i])
            _df0dx_t = transpose(_df0dx[i])
            sn = dot(_dfdx_t, sp) - _df0dx_t
            sig.append(sn)
            s=sn
            j-=1
        # ----


        sig=array(sig)

        return sig[::-1] # Really it is dSigmadalpha




class ParOptProcess(object):
    """This calss corresponds to a optimizational model.
    """

    def __init__(self, model, alpha=1.0, beta=1.0):
        """
        """
        self.model=model
        self.alpha=alpha
        self.beta=beta

    def optimize(self, t, eps=0.001, iters=1000):
        Up=self.model.start_control()
        Xp=self.trajectory(Up)
        Ip=self.model.I(Xp,Up)

        it = 1
        while True:
            beta = self.beta
            while True:
                (Xn, Un) = self.improve(t, Xp, Up, beta=beta)
                In = self.model.I(Xn, Un)

                dI = Ip-In
                if DEBUG>5:
                    print ("DI:", dI)
                if DEBUG>6:
                    print ("Xn:", Xn)
                    print ("Un:", Un)
                if abs(dI)<eps:
                    return In, Xn, Un, it, "Opt"
                if iters<=0:
                    return In, Xn, Un, it, "Nonoptimal"
                iters-=1
                it+=1

                if In>=Ip:
                    beta/=2
                    continue
                else:
                    Xp, Up, Ip = Xn, Un, In
                    break

        raise RuntimeError("this should be not reached")

    def trajectory(self, U):
        x0=self.model.X0
        X = [x0]

        for t, u in enumerate(U): # (0, u0), (1, u1)
            xn=self.model.f(t, X[t], u)
            X.append(xn)

        return array(X)

    def improve(self, t, X, U, **kwargs): # (a,b,c)
        Psi=self.model.Psi(t, X, U, self.alpha)
        _dU=self.dU(t, X, U, Psi=Psi, beta=kwargs['beta'])
        Un = U + _dU
        return self.trajectory(Un), Un

    def dU(self, t, X, U, **kwargs):
        _dHdu=self.model.dHdu(t, X, U, Psi=kwargs['Psi'])
        return _dHdu * kwargs['beta']

class SeconOrderParOptProcess(ParOptProcess):
    """A second order parametric optimization process
    """
    def __init__(self):
        X0=array([0.0, 0.0])
        self.t = arange(2)
        ParOptModel.__init__(self, X0=X0)


    def start_control(self):
        return [
            array([0.1, 0.0]),
            array([0.2, 1.0])
        ]

    def optimize(self, t, eps=0.001, iters=1000):
        Up=self.model.start_control()
        Xp=self.trajectory(Up)
        Ip=self.model.I(Xp,Up)

        it = 1
        while True:
            # Something done with alphas and
            while True:
                (Xn, Un, Psi) = self.improve(t, Xp, Up) #, alpha=)
                In = self.model.I(Xn, Un)
                dI = Ip-In
                if DEBUG>5:
                    print ("DI:", dI)
                if abs(dI)<eps:
                    return In, Xn, Un, it, "Opt"
                if iters<=0:
                    return In, Xn, Un, it, "Nonoptimal"
                iters-=1
                it+=1

                if In>=Ip:
                    beta/=2
                    continue
                else:
                    Xp, Up, Ip = Xn, Un, In
                    break

        raise RuntimeError("this should be not reached")


    def dU(self, t, X, U, **kwargs):
        _dHdu=self.model.dHdu(t, X, U, Psi=kwargs['Psi'])
        return _dHdu * kwargs['beta']

    def Ha(self, t, x, u, Psi, alpha):
        pt=transpose(Psi)
        return pt.dot(self.f(t,x,u))-alpha*self.f0(t,x,u)

# -------------------- tests -------------------------------------------------

class LinModel1(ParOptModel):
    """Possibly simple linear test model.
    """

    def __init__(self):
        X0=array([[1.0]])
        self.h = Ht
        # self.h = 0.001
        # self.h = 0.2
        self.num = int((1.0-0.0) / self.h)
        self.T = linspace(
            start=0.0,
            stop=1.0,
            num=self.num
        )
        self.t = arange(len(self.T))
        ParOptModel.__init__(self, X0=X0, N=1, M=1)

    def start_control(self):
        U = [array([[0.0]]) for t in self.t[:-1]]
        return array(U)

    def F(self, x):
        """ X and U are lists of vectors (arrays)
        """
        return 0.0

    def f(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        ((x0,),)=x
        ((u0,),)=u

        return [[x0+self.h*u0]]


    def f0(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        ((x0,),)=x
        ((u0,),)=u

        return self.h * (x0*x0+u0*u0)


class LinModel2d2du(ParOptModel):
    """Possibly not a simple linear test model.
    """
    def __init__(self):
        X0=array([[1.0],[1.0]])
        self.h = Ht
        self.num = int((1.0-0.0) / self.h)
        self.T = linspace(
            start=0.0,
            stop=1.0,
            num=self.num
        )
        self.t = arange(len(self.T))
        ParOptModel.__init__(self, X0=X0, N=2, M=2)

    def start_control(self):
        U = [array([[0.0],[0.0]]) for t in self.t[:-1]]
        return array(U)

    def F(self, x):
        """ X and U are lists of vectors (arrays)
        """
        return 0.0

    def f(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        ((x0,),(x1,))=x
        ((u0,),(u1,))=u

        return [[x0+self.h*u0],[x1+self.h*u1]]


    def f0(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        ((x0,),(x1,))=x
        ((u0,),(u1,))=u

        return self.h * (x0*x0+x1*x1+u0*u0+u1*u1)


def test2():
    m = LinModel1()

    ip=ParOptProcess(m)
    I, X, U, it, _ = ip.optimize(m.t, eps=0.001, iters=2000)
    print ("")
    print ("X,     U,   ")
    for x,u in zip(X,U):
        print (x,u)
    print ("Result is:", I, "in", it, "iters")
    print (_)

def test2d():
    m = LinModel2d2du()

    ip=ParOptProcess(m)
    I, X, U, it, _ = ip.optimize(m.t, eps=0.001, iters=2000)
    print ("")
    print ("X,     U,   ")
    for x,u in zip(X,U):
        print (x,u)
    print ("Result is:", I, "in", it, "iters")
    print (_)


if __name__=="__main__":
    print ("ok")

    TEST='test2d'
    LOG='restats.log'
    if PROFILE:
        import cProfile, pstats
        cProfile.run(TEST+"()", LOG)
        p=pstats.Stats(LOG)
        p.strip_dirs().sort_stats('time','cumulative').print_stats()
    else:
        eval(TEST+"()")
    quit()
