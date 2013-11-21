import math
import numpy, sympy
from numpy import *
from functools import reduce
from pprint import pprint
import itertools
from sympy import symbols, diff, Symbol

#DEBUG = 20
DEBUG = 0
PROFILE = False # True
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
    for i, x in enumerate(V.x):
        g.setdefault(str(x[0]),xc[i][0])
    for i, u in enumerate(V.u):
        g.setdefault(str(u[0]),uc[i][0])

    rc = eval(f, g)
    return rc

def _vdiff(F, V, num=1): # F is a column

    DF=[]

    for f in F:
        dfs=[]
        for v in V:
            df=diff(f[0], v[0], num)
            dfs.append(df)
        DF.append(dfs)
    return DF

def _mdiff(F_v, V, num=1): # F is not transposed
    F_v_v=[]
    for row in F_v:
        F_v_v.append(_vdiff([[r] for r in row ], V))
    return F_v_v

def _diff(f, V, num=1):

    DF=[]
    for v in V:
        df=diff(f, v[0], num)
        DF.append([df])
    return DF

def _teval(d, T,X,U, V): # d is a code of function to be evaluated at all T time instance.
    rc=[]
    for t in T:
        rc.append(_eval(d, t, X[t], U[t], V))
    return array(rc)


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

        self.v.f=self.f(t, x, u)
        self.v.f0=self.f0(t, x, u)
        self.v.F=self.F(x)

        # -----------------

        self.v.f_x=_vdiff(self.v.f, x)
        self.v.f0_x=_diff(self.v.f0, x)

        self.v.f_x_x=_mdiff(self.v.f_x, x)

        if DEBUG>=5:
            print ("f:", self.v.f)
            print ("f_x:", self.v.f_x)
            print ("f_x_x:", self.v.f_x_x)

            print ("f0:", self.v.f0)
            print ("f0_x:", self.v.f0_x)

        self._f_x=_vcomp(self.v.f_x)
        self._f0_x=_vcomp(self.v.f0_x)

        self.v.F_x=_diff(self.v.F, x)
        if DEBUG>=5:
            print ("F:", self.v.F)
            print ("F_x:", self.v.F_x)
        self._F_x=_vcomp(self.v.F_x)

        self._f_u=_vcomp(_vdiff(self.v.f, u))
        self._f0_u=_vcomp(_diff(self.v.f0, u))

        if DEBUG>=5:
            print (self.v.f, self.v.f0, self._f_x, self._f0_x, sep=", ")

        if DEBUG>=5:
            print (self.v.F, self._F_x, sep=", ")
            print (self.v.f, self.v.f0, self._f_u, self._f0_u, sep=", ")
        self.__cache__={}

    def find_second_order_diffs(self):
        x=self.v.x
        u=self.v.u
        self.v.f_x_x=_mdiff(self.v.f_x, x)
        self.v.F_x_x=_vdiff(self.v.F_x, x)
        self.v.f0_x_x=_vdiff(self.v.f0_x, x)
        if DEBUG>=6:
            print ('f_x_x =', self.v.f_x_x)
            print ('F_x_x =', self.v.f_x_x)
            print ('f0_x_x=', self.v.f_x_x)
        self._f_x_x=_vcomp(self.v.f_x_x)
        self._F_x_x=_vcomp(self.v.F_x_x)
        self._f0_x_x=_vcomp(self.v.f0_x_x)
        # ---- f_u_u, f_u_x, ...
        self.v.f_u_x=_mdiff(self.v.f_u, x)
        self.v.f_u_u=_mdiff(self.v.f_u, u)
        self.v.f0_u_x=_vdiff(self.v.f0_u, x)
        self.v.f0_u_u=_vdiff(self.v.f0_u, u)
        if DEBUG>=6:
            print ('f_u_x =', self.v.f_u_x)
            print ('f_u_u =', self.v.f_u_u)
            print ('f0_u_x=', self.v.f0_u_x)
            print ('f0_u_u=', self.v.f0_u_u)

        self._f_u_x=_vcomp(self.v.f_u_x)
        self._f_u_u=_vcomp(self.v.f_u_u)
        self._f0_u_x=_vcomp(self.v.f0_u_x)
        self._f0_u_u=_vcomp(self.v.f0_u_u)

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

    def F_x(self, xe):
        return eval(self._F_x, {'x':xe, 'array':array})

    def f_x(self, T, X, U, dt=1):
        return _teval(self._f_x,  T,X,U, self.v)

    def f0_x(self, T, X, U, dt=1):
        return _teval(self._f0_x, T,X,U, self.v)

    def f_u(self, T, X, U, dt=1):
        return _teval(self._f_u, T,X,U, self.v)

    def f0_u(self, T, X, U, dt=1):
        return _teval(self._f0_u, T,X,U, self.v)

    def Psi(self, t, X, U, alpha):

        psie = -self.F_x(X[-1])

        psi=[psie]


        _f0_x=self.f0_x(t[:-1], X[:-1], U) # last element is useless
        _f_x =self.f_x (t[:-1], X[:-1], U) # last element is useless


        tt=t[:-1]

        # -----
        # rt.reverse()
        j=len(tt)-1
        p=psie

        while j>=1:
            i=tt[j]
            pp=p
            _f_x_t = transpose(_f_x[i])
            _f0_x_t = transpose(_f0_x[i])

            pn = dot(_f_x_t, pp) - alpha*_f0_x[i] # _t
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

    def H_u(self, t, X, U, Psi):

        _f0_u=self.f0_u(t[:-1], X[:-1], U)
        _f_u =self.f_u (t[:-1], X[:-1], U)
        p=Psi
        if DEBUG>10:
            print ("f_u:", len(_f_u))
            print ("psi:", len(p))
            print ("f0_u", len(_f0_u))
        _s=array([dot(_f_u[i], Psi[i]) for i in range(len(X)-1)]) # Psi[i] is shifted left to 1 step.
        return _s + _f0_u # FIXME Check dot operation.
        #return array(_s) + _f0_u # FIXME Check dot operation.

    def start_control(self):
        raise RuntimeError("should be implemented by subclass")

    #------------ These functions must be defined for Second Order Improvement Process -----

    def krot_d_dd(self, t, X, U):
        psie = -self.F_x(X[-1])
        psi=[psie]
        sige  = -self.F_x_x(X[-1])
        sig=[sige]

        _f0_x=self.f0_x(t, X, U) # last element is useless
        _f_x =self.f_x (t, X, U) # last element is useless

        tt=t[:-1]

        # -----
        # rt.reverse()
        j=len(tt)-1
        p=psie
        s=sige

        while j>=0:
            i=tt[j]
            pp=p
            sp=s
            _f_x_t = transpose(_f_x[i])
            _f0_x_t = transpose(_f0_x[i])
            _f_x_i = _f_x[i]
            _H_x_x_i = H_x_x[i]

            pn = dot(_f_x_t, pp) - _f0_x_t
            sn = dot(dot(_f_x_t, sp), _f_x_i) + _H_x_x_i

            psi.append(pn)
            sig.append(sn)
            p=pn
            s=sn
            j-=1
        # ----

        psi=array(psi)
        sig=array(sig)

        return psi[::-1], sig[::-1] # Really it is dPsidalpha
                                    # Really it is dSigmadalpha

    def f_x_x(self, T, X, U, dt=1):
        return _teval(self._f_x_x, T,X,U, self.v)

    def F_x_x(self, T, X, U, dt=1):
        return _teval(self._F_x_x, T,X,U, self.v)

    def f0_x_x(self, T, X, U, dt=1):
        return _teval(self._f0_x_x, T,X,U, self.v)

    def f_u_x(self, T, X, U, dt=1):
        return _teval(self._f_u_x, T,X,U, self.v)

    def f_u_u(self, T, X, U, dt=1):
        return _teval(self._f_u_u, T,X,U, self.v)

    def f0_u_x(self, T, X, U, dt=1):
        return _teval(self._f0_u_x, T,X,U, self.v)

    def f0_u_u(self, T, X, U, dt=1):
        return _teval(self._f0_u_x, T,X,U, self.v)

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
            Psi=self.model.Psi(t, Xp, Up, self.alpha)
            _H_u=self.model.H_u(t, Xp, Up, Psi=Psi)
            while True:
                _dU=_H_u*beta
                Un = Up + _dU
                Xn = self.trajectory(Un)

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

class SeconOrderParOptProcess(ParOptProcess):
    """A second order parametric optimization process
    """
    def __init__(self):
        X0=array([0.0, 0.0])
        self.t = arange(2)
        ParOptModel.__init__(self, X0=X0)
        self.find_second_order_diffs()


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
        _H_u=self.model.H_u(t, X, U, Psi=kwargs['Psi'])
        return _H_u * kwargs['beta']

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

class LinModel2d2du1(ParOptModel):
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

        return [[x0**5+self.h*u0],[x1**5+self.h*u1]]


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

def test2d1():
    m = LinModel2d2du1()
    return

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
