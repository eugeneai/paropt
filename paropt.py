import math
import numpy, sympy
from numpy import *
from functools import reduce
from pprint import pprint
import itertools
from sympy import symbols, diff

#DEBUG = 20
DEBUG = 5
Ht = 0.01

def compile_method(name, f):
    s="""
def %s(self, t, x, u):
    return %s
    """ % (name, f)

    f=compile(s, "<diffunctions>", 'exec')
    f=exec(f)

    print (s, f)

def _comp(exp):
    return compile(str(exp), "<-%s->" % str(exp), "eval")

def _eval(f, t, X, U):
    rc = eval(f, {'t':t, 'x':X, 'u':U})
    try:
        rc[0]
        return rc
    except TypeError:
        return [array([rc]) for _ in t]

class ParOptModel(object):
    """
    """

    def __init__(self, X0):
        """
        """
        self.X0=X0
        t, x, u = symbols("t, x, u")
        print (t, x, u)
        self._f=self.f(t, x, u)
        self._f0=self.f0(t, x, u)
        self._F=self.F(x)
        self._dfdx=_comp(diff(self._f, x))
        self._df0dx=_comp(diff(self._f0, x))
        self._dFdx=_comp(diff(self._F, x))

        self._dfdu=_comp(diff(self._f, u))
        self._df0du=_comp(diff(self._f0, u))
        #self.dfdu=compile_method("dfdu", self._dfdu)
        if DEBUG>=5:
            print (self._f, self._f0, self._dfdx, self._df0dx, sep=", ")
            print (self._F, self._dFdx, sep=", ")
            print (self._f, self._f0, self._dfdu, self._df0du, sep=", ")




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
        return eval(self._dFdx, {'x':xe})

    def dfdx(self, t, X, U, dt=1):
        return _eval(self._dfdx, t, X, U)

    def df0dx(self, t, X, U, dt=1):
        return _eval(self._df0dx, t, X, U)

    def dfdu(self, t, X, U, dt=1):
        return _eval(self._dfdu, t, X, U)

    def df0du(self, t, X, U, dt=1):
        return _eval(self._df0du, t, X, U)

    def Psi(self, t, X, U, alpha):

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
        return psi[::-1]

    def dHdu(self, t, X, U, Psi):
        _df0du=self.df0du(t, X, U)
        _dfdu =self.dfdu (t, X, U)
        p=Psi
        if DEBUG>10:
            print ("dfdu:", len(_dfdu))
            print ("psi:", len(p))
            print ("df0du", len(_df0du))
        _s=[dot(transpose(Psi[i]),_dfdu[i]) for i in range(len(X))]
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
                print ("DI:", dI)
                #print ("Xn:", Xn)
                #print ("Un:", Un)
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

        for t, u in enumerate(U[:-1]): # (0, u0), (1, u1)
            xn=self.model.f(t, X[t], u)
            X.append(xn)

        # print ("U:", U)
        # print ("X0:", self.model.X0)
        # print ("X:", X)
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

    def __init__(self, model):
        """
        """
        ParOptProcess.__init_(model)


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


# -------------------- tests -------------------------------------------------

class TestModel1(ParOptModel):
    """
    """
    def __init__(self):
        X0=array([0.0, 0.0])
        A=array([0.0, 1.0])
        ParOptModel.__init__(self, X0=X0)

    def start_control(self):
        return [
            array([0.1, 0.0]),
            array([0.2, 1.0])
        ]

    def F(self, x):
        """ X and U are lists of vectors (arrays)
        """
        return 0.01

    def f(self, t, x0, U, dt=1):
        """ X ia a vector of the previous state
        """
        dwdad
        return [x0 for i in t]

    def f0(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        return 0.1

    def dFdx(self, x):
        """ X and U are lists of vectors (arrays)
        """ # [[1..],[.1.],[..1]] # FIXME
        return 1.0

    def dfdx(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return 0.0

    def df0dx(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return 0.0

    def dfdu(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return 0.0

    def df0du(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return 0.0

class LinModel1(ParOptModel):
    """Possibly simple linear test model.
    """
    def __init__(self):
        X0=array([1.0])
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
        ParOptModel.__init__(self, X0=X0)

    def start_control(self):
        U = [array([0.0]) for t in self.t]
        return array(U)

    def F(self, x):
        """ X and U are lists of vectors (arrays)
        """
        return 0.0

    def f(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        return x+self.h*u


    def f0(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        return self.h * (x*x+u*u) #reduce(lambda a, e: a+e[0]**2+e[1]**2, list(zip(X,U))[:-1], 0.0)

    '''
    def dFdx(self, x):
        """ X and U are lists of vectors (arrays)
        """ # [[1..],[.1.],[..1]] # FIXME
        return 0.0

    def dfdx(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return array([array([1.0]) for t in self.t])

    def df0dx(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        if DEBUG:
            print ("!!!")
        r = X*self.h*2.0
        return r

    def dfdu(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return array([array([self.h]) for t in self.t])

    def df0du(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return 2*self.h*U
    '''

class LinModel2(LinModel1):
    """
    """

    def f0(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        return self.h * (x*x) #reduce(lambda a, e: a+e[0]**2+e[1]**2, list(zip(X,U))[:-1], 0.0)

    def dfdu(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return array([array([0.0]) for t in self.t])

    def start_control(self):
        U = [array([0.1]) for t in self.t]
        return array(U)


def test1():
    m = TestModel1()

    ip=ParOptProcess(m)
    print (ip.model.F)
    print ("Result is:", ip.optimize())

def test2():
    m = LinModel1()

    ip=ParOptProcess(m)
    print (ip.model.F)
    I, X, U, it, _ = ip.optimize(m.t, eps=0.001, iters=2000)
    print ("X,     U,   ")
    for x,u in zip(X,U):
        print (x,u)
    print ("Result is:", I, "in", it, "iters")
    print (_)


if __name__=="__main__":
    print ("")
    #import pudb; pu.db
    test2()
    print ("ok")
    quit()
else:
    pass
