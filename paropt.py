import math
import numpy
from numpy import *
from functools import reduce

class ParOptModel(object):
    """
    """

    def __init__(self, X0):
        """
        """
        self.X0=X0

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

    def f0(self, t, X, U, dt=1):
        raise RuntimeError("should be implemented by subclass")

    def dFdx(self, xe):
        """x - ending state of a trajectory
        """
        raise RuntimeError("should be implemented by subclass")

    def dfdx(self, t, X, U, dt=1):
        raise RuntimeError("should be implemented by subclass")

    def df0dx(self, t, X, U, dt=1):
        raise RuntimeError("should be implemented by subclass")

    def dfdu(self, t, X, U, dt=1):
        raise RuntimeError("should be implemented by subclass")

    def df0du(self, t, X, U, dt=1):
        raise RuntimeError("should be implemented by subclass")

    def Psi(self, t, X, U, alpha):
        psie = -self.dFdx(X[-1])
        psi=[psie]
        _df0dx=self.df0dx(t, X, U) # last element is useless
        _dfdx =self.dfdx (t, X, U) # last element is useless
        rt=t[:-1]
        rt.reverse()

        p=psie
        for i in rt:
            p = _dfdx[i].dot(p) - alpha*_df0dx[i]
            psi.append(p)

        psi.reverse()
        return psi

    def dHdu(self, t, X, U, Psi):
        _df0du=self.df0du(t, X, U)
        _dfdu =self.dfdu (t, X, U)
        p=Psi
        return _dfdu.dot(p)+_df0du # FIXME Chack dot operation.

    def start_control(self):
        raise RuntimeError("should be implemented by subclass")


class ParOptProcess(object):
    """This calss corresponds to a optimizational model.
    """

    def __init__(self, model):
        """
        """
        self.model=model
        self.alpha=1.0
        self.beta=1.0

    def optimize(self, eps=0.001):
        Up=self.model.start_control()
        Xp=self.trajectory(Up)
        Ip=self.model.I(Xp,Up)
        while True:
            (Xn, Un) = self.improve(Xp,Up)
            In = self.model.I(Xn, Un)
            dI = Ip-In
            if abs(dI)<eps:
                return In
            Xp, Up, Ip = Xn, Un, In

        raise RuntimeError("this should be not reached")

    def trajectory(self, U):
        X=self.model.X0
        Xs = [X]

        for t, u in enumerate(U): # (0, u0), (1, u1)
            X=self.model.f(t, X, u)
            Xs.append(X)

        print ("U:", U)
        print ("X0:", self.model.X0)
        print ("X:", Xs)

        return Xs

    def improve(self, t, X, U):
        Psi=self.model.Psi(t, X, U, self.alpha)
        dHdu=self.dHdu(t, X, U)
        Un = U + dHdu * self.beta

        return self.model.f(t, X[0], Un), Un



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
        return [x0 for i in t]

    def f0(self, t, X, U, dt=1):
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
        self.h = 0.001
        self.eps = 0.001
        self.num = int((1.0-0.0) / self.h)
        self.T = linspace(
            start=0.0,
            stop=1.0,
            num=self.num
        )
        self.t = arange(len(self.T))
        ParOptModel.__init__(self, X0=X0)

    def start_control(self):
        return [array([0.0]) for t in self.t]

    def F(self, x):
        """ X and U are lists of vectors (arrays)
        """
        return 0.0

    def f(self, t, x0, U, dt=1):
        """ X ia a vector of the previous state
        """
        X=[x0]
        for t in self.t[1:]:
            X.append(X[t-1]+h*U[t-1])
        return X

    def f0(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return self.h * reduce(lambda a, (x,u): a+x*x+u*u, zip(X,U)[:-1], 0.0)

    def dFdx(self, x):
        """ X and U are lists of vectors (arrays)
        """ # [[1..],[.1.],[..1]] # FIXME
        return 0.0

    def dfdx(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return [array([1.0]) for t in self.t]

    def df0dx(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return 2*self.h*X

    def dfdu(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return [array([self.h]) for t in self.t]

    def df0du(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return 2*self.h*U

def test1():
    m = TestModel1()

    ip=ParOptProcess(m)
    print (ip.model.F)
    print ("Result is:", ip.optimize())

def test2():
    m = LinModel1()

    ip=ParOptProcess(m)
    print (ip.model.F)
    print ("Result is:", ip.optimize())


if __name__=="__main__":
    print ("")
    test2()
    print ("ok")
    quit()
else:
    pass
