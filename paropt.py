import math
import numpy
from numpy import *
from functools import reduce
from pprint import pprint
import itertools

#DEBUG = 20
DEBUG = 5
Ht = 0.01

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

    def f0(self, t, x, u, dt=1):
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

        tt=t[:-1]

        # -----
        # rt.reverse()
        j=len(tt)-1
        p=psie

        while j>=0:
            i=tt[j]
            pp=p
            pt = transpose(p)
            pn = dot(pt, _dfdx[i]) - alpha*_df0dx[i]
            psi.append(pn)
            p=pn
            j-=1
        # ----


        #import pdb; pdb.set_trace()
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

        #import pdb; pdb.set_trace()
        it = 1
        while True:
            beta = self.beta
            while True:
                (Xn, Un, Psi) = self.improve(t, Xp, Up)
                In = self.model.I(Xn, Un)
                dI = Ip-In
                print ("DI:", dI)
                #print ("Xn:", Xn)
                #print ("Un:", Un)
                if abs(dI)<eps:
                    return In, Xn, Un, Psi, it, "Opt"
                if iters<=0:
                    return In, Xn, Un, Psi, it, "Nonoptimal"
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

    def improve(self, t, X, U):
        Psi=self.model.Psi(t, X, U, self.alpha)
        _dHdu=self.model.dHdu(t, X, U, Psi)
        Un = U + _dHdu * self.beta

        #import pdb; pdb.set_trace()

        return self.trajectory(Un), Un, Psi



# -------------------- tests -------------------------------------------------

class TestModel1(ParOptModel):
    """
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

    def F(self, x):
        """ X and U are lists of vectors (arrays)
        """
        return 0.01

    def f(self, t, x0, U, dt=1):
        """ X ia a vector of the previous state
        """
        return x0

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
        return [0.0 for i in t]

    def df0dx(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return [0.0 for i in t]

    def dfdu(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return [0.0 for i in t]

    def df0du(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return [0.0 for i in t]

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
    print ("Result is:", ip.optimize(m.t))

def test2():
    m = LinModel1()

    ip=ParOptProcess(m)
    print (ip.model.F)
    I, X, U, Psi, it, _ = ip.optimize(m.t, eps=0.001, iters=2000)
    print ("X,     U,    Psi")
    for x,u,p in zip(X,U,Psi):
        print (x,u,p)
    print ("Result is:", I, "in", it, "iters")
    print (_)


if __name__=="__main__":
    print ("")
    test1()
    print ("ok")
    quit()
else:
    pass
