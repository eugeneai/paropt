import math
import numpy
from numpy import *

class ParOptModel(object):
    """
    """

    def __init__(self, X0, A):
        """
        """
        self.X0=X0
        self.A=A

    def f(self, t, Xp, U, dt=1):
        raise RuntimeError("should be implemented by subclass")

    def F(self, X, U):
        raise RuntimeError("should be implemented by subclass")

    def start_control(self):
        raise RuntimeError("should be implemented by subclass")


class ParOptProcess(object):
    """This calss corresponds to a optimizational model.
    """

    def __init__(self, model):
        """
        """
        self.model=model

    def optimize(self, eps=0.001):
        Up=self.model.start_control()
        Xp=self.trajectory(Up)
        Fp=self.model.F(Xp,Up)
        while True:
            (Xn, Un) = self.improve(Xp,Up)
            Fn = self.model.F(Xn, Un)
            dF = Fp-Fn
            if abs(dF)<eps:
                return Fn
            Xp, Up, Fp = Xn, Un, Fn

        raise RuntimeError("this should be not reached")

    def trajectory(self, U):
        X=self.model.X0
        Xs = [X]

        for t, u in enumerate(U): # (0, u0), (1, u1)
            X=self.model.f(t, X, u)
            Xs.add(X)

        print ("U:", U)
        print ("X0:", self.model.X0)
        print ("Answer:", Xs)

        return Xs

    def improve(self, X, U):
        return (X, U)



# -------------------- tests -------------------------------------------------

class TestModel1(ParOptModel):
    """
    """
    def __init__(self):
        X0=array([0.0, 0.0])
        A=array([0.0, 1.0])
        ParOptModel.__init__(self, X0=X0, A=A)

    def start_control(self):
        return [
            array([0.1, 0.0]),
            array([0.2, 1.0])
        ]

    def f(self, t, X, U, dt=1):
        """ X ia a vector of the previous state
        """
        return X

    def F(self, X, U):
        """ X and U are lists of vectors (arrays)
        """
        return 0.1


def test1():
    m = TestModel1()

    ip=ParOptProcess(m)
    print (ip.model.F)
    print ("Result is:", ip.optimize())


if __name__=="__main__":
    print ("")
    test1()
    print ("ok")
    quit()
else:
    pass
