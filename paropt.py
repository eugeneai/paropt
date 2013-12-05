import math
import numpy, sympy
from numpy import *
from functools import reduce
from pprint import pprint
import itertools
from sympy import symbols, diff, Symbol
import numpy.linalg



#DEBUG = 20
DEBUG = 0
PROFILE = False # True
Ht = 0.2
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
        DF.append(tuple(dfs))
    return tuple(DF)

def _mdiff(F_v, V, num=1): # F is not transposed
    F_v_v=[]
    for row in F_v:
        F_v_v.append(_vdiff([[r] for r in row ], V))
    return tuple(F_v_v)

def _diff(f, V, num=1):

    DF=[]
    for v in V:
        df=diff(f, v[0], num)
        DF.append((df,))
    return tuple(DF)

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
            a.append((Symbol('x%s' % i),))
        self.v.x=tuple(a)
        a=[]
        for i in range(M):
            a.append((Symbol('u%s' % i),))
        self.v.u=tuple(a)

        t, x, u = self.v.t, self.v.x, self.v.u



        self.v.f=self.f(t, x, u)
        self.v.f0=self.f0(t, x, u)
        self.v.F=self.F(x)

        # -----------------

        v=self.v
        f=v.f
        f0=v.f0
        F=v.F
        v.fn={}
        self.v.f_x=_vdiff(self.v.f, x)
        print (f,x)
        v.fn[(f,x)]=v.f_x
        self.v.f0_x=_diff(self.v.f0, x)
        v.fn[(f0,x)]=v.f0_x

        self.v.f_x_x=_mdiff(self.v.f_x, x)
        v.fn[(f,x,x)]=v.f_x_x

        if DEBUG>=5:
            print ("f:", self.v.f)
            print ("f_x:", self.v.f_x)
            print ("f_x_x:", self.v.f_x_x)

            print ("f0:", self.v.f0)
            print ("f0_x:", self.v.f0_x)

        self.c=Helper()
        c=self.c
        c.fn={}
        self._f_x=_vcomp(self.v.f_x)
        c.fn[(f,x)]=self._f_x
        self._f0_x=_vcomp(self.v.f0_x)
        c.fn[(f0,x)]=self._f0_x

        self.v.F_x=_diff(self.v.F, x)
        v.fn[(F,x)]=v.F_x
        if DEBUG>=5:
            print ("F:", self.v.F)
            print ("F_x:", self.v.F_x)
        self._F_x=_vcomp(self.v.F_x)
        c.fn[(F,x)]=self._F_x

        self.v.f_u=_vdiff(self.v.f, u)
        v.fn[(f,u)]=v.f_u
        self.v.f0_u=_diff(self.v.f0, u)
        v.fn[(f0,u)]=v.f0_u
        self._f_u=_vcomp(self.v.f_u)
        c.fn[(f,u)]=self._f_u
        self._f0_u=_vcomp(self.v.f0_u)
        c.fn[(f0,u)]=self._f0_u

        if DEBUG>=5:
            print (self.v.f, self.v.f0, self._f_x, self._f0_x, sep=", ")

        if DEBUG>=5:
            print (self.v.F, self._F_x, sep=", ")
            print (self.v.f, self.v.f0, self._f_u, self._f0_u, sep=", ")
        self.__cache__={}

    def find_second_order_diffs(self):
        x=self.v.x
        u=self.v.u
        c=self.c
        v=self.v
        f=v.f
        f0=v.f0
        F=v.F
        self.v.f_x_x=_mdiff(self.v.f_x, x)
        v.fn[(f,x,x)]=v.f_x_x
        self.v.F_x_x=_vdiff(self.v.F_x, x)
        v.fn[(F,x,x)]=v.F_x_x
        self.v.f0_x_x=_vdiff(self.v.f0_x, x)
        v.fn[(f0,x,x)]=v.f0_x_x
        if DEBUG>=6:
            print ('f_x_x =', self.v.f_x_x)
            print ('F_x_x =', self.v.f_x_x)
            print ('f0_x_x=', self.v.f_x_x)
        self._f_x_x=_vcomp(self.v.f_x_x)
        c.fn[(f,x,x)]=self._f_x_x
        self._F_x_x=_vcomp(self.v.F_x_x)
        c.fn[(F,x,x)]=self._F_x_x
        self._f0_x_x=_vcomp(self.v.f0_x_x)
        c.fn[(f0,x,x)]=self._f0_x_x
        # ---- f_u_u, f_u_x, ...
        self.v.f_u_x=_mdiff(self.v.f_u, x)
        v.fn[(f,u,x)]=self.v.f_u_x
        self.v.f_u_u=_mdiff(self.v.f_u, u)
        v.fn[(f,u,u)]=self.v.f_u_u
        self.v.f0_u_x=_vdiff(self.v.f0_u, x)
        v.fn[(f0,u,x)]=self.v.f0_u_x
        self.v.f0_u_u=_vdiff(self.v.f0_u, u)
        v.fn[(f0,u,u)]=self.v.f0_u_u
        if DEBUG>=6:
            print ('f_u_x =', self.v.f_u_x)
            print ('f_u_u =', self.v.f_u_u)
            print ('f0_u_x=', self.v.f0_u_x)
            print ('f0_u_u=', self.v.f0_u_u)

        self._f_u_x=_vcomp(self.v.f_u_x)
        c.fn[(f,u,x)]=self._f_u_x
        self._f_u_u=_vcomp(self.v.f_u_u)
        c.fn[(f,u,u)]=self._f_u_u
        self._f0_u_x=_vcomp(self.v.f0_u_x)
        c.fn[(f0,u,x)]=self._f0_u_x
        self._f0_u_u=_vcomp(self.v.f0_u_u)
        c.fn[(f0,u,u)]=self._f0_u_u

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

        #for the rest of the interval
        X=X[:-1]
        t=t[:-1]

        _f0_x=self.f0_x(t, X, U) # last element is useless
        _f_x =self.f_x (t, X, U) # last element is useless


        j=len(t)-1
        p=psie

        while j>=1:
            i=t[j]
            pp=p
            _f_x_t = transpose(_f_x[i])
            #_f0_x_t = transpose(_f0_x[i])

            pn = dot(_f_x_t, pp) - alpha*_f0_x[i]

            #if j==1:
            #    print (_f_x_t, pp, _f0_x_t, "=>", pn )
            #    print (self.v[[0.4*x0], [0.4*x1]].f0_x)
            psi.append(pn)
            p=pn
            j-=1

        # ----

        psi=array(psi)

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

    def krot_d_dd(self, t, X, U, alpha=1.0):
        psi=self.Psi(t, X, U, alpha=alpha)
        v=self.v
        H_x_x=self.H((v.x,v.x), t, X, U, Psi, alpha=alpha)


        sige  = -self.fun((v.F,v.x,v.x), 0, X[-1], 0)
        sig=[sige]

        #for the rest of the interval
        X=X[:-1]
        t=t[:-1]

        _f0_x=self.fun((v.f0, v.x), t,X,U)
        _f_x =self.fun((v.f, v,x), t, X, U)
        _f0_x_x=self.fun((v.f0,v.x,v.x), t, X, U)
        _f_x_x=self.fun((v.f,v.x,v.x), t, X, U)


        j=len(t)-1
        sp=sige

        while j>=1:
            i=t[j]
            _f_x_t = transpose(_f_x[i])
            _f_x_i = _f_x[i]

            pn = psi[i]

            H_x_x_i = H_x_x[i]
            sn = dot(dot(_f_x_t, sp), _f_x_i) + H_x_x_i
            sig.append(sn)
            sp=sn
            j-=1

        sig=array(sig)
        return psi, sig[::-1] # Really it is dPsidalpha
                                    # Really it is dSigmadalpha

    def f_x_x(self, T, X, U, dt=1):
        return _teval(self._f_x_x, T,X,U, self.v)

    def F_x_x(self, xe):
        return eval(self._F_x_x, {'x':xe, 'array':array})

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

    def fun(self, vars, T, X, U):
        code=self.c.fn[vars]
        if vars[0]==self.v.F:
            return eval(code, {'x':xe, 'array':array})
        else:
            return _teval(code, T, X, U, self.v)

    def H(self, vars, T, X, U, Psi, alpha = 1.0):
        # calculate H_v_v...(t.
        assert len(vars)>0
        TT=T[:-1]
        XX=X[:-1]

        f=self.fun((self.v.f,)+vars, T, X, U)
        f0=self.fun((self.v.f0,)+vars, T, X, U)

        H = []
        for ipsi, psi in enumerate(Psi):
            H_i = - alpha * f0[ipsi] # !

            for k, p in enumerate(psi):
                # print (_p[0], _f_x_x_i[k])
                H_i += p[0]*f[ipsi][k]

                H.append(H_i)

        H = array(H)
        return H


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
    def __init__(self, model):
        X0=array([0.0, 0.0])
        self.t = arange(2)
        ParOptProcess.__init__(self, model)
        self.model.find_second_order_diffs()

    def optimize(self, t, eps=0.001, iters=1000):
        Up=self.model.start_control()
        Xp=self.trajectory(Up)
        Ip=self.model.I(Xp,Up)
        v=self.model.v


        it = 1
        tc =t[:-1]
        Xpc=Xp[:-1]

        E=numpy.identity(self.M)

        while True:
            # Something done with alphas and
            Psi, Sigma = self.model.krot_d_dd(t, Xp, Up)

            _f_u=self.model.fun((v.f,v.u), tc, Xpc, Up)
            _f_u_u=self.model.fun((v.f,v.u,v.u), tc, Xpc, Up)
            _f_u_t=transpose(_f_u)
            f_part=dot(dot(_f_u_t,Sigma),_f_u)+
                  self.model.H((v.u,v.u), tc, Xpc, Up,  alpha=alpha)-



            f_part=linalg.inv(f_part)


            Un = numpy.copy(Up)
            #dU = numpy.copy(Up)
            Xn = numpy.copy(Xp)
            #dX = numpy.copy(Xp)
            Xn[0]=Xp[0] # In reality it is already copied.

            Xn_p_i=Xp[0]
            for i, psi in enumerate(Psi):
                Xp_i=Xp[i]
                Up_i=Up[i]
                H_u_i=H_u[i]

                dX_i=Xn_p_i-Xp_i
                dU_i=f_part[i] * (H_u_a[i]+dot(s_part[i],dX_i))
                Un_i = Up_i + dU_i

                Xn_i=self.model.f(t[i], Xn_p_i, Un_i)

                Un[i]=Un_i
                Xn[i+1]=Xn_i

                Xn_p_i=Xn_i



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

        return ((x0+self.h*u0,),)


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

        return ((x0+self.h*u0,),(x1+self.h*u1,))


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

        return ((x0**5+self.h*u0,),(x1**5+self.h*u1,))


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

    #ip=ParOptProcess(m)
    ip=SeconOrderParOptProcess(m)
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
