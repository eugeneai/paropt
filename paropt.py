import math
import numpy, sympy
from numpy import *
from functools import reduce
from pprint import pprint
import itertools
from sympy import symbols, diff, Symbol
import numpy.linalg

TupleType=type((1,))


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

def _eprep(t, xc, uc, V):
    g={'t':t, 'array':array}
    for i, x in enumerate(V.x):
        g.setdefault(str(x),xc[i])
    for i, u in enumerate(V.u):
        g.setdefault(str(u),uc[i])
    return g

def _eval(f, t, xc, uc, V, g=None):
    if g == None:
        g=_eprep(t, xc, uc, V)
    rc = eval(f, g)
    return rc

def _rdiff(F, V):
    if type(F) == TupleType:
        return tuple([_rdiff(f, V) for f in F])
    return tuple([diff(F,v) for v in V])

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
            a.append(Symbol('x%s' % i))
        self.v.x=tuple(a)
        a=[]
        for i in range(M):
            a.append(Symbol('u%s' % i))
        self.v.u=tuple(a)

        t, x, u = self.v.t, self.v.x, self.v.u

        self.v.f=self.f(t, x, u)
        self.v.f0=self.f0(t, x, u)
        self.v.F=self.F(x)

        self.c=Helper()
        c=self.c
        c.fn={}
        v=self.v
        v.fn={}
        v.fn[(v.f,)]=v.f
        v.fn[(v.f0,)]=v.f0
        v.fn[(v.F,)]=v.F

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

    def Psi(self, t, X, U, alpha):

        v=self.v

        psie = -self.fun((v.F,v.x), None, X[-1], None)

        psi=[psie]

        #for the rest of the interval
        X=X[:-1]
        t=t[:-1]

        _f0_x=self.fun((v.f0,v.x), t, X, U) # last element is useless
        _f_x =self.fun((v.f,v.x), t, X, U) # last element is useless


        j=len(t)-1
        p=psie

        while j>=1:
            i=t[j]
            pp=p
            _f_x_i = _f_x[i]
            #_f0_x_t = transpose(_f0_x[i])

            pn = dot(pp,_f_x_i) - alpha*_f0_x[i]

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

        XX=X[:-1]
        tt=t[:-1]
        v=self.v
        _f0_u=self.fun((v.f0,v.u), tt, XX, U)
        _f_u =self.fun((v.f,v.u), tt, XX, U)
        p=Psi
        if DEBUG>10:
            print ("f_u:", len(_f_u))
            print ("psi:", len(p))
            print ("f0_u", len(_f0_u))

        _s=array([dot(p,fu) for p,fu in zip(Psi,_f_u)]) # Psi[i] is shifted left to 1 step.
        #_s1=dot(f_u, Psi)
        #_s=dot(Psi,_f_u)

        return _s + _f0_u # FIXME Check dot operation.

    def start_control(self):
        raise RuntimeError("should be implemented by subclass")

    #------------ These functions must be defined for Second Order Improvement Process -----

    def krot_d_dd(self, t, X, U, alpha=1.0):
        Psi=self.Psi(t, X, U, alpha=alpha)
        v=self.v


        sige  = -self.fun((v.F,v.x,v.x), 0, X[-1], 0)
        Sig=[sige]

        #for the rest of the interval
        X=X[:-1]
        t=t[:-1]


        H_x_x=self.H((v.x,v.x), t, X,U, Psi, alpha=alpha)

        _f0_x=self.fun((v.f0, v.x), t,X,U)
        _f_x =self.fun((v.f, v.x), t, X, U)
        _f0_x_x=self.fun((v.f0,v.x,v.x), t, X, U)
        _f_x_x=self.fun((v.f,v.x,v.x), t, X, U)


        j=len(t)-1
        sp=sige


        while j>=1:
            i=t[j]
            _f_x_i = _f_x[i]
            _f_x_t = transpose(_f_x_i)

            pn = Psi[i]

            H_x_x_i = H_x_x[i]
            sn = dot(dot(sp,_f_x_i), _f_x_i) + H_x_x_i
            Sig.append(sn)
            sp=sn
            j-=1

        Sig=array(Sig)
        return Psi, Sig[::-1] # Really it is dPsidalpha
                                    # Really it is dSigmadalpha

    def fun(self, vars, T, X, U):
        code=self.get_code_for(vars)
        if vars[0]==self.v.F:
            return eval(code, {'x':X, 'array':array})
        else:
            return _teval(code, T, X, U, self.v)

    def get_code_for(self, vars):
        if not vars:
            raise ValueError("no variables")
        try:
            return self.c.fn[vars]
        except KeyError:
            pass

        try:
            c=self.c.fn[vars]=_vcomp(self.v.fn[vars])
            return c
        except KeyError:
            pass

        vp=vars[:-1]
        cp=self.get_code_for(vp)
        fp=self.v.fn[vp]
        df=self.v.fn[vars]=_rdiff(fp, vars[-1])
        c=self.c.fn[vars]=_vcomp(df)
        if DEBUG>1:
            print ("A derivative needed for:\n\t", vars, "\nand it is as follows:\n\t", df)
        return c

    def H(self, vars, T, X, U, Psi, alpha = 1.0):
        # calculate H_v_v...(t.
        # print ("Evaluating H:", vars, len(T), len(X), len(U), len(Psi))
        assert len(vars)>0
        assert len(T) == len(X)
        assert len(X)==len(U)
        assert len(U)==len(Psi)


        f=self.fun((self.v.f,)+vars, T, X, U)
        f0=-self.fun((self.v.f0,)+vars, T, X, U)

        H = alpha * f0

        for psi,_H,_f,i in zip(Psi, H, f, range(len(H))):
            _H += dot(psi,_f) # !
            H[i]=_H

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

    def optimize(self, t, eps=0.001, iters=1000, alpha=1.0):
        Up=self.model.start_control()
        Xp=self.trajectory(Up)
        Ip=self.model.I(Xp,Up)
        v=self.model.v


        it = 1
        tc =t[:-1]
        Xpc=Xp[:-1]

        E=numpy.identity(self.model.M)
        _a=alpha

        while True:
            # Something done with alphas and
            Psi, Sigma = self.model.krot_d_dd(t, Xp, Up)


            _f_u=self.model.fun((v.f,v.u), tc, Xpc, Up)

            _f_u_t=_f_u.transpose(0,2,1)

            _f_x=self.model.fun((v.f,v.x), tc, Xpc, Up)
            _f_x_t=_f_x.transpose(0,2,1)

            _f_u_u=self.model.fun((v.f,v.u,v.u), tc, Xpc, Up)
            alpha=_a

            while True:
                f_part=self.model.H((v.u,v.u), tc, Xpc, Up, Psi, alpha=alpha)

                for h, fut, fu, sik, k in zip(f_part, _f_u_t, _f_u, Sigma, range(len(f_part))):
                    h+=dot(dot(fut, alpha*sik),fu)-E
                    f_part[k]=-linalg.inv(h)

                s_part=self.model.H((v.u,v.x), tc, Xpc, Up, Psi, alpha=alpha)

                for h, fxt, sik, fu in zip(s_part, _f_x_t, Sigma, _f_u):
                    s_part[k]+=dot(dot(fxt, alpha*sik),fu)

                H_u=H_u_a=self.model.H((v.u,), tc, Xpc,Up, Psi, alpha=alpha)

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
                    dU_i=dot((H_u_a[i]+dot(dX_i,s_part[i])),f_part[i])
                    Un_i = Up_i + dU_i

                    Xn_i=self.model.f(t[i], Xn_p_i, Un_i)

                    Un[i]=Un_i
                    Xn[i+1]=Xn_i

                    Xn_p_i=Xn_i


                In = self.model.I(Xn, Un)
                dI = Ip-In
                if DEBUG>=0:
                    print ("Ip:", Ip, "In:", In,  "DI:", dI)
                if abs(dI)<eps:
                    return In, Xn, Un, it, "Opt"
                if iters<=0:
                    return In, Xn, Un, it, "Nonoptimal"
                iters-=1
                it+=1

                if In>Ip:
                    alpha*=0.7
                    continue
                else:
                    Xp, Up, Ip = Xn, Un, In
                    break

# -------------------- tests -------------------------------------------------

class LinModel1(ParOptModel):
    """Possibly simple linear test model.
    """

    def __init__(self):
        X0=(1.0,)
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
        U = [(0.0,) for t in self.t[:-1]]
        return array(U)

    def F(self, x):
        """ X and U are lists of vectors (arrays)
        """
        return 0.0

    def f(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        x0=x[0]
        u0=u[0]

        return (x0+self.h*u0,)


    def f0(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        x0=x[0]
        u0=u[0]

        return self.h * (x0*x0+u0*u0)


class LinModel2d2du(ParOptModel):
    """Possibly not a simple linear test model.
    """
    def __init__(self):
        X0=(1.0,1.0)
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
        U = [[0.0,0.0] for t in self.t[:-1]]
        return array(U)

    def F(self, x):
        """ X and U are lists of vectors (arrays)
        """
        return 0.0

    def f(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """
        (x0,x1)=x
        (u0,u1)=u

        return (x0+self.h*u0, x1+self.h*u1)
#        return ((x0+self.h*u0,),(x1+self.h*u1,))


    def f0(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """

        (x0,x1)=x
        (u0,u1)=u
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

def gnuplot(fn, *args):
    # args are I, X, U, it, _ (result: Opt, etc.)
    X=args[0][1]
    o=open(fn+".dat", "w")
    print (args)
    for i, x in enumerate(X[:-1]):
        o.write("\t".join([str(a) for a in x ]))
        o.write("\t")
        for  I, X, U, it, _ in args:
            o.write("\t".join([str(a) for a in U[i]]))
            o.write("\t")
        o.write("\n")
    o.close()

def test2(so=True):
    m = LinModel1()

    if not so:
        ip=ParOptProcess(m)
    else:
        ip=SeconOrderParOptProcess(m)

    I, X, U, it, _ = ip.optimize(m.t, eps=0.001, iters=2000)
    print ("")
    print ("X,     U,   ")
    for x,u in zip(X,U):
        print (x,u)
    print ("Result is:", I, "in", it, "iters")
    print (_)

def test2d(so=True):
    m = LinModel2d2du()

    if not so:
        ip=ParOptProcess(m)
    else:
        ip=SeconOrderParOptProcess(m)
    I, X, U, it, _ = ip.optimize(m.t, eps=0.001, iters=200)
    print ("")
    print ("X,     U,   ")
    for x,u in zip(X,U):
        print (x,u)
    print ("Result is:", I, "in", it, "iters")
    print (_)

def test_with_plot():
    m = LinModel1()
    p1=ParOptProcess(m)
    p2=SeconOrderParOptProcess(m)
    iters=2000
    eps=0.001
    print ("First process:")
    r1=I1, X1, U1, it1, _1 = p1.optimize(m.t, eps=eps, iters=iters)
    print ("Second process:")
    r2=I2, X2, U2, it2, _2 = p2.optimize(m.t, eps=eps, iters=iters)

    gnuplot("test_with_plt", r1,r2)

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

    TEST='test_with_plot'
    LOG='restats.log'
    if PROFILE:
        import cProfile, pstats
        cProfile.run(TEST+"()", LOG)
        p=pstats.Stats(LOG)
        p.strip_dirs().sort_stats('time','cumulative').print_stats()
    else:
        eval(TEST+"()")
    quit()
