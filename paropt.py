import math
import numpy, sympy
from numpy import *
from functools import reduce
from pprint import pprint
import itertools
from sympy import symbols, diff, Symbol
import numpy.linalg
import unittest
from sympy.utilities.lambdify import lambdify
import os

TupleType=type((1,))
ListType=type([])


DEBUG = 20
#DEBUG = 0
PROFILE = False # True
#Ht = 0.01
Ht = 0.2
import time

class Helper():
    pass

class VFCalc(object):
    """ Class automizing differential convertions of
    vector-functions.
    """

    def __init__(self, N,M):
        """
        """
        self.N=N
        self.M=M
        self.v=Helper()
        self.v.x=[Symbol('x'+str(i+1)) for i in range(self.N)]
        self.v.u=[Symbol('u'+str(i+1)) for i in range(self.M)]
        self.v.t=Symbol('t')

    def diff1(self, f, var):# derivation of one vector-variable
        """This is method for figuring out of
        functional differentials of vector-variable lists
        """
        if type(f) in [TupleType,ListType]: #f myght be a tuple(n-ka)
            df=tuple([self.diff1(fi, var) for fi in f])
        else:
            df=tuple([diff(f, vi) for vi in var])
        if len(df)==1:
            df=df[0]    # Original f is not a tuple or a list (!)
        return df

    def diff(self, f, *vars):
        """This is method for figuring out of
        functional differentials of vector-variable lists
        """
        cf=f#current func
        for v in vars:
            cf=self.diff1(cf, v)
        return cf

    def subs(self, f, s):
        if type(f) not in [TupleType,ListType]:
            return f.subs(s) # common subs
        return tuple([self.subs(fi,s) for fi in f])#local subs

    def lambdify(self, f):
        """Compiles function f into a lambda-function
        with args as its arguments.
        """
        l=[self.v.t]
        l.extend(self.v.x)
        l.extend(self.v.u)
        # print (l)
        fl=lambdify(l, f, "numpy")
        """
        def _prepareags(t,X,U):
            if X.ndim>1:
                Xs=[X[:,i:i+1] for i in range(X.shape[1])]
                Us=[U[:,i:i+1] for i in range(U.shape[1])]
                return (t,)+tuple(Xs+Us)
            else:
                return (t,)+tuple(X)+tuple(U)
        def _fgen(t, X, U):
            args=_prepareags(t,X,U)
            ff=fl(*args)
            return ff
        return _fgen
        """
        return fl

    def code(self, f, *vars, debug=False):
        df=self.diff(f, *vars)
        c=self.lambdify(df)
        if debug: return c, df, f
        return c,df


class Model(object):
    """Parametric Optimised model class.
    See docs for individual methods on usage ways.
    """

    def __init__(self, X0, U0):
        self.X0=array(X0)
        self.U0=U0
        self.N=self.X0.shape[0]   # Dimention of X
        self.M=U0.shape[1]   # Dimention of U
        print ('X0:', X0)
        print ('U0:', U0)

    def F(self, x):
        """x - ending state of a trajectory
        """
        return 0.0

    def f(self, t, x, U, dt=1):
        raise RuntimeError("should be implemented by subclass")

    def f0(self, t, x, u, dt=1):
        return 0.0


    #------------ These functions must be defined for Second Order Improvement Process -----


class Process(VFCalc):
    """This calss corresponds to a optimizational model.
    """

    def __init__(self, model, alpha=1.0, beta=1.0):
        """
        """
        VFCalc.__init__(self, model.N, model.M)
        t,x,u=self.v.t,self.v.x,self.v.u
        self.model=model
        self.v.f=model.f(t, x, u)
        self.v.f0=model.f0(t, x, u)
        self.v.F=model.F(x)
        self.alpha=alpha
        self.beta=beta

    def optimize(self, t, eps=0.001, iters=1000):

        Up=self.model.U0

        Xp=self.trajectory(Up)
        Ip=self.I(Xp,Up)

        it = 1

        while True:
            beta = self.beta


            Psi=self.Psi(t, Xp, Up, self.alpha)
            _H_u=self.H((self.v.u,), t[:-1], Xp[:-1], Up, Psi, alpha=1.0)
            while True:
                print (it, beta)
                _dU=_H_u*beta
                Un = Up + _dU
                Xn = self.trajectory(Un)

                In = self.I(Xn, Un)

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
            xn=self.model.f(t, X[t], u) # FIXME: implement via code
            X.append(xn)

        return array(X)

    def I(self, X, U):
        def _add(acc, t):
            return acc + self.model.f0(t, X[t], U[t]) # FIXME: implement via code
        return reduce(_add, range(len(X)-1), self.model.F(X[-1])) # FIXME: implement via code

    def Psi(self, t, X, U, alpha):

        v=self.v

        # import pudb; pu.db

        print (v.F, X)
        print ("-------", t[-1], X[-1], U[-1])


        psie = -self.fun(v.F,(v.x,), t[-1], X[-1], U[-1])

        psi=[psie]

        #for the rest of the interval
        X=X[:-1]
        t=t[:-1]

        _f0_x=self.fun(v.f0, (v.x,), t, X, U) # last element is useless
        _f_x =self.fun(v.f, (v.x,), t, X, U) # last element is useless

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

        #return array([psi[::-1]]).T
        return psi[::-1]

    def krot_d_dd(self, t, X, U, alpha=1.0):
        Psi=self.Psi(t, X, U, alpha=alpha)
        v=self.v


        sige  = -self.fun(v.F, (v.x,v.x), 0, X[-1], 0)
        Sig=[sige]

        #for the rest of the interval
        X=X[:-1]
        t=t[:-1]


        H_x_x=self.H((v.x,v.x), t, X,U, Psi, alpha=alpha)

        _f0_x=self.fun(v.f0, (v.x,), t,X,U)
        _f_x =self.fun(v.f, (v.x,), t, X, U)
        _f0_x_x=self.fun(v.f0, (v.x,v.x), t, X, U)
        _f_x_x=self.fun(v.f, (v.x,v.x), t, X, U)


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

    def fun(self, f, vars, T, X, U):
        # evaluate derivatives of f on vars and substitute t, X,U
        code,df=self.code(f, *vars)
        if X.ndim>1:
            Xs=[X[:,i:i+1] for i in range(X.shape[1])]
            Us=[U[:,i:i+1] for i in range(U.shape[1])]
            m,args=True,(T,)+tuple(Xs+Us)
        else:
            m,args=False,(T,)+tuple(X)+tuple(U)
        rc=code(*args)
        if type(rc) is TupleType:
            rc=array(rc)
            if m:
                rc=rc.T
                rc=rc[0]
        if type(T)==numpy.ndarray:
            dm=False
            if type(rc)!=numpy.ndarray:
                dm=(1,)
            elif T.shape[0]!=rc.shape[0]:
                dm=rc.shape
            if dm:
                nff=numpy.zeros((len(T),)+dm,dtype=float)
                nff[:]=rc
                rc=nff
        return rc

    def H(self, vars, T, X, U, Psi, alpha = 1.0):
        # calculate H_v_v...(t.
        # print ("Evaluating H:", vars, len(T), len(X), len(U), len(Psi))
        assert len(vars)>0
        assert len(T) == len(X)
        assert len(X)==len(U)
        assert len(U)==len(Psi)


        f=self.fun(self.v.f, vars, T, X, U)
        f0=-self.fun(self.v.f0, vars, T, X, U)

        H = alpha * f0

        for psi,_H,_f,i in zip(Psi, H, f, range(len(H))):
            _H += dot(psi,_f) # !
            H[i]=_H

        return H

class SeconOrderProcess(Process):
    """A second order parametric optimization process
    """
    def __init__(self, model, **kwargs):
        X0=array([0.0, 0.0])
        self.t = arange(2)
        Process.__init__(self, model, **kwargs)

    def optimize(self, t, eps=0.001, iters=1000, alpha=None):
        Up=self.model.start_control()
        Xp=self.trajectory(Up)
        Ip=self.I(Xp,Up)
        #v=self.model.v


        it = 1
        tc =t[:-1]
        Xpc=Xp[:-1]

        E=numpy.identity(self.model.M)
        if alpha==None:
            _a=self.alpha
        else:
            _a=alpha

        while True:
            # Something done with alphas and
            Psi, Sigma = self.krot_d_dd(t, Xp, Up)
            _f_u=self.fun(v.f, (v.u,), tc, Xpc, Up)

            _f_u_t=_f_u.transpose(0,2,1)

            _f_x=self.fun(v.f, (v.x,), tc, Xpc, Up)
            _f_x_t=_f_x.transpose(0,2,1)

            _f_u_u=self.fun(v.f,(v.u,v.u), tc, Xpc, Up)
            alpha=_a

            while True:

                #import pudb; pu.db
                f_part=self.H((v.u,v.u), tc, Xpc, Up, Psi, alpha=alpha)

                for h, fut, fu, sik, k in zip(f_part, _f_u_t, _f_u, Sigma, range(len(f_part))):
                    h+=dot(dot(fut, alpha*sik),fu)-E
                    f_part[k]=-linalg.inv(h)

                s_part=self.H((v.u,v.x), tc, Xpc, Up, Psi, alpha=alpha)

                for h, fxt, sik, fu in zip(s_part, _f_x_t, Sigma, _f_u):
                    s_part[k]+=dot(dot(fxt, alpha*sik),fu)

                H_u=H_u_a=self.H((v.u,), tc, Xpc,Up, Psi, alpha=alpha)

                Un = numpy.copy(Up)
                #dU = numpy.copy(Up)
                Xn = numpy.copy(Xp)
                #dX = numpy.copy(Xp)
                Xn[0]=Xp[0] # In reality it is already copied.

                Xn_i=Xp[0]

                for i, psi in enumerate(Psi):
                    Xp_i=Xp[i]
                    Up_i=Up[i]
                    H_u_i=H_u[i]

                    dX_i=Xn_i-Xp_i
                    #dU_i=dot((H_u_a[i]+dot(dX_i,s_part[i])),f_part[i])
                    #dU_i=dot(H_u_a[i],f_part[i])
                    dU_i=H_u_a[i]
                    Un_i = Up_i + dU_i

                    Xnn_i=self.model.f(t[i], Xn_i, Un_i)

                    Un[i]=Un_i
                    Xn[i+1]=Xnn_i

                    Xn_i=Xnn_i


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

class LinModel1(Model):
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
        Model.__init__(self, X0=X0, U0=self.start_control())

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


class LinModel2d2du(Model):
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
        Model.__init__(self, X0=X0, U0=self.start_control())

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

    def f0(self, t, x, u, dt=1):
        """ X ia a vector of the previous state
        """

        (x0,x1)=x
        (u0,u1)=u
        return self.h * (x0*x0+x1*x1+u0*u0+u1*u1)

class LinModel2d2du1(Model):
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
        Model.__init__(self, X0=X0, U0=self.start_control())

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
    for i, (I, X, U, it, _) in enumerate(args):
        o=open(fn+str(i)+".dat", "w")
        lx=len(X)
        for i,x,u in zip(range(lx-1), X[:-1],U):
            t=i/lx
            o.write(str(t)+"\t")
            o.write("\t".join([str(a) for a in x]))
            o.write("\t")
            o.write("\t".join([str(a) for a in u]))
            o.write("\n")
        o.close()
    os.system('gnuplot test_with_plt.p')
    os.system('evince my-plot.ps')

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

def test2d(so=False):
    m = LinModel2d2du()

    if not so:
        ip=Process(m)
    else:
        ip=SeconOrderProcess(m)
    I, X, U, it, _ = ip.optimize(m.t, eps=0.001, iters=200)
    print ("")
    print ("X,     U,   ")
    for x,u in zip(X,U):
        print (x,u)
    print ("Result is:", I, "in", it, "iters")
    print (_)

def test_with_plot():
    m = LinModel1()
    p1=Process(m, beta=1.0)
    p2=SeconOrderProcess(m, alpha=0.001)
    iters=2000
    eps=0.001
    print ("First-order process:")
    print ("-"*20)
    r1=I1, X1, U1, it1, _1 = p1.optimize(m.t, eps=eps, iters=iters)
    print (I1, "iters:", it1)
    print ("Second-order process:")
    print ("-"*20)
    r2=I2, X2, U2, it2, _2 = p2.optimize(m.t, eps=eps, iters=iters)
    print (I2, "iters:", it2)

    gnuplot("test_with_plt", r1,r2)
    #gnuplot("test_with_plt", r1)

def test2d1():
    m = LinModel2d2du1()
    return

    ip=Process(m)
    I, X, U, it, _ = ip.optimize(m.t, eps=0.001, iters=2000)
    print ("")
    print ("X,     U,   ")
    for x,u in zip(X,U):
        print (x,u)
    print ("Result is:", I, "in", it, "iters")
    print (_)

def test_VFCalc():
    d=VFCalc(2,2)
    x1,x2=Symbol('x1'),Symbol('x2')
    u1,u2=Symbol('u1'),Symbol('u2')
    y1=x1**2*u1+x2*u2**2
    y2=x1**2*x2**2*u1**2*u2**2

    print ()
    #print (d.diff1([y1,y2],[x1,x2]))
    #print (d.diff1([y1,y2],[u1,u2]))
    #res=(d.diff([y1,y2], [x1,x2], [u1,u2]))
    #pprint (res)

    X1=ones(10)
    X2=X1
    X2+=1
    X=ndarray(shape=(10,2), dtype=float, order='F')
    X[:,0]=X1
    X[:,1]=X2

    U1=X1+2
    U2=U1+4
    U=ndarray(shape=(10,2), dtype=float, order='F')
    U[:,0]=U1
    U[:,1]=U2

    T=linspace(start=0.0, stop=1.0, num=len(X1))

    fun, f=d.code([y1,y2], [x1,x2],[u1,u2], debug=True)

    print (fun, f)
    print ((T[0],X[0,:],U[0,:]))
    pprint (fun(T[0],X[0,:],U[0,:]))

    quit()


if __name__=="__main__":
    #test_VFCalc()
    #raise SystemExit(0)

    print ("ok")

    TEST='test_with_plot'
    #TEST='test2d'
    LOG='restats.log'
    if PROFILE:
        import cProfile, pstats
        cProfile.run(TEST+"()", LOG)
        p=pstats.Stats(LOG)
        p.strip_dirs().sort_stats('time','cumulative').print_stats()
    else:
        eval(TEST+"()")
    quit()
