CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine h_v(t,vt,st,psi1,psi2,n,m,h1,v1t,
	*c,d,sf0,n2,hv1y,alf,rk)
	real*8 vt,st(m,1),psi1(m,1),psi2,hv1y,alf,h1,t,v1t(n),
	*d(m,m),c,sf0(n2),r(1,m),r1(1,1),psit(1,m),pr,rk
CCCCCC ЯВНАЯ ПРОИЗВОДНАЯ ГАМИЛЬТОНИАНА с АЛЬФА ПО v(t) - скаляр
	
       call transpon(psi1,m,1,psit)

       call umn_matrits(psit,1,m,d,m,r) !PSI'*D=R
       call umn_matrits(r,1,m,st,1,r1) !PSI'*D*ST=R1
	 pr=0.D0
	 do 1 i=1,n
1	  pr=pr+v1t(i) ! СУММА ПО V1T(I=1,7)
	 pr=alf*h1*(sf0(36)*2.D0*(pr-vt)*(-1.D0))
c s0
c	hv1y=psi2-r1(1,1)-alf*h1*s0*c-(1.D0-s0)*pr
c rk
	 hv1y=psi2-r1(1,1)-alf*h1*c-rk*pr/2.D0
	return
	end
CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine h_s(t,vt,st,psi1,psi2,n,m,h1,u3t,A,l,
 	*d,sf0,n2,hs1y,alf,rk)
	real*8 A(m,m),vt,st(m,1),psi1(m,1),psi2,hs1y(m,1),
	*alf,h1,t,u3t,l(m,m),e(m,m),dv(m,m),lu3(m,m),
	*d(m,m),sf0(n2),r(m,1),r2(m,m),pr(m,1),r2t(m,m),rk
CCCCCC ЯВНАЯ ПРОИЗВОДНАЯ ГАМИЛЬТОНИАНА с АЛЬФА ПО st(t)
	e(1:m,1:m)=0.D0
	do 2 i=1,m
2	e(i,i)=1.D0
	dv=vt*d
      lu3=u3t*l
	r2=e+h1*(A-dv-lu3)
      call transpon(r2,m,m,r2t)

	call umn_matrits(r2t,m,m,psi1,1,r)

	do 1 i=1,m
1	pr(i,1)=-2.D0*sf0(i)*(-st(i,1))

c s0
c	hs1y=r+alf*h1*(1.D0-s0)*pr
c rk
	hs1y=r+alf*h1*pr*rk/2.D0
	return
	end
CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine h_u1(t,vt,psi1,psi2,n,m,h1,B1,v1t,w1,
 	*gamma,sf0,n2,hu1y,alf,rk)
	real*8 B1(m,n),gamma(n),vt,psi1(m,1),psi2,hu1y(n,1),
	*alf,h1,t,v1t(n),w1(n),rk,
	*sf0(n2),B1t(n,m),r(n,1),r2(n,1),pr,r3(n,1),r4(n,1),pr1
CCCCCC ЯВНАЯ ПРОИЗВОДНАЯ ГАМИЛЬТОНИАНА с АЛЬФА ПО u1(t)
      call transpon(B1,m,n,B1t)
	call umn_matrits(B1t,n,m,psi1,1,r)
	r=-h1*r

	do 2 i=1,n
c s0
c2	r2(i,1)=alf*h1*s0*(-gamma(i))*w1(i)
2	r2(i,1)=alf*h1*(-gamma(i))*w1(i)

	pr=0.D0
	do 1 i=1,n
1	pr=pr+v1t(i)
	pr1=pr-vt

	do 4 i=1,n
4	r3(i,1)=sf0(36)*2.D0*pr1*w1(i)

	do 5 i=1,n
5	r4(i,1)=sf0(37)*2.D0*pr*w1(i)
c s0
c	hu1y=r-r2-alf*h1*(1.D0-s0)*(r3+r4)
c rk
	hu1y=r-r2-alf*h1*(r3+r4)*rk/2.D0

	return
	end
CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine h_u3(t,st,psi1,psi2,n,m,h1,u3t,l,b,
 	*sf0,n2,hu3y,alf,rk)
	real*8 st(m,1),psi1(m,1),psi2,hu3y,b,
	*alf,h1,t,u3t,l(m,m),
	*sf0(n2),r(m,1),r1(1,1),pr,rt(1,m),rk
CCCCCC ЯВНАЯ ПРОИЗВОДНАЯ ГАМИЛЬТОНИАНА с АЛЬФА ПО u3(t) - скаляр

	call umn_matrits(l,m,m,st,1,r)
      call transpon(r,m,1,rt)

	call umn_matrits(rt,1,m,psi1,1,r1)

CC	pr=sf0(38)*2.D0*u3t	! УБИРАЕМ, БУДЕМ ЯВНО УЧИТЫВАТЬ
c s0
c	hu3y=h1*psi2-r1(1,1)-alf*h1*(s0*b+(1.D0-s0)*pr)
c rk
C	hu3y=h1*psi2-r1(1,1)-alf*h1*(b+rk*pr/2.D0)
	hu3y=h1*psi2-h1*r1(1,1)-alf*h1*b ! БЕЗ ШТРАФА

	return
	end
CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine h_u2(t,psi1,n,m,h1,delta,B2,hu2y,alf,rk)
	real*8 psi1(m,1),hu2y(n,1),B2(m,n),
	*alf,h1,t,delta,B2t(n,m),r1(n,1),del1(n,1),rk
CCCCCC ЯВНАЯ ПРОИЗВОДНАЯ ГАМИЛЬТОНИАНА с АЛЬФА ПО u2(t) - вектор

      call transpon(B2,m,n,B2t)
	call umn_matrits(B2t,n,m,psi1,1,r1)

	r1=h1*r1
	do 1 i=1,n
c s0
c1	del1(i,1)=-alf*h1*s0*delta
1	del1(i,1)=-alf*h1*delta

	hu2y=r1+del1
	return
	end
CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
