	subroutine main1(tn,tk,n,m,n1,n2,eps,h1,s0,s0f,
  	*s,v,A,B1,u1,B2,u2,d,l,u3,v1,gamma,delta,b,c,
     *ffl,nf,alff,alf,fl,fln,iterm,ch)

	implicit real*8 (a-h,o-z)
C   основная процедура
	common /sch/s0f(40)
	common /mfl/sg,beta,imm
	common /lagran/lag
	common /chut/hut(201,20)

	real*8 fl,fln,v(n1),u3(n1),s(n1,m),gamma(n),w(m),w1(n),
     *ffl(nf),alff(nf),B1(m,n),A(m,m),B2(m,n),s0,alf,uu1(201,20),
     *v1(n),u2(n,1),u1(n,1),d(m,m),l(m,m),sf0(n2),pnev0,tn,tk,g(n2),
	*delta
	integer it,ch,k,k1

c в теле основной процедуры
	beta=0.9D0
	imm=0
	if(lag.eq.1) imm=11
	epsne=eps
	if (s0.gt.0.98D0) s0=0.98D0
	sg=(1.D0-s0)/s0
	sgin=sg
	if(iflag.eq.1) goto 99999
	iflag=0
	do 555 i=1,n2
555	sf0(i)=s0f(i)

ccc пересчет исходной системы
	call sin101(tn,tk,n,n1,m,h1,v,s,n2,
 	*u3,A,B1,u1,B2,u2,d,l,v1,gamma,delta,b,c,s0,sf0,fl,fln)
	if(iflag.eq.1) goto 99999

c	stop
CCCC ВЫВОД В ФАЙЛ НАЧАЛЬНОЙ ТРАКТОРИИ И ЗНАЧЕНИЯ ФУНКЦИОНАЛА
116	format(/,15x,'НАЧАЛЬНАЯ ТРАЕКТОРИЯ УПРАВЛЕНИЕ функционал',/)
	write (77,116)
	call pechat_x_u_fl(s,v,m,n1,alf,fl,fln)

c	stop
	alfn=alf
	ffl(1)=fln
	alff(1)=alf
11	flr=fln

CCCCC ВСПОМОГАТЕЛЬНАЯ СИСТЕМА
	write (77,104)
104	format (/,'------------------------------------------------',/)	

	call new(pnev0,n2,tn,tk,n1,n,m,v,s,u3,v1)
	call kowri(sf0,n2,pnev0)
	mmc=0
	it=0
	if(imm.eq.11) write(51,17)
17	format(10x,'модифицированнная функция Лагранжа')
	iflag=0
	bx=alh
	af1=0.5D0*bx
	if(tol.le.0.) tol=0.1*eps
	st=s0
	if(imm.eq.11) st=sg
	izap=1
	af=af1
	it=1
	fp1=1.
	fp2=fp1
ccc сохранение направления спуска
	do 60 i=1,n1
	do 60 j=1,15
60	uu1(i,j)=hut(i,j)

	call sin201(tn,tk,n,n1,m,h1,v,s,n2,
	*u3,A,B1,B2,u2,d,l,v1,w1,gamma,delta,b,c,s0,sf0,				
	*psi1alf,psi2alf,alf)

c	if(ihgr.ne.1) goto 57







	fl01=fl
	fln01=fln
7	flr=fln

16	if(pnev1.lt.0.1*eps) goto 99999

99999 rewind 16
	return
	end

	subroutine sin201(tn,tk,n,n1,m,h1,v,s,n2,
	*u3,A,B1,B2,u2,d,l,v1,w1,gamma,delta,b,c,s0,sf0,				
	*psi1alf,psi2alf,alf)
CCCCCC ПОДСЧЕТ СОПРЯЖЕННОЙ СИСТЕМЫ PSI И SIGMA
CC  !!!!!!!!!!!!! НА ВЫХОД - H_uu,f_u,f'_u,H_u,transpon H_xu	
CCC				
      IMPLICIT REAL*8 (a-h,o-z)
	real*8 v(n1),s(n1,m),u3(n1),fln,tn,tk,
	*h1,t,sf0(n2),v1(n),gamma(n),s0,
	*vt,st(m,1),u3t,u2(n,1),u1(n,1),fl,
	*B1(m,n),A(m,m),B2(m,n),d(m,m),l(m,m),
     *psi1alf(n1,m),psi2alf(n1),w1(n,1),
     *alf,psi1(m,1),psi2,delta,
	*hu1y(n,1),hu3y,hs1y(m,1),hv1y,hu2y(n,1)
CC МАТРИЦЫ ВСПОМОГАТЕЛЬНЫЕ ПРИ СЧЕТЕ
	integer k0,ch,h2,nn

	common /chut/hut(201,20)
	nn=n1-1
	h2=-1
	do 2 i=1,m
2      psi1alf(n1,i)=0.D0
	psi2alf(n1)=0.D0

	q=dfloat(n1-1)
c	q=(tn-tk)/q
	q=(tk-tn)/q
	nr2=int(q/h1)
	if(nr2.le.0.) h1=q
	if(nr2.le.0.) nr2=1


	t=tk
CCCCC ЦЕПОЧКА НА {t1,t1-1,...,t0}
	k0=0
	do 1 nt=nn,1,h2
	t=t-q

CCCCC ЦЕПОЧКА НА {t,t-h}
	vt=v(nt)
	do 33 i=1,m
	psi1(i,1)=psi1alf(nt+1,i)
33	st(i,1)=s(nt,i)
C задаем управление
	u3t=u3(nt)
	psi2=psi2alf(nt+1)

CC ПОДСЧЕТ ПРОИЗВОДНой H_u в момент (t)
	call hu1p1_yavn(t,vt,psi1,psi2,n,m,h1,B1,v1,w1,
  	*gamma,sf0,n2,hu1y,alf,s0)
	call hu3p1_yavn(t,st,psi1,psi2,n,m,h1,u3t,l,b,
  	*sf0,n2,hu3y,alf,s0)
	call hu2p1_yavn(t,psi1,psi2,n,m,h1,delta,B2,
 	*sf0,n2,hu2y,alf,s0)

	do 11 i=1,n
11	hut(nt,i)=hu1y(i,1)
	do 12 i=8,14
12	hut(nt,i)=hu2y(i-n,1)
	hut(nt,15)=hu3y

	
ccc подсчет H_x
	call hvp1_yavn(t,vt,st,psi1,psi2,n,m,h1,v1,
 	*c,d,sf0,n2,hv1y,alf,s0)
	call hsp1_yavn(t,vt,st,psi1,psi2,n,m,h1,u3t,A,l,
  	*d,sf0,n2,hs1y,alf,s0)

	psi2=hv1y
	psi1=hs1y

C PSI И SIGMA В МОМЕНТ nt
	do 5 i=1,m
5	 psi1alf(nt,i)=psi1(i,1)
	psi2alf(nt)=psi2

1	continue
ccc подсчет  H_u в  конце интервала
	do 3 i=1,n
3	hut(n1,i)=hut(n1-1,i)
 
	return
	end
