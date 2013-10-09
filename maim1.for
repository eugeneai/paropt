 	subroutine main2(tn,tk,nn,n,m,n1,n2,eps,h1,sf0,
  	*s,v,A,B1,u1,B2,u2,d,l,u3,v1,gamma,delta,b,c,
     *ffl,nf,alf,fl,fln,iterm,rk,w1)

	implicit real*8 (a-h,o-z)
C   НЯМНБМЮЪ ОПНЖЕДСПЮ
c	common /sch/s0f(40)
c	common /chut/hut(201,20)

	real*8 fl,fln,v(n1),u3(n1),s(n1,m),gamma(n),w1(n),
     *ffl(nf),B1(m,n),A(m,m),B2(m,n),alf,
     *d(m,m),l(m,m),sf0(n2),tn,tk,
	*delta,rk,p,psi1alf(n1,m),psi2alf(n1),
	*umin,hs,u1(n1,n),u2(n1,n),v1(n1,n),
	*un1(n1,n),un2(n1,n),vn1(n1,n),un3(n1),
	*hu1(n1,n),hu2(n1,n),hu3(n1),
	*vx(n1),sx(n1,m)


	integer nn,kol,iter

c Б РЕКЕ НЯМНБМНИ ОПНЖЕДСПШ
c	beta=0.9D0
c	if (s0.gt.0.98D0) s0=0.98D0
c	do 555 i=1,n2
c555	sf0(i)=s0f(i)

CC SF0=1 онйю врн	мювюкэмше оюпюлерпш
C	n3=15
C	rk=1.D0
C	epsp=100.D0
C	alf=1.D0

C	 пюявер мювюкэмни рпюейрнпхх х мювюкэмнцн гмювемхъ тсмйжхнмюкю - flr
	 call traekt_les(tn,nn,n,n1,m,h1,v,s,n2,
 	*u3,A,B1,u1,B2,u2,d,l,v1,gamma,delta,b,c,sf0,fl,flr,rk,p,w1)
C	 оевюрэ мювюкэмшу РПЮЕЙРНПХИ Х сопюбкемхи
	 call pechat_nachsys(n,m,n1,s,v,u1,u2,u3,v1)
C	 бшбнд гмювемхъ тсмйжхнмюкю

CC     нцпюмхвемхъ мю сопюбкемхе U3 ямхгс
	 umin=0.D0
CC     явервхй хрепюжхи
	 kol=0

	 ffl(1)=flr
	 iter=1
C%%%%%%%%%%%%%%%%%%%%%%%
CCC	жхйк лерндю опнейжхх цпюдхемрю
CCC     дн реу онп, онйю бшонкмъеряъ ABS(FLR-FLN)>EPS

CC     пеьемхе янопъфеммни яхярелш х мюунфдемхе опнхгбндмни он соп-ч
3	  call soprhu(tk,nn,n,n1,m,h1,v,s,n2,
 	*u3,A,B1,B2,u2,d,l,v1,w1,gamma,delta,b,c,sf0,				
	*psi1alf,psi2alf,alf,rk,hu1,hu2,hu3)
CC	 оевюрэ янопъфеммни яхярелш
	  call pechat_soprsys(n,n1,hu1,hu2,hu3,iter)

	  hs=1.D0 ! мювюкэмне гмювемхе HS - ьюцю яосяйю
C%%%%%%%%%%%%%%%%%%
CC	 жхйк дпнакемхъ ьюцю яосяйю

C		жхйк он бпелемх	(J)
6	    DO 11 j=1,n1
C	     бшвхякемхе нвепедмнцн опхакхфемхъ сопюбкемхъ
		 DO 10 i=1,n
		  un1(j,i)=u1(j,i)+hs*hu1(j,i) 
10		  un2(j,i)=u2(j,i)+hs*hu2(j,i)
		 un3(j)=u3(j)+hs*hu3(j) 
		  	
C		  бшвхякемхе опнейжхх сопюбкемхъ U3 мю дносярхлне лмнфеярбн
	      IF (un3(j).LT.umin) un3(j)=umin
11	    CONTINUE
C		йнмеж жхйкю он бпелемх (J)

C	 бшвхякемхе vn1, р.й. un1 - МНБНЕ СОПЮБКЕМХЕ 
	   call v1_podschet(un1,w1,n1,n,vn1)
C
C		бшвхякемхе нвепедмнцн опхакхфемхъ тюгнбни рпюейрнпхх (vx,sx)
C	     я мнбшл сопюбкемхел (un1,un2,un3,vn1)
C        х нвепедмнцн опхакхфемхъ лхмхлхгхпселнцн тсмйжхнмюкю (fln) 
		write(77,*)'тсмйжхнмюкш мю', iter,' хрепюжхх'

	    call traekt_les(tn,nn,n,n1,m,h1,vx,sx,n2,
 	*un3,A,B1,un1,B2,un2,d,l,vn1,gamma,delta,b,c,sf0,fl,fln,rk,p,w1)
C		дпнакемхе ьюцю яосяйю
		hs=hs/2.D0
	    kol=kol+1
		IF (fln.GT.flr) GOTO 6
C      йнмеж жхйкю дпнакемхъ ьюцю яосяйю
C%%%%%%%%%%%%%%%%%%
C     гюонлхмюмхе пегскэрюрнб онякедмеи хрепюжхх
	  iter=iter+1
C	  гюонлхмюмхе рпюейрнпхх
	  v=vx
	  s=sx
C	  гюонлхмюмхе сопюбкемхъ
	  u1=un1
	  u2=un2
	  u3=un3
	  v1=vn1
C	  оевюрэ опнлефсрнвмни яхярелш
	  call pechat_promsys(n,m,n1,s,v,u1,u2,u3,v1,iter)
C	  тсмйжхнмюкш
	  FUNK1=fln
	  FUNK0=flr
	  flr=fln
	  ffl(iter)=fln
	  WRITE(*,*)'FFL(',ITER,')=',FLN
	IF (dabs(FUNK1-FUNK0).GT.eps) GOTO 3
CC	 йнмеж жхйкю лерндю опнейжхх цпюдхемрю
C%%%%%%%%%%%%%%%%%%%%%%%
CCC	 оевюрэ пегскэрюрнб







CC	call new(pnev0,n2,tn,tk,n1,n,m,v,s,u3,v1,s0f)
CC	call kowri(sf0,n2,pnev0,s0f)
c16	if(pnev1.lt.0.1*eps) goto 99999
c99999 rewind 16
	return
	end


	subroutine soprhu(tk,nn,n,n1,m,h1,v,s,n2,
	*u3,A,B1,B2,u2,d,l,v1,w1,gamma,delta,b,c,sf0,				
	*psi1alf,psi2alf,alf,rk,hu1,hu2,hu3)
CCCCCC бшвхякемхе янопъфеммни яхярелш х опнхгбндмни H он сопюбкемхч
CCC				
      IMPLICIT REAL*8 (a-h,o-z)
	real*8 v(n1),s(n1,m),u3(n1),tk,
	*h1,t,sf0(n2),v1t(n),gamma(n),
	*vt,st(m,1),u3t,u2t(n,1),u1t(n,1),
	*B1(m,n),A(m,m),B2(m,n),d(m,m),l(m,m),
     *psi1alf(n1,m),psi2alf(n1),w1(n),
     *alf,psi1(m,1),psi2,delta,
	*hu1y(n,1),hu3y,hs1y(m,1),hv1y,hu2y(n,1),rk,
	*u1(n1,n),u2(n1,n),v1(n1,n),
	*hu1(n1,n),hu2(n1,n),hu3(n1)

	integer h2,nn

	do 2 i=1,m
2      psi1alf(n1,i)=0.D0
	psi2alf(n1)=0.D0

	h2=-1
	t=tk ! йнмевмне бпелъ
	call v1_podschet(u1,w1,n1,n,v1)
C	жхйк б напюрмнл бпелемх он (NT)
	do 1 nt=nn,1,h2
	 t=t-h1	 ! рейсыее бпелъ мю ндхм лемэье
C	рейсыюъ рпюейрнпхъ
	 vt=v(nt)
	 psi2=psi2alf(nt+1)
	 do 33 i=1,m
	  psi1(i,1)=psi1alf(nt+1,i)
33	  st(i,1)=s(nt,i)
C	рейсыее сопюбкемхе
	 u3t=u3(nt)	! ПЮЯЬХПЕМХЕ ЛНЫМНЯРХ ОПЕДОПХЪРХЪ
	 do 6 i=1,n
	  u2t(i,1)=u2(nt,i)	! КЕЯНОНЯЮДЙХ Б ЛНЛЕМР t
	  u1t(i,1)=u1(nt,i)	! БШПСАЙХ
6	  v1t(i)=v1(nt,i)	! БШПСАЙХ Б 3-Л РХОЕ КЕЯЮ

C	бшвхякемхе опнхгбндмни H_U б цнд (T)
	 call h_u1(t,vt,psi1,psi2,n,m,h1,B1,v1t,w1,
  	*gamma,sf0,n2,hu1y,alf,rk)
	 call h_u2(t,psi1,n,m,h1,delta,B2,hu2y,alf,rk)

	 call h_u3(t,st,psi1,psi2,n,m,h1,u3t,l,b,
  	*sf0,n2,hu3y,alf,rk)

C	гюонлхмюел H_U б цнд (T)	
	 do 11 i=1,n
	   hu1(nt,i)=hu1y(i,1)
11	   hu2(nt,i)=hu2y(i,1)
	 hu3(nt)=hu3y

C	бшвхякемхе H_X
	 call h_v(t,vt,st,psi1,psi2,n,m,h1,v1t,
 	*c,d,sf0,n2,hv1y,alf,rk)
	 call h_s(t,vt,st,psi1,psi2,n,m,h1,u3t,A,l,
  	*d,sf0,n2,hs1y,alf,rk)

	 psi2=hv1y
	 psi1=hs1y

C PSIALF2 х PSIALF1 б лнлемр (NT)
	 psi2alf(nt)=psi2
	 do 5 i=1,m
5	  psi1alf(nt,i)=psi1(i,1)
1	continue
C	йнмеж жхйкю б напюрмнл бпелемх
 
C	бшвхякемхе H_U б лнлемр (N1)
	 do 3 i=1,n
	  hu1(n1,i)=hu1(n1-1,i)
3	  hu2(n1,i)=hu2(n1-1,i)
	 hu3(n1)=hu3(n1-1)
 
	return
	end


