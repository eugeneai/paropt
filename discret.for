CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine traekt_les(tn,nn,n,n1,m,h1,v,s,n2,
	*u3,A,B1,u1,B2,u2,d,l,v1,gamma,delta,b,c,sf0,fl,fln,rk,p,w1)
CCCCCC ПРОЦЕДУРА НАХОЖДЕНИЯ ТРАЕКТОРИИ ПЛОЩАДЕЙ ЛЕСА НА {t0,t0+1,...,t1}
	implicit real*8 (a-h,o-z)
	real*8 v(n1),s(n1,m),u3(n1),fln,tn,
	*h1,t,sf0(n2),gamma(n),
	*vt,st(m,1),fl,
	*B1(m,n),A(m,m),B2(m,n),d(m,m),l(m,m),rk,p,
	*f(m,1),f1,
     *u1(n1,n),u2(n1,n),v1(n1,n),
	*u3t,u2t(n,1),u1t(n,1),w1(n)
CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CCC ЗАДАЮТСЯ НАЧАЛЬНЫЕ ДАННЫЕ
CC     НАЧАЛЬНОЕ СОСТОЯНИЕ v_0 - мощность Бирюсинского ЛПХ (м3/год)
	v(1)=2465170.D0
CC     НАЧАЛЬНОЕ СОСТОЯНИЕ s_0 - площади разных типов и видов леса (га)
	s(1,1)=328.97D0
	s(1,2)=0.D0
	s(1,3)=779.17D0
	s(1,4)=1753.03D0
	s(1,5)=7235.92D0
	s(1,6)=28217.65D0
	s(1,7)=2231.45D0
	s(1,8)=30213.66D0
	s(1,9)=18896.92D0
	s(1,10)=17434.18D0
	s(1,11)=1359.13D0
	s(1,12)=0.D0
	s(1,13)=9983.1D0
	s(1,14)=16940.43D0
	s(1,15)=55382.56D0
	s(1,16)=899.8D0
	s(1,17)=0.D0
	s(1,18)=12416.21D0
	s(1,19)=17714.94D0
	s(1,20)=50537.55D0
	s(1,21)=703.61D0
	s(1,22)=0.D0
	s(1,23)=8764.82D0
	s(1,24)=5516.96D0
	s(1,25)=22797.24D0
	s(1,26)=16089.77D0
	s(1,27)=0.D0
	s(1,28)=30959.29D0
	s(1,29)=22194.43D0
	s(1,30)=95031.17D0
	s(1,31)=134.D0
	s(1,32)=0.D0
	s(1,33)=1872.5D0
	s(1,34)=1962.04D0
	s(1,35)=32356.19D0
C НАЧАЛЬНОЕ ВРЕМЯ
	t=tn
CCCCC ЦИКЛ ПРЯМОЙ ПО ВРЕМЕНИ
	do 1 nt=1,nn
C	  ВЫЧИСЛЕНИЕ ТРАЕКТОРИИ В ГОД (T)
	  vt=v(nt)		! VT - МОЩНОСТЬ ПРЕДПРИЯТИЯ В ГОД (T)
c	  ПРОВЕРЯЕМ S НА НЕОТРИЦАТЕЛЬНОСТЬ
	  do 3 i=1,m
	    st(i,1)=s(nt,i)	! ST - ПЛОЩАДЬ ТИПОВ ЛЕСА В ГОД (T)
3	 continue
c	  КОНЕЦ ПРОВЕРКИ

C	  ЦИКЛ ПО УПРАВЛЕНИЮ
	  u3t=u3(nt)	! расширение мощности предприятия В ГОД (T)
	  do 6 i=1,n
	  u2t(i,1)=u2(nt,i)	! ЛЕСОПОСАДКИ В ГОД (T)
6	  u1t(i,1)=u1(nt,i)	! ВЫРУБКИ В ГОД (T)
C	  КОНЕЦ ЦИКЛА ПО УПРАВЛЕНИЮ	

CCCCC   ПРАВАЯ ЧАСТЬ В ГОД (T)
	  call fp(n,m,vt,st,h1,u3t,A,B1,u1t,B2,u2t,d,l,f1,f)
	  do 7 i=1,m

7	    s(nt+1,i)=f(i,1)
	  v(nt+1)=f1
	  t=t+h1   ! УВЕЛИЧИВАЕМ ВРЕМЯ НА ОДИН ГОД
1	continue


CCC ПЕЧАТЬ ТРАЕКТОРИЙ
CCCCC ПОДСЧЕТ ФУНКЦИОНАЛА
	call v1_podschet(u1,w1,n1,n,v1)
	call sfunk(tn,nn,n,m,n1,n2,v,s,h1,u3,u2,
	*v1,gamma,delta,b,c,sf0,fl,fln,rk,p)
	return
	end
CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine sfunk(tn,nn,n,m,n1,n2,v,s,h1,u3,u2,v1,
 	*gamma,delta,b,c,sf0,fl,fln,rk,p)
	implicit real*8 (a-h,o-z)

CCCCCC ПРОЦЕДУРА ПОСЧЕТА ФУНКЦИОНАЛА
CCCCCC fln - функционал	основной
C      fl - штрафной
	integer n2
	real*8 t,h1,vt,st(m,1),u3t,gamma(n),delta,b,c,
     *sf0(n2),v1t(n),u2t(n,1),fln,v(n1),s(n1,m),u3(n1),
	*tn,fl,f0,f0s,g(n2),rk,p,razn,
	*u2(n1,n),v1(n1,n)

	t=tn
	fln=0.D0
	fl=0.D0

C ПОДСЧЕТ ПОД ЗНАКОМ СУММЫ
	do 3 nt=1,nn
c формируем траекторию
	 vt=v(nt)
	 do 30 i=1,m
30	 st(i,1)=s(nt,i)
C задаем управление
	 u3t=u3(nt)
       do 6 i=1,n
        u2t(i,1)=u2(nt,i)	! лесопосадки в момент t
6	  v1t(i)=v1(nt,i)	! вырубки в 3-м типе леса

	   call f0p(n,vt,h1,u3t,u2t,v1t,gamma,delta,b,c,f0)
	   call f0ss(n,m,h1,g,n2,vt,st,v1t,sf0,f0s)
c	  print *,'f0s',f0s
C ЗНАЧЕНИЕ ФУНКЦИОНАЛА	
	 fln=fln+f0 ! FLN -ОБЫЧНЫЙ ФУНКЦИОНАЛ
	 fl=fl+f0s ! FL - ШТРАФНОЙ ФУНКЦИОНАЛ
3	 t=t+h1
c s0
c	fln=s0*fln+(1.D0-s0)*fl
	p=fl*rk/2.D0
	fln=fln+p
103	format(/,20x,'vsego fln=',D16.9,'из них  p=',D16.9)
	write(77,103) fln,p
	razn=fln-p
	WRITE(77,'(A20,D16.9)')'РАЗНОСТЬ=',razn
C102	format(/,20x,'f0 i f0s',/)
C	write(77,102)
C101	format(i3,2D16.9)
C	write(77,101) nt,f0,f0s

	return
	end
CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



C	subroutine new(pnev,n2,tn,tk,n1,n,m,v,s,u3,v1,s0f)
C	implicit real*8 (a-h,o-z)
C	integer nn
CC невязки

C	real*8 h,tn,tk,vt,st(m,1),u3t,v1(n),
C	*g(n2),v(n1),s(n1,m),u3(n1),t,s0f(n2),
C	*pnev,sr
 
C	h=(tk-tn)/dfloat(n1-1)
C	pnev=0.D0
C	t=tn
C	nn=n1-1
C	do 201 i=1,n2
C201	 s0f(i)=0.
C	do 202 nt=1,n1
C	 vt=v(nt)
C	 do 203 i=1,m
C203	  st(i,1)=s(nt,i)
C	 u3t=u3(nt)
cC	stop
C	 call gp(t,n,m,g,n2,vt,st,u3t,v1)
C	 do 204 i=1,n2
C	  sr=dmax1(0.D0,g(i))
C	  s0f(i)=s0f(i)+h*(sr**2)
C204	  pnev=pnev+h*sr
C202	 t=t+h
C	do 205 i=1,n2
C205	s0f(i)=dsqrt(s0f(i))
C	return
C	end


C	subroutine kowri(sf0,n2,rnev,s0f)
C	implicit real*8 (a-h,o-z)
C	real*8 sf0(n2),s0f(n2)
CC печать невязок


C100   format(2x,'невязка',4D10.4)
C101   format(2x,'невязки=')
C102	format(2x,'s0f',4D10.4)
C104   format(2x,'штраф=')
C105	format(2x,'sf0=',4D10.4)
C	write(51,100)rnev
C	write(51,101)
C	write(51,102)(s0f(i),i=1,n2)
C	write(51,104)
C	write(51,105)(sf0(i),i=1,n2)
C	return
C	end
