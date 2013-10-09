	program Les
CCC АПРЕЛЬ 2013 - ПЕРЕДЕЛЫВАЕМ В РАБОЧИЙ ВАРИАНТ
	implicit real*8 (a-h,o-z)
	real*8 b,c,delta,s0

      real*8 epsp
	integer n,m,n1,nn,nf,iterm,n2
 	real*8 eps,fln,h1,alf,fl,rk
	real*8,allocatable ::v(:),s(:,:),alff(:),w(:),w1(:),
     *ffl(:),u3(:),A(:,:),B1(:,:),B2(:,:),gamma(:),
     *v1t(:),u2t(:,:),u1t(:,:),A11(:,:),d(:,:),l(:,:),sf0(:),
     *v1(:,:),u2(:,:),u1(:,:)
	integer, allocatable ::i1(:),f(:,:)

c	common /ng/n2
	open(51,file='les1.txt',status='replace')
	open(77,file='les1_traekt_uprav.txt',status='replace')
	open(50,file='MATRIX S.txt',status='replace')
	open(52,file='MATRIX NACH.txt',status='replace')
	open(53,file='SOPR_SYS.txt',status='replace')

ccc   A - матрица А - переходов
ccc   B1 - матрица вырубки


	tn=0.D0
C период планирования 100 лет
	tk=7.D0
C шаг времени 1 год
	h1=1.D0

C n=7 - РАЗМЕРНОСТЬ пород леса, m=i_max=35 - РАЗМЕРНОСТЬ S(t,i)
	n=7
	m=35
C m1 - размерность штрафных коэффициентов	ash(m1)
C n2- кол-во штрафов, s0 -весовой коэффициент штрафа
	n2=m+2
	s0=0.1D0
	ep1=0.1D-10
	iterm=15

C nn - РАЗМЕРНОСТЬ N - КОЛ-ВО ТОЧЕК ВРЕМЕНИ, Т.Е. t ИЗ {0,1,...,N}
C n1 - ФАКТИЧЕСКОЕ КОЛИЧЕСТВО ТОЧЕК (N+1)
	nn=idint(tk)
	n1=nn+1
	nf=17

	allocate (v(n1),u3(n1),s(n1,m),gamma(n),w(m),w1(n),
     *ffl(nf),alff(nf),B1(m,n),A(m,m),B2(m,n),sf0(n2),
     *v1t(n),u2t(n,1),u1t(n,1),A11(m,m),i1(m),f(m,m),d(m,m),l(m,m),
     *u1(n1,n),u2(n1,n),v1(n1,n))

C       МАТРИЦА А ПЕРЕХОДОВ	надо проверить!!!
	A(1:m,1:m)=0.D0
	A(1,1)=-0.0251
	A(2,2)=-0.166
	A(3,1)=0.01934
	A(3,2)=0.166
	A(3,3)=-0.0083
	A(3,15)=0.0002031
	A(3,16)=0.0060175
	A(3,21)=0.02855
	A(3,25)=0.0001
	A(3,26)=0.071527
	A(3,30)=0.00066
      A(3,35)=0.00014
	A(4,3)=0.0083
	A(4,4)=-0.025
	A(4,5)=0.0034
	A(5,4)=0.025
	A(5,5)=-0.01693
	A(6,6)=-0.188907
	A(7,7)=-0.166
	A(8,1)=0.00185
	A(8,6)=0.08059
	A(8,7)=0.166
	A(8,8)=-0.0125
	A(8,10)=0.0008
	A(8,11)=0.0028
      A(8,15)=0.0007769
 	A(8,16)=0.00069
	A(8,21)=0.0074
	A(8,26)=0.0075
	A(8,30)=0.00109
	A(8,31)=0.004999
	A(9,8)=0.0125
	A(9,9)=-0.05
	A(10,9)=0.05
	A(10,10)=-0.003852
	A(11,11)=-0.024249
	A(12,12)=-0.166
	A(13,11)=0.021388
	A(13,12)=0.166
	A(13,13)=-0.0125
	A(13,15)=0.0001
	A(13,16)=0.00233
	A(13,21)=0.0211291
	A(13,26)=0.15648
	A(13,30)=0.00009
	A(13,31)=0.02736
	A(14,13)=0.0125
	A(14,14)=-0.05
	A(15,10)=0.000104
	A(15,14)=0.05
	A(15,15)=-0.0089201
	A(16,16)=-0.0537537
	A(17,17)=-0.166
	A(18,16)=0.042035
	A(18,17)=0.166
	A(18,18)=-0.0125
	A(18,20)=0.0101
	A(18,21)=0.0786
	A(18,26)=0.030409
	A(19,18)=0.0125
	A(19,19)=-0.05
	A(19,30)=0.00696
	A(19,35)=0.0058962
	A(20,10)=0.000036
	A(20,15)=0.0029912
	A(20,19)=0.05
	A(20,20)=-0.0166246
	A(21,21)=-0.32118
	A(22,22)=-0.166
	A(23,5)=0.0135263
	A(23,15)=0.0048489
	A(23,20)=0.0065
	A(23,21)=0.0262
	A(23,22)=0.166
	A(23,23)=-0.0125
	A(23,25)=0.0026
	A(23,26)=0.089278
	A(23,31)=0.0097905
	A(24,23)=0.0125
	A(24,24)=-0.05
	A(24,30)=0.010334
	A(24,35)=0.00745
	A(25,10)=0.0029125
	A(25,24)=0.05
	A(25,25)=-0.0027005
	A(25,35)=0.00121
	A(26,26)=-0.52445
	A(27,27)=-0.166
	A(28,1)=0.0022328
	A(28,6)=0.0614
	A(28,16)=0.00156
	A(28,21)=0.089974
	A(28,26)=0.095776
	A(28,27)=0.166
	A(28,28)=-0.025
	A(28,30)=0.0006
	A(28,31)=0.00175
	A(28,28)=-0.025
	A(29,28)=0.025
	A(29,29)=-0.1
	A(30,29)=0.1
	A(30,30)=-0.019726
	A(31,31)=-0.04526
	A(32,32)=-0.166
	A(33,1)=0.00172
	A(33,6)=0.0470334
	A(33,16)=0.0012016
	A(33,21)=0.0690526
	A(33,26)=0.0735055
	A(33,31)=0.0013529
	A(33,32)=0.166
	A(33,33)=-0.025
	A(34,33)=0.025
	A(34,34)=-0.1
	A(35,34)=0.1
	A(35,35)=-0.01474

c нелегальные коффициенты

c	A(8,10)=0.0009
c	A(23,20)=0.006525
c	A(15,10)=0.001019
C       B1 - МАТРИЦА ПРИ ЛЕСОЗАГОТОВКАХ(ВЫРУБКИ)
	B1(1:m,1:n)=0.D0

	B1(1,1)=-3.805	  ! КЕДР
	B1(5,1)=3.805

	B1(6,2)=-3.1646	  ! СОСНА
	B1(10,2)=3.1646

	B1(11,3)=-3.3036  ! ЛИСТВЕННИЦА
	B1(15,3)=3.3036

	B1(16,4)=-3.932	  ! ЕЛЬ
	B1(20,4)=3.932

	B1(21,5)=-3.362	  ! ПИХТА
	B1(25,5)=3.362	! САМА ЗАПИСАЛА

	B1(26,6)=-4.759	  ! БЕРЕЗА
	B1(30,6)=4.759

	B1(31,7)=-3.2445  ! ОСИНА
	B1(35,7)=3.2445

C       B2 - МАТРИЦА ПРИ ЛЕСОПОСАДКАХ
	B2(1:m,1:n)=0.D0

	B2(1,1)=-1.	  ! КЕДР
	B2(2,1)=1.

	B2(6,2)=-1.	  ! СОСНА
	B2(7,2)=1.

	B2(11,3)=-1.  ! ЛИСТВЕННИЦА
	B2(12,3)=1.

	B2(16,4)=-1.  ! ЕЛЬ
	B2(17,4)=1.

	B2(21,5)=-1.  ! ПИХТА
	B2(22,5)=1.

	B2(26,6)=-1.  ! БЕРЕЗА
	B2(27,6)=1.

	B2(31,7)=-1.  ! ОСИНА
	B2(32,7)=1.
C b - фондоотдача ЛПХ (руб/м3), c - эксплуатац.+непроизводств. расходы (руб/м3)
C gamma(i) - доход от реализации i-го типа леса (руб/м3)
C delta - затраты на посадку 1 га любого типа леса (руб/га)
	b=31.98
	c=6.55
	delta=150.
	do 2 i=1,5
2	gamma(i)=23.1
	gamma(6)=15.7
	gamma(7)=15.7
C w(m) - средние запасы типов леса (м3/га)  набрать в верном порядке
c кедр, сосна, лиственница, ель, пихта, береза, осина
	w(1)=114.6
	w(2)=0.
	w(3)=183.4
	w(4)=222.7
	w(5)=262.8
	w(6)=34.
	w(7)=269.6
	w(8)=206.4
	w(9)=295.9
	w(10)=316.
	w(11)=96.9
	w(12)=0.
	w(13)=218.97
	w(14)=268.6
	w(15)=302.7
	w(16)=50.9
	w(17)=0.
	w(18)=150.98
	w(19)=194.8
	w(20)=254.3
	w(21)=77.6
	w(22)=0.
	w(23)=199.
	w(24)=244.7
	w(25)=297.4
	w(26)=231.7
	w(27)=0.
	w(28)=95.07
	w(29)=143.5
	w(30)=210.13
	w(31)=656.7
	w(32)=0.
	w(33)=160.
	w(34)=210.24
	w(35)=308.21
c w1(n) - средний запас леса 3-го типа(спелые и перестойные) (м3/га),
c        чтобы посчитать вырубку
	do 3 i=1,n
	j=i*5
3	w1(i)=w(j)
c v1(n) - план заготовок древесины в 3-м типе леса (м3/год)
C пока надо взять - Сосна(3-го) типа уходит в - (минус)
	v1t(1)=39130.04
	v1t(2)=975256.4
	v1t(3)=31900.09
	v1t(4)=300679.7
	v1t(5)=146416.6
	v1t(6)=746323.9
	v1t(7)=225465.3

c v1(n) - нормы вырубок в 3-м типе леса (м3/год)
c	v1(1)=39000.
c	v1(2)=312000.
c	v1(3)=31000.
c	v1(4)=301000.
c	v1(5)=146000.
c	v1(6)=253000.
c	v1(7)=232000.

CC     НАЧАЛЬНЫЕ УПРАВЛЕНИЯ u3=U_v - темпы расширения мощности (м3/год)
C   u2(n) - лесопосадки в каждом типе леса, пока не зависит от времени (га/год)
	u2t(1:n,1)=0.D0
c u1(n) - вырубки в 3-м типе леса (га/год)
	 do 4 i=1,n
4	  u1t(i,1)=v1t(i)/w1(i)
CCC ГЛОБАЛЬНЫЕ УПРАВЛЕНИЯ

	 u3(1:n1)=0.D0

	 do	41 j=1,n1
	  do 41 i=1,n
41       u1(j,i)=u1t(i,1)

	 u2(1:n1,1:n)=0.D0
	 call v1_podschet(u1,w1,n1,n,v1)
CC КОНЕЦ ЗАДАНИЯ НАЧАЛЬНЫХ УПРАВЛЕНИЙ

C d(m,m) - сведение лесов в результате пожаров по вине населения (1/(м3/год))
	d(1:m,1:m)=0.D0
	d(1,2)=-0.000000012D0
	d(1,3)=-0.000000012D0
	d(1,4)=-0.000000012D0
	d(1,5)=-0.000000012D0
	d(2,2)=0.000000012D0
	d(3,3)=0.000000012D0
	d(4,4)=0.000000012D0
	d(5,5)=0.000000012D0

	d(6,7)=-0.000000015D0
	d(6,8)=-0.000000015D0
	d(6,9)=-0.000000015D0
	d(6,10)=-0.000000015D0
	d(7,7)=0.000000015D0
	d(8,8)=0.000000015D0
	d(9,9)=0.000000015D0
	d(10,10)=0.000000015D0

	d(11,12)=-0.000000018D0
	d(11,13)=-0.000000018D0
	d(11,14)=-0.000000018D0
	d(11,15)=-0.000000018D0
	d(12,12)=0.000000018D0
	d(13,13)=0.000000018D0
	d(14,14)=0.000000018D0
	d(15,15)=0.000000018D0

	d(16,17)=-0.000000015D0
	d(16,18)=-0.000000015D0
	d(16,19)=-0.000000015D0
	d(16,20)=-0.000000015D0
	d(17,17)=0.000000015D0
	d(18,18)=0.000000015D0
	d(19,19)=0.000000015D0
	d(20,20)=0.000000015D0

	d(21,22)=-0.000000015D0
	d(21,23)=-0.000000015D0
	d(21,24)=-0.000000015D0
	d(21,25)=-0.000000015D0
	d(22,22)=0.000000015D0
	d(23,23)=0.000000015D0
	d(24,24)=0.000000015D0
	d(25,25)=0.000000015D0

	d(26,27)=-0.00000001D0
	d(26,28)=-0.00000001D0
	d(26,29)=-0.00000001D0
	d(26,30)=-0.00000001D0
	d(27,27)=0.00000001D0
	d(28,28)=0.00000001D0
	d(29,29)=0.00000001D0
	d(30,30)=0.00000001D0

	d(31,32)=-0.00000001D0
	d(31,33)=-0.00000001D0
	d(31,34)=-0.00000001D0
	d(31,35)=-0.00000001D0
	d(32,32)=0.00000001D0
	d(33,33)=0.00000001D0
	d(34,34)=0.00000001D0
	d(35,35)=0.00000001D0

c l(m,m) - матрица затрат лесной площади на расширение пр-ва	(1/(м3/год))
	l(1:m,1:m)=0.D0
	do 7 i=1,m
7	l(i,i)=0.0000000623
c A11 - запасная матрица переходов
	A11(1,3)=0.019
	A11(1,8)=0.002
	A11(1,28)=0.002
	A11(1,33)=0.002
	A11(3,4)=0.008
	A11(4,5)=0.025
	A11(5,4)=0.003
	A11(5,5)=-0.017
	A11(5,23)=0.013
	A11(6,8)=0.086
	A11(6,28)=0.061
      A11(6,33)=0.047
	A11(8,9)=0.012
	A11(9,10)=0.050
	A11(10,8)=0.001
	A11(10,10)=-0.004
	A11(9,9)=-0.050
	A11(11,11)=-0.024
	A11(13,13)=-0.012
	A11(14,14)=-0.050
	A11(16,16)=-0.040
	A11(18,18)=-0.0125
	A11(19,19)=-0.050
	A11(21,21)=-0.297
	A11(15,20)=0.003
      A11(15,23)=0.005
 	A11(16,13)=0.002
	A11(16,18)=0.004
	A11(16,28)=0.001
	A11(16,33)=0.001
	A11(18,19)=0.012
	A11(19,20)=0.05
	A11(20,18)=0.01
	A11(20,20)=-0.017
	A11(20,23)=0.006
	A11(21,3)=0.028
	A11(21,8)=0.007
	A11(21,13)=0.021
	A11(1,1)=-0.025
	A11(3,3)=-0.008
	A11(4,4)=-0.025
	A11(6,6)=-0.197
	A11(8,8)=-0.012
	A11(21,18)=0.079
	A11(21,23)=0.003
	A11(21,28)=0.089
	A11(21,33)=0.069
	A11(23,23)=-0.050
	A11(23,24)=0.050
	A11(24,24)=-0.050
	A11(24,25)=0.050
	A11(25,23)=0.002
	A11(25,25)=-0.0027
	A11(26,3)=0.072
	A11(26,8)=0.007
	A11(26,13)=0.156
	A11(26,18)=0.030
	A11(26,23)=0.089
	A11(26,26)=-0.525
	A11(26,28)=0.096
	A11(26,33)=0.073
	A11(28,28)=-0.025
	A11(10,25)=0.003
	A11(11,8)=0.003
	A11(11,13)=0.021
	A11(13,14)=0.012
	A11(15,15)=-0.089
	A11(31,31)=-0.140
	A11(31,33)=0.096
	A11(33,33)=-0.025
	A11(33,34)=0.025
	A11(34,34)=-0.01
	A11(34,35)=0.01
	A11(35,19)=0.006
	A11(35,24)=0.007
	A11(35,25)=0.001
	A11(35,35)=-0.0148
	A11(28,29)=0.025
	A11(30,8)=0.001
	A11(30,19)=0.007
	A11(30,24)=0.001
	A11(30,30)=-0.02
	A11(31,8)=0.005
	A11(31,23)=0.010
	A11(31,28)=0.002

c114	format(5D16.9)

c101	format(/,20x,'matrix A(10) ',/)
c	write (51,101)
c	do 5 i=1,m
c	write (51,114) A(10,i)
c5	continue

	do 6 i=1,n2
6	sf0(i)=1.D0
	rk=1.D0
	epsp=100.D0
	alf=1.D0
	eps=0.001D0 ! ТОЧНОСТЬ ВЫЧИСЛЕНИЯ ФУНКЦИОНАЛА

	call main2(tn,tk,nn,n,m,n1,n2,eps,h1,sf0,
     *s,v,A,B1,u1,B2,u2,d,l,u3,v1,gamma,delta,b,c,
     *ffl,nf,alf,fl,fln,iterm,rk,w1)





c	eps=0.001D0
c	f(1:m,1:m)=0
c	do 8 i=1,m
c	do 8 j=1,m
c	do 8 k=1,m
c	if (dabs(A(i,j)-A11(j,k)).le.eps) f(i,j)=f(i,j)+1
c8	continue

c	do 6 i=1,m
c	write (51,115) (f(i,k),k=1,m)
c115	format(35i3)
c6	continue

	close(77)
	close(51)
	close(50)
	close(52)
	close(53)

	end

	subroutine v1_podschet(u1,w1,n1,n,v1)
	real*8 u1(n1,n),v1(n1,n),w1(n)
CCC  РАСЧЕТ V1=U1*W1 ВО ВСЕ МОМЕНТЫ ВРЕМЕНИ
	do 1 j=1,n1
	 do 1 i=1,n
1	  v1(j,i)=u1(j,i)*w1(i)
	return
	end

	subroutine fp(n,m,vt,st,h1,u3t,A,B1,u1t,B2,u2t,d,l,f1,f)
	implicit real*8 (a-h,o-z)

C ПРАВАЯ ЧАСТЬ ДЛЯ ОДНОЙ ТОЧКИ ДИСКРЕТИЗАЦИИ в момент t
	real*8 vt,st(m,1),f(m,1),h1,u3t,u2t(n,1),u1t(n,1),
     *B1(m,n),A(m,m),B2(m,n),d(m,m),l(m,m),f1,
     *as(m,1),bu1(m,1),ds(m,1),bu2(m,1),dsv(m,1),ls(m,1),lsu3(m,1)

	 call umn_matrits(A,m,m,st,1,as)
	 call umn_matrits(B1,m,n,u1t,1,bu1)
	 call umn_matrits(B2,m,n,u2t,1,bu2)
	 call umn_matrits(d,m,m,st,1,ds)
	 dsv=vt*ds
	 call umn_matrits(l,m,m,st,1,ls)
	 lsu3=u3t*ls
C f1 - скаляр(уравнение на мощность предприятия) в момент t
C f(m,1) - матрица площадей каждого типа древостоев в момент t
	 f1=vt+h1*u3t
	 f=st+h1*(as-bu1+bu2-dsv-lsu3)
	return
	end

	subroutine f0p(n,vt,h1,u3t,u2t,v1t,gamma,delta,b,c,f0)
C ВЫЧИСЛЕНИЕ ЧАСТИ основного ФУНКЦИОНАЛА ПОД ЗНАКОМ СУММЫ В МОМЕНТ t
	implicit real*8 (a-h,o-z)

	real*8 f0,h1,vt,s1,
     *u3t,gamma(n),delta,b,c,v1t(n),u2t(n,1)

C s1 - основной функционал
	s1=0.D0
	do 1 i=1,n
1	s1=s1-gamma(i)*v1t(i)+delta*u2t(i,1)
	s1=s1+b*u3t+c*vt

	f0=h1*s1
	return
	end

	subroutine gp(n,m,g,n2,vt,st,v1t)
	implicit real*8 (a-h,o-z)
C ШТРАФЫ В МОМЕНТ (T)
C ШТРАФ ПО U3T УЧИТЫВАЕМ ПРОЕКЦИЕЙ В ОСНОВНОЙ ПРОЦЕДУРЕ
	real*8 vt,st(m,1),s3,v1t(n),g(n2)

	 s3=0.D0
	 do 1 i=1,n
1	  s3=s3+v1t(i)

	 do 2 i=1,m
2	  g(i)=-st(i,1)

	 g(36)=s3-vt
	 g(37)=-s3
	return
	end

	subroutine f0ss(n,m,h1,g,n2,vt,st,v1t,sf0,f0s)
	implicit real*8 (a-h,o-z)
C ВЫЧИСЛЕНИЕ ЧАСТИ штрафного ФУНКЦИОНАЛА ПОД ЗНАКОМ СУММЫ В МОМЕНТ t
	real*8 f0s,h1,vt,st(m,1),
     *sf0(n2),g(n2),v1t(n)

	 call gp(n,m,g,n2,vt,st,v1t)
	 f0s=0.D0

	 do 1 i=1,n2
	  sr=dmax1(g(i),0.D0)
	  f0s=f0s+sf0(i)*(sr**2)
1	 continue
	 f0s=h1*f0s
	return
	end
