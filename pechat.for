	subroutine pechat_func(p,fln)
CCCC	 оевюрэ тсмйжхнмюкнб
	real*8 fln,p
102	format(/,20x,'гМЮВЕМХЕ ТСМЙЖХНМЮКЮ FLN, p - ЬРПЮТМНЦН')
115	format(2D16.9)
	write (77,102) 
	write (77,115) fln,p
	return
	end

CCCCCCC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	SUBROUTINE PECHAT11(Y,K,L)
	REAL*8 Y(K,L)
CCC оевюрэ люрпхжш я 8-ч ярнкажюлх
106   FORMAT(8(D16.9,' '))
      DO 331 I=1,K
331	 WRITE (50, 106) Y(I,1),Y(I,2),Y(I,3),
     *Y(I,4),Y(I,5),
     *Y(I,6),Y(I,7),Y(I,8)
C     *Y(I,9),Y(I,10),Y(I,11),
C     *Y(I,12),Y(I,13),Y(I,14),Y(I,15),Y(I,16),Y(I,17),
C     *Y(I,18),Y(I,19),Y(I,20),Y(I,21)
	return
	end
CCCCCCC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	SUBROUTINE PECHATnach(Y,K,L)
	REAL*8 Y(K,L)
CCC оевюрэ люрпхжш я 8-ч ярнкажюлх
106   FORMAT(8(D16.9,' '))
      DO 331 I=1,K
331	 WRITE (52, 106) Y(I,1),Y(I,2),Y(I,3),
     *Y(I,4),Y(I,5),
     *Y(I,6),Y(I,7),Y(I,8)
C     *Y(I,9),Y(I,10),Y(I,11),
C     *Y(I,12),Y(I,13),Y(I,14),Y(I,15),Y(I,16),Y(I,17),
C     *Y(I,18),Y(I,19),Y(I,20),Y(I,21)
	return
	end
CCCCCCC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	SUBROUTINE PECHAT_sopr(Y,K,L)
	REAL*8 Y(K,L)
CCC оевюрэ люрпхжш я 8-ч ярнкажюлх
106   FORMAT(8(D16.9,' '))
      DO 331 I=1,K
331	 WRITE (53, 106) Y(I,1),Y(I,2),Y(I,3),
     *Y(I,4),Y(I,5),
     *Y(I,6),Y(I,7),Y(I,8)
C     *Y(I,9),Y(I,10),Y(I,11),
C     *Y(I,12),Y(I,13),Y(I,14),Y(I,15),Y(I,16),Y(I,17),
C     *Y(I,18),Y(I,19),Y(I,20),Y(I,21)
	return
	end
CCCCCCC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine pechat_soprsys(n,n1,hu1,hu2,hu3,iter)
	REAL*8 hu1(n1,n),hu2(n1,n),hu3(n1),
	*hu1tr(n,n1),hu2tr(n,n1)

C	оевюрэ янопъфеммни яхярелш
106    FORMAT(8(D16.9,' '))
         call transpon(hu1,n1,n,hu1tr)
         call transpon(hu2,n1,n,hu2tr)
	  write(53,'(A20,i2)')'H_U1 ХРЕПЮЖХЪ', iter
	   call PECHAT_sopr(hu1tr,n,n1)
	  write(53,'(A20)')'H_U2'
	   call PECHAT_sopr(hu2tr,n,n1)
	  WRITE (53,'(A20)')'H_U3'
	   WRITE (53,106) hu3(1),hu3(2),hu3(3),
     *hu3(4),hu3(5),hu3(6),hu3(7),hu3(8)
	return
	end
CCCCCCC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine pechat_nachsys(n,m,n1,s,v,u1,u2,u3,v1)
	real*8	v(n1),u3(n1),s(n1,m),
	*u1(n1,n),u2(n1,n),v1(n1,n),stt(m,n1),
 	*u1tr(n,n1),u2tr(n,n1),v1tr(n,n1)

C оевюрэ мювюкэмни яхярелш
	   call transpon(s,n1,m,stt)
 	 WRITE (52,'(A20)')'МЮВЮКЭМЮЪ люрпхжю S'
	   call PECHATnach(stt,m,n1)
	 WRITE (52,'(A20)')'МЮВЮКЭМЮЪ люрпхжю V'
106   FORMAT(8(D16.9,' '))
	 WRITE (52, 106) v(1),v(2),v(3),
     *v(4),v(5),v(6),v(7),v(8)
	 write(52,'(A20)')'мювюкэмше сопюбкемхъ'
        call transpon(u1,n1,n,u1tr)
        call transpon(u2,n1,n,u2tr)
	  call transpon(v1,n1,n,v1tr)
	 write(52,'(A20)')'мювюкэмне U1'
        call PECHATnach(u1tr,n,n1)
	 write(52,'(A20)')'мювюкэмне U2'
	  call PECHATnach(u2tr,n,n1)
	 write(52,'(A20)')'мювюкэмне V1'
	  call PECHATnach(v1tr,n,n1) 
	 WRITE (52,'(A20)')'мювюкэмне U3'
	  WRITE (52,106) u3(1),u3(2),u3(3),
     *u3(4),u3(5),u3(6),u3(7),u3(8)

	return
	end
CCCCCCC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CCCCCCC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine pechat_promsys(n,m,n1,s,v,u1,u2,u3,v1,iter)
	real*8	v(n1),u3(n1),s(n1,m),
	*u1(n1,n),u2(n1,n),v1(n1,n),stt(m,n1),
 	*u1tr(n,n1),u2tr(n,n1),v1tr(n,n1)
C оевюрэ опнлефсрнвмни яхярелш
106   FORMAT(8(D16.9,' '))
CC оевюрэ рпюейрнпхх
C	 рпюмяонмхпсел люрпхжс S - stt
       call transpon(s,n1,m,stt)
	 WRITE (50,'(A27,i2)')'люрпхжю S ХРЕПЮЖХЪ',iter
	 call PECHAT11(stt,m,n1)
	 WRITE (50,'(A20,i2)')'люрпхжю V ХРЕПЮЖХЪ',iter
	 WRITE (50, 106) v(1),v(2),v(3),
     *v(4),v(5),v(6),v(7),v(8)
CC йнмеж оевюрх рпюейрнпхх
CC оевюрэ сопюбкемхи
	 write(50,'(A20)')'сопюбкемхъ'
        call transpon(u1,n1,n,u1tr)
        call transpon(u2,n1,n,u2tr)
	  call transpon(v1,n1,n,v1tr)
	 write(50,'(A20)')'U1'
        call PECHAT11(u1tr,n,n1)
	 write(50,'(A20)')'U2'
	  call PECHAT11(u2tr,n,n1)
	 write(50,'(A20)')'V1'
	  call PECHAT11(v1tr,n,n1) 
	 WRITE (50,'(A20,i2)')'U3',iter
	  WRITE (50, 106) u3(1),u3(2),u3(3),
     *u3(4),u3(5),u3(6),u3(7),u3(8)
CC йнмеж оевюрх сопюбкемхи

	return
	end
CCCCCCC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
