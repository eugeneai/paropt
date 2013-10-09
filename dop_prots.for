      SUBROUTINE umn_matrits(A9,N9,M9,B9,L9,C9)
C	 œ–Œ÷≈ƒ”–¿ ”ÃÕŒ∆≈Õ»ﬂ Ã¿“–»÷
      real*8 A9(N9,M9),B9(M9,L9),C9(N9,L9),TEMP
      INTEGER N9,M9,L9
      INTEGER I,J,K
      DO 20 I = 1, N9
         DO 15 J = 1, L9
               TEMP = 0.0D0
            DO 10 K = 1, M9
10               TEMP = A9(I,K) * B9(K,J) + TEMP
15            C9(I,J) = TEMP
20    CONTINUE
      RETURN
      END

	subroutine transpon(a,k,k1,at)
C “–¿Õ—œŒÕ»–Œ¬¿Õ»≈ Ã¿“–»÷€ a, –≈«”À‹“¿“ Ã¿“–»÷¿ - at
	real*8 a(k,k1),at(k1,k)
	integer k,k1
	do 1 i=1,k
	do 1 j=1,k1
1	at(j,i)=a(i,j)
	return
	end
