
!! THE CODE GAMD NUMERICALLY SOLVES A (POSSIBLY STIFF)
!! SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS
!! OR A LINEARLY IMPLICIT DAE
!!
!! Copyright (C)1997-2007
!! Authors: FRANCESCA MAZZIA (mazzia@dm.uniba.it)
!!          FELICE IAVERNARO (felix@dm.uniba.it)
!!
!!
!!This program is free software; you can redistribute it and/or
!!modify it under the terms of the GNU General Public License
!!as published by the Free Software Foundation; either version 2
!!of the License, or (at your option) any later version.
!!
!!This program is distributed in the hope that it will be useful,
!!but WITHOUT ANY WARRANTY; without even the implied warranty of
!!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!GNU General Public License for more details.
!!
!!Licensed under The GNU General Public License, Version 2 or later.
!!    http://www.gnu.org/licenses/info/GPLv2orLater.html
!!
!!You should have received a copy of the GNU General Public License
!!along with this program; if not, write to the Free Software
!!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
!!USA.
!!
!!

!! THIS MODULE IS PART OF THE CODE GAMD
!! THE CODE GAMD NUMERICALLY SOLVES A (POSSIBLY STIFF)
!! SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS
!! OR A LINEARLY IMPLICIT DAE
!!
!! Copyright (C)1997-2007
!! Authors: FRANCESCA MAZZIA (mazzia@dm.uniba.it)
!!          FELICE IAVERNARO (felix@dm.uniba.it)
!!
!!
!!This program is free software; you can redistribute it and/or
!!modify it under the terms of the GNU General Public License
!!as published by the Free Software Foundation; either version 2
!!of the License, or (at your option) any later version.
!!
!!This program is distributed in the hope that it will be useful,
!!but WITHOUT ANY WARRANTY; without even the implied warranty of
!!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!GNU General Public License for more details.
!!
!!Licensed under The GNU General Public License, Version 2 or later.
!!    http://www.gnu.org/licenses/info/GPLv2orLater.html
!!
!!You should have received a copy of the GNU General Public License
!!along with this program; if not, write to the Free Software
!!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
!!USA.
!!
!!

MODULE PRECISION
      IMPLICIT NONE
      INTEGER,PARAMETER     :: PREC = KIND(1.0d0)
      REAL(PREC), PARAMETER :: EPS  = EPSILON(1.d0)
END MODULE PRECISION
MODULE LINALGGAMD
!!-----------------------------------------------------------------------
!!     ADDITIONAL LINEAR ALGEBRA ROUTINES REQUIRED BY GAMD
!!-----------------------------------------------------------------------
!!     VERSION OF MAY 16, 2003
!!-----------------------------------------------------------------------
!!
       USE PRECISION
       IMPLICIT NONE
!! PARAMETER
       INTEGER :: MLLU, MULU, MDIAG, MBDIAG, MBB, MDIFF
  CONTAINS
    !!
    !!     SUBROUTINE DEC
    !!
SUBROUTINE DEC (N, NDIM, A, IP, IER)
  !! VERSION REAL DOUBLE PRECISION
  USE PRECISION
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N,NDIM
  INTEGER, INTENT(IN OUT) :: IP(N), IER
  INTEGER :: NM1,K,KP1,M,I,J
  REAL(PREC), INTENT(IN OUT):: A(NDIM,N)
  REAL(PREC) :: T
  !!-----------------------------------------------------------------------
  !!  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
  !!  INPUT..
  !!     N = ORDER OF MATRIX.
  !!     NDIM = DECLARED DIMENSION OF ARRAY  A .
  !!     A = MATRIX TO BE TRIANGULARIZED.
  !!  OUTPUT..
  !!     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
  !!     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
  !!     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
  !!     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
  !!     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
  !!           SINGULAR AT STAGE K.
  !!  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
  !!  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
  !!  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
  !!
  !!  REFERENCE..
  !!     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
  !!     C.A.C.M. 15 (1972), P. 274.
  !!-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,N
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,N
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN

!!----------------------- END OF SUBROUTINE DEC -------------------------
END SUBROUTINE DEC
                  !!
                  !!     SUBROUTINE SOL
                  !!
SUBROUTINE SOL (N, NDIM, A, B, IP)
  !! VERSION REAL DOUBLE PRECISION
  USE PRECISION
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N,NDIM
  INTEGER, INTENT(IN) :: IP(N)
  INTEGER :: K,NM1,KP1,M,I,KB,KM1
  REAL(PREC), INTENT(IN) :: A(NDIM,N)
  REAL(PREC), INTENT(IN OUT) :: B(N)
  REAL(PREC) :: T
  !!-----------------------------------------------------------------------
  !!  SOLUTION OF LINEAR SYSTEM, A*X = B .
  !!  INPUT..
  !!    N = ORDER OF MATRIX.
  !!    NDIM = DECLARED DIMENSION OF ARRAY  A .
  !!    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
  !!    B = RIGHT HAND SIDE VECTOR.
  !!    IP = PIVOT VECTOR OBTAINED FROM DEC.
  !!  DO NOT USE IF DEC HAS SET IER .NE. 0.
  !!  OUTPUT..
  !!    B = SOLUTION VECTOR, X .
  !!-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        DO 10 I = KP1,N
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
!!----------------------- END OF SUBROUTINE SOL -------------------------
END SUBROUTINE SOL
!!
!!     SUBROUTINE DECB
!!
SUBROUTINE DECB (N, NDIM, A, ML, MU, IP, IER)
  USE PRECISION
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, NDIM, ML, MU
  INTEGER, INTENT(OUT) :: IP(N), IER
  REAL(PREC), INTENT(IN OUT) :: A(NDIM,N)
  REAL(PREC) :: T
  INTEGER :: MD, MD1,NM1,M, JU, I, J, K, KP1, MDL, MM, JK, IJK
  !!-----------------------------------------------------------------------
  !!  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
  !!  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
  !!  INPUT..
  !!     N       ORDER OF THE ORIGINAL MATRIX A.
  !!     NDIM    DECLARED DIMENSION OF ARRAY  A.
  !!     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS
  !!                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
  !!                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
  !!                ML+1 THROUGH 2*ML+MU+1 OF  A.
  !!     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
  !!     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
  !!  OUTPUT..
  !!     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
  !!                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
  !!     IP      INDEX VECTOR OF PIVOT INDICES.
  !!     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
  !!     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
  !!                SINGULAR AT STAGE K.
  !!  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
  !!  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
  !!  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO.
  !!
  !!  REFERENCE..
  !!     THIS IS A MODIFICATION OF
  !!     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
  !!     C.A.C.M. 15 (1972), P. 274.
  !!-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
  5   A(I,J) = 0.D0
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        T = A(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(MD,K)
        A(MD,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = MD1,MDL
 30       A(I,K) = -A(I,K)*T
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          T = A(M,J)
          IF (M .EQ. MM) GO TO 35
          A(M,J) = A(MM,J)
          A(MM,J) = T
 35       CONTINUE
          IF (T .EQ. 0.D0) GO TO 45
          JK = J - K
          DO 40 I = MD1,MDL
            IJK = I - JK
 40         A(IJK,J) = A(IJK,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(MD,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
  !!----------------------- END OF SUBROUTINE DECB ------------------------
END SUBROUTINE DECB
                   !!
                   !!     SUBROUTINE SOLB
                   !!
SUBROUTINE SOLB (N, NDIM, A, ML, MU, B, IP)
  USE PRECISION
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, NDIM, ML, MU
  INTEGER, INTENT(IN) :: IP(N)
  REAL(PREC), INTENT(IN) :: A(NDIM,N)
  REAL(PREC), INTENT(IN OUT) :: B(N)
  REAL(PREC) :: T
  INTEGER :: MD, MD1,NM1,M, JU, I, J, K, KB,MDM,KMD, LM, MDL, IMD

  !!-----------------------------------------------------------------------
  !!  SOLUTION OF LINEAR SYSTEM, A*X = B .
  !!  INPUT..
  !!    N      ORDER OF MATRIX A.
  !!    NDIM   DECLARED DIMENSION OF ARRAY  A .
  !!    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB.
  !!    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
  !!    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
  !!    B      RIGHT HAND SIDE VECTOR.
  !!    IP     PIVOT VECTOR OBTAINED FROM DECB.
  !!  DO NOT USE IF DECB HAS SET IER .NE. 0.
  !!  OUTPUT..
  !!    B      SOLUTION VECTOR, X .
  !!-----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
 10       B(IMD) = B(IMD) + A(I,K)*T
 20     CONTINUE
 25   CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        B(K) = B(K)/A(MD,K)
        T = -B(K)
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
 30       B(IMD) = B(IMD) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(MD,1)
      RETURN
!!----------------------- END OF SUBROUTINE SOLB ------------------------
END SUBROUTINE SOLB





FUNCTION  MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,Y)

  !!-----------------------------------------------------------------------
  !!  MATRIX VECTOR MULTIPLICATION WITH BANDED MATRIX FMAS
  !!  INPUT..
  !!    R      DECLARED DIMENSION OF ARRAY FMAS.
  !!    LDMAS  LEADING DIMENSION OF ARRAY FMAS.
  !!    MLMAS  LOWER BANDWIDTH OF FMAS (DIAGONAL IS NOT COUNTED).
  !!    MUMAS  UPPER BANDWIDTH OF FMAS (DIAGONAL IS NOT COUNTED).
  !!  OUTPUT..
  !!    MATMULB = FMAS*Y
  !!-----------------------------------------------------------------------

  USE PRECISION
  IMPLICIT NONE

  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) :: R, LDMAS, MLMAS, MUMAS
  REAL(PREC), INTENT(IN) ::   FMAS(LDMAS,R), Y(R)
  !!
  !!   OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC) :: MATMULB(R)
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: I, J
  REAL(PREC) ::  S1

  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  !!--------- ONE STEP OF THE ITERATION PROCEDURE

      MATMULB(1:R)=0.D0
      DO I=1,R
         S1=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(R,I+MUMAS)
            S1=S1+FMAS(I-J+MBDIAG,J)*Y(J)
         END DO
         MATMULB(I)=MATMULB(I)+S1
      END DO
END FUNCTION MATMULB


END MODULE LINALGGAMD
MODULE SUBGAMD
!!-----------------------------------------------------------------------
!!     ADDITIONAL ROUTINES REQUIRED BY GAMD
!!-----------------------------------------------------------------------
!!     VERSION OF AUGUST 29, 2003
!!-----------------------------------------------------------------------
!!
       USE PRECISION
       USE LINALGGAMD
       IMPLICIT NONE

!!  ORDER 3

      REAL(PREC), PARAMETER ::    &
     &            B311 = 5d0/12d0, &
     &            B312 = 8d0/12d0, &
     &            B313 = -1d0/12d0

       REAL(PREC), PARAMETER ::  Cp31 = 1d0/24d0
       REAL(PREC), PARAMETER :: &
     &            L321 =  4.012395124208693d-01,&
     &            L331 =  8.819910099032529d-04,&
     &            L341 =  1.728116022258560d-04,&
     &            L332 =  3.680857287181573d-01,&
     &            L342 =  1.635381132422046d-03,&
     &            L343 =  3.688541178419062d-01
!!--------- ORDER 5
!!--------- B3-B5 OF DIMENSION  4
       REAL(PREC), PARAMETER :: &
     &            B3511 = 49d0/720d0, &
     &            B3512 = -83d0/360d0,&
     &            B3513 = 17d0/60d0,  &
     &            B3514 = -53d0/360d0,&
     &            B3515 = 19d0/720d0, &
     &            B3521 = 19d0/720d0, &
     &            B3522 = -23d0/360d0,&
     &            B3523 = 1d0/30d0,   &
     & B3524 = 7d0/360d0,             &
     & B3525 = -11d0/720d0,           &
     & B3531 = -11d0/720d0,           &
     & B3532 = 37d0/360d0,            &
     & B3533 = -13d0/60d0,            &
     & B3534 = 67d0/360d0,            &
     & B3535 = -41d0/720d0,           &
     & B3541 = 19d0/720d0,            &
     & B3542 = -53d0/360d0,           &
     & B3543 =  17d0/60d0,            &
     & B3544 = -83d0/360d0,           &
     & B3545 =  49d0/720d0

!![ 49/720, -83/360,  17/60, -53/360,  19/720]
!![ 19/720, -23/360,   1/30,   7/360, -11/720]
!![-11/720,  37/360, -13/60,  67/360, -41/720]
!![ 19/720, -53/360,  17/60, -83/360,  49/720]

!!--------- A5, B5, B53, B56 Cp5 := MATRICES DEFINING  GAM5

       REAL(PREC), PARAMETER :: &
     & B511 = 251d0/720d0,      &
     & B512 = 323d0/360d0,      &
     & B513 = - 11d0/30d0,      &
     & B514 = 53d0/360d0,       &
     & B515 = -19d0/720d0,      &
     & B521 = -19d0/720d0,      &
     & B522 =  173d0/360d0,     &
     & B523 = 19d0/30d0,        &
     & B524 = -37d0/360d0,      &
     & B525 = 11d0/720d0

       REAL(PREC), PARAMETER:: &
     & B5711 = 1997d0/60480d0, &
     & B5712 = -113d0/630d0,   &
     & B5713 = 1619d0/4032d0,  &
     & B5714 = -715d0/1512d0,  &
     & B5715 = 1241d0/4032d0,  &
     & B5716 = -263d0/2520d0,  &
     & B5717 = 863d0/60480d0,  &
     & B5721 = -733d0/60480d0, &
     & B5722 = 41d0/630d0,     &
     & B5723 = -193d0/1344d0,  &
     & B5724 = 251d0/1512d0,   &
     & B5725 = -425d0/4032d0,  &
     & B5726 = 29d0/840d0,     &
     & B5727 = -271d0/60480d0, &
     & B5731 = -271d0/60480d0, &
     & B5732 = 97d0/5040d0,    &
     & B5733 = -13d0/448d0,    &
     & B5734 = 5d0/378d0,      &
     & B5735 = 37d0/4032d0,    &
     & B5736 = -19d0/1680d0,   &
     & B5737 = 191d0/60480d0,  &
     & B5741 = 191d0/60480d0,  &
     & B5742 = -67d0/2520d0,   &
     & B5743 = 115d0/1344d0,   &
     & B5744 = -211d0/1512d0,  &
     & B5745 = 499d0/4032d0,   &
     & B5746 = -2d0/35d0,      &
     & B5747 = 653d0/60480d0

!!--------- B5 - B7
!!
!![1997/60480,  -113/630, 1619/4032, -715/1512, 1241/4032, -263/2520,  863/60480]
!![-733/60480,    41/630, -193/1344,  251/1512, -425/4032,    29/840, -271/60480]
!![-271/60480,   97/5040,   -13/448,     5/378,   37/4032,  -19/1680,  191/60480]
!![ 191/60480,  -67/2520,  115/1344, -211/1512,  499/4032,     -2/35,  653/60480]
!![-271/60480,    29/840, -425/4032,  251/1512, -193/1344,    41/630, -733/60480]
!![ 863/60480, -263/2520, 1241/4032, -715/1512, 1619/4032,  -113/630, 1997/60480]
!!

       REAL(PREC), PARAMETER :: &
     &  L521 =  3.668340831928216D-01, &
     & L531 =  2.477905683677308D-03,  &
     & L541 = -1.919925047010838D-03,  &
     & L551 =  2.218385581234200D-03,  &
     & L561 = -5.442189351609260D-03,  &
     & L532 =  3.216639533696728D-01,  &
     & L542 =  1.231925763308414D-03,  &
     & L552 =  7.841944627374794D-03,  &
     & L562 =  1.002485104590053D-03,  &
     & L543 =  3.375100828961925D-01,  &
     & L553 = -2.614300734741796D-04,  &
     & L563 =  1.066631182323580D-03,  &
     & L554 =  3.523137378783708D-01,  &
     & L564 = -3.596681121610224D-04,  &
     & L565 =  3.617716171655064D-01,  &
     & CP51 =   3D0/160D0,             &
     & CP52 = -11D0/1440D0

!!--------- A7, B7, B57, B58 Cp7 := MATRICES DEFINING GAM7



        REAL(PREC), PARAMETER :: &
     & B711 = 19087d0/60480d0,   &
     & B712 = 2713d0/2520d0,     &
     & B713 = -15487d0/20160d0,  &
     & B714 = 586d0/945d0,       &
     & B715 = -6737d0/20160d0,   &
     & B716 = 263d0/2520d0,      &
     & B717 = -863d0/60480d0,    &
     & B721 = -863d0/60480d0,    &
     & B722 = 349d0/840d0,       &
     & B723 = 5221d0/6720d0,     &
     & B724 = -254d0/945d0,      &
     & B725 = 811d0/6720d0,      &
     & B726 = -29d0/840d0,       &
     & B727 = 271d0/60480d0,     &
     & B731 = 271d0/60480d0,     &
     & B732 = -23d0/504d0,       &
     & B733 = 10273d0/20160d0,   &
     & B734 = 586d0/945d0,       &
     & B735 = -2257d0/20160d0,   &
     & B736 = 67d0/2520d0,       &
     & B737 = -191d0/60480d0

!!
!!--------- THE LAST THREE ROWS ARE THE REVERSE OF THE FIRST THREE
!!
!! B79 =
!!
!![ 75203/3628800, -280187/1814400, 129781/259200, -238937/259200, 27289/25920, -197687/259200,  88531/259200, -156437/1814400,  33953/3628800]
!![-17827/3628800,   66043/1814400, -30389/259200,   55513/259200, -6281/25920,   44983/259200, -19859/259200,   34453/1814400,  -7297/3628800]
!![  8963/3628800,  -32987/1814400,  15061/259200,  -27257/259200,  3049/25920,  -21527/259200,   9331/259200,  -15797/1814400,   3233/3628800]
!![  3233/3628800,  -10067/1814400,   3601/259200,   -4337/259200,     23/3240,    1393/259200,  -2129/259200,    7123/1814400,  -2497/3628800]
!![ -2497/3628800,   12853/1814400,  -7859/259200,   18583/259200, -2681/25920,   24313/259200, -13589/259200,   30043/1814400,  -8227/3628800]
!![  3233/3628800,  -15797/1814400,   9331/259200,  -21527/259200,  3049/25920,  -27257/259200,  15061/259200,  -32987/1814400,   8963/3628800]
!![ -7297/3628800,   34453/1814400, -19859/259200,   44983/259200, -6281/25920,   55513/259200, -30389/259200,   66043/1814400, -17827/3628800]
!![ 33953/3628800, -156437/1814400,  88531/259200, -197687/259200, 27289/25920, -238937/259200, 129781/259200, -280187/1814400,  75203/3628800]
!!
       REAL(PREC), PARAMETER ::     &
     & B7911 =   75203d0/3628800d0, &
     & B7912 = -280187d0/1814400d0, &
     & B7913 =  129781d0/259200d0,  &
     & B7914 = -238937d0/259200d0,  &
     & B7915 =   27289d0/25920d0,   &
     & B7916 = -197687d0/259200d0,  &
     & B7917 =   88531d0/259200d0,  &
     & B7918 = -156437d0/1814400d0, &
     & B7919 =   33953d0/3628800d0

       REAL(PREC), PARAMETER ::    &
     & B7921 = -17827d0/3628800d0, &
     & B7922 =  66043d0/1814400d0, &
     & B7923 = -30389d0/259200d0,  &
     & B7924 =  55513d0/259200d0,  &
     & B7925 =  -6281d0/25920d0,   &
     & B7926 =  44983d0/259200d0,  &
     & B7927 = -19859d0/259200d0,  &
     & B7928 =  34453d0/1814400d0, &
     & B7929 =  -7297d0/3628800d0

       REAL(PREC), PARAMETER :: &
     & B7931 =   8963d0/3628800d0, &
     & B7932 = -32987d0/1814400d0, &
     & B7933 =  15061d0/259200d0,  &
     & B7934 = -27257d0/259200d0,  &
     & B7935 =   3049d0/25920d0,   &
     & B7936 = -21527d0/259200d0,  &
     & B7937 =   9331d0/259200d0,  &
     & B7938 = -15797d0/1814400d0, &
     & B7939 =   3233d0/3628800d0

       REAL(PREC), PARAMETER :: &
     & B7941 =   3233d0/3628800d0, &
     & B7942 = -10067d0/1814400d0, &
     & B7943 =   3601d0/259200d0,  &
     & B7944 =  -4337d0/259200d0,  &
     & B7945 =     23d0/3240d0,    &
     & B7946 =   1393d0/259200d0,  &
     & B7947 =  -2129d0/259200d0,  &
     & B7948 =   7123d0/1814400d0, &
     & B7949 =  -2497d0/3628800d0


      REAL(PREC),  PARAMETER :: &
     & B7951 =  -2497d0/3628800d0,&
     & B7952 =  12853d0/1814400d0,&
     & B7953 =  -7859d0/259200d0, &
     & B7954 =  18583d0/259200d0, &
     & B7955 =  -2681d0/25920d0,  &
     & B7956 =  24313d0/259200d0, &
     & B7957 = -13589d0/259200d0, &
     & B7958 =  30043d0/1814400d0,&
     & B7959 =  -8227d0/3628800d0

      REAL(PREC),  PARAMETER :: &
     & L721 = 3.023839891568610D-01, &
     & L731 = 3.201698610574002D-05, &
     & L741 = 4.193101163680004D-04, &
     & L751 = 1.686924996069667D-04, &
     & L761 = 4.806043527549464D-05, &
     & L771 = 3.598347048026785D-06, &
     & L781 = 7.892534649789167D-04, &
     & L732 = 2.559868364091398D-01, &
     & L742 = 1.336896192287030D-04, &
     & L752 = 3.080994719931695D-03, &
     & L762 = 1.457177183563680D-04, &
     & L772 = 9.259360509484074D-04, &
     & L782 = 2.397658879381223D-04, &
     & L743 = 2.639734712170458D-01, &
     & L753 = 1.734338929611258D-04, &
     & L763 = 6.704398263264620D-03, &
     & L773 = 4.559927214651730D-05, &
     & L783 = 6.396418554053151D-05, &
     & L754 = 2.817729090368562D-01, &
     & L764 = 2.877761776030408D-04, &
     & L774 = 1.810919475521773D-04, &
     & L784 = 1.009049833235848D-03, &
     & L765 = 2.993040718034231D-01, &
     & L775 = 2.009850887505898D-03, &
     & L785 = 1.748065618845750D-03, &
     & L776 = 3.150349043479135D-01, &
     & L786 = 3.243816792609449D-05, &
     & L787 = 3.271307059448932D-01

      REAL(PREC),  PARAMETER :: &
     & CP71 = 103D0/9061D0,     &
     & CP72 = -13D0/4480D0,     &
     & CP73 =  67D0/42431D0

!!--------- A8, B8, B86, B810 Cp8 := MATRICES DEFINING GAM9

       REAL(PREC), PARAMETER ::   &
     & B911 = 1070017D0/3628800D0,&
     & B912 = 2233547D0/1814400D0,&
     & B913 = -2302297D0/1814400D0,&
     & B914 = 2797679D0/1814400D0, &
     & B915 = -31457D0/22680D0,    &
     & B916 = 1573169D0/1814400D0, &
     & B917 = -645607D0/1814400D0, &
     & B918 = 156437D0/1814400D0,  &
     & B919 = -33953D0/3628800D0,  &
     & B921 = -33953D0/3628800D0,  &
     & B922 = 687797D0/1814400D0,  &
     & B923 =  1622393D0/1814400D0,&
     & B924 = -876271D0/1814400D0, &
     & B925 =   8233D0/22680D0,    &
     & B926 =    -377521D0/1814400D0, &
     & B927 =   147143D0/1814400D0,   &
     & B928 =  -34453D0/1814400D0,    &
     & B929 =   7297D0/3628800D0,     &
     & B931 = 7297D0/3628800D0,       &
     & B932 =  -49813D0/1814400D0,    &
     & B933 =  819143D0/1814400D0,    &
     & B934 =  1315919D0/1814400D0,   &
     & B935 = -5207D0/22680D0,        &
     & B936 =  198929D0/1814400D0,    &
     & B937 =  -71047D0/1814400D0,    &
     & B938 =  15797D0/1814400D0,     &
     & B939 = -3233D0/3628800D0

       REAL(PREC), PARAMETER ::  &
     & B941 = -3233D0/3628800D0,      &
     & B942 = 18197D0/1814400D0,      &
     & B943 =  -108007D0/1814400D0,   &
     & B944 =  954929D0/1814400D0,    &
     & B945 = 13903D0/22680D0,        &
     & B946 = -212881D0/1814400D0,    &
     & B947 =  63143D0/1814400D0,     &
     & B948 =  -12853D0/1814400D0,    &
     & B949 =  2497D0/3628800D0


!!   B910=B9-B10;
!!               prime 4 righe : la 5 e' uguale alla 4;
!!     le ultime 4 uguali alle prime 4 negate (simmetriche);
!!    le prime  5 colonne  : le altre sono simmetriche e negate;
!!
!![ 8183/1036800, -8183/115200,   8183/28800, -57281/86400,  57281/57600, -57281/57600,  57281/86400,  -8183/28800,  8183/115200, -8183/1036800]
!![  -425/290304,    425/32256,    -425/8064,     425/3456,    -425/2304,     425/2304,    -425/3456,     425/8064,   -425/32256,    425/290304]
!![      7/12800,    -63/12800,      63/3200,    -147/3200,     441/6400,    -441/6400,     147/3200,     -63/3200,     63/12800,      -7/12800]
!![-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600]
!![-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600]
!![ 2497/7257600, -2497/806400,  2497/201600,  -2497/86400,   2497/57600,  -2497/57600,   2497/86400, -2497/201600,  2497/806400, -2497/7257600]
!![     -7/12800,     63/12800,     -63/3200,     147/3200,    -441/6400,     441/6400,    -147/3200,      63/3200,    -63/12800,       7/12800]
!![   425/290304,   -425/32256,     425/8064,    -425/3456,     425/2304,    -425/2304,     425/3456,    -425/8064,    425/32256,   -425/290304]
!![-8183/1036800,  8183/115200,  -8183/28800,  57281/86400, -57281/57600,  57281/57600, -57281/86400,   8183/28800, -8183/115200,  8183/1036800]

       REAL(PREC), PARAMETER :: &
     & B91011 =   8183d0/1036800d0,  &
     & B91012 =  -8183d0/115200d0,  &
     & B91013 =   8183d0/28800d0,   &
     & B91014 = -57281d0/86400d0,   &
     & B91015 =  57281d0/57600d0

       REAL(PREC), PARAMETER ::   &
     & B91021 =   -425d0/290304d0,&
     & B91022 =    425d0/32256d0, &
     & B91023 =   -425d0/8064d0,  &
     & B91024 =    425d0/3456d0,  &
     & B91025 =   -425d0/2304d0

       REAL(PREC), PARAMETER :: &
     & B91031 =     7d0/12800d0,&
     & B91032 =   -63d0/12800d0,&
     & B91033 =    63d0/3200d0, &
     & B91034 =  -147d0/3200d0, &
     & B91035 =   441d0/6400d0

       REAL(PREC), PARAMETER ::    &
     & B91041 = -2497d0/7257600d0, &
     & B91042 =  2497d0/806400d0,  &
     & B91043 = -2497d0/201600d0,  &
     & B91044 =  2497d0/86400d0,   &
     & B91045 = -2497d0/57600d0

      REAL(PREC),  PARAMETER :: &
     & Cp91 =  7.892554012345216d-03, &
     & Cp92 = -1.463982583774219d-03, &
     & Cp93 =  5.468749999999983d-04, &
     & Cp94 = -3.440531305114634d-04

!!--------- THE OTHERS ARE THE SAME WITH CHANGED SIGN


       REAL(PREC), PARAMETER :: &
     & L921   =   2.590721934790442d-01,     &
     & L932   =   2.077575545359853d-01,     &
     & L943   =   2.032874698558627d-01,     &
     & L954   =   2.036384888660128d-01,     &
     & L965   =   2.039599505779785d-01,     &
     & L976   =   2.034044409161703d-01,     &
     & L987   =   2.017245408702437d-01,     &
     & L998   =   1.986549276295617d-01


CONTAINS

  SUBROUTINE DECLU(R,JF0,H,LDJAC,LU,LDLU,IPIV,FMAS,LDMAS,MLMAS,MUMAS,ORD,IER,IJOB)
    USE PRECISION
    USE  LINALGGAMD
    IMPLICIT NONE
    !!
    !!   INPUT VARIABLES
    !!------------------------------------
    INTEGER, INTENT(IN):: R, LDJAC, LDLU, LDMAS, MLMAS, MUMAS, ORD, IJOB
    REAL(PREC),INTENT(IN) ::  JF0(LDJAC,R), FMAS(LDMAS,R),  H
    !!
    !!   OUTPUT VARIABLES
    !!------------------------------------
    INTEGER, INTENT(OUT):: IER, IPIV(R)
    REAL(PREC),INTENT(OUT) :: LU(LDLU,R)
    !!
    !!   LOCAL VARIABLES
    !!------------------------------------
    INTEGER :: I, J, IB
    REAL(PREC) :: FAC


    REAL(PREC), PARAMETER :: L31  =  6.411501944628007d-01,  &
                           & L51  =  6.743555662880509D-01,  &
                           & L71  =  7.109158294404152D-01,  &
                           & L91  =  7.440547954061898d-01


    !!
    !!   EXECUTABLE STATEMENTS
    !!---------------------------------
    !!
    SELECT CASE(ORD)

    CASE(1)
       FAC = -L31*H
    CASE(2)
       FAC = -L51*H
    CASE(3)
       FAC = -L71*H
    CASE(4)
       FAC = -L91*H
    END SELECT

    SELECT CASE(IJOB)
    CASE(1)
       !! -------- ODE: JACOBIAN A FULL MATRIX
       DO J=1,R
          LU(1:R,J)= FAC*JF0(1:R,J)
          LU(J,J)=LU(J,J)+1d0
       END DO
       CALL DEC (R,LDLU,LU,IPIV,IER)
       RETURN
    CASE(2)
       !! -------- ODE: JACOBIAN A BAND MATRIX
       DO J=1,R
          LU(MLLU+1:MLLU+MDIAG,J)= FAC*JF0(1:MDIAG,J)
          LU(MDIAG,J)=LU(MDIAG,J)+1d0
       END DO
       CALL DECB (R,LDLU,LU,MLLU,MULU,IPIV,IER)
    CASE(3)
       !! -------- DAE: FMAS A BAND MATRIX, JACOBIAN A FULL MATRIX
       DO J=1,R
          LU(1:R,J)= FAC*JF0(1:R,J)
          DO I=MAX(1,J-MUMAS),MIN(R,J+MLMAS)
             LU(I,J)=LU(I,J)+FMAS(I-J+MBDIAG,J)
          END DO
       END DO
       CALL DEC (R,LDLU,LU,IPIV,IER)
    CASE(4)
       !! -------- DAE: FMAS A BAND MATRIX, JACOBIAN A BAND MATRIX
       DO J=1,R
          LU(MLLU+1:MLLU+MDIAG,J)= FAC*JF0(1:MDIAG,J)
          DO I=1,MBB
            IB=I+MDIFF
            LU(IB,J)=LU(IB,J)+FMAS(I,J)
          END DO
       END DO
       CALL DECB (R,LDLU,LU,MLLU,MULU,IPIV,IER)
    CASE(5)
       !! -------- DAE: FMAS A FULL MATRIX, JACOBIAN A FULL MATRIX
          LU(1:R,1:R)= FMAS(1:R,1:R) + FAC*JF0(1:R,1:R)
       CALL DEC (R,LDLU,LU,IPIV,IER)
    END SELECT
    RETURN

  END SUBROUTINE DECLU

!!
!!  SUBROUTINE SOLLU
!!

SUBROUTINE SOLLU(R,LU,LDLU,F,IPIV,IJOB)
  USE PRECISION
  USE LINALGGAMD
  IMPLICIT NONE

  !!
  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN):: R, LDLU, IPIV(R), IJOB
  REAL(PREC), INTENT(IN) ::  LU(LDLU,R)
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) :: F(R)
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  SELECT CASE(IJOB)
  CASE(1,3,5)
     !! -------- JACOBIAN A FULL MATRIX
     CALL SOL (R,LDLU,LU(1,1),F(1),IPIV(1))
  CASE(2,4)
     !! -------- JACOBIAN A BAND MATRIX
     CALL SOLB (R,LDLU,LU(1,1),MLLU,MULU,F(1),IPIV(1))

  END SELECT

  RETURN

END SUBROUTINE SOLLU


!!
!!  SUBROUTINE NEWTGS
!!

SUBROUTINE NEWTGS(R,DBLK,LU,LDLU,FMAS,LDMAS,MLMAS,MUMAS,IPIV,F,DN,IJOB)
  USE PRECISION
  IMPLICIT NONE
  !!
  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN):: R, LDLU, LDMAS, MLMAS, MUMAS, IPIV(R), IJOB, DBLK
  REAL(PREC), INTENT(IN)::  F(R,DBLK), LU(LDLU,R), FMAS(LDMAS,R)
  !!
  !!   OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(OUT)::  DN(R,DBLK)
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER     ::  J
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!

  SELECT CASE(IJOB)
  CASE(1,2)  !! ODE

     DN(1:R,1) = -F(1:R,1)
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1,1),IPIV(1),IJOB)
     DO J=2,DBLK
        DN(1:R,J) =  -F(1:R,J)+DN(1:R,J-1)
        CALL  SOLLU(R,LU(1,1),LDLU,DN(1,J),IPIV(1),IJOB)
     END DO

  CASE(3,4)  !! DAE: FMAS A BAND MATRIX

     DN(1:R,1) = -F(1:R,1)
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1,1),IPIV(1),IJOB)
     DO J=2,DBLK
        DN(1:R,J) =  -F(1:R,J) + MATMULB(R,FMAS(1,1),LDMAS,MLMAS,MUMAS,DN(1:R,J-1))
        CALL  SOLLU(R,LU(1,1),LDLU,DN(1,J),IPIV(1),IJOB)
     END DO

  CASE(5)    !! DAE: FMAS A FULL MATRIX

     DN(1:R,1) = -F(1:R,1)
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1,1),IPIV(1),IJOB)
     DO J=2,DBLK
        DN(1:R,J) =  -F(1:R,J) + MATMUL(FMAS(1:R,1:R),DN(1:R,J-1))
        CALL  SOLLU(R,LU(1,1),LDLU,DN(1,J),IPIV(1),IJOB)
     END DO



  END SELECT

  RETURN

END SUBROUTINE NEWTGS

!!
!!   FINE SUBROUTINE NEWTGS
!!
SUBROUTINE INTERP(R,TP,YP,T1,F1,NT1,DBLKOLD,DBLK,T0,Y0)
  USE PRECISION
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: R, DBLK, DBLKOLD, NT1
  REAL(PREC), INTENT(IN) :: F1(R,*), T1(*), Y0(R), T0
  REAL(PREC), INTENT(IN OUT) :: YP(R,*), TP(*)
  INTEGER :: I,J, N, IT1, NT2, NTi

  IF (DBLK < DBLKOLD) THEN
     NTi = MAX(5,NT1)
  ELSE
     NTi = MAX(3,NT1)
  END IF
  NT2 = NTi+1
  N = DBLKOLD+1


  DO IT1=2,DBLK+1
     YP(1:R,IT1) = F1(1:R,NTi)
     DO J=NT2,N
        YP(1:R,IT1) = YP(1:R,IT1)*(T1(IT1-1)-TP(J)) + F1(1:R,J)
     END DO
  END DO



  YP(1:R,1) = Y0(1:R)

  TP(1) = T0
  TP(2:DBLK+1) = T1(1:DBLK)

  RETURN
END SUBROUTINE INTERP

!!
!!    SUBROUTINE DIFFDIV
!!

SUBROUTINE DIFFDIV(TP,YP,R,DBLK,NT1)
  USE PRECISION
  IMPLICIT NONE
  !!
  !!   INPUT VARIABLES
  !!-----------------------------------
  INTEGER, INTENT(IN) :: R,  DBLK
  !!
  !!   INPUT/OUTPUT VARIABLES
  !------------------------------------
  INTEGER, INTENT(IN OUT):: NT1
  REAL(PREC), INTENT(IN OUT):: TP(DBLK+1), YP(R,DBLK+1)
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  ::  J, K, N
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!

  N = DBLK+1
  NT1 = 1

  DO J=N-1,NT1,-1
     DO K=1,J
        YP(1:R,K)= ( YP(1:R,K)- YP(1:R,K+1) )/( TP(K)-TP(K+N-J))
     END DO
  END DO

  RETURN
END SUBROUTINE DIFFDIV
!!
!!     FUNCTION  CONTR
!!
FUNCTION CONTR(I,R,T,TP,FF,DBLK,NT1)
  !! ----------------------------------------------------------
  !!     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT. IT PROVIDES AN
  !!     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT T.
  !!     IT GIVES THE VALUE OF THE INTERPOLATION POLYNOMIAL, DEFINED FOR
  !!     THE LAST SUCCESSFULLY COMPUTED STEP.
  !! ----------------------------------------------------------
  USE PRECISION
  IMPLICIT NONE
  !!
  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) ::  R, I, DBLK, NT1
  REAL(PREC) :: CONTR
  REAL(PREC), INTENT(IN) ::  T, TP(DBLK+1), FF(R,DBLK+1)
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER :: J, N, NT2, NTi
  REAL(PREC) :: YP
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  N = DBLK+1
  NTi = MAX(1,NT1)
  NT2=NTi+1
  YP = FF(I,NTi)
  DO J=NT2,N
     YP = YP*(T-TP(J)) + FF(I,J)
  END DO
  CONTR = YP
  RETURN
END FUNCTION  CONTR


    !!
    !!  SUBROUTINE TERMNOT3  (ORDER 3)
    !!
SUBROUTINE  TERMNOT3(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,                    &
                    &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV,            &
                    &  FMAS,LDMAS,MLMAS,MUMAS, SCAL,IJOB,TER,RPAR,IPAR)

  USE PRECISION
  IMPLICIT NONE

  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) :: R,  IJOB, IPIV(R), LDLU, LDMAS, MLMAS, MUMAS, IT, IPAR(1)
  REAL(PREC), INTENT(IN) ::  H, SCAL(R), LU(LDLU,R), FMAS(LDMAS,R), ERRNEWT0,TETAK0,RPAR(1)
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) :: ERRNEWT, YP(R,5), FP(R,5), &
       &                   DN(R), F1(R,5),TP(5)
  INTEGER, INTENT(IN OUT) :: NFCN
  LOGICAL, INTENT(OUT):: TER
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: J, IERR
  REAL(PREC) ::  ERRVJ, SUM, FAC, ZP(R)

  EXTERNAL FCN
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  !!--------- ONE STEP OF THE ITERATION PROCEDURE
  TER = .FALSE.


  SELECT CASE(IJOB)
  CASE(1,2)  !! ODE


      DN(1:R)=(YP(1:R,2)-YP(1:R,1))- &
       &   H*(B311*FP(1:R,1)+B312*FP(1:R,2)+B313*FP(1:R,3))

      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,2)=YP(J,2)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = sqrt(ERRVJ/R)

      ERRNEWT = MAX( ERRNEWT, ERRVJ )

      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR) ! removed ierr

      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
         TER = .TRUE.
         RETURN
      END IF


      DO J=1,R
         SUM = L321*(F1(J,1)-FP(J,2))+B311*FP(J,2) &
              &               +B312*FP(J,3)+ B313*FP(J,4)
         DN(J)=(YP(J,3)-YP(J,2))-H*SUM
      END DO
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,3)=YP(J,3)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = sqrt(ERRVJ/R)

      ERRNEWT = MAX( ERRNEWT, ERRVJ )

      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR) ! removed ierr
      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
         TER = .TRUE.
         RETURN
      END IF

      DO J=1,R
         SUM = L331*(F1(J,1)-FP(J,2))+L332*(F1(J,2)-FP(J,3)) &
              & +B311*FP(J,3)+B312*FP(J,4)+B313*FP(J,5)
         DN(J)=(YP(J,4)-YP(J,3))-H*SUM
      END DO
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,4)=YP(J,4)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = sqrt(ERRVJ/R)

      ERRNEWT = MAX( ERRNEWT, ERRVJ )


      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(4),YP(1,4),F1(1,3),  RPAR,IPAR) ! removed ierr
      IF (IERR .NE.0) THEN
         TER = .TRUE.
         RETURN
      END IF
      NFCN = NFCN + 1


      DO J=1,R
         SUM = L341*(F1(J,1)-FP(J,2))+L342*(F1(J,2)-FP(J,3)) &
              & +L343*(F1(J,3)-FP(J,4))  &
              & +B313*FP(J,3)+B312*FP(J,4)+B311*FP(J,5)
         DN(J)=(YP(J,5)-YP(J,4))-H*SUM
      END DO
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,5)=YP(J,5)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO


      ERRVJ = sqrt(ERRVJ/R)
      ERRNEWT = MAX( ERRNEWT, ERRVJ )



      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(5),YP(1,5),F1(1,4),  RPAR,IPAR) ! removed ierr
      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
          TER = .TRUE.
         RETURN
      END IF


  CASE(3,4)  !! DAE: FMAS A BAND MATRIX

      DN(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,2)-YP(1:R,1))- &
       &   H*(B311*FP(1:R,1)+B312*FP(1:R,2)+B313*FP(1:R,3))

      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,2)=YP(J,2)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = sqrt(ERRVJ/R)

      ERRNEWT = MAX( ERRNEWT, ERRVJ )
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR) ! removed ierr
      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
         TER = .TRUE.
         RETURN
      END IF
      ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,3)-YP(1:R,2))
      DO J=1,R
         SUM = L321*(F1(J,1)-FP(J,2))+B311*FP(J,2) &
              &               +B312*FP(J,3)+ B313*FP(J,4)
         DN(J)=ZP(J)-H*SUM
      END DO
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,3)=YP(J,3)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = sqrt(ERRVJ/R)

      ERRNEWT = MAX( ERRNEWT, ERRVJ )

      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
         TER = .TRUE.
         RETURN
      END IF

      ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,4)-YP(1:R,3))
      DO J=1,R
         SUM = L331*(F1(J,1)-FP(J,2))+L332*(F1(J,2)-FP(J,3)) &
              & +B311*FP(J,3)+B312*FP(J,4)+B313*FP(J,5)
         DN(J)=ZP(J)-H*SUM
      END DO
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,4)=YP(J,4)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = sqrt(ERRVJ/R)

      ERRNEWT = MAX( ERRNEWT, ERRVJ )


      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR) ! removed ierr
      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
         TER = .TRUE.
         RETURN
      END IF



      ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,5)-YP(1:R,4))
      DO J=1,R
         SUM = L341*(F1(J,1)-FP(J,2))+L342*(F1(J,2)-FP(J,3)) &
              & +L343*(F1(J,3)-FP(J,4))  &
              & +B313*FP(J,3)+B312*FP(J,4)+B311*FP(J,5)
         DN(J)=ZP(J)-H*SUM
      END DO
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,5)=YP(J,5)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO


      ERRVJ = sqrt(ERRVJ/R)
      ERRNEWT = MAX( ERRNEWT, ERRVJ )



      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR) ! removed ierr
      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
          TER = .TRUE.
         RETURN
      END IF


  CASE(5)    !! DAE: FMAS A FULL MATRIX

      DN(1:R)=MATMUL(FMAS(1:R,1:R),(YP(1:R,2)-YP(1:R,1)))- &
       &   H*(B311*FP(1:R,1)+B312*FP(1:R,2)+B313*FP(1:R,3))

      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,2)=YP(J,2)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = sqrt(ERRVJ/R)

      ERRNEWT = MAX( ERRNEWT, ERRVJ )
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)  ! removed ierr
      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
         TER = .TRUE.
         RETURN
      END IF

      ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,3)-YP(1:R,2))
      DO J=1,R
         SUM = L321*(F1(J,1)-FP(J,2))+B311*FP(J,2) &
              &               +B312*FP(J,3)+ B313*FP(J,4)
         DN(J)=ZP(J)-H*SUM
      END DO
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,3)=YP(J,3)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = sqrt(ERRVJ/R)

      ERRNEWT = MAX( ERRNEWT, ERRVJ )

      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
         TER = .TRUE.
         RETURN
      END IF

      ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,4)-YP(1:R,3))
      DO J=1,R
         SUM = L331*(F1(J,1)-FP(J,2))+L332*(F1(J,2)-FP(J,3)) &
              & +B311*FP(J,3)+B312*FP(J,4)+B313*FP(J,5)
         DN(J)=ZP(J)-H*SUM
      END DO
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,4)=YP(J,4)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = sqrt(ERRVJ/R)

      ERRNEWT = MAX( ERRNEWT, ERRVJ )


      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR) ! removed ierr
      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
         TER = .TRUE.
         RETURN
      END IF



      ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,5)-YP(1:R,4))
      DO J=1,R
         SUM = L341*(F1(J,1)-FP(J,2))+L342*(F1(J,2)-FP(J,3)) &
              & +L343*(F1(J,3)-FP(J,4))  &
              & +B313*FP(J,3)+B312*FP(J,4)+B311*FP(J,5)
         DN(J)=ZP(J)-H*SUM
      END DO
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,5)=YP(J,5)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO


      ERRVJ = sqrt(ERRVJ/R)
      ERRNEWT = MAX( ERRNEWT, ERRVJ )



      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
         TER = .TRUE.
         RETURN
      END IF

      IERR = 0
      CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR) ! removed ierr
      NFCN = NFCN + 1
      IF (IERR .NE.0) THEN
          TER = .TRUE.
         RETURN
      END IF

  END SELECT

  FP(1:R,2:5) = F1(1:R,1:4)

  RETURN
END SUBROUTINE TERMNOT3
    !!
    !!  SUBROUTINE TERMNOT5  (ORDER 5)
    !!
SUBROUTINE  TERMNOT5(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,                    &
                    &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV,            &
                    &  FMAS,LDMAS,MLMAS,MUMAS, SCAL,IJOB,TER,RPAR,IPAR)

  USE PRECISION
  IMPLICIT NONE

  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) :: R,  IJOB, IPIV(R), LDLU, LDMAS, MLMAS, MUMAS, IT, IPAR(1)
  REAL(PREC), INTENT(IN) ::  H, SCAL(R), LU(LDLU,R), FMAS(LDMAS,R), ERRNEWT0, TETAK0, RPAR(1)
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) :: ERRNEWT, YP(R,7), FP(R,7), &
       &                   DN(R), F1(R,7),TP(7)
  INTEGER, INTENT(IN OUT) :: NFCN
  LOGICAL, INTENT(OUT):: TER
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: J, IERR
  REAL(PREC) ::  ERRVJ, SUM, FAC, ZP(R)

  EXTERNAL FCN
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  !!--------- ONE STEP OF THE ITERATION PROCEDURE



  TER = .FALSE.
  SELECT CASE(IJOB)
  CASE(1,2)  !! ODE

     DO J=1,R
        SUM = B511*FP(J,1)+B512*FP(J,2)+B513*FP(J,3)+B514*FP(J,4)&
             &    +B515*FP(J,5)
        DN(J)=YP(J,2)-YP(J,1)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L521*(F1(J,1)-FP(J,2))+B521*FP(J,1)+B522*FP(J,2)&
             &               +B523*FP(J,3)+B524*FP(J,4)+B525*FP(J,5)
        DN(J)=YP(J,3)-YP(J,2)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )



     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF


     DO J=1,R
        SUM = L531*(F1(J,1)-FP(J,2))+L532*(F1(J,2)-FP(J,3))&
             &+B521*FP(J,2)+B522*FP(J,3)+B523*FP(J,4)+B524*FP(J,5)+B525*FP(J,6)
        DN(J)=YP(J,4)-YP(J,3)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L541*(F1(J,1)-FP(J,2))+L542*(F1(J,2)-FP(J,3))&
             &+L543*(F1(J,3)-FP(J,4))+B521*FP(J,3)+B522*FP(J,4)&
             &+B523*FP(J,5)+B524*FP(J,6)+B525*FP(J,7)
        DN(J)=YP(J,5)-YP(J,4)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )


     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM =  L551*(F1(J,1)-FP(J,2))+L552*(F1(J,2)-FP(J,3))&
             &+L553*(F1(J,3)-FP(J,4))+L554*(F1(J,4)-FP(J,5))&
             &+B525*FP(J,3)+B524*FP(J,4)+B523*FP(J,5)+B522*FP(J,6)+B521*FP(J,7)
        DN(J)=YP(J,6)-YP(J,5)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )


     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM =  L561*(F1(J,1)-FP(J,2))+L562*(F1(J,2)-FP(J,3))&
             &+L563*(F1(J,3)-FP(J,4))+L564*(F1(J,4)-FP(J,5))&
             &+L565*(F1(J,5)-FP(J,6))&
             &+B515*FP(J,3)+B514*FP(J,4)+B513*FP(J,5)+B512*FP(J,6)+B511*FP(J,7)
        DN(J)=YP(J,7)-YP(J,6)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

  CASE(3,4)  !! DAE: FMAS A BAND MATRIX

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,2)-YP(1:R,1))
     DO J=1,R
        SUM = B511*FP(J,1)+B512*FP(J,2)+B513*FP(J,3)+B514*FP(J,4)&
             &    +B515*FP(J,5)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,3)-YP(1:R,2))
     DO J=1,R
        SUM = L521*(F1(J,1)-FP(J,2))+B521*FP(J,1)+B522*FP(J,2)&
             &               +B523*FP(J,3)+B524*FP(J,4)+B525*FP(J,5)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )



     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,4)-YP(1:R,3))
     DO J=1,R
        SUM = L531*(F1(J,1)-FP(J,2))+L532*(F1(J,2)-FP(J,3))&
             &+B521*FP(J,2)+B522*FP(J,3)+B523*FP(J,4)+B524*FP(J,5)+B525*FP(J,6)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,5)-YP(1:R,4))
     DO J=1,R
        SUM = L541*(F1(J,1)-FP(J,2))+L542*(F1(J,2)-FP(J,3))&
             &+L543*(F1(J,3)-FP(J,4))+B521*FP(J,3)+B522*FP(J,4)&
             &+B523*FP(J,5)+B524*FP(J,6)+B525*FP(J,7)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )


     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,6)-YP(1:R,5))
     DO J=1,R
        SUM =  L551*(F1(J,1)-FP(J,2))+L552*(F1(J,2)-FP(J,3))&
             &+L553*(F1(J,3)-FP(J,4))+L554*(F1(J,4)-FP(J,5))&
             &+B525*FP(J,3)+B524*FP(J,4)+B523*FP(J,5)+B522*FP(J,6)+B521*FP(J,7)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )


     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,7)-YP(1:R,6))
     DO J=1,R
        SUM =  L561*(F1(J,1)-FP(J,2))+L562*(F1(J,2)-FP(J,3))&
             &+L563*(F1(J,3)-FP(J,4))+L564*(F1(J,4)-FP(J,5))&
             &+L565*(F1(J,5)-FP(J,6))&
             &+B515*FP(J,3)+B514*FP(J,4)+B513*FP(J,5)+B512*FP(J,6)+B511*FP(J,7)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF




  CASE(5)    !! DAE: FMAS A FULL MATRIX

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,2)-YP(1:R,1))
     DO J=1,R
        SUM = B511*FP(J,1)+B512*FP(J,2)+B513*FP(J,3)+B514*FP(J,4)&
             &    +B515*FP(J,5)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,3)-YP(1:R,2))
     DO J=1,R
        SUM = L521*(F1(J,1)-FP(J,2))+B521*FP(J,1)+B522*FP(J,2)&
             &               +B523*FP(J,3)+B524*FP(J,4)+B525*FP(J,5)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )



     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,4)-YP(1:R,3))
     DO J=1,R
        SUM = L531*(F1(J,1)-FP(J,2))+L532*(F1(J,2)-FP(J,3))&
             &+B521*FP(J,2)+B522*FP(J,3)+B523*FP(J,4)+B524*FP(J,5)+B525*FP(J,6)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,5)-YP(1:R,4))
     DO J=1,R
        SUM = L541*(F1(J,1)-FP(J,2))+L542*(F1(J,2)-FP(J,3))&
             &+L543*(F1(J,3)-FP(J,4))+B521*FP(J,3)+B522*FP(J,4)&
             &+B523*FP(J,5)+B524*FP(J,6)+B525*FP(J,7)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )


     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,6)-YP(1:R,5))
     DO J=1,R
        SUM =  L551*(F1(J,1)-FP(J,2))+L552*(F1(J,2)-FP(J,3))&
             &+L553*(F1(J,3)-FP(J,4))+L554*(F1(J,4)-FP(J,5))&
             &+B525*FP(J,3)+B524*FP(J,4)+B523*FP(J,5)+B522*FP(J,6)+B521*FP(J,7)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )


     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,7)-YP(1:R,6))
     DO J=1,R
        SUM =  L561*(F1(J,1)-FP(J,2))+L562*(F1(J,2)-FP(J,3))&
             &+L563*(F1(J,3)-FP(J,4))+L564*(F1(J,4)-FP(J,5))&
             &+L565*(F1(J,5)-FP(J,6))&
             &+B515*FP(J,3)+B514*FP(J,4)+B513*FP(J,5)+B512*FP(J,6)+B511*FP(J,7)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF


END SELECT

  FP(1:R,2:7) = F1(1:R,1:6)

  RETURN
END SUBROUTINE TERMNOT5

    !!
    !!  SUBROUTINE TERMNOT7  (ORDER 7)
    !!
SUBROUTINE  TERMNOT7(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,                    &
                    &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV,            &
                    &  FMAS,LDMAS,MLMAS,MUMAS, SCAL,IJOB,TER,RPAR,IPAR)

  USE PRECISION
  IMPLICIT NONE

  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) :: R,  IJOB, IPIV(R), LDLU, LDMAS, MLMAS, MUMAS, IT, IPAR(1)
  REAL(PREC), INTENT(IN) ::  H, SCAL(R), LU(LDLU,R), FMAS(LDMAS,R), ERRNEWT0, TETAK0, RPAR(1)
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) :: ERRNEWT, YP(R,9), FP(R,9), &
       &                   DN(R), F1(R,9),TP(9)
  INTEGER, INTENT(IN OUT) :: NFCN
  LOGICAL, INTENT(OUT):: TER
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: J, IERR
  REAL(PREC) ::  ERRVJ, SUM, FAC, ZP(R)

  EXTERNAL FCN
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  !!--------- ONE STEP OF THE ITERATION PROCEDURE


  TER = .FALSE.

  SELECT CASE(IJOB)
  CASE(1,2)  !! ODE

     DO J=1,R
        SUM= B711*FP(J,1)+B712*FP(J,2)+B713*FP(J,3)&
             &          +B714*FP(J,4)+B715*FP(J,5)+B716*FP(J,6)+B717*FP(J,7)
        DN(J) =YP(J,2)-YP(J,1)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2) - DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )


     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L721*(F1(J,1)-FP(J,2))+B721*FP(J,1)+B722*FP(J,2)+&
             &B723*FP(J,3)+B724*FP(J,4)+B725*FP(J,5)+B726*FP(J,6)+B727*FP(J,7)
        DN(J) = YP(J,3)-YP(J,2)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = +L731*(F1(J,1)-FP(J,2))+L732*(F1(J,2) -FP(J,3))&
             & +B731*FP(J,1)+B732*FP(J,2)+B733*FP(J,3)&
             & +B734*FP(J,4)+B735*FP(J,5)+B736*FP(J,6)+B737*FP(J,7)
        DN(J) =YP(J,4)-YP(J,3)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM =L741*(F1(J,1)-FP(J,2))+L742*(F1(J,2)-FP(J,3))&
             & +L743*(F1(J,3)-FP(J,4))&
             & +B731*FP(J,2)+B732*FP(J,3)+B733*FP(J,4)&
             & +B734*FP(J,5)+B735*FP(J,6)+B736*FP(J,7)+B737*FP(J,8)
        DN(J) =YP(J,5)-YP(J,4)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4),  RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L751*(F1(J,1)-FP(J,2))+L752*(F1(J,2)-FP(J,3))&
             &+L753*(F1(J,3)-FP(J,4))+L754*(F1(J,4)-FP(J,5))&
             & +B731*FP(J,3)+B732*FP(J,4)+B733*FP(J,5)&
             & +B734*FP(J,6)+B735*FP(J,7)+B736*FP(J,8)+B737*FP(J,9)
        DN(J)=YP(J,6)-YP(J,5)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L761*(F1(J,1)-FP(J,2))+L762*(F1(J,2)-FP(J,3))&
             &+L763*(F1(J,3)-FP(J,4))+L764*(F1(J,4)-FP(J,5))&
             &+L765*(F1(J,5)-FP(J,6))&
             & +B737*FP(J,3)+B736*FP(J,4)+B735*FP(J,5)&
             & +B734*FP(J,6)+B733*FP(J,7)+B732*FP(J,8)+B731*FP(J,9)
        DN(J) = YP(J,7)-YP(J,6)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L771*(F1(J,1)-FP(J,2))+L772*(F1(J,2)-FP(J,3))&
        &+ L773*(F1(J,3)-FP(J,4))+L774*(F1(J,4)-FP(J,5))&
        &+ L775*(F1(J,5)-FP(J,6))+L776*(F1(J,6)-FP(J,7))&
        & +B727*FP(J,3)+B726*FP(J,4)+B725*FP(J,5)&
        & +B724*FP(J,6)+B723*FP(J,7)+B722*FP(J,8)+B721*FP(J,9)
        DN(J) = YP(J,8)-YP(J,7)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,8)=YP(J,8)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(8),YP(1,8),F1(1,7),  RPAR,IPAR)   ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF


     DO J=1,R
        SUM = L781*(F1(J,1)-FP(J,2))+L782*(F1(J,2)-FP(J,3))&
        &+ L783*(F1(J,3)-FP(J,4))+L784*(F1(J,4)-FP(J,5))&
        &+ L785*(F1(J,5)-FP(J,6))+L786*(F1(J,6)-FP(J,7))&
        &+ L787*(F1(J,7)-FP(J,8))&
        & +B717*FP(J,3)+B716*FP(J,4)+B715*FP(J,5)&
        & +B714*FP(J,6)+B713*FP(J,7)+B712*FP(J,8)+B711*FP(J,9)
        DN(J) = YP(J,9)-YP(J,8)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,9)=YP(J,9)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF


  CASE(3,4)  !! DAE: FMAS A BAND MATRIX

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,2)-YP(1:R,1))
     DO J=1,R
        SUM= B711*FP(J,1)+B712*FP(J,2)+B713*FP(J,3)&
             &          +B714*FP(J,4)+B715*FP(J,5)+B716*FP(J,6)+B717*FP(J,7)
        DN(J) =ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2) - DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )


     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,3)-YP(1:R,2))
     DO J=1,R
        SUM = L721*(F1(J,1)-FP(J,2))+B721*FP(J,1)+B722*FP(J,2)+&
             &B723*FP(J,3)+B724*FP(J,4)+B725*FP(J,5)+B726*FP(J,6)+B727*FP(J,7)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,4)-YP(1:R,3))
     DO J=1,R
        SUM = +L731*(F1(J,1)-FP(J,2))+L732*(F1(J,2) -FP(J,3))&
             & +B731*FP(J,1)+B732*FP(J,2)+B733*FP(J,3)&
             & +B734*FP(J,4)+B735*FP(J,5)+B736*FP(J,6)+B737*FP(J,7)
        DN(J) =ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,5)-YP(1:R,4))
     DO J=1,R
        SUM =L741*(F1(J,1)-FP(J,2))+L742*(F1(J,2)-FP(J,3))&
             & +L743*(F1(J,3)-FP(J,4))&
             & +B731*FP(J,2)+B732*FP(J,3)+B733*FP(J,4)&
             & +B734*FP(J,5)+B735*FP(J,6)+B736*FP(J,7)+B737*FP(J,8)
        DN(J) =ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,6)-YP(1:R,5))
     DO J=1,R
        SUM = L751*(F1(J,1)-FP(J,2))+L752*(F1(J,2)-FP(J,3))&
             &+L753*(F1(J,3)-FP(J,4))+L754*(F1(J,4)-FP(J,5))&
             & +B731*FP(J,3)+B732*FP(J,4)+B733*FP(J,5)&
             & +B734*FP(J,6)+B735*FP(J,7)+B736*FP(J,8)+B737*FP(J,9)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,7)-YP(1:R,6))
     DO J=1,R
        SUM = L761*(F1(J,1)-FP(J,2))+L762*(F1(J,2)-FP(J,3))&
             &+L763*(F1(J,3)-FP(J,4))+L764*(F1(J,4)-FP(J,5))&
             &+L765*(F1(J,5)-FP(J,6))&
             & +B737*FP(J,3)+B736*FP(J,4)+B735*FP(J,5)&
             & +B734*FP(J,6)+B733*FP(J,7)+B732*FP(J,8)+B731*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,8)-YP(1:R,7))
     DO J=1,R
        SUM = L771*(F1(J,1)-FP(J,2))+L772*(F1(J,2)-FP(J,3))&
        &+ L773*(F1(J,3)-FP(J,4))+L774*(F1(J,4)-FP(J,5))&
        &+ L775*(F1(J,5)-FP(J,6))+L776*(F1(J,6)-FP(J,7))&
        & +B727*FP(J,3)+B726*FP(J,4)+B725*FP(J,5)&
        & +B724*FP(J,6)+B723*FP(J,7)+B722*FP(J,8)+B721*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,8)=YP(J,8)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF


     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,9)-YP(1:R,8))
     DO J=1,R
        SUM = L781*(F1(J,1)-FP(J,2))+L782*(F1(J,2)-FP(J,3))&
        &+ L783*(F1(J,3)-FP(J,4))+L784*(F1(J,4)-FP(J,5))&
        &+ L785*(F1(J,5)-FP(J,6))+L786*(F1(J,6)-FP(J,7))&
        &+ L787*(F1(J,7)-FP(J,8))&
        & +B717*FP(J,3)+B716*FP(J,4)+B715*FP(J,5)&
        & +B714*FP(J,6)+B713*FP(J,7)+B712*FP(J,8)+B711*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,9)=YP(J,9)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

  CASE(5)    !! DAE: FMAS A FULL MATRIX

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,2)-YP(1:R,1))
     DO J=1,R
        SUM= B711*FP(J,1)+B712*FP(J,2)+B713*FP(J,3)&
             &          +B714*FP(J,4)+B715*FP(J,5)+B716*FP(J,6)+B717*FP(J,7)
        DN(J) =ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2) - DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )


     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,3)-YP(1:R,2))
     DO J=1,R
        SUM = L721*(F1(J,1)-FP(J,2))+B721*FP(J,1)+B722*FP(J,2)+&
             &B723*FP(J,3)+B724*FP(J,4)+B725*FP(J,5)+B726*FP(J,6)+B727*FP(J,7)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,4)-YP(1:R,3))
     DO J=1,R
        SUM = +L731*(F1(J,1)-FP(J,2))+L732*(F1(J,2) -FP(J,3))&
             & +B731*FP(J,1)+B732*FP(J,2)+B733*FP(J,3)&
             & +B734*FP(J,4)+B735*FP(J,5)+B736*FP(J,6)+B737*FP(J,7)
        DN(J) =ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,5)-YP(1:R,4))
     DO J=1,R
        SUM =L741*(F1(J,1)-FP(J,2))+L742*(F1(J,2)-FP(J,3))&
             & +L743*(F1(J,3)-FP(J,4))&
             & +B731*FP(J,2)+B732*FP(J,3)+B733*FP(J,4)&
             & +B734*FP(J,5)+B735*FP(J,6)+B736*FP(J,7)+B737*FP(J,8)
        DN(J) =ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,6)-YP(1:R,5))
     DO J=1,R
        SUM = L751*(F1(J,1)-FP(J,2))+L752*(F1(J,2)-FP(J,3))&
             &+L753*(F1(J,3)-FP(J,4))+L754*(F1(J,4)-FP(J,5))&
             & +B731*FP(J,3)+B732*FP(J,4)+B733*FP(J,5)&
             & +B734*FP(J,6)+B735*FP(J,7)+B736*FP(J,8)+B737*FP(J,9)
        DN(J)=ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,7)-YP(1:R,6))
     DO J=1,R
        SUM = L761*(F1(J,1)-FP(J,2))+L762*(F1(J,2)-FP(J,3))&
             &+L763*(F1(J,3)-FP(J,4))+L764*(F1(J,4)-FP(J,5))&
             &+L765*(F1(J,5)-FP(J,6))&
             & +B737*FP(J,3)+B736*FP(J,4)+B735*FP(J,5)&
             & +B734*FP(J,6)+B733*FP(J,7)+B732*FP(J,8)+B731*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6),  RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,8)-YP(1:R,7))
     DO J=1,R
        SUM = L771*(F1(J,1)-FP(J,2))+L772*(F1(J,2)-FP(J,3))&
        &+ L773*(F1(J,3)-FP(J,4))+L774*(F1(J,4)-FP(J,5))&
        &+ L775*(F1(J,5)-FP(J,6))+L776*(F1(J,6)-FP(J,7))&
        & +B727*FP(J,3)+B726*FP(J,4)+B725*FP(J,5)&
        & +B724*FP(J,6)+B723*FP(J,7)+B722*FP(J,8)+B721*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,8)=YP(J,8)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF


     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,9)-YP(1:R,8))
     DO J=1,R
        SUM = L781*(F1(J,1)-FP(J,2))+L782*(F1(J,2)-FP(J,3))&
        &+ L783*(F1(J,3)-FP(J,4))+L784*(F1(J,4)-FP(J,5))&
        &+ L785*(F1(J,5)-FP(J,6))+L786*(F1(J,6)-FP(J,7))&
        &+ L787*(F1(J,7)-FP(J,8))&
        & +B717*FP(J,3)+B716*FP(J,4)+B715*FP(J,5)&
        & +B714*FP(J,6)+B713*FP(J,7)+B712*FP(J,8)+B711*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,9)=YP(J,9)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF




  END SELECT

  FP(1:R,2:9) = F1(1:R,1:8)



  RETURN


END SUBROUTINE TERMNOT7
    !!
    !!  SUBROUTINE TERMNOT9  (ORDER 9)
    !!
SUBROUTINE  TERMNOT9(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,                    &
                    &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV,            &
                    &  FMAS,LDMAS,MLMAS,MUMAS, SCAL,IJOB,TER,RPAR,IPAR)

  USE PRECISION
  IMPLICIT NONE

  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) :: R,  IJOB, IPIV(R), LDLU, LDMAS, MLMAS, MUMAS, IT, IPAR(1)
  REAL(PREC), INTENT(IN) ::  H, SCAL(R), LU(LDLU,R), FMAS(LDMAS,R), ERRNEWT0, TETAK0, RPAR(1)
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) :: ERRNEWT, YP(R,10), FP(R,10), &
       &                   DN(R), F1(R,10),TP(10)
  INTEGER, INTENT(IN OUT) :: NFCN
  LOGICAL, INTENT(OUT):: TER
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: J, IERR
  REAL(PREC) ::  ERRVJ, SUM, FAC, ZP(R)

  EXTERNAL FCN
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  !!--------- ONE STEP OF THE ITERATION PROCEDURE


  TER = .FALSE.
  SELECT CASE(IJOB)
  CASE(1,2)  !! ODE

     DO J=1,R
        SUM = B911*FP(J,1)+B912*FP(J,2)+B913*FP(J,3)+B914*FP(J,4)&
             &+B915*FP(J,5)+B916*FP(J,6)+B917*FP(J,7)+B918*FP(J,8)+B919*FP(J,9)
        DN(J) = YP(J,2)-YP(J,1)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L921*(F1(J,1)-FP(J,2))+B921*FP(J,1)+B922*FP(J,2)&
             &   +B923*FP(J,3)+B924*FP(J,4)+B925*FP(J,5)&
             &   +B926*FP(J,6)+B927*FP(J,7)+B928*FP(J,8)+B929*FP(J,9)
        DN(J) = YP(J,3)-YP(J,2)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L932*(F1(J,2)-FP(J,3))&
             &+B931*FP(J,1)+B932*FP(J,2)+B933*FP(J,3)+B934*FP(J,4)+B935*FP(J,5)&
             &+B936*FP(J,6)+B937*FP(J,7)+B938*FP(J,8)+B939*FP(J,9)
             DN(J) = YP(J,4)-YP(J,3)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L943*(F1(J,3)-FP(J,4))+B941*FP(J,1)+B942*FP(J,2)&
             & +B943*FP(J,3)+B944*FP(J,4)+B945*FP(J,5)+B946*FP(J,6)&
             & +B947*FP(J,7)+B948*FP(J,8)+B949*FP(J,9)
        DN(J) =YP(J,5)-YP(J,4)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L954*(F1(J,4)-FP(J,5))+B941*FP(J,2)+B942*FP(J,3)&
             & +B943*FP(J,4)+B944*FP(J,5)+B945*FP(J,6)+B946*FP(J,7)&
             & +B947*FP(J,8)+B948*FP(J,9)+B949*FP(J,10)
        DN(J) =YP(J,6)-YP(J,5)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
           SUM = L965*(F1(J,5)-FP(J,6))+B949*FP(J,2)+B948*FP(J,3)&
             & +B947*FP(J,4)+B946*FP(J,5)+B945*FP(J,6)+B944*FP(J,7)&
             & +B943*FP(J,8)+B942*FP(J,9)+B941*FP(J,10)
        DN(J) =YP(J,7)-YP(J,6)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L976*(F1(J,6)-FP(J,7))+B939*FP(J,2)+B938*FP(J,3)&
             & +B937*FP(J,4)+B936*FP(J,5)+B935*FP(J,6)+B934*FP(J,7)&
             & +B933*FP(J,8)+B932*FP(J,9)+B931*FP(J,10)
        DN(J) =YP(J,8)-YP(J,7)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,8)=YP(J,8)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L987*(F1(J,7)-FP(J,8))+B929*FP(J,2)+B928*FP(J,3)&
             & +B927*FP(J,4)+B926*FP(J,5)+B925*FP(J,6)+B924*FP(J,7)&
             & +B923*FP(J,8)+B922*FP(J,9)+B921*FP(J,10)
        DN(J) =YP(J,9)-YP(J,8)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,9)=YP(J,9)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L998*(F1(J,8)-FP(J,9))+B919*FP(J,2)+B918*FP(J,3)&
             & +B917*FP(J,4)+B916*FP(J,5)+B915*FP(J,6)+B914*FP(J,7)&
             & +B913*FP(J,8)+B912*FP(J,9)+B911*FP(J,10)
        DN(J) =YP(J,10)-YP(J,9)-H*SUM
     END DO

     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,10)=YP(J,10)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(10),YP(1,10),F1(1,9),  RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF


  CASE(3,4)  !! DAE: FMAS A BAND MATRIX

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,2)-YP(1:R,1))
     DO J=1,R
        SUM = B911*FP(J,1)+B912*FP(J,2)+B913*FP(J,3)+B914*FP(J,4)&
             &+B915*FP(J,5)+B916*FP(J,6)+B917*FP(J,7)+B918*FP(J,8)+B919*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,3)-YP(1:R,2))
     DO J=1,R
        SUM = L921*(F1(J,1)-FP(J,2))+B921*FP(J,1)+B922*FP(J,2)&
             &   +B923*FP(J,3)+B924*FP(J,4)+B925*FP(J,5)&
             &   +B926*FP(J,6)+B927*FP(J,7)+B928*FP(J,8)+B929*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,4)-YP(1:R,3))
     DO J=1,R
        SUM = L932*(F1(J,2)-FP(J,3))&
             &+B931*FP(J,1)+B932*FP(J,2)+B933*FP(J,3)+B934*FP(J,4)+B935*FP(J,5)&
             &+B936*FP(J,6)+B937*FP(J,7)+B938*FP(J,8)+B939*FP(J,9)
             DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,5)-YP(1:R,4))
     DO J=1,R
        SUM = L943*(F1(J,3)-FP(J,4))+B941*FP(J,1)+B942*FP(J,2)&
             & +B943*FP(J,3)+B944*FP(J,4)+B945*FP(J,5)+B946*FP(J,6)&
             & +B947*FP(J,7)+B948*FP(J,8)+B949*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,6)-YP(1:R,5))
     DO J=1,R
        SUM = L954*(F1(J,4)-FP(J,5))+B941*FP(J,2)+B942*FP(J,3)&
             & +B943*FP(J,4)+B944*FP(J,5)+B945*FP(J,6)+B946*FP(J,7)&
             & +B947*FP(J,8)+B948*FP(J,9)+B949*FP(J,10)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,7)-YP(1:R,6))
     DO J=1,R
           SUM = L965*(F1(J,5)-FP(J,6))+B949*FP(J,2)+B948*FP(J,3)&
             & +B947*FP(J,4)+B946*FP(J,5)+B945*FP(J,6)+B944*FP(J,7)&
             & +B943*FP(J,8)+B942*FP(J,9)+B941*FP(J,10)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,8)-YP(1:R,7))
     DO J=1,R
        SUM = L976*(F1(J,6)-FP(J,7))+B939*FP(J,2)+B938*FP(J,3)&
             & +B937*FP(J,4)+B936*FP(J,5)+B935*FP(J,6)+B934*FP(J,7)&
             & +B933*FP(J,8)+B932*FP(J,9)+B931*FP(J,10)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,8)=YP(J,8)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,9)-YP(1:R,8))
     DO J=1,R
        SUM = L987*(F1(J,7)-FP(J,8))+B929*FP(J,2)+B928*FP(J,3)&
             & +B927*FP(J,4)+B926*FP(J,5)+B925*FP(J,6)+B924*FP(J,7)&
             & +B923*FP(J,8)+B922*FP(J,9)+B921*FP(J,10)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,9)=YP(J,9)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,10)-YP(1:R,9))
     DO J=1,R
        SUM = L998*(F1(J,8)-FP(J,9))+B919*FP(J,2)+B918*FP(J,3)&
             & +B917*FP(J,4)+B916*FP(J,5)+B915*FP(J,6)+B914*FP(J,7)&
             & +B913*FP(J,8)+B912*FP(J,9)+B911*FP(J,10)
        DN(J) = ZP(J)-H*SUM
     END DO

     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,10)=YP(J,10)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(10),YP(1,10),F1(1,9),  RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

  CASE(5)    !! DAE: FMAS A FULL MATRIX

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,2)-YP(1:R,1))
     DO J=1,R
        SUM = B911*FP(J,1)+B912*FP(J,2)+B913*FP(J,3)+B914*FP(J,4)&
             &+B915*FP(J,5)+B916*FP(J,6)+B917*FP(J,7)+B918*FP(J,8)+B919*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,3)-YP(1:R,2))
     DO J=1,R
        SUM = L921*(F1(J,1)-FP(J,2))+B921*FP(J,1)+B922*FP(J,2)&
             &   +B923*FP(J,3)+B924*FP(J,4)+B925*FP(J,5)&
             &   +B926*FP(J,6)+B927*FP(J,7)+B928*FP(J,8)+B929*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,4)-YP(1:R,3))
     DO J=1,R
        SUM = L932*(F1(J,2)-FP(J,3))&
             &+B931*FP(J,1)+B932*FP(J,2)+B933*FP(J,3)+B934*FP(J,4)+B935*FP(J,5)&
             &+B936*FP(J,6)+B937*FP(J,7)+B938*FP(J,8)+B939*FP(J,9)
             DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,5)-YP(1:R,4))
     DO J=1,R
        SUM = L943*(F1(J,3)-FP(J,4))+B941*FP(J,1)+B942*FP(J,2)&
             & +B943*FP(J,3)+B944*FP(J,4)+B945*FP(J,5)+B946*FP(J,6)&
             & +B947*FP(J,7)+B948*FP(J,8)+B949*FP(J,9)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,6)-YP(1:R,5))
     DO J=1,R
        SUM = L954*(F1(J,4)-FP(J,5))+B941*FP(J,2)+B942*FP(J,3)&
             & +B943*FP(J,4)+B944*FP(J,5)+B945*FP(J,6)+B946*FP(J,7)&
             & +B947*FP(J,8)+B948*FP(J,9)+B949*FP(J,10)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,7)-YP(1:R,6))
     DO J=1,R
           SUM = L965*(F1(J,5)-FP(J,6))+B949*FP(J,2)+B948*FP(J,3)&
             & +B947*FP(J,4)+B946*FP(J,5)+B945*FP(J,6)+B944*FP(J,7)&
             & +B943*FP(J,8)+B942*FP(J,9)+B941*FP(J,10)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,8)-YP(1:R,7))
     DO J=1,R
        SUM = L976*(F1(J,6)-FP(J,7))+B939*FP(J,2)+B938*FP(J,3)&
             & +B937*FP(J,4)+B936*FP(J,5)+B935*FP(J,6)+B934*FP(J,7)&
             & +B933*FP(J,8)+B932*FP(J,9)+B931*FP(J,10)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,8)=YP(J,8)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,9)-YP(1:R,8))
     DO J=1,R
        SUM = L987*(F1(J,7)-FP(J,8))+B929*FP(J,2)+B928*FP(J,3)&
             & +B927*FP(J,4)+B926*FP(J,5)+B925*FP(J,6)+B924*FP(J,7)&
             & +B923*FP(J,8)+B922*FP(J,9)+B921*FP(J,10)
        DN(J) = ZP(J)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,9)=YP(J,9)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,10)-YP(1:R,9))
     DO J=1,R
        SUM = L998*(F1(J,8)-FP(J,9))+B919*FP(J,2)+B918*FP(J,3)&
             & +B917*FP(J,4)+B916*FP(J,5)+B915*FP(J,6)+B914*FP(J,7)&
             & +B913*FP(J,8)+B912*FP(J,9)+B911*FP(J,10)
        DN(J) = ZP(J)-H*SUM
     END DO

     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,10)=YP(J,10)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(10),YP(1,10),F1(1,9), RPAR,IPAR)  ! removed ierr
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF



  END SELECT


  FP(1:R,2:10) = F1(1:R,1:9)

  RETURN

  RETURN
END SUBROUTINE TERMNOT9


    !!
    !!  SUBROUTINE ESTERR
    !! ERRORS ESTIMATION.    ERRSAME: THE CURRENT ORDER
    !!                         ERRUP: GREATER ORDER (THAN ERRSAME)
    !!                       ERRDOWN: LOWER ORDER
    !!

SUBROUTINE  ESTERR(ERRV, ERRSAME, ERRUP, ERRDOWN, FP,            &
     &     R, H, ORD, DBLK, LU, LDLU, FMAS, LDMAS, MLMAS, MUMAS, &
     &     IPIV, F, DN,SCAL,ORDMAX,ORDMIN,IJOB)


  !!
  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) :: R, ORD,ORDMIN, ORDMAX,DBLK,IJOB,IPIV(R),LDLU, LDMAS, MLMAS, MUMAS
  REAL(PREC), INTENT(IN) :: H, SCAL(R),  LU(LDLU,R), FMAS(LDMAS,R)
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) ::  ERRV(DBLK), ERRSAME, ERRUP, ERRDOWN,   &
       &                  FP(R,11), F(R,11), DN(R,11)
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: I, J
  REAL(PREC) ::  ERRVJ,   &
       &              FP1,FP2,FP3,FP4,FP5,FP6,FP7,FP8,FP9,FP10
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!--------- ERRSAME ESTIMATION
  SELECT CASE(ORD)
  CASE(1)

     DO I=1,R
        FP1= FP(I,1)
        FP2= FP(I,2)
        FP3= FP(I,3)
        FP4= FP(I,4)
        FP5= FP(I,5)
        F(I,1) = H*(B3511*FP1+B3512*FP2+B3513*FP3+B3514*FP4+B3515*FP5)
        F(I,2) = H*(B3521*FP1+B3522*FP2+B3523*FP3+B3524*FP4+B3525*FP5)
        F(I,3) = H*(B3531*FP1+B3532*FP2+B3533*FP3+B3534*FP4+B3535*FP5)
        F(I,4) = H*(B3541*FP1+B3542*FP2+B3543*FP3+B3544*FP4+B3545*FP5)
       END DO

  CASE(2)
     DO J=1,R
        FP1= FP(J,1)
        FP2= FP(J,2)
        FP3= FP(J,3)
        FP4= FP(J,4)
        FP5= FP(J,5)
        FP6= FP(J,6)
        FP7= FP(J,7)
        F(J,1) = H*(B5711*FP1+B5712*FP2+B5713*FP3&
             &            +B5714*FP4+B5715*FP5+B5716*FP6+B5717*FP7)
        F(J,2) = H*(B5721*FP1+B5722*FP2+B5723*FP3&
             &            +B5724*FP4+B5725*FP5+B5726*FP6+B5727*FP7)
        F(J,3) = H*(B5731*FP1+B5732*FP2+B5733*FP3&
             &            +B5734*FP4+B5735*FP5+B5736*FP6+B5737*FP7)
        F(J,4) = H*(B5741*FP1+B5742*FP2+B5743*FP3&
             &            +B5744*FP4+B5745*FP5+B5746*FP6+B5747*FP7)
        F(J,5) = H*(B5727*FP1+B5726*FP2+B5725*FP3&
             &            +B5724*FP4+B5723*FP5+B5722*FP6+B5721*FP7)
        F(J,6) = H*(B5717*FP1+B5716*FP2+B5715*FP3&
             &            +B5714*FP4+B5713*FP5+B5712*FP6+B5711*FP7)
     END DO

  CASE(3)
     DO J=1,R
        FP1= FP(J,1)
        FP2= FP(J,2)
        FP3= FP(J,3)
        FP4= FP(J,4)
        FP5= FP(J,5)
        FP6= FP(J,6)
        FP7= FP(J,7)
        FP8= FP(J,8)
        FP9= FP(J,9)
        F(J,1) = H*(B7911*FP1+B7912*FP2+B7913*FP3+B7914*FP4&
             & +B7915*FP5+B7916*FP6+B7917*FP7+B7918*FP8+B7919*FP9)
        F(J,2) = H*(B7921*FP1+B7922*FP2+B7923*FP3+B7924*FP4&
             & +B7925*FP5+B7926*FP6+B7927*FP7+B7928*FP8+B7929*FP9)
        F(J,3) = H*(B7931*FP1+B7932*FP2+B7933*FP3+B7934*FP4&
             & +B7935*FP5+B7936*FP6+B7937*FP7+B7938*FP8+B7939*FP9)
        F(J,4) = H*(B7941*FP1+B7942*FP2+B7943*FP3+B7944*FP4&
             & +B7945*FP5+B7946*FP6+B7947*FP7+B7948*FP8+B7949*FP9)
        F(J,5) = H*(B7951*FP1+B7952*FP2+B7953*FP3+B7954*FP4&
             & +B7955*FP5+B7956*FP6+B7957*FP7+B7958*FP8+B7959*FP9)
        F(J,6) = H*(B7939*FP1+B7938*FP2+B7937*FP3+B7936*FP4&
             & +B7935*FP5+B7934*FP6+B7933*FP7+B7932*FP8+B7931*FP9)
        F(J,7) = H*(B7929*FP1+B7928*FP2+B7927*FP3+B7926*FP4&
             & +B7925*FP5+B7924*FP6+B7923*FP7+B7922*FP8+B7921*FP9)
        F(J,8) = H*(B7919*FP1+B7918*FP2+B7917*FP3+B7916*FP4&
             & +B7915*FP5+B7914*FP6+B7913*FP7+B7912*FP8+B7911*FP9)
     END DO

  CASE(4)
     DO J=1,R
        FP1= FP(J,1)
        FP2= FP(J,2)
        FP3= FP(J,3)
        FP4= FP(J,4)
        FP5= FP(J,5)
        FP6= FP(J,6)
        FP7= FP(J,7)
        FP8= FP(J,8)
        FP9= FP(J,9)
        FP10= FP(J,10)
        F(J,1) = H*(B91011*(FP1-FP10)+B91012*(FP2-FP9)+B91013*(FP3-FP8)&
             & +B91014*(FP4-FP7)+B91015*(FP5-FP6) )
        F(J,2) = H*(B91021*(FP1-FP10)+B91022*(FP2-FP9)+B91023*(FP3-FP8)&
             & +B91024*(FP4-FP7)+B91025*(FP5-FP6) )
        F(J,3) = H*(B91031*(FP1-FP10)+B91032*(FP2-FP9)+B91033*(FP3-FP8)&
             & +B91034*(FP4-FP7)+B91035*(FP5-FP6) )
        F(J,4) = H*(B91041*(FP1-FP10)+B91042*(FP2-FP9)+B91043*(FP3-FP8)&
             & +B91044*(FP4-FP7)+B91045*(FP5-FP6) )
        F(J,5) =  F(J,4)
        F(J,6) = -F(J,4)
        F(J,7) = -F(J,3)
        F(J,8) = -F(J,2)
        F(J,9) = -F(J,1)
     END DO

  END SELECT

  !!--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRSAME

  CALL NEWTGS(R,DBLK,LU(1,1),LDLU,FMAS(1,1),LDMAS,MLMAS,MUMAS,&
     &      IPIV(1),F(1,1),DN(1,1),IJOB)


  !!--------- COMPUTE  ERRSAME AND ERRV (VECTOR ERROR)

  ERRSAME = 0D0
  DO J=1,DBLK
     ERRV(J) = 0D0
     DO I=1,R
        FP1 = (DN(I,J)/SCAL(I) )
        ERRV(J) =  ERRV(J)+ FP1*FP1
     END DO
     ERRV(J) = SQRT(ERRV(J)/R)
     ERRSAME = MAX( ERRSAME, ERRV(J) )
  END DO
  ERRSAME = MAX(ERRSAME, 1d-15)
  ERRDOWN = 0D0
  ERRUP   = 0D0
  IF ( ERRSAME .LE. 1d0) THEN

     IF (ORD .LT. ORDMAX) THEN
        !!--------- ERRUP ESTIMATION
        SELECT CASE(ORD)
        CASE(1)
           DO I=1,R
              FP1 =  F(I,1)/CP31
              FP2 =  F(I,2)/CP31
              FP3 =  F(I,3)/CP31
              FP4 = -F(I,4)/CP31
              F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP51
              F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP52
              F(I,3) = -(-FP2 + 2d0*FP3 - FP4)*CP52
              F(I,4) = -(-FP2 + 2d0*FP3 - FP4)*CP51
           END DO

        CASE(2)
           DO I = 1, R
              FP1 =  F(I,1)/CP51
              FP2 =  F(I,2)/CP52
              FP3 =  F(I,3)/CP52
              FP4 =  F(I,4)/CP52
              FP5 = -F(I,5)/CP52
              FP6 = -F(I,6)/CP51
              F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP71
              F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP72
              F(I,3) =  (-FP2 + 2d0*FP3 - FP4)*CP73
              F(I,4) = -(-FP3 + 2d0*FP4 - FP5)*CP73
              F(I,5) = -(-FP4 + 2d0*FP5 - FP6)*CP72
              F(I,6) = -(-FP4 + 2d0*FP5 - FP6)*CP71

           END DO
        CASE(3)

           DO I = 1, R
              FP1 =  F(I,1)/CP71
              FP2 =  F(I,2)/CP72
              FP3 =  F(I,3)/CP73
              FP4 =  F(I,4)/CP73
              FP5 =  F(I,5)/CP73
              FP6 = -F(I,6)/CP73
              FP7 = -F(I,7)/CP72
              FP8 = -F(I,8)/CP71
              F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP91
              F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP92
              F(I,3) =  (-FP2 + 2d0*FP3 - FP4)*CP93
              F(I,4) =  (-FP3 + 2d0*FP4 - FP5)*CP94
              F(I,5) = -(-FP4 + 2d0*FP5 - FP6)*CP94
              F(I,6) = -(-FP5 + 2d0*FP6 - FP7)*CP93
              F(I,7) = -(-FP6 + 2d0*FP7 - FP8)*CP92
              F(I,8) = -(-FP6 + 2d0*FP7 - FP8)*CP91
           END DO

        END SELECT

        !!--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRUP


        CALL NEWTGS(R,DBLK,LU(1,1),LDLU,FMAS(1,1),LDMAS,MLMAS,MUMAS,IPIV(1),F(1,1),DN(1,1),IJOB)


        !!-------- COMPUTE ERRUP
        ERRUP = 0D0
        DO J=1,DBLK
           ERRVJ = 0D0
           DO I=1,R
              FP1 = (DN(I,J)/SCAL(I) )
              ERRVJ =  ERRVJ + FP1*FP1
           END DO
           ERRVJ = SQRT(ERRVJ/R)
           ERRUP = MAX( ERRUP, ERRVJ )
        END DO
        ERRUP = max(ERRUP, 1d-15)

     END IF

     IF (ORD .GT. ORDMIN) THEN
        !!--------- ERRDOWN ESTIMATION
        SELECT CASE(ORD)
        CASE(2)

           DO J=1,R
              FP1= FP(J,1)
              FP2= FP(J,2)
              FP3= FP(J,3)
              FP4= FP(J,4)
              FP5= FP(J,5)
              FP6= FP(J,6)
              FP7= FP(J,7)
              F(J,1) = H*(B3511*FP1+B3512*FP2+B3513*FP3+B3514*FP4+B3515*FP5)
              F(J,2) = H*(B3521*FP1+B3522*FP2+B3523*FP3+B3524*FP4+B3525*FP5)
              F(J,3) = H*(B3521*FP2+B3522*FP3+B3523*FP4+B3524*FP5+B3525*FP6)
              F(J,4) = H*(B3521*FP3+B3522*FP4+B3523*FP5+B3524*FP6+B3525*FP7)
              F(J,5) = H*(B3531*FP3+B3532*FP4+B3533*FP5+B3534*FP6+B3535*FP7)
              F(J,6) = H*(B3541*FP3+B3542*FP4+B3543*FP5+B3544*FP6+B3545*FP7)
           END DO
        CASE(3)

           DO J=1,R
              FP1= FP(J,1)
              FP2= FP(J,2)
              FP3= FP(J,3)
              FP4= FP(J,4)
              FP5= FP(J,5)
              FP6= FP(J,6)
              FP7= FP(J,7)
              FP8= FP(J,8)
              FP9= FP(J,9)
              F(J,1) = H*(B5711*FP1+B5712*FP2+B5713*FP3   &
                   &            +B5714*FP4+B5715*FP5+B5716*FP6+B5717*FP7)
              F(J,2) = H*(B5721*FP1+B5722*FP2+B5723*FP3   &
                   &            +B5724*FP4+B5725*FP5+B5726*FP6+B5727*FP7)
              F(J,3) = H*(B5731*FP1+B5732*FP2+B5733*FP3   &
                   &            +B5734*FP4+B5735*FP5+B5736*FP6+B5737*FP7)
              F(J,4) = H*(B5731*FP2+B5732*FP3+B5733*FP4   &
                   &            +B5734*FP5+B5735*FP6+B5736*FP7+B5737*FP8)
              F(J,5) = H*(B5731*FP3+B5732*FP4+B5733*FP5   &
                   &            +B5734*FP6+B5735*FP7+B5736*FP8+B5737*FP9)
              F(J,6) = H*(B5741*FP3+B5742*FP4+B5743*FP5   &
                   &            +B5744*FP6+B5745*FP7+B5746*FP8+B5747*FP9)
              F(J,7) = H*(B5727*FP3+B5726*FP4+B5725*FP5   &
                   &            +B5724*FP6+B5723*FP7+B5722*FP8+B5721*FP9)
              F(J,8) = H*(B5717*FP3+B5716*FP4+B5715*FP5   &
                   &            +B5714*FP6+B5713*FP7+B5712*FP8+B5711*FP9)

           END DO
        CASE(4)


           DO J=1,R
              FP1= FP(J,1)
              FP2= FP(J,2)
              FP3= FP(J,3)
              FP4= FP(J,4)
              FP5= FP(J,5)
              FP6= FP(J,6)
              FP7= FP(J,7)
              FP8= FP(J,8)
              FP9= FP(J,9)
              FP10= FP(J,10)
              F(J,1) = H*(B7911*FP1+B7912*FP2+B7913*FP3+B7914*FP4&
                   & +B7915*FP5+B7916*FP6+B7917*FP7+B7918*FP8+B7919*FP9)

              F(J,2) = H*(B7921*FP1+B7922*FP2+B7923*FP3+B7924*FP4&
                   & +B7925*FP5+B7926*FP6+B7927*FP7+B7928*FP8+B7929*FP9)

              F(J,3) = H*(B7931*FP1+B7932*FP2+B7933*FP3+B7934*FP4&
                   & +B7935*FP5+B7936*FP6+B7937*FP7+B7938*FP8+B7939*FP9)

              F(J,4) = H*(B7941*FP1+B7942*FP2+B7943*FP3+B7944*FP4&
                   & +B7945*FP5+B7946*FP6+B7947*FP7+B7948*FP8+B7949*FP9)

              F(J,5) = -H*(B7941*FP2+B7942*FP3+B7943*FP4+B7944*FP5&
                   & +B7945*FP6+B7946*FP7+B7947*FP8+B7948*FP9+B7949*FP10)


              F(J,6) = -H*(B7951*FP2+B7952*FP3+B7953*FP4+B7954*FP5&
                   & +B7955*FP6+B7956*FP7+B7957*FP8+B7958*FP9+B7959*FP10)

              F(J,7) = -H*(B7939*FP2+B7938*FP3+B7937*FP4+B7936*FP5&
                   & +B7935*FP6+B7934*FP7+B7933*FP8+B7932*FP9+B7931*FP10)

              F(J,8) = -H*(B7929*FP2+B7928*FP3+B7927*FP4+B7926*FP5&
                   & +B7925*FP6+B7924*FP7+B7923*FP8+B7922*FP9+B7921*FP10)

              F(J,9) = -H*(B7919*FP2+B7918*FP3+B7917*FP4+B7916*FP5&
                   & +B7915*FP6+B7914*FP7+B7913*FP8+B7912*FP9+B7911*FP10)

           END DO

        END SELECT
        !!--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRDOWN

        CALL NEWTGS(R,DBLK,LU(1,1),LDLU,FMAS(1,1),LDMAS,MLMAS,MUMAS,IPIV(1),F(1,1),DN(1,1),IJOB)


        !!--------- COMPUTE ERRDOWN
        ERRDOWN = 0D0
        DO J=1,DBLK
           ERRVJ = 0D0
           DO I=1,R
              FP1 = (DN(I,J)/SCAL(I) )
              ERRVJ =  ERRVJ + FP1*FP1
           END DO
           ERRVJ = SQRT(ERRVJ/R)
           ERRDOWN = MAX(ERRDOWN, ERRVJ )
        END DO
        ERRDOWN = max(ERRDOWN, 1d-15)
     END IF

  END IF
  RETURN
END SUBROUTINE ESTERR


END MODULE SUBGAMD



SUBROUTINE   GAMD(R,FCN,T0,Y0,TEND,H,            &
       &                  RTOL,ATOL,ITOL,        &
       &                  JAC ,IJAC,MLJAC,MUJAC, &
       &                  MAS ,IMAS,MLMAS,MUMAS, &
       &                  SOLOUT,IOUT,           &
       &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

    !!
    !!     PURPOSE: THE CODE GAMD NUMERICALLY SOLVES A (POSSIBLY STIFF)
    !!              SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS
    !!              IN THE FORM  Y'=F(T,Y), OR A LINEARLY IMPLICIT DAE
    !!              MY'=F(T,Y) WITH CONSTANT MASS (FULL OR BANDED) MATRIX M,
    !!              WITH A GIVEN INITIAL CONDITION. IT IS INTENDED AS AN
    !!              EXTENSION OF THE CODE GAM TO DAEs (SEE THE REVISION
    !!              HISTORY FOR DETAILS)
    !!
    !!     AUTHORS: F. IAVERNARO AND F. MAZZIA
    !!              UNIVERSITA' DEGLI STUDI DI BARI,
    !!              DIPARTIMENTO DI MATEMATICA
    !!              VIA ORABONA 4, 70125 BARI, ITALY
    !!              E-MAIL:  MAZZIA@DM.UNIBA.IT
    !!                       FELIX@DM.UNIBA.IT
    !!
    !!     METHODS: THE METHODS USED ARE IN THE CLASS OF BOUNDARY VALUE
    !!              METHODS (BVMs), NAMELY THE GENERALIZED ADAMS METHODS
    !!              (GAMs) OF ORDER 3-5-7-9 WITH STEP SIZE CONTROL.
    !!
    !!  REFERENCES: L.BRUGNANO, D.TRIGIANTE,  Solving Differential Problems
    !!              by Multistep Initial and Boundary Value Methods,
    !!              Gordon & Breach,1998.
    !!
    !!              F.IAVERNARO, F.MAZZIA,  Block-Boundary Value Methods
    !!              for the solution of Ordinary Differential Equations,
    !!              Siam J. Sci. Comput. 21 (1) (1999) 323--339.
    !!
    !!              F.IAVERNARO, F.MAZZIA,  Solving Ordinary Differential
    !!              Equations by Generalized Adams Methods: properties and
    !!              implementation techniques,
    !!              Appl. Num. Math. 28 (2-4) (1998) 107-126.
    !!
    !! DESCRIPTION: THE CODE GAMD CONSISTS OF TWO parts:
    !!              - FROM HERE ON: THE MAIN SUBROUTINES THAT IMPLEMENT THE
    !!                           INTEGRATION PROCEDURE;
    !!
    !!              - BEFOR:  THE ADDITIONAL LINEAR ALGEBRA ROUTINES
    !!                        REQUIRED BY gamd.f90 PLUS SOME OTHER SUBROUTINES
    !!                        PROPER OF THE USED METHODS;
    !!
    !!    COMMENTS: THE PHILOSOFY AND THE STYLE USED IN WRITING THE CODE ARE VERY
    !!              SIMILAR TO THOSE CHARACTERIZING THE CODE  RADAU5.
    !!              INDEED THE AUTHORS IMPORTED FROM RADAU5 SOME SUBROUTINES,
    !!              COMMENTS AND IMPLEMENTATION TECHNIQUES LEAVING UNCHANGED
    !!              THE NAME AND THE MEANING OF  A NUMBER OF VARIABLES.
    !!              THE AUTHORS ARE VERY GRATEFUL TO ANYONE USES THE CODE AND
    !!              WOULD APPRECIATE ANY CRITICISM AND REMARKS ON HOW IT PERFORMS.
    !!
    !! REVISION HISTORY (YYYY/MM/DD):
    !!
    !!                   1997/20/08
    !!                      - FIRST VERSION OF GAM
    !!
    !!                   1999/11/25
    !!                      - CORRECTED OUTPUT LAST STEPSIZE
    !!                      - CORRECTED INPUT  H
    !!
    !!                   2003/23/08
    !!                      - FIRST VERSION OF GAMD
    !!                      - REWRITTEN IN FORTRAN 90
    !!                      - EXTENDED TO LINEARLY IMPLICIT DAES  MY'=F(T,Y)
    !!                        WITH CONSTANT MASS MATRIX M
    !!                      - ADDED IERR IN FCN
    !!
    !!                   2006/24/01
    !!                      - CORRECTED THE VALUE CALJAC TO AVOID
    !!                        THE COMPUTATION OF JACOBIAN MORE THEN ONE
    !!                        TIME PER STEPS
    !!                      - CHANGED THE DEFINITION OF TP AND T1 IN ALL
    !!                        THE SUBROUTINE TP(DBLK+1) --> TP(*)
    !!
    !!                   2006/15/02
    !!                      - CORRECTED ALL THE RUN TIM ERROR FOR
    !!                        SALFORD FTN95 FORTRAN
    !!
    !!                   2007/15/03
    !!                       - CORRECTED INSTRUCTION FOR WORK INDEX 5:7
    !!                   2007/24/05
    !!                       - ADDED THE GNU General Public License
    !! -------------------------------------------------------------------------
    !!     INPUT PARAMETERS
    !! -------------------------------------------------------------------------
    !!     R           DIMENSION OF THE SYSTEM
    !!
    !!     FCN         NAME (EXTERNAL) OF THE SUBROUTINE COMPUTING THE
    !!                 VALUE OF F(T,Y):
    !!                    SUBROUTINE FCN(R,T,Y,F,IERR,RPAR,IPAR)
    !!                    DOUBLE PRECISION X,Y(R),F(R)
    !!                    INTEGER IERR
    !!                    F(1)=...   ETC.
    !!                  IERR = -1 to prevent OVERFLOW
    !!                 (RPAR, IPAR    SEE BELOW)
    !!
    !!     T0          INITIAL T-VALUE
    !!
    !!     Y0          INITIAL VALUES FOR Y
    !!
    !!     TEND        FINAL T-VALUE (TEND-T MUST BE POSITIVE)
    !!
    !!     H           INITIAL STEP SIZE GUESS;
    !!                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
    !!                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD.
    !!                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS
    !!                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6).
    !!
    !!     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
    !!                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH R.
    !!
    !!     ITOL        SWITCH FOR RTOL AND ATOL:
    !!                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
    !!                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
    !!                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
    !!                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
    !!                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
    !!                     RTOL(I)*ABS(Y(I))+ATOL(I).
    !!
    !! Modified to be used in R
    !!     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
    !!                 THE PARTIAL DERIVATIVES OF F(T,Y) WITH RESPECT TO Y
    !!                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
    !!                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
    !!                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
    !!                    SUBROUTINE JAC(R,T,ML,MU,Y,DFY,LDFY,RPAR,IPAR)
    !!                    DOUBLE PRECISION T,Y(R),DFY(LDFY,R)
    !!                    DFY(1,1)= ...
    !!                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS
    !!                 FURNISHED BY THE CALLING PROGRAM.
    !!                 IF (MLJAC.EQ.R) THE JACOBIAN IS SUPPOSED TO
    !!                    BE FULL AND THE PARTIAL DERIVATIVES ARE
    !!                    STORED IN DFY AS
    !!                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
    !!                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
    !!                    THE PARTIAL DERIVATIVES ARE STORED
    !!                    DIAGONAL-WISE AS
    !!                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
    !!
    !!     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
    !!                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
    !!                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
    !!                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
    !!
    !!     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
    !!                    MLJAC=R: JACOBIAN IS A FULL MATRIX. THE LINEAR
    !!                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
    !!                       0<=MLJAC<R: MLJAC IS THE LOWER BANDWITH OF JACOBIAN
    !!                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
    !!                       THE MAIN DIAGONAL).
    !!
    !!     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
    !!                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
    !!                 NEED NOT BE DEFINED IF MLJAC=R.
    !!
    !!     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
    !!                 MATRIX M.
    !!                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
    !!                 MATRIX AND NEEDS NOT TO BE DEFINED;
    !!                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
    !!                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
    !!                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR)
    !!                    DOUBLE PRECISION AM(LMAS,N)
    !!                    AM(1,1)= ....
    !!                  IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
    !!                  AS FULL MATRIX LIKE
    !!                         AM(I,J) = M(I,J)
    !!                  ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
    !!                  DIAGONAL-WISE AS
    !!                         AM(I-J+MUMAS+1,J) = M(I,J).
    !!
    !!     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
    !!                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
    !!                       MATRIX, MAS IS NEVER CALLED.
    !!                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
    !!
    !!    MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
    !!                       MLMAS=N: THE FULL MATRIX CASE. THE LINEAR ALGEBRA
    !!                                IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
    !!                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE MATRIX
    !!                                (>= NUMBER OF NON-ZERO DIAGONALS BELOW
    !!                                THE MAIN DIAGONAL).
    !!                 MLMAS IS SUPPOSED TO BE .LE. MLJAC.
    !!
    !!     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON-
    !!                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
    !!                 NEED NOT BE DEFINED IF MLMAS=N.
    !!                 MUMAS IS SUPPOSED TO BE .LE. MUJAC.
    !!
    !!     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
    !!                 NUMERICAL SOLUTION DURING INTEGRATION.
    !!                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
    !!                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
    !!                 IT MUST HAVE THE FORM
    !!                    SUBROUTINE SOLOUT(R,TP,YP,FF,NT,DBLK,ORD,RPAR,IPAR,IRTRN)
    !!                    INTEGER R, DBLK, ORD, IPAR(*), IRTRN, NT
    !!                    DOUBLE PRECISION TP(*), YP(R,*), RPAR(*), FF(R,*)
    !!                    ....
    !!                 SOLOUT FURNISHES THE SOLUTION "YP" AT THE
    !!                    GRID-POINTS "TP(*)".
    !!                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
    !!                    IS SET <0, GAM  RETURNS TO THE CALLING PROGRAM.
    !!
    !!                 CONTINUOUS OUTPUT:
    !!                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
    !!                 FOR THE INTERVAL [TP(1),TP(DBLK+1)] IS AVAILABLE THROUGH
    !!                 THE FUNCTION
    !!                        >>>   CONTR(I,R,T,TP,FF,DBLK,NT)   <<<
    !!                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
    !!                 COMPONENT OF THE SOLUTION AT THE POINT T. THE VALUE
    !!                 T SHOULD LIE IN THE INTERVAL [TP(1),TP(DBLK+1)] ON
    !!                 WHICH THE SOLUTION IS COMPUTED AT CURRENT STEP.
    !!                 DO NOT CHANGE THE ENTRIES OF FF(R,*) and NT, IF THE
    !!                 DENSE OUTPUT FUNCTION IS USED.
    !!
    !!     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
    !!                    IOUT=0: SUBROUTINE IS NEVER CALLED
    !!                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
    !!
    !!     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
    !!                 WORK(1), WORK(2),.., WORK(LWORK) SERVE AS PARAMETERS
    !!                 FOR THE CODE. FOR STANDARD USE OF THE CODE
    !!                 WORK(1),..,WORK(LWORK) MUST BE SET TO ZERO BEFORE
    !!                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE.
    !!
    !!     LWORK       DECLARED LENGTH OF ARRAY "WORK". IN THIS VERSION SET
    !!                 LWORK=21
    !!
    !!     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
    !!                 IWORK(1),IWORK(2),...,IWORK(LIWORK) SERVE
    !!                 AS PARAMETERS FOR THE CODE. FOR STANDARD USE,
    !!                 SET IWORK(1),..,IWORK(9) TO ZERO BEFORE CALLING.
    !!                 - IWORK(10),...,IWORK(24) SERVE AS WORKING AREA.
    !!                 - IWORK(25),IWORK(26),IWORK(27) CONTAIN THE DIMENSION OF
    !!                   THE INDEX 1, INDEX2, INDEX3, VARIABLES RESPECTIVELY.
    !!                   THEY MUST BE PASSED IN GAMD AS INPUT VARIABLES WHEN
    !!                   SOLVING A DAE; THEY MAY BE SET EQUAL TO ZERO WHEN
    !!                   SOLVING AN ODE.
    !!
    !!
    !!     LIWORK      DECLARED LENGTH OF ARRAY "IWORK". IN THIS VERSION SET
    !!                 LIWORK=27
    !!
    !!     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
    !!                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
    !!                 PROGRAM AND THE FCN, JAC SUBROUTINES.
    !!
    !! -------------------------------------------------------------------------
    !!     SOPHISTICATED SETTING OF PARAMETERS
    !! -------------------------------------------------------------------------
    !!              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
    !!              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),...
    !!              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO.
    !!              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
    !!
    !!    IWORK(1)  NOT USED
    !!
    !!    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
    !!              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
    !!
    !!    IWORK(3)  ORDMIN, 3 <= ORDMIN <= 9,
    !!
    !!    IWORK(4)  ORDMAX, ORDMIN <= ORDMAX <= 9
    !!
    !!    IWORK(5)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATIONS FOR THE
    !!              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP FOR ORDER 3.
    !!              THE DEFAULT VALUE (FOR IWORK(5)=0) IS 10.
    !!
    !!    IWORK(6)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
    !!              ORDER 5, THE DEFAULT VALUE (FOR IWORK(6)=0) IS 18.
    !!
    !!    IWORK(7)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
    !!              ORDER 7, THE DEFAULT VALUE (FOR IWORK(7)=0) IS 26.
    !!
    !!    IWORK(8)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
    !!              ORDER 9, THE DEFAULT VALUE (FOR IWORK(5)=0) IS 36.
    !!
    !!    IWORK(9)  NOT YET USED
    !!
    !!    IWORK(10:24) USED AS OUTPUT PARAMETERS (SEE BELOW)
    !!
    !!       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR
    !!       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1.
    !!       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT
    !!       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER.
    !!       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE
    !!       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2.
    !!
    !!    IWORK(25)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR
    !!               ODE'S THIS EQUALS THE DIMENSION OF THE SYSTEM.
    !!               DEFAULT IWORK(25)=N.
    !!
    !!    IWORK(26)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(26)=0.
    !!
    !!    IWORK(27)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(27)=0.
    !!
    !!
    !!    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
    !!
    !!    WORK(2)   HMAX  MAXIMAL STEP SIZE, DEFAULT TEND-T0.
    !!
    !!    WORK(3)   THET DECIDE WHETHER THE JACOBIAN SHOULD BE RECOMPUTED
    !!
    !!    WORK(4)   FACNEWT:  stopping criterion for splitting-Newton method
    !!                 for small values of min(abs(y_i)) and min(abs(f_j)).
    !!
    !!    WORK(5)   TETAK0(1) stopping criterium for the order 3
    !!              splitting-Newton method:
    !!              the iterates must be decreasing by a factor tetak0(1)
    !!
    !!    WORK(6)   TETAK0(2) stopping criterium for the order 5
    !!              splitting-Newton method:
    !!              the iterates must be decreasing by a factor tetak0(2)
    !!
    !!    WORK(7)   TETAK0(3) stopping criterium for the order 7
    !!              splitting-Newton method:
    !!              the iterates must be decreasing by a factor tetak0(3)
    !!
    !!    WORK(8)   TETAK0(4) stopping criterium for the order 9
    !!              splitting-Newton method:
    !!              the iterates must be decreasing by a factor tetak0(4)
    !!
    !!    WORK(9)   CS(2): EMPIRICAL COMPUTATIONAL COST FOR ORDER  5 METHOD
    !!              USED IN THE ORDER VARIATION STRATEGY
    !!              (DEFAULT WORK(6) = 2.4D0)
    !!
    !!    WORK(10)  CS(3): EMPIRICAL COMPUTATIONAL COST FOR ORDER  7 METHOD
    !!              USED IN THE ORDER VARIATION STRATEGY
    !!              (DEFAULT WORK(6) = 4.0D0)
    !!
    !!    WORK(11)  CS(4): EMPIRICAL COMPUTATIONAL COST FOR ORDER  9 METHOD
    !!              USED IN THE ORDER VARIATION STRATEGY
    !!              (DEFAULT WORK(6) = 7.2D0)
    !!
    !!    WORK(12)-WORK(13)   FACL-FACR: PARAMETERS FOR STEP SIZE SELECTION
    !!               THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
    !!               FACL  <=  HNEW/HOLD <= FACR
    !!               (DEFAULT WORK(9) = 0.12, WORK(10) = 10 )
    !!
    !!    WORK(14)  SFDOWN:SAFETY FACTOR IN STEP SIZE PREDICTION
    !!                  USED FOR THE LOWER ORDER METHOD
    !!                  (DEFAULT WORK(11) = 20D0)
    !!
    !!    WORK(15)  SFUP:SAFETY FACTOR IN STEP SIZE PREDICTION
    !!                  USED FOR THE UPPER ORDER METHOD
    !!                  (DEFAULT WORK(12) = 40D0)
    !!
    !!    WORK(16)  SFSAME: SAFETY FACTOR IN STEP SIZE PREDICTION
    !!                  USED FOR THE CURRENT ORDER METHOD
    !!                  (DEFAULT WORK(13) = 18D0)
    !!
    !!    WORK(17)  SF: SAFETY FACTOR IN STEP SIZE PREDICTION
    !!                  USED FOR THE CURRENT ORDER METHOD WHEN IS
    !!                  FAILED THE ERROR CONTROL TEST (DEFAULT WORK(14) = 15D0)
    !!
    !!
    !!
    !!    WORK(18)  FACNEWT stopping criterion for splitting-Newton method ORDER 3
    !!                  (DEFAULT WORK(15) = 1.0D-3)
    !!
    !!    WORK(19)  FACNEWT stopping criterion for splitting-Newton method ORDER 5
    !!                  (DEFAULT WORK(16) = 9.0D-2)
    !!
    !!    WORK(20)  FACNEWT stopping criterion for splitting-Newton method ORDER 7
    !!                  (DEFAULT WORK(17) = 9.0D-1)
    !!
    !!    WORK(21)  FACNEWT stopping criterion for splitting-Newton method ORDER 9
    !!                  (DEFAULT WORK(18) = 9.9D-1)
    !!
    !! -------------------------------------------------------------------------
    !!     OUTPUT PARAMETERS
    !! -------------------------------------------------------------------------
    !!     T0          T-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
    !!                 (AFTER SUCCESSFUL RETURN T0=TEND).
    !!
    !!     Y(N)        NUMERICAL SOLUTION AT T0
    !!
    !!     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
    !!
    !!     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
    !!                   IDID= 1  COMPUTATION SUCCESSFUL,
    !!                   IDID=-1  INPUT IS NOT CONSISTENT,
    !!                   IDID=-2  LARGER NMAX IS NEEDED,
    !!                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
    !!                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
    !!                   IDID=-5  GAM CANNOT HANDLES IERR=-1.
    !!
    !!   IWORK(10)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
    !!                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)
    !!   IWORK(11)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
    !!                      OR NUMERICALLY)
    !!   IWORK(12)  NSTEP(1)  NUMBER OF COMPUTED STEPS   ORD 3
    !!   IWORK(13)  NSTEP(2)  NUMBER OF COMPUTED STEPS   ORD 5
    !!   IWORK(14)  NSTEP(3)  NUMBER OF COMPUTED STEPS   ORD 7
    !!   IWORK(15)  NSTEP(4)  NUMBER OF COMPUTED STEPS   ORD 9
    !!   IWORK(16)  NNEWT(1)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 3
    !!   IWORK(17)  NNEWT(2)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 5
    !!   IWORK(18)  NNEWT(3)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 7
    !!   IWORK(19)  NNEWT(4)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 9
    !!   IWORK(20)  NERR(1)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 3
    !!   IWORK(21)  NERR(2)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 5
    !!   IWORK(22)  NERR(3)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 7
    !!   IWORK(23)  NERR(4)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 9
    !!   IWORK(24)  NDEC      NUMBER OF LU-DECOMPOSITIONS
    !!-----------------------------------------------------------------------
    !!     DECLARATIONS
    !! -------------------------------------------------------------------------
    USE PRECISION
    IMPLICIT NONE
    !!
    !!   INPUT VARIABLES
    !!------------------------------------
    INTEGER, INTENT(IN) :: R, ITOL, IJAC, IMAS, IOUT, IPAR(*), LIWORK, LWORK

    REAL(PREC), INTENT(IN) :: TEND, ATOL(*), RTOL(*), RPAR(*)
    !!
    !!   INPUT/OUTPUT VARIABLES
    !!------------------------------------
    INTEGER, INTENT(IN OUT)    ::  IWORK(LIWORK), IDID,  &
         &  MLJAC, MUJAC, MLMAS, MUMAS
    REAL(PREC), INTENT(IN OUT) :: T0, Y0(R), H, WORK(LWORK)

    !!
    !!   LOCAL VARIABLES
    !!------------------------------------
    REAL(PREC) :: FACNORD(4), HMAX, THET, FACNEWT, TETAK0(4), CS(4), FACL, FACR,       &
         &                  SFDOWN, SFUP, SFSAME, SF, UROUND
    INTEGER :: ORDMIN, ORDMAX, ITINT(4), ITMAX, IJOB, NMAX, LDJAC, LDLU, LDMAS
    INTEGER ::  NDEC, NFCN, NJAC, NSTEP(4), NNEWT(4), NERR(4)
    INTEGER ::  IEYP, IEFP, IEDN,IEF,IEF1, IEJF0,IELU, ISTORE, &
         &   I,IEIPIV, IESC, NIND1, NIND2, NIND3
    LOGICAL ::  ARRET, JBAND, IMPLCT

    !!
    !!   EXTERNAL FUNCTIONS
    !!------------------------------------
    !!EXTERNAL FCN, JAC, MAS, SOLOUT
    INTERFACE
       !!-----------------------------------------------------------------------
       SUBROUTINE fcn(neqn,t,y,dy,rpar,ipar)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: neqn, ipar(*)
         REAL(PREC), INTENT(IN) :: t,y(neqn),rpar(*)
         REAL(PREC), INTENT(OUT) :: dy(neqn)
       END SUBROUTINE fcn
       !!-----------------------------------------------------------------------
       SUBROUTINE jac(neqn,t,y,mu,ml,jacob,ldim,rpar,ipar)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: neqn,ldim,ipar(*),mu,ml
         REAL(PREC), INTENT(IN) :: t,y(neqn),rpar(*)
         REAL(PREC), INTENT(OUT) :: jacob(ldim,neqn)
       END SUBROUTINE jac
       !!-----------------------------------------------------------------------
       SUBROUTINE mas(neqn,am,ldim,rpar,ipar)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: neqn, ldim, ipar(*)
         REAL(PREC), INTENT(IN) :: rpar(*)
         REAL(PREC), INTENT(OUT) :: am(ldim,neqn)
        END SUBROUTINE mas
       !!-----------------------------------------------------------------------
       SUBROUTINE SOLOUT(R,TP,YP,F1,NT1,DBLK,ORD,RPAR,IPAR,IRTRN)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: R, DBLK, ORD, IPAR(*),NT1
         INTEGER, INTENT(OUT) :: IRTRN
         REAL(PREC), INTENT(IN) :: TP(*),RPAR(*),F1(R,*)
         REAL(PREC), INTENT(IN) :: YP(R,*)
       END SUBROUTINE solout
    END INTERFACE
    !! -------------------------------------------------------------------------
    !!     SETTING THE PARAMETERS
    !! -------------------------------------------------------------------------


    ARRET   = .FALSE.
    !! -------- NMAX := THE MAXIMAL NUMBER OF STEPS -----
    IF (IWORK(2).EQ.0) THEN
       NMAX=100000
    ELSE
       NMAX=IWORK(2)
       IF (NMAX.LE.0) THEN
      CALL Rprinti1('Wrong input iwork(2) = ',IWORK(2))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- ORDMIN  :=  MINIMAL ORDER
    IF (IWORK(3).EQ.0) THEN
       ORDMIN = 1
    ELSE
       ORDMIN=(IWORK(3)-1)/2
        IF ((ORDMIN.LE.0).OR.(ORDMIN.GT.4)) THEN
      CALL Rprinti1('Curious input iwork(3) = ',IWORK(3))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- ORDMAX :=  MAXIMAL ORDER
    IF (IWORK(4).EQ.0) THEN
       ORDMAX = 4
    ELSE
       ORDMAX=(IWORK(4)-1)/2
       IF ((ORDMAX.LE.0).OR.(ORDMAX.GT.4).OR.(ORDMAX.LT.ORDMIN)) THEN
      CALL Rprinti1('Curious input iwork(4) = ',IWORK(4))
          ARRET=.TRUE.
       END IF
    END IF
    !! -------- ITINT(1) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 3
    IF (IWORK(5).EQ.0) THEN
       ITINT(1)=12
    ELSE
       ITINT(1)=IWORK(5)
       IF (ITINT(1).LE.0) THEN
      CALL Rprinti1('Curious input iwork(5) = ',IWORK(5))
          ARRET=.TRUE.
       END IF
    END IF
    ITMAX = ITINT(1)
    !! -------- ITINT(2) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 5
    IF (IWORK(6).EQ.0) THEN
       ITINT(2)=18
    ELSE
       ITINT(2)=IWORK(6)
       IF (ITINT(2).LT.0) THEN
      CALL Rprinti1('Curious input iwork(6) = ',IWORK(6))
          ARRET=.TRUE.
       END IF
    END IF
    !! -------- ITINT(3) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 7
    IF (IWORK(7).EQ.0) THEN
       ITINT(3)= 26
    ELSE
       ITINT(3)=IWORK(7)
       IF (ITINT(3).LT.0) THEN
      CALL Rprinti1('Curious input iwork(7) = ',IWORK(7))
          ARRET=.TRUE.
       END IF
    END IF
    !! -------- PARAMETER FOR DIFFERENTIAL-ALGEBRAIC COMPONENTS
      NIND1=IWORK(25)
      NIND2=IWORK(26)
      NIND3=IWORK(27)

      IF (NIND1.EQ.0) NIND1=R
      IF (NIND1+NIND2+NIND3.NE.R) THEN
      CALL Rprinti3('Curious input iwork(25,26,27) ',NIND1,NIND2,NIND3)
         ARRET=.TRUE.
      END IF

    !! -------- ITINT(4) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 9
    IF (IWORK(8).EQ.0) THEN
       ITINT(4)= 36
    ELSE
       ITINT(4)=IWORK(8)
       IF (ITINT(4).LT.0) THEN
      CALL Rprinti1('Curious input iwork(8) = ',IWORK(8))
          ARRET=.TRUE.
       END IF
    END IF

    !! -------- UROUND :=  SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0
    IF (WORK(1).EQ.0.0D0) THEN
       UROUND=2.30D-16
    ELSE
       UROUND=WORK(1)
       IF (UROUND.LE.1.0D-19.OR.UROUND.GE.1.0D0) THEN
      CALL Rprintd1('Coefficients have 20 digits, Uround = ',WORK(1))
          ARRET=.TRUE.
       END IF
    END IF
    !! -------- HMAX := MAXIMAL STEP SIZE
    IF (WORK(2).EQ.0.D0) THEN
       HMAX=TEND-T0
    ELSE
       HMAX=WORK(2)
       IF (HMAX.GT.TEND-T0) THEN
          HMAX=TEND-T0
       END IF
    END IF
    !! -------- THET:  DECIDE WHETHER THE JACOBIAN SHOULD BE RECOMPUTED
    IF (WORK(3).EQ.0.D0) THEN
       THET = 0.005
    ELSE
       THET=WORK(3)
       IF (THET .GT. 1d0) THEN
      CALL Rprintd1('Curious input work(3) = ',WORK(3))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACNEWT: STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
    !!--------           FOR SMALL VALUES OF min(abs(y_i)) and min(abs(f_j))
    IF (WORK(4).EQ.0.D0) THEN
       FACNEWT= 1d-2
       FACNEWT=MAX(FACNEWT,EPS/RTOL(1) )
    ELSE
       FACNEWT=WORK(4)
       FACNEWT=MAX(FACNEWT,EPS/RTOL(1) )
       IF (FACNEWT.GE.1.0D0) THEN
      CALL Rprintd1('Wrong input for work(4) ',WORK(4))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACNORD(1): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
    !!---------             ORDER 3
    IF (WORK(15).EQ.0.D0) THEN
       SELECT CASE (IMAS)
       CASE(0)    ! ODE CASE
          FACNORD(1) = 1D-3
       CASE(1)    ! DAE CASE
          FACNORD(1) = 5d-3
       END SELECT
       FACNORD(1) = MAX(FACNORD(1) ,EPS/RTOL(1) )
    ELSE
       FACNORD(1) = WORK(15)
       FACNORD(1) = MAX(FACNORD(1) ,EPS/RTOL(1) )
       IF (FACNEWT.GE.1.0D0) THEN
      CALL Rprintd1('Wrong input for work(15) ',WORK(15))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACNORD(2): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
    !!--------              ORDER 5
    IF (WORK(16).EQ.0.D0) THEN
       SELECT CASE (IMAS)
       CASE(0)    ! ODE CASE
          FACNORD(2) = 5d-2
       CASE(1)    ! DAE CASE
          FACNORD(2) = 5D-2
       END SELECT
       FACNORD(2) = MAX(FACNORD(2) ,EPS/RTOL(1) )
    ELSE
       FACNORD(2) = WORK(16)
       FACNORD(2) = MAX(FACNORD(2) ,EPS/RTOL(1) )
       IF (FACNORD(2).GE.1.0D0) THEN
      CALL Rprintd1('Wrong input for work(16) ',WORK(16))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACNORD(3): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
    !!---------             ORDER 7
    IF (WORK(17).EQ.0.D0) THEN
       SELECT CASE (IMAS)
       CASE(0)    ! ODE CASE
          FACNORD(3) = 1.0d-1
       CASE(1)    ! DAE CASE
          FACNORD(3) = 5.0d-2
       END SELECT
       FACNORD(3) = MAX(FACNORD(3), EPS/RTOL(1) )
    ELSE
       FACNORD(3) = WORK(17)
       FACNORD(3) = MAX(FACNORD(3), EPS/RTOL(1) )
       IF (FACNORD(3).GE.1.0D0) THEN
      CALL Rprintd1('Wrong input for work(17) ',WORK(17))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACNORD(4): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
    !!---------             ORDER 9
    IF (WORK(18).EQ.0.D0) THEN
       SELECT CASE (IMAS)
       CASE(0)    ! ODE CASE
          FACNORD(4) = 1.0d-1
       CASE(1)    ! DAE CASE
          FACNORD(4) = 5.0d-2
       END SELECT
       FACNORD(4) = MAX(FACNORD(4),EPS/RTOL(1) )
    ELSE
       FACNORD(4) = WORK(18)
       FACNORD(4) = MAX(FACNORD(4),EPS/RTOL(1) )
       IF (FACNORD(4).GE.1.0D0) THEN
      CALL Rprintd1('Wrong input for work(18) ',WORK(18))
          ARRET=.TRUE.
       END IF
    END IF

    !!--------- TETAK0(1:4): STOPPING CRITERIUM FOR THE SPLITTING-NEWTON METHOD
    !!---------              THE ERROR IN THE ITERATES MUST BE DECREASING BY A FACTOR
    !!---------              TETAK0(I), i=1,..,4 ACCORDING TO THE SELECTED ORDER
    IF (IMAS==0) THEN
       TETAK0(1:4) = 0.9D0                         ! ODE CASE
    ELSE
       TETAK0(1:4) = (/ 0.9D0, 1.5d0, 20D0, 20D0 /)    ! DAE CASE
    END IF
    DO I=1,4
       IF (WORK(4+I).NE.0.D0) THEN
          TETAK0(I) = WORK(4+I)
       END IF
    END DO
    IF (MINVAL(TETAK0).LE.0.0D0) THEN
      CALL Rprintd4('Wrong input for work(5:8) ',WORK(5:8))
       ARRET=.TRUE.
    END IF

    CS(1) = 1.0D0
    !!--------- CS(2): EMPIRICAL COMPUTATIONAL COST FOR ORDER 5
    !!---------        USED IN THE ORDER VARIATION STRATEGY.
    IF (WORK(9).EQ.0.D0) THEN
       CS(2) = 2.4D0
    !   IF ((NIND2.NE.0).OR. (NIND3.NE.0)) CS(2) = 2.8d0;
    ELSE
       CS(2) = WORK(9)
       IF (CS(2).LE.0.0D0) THEN
      CALL Rprintd1('Wrong input for work(9) ',WORK(9))
          ARRET=.TRUE.
       END IF
    END IF
    !!---------  CS(3): EMPIRICAL COMPUTATIONAL COST FOR ORDER 7
    !!---------         USED IN THE ORDER VARIATION STRATEGY.
    IF (WORK(10).EQ.0.D0) THEN
       CS(3) = 4.0D0
    !   IF  ((NIND2.NE.0).OR. (NIND3.NE.0)) CS(2) = 4.5d0;
    ELSE
       CS(3) = WORK(10)
       IF (CS(3).LE.0.0D0) THEN
      CALL Rprintd1('Wrong input for work(10) ',WORK(10))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- CS(4): EMPIRICAL COMPUTATIONAL COST FOR ORDER 9
    !!---------        USED IN THE ORDER VARIATION STRATEGY.
    IF (WORK(11).EQ.0.D0) THEN
       CS(4) =7.2D0
    !  IF ((NIND2.NE.0).OR. (NIND3.NE.0)) CS(2) = 8d0;
    ELSE
       CS(4) = WORK(11)
       IF (CS(4).LE.0.0D0) THEN
      CALL Rprintd1('Wrong input for work(11) ',WORK(11))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACL: PARAMETER FOR STEP SIZE SELECTION
    !!---------       THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
    !!---------       FACL <= HNEW/HOLD
    IF (WORK(12).EQ.0.D0) THEN
       FACL = 0.12D0
    ELSE
       FACL = WORK(12)
       IF (FACL.LE.0.0D0) THEN
      CALL Rprintd1('Wrong input for work(12) ',WORK(12))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACR: PARAMETER FOR STEP SIZE SELECTION
    !!---------       THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
    !!---------       HNEW/HOLD <= FACR
    IF (WORK(13).EQ.0.D0) THEN
       FACR = 10D0
    ELSE
       FACR = WORK(13)
       IF (FACR.LE.0.0D0) THEN
      CALL Rprintd1('Wrong input for work(13) ',WORK(13))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- SFDOWN: SAFETY FACTOR IN STEP SIZE PREDICTION
    !!---------         USED FOR THE LOWER ORDER METHOD
    IF (WORK(14).EQ.0.D0) THEN
       SFDOWN = 20.0D0
    ELSE
       SFDOWN = WORK(14)
       IF (SFDOWN.LE.0.0D0) THEN
      CALL Rprintd1('Wrong input for work(14) ',WORK(14))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- SFUP:  SAFETY FACTOR IN STEP SIZE PREDICTION
    !!---------        USED FOR THE UPPER ORDER METHOD
    IF (WORK(15).EQ.0.D0) THEN
       SFUP = 40.0D0
    ELSE
       SFUP = WORK(15)
       IF (SFUP.LE.0.0D0) THEN
      CALL Rprintd1('Wrong input for work(15) ',WORK(15))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- SFSAME: SAFETY FACTOR IN STEP SIZE PREDICTION
    !!---------         USED FOR THE CURRENT ORDER METHOD
    IF (WORK(16).EQ.0.D0) THEN
       SFSAME = 18.0D0
    ELSE
       SFSAME = WORK(16)
       IF (SFSAME.LE.0.0D0) THEN
      CALL Rprintd1('Wrong input for work(16) ',WORK(16))
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- SF: SAFETY FACTOR IN STEP SIZE PREDICTION
    !!---------     USED FOR THE CURRENT ORDER METHOD WHEN
    !!---------     THE ERROR CONTROL TEST fails
    IF (WORK(17).EQ.0.D0) THEN
       SF = 15.0D0
    ELSE
       SF = WORK(17)
       IF (SF.LE.0.0D0) THEN
      CALL Rprintd1('Wrong input for work(17) ',WORK(17))
          ARRET=.TRUE.
       END IF
    END IF
    !! -------- CHECK  THE TOLERANCES
    IF (ITOL.EQ.0) THEN
       IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE. EPS) THEN

      CALL Rprint('Tolerances are too small')
      
          ARRET=.TRUE.
       END IF
    ELSE
       DO I=1,R
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE. EPS) THEN
      CALL Rprint('Tolerances(',I,') are too small')
             ARRET=.TRUE.
          END IF
       END DO
    END IF
    !! -------------------------------------------------------------------------
    !!     COMPUTATION OF ARRAY ENTRIES
    !! -------------------------------------------------------------------------
    !! -------- BANDED OR NOT, IMPLICIT ?
    JBAND=MLJAC.LT.R
    IMPLCT=IMAS.NE.0
    !! -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS
    IF (JBAND) THEN
       LDJAC = MLJAC+MUJAC+1
       LDLU  = MLJAC+LDJAC
    ELSE
       LDJAC = R
       LDLU  = R
    END IF
!    IF (JBAND) THEN
!       IJOB=2
!    ELSE
!       IJOB=1
!    END IF
    !! --------  MASS MATRIX
      IF (IMPLCT) THEN
          IF (MLMAS.NE.R) THEN
              LDMAS=MLMAS+MUMAS+1
              IF (JBAND) THEN
                 IJOB=4
              ELSE
                 IJOB=3
              END IF
          ELSE
              MUMAS=R
              LDMAS=R
              IJOB=5
          END IF
    !! ------   BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF "JAC"
          IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN
      CALL Rprint('Bandwidth of MAS not smaller than bandwidth of JAC')
            ARRET=.TRUE.
          END IF
      ELSE
          LDMAS=0
          IF (JBAND) THEN
             IJOB=2
          ELSE
             IJOB=1
          END IF
      END IF
      LDMAS=MAX(1,LDMAS)


!! -------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK_IN
!! -------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK
!!      IEYP  = 19
!!      IEFP  = IEYP  + 10*R
!!      IEDN  = IEFP  + 10*R
!!      IEF   = IEDN  + R
!!      IEF1  = IEF   + 10*R
!!      IESC  = IEF1  + 10*R
!!      IEJF0 = IESC  + R
!!      IELU  = IEJF0 + R*LDJAC
!!--------- TOTAL STORAGE REQUIREMENT
!!      ISTORE = IELU + R*LDLU - 1
      ISTORE = 21
      IF(ISTORE.GT.LWORK)THEN
      CALL Rprinti1('Insufficient storage for work, min. is ',ISTORE)
         ARRET=.TRUE.
      END IF
!! -------- ENTRY POINTS FOR INTEGER WORKSPACE
!!      IEIPIV=25
!! -------- TOTAL REQUIREMENT
!!      ISTORE=IEIPIV+R-1
      ISTORE = 27
      IF (ISTORE.GT.LIWORK) THEN
      CALL Rprinti1('Insuff. storage for iwork, min.  = ',ISTORE)
         ARRET=.TRUE.
      END IF
    !! -------- WHEN A FAIL HAS OCCURED, GAM RETURNs WITH IDID=-1
    IF (ARRET) THEN
       IDID=-1
       RETURN
    END IF


    !! -------------------------------------------------------------------------
    !!     CALL TO CORE INTEGRATOR
    !! -------------------------------------------------------------------------
    CALL ETRO(R, FCN, T0, Y0, TEND, HMAX, H, RTOL, ATOL, ITOL,                       &
         &   JAC, IJAC, MLJAC, MUJAC, MAS, MLMAS, MUMAS, SOLOUT, IOUT, IDID, NMAX,   &
         &   UROUND, THET, FACNEWT, FACNORD, TETAK0, CS, FACL, FACR, SFDOWN,         &
         &   SFUP, SFSAME, SF, ORDMIN, ORDMAX, ITINT, ITMAX,                         &
         &   IMPLCT, JBAND, IJOB, LDJAC, LDLU, LDMAS, NIND1, NIND2, NIND3,           &
         &   NFCN, NJAC, NSTEP, NNEWT, NERR, NDEC, RPAR, IPAR)

    IWORK(10)= NFCN
    IWORK(11)= NJAC
    IWORK(12)= NSTEP(1)
    IWORK(13)= NSTEP(2)
    IWORK(14)= NSTEP(3)
    IWORK(15)= NSTEP(4)
    IWORK(16)= NNEWT(1)
    IWORK(17)= NNEWT(2)
    IWORK(18)= NNEWT(3)
    IWORK(19)= NNEWT(4)
    IWORK(20)= NERR(1)
    IWORK(21)= NERR(2)
    IWORK(22)= NERR(3)
    IWORK(23)= NERR(4)
    IWORK(24)= NDEC

    RETURN
  CONTAINS
    !!END SUBROUTINE GAM
    !!
    !!--------- END OF SUBROUTINE GAM
    !!
    !! -------------------------------------------------------------------------
    !!     SUBROUTINE  ETRO (Extended trapezoidal Rules of Odd order,
    !!                       that is GAMs)
    !! -------------------------------------------------------------------------
    SUBROUTINE  ETRO(R,FCN,T0,Y0,TEND,HMAX,H,RTOL,ATOL,ITOL,                  &
         &   JAC,IJAC,MLJAC,MUJAC,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,           &
         &   NMAX,UROUND,THET,FACNEWT,FACNORD,TETAK0,CS,FACL,FACR,SFDOWN,     &
         &   SFUP,SFSAME,SF, ORDMIN,ORDMAX,ITINT,ITMAX,                       &
         &   IMPLCT,JBAND,IJOB,LDJAC,LDLU,LDMAS,NIND1,NIND2,NIND3,            &
         &   NFCN,NJAC,NSTEP,NNEWT,NERR,NDEC,RPAR,IPAR)
      !! -------------------------------------------------------------------------
      !!     CORE INTEGRATOR FOR GAM
      !!     PARAMETERS SAME AS IN GAM WITH ADDED WORKSPACE
      !! -------------------------------------------------------------------------
      !!     DECLARATIONS
      !! ----------------------------------------------------------

      USE LINALGGAMD
      USE SUBGAMD
      IMPLICIT NONE
      !!
      !!   COMMON
      !!------------------------------------
      !!COMMON/LINAL/MLLU,MULU,MDIAG
      !!
      !!   INPUT VARIABLES
      !!------------------------------------
      INTEGER, INTENT(IN) ::  R, ORDMIN, ORDMAX,  ITOL, IJAC, MLJAC, MUJAC, MLMAS, MUMAS, IOUT, &
           &         IPAR(*), ITINT(4),  IJOB, NMAX, LDJAC, LDLU, LDMAS, NIND1, NIND2, NIND3

      REAL(PREC), INTENT(IN)  :: TEND, ATOL(*), RTOL(*), RPAR(*), FACNORD(4),&
           &                  HMAX, THET, FACNEWT, TETAK0(4), CS(4), FACL, FACR,       &
           &                  SFDOWN, SFUP, SFSAME, SF, UROUND


      LOGICAL, INTENT(IN) :: IMPLCT
      !!
      !!   OUTPUT VARIABLES
      !!------------------------------------
      INTEGER, INTENT(OUT) :: NDEC, NFCN, NJAC, NSTEP(4), NNEWT(4), NERR(4)
      !!
      !!   INPUT/OUTPUT VARIABLES
      !!----
      REAL(PREC), INTENT(IN OUT) :: T0, Y0(R), H

      INTEGER, INTENT(IN OUT) ::  IDID, ITMAX
      !!
      !!   LOCAL VARIABLES
      !!------------------------------------
      REAL(PREC), ALLOCATABLE :: SCAL(:), YP(:,:), FP(:,:),        &
           &       F(:,:),DN(:), F1(:,:), JF0(:,:), LU(:,:), FMAS(:,:)

      INTEGER, ALLOCATABLE :: IPIV(:)



      INTEGER :: I, J, NSING,                                      &
           &          FAILNI, FAILEI, ORDOLD, ORD, ORD2, ORDN,     &
           &          IT, DBL(4), DBLK, DBLKOLD, ORDDEC,           &
           &          NSTEPS, IRTRN, NT1, IER

      REAL(PREC) :: ERRV(10), TP(11), T1(11), YSAFE, DELT,         &
           &                  THETA, TETAK, TETAKOLD,THETAPREC,    &
           &                  HOLD, HDEC,  ERRNEWT, ERRNEWT0,ESP,  &
           &                  ERRUP, ERRSAME, ERRDOWN, RR, RRN, TH, THN, FACN
      REAL(PREC) ::    DYTH,QNEWT, HHFAC, HACC, R0, HEXTRAP, CPORD(4), ERRSAMEOLD

      INTEGER ::  IFAC,  IERR
      LOGICAL :: JBAND, CALJAC, NEWJAC, JVAI, TER, EXTRAP
      LOGICAL :: VAR, INDEX1, INDEX2, INDEX3
      !!
      !!   EXTERNAL FUNCTIONS
      !!------------------------------------
      !!EXTERNAL FCN,JAC,MAS,SOLOUT
       INTERFACE
       !!-----------------------------------------------------------------------
       SUBROUTINE fcn(neqn,t,y,dy,rpar,ipar)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: neqn, ipar(*)
         REAL(PREC), INTENT(IN) :: t,y(neqn),rpar(*)
         REAL(PREC), INTENT(OUT) :: dy(neqn)
       END SUBROUTINE fcn
       !!-----------------------------------------------------------------------
       SUBROUTINE jac(neqn,t,y,mu,ml,jacob,ldim,rpar,ipar)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: neqn,ldim,ipar(*),mu,ml
         REAL(PREC), INTENT(IN) :: t,y(neqn),rpar(*)
         REAL(PREC), INTENT(OUT) :: jacob(ldim,neqn)
       END SUBROUTINE jac
       !!-----------------------------------------------------------------------
       SUBROUTINE mas(neqn,am,ldim,rpar,ipar)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: neqn, ldim, ipar(*)
         REAL(PREC), INTENT(IN) :: rpar(*)
         REAL(PREC), INTENT(OUT) :: am(ldim,neqn)
        END SUBROUTINE mas
       !!-----------------------------------------------------------------------
       SUBROUTINE SOLOUT(R,TP,YP,F1,NT1,DBLK,ORD,RPAR,IPAR,IRTRN)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: R, DBLK, ORD, IPAR(*),NT1
         INTEGER, INTENT(OUT) :: IRTRN
         REAL(PREC), INTENT(IN) :: TP(*),RPAR(*),F1(R,*)
         REAL(PREC), INTENT(IN) :: YP(R,*)
       END SUBROUTINE solout
    END INTERFACE
      !! -------- CONSTANTS


      ALLOCATE( SCAL(R), YP(R,11), FP(R,11),                            &
           &       F(R,11), DN(R), F1(R,11), JF0(LDJAC,R), LU(LDLU,R),  &
           &       FMAS(LDMAS,R), IPIV(R)  )


      MLLU=MLJAC
      MULU=MUJAC
      MDIAG=MLLU + MULU +1
      !! ------- CHECK THE INDEX OF THE PROBLEM -----
      INDEX1=NIND1.NE.0
      INDEX2=NIND2.NE.0
      INDEX3=NIND3.NE.0
      !! ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------

      IF (IMPLCT) THEN
        CALL MAS(R,FMAS,LDMAS,RPAR,IPAR)
        MBDIAG=MUMAS+1
        MBB=MLMAS+MUMAS+1
        MDIFF=MLJAC+MUJAC-MUMAS
      END IF


      !!--------- DBL(1:4) := SIZE OF THE COEFFICIENT MATRICES DEFINING THE GAMs

      DBL(1) = 4
      DBL(2) = 6
      DBL(3) = 8
      DBL(4) = 9

      ORD  = ORDMIN
      NSTEPS = 0
      NFCN = 0
      NJAC = 0
      NDEC = 0
      NSTEP(1:4)=0
      NNEWT(1:4)=0
      NERR(1:4)=0


      !!--------- STARTING VALUES FOR NEWTON ITERATION
      DBLK = DBL(ORD)
      DBLKOLD = DBLK
      H    = MIN( H, ABS(TEND-T0)/DBLK )
      IERR = 0
      CALL FCN(R,T0,Y0,FP(1,1), RPAR,IPAR)  ! removed ierr
      NFCN = NFCN + 1
      IF (IERR.NE.0) THEN
         IDID=-5
      CALL Rprint('GAM: ERROR: cannot handle IERR = -1 in (T0,Y0)')
          RETURN
      ENDIF

      !! SCALING FACTOR OF THE STEPSIZE FOR THE HIGHER INDEX PROBLEMS
      !! USED TO DEFINE THE SCALING FACTOR FOR THE COMPUTATION OF THE ERROR

          CPORD(1) = 1d0
          CPORD(2) = 1d-1
          CPORD(3) = 1d-2
          CPORD(4) = 1d-3



      !! -------- NUMBER OF FAILURES IN THE SPLITTING-NEWTON SCHEME
      FAILNI = 0
      !! -------- NUMBER OF FAILURES DUE TO THE ERROR TEST
      FAILEI = 0
      NSING  = 0
      ORDOLD = 2*ORD
      ORDDEC = 2*ORD
      HOLD   = 2*H
      HACC   = 2*H
      HDEC   = 2*H
      CALJAC = .TRUE.
      EXTRAP = .FALSE.
      ITMAX  = ITINT(ORD)

      CALL FCN(R,T0,Y0,DN(1),RPAR,IPAR)  ! removed ierr


      !!--------- MAIN LOOP (ADVANCING IN TIME)
      LOOP_TIME : DO
         !! 100   CONTINUE
         !!--------- (EVENTUALLY) COMPUTE THE JACOBIAN MATRIX NUMERICALLY
         NEWJAC = .FALSE.
         ERRSAME = 0d0
         IF (CALJAC) THEN
            IF (IJAC.EQ.0) THEN
               DO I=1,R
                  YSAFE=Y0(I)
                  DELT=SQRT(EPS*MAX(1.D-5,ABS(YSAFE)))
                  Y0(I)=YSAFE+DELT
                  IERR = 0
                  CALL FCN(R,T0,Y0,DN,RPAR,IPAR)  ! removed ierr
                  IF (IERR.NE.0) THEN
                    IDID=-5
       CALL Rprint('ERROR:  IERR = -1 in numerical Jacobians')
                    RETURN
                  ENDIF
                  IF (JBAND) THEN
                     DO J=MAX(1,I-MUJAC),MIN(R,I+MLJAC)
                        JF0(J-I+MUJAC+1,I) = (DN(J)-FP(J,1))/DELT
                     END DO
                  ELSE
                    JF0(1:R,I)=(DN(1:R)-FP(1:R,1))/DELT
                  END IF
                  Y0(I)=YSAFE
               END DO
            ELSE
               !! -------- COMPUTE JACOBIAN MATRIX ANALYTICALLY

               CALL JAC(R,T0,Y0(1),MLJAC,MUJAC,JF0(1,1),LDJAC,RPAR,IPAR)

            END IF
            NJAC = NJAC + 1
            NEWJAC = .TRUE.
         END IF




         !!--------- DEFINE SCAL

         THN = 1d0
         J    = 0
         IF (ITOL.EQ.0) THEN
            DO I=1,R
               RRN = ABS(Y0(I))
               SCAL(I)=ATOL(1)+RTOL(1)*RRN
               IF (RRN .LT. THN) THEN
                  J = I
                  THN = RRN
               ENDIF
            END DO


         ELSE
            DO I=1,R
               RRN = ABS( Y0(I) )
               SCAL(I)=ATOL(I)+RTOL(I)*RRN
               IF (RRN .LT. THN) THEN
                  J = I
                  THN = RRN
               ENDIF
            END DO

         END IF



         !!--------- DEFINE  FACN
         FACN = FACNORD(ORD)
         IF (THN .LT. 1d-1) THEN
            IF (ABS(FP(J,1)) .LT. 1d-5) THEN
               FACN = MIN(FACNEWT, FACNORD(ORD) )
            END IF
         END IF





         !!--------- DEFINE TP AND YP
         IER = 1
         DO WHILE (IER .NE. 0)
          IF (EXTRAP) THEN
            T1(1) = T0+H
            DO I=2,DBLK+1
               T1(I) = T1(I-1)+H
            END DO
            CALL INTERP(R,TP(1),YP(1,1),T1(1),F1(1,1),NT1,DBLKOLD,DBLK,T0,Y0(1))
          ELSE
            TP(1) = T0
            YP(1:R,1) = Y0(1:R)
            DO I=2,DBLK+1
               YP(1:R,I) = Y0(1:R)
               TP(I) = TP(I-1)+H
            END DO
          END IF

          IER = 0
          DO J = 2, DBLK+1
             IERR = 0
             CALL FCN(R,TP(J),YP(1,J),FP(1,J), RPAR,IPAR)  ! removed ierr
             IER = MIN( IER, IERR)
          END DO
          NFCN = NFCN + DBLK
          IF (IER .NE. 0) THEN
            H = MIN(HACC,H/2d0)
          END IF
         END DO
         HEXTRAP  = H

         !!--------- FACTORIZE THE ITERATION MATRIX
         IF ((ORDDEC.NE.ORD) .OR. (HDEC.NE.H).OR.(NEWJAC) ) THEN
            IER = 1
            DO WHILE ( IER .NE. 0)
               CALL DECLU(R,JF0(1,1),H,LDJAC,LU(1,1),LDLU,IPIV,FMAS,LDMAS,MLMAS,MUMAS,ORD,IER,IJOB)
               NDEC = NDEC + 1
               IF (IER.NE.0) THEN
                  NSING = NSING + 1
                  IF (NSING.GT.5) THEN
      CALL Rprinti1('Matrix is repeatedly singular, IER= ',IER)
      
      CALL Rprintd1(' exit of GAM at t = ',T0)
      
                     IDID=-4
                     GOTO 800
                  ELSE
                     H = MIN(HACC,H/2D0)
                  END IF
               END IF
            END DO
            HDEC = H
            ORDDEC = ORD
         END IF
         IF (H .NE. HEXTRAP) THEN

         !!--------- DEFINE AGAIN  TP AND YP
            IER = 1
            DO WHILE (IER .NE. 0)
            IF (EXTRAP) THEN
              T1(1) = T0+H
              DO I=2,DBLK+1
                 T1(I) = T1(I-1)+H
              END DO
              CALL INTERP(R,TP(1),YP(1,1),T1(1),F1(1,1),NT1,DBLKOLD,DBLK,T0,Y0(1))
            ELSE
              TP(1) = T0
              YP(1:R,1) = Y0(1:R)
              DO I=2,DBLK+1
                YP(1:R,I) = Y0(1:R)
                TP(I) = TP(I-1)+H
              END DO
            END IF
            IER = 0
            DO J = 2, DBLK+1
               IERR = 0
               CALL FCN(R,TP(J),YP(1,J),FP(1,J), RPAR,IPAR)  ! removed ierr
               IER = MIN(IERR, IER)
            END DO
            NFCN = NFCN + DBLK
 !           IF (IER .NE. 0) THEN
 !             H = MIN(HACC,H/2d0)
 !           END IF

           END DO
         END IF



         !! DAE:  SCALE THE TOLERANCES ACCORDING TO THE INDEX OF THE VARIABLES
          IF (INDEX2) THEN
             DO I=NIND1+1,NIND1+NIND2
                SCAL(I)=SCAL(I)/(MIN(1D0,H*CPORD(ORD)))
             END DO
          END IF
          IF (INDEX3) THEN
             DO I=NIND1+NIND2+1,NIND1+NIND2+NIND3
                SCAL(I)=SCAL(I)/(MIN(1D0,H*H*CPORD(ORD)**2))
             END DO
          END IF

         !!---------- COMPUTE THE NUMERICAL SOLUTION AT TIMES T1(1)...T1(DBLK)
         !!---------- DEFINE VARIABLES NEEDED IN THE ITERATION
         ERRNEWT  = FACN+1d0
         ERRNEWT0 = FACN+1d0
         !!----------
         TETAK    = 1.0D0
         THETA    = 1.0D0
         TETAKOLD = 1.0D0
         THETAPREC = THETA
         ITMAX    = ITINT(ORD)
         IT  = 0

         !!
         !!-------- SPLITTING NEWTON LOOP
         !!
         NEWT_LOOP : DO
            !! 300    CONTINUE
            ERRNEWT0 = ERRNEWT
            ERRNEWT  = 0D0

            !!--------- COMPUTE ONE ITERATION FOR THE SELECTED ORDER
            ORD_NEWT : SELECT CASE(ORD)
            CASE(1)

               CALL      TERMNOT3(R,FCN,H,IT,DN(1), F(1,1),FP(1,1),YP(1,1),TP(1),NFCN,   &
                    &  ERRNEWT,ERRNEWT0,TETAK0(1),LU(1,1), LDLU,IPIV(1),            &
                    &  FMAS(1,1),LDMAS,MLMAS,MUMAS, SCAL(1),IJOB,TER,RPAR,IPAR)

            CASE(2)
               CALL      TERMNOT5(R,FCN,H,IT,DN(1), F(1,1),FP(1,1),YP(1,1),TP(1),NFCN,        &
                    &  ERRNEWT,ERRNEWT0,TETAK0(2),LU(1,1), LDLU,IPIV(1),            &
                    &  FMAS(1,1),LDMAS,MLMAS,MUMAS, SCAL(1),IJOB,TER,RPAR,IPAR)

            CASE(3)
               CALL      TERMNOT7(R,FCN,H,IT,DN(1), F(1,1),FP(1,1),YP(1,1),TP(1),NFCN,        &
                    &  ERRNEWT,ERRNEWT0,TETAK0(3),LU(1,1), LDLU,IPIV(1),            &
                    &  FMAS(1,1),LDMAS,MLMAS,MUMAS, SCAL(1),IJOB,TER,RPAR,IPAR)

            CASE(4)
               CALL      TERMNOT9(R,FCN,H,IT,DN(1), F(1,1),FP(1,1),YP(1,1),TP(1),NFCN,        &
                    &  ERRNEWT,ERRNEWT0,TETAK0(4),LU(1,1), LDLU,IPIV(1),            &
                    &  FMAS(1,1),LDMAS,MLMAS,MUMAS, SCAL(1),IJOB,TER,RPAR,IPAR)

            END SELECT ORD_NEWT

            IF (TER .OR. .NOT.(ERRNEWT .GT. 0.0d0))  THEN
               ERRNEWT = (FACN + 1)
               EXIT NEWT_LOOP
            END IF


            !!--------- COMPUTE TETAK, ETAK

            TETAKOLD = TETAK
            TETAK    = ERRNEWT/ERRNEWT0

            IF (IT.LT.2) THEN
               THETA=TETAK0(ORD)/2
            ELSE IF (IT .EQ. 2) THEN
               THETA = TETAK
            ELSE IF (IT .GT. 2) THEN
               THETA = SQRT(TETAK*TETAKOLD)
            END IF
            IT = IT+1

            JVAI = (IT .LE. ITMAX).AND.(ERRNEWT.GT.FACN) .AND. &
                &((THETA.LT.TETAK0(ORD)).OR.(IT.LE.2)) .AND. &
                &  (ERRNEWT .GT.0d0)
          !   JVAI = (IT .LE. ITMAX).AND.(ERRNEWT.GT.FACN)
            IF (.NOT. (JVAI)) EXIT  NEWT_LOOP
            !!
            !!--------- END OF NEWTON LOOP
            !!
         END DO NEWT_LOOP
         !! 999   CONTINUE



            IF (ERRNEWT.GT.FACN) THEN
            !!--------- THE ITERATION DOES NOT CONVERGE
            FAILNI = FAILNI + 1
            NNEWT(ORD) = NNEWT(ORD)+1
            !!--------- CHOICE OF THE NEW STEPSIZE
            HOLD = H
            H=MIN(HACC,H/2d0)
            DBLKOLD = DBLK
            EXTRAP = .FALSE.
            IF (FAILNI .EQ. 1) THEN
               CALJAC = .NOT. NEWJAC
            ELSE
               CALJAC = .FALSE.
            END IF
            !!-------- RETURN TO THE MAIN LOOP
         ELSE

            !!--------- THE ITERATION CONVERGES
            !!--------- ERROR ESTIMATION

            ERRSAMEOLD = ERRSAME

            CALL  ESTERR(ERRV, ERRSAME, ERRUP, ERRDOWN, FP, &
                 &     R, H, ORD, DBLK, LU(1,1), LDLU, FMAS(1,1), LDMAS, MLMAS, MUMAS, &
                 &     IPIV(1), F(1,1), F1(1,1), SCAL(1), ORDMAX,ORDMIN,IJOB)


            IF (  FAILEI > 5) THEN
                IF ( ABS(ERRSAME-ERRSAMEOLD) < 1d-4) THEN
                   ERRSAME = 0.8d0
                END IF
            END IF




            IF ( ERRSAME .GT. 0.8d0  ) THEN
               FAILEI = FAILEI + 1
               NERR(ORD) = NERR(ORD) + 1
               IF (FAILEI .EQ. 1) THEN
                 CALJAC = (THETA .GT. THET) .and. (.NOT. NEWJAC)
               ELSE
                 CALJAC=.FALSE.
               END IF
               !!           IF (IT .LE. DBLK+1 ) CALJAC = .FALSE.
               !!--------- NEW STEPSIZE SELECTION
               ORD2 = 2*ORD
               ESP = 1D0/(ORD2+1D0)
               RRN=MAX(FACL,MIN(FACR,(SF*ERRSAME)**ESP))
               HOLD = H
               H = H/RRN
               DBLKOLD = DBLK
               F1(1:R,1:DBLKOLD+1) = YP(1:R,1:DBLKOLD+1)
               CALL DIFFDIV(TP,F1,R,DBLK,NT1)
               EXTRAP = .TRUE.

               !!--------- RETURN TO THE MAIN LOOP
            ELSE

               !!--------- THE STEPSIZE IS ACCEPTED

               NSTEP(ORD) = NSTEP(ORD)+1
               T0 = TP(DBLK+1)
               Y0(1:R) = YP(1:R,DBLK+1)
               FP(1:R,1) = FP(1:R,DBLK+1)
               !!--------- NEW STEPSIZE SELECTION
               ORD2 = 2*ORD
               ESP = 1D0/(ORD2+1d0)
               RRN=MAX(FACL,MIN(FACR,(SFSAME*ERRSAME)**ESP))
               THN=DBL(ORD)/(CS(ORD)*RRN)
               ORDN = ORD

               IF  (ORD.LT.ORDMAX) THEN
                     ESP = 1D0/(ORD2+3D0)
                     RR=MAX(FACL,MIN(FACR,(SFUP*ERRUP)**ESP))

                     TH=DBL(ORD+1)/(CS(ORD+1)*RR )
                     IF (TH .GT. THN ) THEN
                        ORDN = ORD + 1
                        RRN  = RR
                        THN  = TH
                     END IF
               END IF

               IF ( ORD.GT.ORDMIN)  THEN
                  ESP = 1D0/(ORD2-1d0)
                  RR=MAX(FACL,MIN(FACR,(SFDOWN*ERRDOWN)**ESP))
                  TH=DBL(ORD-1)/(CS(ORD-1)*RR )
                  IF ( (TH .GT. THN )) THEN
                     ORDN = ORD - 1
                     RRN  = RR
                  END IF
               END IF

               HOLD = H
               HACC = H

               IF (ORDN.GT.ORD) THEN
                  H = MIN(H/RRN,HOLD)
               ELSE
                   H = H/RRN
               END IF


               ORDOLD = ORD
               ORD = ORDN
               DBLKOLD = DBLK
               DBLK = DBL(ORD)

               CALJAC = (THETA .GT. THET)
 !               write(6,*) 'THETA ', THETA,  TETAK0(ORD)/2, THET

               IF ((FAILNI.NE.0).OR.(FAILEI.NE.0)) THEN
                  H = MIN( H, HOLD)
               END IF
               IF  (.NOT. CALJAC) THEN
                  IF ((H/HOLD.LE.1.1D0 ).AND.(H/HOLD.GE.0.9D0))THEN
                     H = HOLD
                  END IF
               END IF
               H = MIN( H, MIN(HMAX, (TEND-T0)/DBLK) )



               F1(1:R,1:DBLKOLD+1) = YP(1:R,1:DBLKOLD+1)
               CALL DIFFDIV(TP,F1,R,DBLKOLD,NT1)
               EXTRAP = .TRUE.

               IF (IOUT.NE.0) THEN
                  !!--------- CALL SOLOUT
                  IRTRN = 0

                  CALL SOLOUT(R,TP,YP,F1 &
                             &      ,NT1,DBLKOLD,ORDOLD,RPAR,IPAR,IRTRN)
                  IF (IRTRN.LT.0) GOTO 800
               END IF
               IF (NSTEPS .EQ. 0) THEN
                  FAILNI = 0
                  FAILEI = 0
               ELSE
                  FAILNI = MAX(FAILNI-1,0)
                  FAILEI = MAX(FAILEI-1,0)

               END IF
               NSING  = 0
            END IF
            !!--------- END IF ERRSAME > 1
         END IF
         !!--------- END IF ERRNEWT > 1
         NSTEPS = NSTEPS + 1
         IF (0.1d0*ABS(T0-TEND)/DBLK .GT. ABS(T0)*EPS ) THEN

            IF (0.1d0*ABS(H) .LE. ABS(T0+EPS)*EPS) THEN
      CALL Rprintd1(' Stepsize too small, h = ',H)
      CALL Rprintd1(' exit of GAM at t = ',T0)
               IDID=-3
               !!  GOTO 800
               EXIT LOOP_TIME
            END IF
            IF (NSTEPS.GT.NMAX) THEN
      CALL Rprinti1('More than nmax = Steps are needed', NMAX)
      CALL Rprintd1(' exit of GAM at t = ',T0)
               IDID=-2
               !!  GOTO 800
               EXIT LOOP_TIME
            END IF
            CYCLE LOOP_TIME
            !! GOTO 100
            !!
            !!---------- END WHILE T0 < T
            !!
         ELSE
            H    = HACC
            IDID = 1
            EXIT LOOP_TIME
         END IF
    END DO LOOP_TIME
    DEALLOCATE( SCAL, YP, FP, F,DN, F1, JF0, LU, FMAS, IPIV)

800      RETURN
     END SUBROUTINE ETRO

END SUBROUTINE GAMD


!
! Karline: made subroutine from CONTR called CONTOUT
!

SUBROUTINE contout(R,T,TP,FF,DBLK,NT1, CONTR)
! ----------------------------------------------------------
!     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT. IT PROVIDES AN
!     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT T.
!     IT GIVES THE VALUE OF THE INTERPOLATION POLYNOMIAL, DEFINED FOR
!     THE LAST SUCCESSFULLY COMPUTED STEP.
! ----------------------------------------------------------
  USE PRECISION
  IMPLICIT NONE
  !!
  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) ::  R, DBLK, NT1
  REAL(PREC), INTENT(OUT) :: CONTR(R)
  REAL(PREC), INTENT(IN) ::  T, TP(DBLK+1), FF(R,DBLK+1)
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER :: I, J, N, NT2, NTi
  REAL(PREC) :: YP
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  N = DBLK+1
  NTi = MAX(1,NT1)
  NT2 = NTi+1
  DO I = 1, R
    YP = FF(I,NTi)
    DO J=NT2,N
       YP = YP*(T-TP(J)) + FF(I,J)
    ENDDO
    CONTR(I) = YP
  ENDDO

  RETURN
END SUBROUTINE contout

