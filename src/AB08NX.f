      SUBROUTINE AB08NX( N, M, P, RO, SIGMA, SVLMAX, ABCD, LDABCD,
     $                   NINFZ, INFZ, KRONL, MU, NU, NKROL, TOL, IWORK,
     $                   DWORK, LDWORK, INFO )
C
C     RELEASE 4.0, WGS COPYRIGHT 1999.
C
C     PURPOSE
C
C     To extract from the (N+P)-by-(M+N) system
C                  ( B  A )
C                  ( D  C )
C     an (NU+MU)-by-(M+NU) "reduced" system
C                  ( B' A')
C                  ( D' C')
C     having the same transmission zeros but with D' of full row rank.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of state variables.  N >= 0. 
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0. 
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0. 
C
C     RO      (input/output) INTEGER
C             On entry,
C             = P     for the original system;
C             = MAX(P-M, 0) for the pertransposed system.
C             On exit, RO contains the last computed rank.
C     
C     SIGMA   (input/output) INTEGER
C             On entry,
C             = 0  for the original system;
C             = M  for the pertransposed system.
C             On exit, SIGMA contains the last computed value sigma in
C             the algorithm.
C     
C     SVLMAX  (input) DOUBLE PRECISION
C             An estimate of the largest singular value of the original
C             matrix ABCD (for instance, the Frobenius norm of ABCD).
C             SVLMAX >= 0.
C            
C     ABCD    (input/output) DOUBLE PRECISION array, dimension 
C             (LDABCD,M+N)
C             On entry, the leading (N+P)-by-(M+N) part of this array
C             must contain the compound input matrix of the system.
C             On exit, the leading (NU+MU)-by-(M+NU) part of this array
C             contains the reduced compound input matrix of the system.
C     
C     LDABCD  INTEGER
C             The leading dimension of array ABCD.
C             LDABCD >= MAX(1,N+P).
C     
C     NINFZ   (input/output) INTEGER
C             On entry, the currently computed number of infinite zeros.
C             It should be initialized to zero on the first call.
C             NINFZ >= 0.
C             On exit, the number of infinite zeros.
C
C     INFZ    (input/output) INTEGER array, dimension (N)
C             On entry, INFZ(i) must contain the current number of
C             infinite zeros of degree i, where i = 1,2,...,N, found in
C             the previous call(s) of the routine. It should be
C             initialized to zero on the first call.
C             On exit, INFZ(i) contains the number of infinite zeros of
C             degree i, where i = 1,2,...,N.
C
C     KRONL   (input/output) INTEGER array, dimension (N+1)
C             On entry, this array must contain the currently computed
C             left Kronecker (row) indices found in the previous call(s)
C             of the routine. It should be initialized to zero on the
C             first call.
C             On exit, the leading NKROL elements of this array contain
C             the left Kronecker (row) indices.
C
C     MU      (output) INTEGER
C             The normal rank of the transfer function matrix of the
C             original system.
C
C     NU      (output) INTEGER
C             The dimension of the reduced system matrix and the number
C             of (finite) invariant zeros if D' is invertible.
C
C     NKROL   (output) INTEGER
C             The number of left Kronecker indices.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance used in rank decisions to determine the
C             effective rank, which is defined as the order of the
C             largest leading (or trailing) triangular submatrix in the
C             QR (or RQ) factorization with column (or row) pivoting
C             whose estimated condition number is less than 1/TOL.
C             NOTE that when SVLMAX > 0, the estimated ranks could be
C             less than those defined above (see SVLMAX).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (MAX(M,P+1))
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  
C             LDWORK >= MAX( MIN(P+1,M) + MAX(3*M,N),
C                            MIN(P+1,N) + MAX(3*(P+1),N+M) ).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     REFERENCES
C
C     [1] Svaricek, F.
C         Computation of the Structural Invariants of Linear
C         Multivariable Systems with an Extended Version of
C         the Program ZEROS.
C         System & Control Letters, 6, pp. 261-266, 1985.
C
C     [2] Emami-Naeini, A. and Van Dooren, P.
C         Computation of Zeros of Linear Multivariable Systems.
C         Automatica, 18, pp. 415-430, 1982.
C
C     NUMERICAL ASPECTS
C                            
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
C     Supersedes Release 2.0 routine AB08BZ by F. Svaricek.
C
C     REVISIONS
C
C     V. Sima, Oct. 1997, Feb. 1998.
C     A. Varga, May 1999.
C
C     KEYWORDS
C
C     Generalized eigenvalue problem, Kronecker indices, multivariable
C     system, orthogonal transformation, structural invariant.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDABCD, LDWORK, M, MU, N, NINFZ, NKROL,
     $                  NU, P, RO, SIGMA
      DOUBLE PRECISION  SVLMAX, TOL
C     .. Array Arguments ..
      INTEGER           INFZ(*), IWORK(*), KRONL(*)
      DOUBLE PRECISION  ABCD(LDABCD,*), DWORK(*)
C     .. Local Scalars ..
      INTEGER           I1, IK, IROW, ITAU, IZ, JWORK, M1, MM1, MNTAU,
     $                  MNU, N1, RANK, RO1, TAU, WRKOPT
      DOUBLE PRECISION  AII
C     .. Local Arrays ..
      DOUBLE PRECISION  SVAL(3)
C     .. External Subroutines ..
      EXTERNAL          DLAPMT, DLARF, DLASET, DORM2R, MB03OY, MB03PY,
     $                  MB04ID, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( RO.NE.P .AND. RO.NE.MAX( P-M, 0 ) ) THEN
         INFO = -4
      ELSE IF( SIGMA.NE.0 .AND. SIGMA.NE.M ) THEN
         INFO = -5
      ELSE IF( SVLMAX.LT.ZERO ) THEN
         INFO = -6
      ELSE IF( LDABCD.LT.MAX( 1, N+P ) ) THEN
         INFO = -8
      ELSE IF( NINFZ.LT.0 ) THEN
         INFO = -9
      ELSE IF( LDWORK.LT.MAX( MIN( P+1, M ) + MAX( 3*M, N ),
     $                        MIN( P+1, N ) + MAX( 3*(P+1), N+M ) ) )
     $      THEN
         INFO = -18
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB08NX', -INFO )
         RETURN
      END IF
C
      MU = P
      NU = N
C
      IZ = 0
      IK = 1
      MM1 = M + 1
      NKROL = 0
      WRKOPT = 1
C
   20 IF ( MU.EQ.0 )
     $   GO TO 80
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.)
C
      RO1 = RO
      MNU = M + NU
      ITAU = 1
      IF ( M.GT.0 ) THEN
         IROW = NU + 1
         IF ( SIGMA.EQ.0 ) THEN
            SIGMA = 1
            M1 = 0
         ELSE
            RO1 = RO1 + 1
            M1  = SIGMA - 1
            IF ( M1.GT.0 ) THEN
               JWORK = ITAU + M1
C
C              Compress rows of D.  First exploit triangular shape.
C              Workspace: need   max(SIGMA,2*SIGMA-3,M+N);
C                         prefer larger.
C
               CALL MB04ID( RO+M1, M1, M1-1, MNU-M1, ABCD(IROW,1),
     $                      LDABCD, ABCD(IROW,SIGMA), LDABCD,
     $                      DWORK(ITAU), DWORK(JWORK), LDWORK-JWORK+1,
     $                      INFO )
               CALL DLASET( 'Lower', RO+M1, M1, ZERO, ZERO,
     $                      ABCD(IROW+1,1), LDABCD )
               WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
            END IF
         END IF
C
C        Continue with Householder with column pivoting.
C        Workspace: need   min(RO1,M-SIGMA+1) + max(1,3*(M-SIGMA+1)).
C
         JWORK = ITAU + MIN( RO1, M-M1 )
C
         CALL MB03OY( RO1, M-M1, ABCD(NU+SIGMA,SIGMA), LDABCD, TOL,
     $                SVLMAX, RANK, SVAL, IWORK, DWORK(ITAU),
     $                DWORK(JWORK), INFO )
         WRKOPT = MAX( WRKOPT, JWORK + 3*( M-M1 ) - 1 )
C
C        Apply the column permutations to submatrices B and part of D.
C
         CALL DLAPMT( .TRUE., NU+M1, M-M1, ABCD(1,SIGMA), LDABCD, 
     $                IWORK )
C
         IF ( RANK.GT.0 ) THEN
C
C           Apply the Householder transformations to the submatrix C.
C           Workspace: need   min(RO1,M-SIGMA+1) + NU.
C           
            CALL DORM2R( 'Left', 'Transpose', RO1, NU, RANK,
     $                   ABCD(NU+SIGMA,SIGMA), LDABCD, DWORK(ITAU),
     $                   ABCD(NU+SIGMA,MM1), LDABCD, DWORK(JWORK), 
     $                   INFO )
            WRKOPT = MAX( WRKOPT, MIN( RO1, M-M1 ) + NU )
            CALL DLASET( 'Lower', RO1-1, MIN( RO1-1, RANK ), ZERO, ZERO,
     $                   ABCD(NU+SIGMA+1,SIGMA), LDABCD )
            RO1 = RO1 - RANK
         END IF
      END IF
C
      TAU = RO1
      SIGMA = MU - TAU
C
C     Determination of the orders of the infinite zeros.
C
      IF ( IZ.GT.0 ) THEN
         INFZ(IZ) = INFZ(IZ) + RO - TAU
         NINFZ = NINFZ + IZ*( RO - TAU )
      END IF
      IF ( RO1.EQ.0 )
     $   GO TO 80
      IZ = IZ + 1
C
      IF ( NU.LE.0 ) THEN
         MU = SIGMA
         NU = 0
         RO = 0
      ELSE
C
C        Compress the columns of C using RQ factorization with row
C        pivoting, P * C = R * Q.
C        Workspace: need   min(P+1,NU) + 3*(P+1).
C
         I1 = NU + SIGMA
         MNTAU = MIN( TAU, NU )
         JWORK = ITAU + MNTAU
C
         CALL MB03PY( TAU, NU, ABCD(I1+1,MM1), LDABCD, TOL, SVLMAX,
     $                RANK, SVAL, IWORK, DWORK(ITAU), DWORK(JWORK),
     $                INFO )
         WRKOPT = MAX( WRKOPT, JWORK + 3*TAU - 1 )
         IF ( RANK.GT.0 ) THEN
            IROW = I1 + TAU
            N1 = NU
C
            DO 40 ITAU = MNTAU, MNTAU - RANK + 1, -1
C
C              Apply Q(itau) to the first N1 columns of A from the
C              right.
C              Workspace: need   min(P+1,NU) + NU + SIGMA.
C
               AII = ABCD(IROW,M+N1)
               ABCD(IROW,M+N1) = ONE
               CALL DLARF( 'Right', I1, N1, ABCD(IROW,MM1), LDABCD,
     $                     DWORK(ITAU), ABCD(1,MM1), LDABCD,
     $                     DWORK(JWORK) )
               ABCD(IROW,M+N1) = AII
               IROW = IROW - 1
               N1 = N1 - 1
   40       CONTINUE
C
            WRKOPT = MAX( WRKOPT, I1 + JWORK - 1 )
            IROW = I1 + TAU
            N1 = NU
C
            DO 60 ITAU = MNTAU, MNTAU - RANK + 1, -1
C
C              Apply Q(itau) to the first N1 rows and M + N1 columns of
C              [ B  A ] from the left.
C              (LAPACK routine DORMR2 cannot be used efficiently.)
C              Workspace: need   min(P+1,NU) + M + NU.
C
               AII = ABCD(IROW,M+N1)
               ABCD(IROW,M+N1) = ONE
               CALL DLARF( 'Left', N1, M+N1, ABCD(IROW,MM1), LDABCD,
     $                     DWORK(ITAU), ABCD, LDABCD, DWORK(JWORK) )
               ABCD(IROW,M+N1) = AII
               IROW = IROW - 1
               N1 = N1 - 1
   60       CONTINUE
C
            WRKOPT = MAX( WRKOPT, MNU + JWORK - 1 )
C
            CALL DLASET( 'Full', RANK, NU-RANK, ZERO, ZERO,
     $                   ABCD(I1+TAU-RANK+1,MM1), LDABCD )
            IF ( RANK.GT.1 )
     $         CALL DLASET( 'Lower', RANK-1, RANK-1, ZERO, ZERO,
     $                      ABCD(I1+TAU-RANK+2,MM1+NU-RANK), LDABCD )
         END IF
C
         RO = RANK
      END IF
C
C     Determine the left Kronecker indices (row indices).
C
      KRONL(IK) = KRONL(IK) + TAU - RO
      NKROL = NKROL + KRONL(IK)
      IK = IK + 1
C
      NU = NU - RO
      MU = SIGMA + RO
      IF ( RO.NE.0 )
     $   GO TO 20
C
   80 CONTINUE
      DWORK(1) = WRKOPT
      RETURN
C *** Last line of AB08NX ***
      END
