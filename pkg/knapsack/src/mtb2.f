      SUBROUTINE MTB2(N,P,W,B,C,Z,X,JDIM1,JDIM2,JFO,JFS,JCK,JUB,
     1                ID1,ID2,ID3,ID4,ID5,ID6,ID7,RD8)
C
C THIS SUBROUTINE SOLVES THE BOUNDED SINGLE KNAPSACK PROBLEM
C
C MAXIMIZE  Z = P(1)*X(1) + ... + P(N)*X(N)
C
C SUBJECT TO:   W(1)*X(1) + ... + W(N)*X(N) .LE. C ,
C               0 .LE. X(J) .LE. B(J) FOR J=1,...,N,
C               X(J)  INTEGER         FOR J=1,...,N.
C
C THE PROGRAM IS INCLUDED IN THE VOLUME
C   S. MARTELLO, P. TOTH, "KNAPSACK PROBLEMS: ALGORITHMS
C   AND COMPUTER IMPLEMENTATIONS", JOHN WILEY, 1990
C AND IMPLEMENTS THE TRANSFORMATION METHOD DESCRIBED IN
C SECTION  3.2 .
C THE PROBLEM IS TRANSFORMED INTO AN EQUIVALENT 0-1 KNAPSACK
C PROBLEM AND THEN SOLVED THROUGH SUBROUTINE MT2. THE USER
C MUST LINK MT2 AND ITS SUBROUTINES TO THIS PROGRAM.
C
C THE INPUT PROBLEM MUST SATISFY THE CONDITIONS
C
C   1) 2 .LE. N .LE. JDIM1 - 1 ;
C   2) P(J), W(J), B(J), C  POSITIVE INTEGERS;
C   3) MAX (B(J)*W(J)) .LE. C ;
C   4) B(1)*W(1) + ... + B(N)*W(N) .GT. C ;
C   5) 2 .LE. N + ( LOG2(B(1)) + ... + LOG2(B(N)) ) .LE. JDIM2 - 3 ;
C
C AND, IF  JFS = 1 ,
C
C   6) P(J)/W(J) .GE. P(J+1)/W(J+1) FOR J=1,...,N-1.
C
C MTB2 CALLS  4  PROCEDURES: CHMTB2, SOL, TRANS AND MT2 (EXTERNAL).
C
C COMMUNICATION TO THE PROGRAM IS ACHIEVED SOLELY THROUGH THE PARAMETER
C LIST OF MTB2.
C NO MACHINE-DEPENDENT CONSTANT IS USED.
C THE PROGRAM IS WRITTEN IN 1967 AMERICAN NATIONAL STANDARD FORTRAN
C AND IS ACCEPTED BY THE PFORT VERIFIER (PFORT IS THE PORTABLE
C SUBSET OF ANSI DEFINED BY THE ASSOCIATION FOR COMPUTING MACHINERY).
C THE PROGRAM HAS BEEN TESTED ON A DIGITAL VAX 11/780 AND AN H.P.
C 9000/840.
C
C MTB2 NEEDS
C   4  ARRAYS ( P ,  W ,  B  AND  X ) OF LENGTH AT LEAST  JDIM1 ;
C   8  ARRAYS ( ID1 ,  ID2 ,  ID3 ,  ID4 ,  ID5 ,  ID6 ,  ID7  AND
C               RD8 ) OF LENGTH AT LEAST  JDIM2 .
C
C MEANING OF THE INPUT PARAMETERS:
C N     = NUMBER OF ITEM TYPES;
C P(J)  = PROFIT OF EACH ITEM OF TYPE  J  (J=1,...,N);
C W(J)  = WEIGHT OF EACH ITEM OF TYPE  J  (J=1,...,N);
C B(J)  = NUMBER OF ITEMS OF TYPE  J  AVAILABLE  (J=1,...,N);
C C     = CAPACITY OF THE KNAPSACK;
C JDIM1 = DIMENSION OF P, W, B, X;
C JDIM2 = DIMENSION OF ID1, ID2, ID3, ID4, ID5, ID6, ID7, RD8.
C JFO   = 1 IF OPTIMAL SOLUTION IS REQUIRED,
C       = 0 IF APPROXIMATE SOLUTION IS REQUIRED;
C JFS   = 1 IF THE ITEMS ARE ALREADY SORTED ACCORDING TO
C           DECREASING PROFIT PER UNIT WEIGHT (SUGGESTED
C           FOR LARGE  B(J)  VALUES),
C       = 0 OTHERWISE;
C JCK   = 1 IF CHECK ON THE INPUT DATA IS DESIRED,
C       = 0 OTHERWISE;
C
C MEANING OF THE OUTPUT PARAMETERS:
C Z    = VALUE OF THE SOLUTION FOUND IF  Z .GT. 0 ,
C      = ERROR IN THE INPUT DATA (WHEN JCK=1) IF Z .LT. 0 : CONDI-
C        TION  - Z  IS VIOLATED;
C X(J) = NUMBER OF ITEMS OF TYPE  J  IN THE SOLUTION FOUND;
C JUB  = UPPER BOUND ON THE OPTIMAL SOLUTION VALUE (TO EVALUATE Z
C        WHEN JFO=0).
C
C ARRAYS ID1, ID2, ID3, ID4, ID5, ID6, ID7 AND RD8 ARE DUMMY.
C
C ALL THE PARAMETERS BUT RD8 ARE INTEGER. ON RETURN OF MTB2 ALL THE
C INPUT PARAMETERS ARE UNCHANGED.
C
      INTEGER P(JDIM1),W(JDIM1),B(JDIM1),X(JDIM1),C,Z
      INTEGER ID1(JDIM2),ID2(JDIM2),ID3(JDIM2),ID4(JDIM2),ID5(JDIM2),
     1        ID6(JDIM2),ID7(JDIM2)
      REAL    RD8(JDIM2)
      Z = 0
      IF ( JCK .EQ. 1 ) CALL CHMTB2(N,P,W,B,C,JFS,Z,JDIM1)
      IF ( Z .LT. 0 ) RETURN
C
C TRANSFORM THE BOUNDED KNAPSACK PROBLEM INTO AN EQUIVALENT
C 0-1 KNAPSACK PROBLEM.
C
      CALL TRANS(N,P,W,B,JDIM1,JDIM2,NT,ID1,ID2)
      IF ( NT .GT. 0 ) GO TO 10
      Z = - 5
      RETURN
C SOLVE THE EQUIVALENT 0-1 KNAPSACK PROBLEM.
   10 CALL MT2(NT,ID1,ID2,C,Z,ID3,JDIM2,JFO,JFS,0,JUB,
     1         ID4,ID5,ID6,ID7,RD8)
C DETERMINE THE SOLUTION VECTOR FOR THE ORIGINAL PROBLEM.
      CALL SOL(N,B,ID3,JDIM1,JDIM2,X)
      RETURN
      END
      SUBROUTINE CHMTB2(N,P,W,B,C,JFS,Z,JDIM1)
C
C CHECK THE INPUT DATA.
C
      INTEGER P(JDIM1),W(JDIM1),B(JDIM1),C,Z
      IF ( N .GT. 1 .AND. N .LE. JDIM1 - 1 ) GO TO 10
      Z = - 1
      RETURN
   10 IF ( C .GT. 0 ) GO TO 30
   20 Z = - 2
      RETURN
   30 JSW = 0
      DO 40 J=1,N
        IF ( P(J) .LE. 0 ) GO TO 20
        IF ( W(J) .LE. 0 ) GO TO 20
        IF ( B(J) .LE. 0 ) GO TO 20
        JSW = JSW + B(J)*W(J)
        IF ( B(J)*W(J) .GT. C ) GO TO 50
   40 CONTINUE
      IF ( JSW .GT. C ) GO TO 60
      Z = - 4
      RETURN
   50 Z = - 3
      RETURN
   60 IF ( JFS .EQ. 0 ) RETURN
      RR = FLOAT(P(1))/FLOAT(W(1))
      DO 70 J=2,N
        R = RR
        RR = FLOAT(P(J))/FLOAT(W(J))
        IF ( RR .GT. R ) GO TO 80
   70 CONTINUE
      RETURN
   80 Z = - 6
      RETURN
      END
      SUBROUTINE SOL(N,B,XT,JDIM1,JDIM2,X)
C
C DETERMINE THE SOLUTION VECTOR X FOR THE ORIGINAL PROBLEM.
C
      INTEGER B(JDIM1),XT(JDIM2),X(JDIM1)
      NT = 0
      DO 20 J=1,N
        ISUM = 0
        ID = 1
        X(J) = 0
   10   NT = NT + 1
        ISUM = ISUM + ID
        X(J) = X(J) + ID*XT(NT)
        ID = ID*2
        IF ( ID + ISUM .LE. B(J) ) GO TO 10
        IF ( ISUM .EQ. B(J) ) GO TO 20
        ID = B(J) - ISUM
        NT = NT + 1
        X(J) = X(J) + ID*XT(NT)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE TRANS(N,P,W,B,JDIM1,JDIM2,NT,PT,WT)
C
C TRANSFORM A BOUNDED KNAPSACK PROBLEM (N, P, W, B) INTO
C A 0-1 KNAPSACK PROBLEM (NT, PT, WT ).
C
      INTEGER P(JDIM1),W(JDIM1),B(JDIM1),PT(JDIM2),WT(JDIM2)
      JDMAX = JDIM2 - 3
      NT = 0
      DO 20 J=1,N
        ISUM = 0
        ID = 1
   10   NT = NT + 1
        IF ( NT .GT. JDMAX ) GO TO 30
        PT(NT) = P(J)*ID
        WT(NT) = W(J)*ID
        ISUM = ISUM + ID
        ID = ID*2
        IF ( ID + ISUM .LE. B(J) ) GO TO 10
        IF ( ISUM .EQ. B(J) ) GO TO 20
        ID = B(J) - ISUM
        NT = NT + 1
        IF ( NT .GT. JDMAX ) GO TO 30
        PT(NT) = P(J)*ID
        WT(NT) = W(J)*ID
   20 CONTINUE
      RETURN
   30 NT = - 5
      RETURN
      END