C------------------------------------------------------------------
C       Initial nuclear Hessian guess per Fischer, Almlof (J. Phys.
C       Chem. 1992, 96,9768-9774), with (NXM6) int. coords.
C
C       Does not yet include force constants for out-of-plane modes
C------------------------------------------------------------------

            NATOMS=NX/3
            NATOMSQ=NATOMS**2
            CALL GETREC(20, 'JOBARC', 'IBONDTO ', NATOMSQ, IBNDTO)
C
C       Calculation of no. of bonds attached to each atom, L'
C
            DO I=1,NATOMS
               LP(I)=0
               DO J=1,NATOMS
                  LP(I)=LP(I)+IBNDTO(I,J)
               END DO
            END DO

            CALL GETREC(20, 'JOBARC', 'CONCTVTY', NX, NCON)
            CALL GETREC(20, 'JOBARC', 'COORD', NX, CARTCOORD)
            CALL GETREC(20, 'JOBARC', 'ATOMCHRG', NATOMS, ATCNM)

            CALL GETCOVLRADI(ATCNM, SMOFCOVRADI, NATOMS)

            IJ=1
            DO I=NATOMS,2
C
C       Force constants for stretching modes
C
               J=3*I
               ATNM=3*I-2
               K=3*NCON(ATNM)

               CARTXSQ=(CARTCOORD(J-2)-CARTCOORD(K-2))**2
               CARTYSQ=(CARTCOORD(J-1)-CARTCOORD(K-1))**2
               CARTZSQ=(CARTCOORD(J)-CARTCOORD(K))**2
               DISTJK=SQRT(CARTXSQ+CARTYSQ+CARTZSQ)
               DUMMY1=DISTJK-SMOFCOVRADI(I,NCON(ATNM))

               HI(IJ,IJ)=0.3601*EXP(-1.944*DUMMY1)
               IJ=IJ+1
C
C       Force constants for bending modes
C
               IF (NATOMS.GT.2) THEN
                  ATNM2=3*I-1
                  L=3*NCON(ATNM2)

                  CARTXSQ=(CARTCOORD(K-2)-CARTCOORD(L-2))**2
                  CARTYSQ=(CARTCOORD(K-1)-CARTCOORD(L-1))**2
                  CARTZSQ=(CARTCOORD(K)-CARTCOORD(L))**2
                  DISTKL=SQRT(CARTXSQ+CARTYSQ+CARTZSQ)
                  DUMMY2=DISTKL+DISTJK-SMOFCOVRADI(I,NCON(ATNM))
                &        -SMOFCOVRADI(NCON(ATNM),NCON(ATNM2))

                  HI(IJ,IJ)=0.089+0.11*EXP(-0.44*DUMMY2)
                  IJ=IJ+1
C
C       Force constants for torsional modes
C
                  IF (NATOMS.GT.3) THEN
                     NUMBND=LP(NCON(ATNM))+LP(NCON(ATNM2))-2
                     DUMMY3=DUMMY2-DUMMY1
                     DUMMY4=DISTKL
                    &       *SMOFCOVRADI(NCON(ATNM),NCON(ATNM2))

                     HI(IJ,IJ)=0.0015+(14.0*(NUMBND**0.57)
                    &          *EXP(-2.85*DUMMY3))/(DUMMY4**4.0)
                     IJ=IJ+1
                  END IF
               END IF
            END DO
C------------------------------------------------------------------
C       End of Fischer-Almlof Hessian guess code
C------------------------------------------------------------------
