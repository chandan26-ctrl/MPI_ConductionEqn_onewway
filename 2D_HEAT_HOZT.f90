!****************************************************************************
! *   FILE         = 2D_HEAT_HOZT.f90                                       *
! *   DESCRIPTION  = 1D ROW DECOMPOSITION                                   *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *
!****************************************************************************
!*                                                                          *
!*          T=100.                                                          *
!*         _________                    ___________________                 *
!*         |       |                   |       1           |                *
!*         |       |                 --|-------------------|--              *
!* T=400.  |T=100.0| T =600.           |       2           |                *
!*         |       |                 --|-------------------|--              *
!*         |_______|                   |       3           |                *
!*          T=1000.0                 --|-------------------|--              *
!                                      |__ ____4_____ _____|                *
!*                                                                          *
!****************************************************************************
     PROGRAM TWOD_HEAT_HOZT
     IMPLICIT NONE
     INCLUDE 'mpif.h'
     INTEGER, PARAMETER :: NX=1000,NY=1000
     REAL, ALLOCATABLE :: T(:,:),TN(:,:), X(:), Y(:)
     REAL, ALLOCATABLE :: TSENDT(:),TSENDB(:),TRECVB(:),TRECVT(:)
     REAL :: DX,DY,S,RMS,G_RMS, YMAX, Y0
     INTEGER :: I,J,K,M
     INTEGER :: NUMTASKS,MYID,IERR
     INTEGER :: NYMAX, NSTART,NEND,NREM,ICT,ITA
     DOUBLE PRECISION :: T1, T2
     INTEGER :: STATUS(MPI_STATUS_SIZE)

     CHARACTER(LEN=20) :: FILE_ID
     CHARACTER(LEN=50) :: FILE_NAME

     CALL MPI_INIT(IERR)
     CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMTASKS, IERR)
   
     
      DX=5.0/(NX-1)
      DY=5.0/(NY-1-(NUMTASKS-1))
     
      NYMAX=NY/NUMTASKS
      NREM=NY-NYMAX*NUMTASKS
      IF(MYID.LT.NREM) THEN
       NYMAX=NYMAX+1
      END IF
      
      PRINT*,'PROCESSOR ID:',MYID,'WORKING ON', NYMAX, 'NODES'
     
      IF (NREM.EQ.0) THEN
       NSTART=1+NY/NUMTASKS*MYID
      ELSE
       IF(MYID.LT.NREM) NSTART=1+(NY/NUMTASKS+1)*MYID
       IF(MYID.GT.NREM) NSTART=1+NY/NUMTASKS*MYID+NREM
      ENDIF
       NEND=NSTART+NYMAX

      IF(MYID.NE.0) NYMAX=NYMAX+1
      IF(MYID.NE.NUMTASKS-1) NYMAX=NYMAX+1

 
        ALLOCATE(T(NX,NYMAX))
        ALLOCATE(TN(NX,NYMAX))
        ALLOCATE(TSENDT(NX),TRECVT(NX),TSENDB(NX),TRECVB(NX))
        ALLOCATE ( X(NX), Y(NYMAX))

        X(1) = 0.0
        DO I = 2, NX
           X(I) = X(I-1) +DX
        END DO

       IF (MYID .EQ. 0) THEN
          Y(1) = 0.0
        DO I = 2, NYMAX
           Y(I) = Y(I-1) +DY
        END DO
      END IF

    DO I = 1, NUMTASKS-1
       IF (MYID .EQ. I-1) THEN
          YMAX = Y(NYMAX)
         CALL MPI_SEND ( YMAX, 1, MPI_REAL, MYID+1,10, MPI_COMM_WORLD, IERR)
       END IF
       IF (MYID .EQ. I) THEN
         CALL MPI_RECV(Y0, 1, MPI_REAL, MYID-1,10, MPI_COMM_WORLD, STATUS, IERR)
         Y(1) = Y0-2.0*DY
         DO J = 2, NYMAX
            Y(J) = Y(J-1) +DY
         END DO
       END IF
   END DO
       



 !INITIAL CONDITION
      DO I=1,NX
       DO J=1,NYMAX
          T(I,J)=100.0
       END DO
      END DO

   !BOUNDARY CONDTION
     DO J=1,NYMAX
        T(1,J)=400.0 !LEFT
        T(NX,J)=600.0 !RIGHT
     END DO
     
     IF (MYID.EQ.0) THEN
       DO I=1,NX
         T(I,1)=1000.0! BOTTOM
       END DO
     END IF
     IF(MYID.EQ.NUMTASKS-1) THEN
        DO I=1,NX
           T(I,NYMAX)=100.0 !TOP
        END DO
     END IF

      T1 = MPI_Wtime()

      ITA=0
        
71    ITA=ITA+1
      DO J=2,NYMAX-1
      DO I=2,NX-1
         TN(I,J) = (DY**2*(T(I+1,J)+T(I-1,J)) + DX**2*(T(I,J+1)+T(I,J-1)))/(2.0*DY**2+2.0*DX**2)
      END DO
      END DO
     
      IF (MYID.NE.NUMTASKS-1) THEN
         ICT=1
         DO I=2,NX-1
          TSENDT(ICT)=TN(I,NYMAX-1)
          ICT=ICT+1
         END DO
     END IF
     IF(MYID.NE.0) THEN
        ICT=1
         DO I=2,NX-1
         TSENDB(ICT)=TN(I,2)
         ICT=ICT+1
        END DO
     END IF
     

    IF (MYID.LT.NUMTASKS-1) THEN
      CALL MPI_SEND(TSENDT,NX-2,MPI_REAL,MYID+1,10, MPI_COMM_WORLD,IERR)
    END IF
    IF(MYID.GT.0) THEN
      CALL MPI_RECV(TRECVT,NX-2,MPI_REAL,MYID-1,10,MPI_COMM_WORLD,STATUS,IERR)
    END IF

    

    IF(MYID.GT.0) THEN
      CALL MPI_SEND(TSENDB,NX-2,MPI_REAL,MYID-1,20,MPI_COMM_WORLD,IERR)
    END IF
    IF(MYID.LT.NUMTASKS-1) THEN
      CALL MPI_RECV(TRECVB,NX-2,MPI_REAL,MYID+1,20,MPI_COMM_WORLD,STATUS,IERR)
    END IF

   


    IF (MYID.LT.NUMTASKS-1) THEN
        ICT=1
         DO I=2,NX-1
         TN(I,NYMAX)=TRECVB(ICT)
        ICT=ICT+1
        END DO
    END IF

   IF(MYID.GT.0) THEN
      ICT=1
      DO I=2,NX-1
       TN(I,1)=TRECVT(ICT)
      ICT=ICT+1
      END DO
   END IF
   
   S=0.0
   DO I=2,NX-1
      DO J=2,NYMAX-1
      S=S+(TN(I,J)-T(I,J))**2
      END DO
   END DO
    S=S/((NX-2)*(NYMAX-2))
    RMS=SQRT(S)
    
   IF (MYID .GT. 0 .AND. MYID .LT.NUMTASKS-1) THEN
    DO I=2,NX-1
     DO J=1,NYMAX
       T(I,J)=TN(I,J)
     END DO
    END DO
   ELSE
    IF(MYID.EQ.0) THEN
      DO I=2,NX-1
        DO J=2,NYMAX
         T(I,J)=TN(I,J)
        END DO
      END DO
    END IF
    IF(MYID.EQ.NUMTASKS-1) THEN
       DO I=2,NX-1
        DO J=1,NYMAX-1
          T(I,J)=TN(I,J)
        END DO
       END DO
    END IF
   END IF

   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(RMS,G_RMS,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,IERR)

   IF(MYID.EQ.0) THEN
    PRINT*, 'ERROR=', G_RMS,'ITERATION',ITA
   END IF

   IF(G_RMS.GT.1.0E-5) GOTO 71
   IF(MYID.EQ.0) PRINT*,'CONVERGED'
    T2 = MPI_Wtime()
   IF (MYID .EQ. 0) THEN
    PRINT*, 'ELAPSED TIME:', T2-T1     
   END IF

! WRITE FILES FOR TECPLOT FROM EACH PROCESSOR
   WRITE(FILE_ID, '(I6)') MYID+1
   FILE_NAME = 'PRINT' // TRIM(ADJUSTL(FILE_ID)) // '.dat'

   OPEN(FILE=TRIM(FILE_NAME), UNIT = 4)
   

   WRITE(4,*) " TITLE = ""FLOW-FIELD 2D"" "
   WRITE(4,*)  "VARIABLES =  X, Y, T"
   
   IF(MYID.GT.0 .AND. MYID.LT.NUMTASKS-1)THEN
   WRITE(4,*) "ZONE T=  ""N= 0"", I= ", NX, ", J= ", NYMAX, ",F=POINT"
        DO J=1,NYMAX 
        DO I=1,NX
           WRITE(4,*)X(I),Y(J),T(I,J)
         END DO
       END DO
   ELSE 
     IF(MYID.EQ.0) THEN
      WRITE(4,*) "ZONE T=  ""N= 0"", I= ", NX, ", J= ", NYMAX, ",F=POINT"
         DO J=1,NYMAX 
          DO I=1,NX
           WRITE(4,*)X(I),Y(J),T(I,J)
         END DO
        END DO
      ENDIF
      IF(MYID.EQ.NUMTASKS-1) THEN
       WRITE(4,*) "ZONE T=  ""N= 0"", I= ", NX, ", J= ", NYMAX, ",F=POINT"
         DO J=1,NYMAX
          DO I=1,NX
           WRITE(4,*)X(I),Y(J),T(I,J)
          END DO
         END DO
      END IF
    END IF
    CLOSE (4)
 
   CALL MPI_FINALIZE(IERR)

   
   END PROGRAM 

    
    
    
     

