!****************************************************************************
! *   FILE         = 2D_HEAT_VERT.f90                                       *
! *   DESCRIPTION  = 1D COLUMN DECOMPOSITION                                *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *
!****************************************************************************
!*                                                                          *
!*          T=100.                                                          *
!*         _________                    ____|____|____|____    Y            *
!*         |       |                   |    |    |    |    |   |            *
!*         |       |                   |    |    |    |    |   |            *
!* T=400.  |T=100.0| T =600.           | 1  | 2  | 3  | 4  |   |            *
!*         |       |                   |    |    |    |    |   |            *
!*         |_______|                   |    |    |    |    |   |_______X    *
!*          T=1000.0                   |____|____|____|____|                *
!*                                          |    |    |                     *
!****************************************************************************
     PROGRAM TWOD_HEAT_VERT
     IMPLICIT NONE
     INCLUDE 'mpif.h'
     INTEGER, PARAMETER :: NX=1000,NY=1000
     REAL, ALLOCATABLE :: T(:,:),TN(:,:), X(:), Y(:)
     REAL, ALLOCATABLE :: TSENDL(:),TSENDR(:),TRECVR(:),TRECVL(:)
     REAL :: DX,DY,S,RMS,G_RMS, XMAX, X0
     INTEGER :: I,J,K,M
     INTEGER :: NUMTASKS,MYID,IERR
     INTEGER :: NMAX, NSTART,NEND,NREM,ICT,ITA
     DOUBLE PRECISION :: T1, T2
     INTEGER :: STATUS(MPI_STATUS_SIZE)

     CHARACTER(LEN=20) :: FILE_ID
     CHARACTER(LEN=50) :: FILE_NAME

     CALL MPI_INIT(IERR)
     CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMTASKS, IERR)
   
     
      DX=5.0/(NX-1-(NUMTASKS-1))
      DY=5.0/(NY-1)
     
      NMAX=NX/NUMTASKS
      NREM=NX-NMAX*NUMTASKS
      IF(MYID.LT.NREM) THEN
       NMAX=NMAX+1
      END IF
      
      
      PRINT*,'PROCESSOR ID:',MYID,'WORKING ON', NMAX, 'NODES'
     
      IF (NREM.EQ.0) THEN
       NSTART=1+NX/NUMTASKS*MYID
      ELSE
       IF(MYID.LT.NREM) NSTART=1+(NX/NUMTASKS+1)*MYID
       IF(MYID.GT.NREM) NSTART=1+NX/NUMTASKS*MYID+NREM
      ENDIF
       NEND=NSTART+NMAX

      IF(MYID.NE.0) NMAX=NMAX+1
      IF(MYID.NE.NUMTASKS-1) NMAX=NMAX+1

 
        ALLOCATE(T(NMAX,NY))
        ALLOCATE(TN(NMAX,NY))
        ALLOCATE(TSENDL(NY),TRECVL(NY),TSENDR(NY),TRECVR(NY))
        ALLOCATE(X(NMAX), Y(NY))

!     GENERATE GRID
       Y(1) = 0.0
      DO  J =2, NY
          Y(J) = Y(J-1) + DY
      END DO

       IF (MYID .EQ. 0) THEN
          X(1) = 0.0
        DO I = 2, NMAX
           X(I) = X(I-1) +DX
        END DO
      END IF

    DO I = 1, NUMTASKS-1
       IF (MYID .EQ. I-1) THEN
          XMAX = X(NMAX)
         CALL MPI_SEND ( XMAX, 1, MPI_REAL, MYID+1,10, MPI_COMM_WORLD, IERR)
       END IF
       IF (MYID .EQ. I) THEN
         CALL MPI_RECV(X0, 1, MPI_REAL, MYID-1,10, MPI_COMM_WORLD, STATUS, IERR)
         X(1) = X0-2.0*DX
         DO J = 2, NMAX
            X(J) = X(J-1) +DX
         END DO
       END IF
   END DO


 !INITIAL CONDITION
      DO I=1,NMAX
       DO J=1,NY
          T(I,J)=100.0
       END DO
      END DO

   !BOUNDARY CONDTION
     DO I=1,NMAX
        T(I,1)=1000.0 !BOTTOM
        T(I,NY)=100.0 !TOP
     END DO
     
     IF (MYID.EQ.0) THEN
       DO J=1,NY
         T(1,J)=400.0 ! LEFT BOUNDARY CONDITON
       END DO
     END IF
     IF(MYID.EQ.NUMTASKS-1) THEN
        DO J=1,NY
           T(NMAX,J)=600.0 !RIGHT BOUNDARY CONDITION
        END DO
     END IF

      T1 = MPI_Wtime()
      ITA=0
        
71    ITA=ITA+1
      DO J=2,NY-1
      DO I=2,NMAX-1
         TN(I,J) = (DY**2*(T(I+1,J)+T(I-1,J)) + DX**2*(T(I,J+1)+T(I,J-1)))/(2.0*DY**2+2.0*DX**2)
      END DO
      END DO
     
      IF (MYID.NE.NUMTASKS-1) THEN
         ICT=1
         DO J=2,NY-1
          TSENDR(ICT)=TN(NMAX-1,J)
          ICT=ICT+1
         END DO
     END IF
     IF(MYID.NE.0) THEN
        ICT=1
        DO J=2,NY-1
         TSENDL(ICT)=TN(2,J)
         ICT=ICT+1
        END DO
     END IF
     

    IF (MYID.LT.NUMTASKS-1) THEN
      CALL MPI_SEND(TSENDR,NY-2,MPI_REAL,MYID+1,10, MPI_COMM_WORLD,IERR)
    END IF
    IF(MYID.GT.0) THEN
      CALL MPI_RECV(TRECVR,NY-2,MPI_REAL,MYID-1,10,MPI_COMM_WORLD,STATUS,IERR)
    END IF 

    IF(MYID.GT.0) THEN
      CALL MPI_SEND(TSENDL,NY-2,MPI_REAL,MYID-1,20,MPI_COMM_WORLD,IERR)
    END IF
    IF(MYID.LT.NUMTASKS-1) THEN
      CALL MPI_RECV(TRECVL,NY-2,MPI_REAL,MYID+1,20,MPI_COMM_WORLD,STATUS,IERR)
    END IF



    IF (MYID.LT.NUMTASKS-1) THEN
        ICT=1
        DO J=2,NY-1
         TN(NMAX,J)=TRECVL(ICT)
        ICT=ICT+1
        END DO
    END IF

   IF(MYID.GT.0) THEN
      ICT=1
      DO J=2,NY-1
       TN(1,J)=TRECVR(ICT)
      ICT=ICT+1
      END DO
   END IF
   
   S=0.0
   DO I=2,NMAX-1
      DO J=2,NY-1
      S=S+(TN(I,J)-T(I,J))**2
      END DO
   END DO
    S=S/((NMAX-2)*(NY-2))
    RMS=SQRT(S)
    
   IF (MYID .GT. 0 .AND. MYID .LT.NUMTASKS-1) THEN
    DO I=1,NMAX
     DO J=2,NY-1
       T(I,J)=TN(I,J)
     END DO
    END DO
   ELSE
    IF(MYID.EQ.0) THEN
      DO I=2,NMAX
        DO J=2,NY-1
         T(I,J)=TN(I,J)
        END DO
      END DO
    END IF
    IF(MYID.EQ.NUMTASKS-1) THEN
       DO I=1,NMAX-1
        DO J=2,NY-1
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
   WRITE(4,*) "ZONE T=  ""N= 0"", I= ", NMAX, ", J= ", NY, ",F=POINT"
        DO J=1,NY 
        DO I=1,NMAX
           WRITE(4,*)X(I),Y(J),T(I,J)
         END DO
       END DO
   ELSE 
     IF(MYID.EQ.0) THEN
      WRITE(4,*) "ZONE T=  ""N= 0"", I= ", NMAX, ", J= ", NY, ",F=POINT"
         DO J=1,NY 
          DO I=1,NMAX
           WRITE(4,*)X(I),Y(J),T(I,J)
         END DO
        END DO
      ENDIF
      IF(MYID.EQ.NUMTASKS-1) THEN
       WRITE(4,*) "ZONE T=  ""N= 0"", I= ", NMAX, ", J= ", NY, ",F=POINT"
         DO J=1,NY 
          DO I=1,NMAX
           WRITE(4,*)X(I),Y(J),T(I,J)
          END DO
         END DO
      END IF
    END IF
    CLOSE (4)
 
   CALL MPI_FINALIZE(IERR)

   
   END PROGRAM 
    
    
    
     

