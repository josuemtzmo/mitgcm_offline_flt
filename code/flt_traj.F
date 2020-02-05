#include "FLT_OPTIONS.h"

CBOP 0
C !ROUTINE: FLT_TRAJ

C !INTERFACE:
      SUBROUTINE FLT_TRAJ (
     I                      myTime, myIter, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE FLT_TRAJ
C     | o This routine samples the model state at float position
C     |   every flt_int_traj time steps and writes output.
C     *==========================================================*

C     !USES:
      IMPLICIT NONE
C     == global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "FLT_SIZE.h"
#include "FLT.h"
#include "FLT_BUFF.h"
#ifdef ALLOW_EXCH2
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"
#endif

C     !INPUT PARAMETERS:
C     myTime :: current time in simulation
C     myIter :: current iteration number
C     myThid :: my Thread Id number
      _RL myTime
      INTEGER myIter, myThid

C     !FUNCTIONS:
      _RL FLT_MAP_K2R
      EXTERNAL FLT_MAP_K2R

C     !LOCAL VARIABLES:
      INTEGER bi, bj, nFlds
      INTEGER ip, kp, ii, i, j
      _RL ix, jy, i0x, j0y, xx, yy, zz
      _RL uu, vv, tt, ss, pp, ttke, vort
      _RL KE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL vort3(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER imax, iG, jG
      PARAMETER (imax=16)
      _RL tmp(imax), tmp2(imax)
      _RL npart_read, npart_times
      _RS dummyRS(1)
      INTEGER fp, ioUnit, irecord, pirecord
      CHARACTER*(MAX_LEN_FNAM) fn
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifdef ALLOW_EXCH2
      INTEGER nT
#endif
CEOP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   set number of fields to write
      nFlds = 0
      IF ( flt_selectTrajOutp.GE.1 ) nFlds = nFlds + 8
      IF ( flt_selectTrajOutp.GE.2 ) nFlds = nFlds + 5
C--   Modified section to write KE and vorticity along the particle track.
      IF ( flt_selectTrajOutp.GE.3 ) THEN 
        nFlds = nFlds + 3
        KE = 1/2*((uVel*uVel) + (vVel*vVel))
        DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)
            DO j=2,sNy
              DO i=2,sNx
                vort3(i,j) = recip_rAz(i,j,bi,bj)*(
     &               vVel( i ,j,1,bi,bj)*dyC( i ,j,bi,bj)
     &               -vVel(i-1,j,1,bi,bj)*dyC(i-1,j,bi,bj)
     &               -uVel(i, j ,1,bi,bj)*dxC(i, j ,bi,bj)
     &               +uVel(i,j-1,1,bi,bj)*dxC(i,j-1,bi,bj)
     &                                                    )
C               Look at file mon_vort3.F
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C--   check buffer size
      IF ( nFlds.GT.fltBufDim ) THEN
         _BEGIN_MASTER(myThid)
         WRITE(msgBuf,'(3(A,I4))') ' FLT_TRAJ: fltBufDim=', fltBufDim,
     &                             ' too small (<', nFlds, ' )'
         CALL PRINT_ERROR( msgBuf, myThid )
         WRITE(msgBuf,'(2A)')     ' FLT_TRAJ: => increase fltBufDim',
     &                            ' in "FLT_SIZE.h" & recompile'
         CALL PRINT_ERROR( msgBuf, myThid )
         _END_MASTER(myThid)
         CALL ALL_PROC_DIE( myThid )
         STOP 'ABNORMAL END: S/R FLT_TRAJ'
      ENDIF

      IF ( myIter.EQ.nIter0 .OR. flt_selectTrajOutp.LE.0 ) RETURN
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Calculate position + other fields at float position and fill up IO-buffer

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

#ifdef ALLOW_EXCH2
         nT = W2_myTileList(bi,bj)
         i0x = DFLOAT( exch2_txGlobalo(nT) - 1 )
         j0y = DFLOAT( exch2_tyGlobalo(nT) - 1 )
#else
         i0x = DFLOAT( myXGlobalLo-1 + (bi-1)*sNx )
         j0y = DFLOAT( myYGlobalLo-1 + (bj-1)*sNy )
#endif
         DO ip=1,npart_tile(bi,bj)

            ix = ipart(ip,bi,bj)
            jy = jpart(ip,bi,bj)
            print*,'Particle',ip
            print*,ix,jy
            print*,bi,bj
            CALL FLT_MAP_IJLOCAL2XY( xx, yy,
     I                               ix, jy, bi,bj, myThid )
            print*,xx,yy
            zz = FLT_MAP_K2R( kpart(ip,bi,bj),bi,bj,myThid )
            kp = NINT(kpart(ip,bi,bj))
            tmp(1) = npart(ip,bi,bj)
            tmp(2) = myTime
            tmp(3) = xx
            tmp(4) = yy
            tmp(5) = zz
            tmp(6) = ix + i0x
            tmp(7) = jy + j0y
            tmp(8) = kpart(ip,bi,bj)

            IF ( ( flt_selectTrajOutp.GE.2 )   .AND.
     &           ( myTime.GE.tstart(ip,bi,bj)) .AND.
     &           ( tend(ip,bi,bj).EQ.-1. .OR. myTime.LE.tend(ip,bi,bj))
     &         ) THEN
              IF ( kp.LT.1 .OR. kp.GT.Nr ) THEN
                WRITE(msgBuf,'(2A,I8)') '** WARNING ** FLT_TRAJ: ',
     &            ' illegal value for kp=',kp
                CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                              SQUEEZE_RIGHT, myThid )
                WRITE(msgBuf,'(A,1P5E20.13)')
     &            ' FLT_TRAJ: ', (flt_io_buff(ii,ip,bi,bj),ii=1,5)
                CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                              SQUEEZE_RIGHT, myThid )
c               CALL PRINT_ERROR( msgBuf, myThid )
c               STOP 'ABNORMAL END: S/R FLT_TRAJ'
C-- jmc: not sure if this is right but added to avoid Pb in FLT_BILINEAR:
                kp = MIN( MAX(kp,1), Nr)
              ENDIF
              CALL FLT_BILINEAR  (ix,jy,uu,uVel,  kp,1,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,vv,vVel,  kp,2,bi,bj,myThid)
              CALL FLT_BILINEAR2D (ix,jy,pp,etaN,     0,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,tt,theta, kp,0,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,ss,salt,  kp,0,bi,bj,myThid)
              tmp( 9) = pp
              tmp(10) = uu
              tmp(11) = vv
              
              tmp(12) = tt
              tmp(13) = ss
              
              IF ( ( flt_selectTrajOutp.GE.3 ) ) THEN
                CALL FLT_BILINEAR  (ix,jy,vort,vort3,
     &                               kp,0,bi,bj,myThid)
                CALL FLT_BILINEAR  (ix,jy,ttke, KE, 
     &                               kp,1,bi,bj,myThid)
                tmp(14) = vort
                tmp(15) = ttke
                tmp(16) = 1
                CLOSE( ioUnit )
              ENDIF
C DELETE NEXT 2 LINES
              tmp(12) = vort
              tmp(13) = ttke

            ELSEIF ( flt_selectTrajOutp.GE.2 ) THEN
              tmp( 9) = flt_nan
              tmp(10) = flt_nan
              tmp(11) = flt_nan
              tmp(12) = flt_nan
              tmp(13) = flt_nan
            ENDIF

            DO ii=1,nFlds
              flt_io_buff(ii,ip,bi,bj) = tmp(ii)
            ENDDO

         ENDDO

       ENDDO
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Write shared buffer to file

      _BARRIER
      _BEGIN_MASTER(myThid)

      fn = 'float_trajectories'
      fp = writeBinaryPrec

      DO bj=1,nSy
       DO bi=1,nSx

C (1) read actual number floats from file (if exists)
         ioUnit = -2
         CALL MDS_READVEC_LOC(  fn, fp, ioUnit, 'RL', nFlds,
     O                          tmp, dummyRS,
     I                          bi, bj, 1, myThid )
         IF ( ioUnit.GT. 0 ) THEN
            npart_read  = tmp(1)
            npart_times = tmp(5)
            ii = NINT(tmp(7))
C-       for backward compatibility with old trajectory files:
            IF ( ii.EQ.0 ) ii = 15
            IF ( ii.NE.nFlds ) THEN
              WRITE(msgBuf,'(A,I4,A)')
     &            'FLT_TRAJ: nFlds=', nFlds,' different from'
              CALL PRINT_ERROR( msgBuf, myThid )
              WRITE(msgBuf,'(3A,I4,A)')
     &            'previous file (',fn(1:18),') value =',ii
              CALL PRINT_ERROR( msgBuf, myThid )
              CALL ALL_PROC_DIE( 0 )
              STOP 'ABNORMAL END: S/R FLT_TRAJ'
            ENDIF
C-       close the read-unit (safer to use a different unit for writing)
            CLOSE( ioUnit )
         ELSE
            npart_read  = 0.
            npart_times = 0.
            tmp(2)      = myTime
         ENDIF

C (2) write new actual number floats and time axis into file
C-    the standard routine mds_writevec_loc can be used here

C     total number of records in this file
         tmp(1) = DBLE(npart_tile(bi,bj))+npart_read
C     first time of writing floats (do not change when written)
c        tmp(2) = tmp(2)
C     current time
         tmp(3) = myTime
C     timestep
         tmp(4) = flt_int_traj
C     total number of timesteps
         tmp(5) = npart_times + 1.
C     total number of floats
         tmp(6) = max_npart
C     total number of fields
         tmp(7) = nFlds
         DO ii=8,nFlds
           tmp(ii) = 0.
         ENDDO
         ioUnit = -1
         CALL MDS_WRITEVEC_LOC( fn, fp, ioUnit, 'RL', nFlds,
     &                          tmp, dummyRS,
     &                          bi, bj, -1, myIter, myThid )

         iG=bi+(myXGlobalLo-1)/sNx
         jG=bj+(myYGlobalLo-1)/sNy

         DO ip=1,npart_tile(bi,bj)
C            irecord = npart_read+ip+1
            irecord = ip+1
            IF ( ( flt_selectTrajOutp.GE.3 ) ) THEN
              IF ( npart_read.NE.0 ) THEN
C              pirecord=irecord-npart_tile(bi,bj)
              pirecord=irecord
              CALL MDS_READVEC_LOC(  fn, fp, ioUnit, 'RL', nFlds,
     O                               tmp2, dummyRS,
     I                               bi, bj, pirecord, myThid )
              DO ii=14,nFlds
                flt_io_buff(ii,ip,bi,bj) = (flt_io_buff(ii,ip,bi,bj)
     &                                      + tmp2(ii))
              ENDDO
              ENDIF
            ENDIF
            IF ( ip.NE.npart_tile(bi,bj) ) irecord = -irecord
            
            CALL MDS_WRITEVEC_LOC( fn, fp, ioUnit, 'RL', nFlds,
     I                             flt_io_buff(1,ip,bi,bj), dummyRS,
     I                             bi, bj, irecord, myIter, myThid )
            
         ENDDO
         CLOSE( ioUnit )

       ENDDO
      ENDDO

      _END_MASTER(myThid)
      _BARRIER

      RETURN
      END
