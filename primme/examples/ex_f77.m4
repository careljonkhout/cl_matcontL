C*******************************************************************************
C  Copyright (c) 2018, College of William & Mary                                   
C  All rights reserved.                                                            
C                                                                                  
C  Redistribution and use in source and binary forms, with or without
C  modification, are permitted provided that the following conditions are met:     
C      * Redistributions of source code must retain the above copyright
C        notice, this list of conditions and the following disclaimer.             
C      * Redistributions in binary form must reproduce the above copyright         
C        notice, this list of conditions and the following disclaimer in the       
C        documentation and/or other materials provided with the distribution.      
C      * Neither the name of the College of William & Mary nor the
C        names of its contributors may be used to endorse or promote products      
C        derived from this software without specific prior written permission.     
C                                                                                  
C  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
C  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
C  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          
C  DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY       
C  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      
C  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
C  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
C  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
C  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C  
C  PRIMME: https://github.com/primme/primme
C  Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
*******************************************************************************
*
*  Example to compute the k largest eigenvalues in a 1-D Laplacian matrix.
*
*******************************************************************************
define(`PRIMME_SCALAR', ifdef(`USE_PETSC', `PetscScalar', ifdef(`USE_COMPLEX', `complex*16', `real*8')))dnl
define(`PRIMME_REAL', ifdef(`USE_PETSC', `PetscReal', `real*8'))dnl

        Program primmeF77Example
!-----------------------------------------------------------------------
ifdef(`USE_PETSC', ``#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscmat.h>
        use petscvec
        use petscmat
        use petscpc
'')dnl
        implicit none
        include 'primme_f77.h'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Pointer to the PRIMME data structure used internally by PRIMME
!
!       Note that for 64 bit systems, pointers are 8 bytes so use:
        integer*8 primme
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Problem setup
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        ! Solver Parameters
        integer*8 n,NUMEmax,BASISmax,BLOCKmax,maxMatvecs,
     :          printLevel, method, whichEvals, numTargetShifts
        real*8 ETOL

        parameter (
     :            n               = 100,
     :            BASISmax        = 12,
     :            NUMEmax         = 5,
     :            BLOCKmax        = 1,
     :            maxMatvecs      = 300000,
     :            ETOL            = ifdef(`USE_PETSC', 1.0D-5, 1.0D-12),
     :            printLevel      = 5,
     :            whichEvals      = primme_smallest,
     :            numTargetShifts = 2,
     :            method          = PRIMME_DYNAMIC
     :  )
        real*8 TargetShifts(numTargetShifts)
        data TargetShifts /3.0D0, 5.1D0/
ifdef(`USE_PETSC', `
        external generateLaplacian1D, PETScMatvec, ApplyPCPrecPETSC,
     :           par_GlobalSum
', `
        external MV, ApplyPrecon
')dnl

!       Eigenvalues, eigenvectors, and their residual norms
!
        PRIMME_REAL evals(NUMEmax)
        PRIMME_REAL rnorms(NUMEmax)
        PRIMME_SCALAR evecs(n*NUMEmax)

!       Other vars
!
ifdef(`USE_PETSC', `ifdef(`USE_POINTER',
`        Mat, target :: A
        PC, target :: pc
        MPI_Comm, target :: comm
', `        Mat A
        PC pc
        COMMON A, pc
')dnl
        PetscErrorCode ierr
        PetscInt nLocal
        integer*8 numProcs,procID
        integer i,numProcsMpi,procIDMpi
', `        integer i,ierr
')dnl
        real*8  epsil, aNorm
        integer*8 numIts, numMatvecs

!-----------------------------------------------------------------------
!       Start executable 
!-----------------------------------------------------------------------
!
ifdef(`USE_PETSC', `        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
')dnl
!       ----------------------------------------------------------------
!       Initialize PRIMME
!       ----------------------------------------------------------------
!
        call primme_initialize_f77(primme)

!       Set a few basic solver parameters
        call primme_set_member_f77(primme, PRIMME_n, n, ierr)
        call primme_set_member_f77(primme, PRIMME_numEvals, NUMEmax,
     :                                                             ierr)
        call primme_set_member_f77(primme, PRIMME_eps, ETOL, ierr)
        call primme_set_member_f77(primme, PRIMME_target,
     :                                                 whichEvals, ierr)
        call primme_set_member_f77(primme, PRIMME_numTargetShifts, 
     :                                            numTargetShifts, ierr)
        call primme_set_member_f77(primme, PRIMME_targetShifts, 
     :                                               TargetShifts, ierr)

!       Set matvec 
ifdef(`USE_PETSC', `        call generateLaplacian1D(n, A, ierr)
ifdef(`USE_POINTER', `        call primme_set_member_f77(primme, PRIMME_matrix, A, ierr)
')dnl
        call primme_set_member_f77(primme, PRIMME_matrixMatvec,
     :                                                PETScMatvec, ierr)
', `        call primme_set_member_f77(primme, PRIMME_matrixMatvec,
     :                                                         MV, ierr)
')dnl
        
ifdef(`USE_PETSC', `!       Set parallel parameters
        call MatGetLocalSize(A, nLocal, PETSC_NULL_INTEGER, ierr)
        call primme_set_member_f77(primme, PRIMME_nLocal, nLocal, ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD, numProcsMpi, ierr)
        numProcs = numProcsMpi
        call primme_set_member_f77(primme, PRIMME_numProcs, numProcs,
     :                                                             ierr)
        call MPI_Comm_rank(PETSC_COMM_WORLD, procIDMpi, ierr);
        procID = procIDMpi
        call primme_set_member_f77(primme, PRIMME_procID,procID,ierr)
ifdef(`USE_POINTER', `        comm = PETSC_COMM_WORLD
        call primme_set_member_f77(primme, PRIMME_commInfo,comm,ierr)
')dnl
        call primme_set_member_f77(primme, PRIMME_globalSumReal,
     :                                        par_GlobalSum, ierr)
')dnl

!       Set preconditioner  (optional)
ifdef(`USE_PETSC', `        pc = tPC(0)
        call PCCreate(PETSC_COMM_WORLD, pc, ierr)
        call PCSetType(pc, PCJACOBI, ierr)
        call PCSetOperators(pc, A, A, ierr)
        call PCSetFromOptions(pc, ierr)
        call PCSetUp(pc, ierr)
ifdef(`USE_POINTER', `        call primme_set_member_f77(primme, 
     :       PRIMME_preconditioner, pc, ierr)
')dnl
        call primme_set_member_f77(primme, 
     :       PRIMME_applyPreconditioner, ApplyPCPrecPETSC, ierr)
', `        call primme_set_member_f77(primme, 
     :       PRIMME_applyPreconditioner, ApplyPrecon, ierr)
')dnl
        call primme_set_member_f77(primme, 
     :       PRIMME_correctionParams_precondition, 0, ierr)
!
!       Set a few other solver parameters (optional) 
!
        call primme_set_member_f77(primme, PRIMME_maxBasisSize, 
     :                                                   BASISmax, ierr)
        call primme_set_member_f77(primme, PRIMME_maxBlockSize,
     :                                                   BLOCKmax, ierr)
        call primme_set_member_f77(primme, PRIMME_printLevel, 
     :                                                 printLevel, ierr)
        call primme_set_member_f77(primme, PRIMME_maxMatvecs,
     :                                                 maxMatvecs, ierr)
!
!       Set the method to be used (after n, numEvals, and precondition have
!       been set. Also after basisSize is set, if desired.)
        call primme_set_method_f77(primme, method, ierr)

        if (ierr .lt. 0) 
     :     write(*,*) 'No preset method. Using custom settings'

!       ----------------------------------------------------------------
!       Display what parameters are used
!       ----------------------------------------------------------------

       ifdef(`USE_PETSC', ` if (procID.eq.0)') call primme_display_params_f77(primme)

!       ----------------------------------------------------------------
!       Calling the PRIMME solver
!       ----------------------------------------------------------------
ifdef(`USE_PETSC', ``
#if defined(PETSC_USE_COMPLEX) && defined(PETSC_USE_REAL_SINGLE)
        call cprimme_f77(evals, evecs, rnorms, primme, ierr)
#elif defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_REAL_SINGLE)
        call zprimme_f77(evals, evecs, rnorms, primme, ierr)
#elif !defined(PETSC_USE_COMPLEX) && defined(PETSC_USE_REAL_SINGLE)
        call sprimme_f77(evals, evecs, rnorms, primme, ierr)
#elif !defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_REAL_SINGLE)
        call dprimme_f77(evals, evecs, rnorms, primme, ierr)
#endif
'', `
        call ifdef(`USE_COMPLEX', ifdef(`USE_FLOAT',`c',`z'), ifdef(`USE_FLOAT',`s',`d'))primme_f77(evals, evecs, rnorms, primme, ierr)
')dnl
        if (ierr.ne.0) then
          stop 1
        endif

!       ----------------------------------------------------------------
!       Reporting results

ifdef(`USE_PETSC', ``        if (procID.eq.0) then
' define(sp, `   ')', `define(sp, `')')dnl
        sp()if (ierr.eq.0) then
        sp()   print *, 'PRIMME has returned successfully'
        sp()else 
        sp()   print *, 'PRIMME returned with error: ', ierr
        sp()endif

!       sp()
!       sp()Example of obtaining primme members from the driver:
!       sp()NOTE: don't use primme_get_member_f77, which can only be used in a callback
!
        sp()call primmetop_get_member_f77(primme, PRIMME_eps, epsil,
     :  sp()                                                      ierr)
        sp()call primmetop_get_member_f77(primme, PRIMME_aNorm,
     :  sp()                                              aNorm, ierr)
        sp()call primmetop_get_member_f77(primme,
     :  sp()         PRIMME_stats_numOuterIterations, numIts, ierr)
        sp()call primmetop_get_member_f77(primme,
     :  sp()              PRIMME_stats_numMatvecs, numMatvecs, ierr)
        sp()print '(A,E8.2,/,A,e12.5,/,A,I8,/,A,I8)',
     :  sp()                           'Tolerance used:   ',epsil,
     :  sp()                           'Estimated norm(A):',aNorm,
     :  sp()                           'Iterations:       ',numIts,
     :  sp()                           'Matvecs:          ',numMatvecs
!
!       sp()Reporting of evals and residuals
!
        sp()do i = 1, numemax
        sp()   write (*, 9000) i, evals(i),rnorms(i)
        sp()enddo
 9000   sp()FORMAT (1x,'E(',i1,') = ',G24.16,4x,
     &  sp()       'residual norm =', E12.4)
        call primme_free_f77(primme)
ifdef(`USE_PETSC',`        endif

        call PetscFinalize(ierr)
')dnl
        stop
        write(0,*) 'ERROR! No data in the file'
        stop
        end
!-----------------------------------------------------------------------
! Supporting subroutines
!-----------------------------------------------------------------------
!       ----------------------------------------------------------------
changequote(`[',`]')
ifdef([USE_PETSC], [
        subroutine generateLaplacian1D(n0,A,ierr)
!       ----------------------------------------------------------------
        use petscmat
        implicit none
        integer*8 n0
        PetscInt n, one, two, three
        Mat A

        PetscScalar value(3)
        PetscInt i, Istart,Iend,col(3)
        PetscErrorCode ierr

        A = tMat(0)
        call MatCreate(PETSC_COMM_WORLD, A, ierr)
        n = n0
        call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
        call MatSetUp(A, ierr)
        call MatGetOwnershipRange(A, Istart, Iend, ierr)
        one = 1
        two = 2
        three = 3
        do i=Istart,Iend-1
           if (i.eq.0) then
              col = (/0, 1, 0/)
              value = (/2.0, -1.0, 0.0/)
              call MatSetValues(A, one, i, two, col, value,
     :                                             INSERT_VALUES, ierr)
           else if (i.eq.n-1) then
              col = (/n-2, n-1, 0/)
              value = (/-1.0, 2.0, 0.0/)
              call MatSetValues(A, one, i, two, col, value,
     :                                             INSERT_VALUES, ierr)
           else
              col = (/i-1, i, i+1/)
              value = (/-1.0, 2.0, -1.0/)
              call MatSetValues(A, one, i, three, col, value,
     :                                             INSERT_VALUES, ierr)
           endif
           call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
           call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
        enddo
        end
        subroutine PETScMatvec(x,ldx,y,ldy,k,primme,err)
!       ----------------------------------------------------------------
ifdef([USE_POINTER], [        use iso_c_binding
])dnl
        use petscvec
        use petscmat
        implicit none
        include 'primme_f77.h'
        integer*8 ldx,ldy
        PRIMME_SCALAR x(ldx,*), y(ldy,*)
        integer*8 primme
        integer k,err,j
ifdef([USE_POINTER], [        Mat, pointer :: A
        type(c_ptr) :: pA
], [        Mat A
        COMMON A
])dnl
        Vec xvec,yvec
        PetscErrorCode ierr

ifdef([USE_POINTER], [        call primme_get_member_f77(primme, PRIMME_matrix, pA, err)
        call c_f_pointer(pA, A)
])
        xvec = tVec(0)
        yvec = tVec(0)
        call MatCreateVecs(A, xvec, yvec, ierr)
        do j=1,k
           call VecPlaceArray(xvec, x(1,j), ierr)
           call VecPlaceArray(yvec, y(1,j), ierr)
           call MatMult(A, xvec, yvec, ierr)
           call VecResetArray(xvec, ierr)
           call VecResetArray(yvec, ierr)
        enddo
        call VecDestroy(xvec, ierr)
        call VecDestroy(yvec, ierr)
        err = 0
        end
        subroutine ApplyPCPrecPETSc(x,ldx,y,ldy,k,primme,err)
!       ----------------------------------------------------------------
ifdef([USE_POINTER], [        use iso_c_binding
])dnl
        use petscvec
        use petscmat
        use petscpc
        implicit none
        include 'primme_f77.h'
        integer*8 ldx,ldy
        PRIMME_SCALAR x(ldx,*), y(ldy,*)
        integer*8 primme
        integer k,err,j
ifdef([USE_POINTER], [        Mat, pointer :: A
        PC, pointer :: pc
        type(c_ptr) :: pA, ppc
], [        Mat A
        PC pc
        COMMON A, pc
])dnl
        Vec xvec,yvec
        PetscErrorCode ierr

ifdef([USE_POINTER], [        call primme_get_member_f77(primme, PRIMME_matrix, pA, err)
        call primme_get_member_f77(primme, PRIMME_preconditioner,
     :                                                        ppc, err)
        call c_f_pointer(pA, A)
        call c_f_pointer(ppc, pc)
])
        xvec = tVec(0)
        yvec = tVec(0)
        call MatCreateVecs(A, xvec, yvec, ierr)
        do j=1,k
           call VecPlaceArray(xvec, x(1,j), ierr)
           call VecPlaceArray(yvec, y(1,j), ierr)
           call PCApply(pc, xvec, yvec, ierr)
           call VecResetArray(xvec, ierr)
           call VecResetArray(yvec, ierr)
        enddo
        call VecDestroy(xvec, ierr)
        call VecDestroy(yvec, ierr)
        err = 0
        end
        subroutine par_GlobalSum(x,y,k,primme,ierr)
!       ----------------------------------------------------------------
        use iso_c_binding
        use petscsys
        implicit none
        include 'primme_f77.h'
        real*8, target :: x(k), y(k)
        integer*8 primme
        integer k,ierr
ifdef([USE_POINTER], [        MPI_Comm, pointer :: comm
        type(c_ptr) :: pcomm

        call primme_get_member_f77(primme, PRIMME_commInfo, pcomm, ierr)
        call c_f_pointer(pcomm, comm)
])dnl
        if (c_associated(c_loc(x),c_loc(y))) then
          call MPI_Allreduce(MPI_IN_PLACE, y, k, MPIU_REAL, MPIU_SUM,
     :                                 ifdef([USE_POINTER], [comm], [PETSC_COMM_WORLD]), ierr)
        else
          call MPI_Allreduce(x, y, k, MPIU_REAL, MPIU_SUM,
     :                                 ifdef([USE_POINTER], [comm], [PETSC_COMM_WORLD]), ierr)
        endif
        end
], [
!       1-D Laplacian block matrix-vector product, Y = A * X, where
!      
!       - X, input dense matrix of size primme.n x blockSize;
!       - Y, output dense matrix of size primme.n x blockSize;
!       - A, tridiagonal square matrix of dimension primme.n with this form:
!      
!            [ 2 -1  0  0  0 ... ]
!            [-1  2 -1  0  0 ... ]
!            [ 0 -1  2 -1  0 ... ]
!             ...
!      
        subroutine MV(x,ldx,y,ldy,k,primme,ierr)
!       ----------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
        integer*8 ldx, ldy
        PRIMME_SCALAR x(ldx,*), y(ldy,*)
        integer*8 primme
        integer*8 n, i
        integer k,ierr,j
        call primme_get_member_f77(primme, PRIMME_n, n, ierr)
        do j=1,k
           do i=1,n
              y(i,j) = 0
              if (i.ge.2) then
                 y(i,j) = y(i,j) - x(i-1,j)
              endif
              y(i,j) = y(i,j) + 2.*x(i,j)
              if (i.le.n-1) then
                 y(i,j) = y(i,j) - x(i+1,j)
              endif
           enddo
        enddo
        ierr = 0
        end

!       This performs Y = M^{-1} * X, where
!      
!       - X, input dense matrix of size primme.n x blockSize;
!       - Y, output dense matrix of size primme.n x blockSize;
!       - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
!      
        subroutine ApplyPrecon(x,ldx,y,ldy,k,primme, ierr)
!       ----------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
        integer*8 ldx, ldy
        PRIMME_SCALAR x(ldx,*), y(ldy,*)
        integer*8 primme
        integer*8 n, i
        integer k,ierr,j
        call primme_get_member_f77(primme, PRIMME_n, n, ierr)
        do j=1,k
           do i=1,n
              y(i,j) = x(i,j)/2.0
           enddo
        enddo
        ierr = 0
        end
])dnl
