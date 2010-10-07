#include "fintrf.h"
C======================================================================
#if 0
C     
C     bestparaboloid.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     bestparaboloid.f
C
C     Computational function that takes a scalar and doubles it.
      
C     This is a MEX-file for MATLAB.
C     Copyright 2010 Andres Gonzalez-Mancera
C======================================================================
C     Gateway subroutine
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwSize mxGetM, mxGetN

C     Pointers to input/output mxArrays:
      mwPointer pr_nb,pr_normal,pr_tol,pr_itmax, pr_out

C     Array information:
      mwSize m_nb, n_nb, size_nb
      mwSize m_normal, n_normal, size_normal
      
C     Arguments for computational routine:
      real*8  tol, cmean, kg
      integer maxit
      real*8 normalin(3), normalout(3)

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 4) then
         call mexErrMsgIdAndTxt ('MATLAB:bestparaboloid:nInput',
     +                           'Four input required.')
      elseif(nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:bestparaboloid:nOutput',
     +                           'Too many output arguments.')
      endif

C     Get the size of the input arrays.
      m_nb = mxGetM(prhs(1))
      n_nb = mxGetN(prhs(1))
      size_nb = m_nb * n_nb
      m_normal = mxGetM(prhs(2))
      n_normal = mxGetN(prhs(2))
      size_normal = m_normal * n_normal

C     Create Fortran array from the input arguments.
      pr_nb = mxGetPr(prhs(1))
      pr_normal = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(pr_normal,normal_in,size_normal)
      pr_tol = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(pr_tol,tol,1)
      pr_itmax = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(pr_itmax,itmax,1)

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(m_normal, n_normal, 0)
      out_normal = mxGetPr(plhs(1))
      plhs(2) = mxCreateDoubleMatrix(1, 1, 0)
      out_cmean = mxGetPr(plhs(2))
      plhs(3) = mxCreateDoubleMatrix(1, 1, 0)
      out_kg = mxGetPr(plhs(3))

C     Call the computational routine.
      call bestparaboloid(%VAL(pr_nb),normal_in,tol,itmax,m_nb,n_nb,
     +                       normalout,cmean,kg)
      
C     Load the data into pointers, which is the output to MATLAB.
      call mxCopyReal8ToPtr(normalout,out_normal,size_normal)
      call mxCopyReal8ToPtr(cmean,out_cmean,1)
      call mxCopyReal8ToPtr(kg,out_kg,1)       

      return
      end

C-----------------------------------------------------------------------
C     Computational routine

      subroutine bestparaboloid(%VAL(pr_nb),normal_in,tol,itmax,
     +                          m_nb,n_nb,normalout,cmean,kg)
      real*8 x_input, y_output

      y_output = 2.0 * x_input
      return
      end
