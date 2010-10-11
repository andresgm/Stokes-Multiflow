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
      mwPointer out_normal, out_cmean, out_kg

C     Array information:
      mwSize m_nb, n_nb, size_nb
      mwSize m_normal, n_normal, size_normal
      
C     Arguments for computational routine:
      real*8  tol, cmean, kg
      integer itmax
      real*8 normalin(3), normalout(3), nb(3,6)

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
      call mxCopyPtrToReal8(pr_nb,nb,size_nb)
      pr_normal = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(pr_normal,normalin,size_normal)
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
      call bestparaboloid(nb,normalin,tol,itmax,m_nb,n_nb,
     +                       normalout,cmean,kg)
      
C     Load the data into pointers, which is the output to MATLAB.
      call mxCopyReal8ToPtr(normalout,out_normal,size_normal)
      call mxCopyReal8ToPtr(cmean,out_cmean,1)
      call mxCopyReal8ToPtr(kg,out_kg,1)       

      return
      end

C-----------------------------------------------------------------------
C     Computational routine

      subroutine bestparaboloid(nb,normal_in,tol,itmax,
     +                          m_nb,n_nb,normalout,cmean,kg)
      mwSize m_nb,n_nb,it
      real*8 normalin(3), normalout(3), nb(m_nb,n_nb)
      real*8  tol, cmean, kg
      integer maxit
      real*8 rhs(5), am(5,5)
      integer n,nrhs,lda,ldb,ipiv(5),info
      
      real*8 normalant(3), normalnew(3), normalnewl(3)
      real*8 r1(3),r2(3),r3(3)
      real*8 mr(3,3), eye(3,3), pga(3,3)
      real*8 xlocnodes(3,n_nb), mata(5,5), vectb(5), paramvec(5)
      real*8 xi,yi,zi,dist,a,b,c,d,e
      
      n= 5
      nrhs = 1
      lda = 5
      ldb = 5
      
      eye(:,1) = [1,0,0]
      eye(:,2) = [0,1,0]
      eye(:,3) = [0,0,1]

      normalant = [0,0,0]
      normalnew = normal_in
      it = 0
      do while (normesp(normalnew - normalant) > tol)
      it = it + 1
C     matriz de rotaciones
         r3 = normalnew
         pga = eye-nn(normalnew)
         r1 = matmul(pga,nb(:,1))
         r1 = univect(r1)
         r2 = crossv(r3,r1)
      
         mr(:,1) = r1
         mr(:,2) = r2
         mr(:,3) = r3
         
         mr = transpose(mr)
      
         
         do k = 1,n_nb
             xlocnodes(:,k) = matmul(mr,nb(:,k))
         enddo
         
         mata = 0
         vectb = 0
         
         do i = 1,n_nb
            xi = xlocnodes(1,i)
            yi = xlocnodes(2,i)
            zi = xlocnodes(3,i)
            dist = xi**2+yi**2+zi**2
            mata(1,1) = mata(1,1) + xi**4/dist
            mata(1,2) = mata(1,2) + xi**3*yi/dist
            mata(1,3) = mata(1,3) + yi**2*xi**2/dist
            mata(1,4) = mata(1,4) + xi**3/dist
            mata(1,5) = mata(1,5) + xi**2*yi/dist
            mata(2,1) = mata(2,1) + xi**3*yi/dist
            mata(2,2) = mata(2,2) + xi**2*yi**2/dist
            mata(2,3) = mata(2,3) + xi*yi**3/dist
            mata(2,4) = mata(2,4) + xi**2*yi/dist
            mata(2,5) = mata(2,5) + xi*yi**2/dist
            mata(3,1) = mata(3,1) + xi**2*yi**2/dist
            mata(3,2) = mata(3,2) + xi*yi**3/dist
            mata(3,3) = mata(3,3) + yi**4/dist
            mata(3,4) = mata(3,4) + xi*yi**2/dist
            mata(3,5) = mata(3,5) + yi**3/dist
            mata(4,1) = mata(4,1) + xi**3/dist
            mata(4,2) = mata(4,2) + xi**2*yi/dist
            mata(4,3) = mata(4,3) + xi*yi**2/dist
            mata(4,4) = mata(4,4) + xi**2/dist
            mata(4,5) = mata(4,5) + xi*yi/dist
            mata(5,1) = mata(5,1) + xi**2*yi/dist
            mata(5,2) = mata(5,2) + xi*yi**2/dist
            mata(5,3) = mata(5,3) + yi**3/dist
            mata(5,4) = mata(5,4) + xi*yi/dist
            mata(5,5) = mata(5,5) + yi**2/dist
            vectb(1) = vectb(1) + xi**2*zi/dist
            vectb(2) = vectb(2) + xi*yi*zi/dist
            vectb(3) = vectb(3) + yi**2*zi/dist
            vectb(4) = vectb(4) + xi*zi/dist
            vectb(5) = vectb(5) + yi*zi/dist
         end do
         
         am = mata
         rhs = vectb
         call dgesv(n,nrhs,am,lda,ipiv,rhs,ldb,info)
         paramvec = rhs
         a = paramvec(1)
         b = paramvec(2)
         c = paramvec(3)
         d = paramvec(4)
         e = paramvec(5)
         normalant = normalnew
         normalnewl = [-d,-e,1.d0]/(sqrt(d**2 + e**2 + 1))
C          transforme la normal a la base Global
         normalnew = matmul(mr,normalnewl)
C          curvatura media
         cmean = -(a+c+a*e**2+c*d**2-b*d*e)/(1 + d**2 + e**2)**1.5
C          curvatura gaussiana
         Kg= (4*a*c-b**2)/(1 + d**2 + e**2)**2
         if (it >= itmax) then
            print*,'Calculo fallido en la curvatura, fin de simulacion!'
            return
         endif
      enddo
      normalout = normalnew
      
      return
      
      contains

      function nn(n)
          real*8 nn(3,3), n(3)

          do i=1,3
              do j=1,3
                  nn(i,j) = n(i)*n(j)
              end do
          end do
      end function
      
      function normesp(x)
        real*8 x(3)
        real*8 normesp
        
        normesp = sqrt(x(1)**2+x(2)**2+x(3)**2)
      end function
      
      function univect(x)
        real*8 x(3), univect(3)
        
        univect=x/normesp(x)
      end function
      
      function crossv(a,b)
        real*8 a(3),b(3),crossv(3)
        
        c1 = a(2)*b(3)-a(3)*b(2)
        c2 = a(3)*b(1)-a(1)*b(3)
        c3 = a(1)*b(2)-a(2)-b(1)
        
        crossv = [c1,c2,c3]
      end function
      
      end subroutine
      
      
