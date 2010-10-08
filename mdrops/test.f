      program test
          real*8 x1(3),x2(3),x3(3), xx(3,3), eye(3,3)
          real*8 mr(3,3)
          
          mr = 0
          
          print*, mr
      
          x1 = [1,2,3]
          x2 = [4,5,6]
          x3 = [7,8,9]
          
          mr(:,1) = x1
          mr(:,2) = x2
          mr(:,3) = x3
          
          print*, mr
          mr = transpose(mr)
          print*, mr
          
          xx = nn(x1)
      
          eye(:,1) = [1,0,0]
          eye(:,2) = [0,1,0]
          eye(:,3) = [0,0,1]
      
          print*, x1/normesp(x1), normesp(x1/normesp(x1))
      
      contains
      
      function nn(n)
          real*8 nn(3,3), n(:)
      
          do i=1,3
              do j=1,3
                  nn(i,j) = n(i)*n(j)
              end do
          end do
      end function
      
      function normesp(x)
          real*8 x(:)
          real*8 normesp
      
          normesp = sqrt(x(1)**2+x(2)**2+x(3)**2)
      end function
      
      end program

