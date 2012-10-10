      program test
!     Test space discretization module
      use prec
      use space_discretization
      implicit none

      integer, parameter :: n = 4
      real(dp), dimension(1:n) :: xs
      real(dp), dimension(0:n) :: xs05

      integer :: i

      call gauss_nodes(n,xs)
      write (*,*) (xs(i),i=1,n)

      call gl_nodes(n,xs05)
      write (*,*) (xs05(i),i=0,n)

      end program test
