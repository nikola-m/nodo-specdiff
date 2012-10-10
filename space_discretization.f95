      module space_discretization
!     Discretization for Spectral Difference method.
      use prec
      implicit none
   
      real(dp), parameter :: pi = 4.*atan(1.0_dp)

      contains

      subroutine gauss_nodes(n,xs)
!     Solution points - a 1d array of Gauss nodes
      implicit none
      integer, intent(in) :: n ! If polynomial is of order (n-1) there are n nodes
      real(dp),dimension(1:n),intent(out) :: xs ! Gauss nodes

!     Locals
      integer :: j
   
      do j=1,n
        xs(j) = 0.5*( 1.-cos( (2*j-1)/dble(2*n) *pi) )
      enddo     

      end subroutine gauss_nodes

      subroutine gl_nodes(n,xs05)
!     Flux points - a 1d array of Gauss-Lobatto nodes
      implicit none
      integer, intent(in) :: n ! If polynomial is of order (n-1) there are n nodes
      real(dp),dimension(0:n),intent(out) :: xs05 ! Gauss-Lobatto nodes

!     Locals
      integer :: j
   
      do j=0,n
        xs05(j) = 0.5*( 1.-cos( j/dble(n) * pi ) )
      enddo     

      end subroutine gl_nodes
      

      end module space_discretization
