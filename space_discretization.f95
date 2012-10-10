      module space_discretization
!     Discretization for Spectral Difference method.
      use prec
      implicit none
   
      real(dp), parameter :: pi = 4.*atan(1.0_dp)

      contains

      subroutine gauss_nodes(n,xs)
!     Solution points - a 1d array of Chebyshev-Gauss nodes
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
!     Flux points - a 1d array of Chebyshev-Gauss-Lobatto nodes
      implicit none
      integer, intent(in) :: n ! If polynomial is of order (n-1) there are n nodes
      real(dp),dimension(0:n),intent(out) :: xs05 ! Gauss-Lobatto nodes

!     Locals
      integer :: j
   
      do j=0,n
        xs05(j) = 0.5*( 1.-cos( j/dble(n) * pi ) )
      enddo     

      end subroutine gl_nodes
      
      function lagrange_basis_h(i,x,xs) result(h)
!     i-th Lagrange basis function for reconstruction of solution.
!     We use it to construct a polynomial of order N-1,
!     using the values at N solution points.
!     Arguments:
!     i  - which basis func. (integer)
!     x - at what point we evaluete Lagrange basis func (real)
!     xs - Chebyshev-Gauss nodes (1d array of reals, array size, N)
      implicit none
      real(dp) :: h
!
      integer, intent(in) :: i 
      real(dp), intent(in) :: x
      real(dp),dimension(1:n),intent(in) :: xs

      x_i = xs(i)

      h = 0.0_dp
      do j=1,n
      if (j==i) cycle
        h = h * (x - xs(j))/(x_i-xs(j))
      end do
      
      end function lagrange_basis_h

      function lagrange_basis_l(i,x,xs05) result(l05)
!     i-th Lagrange basis function for reconstruction of the flux.
!     We use it to construct a polynomial of order N,
!     using the values at N+1 flux points.
!     Arguments:
!     i  - which basis func. (integer)
!     x - at what point we evaluete Lagrange basis func (real)
!     xs05 - Chebyshev-Gauss-Lobatto nodes (1d array of reals, array size, N+1)
      implicit none
      real(dp) :: l05
!
      integer, intent(in) :: i 
      real(dp), intent(in) :: x
      real(dp),dimension(0:n),intent(in) :: xs05

      x_i05 = xs05(i)

      h = 0.0_dp
      do j=0,n
      if (j==i) cycle
        l05 = l05 * (x - xs05(j))/(x_i05-xs05(j))
      end do
      
      end function lagrange_basis_l

      end module space_discretization
