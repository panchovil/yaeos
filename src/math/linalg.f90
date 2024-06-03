module yaeos__math_linalg
   !! Wrapper module around LAPACK's `dgesv`
   use yaeos__constants, only: pr
   implicit none

   abstract interface
      subroutine f_1d(x, f, df)
         import pr
         real(pr), intent(in) :: x
         real(pr), intent(out) :: f
         real(pr), intent(out) :: df
      end subroutine
   end interface

   interface newton
      module procedure :: newton_1d
   end interface

contains
   function solve_system(a, b) result(x)
      !! Solve a linear sytem AX = b
      real(pr), intent(in) :: b(:)
      real(pr), intent(in) :: a(size(b), size(b))
      integer, parameter :: dp = selected_real_kind(15)

      real(pr) :: x(size(b))

      real(dp) :: a_lapack(size(b), size(b)), b_lapack(size(b))
      integer :: n, nrhs, lda, ipiv(size(b)), ldb, info

      interface
         subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            integer :: n
            integer :: nrhs
            real(dp) :: a(n, n)
            integer :: lda
            integer :: ipiv(n)
            real(dp) :: b(n)
            integer :: ldb
            integer :: info
         end subroutine
      end interface

      n = size(a, dim=1)
      nrhs = 1
      lda = n
      ldb = n

      a_lapack = a
      b_lapack = b
      call dgesv(n, nrhs, a_lapack, lda, ipiv, b_lapack, ldb, info)

      if (info > 0) error stop 1

      x = b_lapack
   end function solve_system

   subroutine newton_1d(f, x, tol, max_iters)
      procedure(f_1d) :: f
      real(pr), intent(in out) :: x
      real(pr), intent(in) :: tol
      integer, intent(in) :: max_iters

      integer :: i
      real(pr) :: error, fval, df, step

      error = 10

      step = 10

      do i=1, max_iters
         if (abs(error) < tol .or. abs(step) < tol)  exit
         call f(x, fval, df)

         step = fval/df

         do while (abs(step) > 0.5 * abs(x))
            step = step/2
         end do

         x = x - step
      end do
   end subroutine
end module