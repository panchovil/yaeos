module yaeos__fitting
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel, CubicEoS
   use yaeos__equilibria, only: &
      EquilibriaState, saturation_pressure, saturation_temperature, flash
   implicit none

   type :: FittingProblem
      !! Fitting problem setting
      !!
      !! # Description
      !! This derived type holds all the relevant information for a parameter
      !! optimization problem. It keeps the base model structure that will be
      !! optimized and a procedure `get_model_from_X` that should reconstruct
      !! the model with the desired parameters to optimize.
      real(pr) :: solver_tolerance = 1e-5_pr
      real(pr) :: parameter_step = 0.1_pr

      class(ArModel), allocatable :: initial_model

      type(EquilibriaState), allocatable :: experimental_points(:)
      procedure(model_from_X), pointer :: get_model_from_X
   end type FittingProblem

   abstract interface
      subroutine model_from_X(self, X, model)
         import ArModel, FittingProblem, pr
         real(pr), intent(in) :: X(:)
         class(FittingProblem), intent(in) :: self
         class(ArModel), intent(out) :: model
      end subroutine model_from_X
   end interface

   type(CubicEoS) :: model

contains

   real(pr) function optimize(X, func_data) result(y)
      use nlopt_wrap, only: create, nlopt_opt, nlopt_algorithm_enum
      use nlopt_callback, only: nlopt_func, create_nlopt_func

      real(pr), intent(in out) :: X(:) !! Vector of parameters to fit
      type(FittingProblem) :: func_data !! Parametrization details

      real(pr) :: dx(size(X))

      type(nlopt_opt) :: opt !! Optimizer
      type(nlopt_func) :: f !! Function to optimize
      integer :: stat

      ! opt = nlopt_opt(nlopt_algorithm_enum%LN_NELDERMEAD, size(X))
      opt = nlopt_opt(nlopt_algorithm_enum%LN_PRAXIS, size(X))
      ! opt = nlopt_opt(nlopt_algorithm_enum%LN_COBYLA, size(X))

      f = create_nlopt_func(fobj, f_data=func_data)

      dx = func_data%parameter_step
      call opt%set_ftol_rel(func_data%solver_tolerance)

      call opt%set_initial_step(dx)
      call opt%set_min_objective(f)
      call opt%optimize(x, y, stat)
   end function optimize

   real(pr) function fobj(x, gradient, func_data)
      !! Objective function to fit phase-equilibria points.
      !!
      !! # Description
      !! ...
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  ...
      !! ```
      use yaeos__math, only: sq_error
      real(pr), intent(in) :: x(:)
      real(pr), optional, intent(in out) :: gradient(:)
      class(*), optional, intent(in) :: func_data

      !! Thermodynamic model to make calculations inside
      type(EquilibriaState) :: model_point !! Each solved point
      type(EquilibriaState), pointer :: exp_point

      integer :: i

      real(pr) :: p_exp, t_exp

      fobj = 0
      select type(func_data)
       type is(FittingProblem)

         call func_data%get_model_from_X(X, model)

         do i=1, size(func_data%experimental_points)
            exp_point => func_data%experimental_points(i)
            
            select case(exp_point%kind)
             case("bubble")
               model_point = saturation_pressure(&
                  model, exp_point%x, exp_point%t, kind="bubble", &
                  p0=exp_point%p, y0=exp_point%y &
               )
               fobj = fobj + sq_error(exp_point%p, model_point%p)
            end select
         end do
      end select
      write(1, *) X, fobj
   end function fobj
end module yaeos__fitting
