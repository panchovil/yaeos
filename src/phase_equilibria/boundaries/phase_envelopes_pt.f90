module yaeos__phase_equilibria_boundaries_phase_envelopes_pt
   !! Phase boundaries line on the \(PT\) plane calculation procedures.
   use yaeos_constants, only: pr
   use yaeos_models, only: ArModel
   use yaeos_equilibria_equilibria_state, only: EquilibriaState
   use yaeos_thermoprops, only: fugacity_tp
   use yaeos__math_continuation, only: &
      continuation, continuation_solver, continuation_stopper
   implicit none

   type :: CriticalPoint
      !! Critical point
      real(pr) :: T !! Temperature [K]
      real(pr) :: P !! Pressure [bar]
   end type CriticalPoint

   type :: PTEnvel2
      !! Two-phase isopleth.
      !! Phase boundary line of a fluid at constant composition.
      type(EquilibriaState), allocatable :: points(:)
      !! Each point through the line.
      type(CriticalPoint), allocatable :: cps(:)
      !! Critical points found along the line.
   contains
      procedure, pass :: write =>  write_PTEnvel2
      generic, public :: write (FORMATTED) => write
   end type PTEnvel2

contains

   function pt_envelope_2ph(&
      model, z, first_point, &
      points, iterations, delta_0, specified_variable_0, &
      solver, stop_conditions &
      ) result(envelopes)
      !! PT two-phase envelope calculation procedure.
      !!
      !! Phase envelope calculation using the continuation method.
      !! Defaults to solving the saturation temperature and continues with
      !! an increment in it. The variable to specify can be changed by modifying
      !! `specified_variable_0` with the corresponding variable number.
      ! ========================================================================
      use stdlib_optval, only: optval
      class(ArModel), intent(in) :: model
      !! Thermodyanmic model
      real(pr), intent(in) :: z(:)
      !! Vector of molar fractions
      type(EquilibriaState) :: first_point
      integer, optional, intent(in) :: points
      !! Maxmimum number of points, defaults to 500
      integer, optional, intent(in) :: iterations
      !! Point solver maximum iterations, defaults to 100
      real(pr), optional, intent(in) :: delta_0
      !! Initial extrapolation \(\Delta\)
      integer, optional, intent(in) :: specified_variable_0
      !! Position of specified variable, since the vector of variables is
      !! \(X = [lnK_i, \dots, lnT, lnP]\) the values for specification
      !! will be \([1 \dots nc]\) for the equilibria constants, \(nc+1\) for
      !! \(lnT\) and \(nc + 2\) for \(lnT\).
      procedure(continuation_solver), optional :: solver
      !! Specify solver for each point, defaults to a full newton procedure
      procedure(continuation_stopper), optional :: stop_conditions
      !! Function that returns true if the continuation method should stop
      type(PTEnvel2) :: envelopes
      ! ------------------------------------------------------------------------

      integer :: nc !! Number of components
      integer :: ns !! Number of specified variable
      real(pr) :: dS0 !! Initial specification step
      real(pr) :: S0 !! Initial specification value

      integer :: max_points !! Maximum number of points
      integer :: max_iterations

      real(pr) :: X(size(z) + 2)
      real(pr), allocatable :: XS(:, :)

      character(len=14) :: kind

      ! ========================================================================
      ! Handle input
      ! ------------------------------------------------------------------------
      kind = first_point%kind
      nc = size(z)
      max_points = optval(points, 500)
      max_iterations = optval(iterations, 100)
      ns = optval(specified_variable_0, nc+1)
      dS0 = optval(delta_0, 0.1_pr)


      ! Correctly define the K-values based on the provided incipient point.
      select case(first_point%kind)
       case("bubble", "liquid-liquid")
         X(:nc) = log(first_point%y/z)
       case("dew")
         X(:nc) = 1 - log(first_point%y/z)
      end select

      X(nc+1) = log(first_point%T)
      X(nc+2) = log(first_point%P)
      S0 = X(ns)

      allocate(envelopes%points(0), envelopes%cps(0))
      ! ========================================================================
      ! Trace the line using the continuation method.
      ! ------------------------------------------------------------------------
      XS = continuation(&
         foo, X, ns0=ns, S0=S0, &
         dS0=dS0, max_points=max_points, solver_tol=1.e-9_pr, &
         update_specification=update_specification, &
         stop=stop_conditions &
         )
   contains
      subroutine foo(X, ns, S, F, dF, dFdS)
         !! Function that needs to be solved at each envelope point
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: ns
         real(pr), intent(in) :: S

         real(pr), intent(out) :: F(:)
         real(pr), intent(out) :: dF(:, :)
         real(pr), intent(out) :: dFdS(:)

         real(pr) :: y(nc)
         real(pr) :: Vz, Vy, lnphip_z(nc), lnphip_y(nc)
         real(pr) :: dlnphi_dt_z(nc), dlnphi_dt_y(nc)
         real(pr) :: dlnphi_dp_z(nc), dlnphi_dp_y(nc)
         real(pr) :: dlnphi_dn_z(nc, nc), dlnphi_dn_y(nc, nc)

         real(pr) :: T, P, K(nc)

         integer :: i, j

         F = 0
         dF = 0

         K = exp(X(:nc))
         T = exp(X(nc+1))
         P = exp(X(nc+2))

         y = K*z

         call fugacity_tp(&
            model, z, T, P, V=Vz, root_type="stable", &
            lnphip=lnphip_z, dlnPhidt=dlnphi_dt_z, &
            dlnPhidp=dlnphi_dp_z, dlnphidn=dlnphi_dn_z &
            )
         call fugacity_tp(&
            model, y, T, P, V=Vy, root_type="stable", &
            lnphip=lnphip_y, dlnPhidt=dlnphi_dt_y, &
            dlnPhidp=dlnphi_dp_y, dlnphidn=dlnphi_dn_y &
            )

         F(:nc) = X(:nc) + lnphip_y - lnphip_z
         F(nc + 1) = sum(y - z)
         F(nc + 2) = X(ns) - S

         ! Jacobian Matrix
         do j=1,nc
            df(:nc, j) = dlnphi_dn_y(:, j) * y(j)
            df(j, j) = dF(j, j) + 1
         end do

         df(:nc, nc + 1) = T * (dlnphi_dt_y - dlnphi_dt_z)
         df(:nc, nc + 2) = P * (dlnphi_dp_y - dlnphi_dp_z)

         df(nc + 1, :nc) = y

         df(nc + 2, :) = 0
         df(nc + 2, ns) = 1

         dFdS = 0
         dFdS(nc+2) = -1
      end subroutine foo

      subroutine update_specification(X, ns, S, dS, dXdS, iterations)
         !! Update the specification during continuation.
         real(pr), intent(in out) :: X(:)
         !! Vector of variables \([lnK_i \dots , lnT, lnP]\)
         integer, intent(in out) :: ns
         !! Number of specified variable in the vector
         real(pr), intent(in out) :: S
         !! Variable specification value
         real(pr), intent(in out) :: dS
         !! Step in specification
         real(pr), intent(in out) :: dXdS(:)
         !! Variation of variables with respect to specification
         integer, intent(in) :: iterations
         !! Iterations used in the solver

         real(pr) :: Xnew(nc+2), Xold(nc+2), Xc(nc+2)

         Xold = X

         if (maxval(abs(X(:nc))) < 0.1_pr) then
            ns = maxloc(abs(dXdS(:nc)), dim=1)
         else
            ns = maxloc(abs(dXdS), dim=1)
         end if

         dS = dXdS(ns) * dS
         dXdS = dXdS/dXdS(ns)

         dS = sign(1.0_pr, dS) * minval([ &
            max(sqrt(abs(X(ns))/10._pr), 0.1_pr), &
            abs(dS)*3/iterations &
            ] &
            )

         dS = sign(1.0_pr, dS) * maxval([abs(dS), 0.01_pr])

         ! Save the point
         envelopes%points = [&
            envelopes%points, &
            EquilibriaState(&
            kind=kind, &
            x=z, Vx=0._pr, y=exp(X(:nc))*z, Vy=0, &
            T=exp(X(nc+1)), P=exp(X(nc+2)), beta=0._pr, iters=0) &
            ]
         
         ! Jump over critical point
         do while (maxval(abs(X(:nc))) < 0.05)
            S = S + dS
            X = X + dXdS*dS
         end do

         Xnew = X + dXdS*dS

         if (all(Xold(:nc) * (Xnew(:nc)) < 0)) then
            select case(kind)
             case("dew")
               kind = "bubble"
             case("bubble")
               kind = "dew"
            end select
            Xc = (Xnew - Xold)/(Xnew(ns) - Xold(ns)) * (Xnew(ns)) + (Xold)
            envelopes%cps = [&
               envelopes%cps, CriticalPoint(T=exp(Xc(nc+1)), P=exp(X(nc+2))) &
               ]
         end if
      end subroutine update_specification

   end function pt_envelope_2ph

   subroutine write_PTEnvel2(pt2, unit, iotype, v_list, iostat, iomsg)
      class(PTEnvel2), intent(in) :: pt2
      integer, intent(in) :: unit
      character(*), optional, intent(in) :: iotype
      integer, optional, intent(in)  :: v_list(:)
      integer, intent(out) :: iostat
      character(*), optional, intent(inout) :: iomsg

      integer, allocatable :: cps(:)
      integer :: cp
      integer :: i, nc


      allocate(cps(0))
      do i=1,size(pt2%cps)
         cp = minloc(&
            (pt2%points%T - pt2%cps(i)%T)**2 &
            + (pt2%points%P - pt2%cps(i)%P)**2, dim=1&
            )
         cps = [cps, cp]
      end do

      write(unit, "(A, /, /)") "#PTEnvel2"

      write(unit, "(A, /)") "#" // pt2%points(1)%kind
      
      do i=1, size(pt2%points)
         write(unit, *) pt2%points(i)
         write(unit, "(/)")
         if (any(cps - i == 0)) then
            write(unit, "(/, /)")
            write(unit, "(A, /)") "#" // pt2%points(i+1)%kind
         end if
      end do
      
      write(unit, "(/, /, A, /)") "#Critical"
      do cp = 1, size(cps)
         write(unit, *) pt2%cps(cp)%T, pt2%cps(cp)%P
      end do
   end subroutine write_PTEnvel2

end module yaeos__phase_equilibria_boundaries_phase_envelopes_pt
