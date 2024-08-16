module yaeos__equilibria_boundaries_nano_phase_envelopes_pt
    !! Phase boundaries line on the \(PT\) plane calculation procedures.
    use yaeos__constants, only: pr
    use yaeos__models, only: ArModel
    use yaeos__equilibria_equilibrium_state, only: NanoEquilibriumState
    use yaeos__math_continuation, only: &
       continuation, continuation_solver, continuation_stopper
    implicit none
 
    type :: NanoCriticalPoint
       !! Critical point
       real(pr) :: T !! Temperature [K]
       real(pr) :: P !! Pressure [bar]
    end type NanoCriticalPoint
 
    type :: NanoPTEnvel2
       !! Two-phase isopleth.
       !! Phase boundary line of a fluid at constant composition.
       type(NanoEquilibriumState), allocatable :: points(:)
       !! Each point through the line.
       type(NanoCriticalPoint), allocatable :: cps(:)
       !! Critical points found along the line.
    contains
       procedure, pass :: write =>  nano_pt_envelope_2ph
       generic, public :: write (FORMATTED) => write
    end type NanoPTEnvel2
 
    ! Saved volume values
    real(pr), private :: Vz
    real(pr), private :: Vy
    ! Saved Laplace values
    real(pr), private :: Pcap, IFT

 
contains
    function nano_pt_envelope_2ph(&
        model, z, r_poro, ang_cont, Par, first_point, &
        points, iterations, delta_0, specified_variable_0, &
        solver, stop_conditions &
        ) result(nano_envelopes)
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
        real(pr), intent(in) :: r_poro, ang_cont, Par(:)
        ! Capilar variables
        type(NanoEquilibriumState) :: first_point
        integer, optional, intent(in) :: points
        !! Maxmimum number of points, defaults to 500
        integer, optional, intent(in) :: iterations
        !! Point solver maximum iterations, defaults to 100
        real(pr), optional, intent(in) :: delta_0
        !! Initial extrapolation \(\Delta\)
        integer, optional, intent(in) :: specified_variable_0
        !! Position of specified variable, since the vector of variables is
        !! \(X = [lnK_i, \dots, lnT, Px, Py]\) the values for specification
        !! will be \([1 \dots nc]\) for the equilibria constants, \(nc+1\) for
        !! \(lnT\), \(nc + 2\) for \(Px\) and \(nc + 3\) for \(Py\).
        procedure(continuation_solver), optional :: solver
        !! Specify solver for each point, defaults to a full newton procedure
        procedure(continuation_stopper), optional :: stop_conditions
        !! Function that returns true if the continuation method should stop
        type(NanoPTEnvel2) :: nano_envelopes
        ! ------------------------------------------------------------------------

        integer :: nc !! Number of components
        integer :: ns !! Number of specified variable
        real(pr) :: dS0 !! Initial specification step
        real(pr) :: S0 !! Initial specification value

        integer :: max_points !! Maximum number of points
        integer :: max_iterations !! Maximum number of iterations

        real(pr) :: X(size(z) + 3)
        !! Vector of variables used in the continuation method
        real(pr), allocatable :: XS(:, :)
        !! All the calculated variables that are returned on the continuation
        !! method procedure (unused since each point is saved on the fly)

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
        X(:nc) = log(first_point%x/z)
        end select

        X(nc+1) = log(first_point%T)
        X(nc+2) = first_point%Px
        X(nc+3) = first_point%Py
        S0 = X(ns)

        allocate(nano_envelopes%points(0), nano_envelopes%cps(0))
        ! ========================================================================
        ! Trace the line using the continuation method.
        ! ------------------------------------------------------------------------
        XS = continuation(&
        foo, X, ns0=ns, S0=S0, &
        dS0=dS0, max_points=max_points, solver_tol=1.e-9_pr, &
        update_specification=update_spec, &
        solver=solver, stop=stop_conditions &
        )

    contains
        subroutine Laplace(y_in, IFT_out, Pcap_out)
            !real(pr), intent(in) :: r_poro_in, ang_cont_in, Par_in(:)
            !real(pr), intent(in) :: Vx_in, Vy_in, y_in(:)
            real, intent(in) :: y_in(:)
            real(pr), intent(out) :: IFT_out, Pcap_out
            integer :: i
            !Par [cm^3/mol * (mN/m)^1/4], r_poro [m], ang_cont [rad]
            !IFT [(mN/m)^1/4]
            do i=1,nc
                IFT_out=IFT_out+((Par(i)/1000.0)*(z(i)/Vx-y_in(i)/Vy))
            end do
            !Pcap [bar]
            Pcap_out = (0.00000001*2.0*(IFT_out**4)*cos(ang_cont))/r_poro !E=4
        end subroutine Laplace
        
        subroutine foo(X, ns, S, F, dF, dFdS)   
            !! Function that needs to be solved at each envelope point
            real(pr), intent(in) :: X(:)
            integer, intent(in) :: ns
            real(pr), intent(in) :: S

            real(pr), intent(out) :: F(:)
            real(pr), intent(out) :: dF(:, :)
            real(pr), intent(out) :: dFdS(:)

            character(len=14) :: kind_z, kind_y

            real(pr) :: y(nc)
            real(pr) :: lnFug_z(nc), lnFug_y(nc) 
            real(pr) :: dlnphi_dt_z(nc), dlnphi_dt_y(nc)
            real(pr) :: dlnphi_dp_z(nc), dlnphi_dp_y(nc)
            real(pr) :: dlnphi_dn_z(nc, nc), dlnphi_dn_y(nc, nc)

            real(pr) :: T, Py, Pz, K(nc)

            integer :: i, j

            !! Capilar Preasure variables
            !real(pr), intent(in) :: IFT, r_poro, ang_cont, Pcap, Par(nc) 
            real(pr) :: dPdV_y, dPdV_z, dVdT_y, dVdT_z
            real(pr) :: dVdn_y(nc)
        
            !! Capilar jacobian intermediate variables
            real(pr) :: var_dFn2, var_dFn2_dK, var_dFn2_dT, var_dFn2_dPz, var_dFn2_dPy

            F = 0
            dF = 0
            
            K = exp(X(:nc))
            T = exp(X(nc+1))
            Pz = X(nc+2)
            Py = X(nc+3)
            
            y = K*z
            select case(kind)
            case ("bubble")
                kind_z = "liquid"
                kind_y = "vapor"
            case ("dew")
                kind_z = "vapor"
                kind_y = "liquid"
            case default
                kind_z = "stable"
                kind_y = "stable"
            end select   

            call model%lnphi_pt(&
                z, P=Pz, T, V=Vz, root_type=kind_z, &
                lnFug=lnFug_z, dlnPhidt=dlnphi_dt_z, &
                dlnPhidp=dlnphi_dp_z, dlnphidn=dlnphi_dn_z, &
                dPdV=dPdV_z, dVdT=dVdT_z)
            call model%lnphi_pt(&
                y, P=Py, T, V=Vy, root_type=kind_y, &
                lnFug=lnFug_y, dlnPhidt=dlnphi_dt_y, &
                dlnPhidp=dlnphi_dp_y, dlnphidn=dlnphi_dn_y, &
                dPdV=dPdV_y, dVdT=dVdT_y, dVdn=dVdn_y)
            call Laplace(y_in=y, IFT_out=IFT, Pcap_out= Pcap)

            F(:nc) = X(:nc) + lnFug_y - lnFug_z
            F(nc + 1) = sum(y - z)
            F(nc + 2) = Pz - Py + Pcap
            F(nc + 3) = X(ns) - S
            
            !! Jacobian intermediate variables
            var_dFn2 = 1.0E-11*(8._pr*cos(ang_cont)/r_poro)*(IFT**3) !! 1.0E-11 is an unit conversion 
            do i=1,nc
                var_dFn2_dK = var_dFn2_dK+(Par(i)*((y(i)*dVdn_y(i)/(Vy**2))-(1._pr/Vy)))
                var_dFn2_dT = var_dFn2_dT+(Par(i)*((y(i)*dVdT_y/(Vy**2))-(z(i)*dVdT_z/(Vz**2))))
                var_dFn2_dPz = var_dFn2_dPz + (-Par(i)*z(i))
                var_dFn2_dPy = var_dFn2_dPy + (Par(i)*y(i))
            end do
            
            !! Jacobian Matrix
            do j=1,nc
                df(:nc, j) = dlnphi_dn_y(:, j) * y(j)
                df(j, j) = dF(j, j) + 1._pr
            end do

            df(:nc, nc + 1) = T * (dlnphi_dt_y - dlnphi_dt_z)
            df(:nc, nc + 2) = -(1._pr/Pz) - (dlnphi_dp_z)
            df(:nc, nc + 3) = (1._pr/Py) + (dlnphi_dp_y)
      
            df(nc + 1, :nc) = y
            
            df(nc + 2, :nc) = y*var_dFn2*var_dFn2_dK
            df(nc + 2, nc + 1) = var_dFn2*T*var_dFn2_dT
            df(nc + 2, nc + 2) = 1._pr + var_dFn2*var_dFn2_dPz*(1._pr/(dPdV_z*(Vz**2))) 
            df(nc + 2, nc + 3) = - 1._pr + var_dFn2*var_dFn2_dPy*(1._pr/(dPdV_y*(Vy**2)))
    
            df(nc + 3, :) = 0._pr
            df(nc + 3, ns) = 1._pr

            dFdS = 0._pr
            dFdS(nc + 2) = -1._pr
        end subroutine foo
        
        subroutine update_spec(X, ns, S, dS, dXdS, step_iters)
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
            integer, intent(in) :: step_iters
            !! Iterations used in the solver
   
            real(pr) :: maxdS
   
            ! =====================================================================
            ! Update specification
            ! - Dont select T or P near critical points
            ! - Update dS wrt specification units
            ! - Set step
            ! ---------------------------------------------------------------------
            if (maxval(abs(X(:nc))) < 0.1_pr) then
               ns = maxloc(abs(dXdS(:nc)), dim=1)
               maxdS=0.01_pr
            else
               ns = maxloc(abs(dXdS), dim=1)
               maxdS = 0.05_pr
            end if
   
            dS = dXdS(ns) * dS
            dXdS = dXdS/dXdS(ns)
   
            dS = sign(1.0_pr, dS) * minval([ &
               max(sqrt(abs(X(ns))/10._pr), 0.1_pr), &
               abs(dS)*3/step_iters &
               ] &
               )
   
            dS = sign(1.0_pr, dS) * maxval([abs(dS), maxdS])
   
            call save_point(X, step_iters)
            call detect_critical(X, dXdS, ns, S, dS)
        end subroutine update_spec 
        
        subroutine save_point(X, iters)
            !! Save the converged point
            real(pr), intent(in) :: X(:)
            integer, intent(in) :: iters
            type(NanoEquilibriumState) :: point
   
            real(pr) :: y(nc), T, Pz, Py
            
            y = exp(X(:nc))*z
            T = exp(X(nc + 1))
            Pz = (X(nc + 2))
            Py = (X(nc + 3))

            select case(kind)
                case("bubble")
                point = EquilibriumState(&
                    kind="bubble", x=z, Vx=Vz, y=y, Vy=Vy, &
                    T=T, Px=Pz, Py=Py, Pcap=Pcap, beta=0._pr, iters=iters &
                    )
                case("dew")
                point = EquilibriumState(&
                    kind="dew", x=y, Vx=Vy, y=z, Vy=Vz, &
                    T=T, Px=Py, Py=Pz, Pcap=Pcap, beta=1._pr, iters=iters &
                    )
                case default
                point = EquilibriumState(&
                    kind=kind, x=z, Vx=Vz, y=y, Vy=Vy, &
                    T=T, Px=Pz, Py=Py, Pcap=Pcap, beta=0._pr, iters=iters &
                    )
            end select
   
            nano_envelopes%points = [nano_envelopes%points, point]
        end subroutine save_point
        
        subroutine detect_critical(X, dXdS, ns, S, dS)
            !! # `detect_critical`
            !! Critical point detection
            !!
            !! # Description
            !! If the values of lnK (X[:nc]) change sign then a critical point
            !! Has passed, since for this to happen all variables should pass
            !! through zero. Near critical points (lnK < 0.05) points are harder
            !! to converge, so more steps in the extrapolation vector are made to
            !! jump over the critical point.
            !! If the critical point is detected then the kind of the point is
            !! changed and the point is saved using an interpolation knowing that
            !!
            !! \[
            !!   X_c = a * X + (1-a)*X_{new}
            !! \]
            !!
            !! With \(X_c\) is the variables at the critical point, \(X_{new}\)
            !! is the new initialization point of the method and \(a\) is the
            !! parameter to interpolate the values. This subroutine finds the
            !! value of  \(a\) to obtain \(X_c\).
            real(pr), intent(in out) :: X(:) !! Vector of variables
            real(pr), intent(in out) :: dXdS(:) !! Variation of variables wrt S
            integer, intent(in out) :: ns !! Number of specified variable
            real(pr), intent(in out) :: S !! Specification value
            real(pr), intent(in out) :: dS !! Step in specification
            real(pr) :: Xc(nc+3) !! Value at (near) critical point
            real(pr) :: a !! Parameter for interpolation
   
            real(pr) :: Xold(size(X)) !! Old value of X
            real(pr) :: Xnew(size(X)) !! Value of the next initialization
   
            Xold = X
   
            do while (maxval(abs(X(:nc))) < 0.05)
               ! If near a critical point, jump over it
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
                case default
                  kind = "liquid-liquid"
               end select
   
               ! 0 = a*X(ns) + (1-a)*Xnew(ns) < Interpolation equation to get X(ns) = 0
               a = -Xnew(ns)/(X(ns) - Xnew(ns))
               Xc = a * X + (1-a)*Xnew
   
               nano_envelopes%cps = [&
                  nano_envelopes%cps, NanoCriticalPoint(T=exp(Xc(nc+1)), P=Xc(nc+3)) &
                  ]
               X = Xc + dXdS*dS
            end if
        end subroutine detect_critical
    end function nano_pt_envelope_2ph
    subroutine write_NanoPTEnvel2(pt2, unit, iotype, v_list, iostat, iomsg)
        class(NanoPTEnvel2), intent(in) :: pt2
        integer, intent(in) :: unit
        character(*),  intent(in) :: iotype
        integer, intent(in)  :: v_list(:)
        integer, intent(out) :: iostat
        character(*), intent(inout) :: iomsg
  
        integer, allocatable :: cps(:)
        integer :: cp
        integer :: i, nc
  
  
        if (size(pt2%points) == 0) return
        allocate(cps(0))
        do i=1,size(pt2%cps)
           cp = minloc(&
              (pt2%points%T - pt2%cps(i)%T)**2 &
              + (pt2%points%Py - pt2%cps(i)%Py)**2, dim=1&
              )
           cps = [cps, cp]
        end do
  
        write(unit,  "(A, /, /)", iostat=iostat) "#NanoPTEnvel2"
  
        write(unit, "(A, /)") "#" // pt2%points(1)%kind
  
        do i=1, size(pt2%points)-1
           ! Change label if passed a critical point
           if (any(cps - i == 0) .and. i < size(pt2%points)) then
              write(unit, "(/, /)")
              write(unit, "(A, /)") "#" // pt2%points(i+1)%kind
           end if
  
           write(unit, *) pt2%points(i)
           write(unit, "(/)")
        end do
  
        write(unit, "(/, /, A, /)") "#Critical"
        do cp = 1, size(cps)
           write(unit, *) pt2%cps(cp)%T, pt2%cps(cp)%P
        end do
    end subroutine write_NanoPTEnvel2
end module yaeos__equilibria_boundaries_nano_phase_envelopes_pt