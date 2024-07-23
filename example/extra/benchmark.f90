module bench
    use yaeos, only: pr, R, Substances, AlphaSoave, CubicEoS, &
        QMR, PengRobinson76, ArModel
    use hyperdual_pr76, only: PR76, setup_adiff_pr76 => setup
    ! use TapeRobinson, only: setup_taperobinson => setup_model
    implicit none

    real(pr), allocatable :: z(:), tc(:), pc(:), w(:), kij(:,:), lij(:,:)

contains

    subroutine set_vars(n)
        integer, intent(in) :: n
        integer :: i
        z = [(i, i=1,n)]
        z = z/sum(z)
        tc = z * 2
        pc = z * 100
        w  = z/10
        kij = reshape([(i,i=1,n**2)], [n,n])
        lij = kij/2
    end subroutine

    subroutine yaeos__run(n, dn, f_p, model_name)
        integer :: n
        logical :: dn
        logical :: f_p
        character(len=*), intent(in) :: model_name
        class(ArModel), allocatable :: model
        real(pr) :: lnfug(n), dlnphidp(n), dlnphidt(n), dlnphidn(n,n)


        real(pr) :: v, t, p

        call set_vars(n)

        select case(model_name)
        case ("Analytic PR76")
            model = PengRobinson76(tc, pc, w, kij, lij)
        case ("Adiff PR76")
            model = setup_adiff_pr76(tc, pc, w, kij, lij)
        case ("Tape PR76")
            ! model = setup_taperobinson(tc, pc, w, kij, lij)
        end select

        v = 1.0_pr
        t = 150._pr
        p = 15

        if (dn) then
            if (f_p) then
                call model%lnphi_pt(&
                    z, P, T, root_type="stable", &
                    lnPhi=lnfug, dlnphidp=dlnphidp, dlnphidn=dlnphidn)
            else
                call model%lnphi_vt(z, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnphidn)
            end if
        else
            call model%lnphi_vt(z, V, T, lnPhi=lnfug)
        end if
    end subroutine

    subroutine run_bench(nmax, all_derivs, f_p, eos)
        integer, parameter :: nevals=1e3
        integer :: nmax
        logical :: all_derivs
        logical :: f_p
        character(len=*), intent(in) :: eos
        real(pr) :: time, std, mean
        real(pr) :: et, st
        integer :: i, n
        time = 0
        std = 0
        mean = 0

        print *, "running: ", eos, "all derivs: ", all_derivs
        do n=1,nmax
            do i=1,nevals
                call cpu_time(st)
                    call yaeos__run(n, all_derivs, f_p, eos)
                call cpu_time(et)

                time = (et-st)*1e6
                mean = mean + time
                std = std + time**2
            end do
            mean = mean/nevals
            std = sqrt(1._pr/(nevals-1) * abs(std - nevals*mean**2))

            print *, n, mean, std
        end do
    end subroutine
    
    subroutine main()
        use iso_fortran_env, only: compiler_version, compiler_options
        integer :: n=20
        logical :: allderivs=.false.
        logical :: fug_p = .false.

        print *, "============================================================="
        print *, "BENCHMARKING YAEOS"
        print *, "-------------------------------------------------------------"

        print *, "Compiler: ", compiler_version()
        print *, "Compiler options: ", compiler_options()

        print *, "--------------------------------------------------------------"
        print *, "Helmholtz Function valuation"
        call run_bench(n, allderivs, fug_p, "Analytic PR76")
        ! call run_bench(n, allderivs, fug_p, "Tape PR76")
        ! call run_bench(n, allderivs, fug_p, "Adiff PR76")

        allderivs = .true.
        call run_bench(n, allderivs, fug_p, "Analytic PR76")
        ! call run_bench(n, allderivs, fug_p, "Tape PR76")
        ! call run_bench(n, allderivs, fug_p, "Adiff PR76")

        print *, "--------------------------------------------------------------"
        print *, "Flashes"
        flash: block
        use yaeos, only: pr, CubicEoS, PengRobinson76, EquilibriumState, flash
        real(pr) :: tc(2), pc(2), w(2), n(2), T, P
        type(CubicEoS) :: model
        type(EquilibriumState) :: flash_result
        integer :: iter

        ! Methane/ Butane mixture
        n = [0.4, 0.6]                      ! Composition
        tc = [190.564, 425.12]              ! Critical temperatures
        pc = [45.99, 37.96]                 ! Critical pressures
        w = [0.0115478, 0.200164]           ! Acentric factors

        ! Use the PengRobinson76 model
        model = PengRobinson76(tc, pc, w)

        ! Set pressure and temperatures
        P = 60
        T = 294

        ! Calculate flashes
        flash_result = flash(model, n, t=t, p_spec=p, iters=iter)

        end block flash


        print *, "============================================================="
    end subroutine

end module
