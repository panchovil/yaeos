program capillary_example_sandoval
    
    use yaeos, only: pr, &
        SoaveRedlichKwong, PengRobinson76, PengRobinson78, RKPR, &
        EquilibriumState, ArModel, PTEnvel2, &
        pt_envelope_2ph, saturation_pressure, saturation_temperature, k_wilson &  
        , NanoEquilibriumState, nano_pt_envelope_2ph, NanoPTEnvel2
    !use yaeos__equilibria_equilibrium_state
    !use yaeos__equilibria_boundaries_nano_phase_envelopes_pt
        implicit none
    integer, parameter :: nc=7
    class(ArModel), allocatable :: model ! Thermodynamic model to be used
    type(EquilibriumState) :: sat_point   ! Init bulk
    type(NanoEquilibriumState) :: init_point ! Init Nano

    type(PTEnvel2) :: envelope           ! PT Phase envelope
    type(NanoPTEnvel2) :: nano_envelope  ! PT Phase nano envelope
    real(pr) :: tc(nc), pc(nc), w(nc), z(nc), kij(nc,nc), lij(nc,nc) ! Component's values
    ! Capilar's values
    real(pr) :: r_poro, ang_cont
    real(pr), allocatable :: Parachor(:)
    real(pr) :: Pcap_init, IFT_init, Py_init, Px_init
    character(len=99) :: cadena

    ! get values
    call values(z, tc, pc, w, kij, lij, r_poro, ang_cont)
    
    ! model definition
    model = SoaveRedlichKwong(tc, pc, w, kij, lij)
    
    !! ------------------------- BUBBLE ENVELOPE ------------------
    ! initialize the phase envelope
    sat_point = saturation_pressure(model, z, P0=50.5_pr, kind="bubble", t=300._pr)
    !sat_point = saturation_temperature(model, z, P=50.5_pr, kind="bubble", t0=300._pr)
    ! Calculate 1 point of bulk envelope
    envelope = pt_envelope_2ph(model, z, sat_point, points=350, iterations=500)
    write(1,*) envelope

    ! Capillay envelope
    allocate(Parachor(nc))
    call Parachor_values(tc, pc, w, Parachor)
    call Laplace_init(y_in=sat_point%y, x_in=sat_point%x,& 
                    Vx_in=sat_point%Vx, Vy_in=sat_point%Vy,&
                     Par_in=Parachor, IFT_out=IFT_init, Pcap_out=Pcap_init)
    
    Py_init=sat_point%P                 
    Px_init=Py_init-Pcap_init
    !write(1,*) "-----------------Bulk------------------------"
    !write(1,*) "iteraciones bulk ",sat_point%iters
    !write(1,*) "T inicial ",sat_point%T
    !write(1,*) "Vx inicial ",sat_point%Vx
    !write(1,*) "Vy inicial ",sat_point%Vy
    !write(1,*) "K inicial ",(sat_point%y/sat_point%x)
    !write(1,*) "Pcap inicial ",Pcap_init
    !write(1,*) "Px inicial ",Px_init
    !write(1,*) "Py inicial ",Py_init

    init_point%kind=sat_point%kind
    init_point%iters=sat_point%iters
    init_point%y=sat_point%y
    init_point%x=sat_point%x
    init_point%Vy=sat_point%Vy
    init_point%Vx=sat_point%Vx
    init_point%T=sat_point%T
    init_point%Pcap=Pcap_init
    init_point%Py=Py_init
    init_point%Px=Px_init
    init_point%beta=sat_point%beta

    !nano_envelope = nano_pt_envelope_2ph(model, z, r_poro, ang_cont,&
    !                 Parachor, init_point, points=3000)
    !write(2,*) nano_envelope
    !write(1,*) "-----------------Nano------------------------"
    !write(1,*) "iteraciones Nano ",nano_envelope%points(1)%iters
    !write(1,*) "T final ", nano_envelope%points(1)%T
    !write(1,*) "Vx final ", nano_envelope%points(1)%Vx
    !write(1,*) "Vy final ", nano_envelope%points(1)%Vy
    !write(1,*) "K final ", (nano_envelope%points(1)%y/nano_envelope%points(1)%x)
    !write(1,*) "Pcap final ", nano_envelope%points(1)%Pcap
    !write(1,*) "Px final ", nano_envelope%points(1)%Px
    !write(1,*) "Py final ", nano_envelope%points(1)%Py
    !write(cadena,"(F16.8)") nano_envelope%points(1)%Py
    !write(*,*) "Py final ", cadena

    !! -------------------------- DEW ENVELOPE ---------------------------

    ! sat_point = saturation_temperature(model, z, P=0.5_pr, kind="dew", t0=200._pr)
    ! !sat_point = saturation_pressure(model, z, P0=0.5_pr, kind="dew", t=300._pr)


    ! !! Calculate 1 point of bulk envelope
    ! envelope = pt_envelope_2ph(model, z, sat_point, points=100)
    ! write(3,*) envelope

    ! call Laplace_init(y_in=sat_point%y, x_in=sat_point%x,& 
    ! Vx_in=sat_point%Vx, Vy_in=sat_point%Vy,&
    !  Par_in=Parachor, IFT_out=IFT_init, Pcap_out=Pcap_init)

    ! Py_init=sat_point%P                 
    ! Px_init=Py_init-Pcap_init
    ! init_point%kind=sat_point%kind
    ! init_point%iters=sat_point%iters
    ! init_point%y=sat_point%y
    ! init_point%x=sat_point%x
    ! init_point%Vy=sat_point%Vy
    ! init_point%Vx=sat_point%Vx
    ! init_point%T=sat_point%T
    ! init_point%Pcap=Pcap_init
    ! init_point%Py=Py_init
    ! init_point%Px=Px_init
    ! init_point%beta=sat_point%beta
    ! nano_envelope = nano_pt_envelope_2ph(model, z, r_poro, ang_cont,&
    !                  Parachor, init_point, points=500)
    ! write(4,*) nano_envelope
    ! print*, 1e-2, 0.01

    ! !sat_point = saturation_pressure(model, z, P0=900._pr, kind="liquid-liquid", t=500._pr)




contains
    subroutine values(z_in, tc_in, pc_in, w_in, kij_in, lij_in, r_poro_in, ang_cont_in)
        real(pr), intent(out) :: tc_in(:), pc_in(:), w_in(:),&
                                z_in(:), kij_in(:,:), lij_in(:,:),&
                                r_poro_in, ang_cont_in
        
        !! names="N2" "CH4" "C2H6" "C3H8" "nC4" "nC5" "nC6"
        !! composition vector                        
        z_in = (/0.014, 0.943, 0.027, 0.0074, 0.0049, 0.0027, 0.001/)
        !! Critical Temperature
        tc_in = (/34.0, 45.99, 48.72, 42.48, 37.96, 33.70, 30.25/)
        !! Critical Pressure
        pc_in = (/126.20, 190.56, 305.32, 369.83, 425.12, 469.70, 507.60/)
        !! Acentric Factor
        w_in = (/0.0377, 0.0115, 0.0995, 0.1523, 0.2002, 0.2515, 0.3013/)
        lij_in = 0
        kij_in = 0
        Kij_in(1,2) = 0.0278
        Kij_in(1,3) = 0.4070
        Kij_in(1,4) = 0.0763
        Kij_in(1,5) = 0.0700
        Kij_in(1,6) = 0.0787
        Kij_in(1,7) = 0.1496

        Kij_in(2,1) = 0.0278
        Kij_in(2,3) = -0.0078
        Kij_in(2,4) = 0.0090
        Kij_in(2,5) = 0.0056
        Kij_in(2,6) = 0.0190
        Kij_in(2,7) = 0.0374



        
        !Capillary pressure variables
        r_poro_in=1E-8 !radio cualquiera de 10 nm
        ang_cont_in=0 !angulo cualquiera de 0º en radianes
        !ang_cont=ang_cont*3.14/180.0 !la variable esta en º y se necesita en radianes

    end subroutine values
    subroutine Parachor_values(tc_in, pc_in, w_in, Parachor_out)
        real(pr), intent(out) :: Parachor_out(:)
        real(pr), intent(in) :: tc_in(:), pc_in(:), w_in(:)
        !Parachor :: [cm^3/mol * (mN/m)^1/4]
        !Parachor_out= 40.1684*(0.151-0.0464*w_in)*(tc_in**(13._pr/12._pr))/(pc_in**(5._pr/6._pr))  ! https://doi.org/10.1002/cjce.5450750617 // eq(12)
        Parachor_out = (/61.12, 74.05, 112.91, 154.03, 193.90, 236.00, 276.71/)
    end subroutine Parachor_values
    subroutine Laplace_init(y_in, x_in, Vx_in, Vy_in, Par_in, IFT_out, Pcap_out)
        !real(pr), intent(in) :: r_poro_in, ang_cont_in, Par_in(:)
        !real(pr), intent(in) :: Vx_in, Vy_in, y_in(:)
        real(pr), intent(in) :: y_in(:), x_in(:), Vx_in, Vy_in, Par_in(:)
        real(pr), intent(out) :: IFT_out, Pcap_out
        integer :: i
        !! Par [cm^3/mol * (mN/m)^1/4], r_poro [m], ang_cont [rad]
        !! IFT [(mN/m)^1/4]
        IFT_out = sum((Par_in/1E3)*(x_in/Vx_in-y_in/Vy_in))
        !! Pcap [bar]
        Pcap_out = (1E-8*2._pr*(IFT_out**4)*cos(ang_cont))/r_poro !E=4
        
    end subroutine Laplace_init

end program capillary_example_sandoval
