program ejemplo_gera
    
    use yaeos, only: pr, &
        SoaveRedlichKwong, PengRobinson76, PengRobinson78, RKPR, &
        EquilibriumState, ArModel, PTEnvel2, &
        pt_envelope_2ph, saturation_pressure, saturation_temperature, k_wilson &  
        , NanoEquilibriumState, nano_pt_envelope_2ph, NanoPTEnvel2
    !use yaeos__equilibria_equilibrium_state
    !use yaeos__equilibria_boundaries_nano_phase_envelopes_pt
        implicit none
    integer, parameter :: nc=3
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
    sat_point = saturation_pressure(model, z, P0=700._pr, kind="bubble", t=145._pr)

    ! Calculate 1 point of bulk envelope
    envelope = pt_envelope_2ph(model, z, sat_point, points=200)
    write(1,*) envelope

    ! Capillay envelope
    allocate(Parachor(nc))
    call Parachor_values(tc, pc, w, Parachor)
    call Laplace_init(y_in=sat_point%y, x_in=sat_point%x,& 
                    Vx_in=sat_point%Vx, Vy_in=sat_point%Vy,&
                     Par_in=Parachor, IFT_out=IFT_init, Pcap_out=Pcap_init)
    
    Py_init=sat_point%P                 
    Px_init=Py_init-Pcap_init
    write(*,*) "-----------------Bulk------------------------"
    write(*,*) "iteraciones bulk ",sat_point%iters
    write(*,*) "kind bulk ",sat_point%kind
    write(*,*) "T inicial ",sat_point%T
    write(*,*) "Vx inicial ",sat_point%Vx
    write(*,*) "Vy inicial ",sat_point%Vy
    write(*,*) "K inicial ",(sat_point%y/sat_point%x)
    write(*,*) "Pcap inicial ",Pcap_init
    write(*,*) "Px inicial ",Px_init
    write(*,*) "Py inicial ",Py_init

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

    write(3,*) "nano"
    nano_envelope = nano_pt_envelope_2ph(model, z, r_poro, ang_cont,&
                     Parachor, init_point, points=2000, iterations=1000)
    write(2,*) nano_envelope
    write(*,*) "-----------------Nano------------------------"
    write(*,*) "iteraciones Nano ",nano_envelope%points(1)%iters
    write(*,*) "T final ", nano_envelope%points(1)%T
    write(*,*) "Vx final ", nano_envelope%points(1)%Vx
    write(*,*) "Vy final ", nano_envelope%points(1)%Vy
    write(*,*) "K final ", (nano_envelope%points(1)%y/nano_envelope%points(1)%x)
    write(*,*) "Pcap final ", nano_envelope%points(1)%Pcap
    write(*,*) "Px final ", nano_envelope%points(1)%Px
    write(*,*) "Py final ", nano_envelope%points(2)%Py
    !write(cadena,"(F16.8)") nano_envelope%points(1)%Py
    !write(*,*) "Py final ", cadena

   

    !sat_point = saturation_pressure(model, z, P0=900._pr, kind="liquid-liquid", t=500._pr)




contains
    subroutine values(z_in, tc_in, pc_in, w_in, kij_in, lij_in, r_poro_in, ang_cont_in)
        real(pr), intent(out) :: tc_in(:), pc_in(:), w_in(:),&
                                z_in(:), kij_in(:,:), lij_in(:,:),&
                                r_poro_in, ang_cont_in
        
        !! names="CH4" "C50-C60" "Asph"
        !! composition vector                        
        z_in = (/0.621476711621686, 0.294896838294686, 0.0836264500836264/)
        !! Critical Temperature
        tc_in = (/190.56, 956.95, 1286.75/)
        !! Critical Pressure
        pc_in = (/45.99, 13.72, 18.11/)
        !! Acentric Factor
        w_in = (/0.0115, 1.313, 1.274/)
        lij_in = 0
        kij_in = 0
        Kij_in(1,3) = 0.017
        Kij_in(3,1) = 0.017

        
        !Capillary pressure variables
        r_poro_in=1E-7 !radio cualquiera de 100 nm
        ang_cont_in=1.0472 !angulo cualquiera de 60º en radianes
        !ang_cont=ang_cont*3.14/180.0 !la variable esta en º y se necesita en radianes

    end subroutine values
    subroutine Parachor_values(tc_in, pc_in, w_in, Parachor_out)
        real(pr), intent(out) :: Parachor_out(:)
        real(pr), intent(in) :: tc_in(:), pc_in(:), w_in(:)
        !Parachor :: [cm^3/mol * (mN/m)^1/4]
        !Parachor_out= 40.1684*(0.151-0.0464*w_in)*(tc_in**(13._pr/12._pr))/(pc_in**(5._pr/6._pr))  ! https://doi.org/10.1002/cjce.5450750617 // eq(12)
        Parachor_out = (/0.07405, 2.76, 3.2/)*1E3

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

end program ejemplo_gera
