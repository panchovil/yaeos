program capilar_envelope
    
    use yaeos, only: pr, &
        SoaveRedlichKwong, PengRobinson76, PengRobinson78, RKPR, &
        EquilibriumState, ArModel, PTEnvel2, &
        pt_envelope_2ph, saturation_pressure, saturation_temperature, k_wilson &  
        , NanoEquilibriumState, nano_pt_envelope_2ph, NanoPTEnvel2
    !use yaeos__equilibria_equilibrium_state
    !use yaeos__equilibria_boundaries_nano_phase_envelopes_pt
        implicit none
    integer, parameter :: nc=8
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

    ! get values
    call values(z, tc, pc, w, kij, lij, r_poro, ang_cont)
    
    ! model definition
    model = SoaveRedlichKwong(tc, pc, w, kij, lij)
    
    ! Calculate a dew point at low pressure to later 
    ! initialize the phase envelope
    sat_point = saturation_temperature(model, z, P=1._pr, kind="dew", t0=150._pr)
    ! Calculate 1 point of bulk envelope
    envelope = pt_envelope_2ph(model, z, sat_point, points=2)
    !write(1,*) envelope%Points(1)
    !write(*,*) envelope%Points(1)%y

    ! Capillay envelope
    allocate(Parachor(nc))
    call Parachor_values(tc, pc, w, Parachor)
    !print*, "Par", Parachor
    call Laplace_init(envelope%points(2)%y, envelope%points(2)%Vx, envelope%points(2)%Vy,&
                     Parachor, IFT_init, Pcap_init)
    
    print*, Pcap_init

    Py_init=envelope%points(2)%P                 
    Px_init=Py_init-Pcap_init

    init_point%kind=envelope%points(2)%kind
    init_point%iters=envelope%points(2)%iters
    init_point%y=envelope%points(2)%y
    init_point%x=envelope%points(2)%x
    init_point%Vy=envelope%points(2)%Vy
    init_point%Vx=envelope%points(2)%Vx
    init_point%T=envelope%points(2)%T
    init_point%Pcap=Pcap_init
    init_point%Py=Py_init
    init_point%Px=Px_init
    init_point%beta=envelope%points(2)%beta
    !print*, init_point
    nano_envelope = nano_pt_envelope_2ph(model, z, r_poro, ang_cont, Parachor, init_point, points=50, iterations=500)
    write(*,*) nano_envelope%points(:)%iters
    write(1,*) nano_envelope

    !print*, Parachor
contains
    subroutine values(z_in, tc_in, pc_in, w_in, kij_in, lij_in, r_poro_in, ang_cont_in)
        real(pr), intent(out) :: tc_in(:), pc_in(:), w_in(:),&
                                z_in(:), kij_in(:,:), lij_in(:,:),&
                                r_poro_in, ang_cont_in
        
        !! names="CO2" "C1-N2" "C2-C3" "C4" "C5" "C6" "C7+n" "Asph"
        !! composition vector                        
        z_in = (/0.0246,0.3694,0.0752,0.0193,0.0157,0.0162,0.47145,0.00815/)
        !! Critical Temperature
        tc_in = (/304.0390,189.4280,339.8720,419.8170,465.0940,507.3170,860.3720,1424.8170/)
        !! Critical Pressure
        pc_in = (/73.7900,45.8300,45.4100,37.5400,33.8000,32.9000,12.4600,12.2900/)
        !! Acentric Factor
        w_in = (/0.225000,0.008500,0.127100,0.187800,0.239700,0.275000,1.022000,1.441000/)
        lij_in = 0
        kij_in = 0
        Kij_in(2,7) = 0.053
        Kij_in(2:5,8) = 0.135
        Kij_in(7,2) = 0.053
        Kij_in(8,2:5) = 0.135
        
        !Capillary pressure variables
        r_poro_in=0.0000001 !radio cualquiera de 100 nm
        ang_cont_in=1.0472 !angulo cualquiera de 60º en radianes
        !ang_cont=ang_cont*3.14/180.0 !la variable esta en º y se necesita en radianes

    end subroutine values
    subroutine Parachor_values(tc_in, pc_in, w_in, Parachor_out)
        real(pr), intent(out) :: Parachor_out(:)
        real(pr), intent(in) :: tc_in(:), pc_in(:), w_in(:)
        !Parachor :: cm^3/mol*(mN/m)^1/4
        Parachor_out= 40.1684*(0.151-0.0464*w_in)*(tc_in**(13.0/12.0))/(pc_in**(5.0/6.0))  ! https://doi.org/10.1002/cjce.5450750617 // eq(12)
    end subroutine Parachor_values
    subroutine Laplace_init(y_in, Vz_in, Vy_in, Par_in, IFT_out, Pcap_out)
        !real(pr), intent(in) :: r_poro_in, ang_cont_in, Par_in(:)
        !real(pr), intent(in) :: Vx_in, Vy_in, y_in(:)
        real(pr), intent(in) :: y_in(:), Vz_in, Vy_in, Par_in(:)
        real(pr), intent(out) :: IFT_out, Pcap_out
        integer :: i
        !Par [cm^3/mol * (mN/m)^1/4], r_poro [m], ang_cont [rad]
        !IFT [(mN/m)^1/4]
        !do i=1,nc
        !    IFT_out=IFT_out+((Par_in(i)/1000.0)*(z(i)/Vz_in-y_in(i)/Vy_in))
            !print*, (Par_in(i)/1000.0)*(z(i)/Vz_in-y_in(i)/Vy_in)
        !end do
        !print*, 0.1*cos(ang_cont)*2.0*((sum((Par_in/1000.0)*(z/Vz_in-y_in/Vy_in)))**4)
        IFT_out = sum((Par_in/1000.0)*(z/Vz_in-y_in/Vy_in))
        !Pcap [bar]
        Pcap_out = (0.00000001*2.0*(IFT_out**4)*cos(ang_cont))/r_poro !E=4
        
    end subroutine Laplace_init

end program capilar_envelope

