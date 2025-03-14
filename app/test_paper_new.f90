program test_paper_new
   !! Program for calculation of phase diagrams. 
   !use forsus, only: Substance, forsus_dir, forsus_default_dir
   use yaeos, only: pr, R, &
      SoaveRedlichKwong, PengRobinson76, PengRobinson78, RKPR, PengRobinson78Nano,&
      EquilibriumState, ArModel, PTEnvel2, &
      pt_envelope_2ph, saturation_pressure, saturation_temperature, k_wilson ,&
      NanoEquilibriumState, nano_pt_envelope_2ph, NanoPTEnvel2
   use hyperdual_pr78_nano
   implicit none

   ! ===========================================================================
   ! Variables definition
   ! ---------------------------------------------------------------------------
   integer, parameter :: nc=5  
   integer :: i
   class(ArModel), allocatable :: model, model_nano ! Thermodynamic model to be used
   type(EquilibriumState) :: sat_point, sat_point_nano!, sat_point_auto_eos            ! Init
   type(NanoEquilibriumState) :: init_point_cap, init_point_nano_cap
   type(PTEnvel2) :: envelope, envelope_nano!, envelope_auto_eos  
   type(NanoPTEnvel2) :: envelope_cap ,envelope_nano_cap
   !class(PR78_nano_autodiff), allocatable :: auto_eos 

                       ! PT Phase envelope
   real(pr) :: tc(nc), pc(nc), w(nc)                     ! Component's critical constants
   real(pr) :: z(nc), kij(nc,nc), lij(nc,nc)             ! Termodynamic variables
   real(pr) :: rp, LJ_par(nc), ang_cont, Parachor(nc)    ! Nano parameters
   real(pr) :: tcc(nc), pcc(nc), delP(nc), delT(nc)      ! Critical shift parameters

   character(len=500) :: header, header_cap, names(nc)  ! Cadena para almacenar el encabezado
   character(len=1000) :: fmt, fmt_cap
   
   ! ===========================================================================
   ! Compound definition
   ! ---------------------------------------------------------------------------
   !! names="C1" "C2-C4" "C5-C7" "C8-C9" "C10+"   
   !! composition vector                        
   z = (/0.2506_pr, 0.22_pr, 0.20_pr, 0.13_pr, 0.1994_pr/)

   z = z/sum(z)
   print*, z
   print*, sum(z)
   !! Critical Pressure
   pc = (/45.40_pr, 42.54_pr, 33.76_pr, 30.91_pr, 21.58_pr/)

   !! Critical Temperature
   tc = (/190.60_pr, 363.30_pr, 511.56_pr, 579.34_pr, 788.74_pr/)

   !! Acentric Factor
   w = (/0.008_pr, 0.1432_pr, 0.2474_pr, 0.2861_pr, 0.6869_pr/)


   lij = 0.0_pr
   kij = 0.0_pr
   kij = reshape([ &
   0.0_pr, 0.0078_pr, 0.0242_pr, 0.0324_pr, 0.0779_pr, &  ! Fila 1
   0.0078_pr, 0.0_pr, 0.0046_pr, 0.0087_pr, 0.0384_pr, &  ! Fila 2
   0.0242_pr, 0.0046_pr, 0.0_pr, 0.0006_pr, 0.0169_pr, &  ! Fila 3
   0.0324_pr, 0.0087_pr, 0.0006_pr, 0.0_pr, 0.0111_pr, &  ! Fila 4
   0.0779_pr, 0.0384_pr, 0.0169_pr, 0.0111_pr, 0.0_pr  &  ! Fila 5
   ], shape=[nc, nc])

   
   rp = 20.0_pr !nm
   ang_cont = 1.0472_pr !rad = 60º
   Parachor = (/77.0_pr, 145.2_pr, 250.0_pr, 306.0_pr, 686.3_pr/)
   !! correlacion con respecto a propiedades criticas
   LJ_par = 0.244_pr*((Tc/Pc)**(1.0/3.0)) !nm

   !! changing tc and pc to tcc and pcc 
   delP = 0.9793_pr*((rp/LJ_par)**(-0.6366_pr))
   delT = 0.7597_pr*((rp/LJ_par)**(-0.7708_pr))
   Pcc = Pc - delP * Pc
   Tcc = Tc - delT * Tc

   ! Model definition
   model = PengRobinson78(tc, pc, w, kij, lij)
   model_nano = PengRobinson78Nano(LJ_par, rp, tc, pc, w, kij, lij)
   !auto_eos = setup(LJ_par=LJ_par, rp=rp, tc=tc, pc=pc, w=w, kij=kij, lij=lij)

   !! Envelope's format
   ! Inicia el encabezado con "T, P"
   names(1) = "C1"
   names(2) = "C2-C4"
   names(3) = "C5-C7"
   names(4) = "C8-C9"
   names(5) = "nC10+"
   header = "         T                     P"
   header_cap = "         T                     Py                     Px                     Pcap"
   ! Añade "lnK1", "lnK2", ..., "lnKn" dinámicamente
   do i = 1, nc
      header = trim(adjustl(header)) // "          lnK" // trim(adjustl(itoa(i))) // "(" // trim(adjustl(names(i))) // ")"
      header_cap = trim(adjustl(header_cap)) // "          lnK" // trim(adjustl(itoa(i))) // "(" // trim(adjustl(names(i))) // ")"
   end do
   ! Define el formato fijo para los datos
   fmt = "(F25.16, F25.16"  ! Espacio fijo para "T" y "P"
   fmt_cap = "(F25.16, F25.16, F25.16, F25.16" ! Espacio fijo para "T" , "Py", "Px" y "Pcap"
   do i = 1, nc
       fmt = trim(fmt) // ", F25.16"  ! Espacio fijo para cada "lnK"
       fmt_cap = trim(fmt_cap) // ", F25.16"  ! Espacio fijo para cada "lnK"
   end do
   fmt = trim(fmt) // ")"
   fmt_cap = trim(fmt_cap) // ")"

   ! ------------------------- BUBBLE ENVELOPE ------------------

   sat_point = saturation_temperature(model, z, P=90._pr, kind="bubble", t0=310._pr)
   sat_point_nano = saturation_temperature(model_nano, z, P=90._pr, kind="bubble", t0=310._pr)
   !sat_point_auto_eos = saturation_temperature(auto_eos, z, P=90._pr, kind="bubble", t0=310._pr)

   ! Laplace's subroutine needs the pore radius in meters
   call Laplace_init(sat_point, Parachor, ang_cont, rp/1E9, init_point_cap)
   call Laplace_init(sat_point_nano, Parachor, ang_cont, rp/1E9, init_point_nano_cap)

   print*, "sat point bulk bubble", sat_point%iters
   print*, "sat point nano bubble", sat_point_nano%iters
   !print*, "sat point nano auto bubble", sat_point_auto_eos%iters

   ! Calculate phase envelope
   envelope = pt_envelope_2ph(model, z, sat_point)
   envelope_nano = pt_envelope_2ph(model_nano, z, sat_point_nano)
   ! pore radius in meters
   !envelope_cap = nano_pt_envelope_2ph(model, z, rp/1E9, ang_cont, Parachor, init_point_cap)
   envelope_nano_cap = nano_pt_envelope_2ph(model_nano, z, rp/1E9, ang_cont, Parachor, init_point_nano_cap, 1000, 1000)
   !envelope_auto_eos = pt_envelope_2ph(auto_eos, z, sat_point_auto_eos)

   !! --------------------------------------------------------------------------
   
   ! write(1, "(A)") trim(header)
   ! do i=1,size(envelope%points)
   !    write(1,fmt) (((envelope%points(i)%T)-273.15_pr)*&
   !    (9.0_pr/5.0_pr))+32.0_pr , 14.5038_pr*envelope%points(i)%P, &
   !    log(envelope%points(i)%y/envelope%points(i)%x)
   ! end do
   do i=1,size(envelope%points)

      write (1,*) envelope%points(i), envelope%points(i)%kind
   end do
   write(*,*)"bulk bubble cps:", envelope%cps
   
   !! --------------------------------------------------------------------------

   ! write(1, "(A)") trim(header)
   ! do i=1,size(envelope_nano%points)
   !    write(1,fmt) (((envelope_nano%points(i)%T)-273.15_pr)*&
   !    (9.0_pr/5.0_pr))+32.0_pr , 14.5038_pr*envelope_nano%points(i)%P, &
   !    log(envelope_nano%points(i)%y/envelope_nano%points(i)%x)
   ! end do
   ! write(*,*)"nano bubble cps:", envelope_nano%cps

   !! --------------------------------------------------------------------------

   ! write(3, "(A)") trim(header)
   ! do i=1,size(envelope_auto_eos%points)
   !    write(3,fmt) (((envelope_auto_eos%points(i)%T)-273.15_pr)*&
   !    (9.0_pr/5.0_pr))+32.0_pr , 14.5038_pr*envelope_auto_eos%points(i)%P, &
   !    log(envelope_auto_eos%points(i)%y/envelope_auto_eos%points(i)%x)
   ! end do
   ! write(*,*)"nano auto bubble cps:", envelope_auto_eos%cps

   !! --------------------------------------------------------------------------
   
   ! write(3, "(A)") trim(header_cap)
   ! do i=1,size(envelope_cap%points)
   !    write(3,fmt_cap) (((envelope_cap%points(i)%T)-273.15_pr)*&
   !    (9.0_pr/5.0_pr))+32.0_pr , 14.5038_pr*envelope_cap%points(i)%Py, &
   !    14.5038_pr*envelope_cap%points(i)%Px, 14.5038_pr*envelope_cap%points(i)%Pcap,&
   !    log(envelope_cap%points(i)%y/envelope_cap%points(i)%x)
   ! end do
   ! write(*,*)"cap bubble cps:", envelope_cap%cps

   !! --------------------------------------------------------------------------

   ! write(3, "(A)") trim(header_cap)
   ! do i=1,size(envelope_nano_cap%points)
   !    write(3,fmt_cap) (((envelope_nano_cap%points(i)%T)-273.15_pr)*&
   !    (9.0_pr/5.0_pr))+32.0_pr , 14.5038_pr*envelope_nano_cap%points(i)%Py, &
   !    14.5038_pr*envelope_nano_cap%points(i)%Px, 14.5038_pr*envelope_nano_cap%points(i)%Pcap,&
   !    log(envelope_nano_cap%points(i)%y/envelope_nano_cap%points(i)%x)
   ! end do
   ! write(*,*)"nano cap bubble cps:", envelope_nano_cap%cps

   write(3,*) envelope_nano_cap


   !! --------------------------------------------------------------------------

   !! -------------------------- DEW ENVELOPE ---------------------------

   sat_point = saturation_temperature(model, z, P=0.5_pr, kind="dew", t0=500._pr)
   sat_point_nano = saturation_temperature(model_nano, z, P=0.5_pr, kind="dew", t0=500._pr)
   !sat_point_auto_eos = saturation_temperature(auto_eos, z, P=0.5_pr, kind="dew", t0=500._pr)
   
   ! Laplace's subroutine needs the pore radius in meters
   call Laplace_init(sat_point, Parachor, ang_cont, rp/1E9, init_point_cap)
   call Laplace_init(sat_point_nano, Parachor, ang_cont, rp/1E9, init_point_nano_cap)

   print*, "sat point bulk dew", sat_point%iters
   print*, "sat point nano dew", sat_point_nano%iters
   !print*, "sat point nano auto dew", sat_point_auto_eos%iters


   envelope = pt_envelope_2ph(model, z, sat_point)
   envelope_nano = pt_envelope_2ph(model_nano, z, sat_point_nano)
   ! pore radius in meters
   !envelope_cap = nano_pt_envelope_2ph(model, z, rp/1E9, ang_cont, Parachor, init_point_cap, 5000)
   envelope_nano_cap = nano_pt_envelope_2ph(model_nano, z, rp/1E9, ang_cont, Parachor, init_point_nano_cap, 900)   
   !envelope_auto_eos = pt_envelope_2ph(auto_eos, z, sat_point)



   !! --------------------------------------------------------------------------   
   
   ! write(2, "(A)") trim(header)
   ! do i=1,size(envelope%points)
   !    write(2,fmt) (((envelope%points(i)%T)-273.15_pr)*&
   !    (9.0_pr/5.0_pr))+32.0_pr , 14.5038_pr*envelope%points(i)%P, &
   !    log(envelope%points(i)%y/envelope%points(i)%x)
   ! end do
   do i=1,size(envelope%points)

      write (2,*) envelope%points(i), envelope%points(i)%kind
   end do
   write(*,*)"bulk dew cps:", envelope%cps
   
   !! --------------------------------------------------------------------------

   ! write(2, "(A)") trim(header)
   ! do i=1,size(envelope_nano%points)
   !    write(2,fmt) (((envelope_nano%points(i)%T)-273.15_pr)*&
   !    (9.0_pr/5.0_pr))+32.0_pr , 14.5038_pr*envelope_nano%points(i)%P, &
   !    log(envelope_nano%points(i)%y/envelope_nano%points(i)%x)
   ! end do
   ! write(*,*)"nano dew cps:", envelope_nano%cps

   !----------------------------------------------------------------------------

   ! write(4, "(A)") trim(header)
   ! do i=1,size(envelope_auto_eos%points)
   !    write(4,fmt) (((envelope_auto_eos%points(i)%T)-273.15_pr)*&
   !    (9.0_pr/5.0_pr))+32.0_pr , 14.5038_pr*envelope_auto_eos%points(i)%P, &
   !    log(envelope_auto_eos%points(i)%y/envelope_auto_eos%points(i)%x)
   ! end do
   ! write(*,*)"nano auto dew cps:", envelope_auto_eos%cps

   !! --------------------------------------------------------------------------
   
   ! write(4, "(A)") trim(header_cap)
   ! do i=1,size(envelope_cap%points)
   !    write(4,fmt_cap) (((envelope_cap%points(i)%T)-273.15_pr)*&
   !    (9.0_pr/5.0_pr))+32.0_pr , 14.5038_pr*envelope_cap%points(i)%Py, &
   !    14.5038_pr*envelope_cap%points(i)%Px, 14.5038_pr*envelope_cap%points(i)%Pcap,&
   !    log(envelope_cap%points(i)%y/envelope_cap%points(i)%x)
   ! end do
   ! write(*,*)"cap bubble cps:", envelope_cap%cps

   !! --------------------------------------------------------------------------

   ! write(4, "(A)") trim(header_cap)
   ! do i=1,size(envelope_nano_cap%points)
   !    write(4,fmt_cap) (((envelope_nano_cap%points(i)%T)-273.15_pr)*&
   !    (9.0_pr/5.0_pr))+32.0_pr , 14.5038_pr*envelope_nano_cap%points(i)%Py, &
   !    14.5038_pr*envelope_nano_cap%points(i)%Px, 14.5038_pr*envelope_nano_cap%points(i)%Pcap,&
   !    log(envelope_nano_cap%points(i)%y/envelope_nano_cap%points(i)%x)
   ! end do
   write(*,*)"nano cap dew cps:", envelope_nano_cap%cps
   write(4,*) envelope_nano_cap
   !! --------------------------------------------------------------------------
   write(*,*) "          ","C1","                       ", "C2-C4","                     ",&
    "C5-C7","                     ","C8-C9","                     ","C10+"   
   write(*,*) LJ_par

   contains

   ! Función para convertir números enteros a cadenas (para el índice i)
   character(len=10) function itoa(i)
     integer, intent(in) :: i
     write(itoa, '(I0)') i
   end function itoa
   ! interface
   subroutine Laplace_init(sat_point_in, Par_in, ang_cont_in, rp_in, init_point_out)
      ! use yaeos, only: pr, EquilibriumState, NanoEquilibriumState
      ! implicit none
      type(EquilibriumState), intent(in) :: sat_point_in
      type(NanoEquilibriumState), intent(out) :: init_point_out
      real(pr), intent(in) :: rp_in, ang_cont_in, Par_in(nc)
      real(pr) :: IFT, Pcap
      !! Calculation of first capillary pressure
      !! Par [cm^3/mol * (mN/m)^1/4], r_poro [m], ang_cont [rad]
      !! IFT [(mN/m)^1/4]
      IFT = sum((Par_in/1E3)*(sat_point_in%x/sat_point_in%Vx-sat_point_in%y/sat_point_in%Vy))
      !! Pcap [bar]
      Pcap = (1E-8*2._pr*(IFT**4)*cos(ang_cont_in))/(rp_in) !E=4
      !! Load the capillary saturation point with the data calculated until now
      init_point_out%kind = sat_point_in%kind
      init_point_out%iters = sat_point_in%iters
      init_point_out%beta = sat_point_in%beta
      init_point_out%ns = sat_point_in%ns
      init_point_out%y = sat_point_in%y
      init_point_out%x = sat_point_in%x
      init_point_out%Vy = sat_point_in%Vy
      init_point_out%Vx = sat_point_in%Vx
      init_point_out%T = sat_point_in%T
      init_point_out%Pcap = Pcap
      init_point_out%Py = sat_point_in%P
      init_point_out%Px = sat_point_in%P-Pcap

   end subroutine Laplace_init
   ! end interface
end program test_paper_new