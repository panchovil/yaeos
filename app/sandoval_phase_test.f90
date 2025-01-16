program sandoval_phase_test
   !! Program for calculation of phase diagrams. 
   !use forsus, only: Substance, forsus_dir, forsus_default_dir
   use yaeos, only: pr, &
      SoaveRedlichKwong, PengRobinson76, PengRobinson78, RKPR, PengRobinson78Nano,&
      EquilibriumState, ArModel, PTEnvel2, &
      pt_envelope_2ph, saturation_pressure, saturation_temperature, k_wilson 
   use hyperdual_pr78_nano
   implicit none

   ! ===========================================================================
   ! Variables definition
   ! ---------------------------------------------------------------------------
   integer, parameter :: nc=7  
   integer :: i
   class(ArModel), allocatable :: model, model_nano ! Thermodynamic model to be used
   type(EquilibriumState) :: sat_point, sat_point_nano, sat_point_auto_eos              ! Init
   type(PTEnvel2) :: envelope, envelope_nano, envelope_auto_eos  
   class(PR78_nano_autodiff), allocatable :: auto_eos 
                       ! PT Phase envelope
   real(pr) :: tc(nc), pc(nc), w(nc)                ! Component's critical constants
   real(pr) :: z(nc), kij(nc,nc), lij(nc,nc)        ! Termodynamic variables
   real(pr) :: rp, LJ_par(nc)                       ! Nano parameters
   real(pr) :: pcc(nc), delP(nc), tcc(nc), delT(nc)

   character(len=500) :: header, names(nc)  ! Cadena para almacenar el encabezado
   character(len=1000) :: fmt
   
   ! ===========================================================================
   ! Compound definition
   ! ---------------------------------------------------------------------------
   !! names="N2" "CH4" "C2H6" "C3H8" "nC4" "nC5" "nC6"
   
   !! composition vector                        
   z = (/0.0014_pr, 0.943_pr, 0.027_pr, 0.0074_pr, 0.0049_pr, 0.0027_pr, 0.001_pr/)
   !z = (/0.0014_pr, 0.943_pr, 0.027_pr, 0.0074_pr, 0.0049_pr/)

   z = z/sum(z)
   print*, z
   print*, sum(z)
   !! Critical Pressure
   pc = (/34.0_pr, 45.99_pr, 48.72_pr, 42.48_pr, 37.96_pr, 33.70_pr, 30.25_pr/)
   !pc = (/34.0_pr, 45.99_pr, 48.72_pr, 42.48_pr, 37.96_pr/)

   !! Critical Temperature
   tc = (/126.20_pr, 190.56_pr, 305.32_pr, 369.83_pr, 425.12_pr, 469.70_pr, 507.60_pr/)
   !tc = (/126.20_pr, 190.56_pr, 305.32_pr, 369.83_pr, 425.12_pr/)

   !! Acentric Factor
   w = (/0.0377_pr, 0.0115_pr, 0.0995_pr, 0.1523_pr, 0.2002_pr, 0.2515_pr, 0.3013_pr/)
   !w = (/0.0377_pr, 0.0115_pr, 0.0995_pr, 0.1523_pr, 0.2002_pr/)

   lij = 0.0_pr
   kij = 0.0_pr
   Kij(1,2) = 0.0278_pr
   Kij(1,3) = 0.4070_pr
   Kij(1,4) = 0.0763_pr
   Kij(1,5) = 0.0700_pr
   Kij(1,6) = 0.0787_pr
   Kij(1,7) = 0.1496_pr

   Kij(2,1) = 0.0278_pr
   Kij(2,3) = -0.0078_pr
   Kij(2,4) = 0.0090_pr
   Kij(2,5) = 0.0056_pr
   Kij(2,6) = 0.0190_pr
   Kij(2,7) = 0.0374_pr
   
   rp = 50.0_pr !nm
   LJ_par = (/0.364_pr, 0.375_pr, 0.443_pr, 0.472_pr, 0.506_pr, 0.529_pr, 0.550_pr/) 
!   LJ_par = (/0.364_pr, 0.375_pr, 0.443_pr, 0.472_pr, 0.506_pr/) 

   ! !! Critical Pressure shift
   ! delP = 0.9793*((rp/LJ_par)**-0.6366)
   ! !! Critical Pressure modified
   ! pcc = -((delP*pc)-pc)
   ! !! Critical Temperature shift
   ! delT = 0.7597*((rp/LJ_par)**-0.7708)
   ! !! Critical Pressure modified
   ! tcc = -((delT*tc)-tc)

   ! Model definition
   model = PengRobinson78(tc, pc, w, kij, lij)
   model_nano = PengRobinson78Nano(LJ_par, rp, tc, pc, w, kij, lij)
   auto_eos = setup(LJ_par=LJ_par, rp=rp, tc=tc, pc=pc, w=w, kij=kij, lij=lij)
   
   ! Calculate a dew point at low pressure to later 
   ! initialize the phase envelope
   ! sat_point = saturation_temperature(model, z, P=1._pr, kind="dew", t0=300._pr)
   ! sat_point_nano = saturation_temperature(model_nano, z, P=1._pr, kind="dew", t0=300._pr)
   ! print*, "pc", pc
   ! print*, "pcc", pcc
   ! print*, "tc", tc
   ! print*, "tcc", tcc


   ! Inicia el encabezado con "T, P"
   names(1) = "N2"
   names(2) = "CH4"
   names(3) = "C2H6"
   names(4) = "C3H8"
   names(5) = "nC4"
   names(6) = "nC5"
   names(7) = "nC6"

   !names=(/"N2", "CH4", "C2H6", "C3H8", "nC4", "nC5", "nC6"/)
   header = "         T                     P"

   ! Añade "lnK1", "lnK2", ..., "lnKn" dinámicamente
   do i = 1, nc
      header = trim(adjustl(header)) // "          lnK" // trim(adjustl(itoa(i))) // "(" // trim(adjustl(names(i))) // ")"

      !write(header(len_trim(header)+1:), "(A, '          lnK', I0)") " ", i,"(",names(i),")"
   end do

   ! Define el formato fijo para los datos
   fmt = "(F25.16, F25.16"  ! Espacio fijo para "T" y "P"
   do i = 1, nc
       fmt = trim(fmt) // ", F25.16"  ! Espacio fijo para cada "lnK"
   end do
   fmt = trim(fmt) // ")"
   ! ------------------------- BUBBLE ENVELOPE ------------------

   sat_point = saturation_temperature(model, z, P=1.5_pr, kind="bubble", t0=115._pr)
   sat_point_nano = saturation_temperature(model_nano, z, P=1.5_pr, kind="bubble", t0=115._pr)
   sat_point_auto_eos = saturation_temperature(auto_eos, z, P=1.5_pr, kind="bubble", t0=115._pr)
   ! Calculate phase envelope
   envelope = pt_envelope_2ph(model, z, sat_point)
   envelope_nano = pt_envelope_2ph(model_nano, z, sat_point_nano, points=500, iterations=1000)
   envelope_auto_eos = pt_envelope_2ph(auto_eos, z, sat_point_auto_eos)

   !write(1, "(*(A,2x))") "kind","T", "P", "beta","x", "y", "Vx", "Vy"
   !write(*, *) envelope%points(1)
   write(1, "(A)") trim(header)
   ! do i=1,size(envelope%points)
   !    write(1,fmt) envelope%points(i)%T ,envelope%points(i)%P,log(envelope%points(i)%y/envelope%points(i)%x)
   !    !write(1,*) log(envelope%points(i)%y/envelope%points(i)%x)
   ! end do
   write(1,*) envelope
   ! write(1,*), envelope%cps
   !write(2, "(*(A,2x))") "kind","T", "P", "beta","x", "y", "Vx", "Vy"
   !write(*, *) envelope_nano%points(28)
   write(3, "(A)") trim(header)
   ! do i=1,size(envelope_nano%points)
   !    !write(1,*) envelope_nano%points(i)%T ,envelope_nano%points(i)%P
   !    write(3,fmt) envelope_nano%points(i)%T ,envelope_nano%points(i)%P, log(envelope_nano%points(i)%y/envelope_nano%points(i)%x)
   ! end do
   write(3,*) envelope_nano
   ! do i=1,size(envelope_auto_eos%points)
   !    !write(3,*) envelope_auto_eos%points(i)%T ,envelope_auto_eos%points(i)%P
   !    write(3,*) envelope_auto_eos%points(i)%T ,envelope_auto_eos%points(i)%P, &
   !    log(envelope_auto_eos%points(i)%y/envelope_auto_eos%points(i)%x)
   ! end do
   !write(3,*) envelope_auto_eos
   !write(2,*) envelope_nano
   !print*, size(envelope_nano%points)

   !! -------------------------- DEW ENVELOPE ---------------------------

   sat_point = saturation_temperature(model, z, P=0.5_pr, kind="dew", t0=200._pr)
   sat_point_nano = saturation_temperature(model_nano, z, P=0.5_pr, kind="dew", t0=200._pr)
   sat_point_auto_eos = saturation_temperature(auto_eos, z, P=0.5_pr, kind="dew", t0=200._pr)

   envelope = pt_envelope_2ph(model, z, sat_point)
   envelope_nano = pt_envelope_2ph(model_nano, z, sat_point_nano, points=500, iterations=500)
   envelope_auto_eos = pt_envelope_2ph(auto_eos, z, sat_point)

   ! write(3, "(*(A,2x))") "kind","T", "P", "beta","x", "y", "Vx", "Vy"
   ! write(*, *) envelope%points(1)
   ! write(3,*) envelope
   
   ! write(4, "(*(A,2x))") "kind","T", "P", "beta","x", "y", "Vx", "Vy"
   ! write(*, *) envelope_nano%points(1)
   ! write(4,*) envelope_nano
   write(2, "(A)") trim(header)

   do i=1,size(envelope%points)
      !write(2,*) envelope%points(i)%T ,envelope%points(i)%P
      write(2,fmt) envelope%points(i)%T ,envelope%points(i)%P, log(envelope%points(i)%y/envelope%points(i)%x)
   end do
   
   write(4, "(A)") trim(header)
   do i=1,size(envelope_nano%points)
      !write(2,*) envelope_nano%points(i)%T ,envelope_nano%points(i)%P
      write(4,fmt) envelope_nano%points(i)%T ,envelope_nano%points(i)%P, log(envelope_nano%points(i)%y/envelope_nano%points(i)%x)
   end do
   ! do i=1,size(envelope_auto_eos%points)
   !    !write(4,*) envelope_auto_eos%points(i)%T ,envelope_auto_eos%points(i)%P
   !    write(4,*) envelope_auto_eos%points(i)%T ,envelope_auto_eos%points(i)%P, &
   !    log(envelope_auto_eos%points(i)%y/envelope_auto_eos%points(i)%x)
   ! end do
   ! !write(4,*) envelope_auto_eos
   ! write(6,*) "asd"

   contains

   ! Función para convertir números enteros a cadenas (para el índice i)
   character(len=10) function itoa(i)
     integer, intent(in) :: i
     write(itoa, '(I0)') i
   end function itoa

end program sandoval_phase_test