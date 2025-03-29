program sandoval_fluid
   use yaeos, only: pr, ArModel, PengRobinson78, PengRobinson78Nano,&
               saturation_pressure, saturation_temperature, k_wilson ,&
               EquilibriumState, PTEnvel2, pt_envelope_2ph,&
               NanoEquilibriumState, NanoPTEnvel2, nano_pt_envelope_2ph 

   use Capillary_Initializer, only: capillary_init_point               
   implicit none

   ! ===========================================================================
   ! Variables definition
   ! ---------------------------------------------------------------------------
   integer, parameter :: nc=7  
   integer :: i
   class(ArModel), allocatable :: model_bulk, model_conf                        ! Thermodynamic model to be used
   type(EquilibriumState) :: sat_point_bulk, sat_point_conf                     ! Init without Capillary
   type(NanoEquilibriumState) :: init_point_bulk_cap, init_point_conf_cap       ! Init with Capillary
   type(PTEnvel2) :: envelope_bulk, envelope_conf                               ! PT Phase envelope without Capillary
   type(NanoPTEnvel2) :: envelope_bulk_cap ,envelope_conf_cap                   ! PT Phase envelope with Capillary

   real(pr) :: tc(nc), pc(nc), w(nc)                                            ! Component's critical constants
   real(pr) :: z(nc), kij(nc,nc), lij(nc,nc)                                    ! Termodynamic variables
   real(pr) :: rp, LJ_par(nc), ang_cont, Parachor(nc)  
   real(pr) :: conv_bar_kPa                         ! Nano parameters

   character(len=500) :: header, header_cap, names(nc)                          ! Format of envelopes
   character(len=1000) :: fmt, fmt_cap                                          ! Format of envelopes

   conv_bar_kPa = 100._pr
   ! ===========================================================================
   ! Format
   ! ---------------------------------------------------------------------------
   names(1) = "N2"
   names(2) = "CH4"
   names(3) = "C2H6"
   names(4) = "C3H8"
   names(5) = "nC4"
   names(6) = "nC5"
   names(7) = "nC6"
   header = "         T                     P"
   header_cap = "         T                     Py                     Px                     Pcap"

   ! Add lnK1", "lnK2", ..., "lnKn" dynamically
   do i = 1, nc
      header = trim(adjustl(header)) // "          lnK" // trim(adjustl(itoa(i))) // "(" // trim(adjustl(names(i))) // ")"
      header_cap = trim(adjustl(header_cap)) // "          lnK" // trim(adjustl(itoa(i))) // "(" // trim(adjustl(names(i))) // ")"
   end do

   ! Define fix format for the data 
   fmt = "(F25.16, F25.16"                      ! Fix space for "T" y "P"
   fmt_cap = "(F25.16, F25.16, F25.16, F25.16"  ! Fix space for "T" , "Py", "Px" y "Pcap"
   do i = 1, nc
       fmt = trim(fmt) // ", F25.16"            ! Fix space for each "lnK"
       fmt_cap = trim(fmt_cap) // ", F25.16"    ! Fix space for each "lnK"
   end do
   fmt = trim(fmt) // ")"
   fmt_cap = trim(fmt_cap) // ")"

   ! ===========================================================================
   ! Values loading
   ! ---------------------------------------------------------------------------

   !! names="N2" "CH4" "C2H6" "C3H8" "nC4" "nC5" "nC6"
   
   !! composition vector                        
   z = (/0.0014_pr, 0.943_pr, 0.027_pr, 0.0074_pr, 0.0049_pr, 0.0027_pr, 0.001_pr/)
   !! Critical Pressure
   pc = (/34.0_pr, 45.99_pr, 48.72_pr, 42.48_pr, 37.96_pr, 33.70_pr, 30.25_pr/)
   !! Critical Temperature
   tc = (/126.20_pr, 190.56_pr, 305.32_pr, 369.83_pr, 425.12_pr, 469.70_pr, 507.60_pr/)
   !! Acentric Factor
   w = (/0.0377_pr, 0.0115_pr, 0.0995_pr, 0.1523_pr, 0.2002_pr, 0.2515_pr, 0.3013_pr/)
   !! lij and kij matrix
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
      
   !! Capillary pressure variables
   rp = 10 !! [nm]
   rp = rp/1E9 !! [m]
   ang_cont = 1.0472 !! angulo cualquiera de 60º en radianes
   LJ_par = (/0.364_pr, 0.375_pr, 0.443_pr, 0.472_pr, 0.506_pr, 0.529_pr, 0.550_pr/) !! [nm]
   Parachor = (/61.12, 74.05, 112.91, 154.03, 193.90, 236.00, 276.71/)!! [cm^3/mol * (mN/m)^1/4]
   
   !! Model definition
   model_bulk = PengRobinson78(tc, pc, w, kij, lij)
   model_conf = PengRobinson78Nano(LJ_par, rp*1E9, tc, pc, w, kij, lij)  ! pore radius in nanometers

   ! ===========================================================================
   ! BUBBLE ENVELOPES
   ! ---------------------------------------------------------------------------
   !! initialize the phase envelope
   sat_point_bulk = saturation_temperature(model_bulk, z, P=1.5_pr, kind="bubble", t0=115._pr)
   ! sat_point_conf = saturation_temperature(model_conf, z, P=1.5_pr, kind="bubble", t0=115._pr)
   ! init_point_bulk_cap = capillary_init_point(model_bulk, z, Parachor, ang_cont, rp, "bubble",&
   ! T0_in=115._pr, P_in=1.5_pr)
   init_point_conf_cap = capillary_init_point(model_conf, z, Parachor, ang_cont, rp, "bubble",&
   T0_in=115._pr, P_in=1.5_pr)

   !! Calculate the phase envelope
   envelope_bulk = pt_envelope_2ph(model_bulk, z, sat_point_bulk)
   ! envelope_conf = pt_envelope_2ph(model_conf, z, sat_point_conf, iterations=1000 )
   ! envelope_bulk_cap = nano_pt_envelope_2ph(model_bulk, z, rp, ang_cont, Parachor, init_point_bulk_cap, 5000)
   envelope_conf_cap = nano_pt_envelope_2ph(model_conf, z, rp, ang_cont, Parachor, init_point_conf_cap, 1800)

   !! writing envelope
   !! --------------------------------------------------------------------------
   !! Bulk envelope
   write(1, "(A)") trim(header)
   do i=1,size(envelope_bulk%points)
      write(1,fmt) envelope_bulk%points(i)%T, envelope_bulk%points(i)%P*conv_bar_kPa, &
      log(envelope_bulk%points(i)%y/envelope_bulk%points(i)%x)
   end do
   write(1,*) envelope_bulk%cps%T, envelope_bulk%cps%P*conv_bar_kPa
   
   !! --------------------------------------------------------------------------
   ! Conf envelope
   ! write(1, "(A)") trim(header)
   ! do i=1,size(envelope_conf%points)
   !    write(1,fmt) envelope_conf%points(i)%T, envelope_conf%points(i)%P*conv_bar_kPa, &
   !    log(envelope_conf%points(i)%y/envelope_conf%points(i)%x)
   ! end do
   ! write(1,*) envelope_conf%cps%T, envelope_conf%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Cap envelopre
   ! write(1, "(A)") trim(header_cap)
   ! do i=1,size(envelope_bulk_cap%points)
   !    write(1,fmt_cap) envelope_bulk_cap%points(i)%T, envelope_bulk_cap%points(i)%Py*conv_bar_kPa, &
   !    envelope_bulk_cap%points(i)%Px*conv_bar_kPa, envelope_bulk_cap%points(i)%Pcap*conv_bar_kPa,&
   !    log(envelope_bulk_cap%points(i)%y/envelope_bulk_cap%points(i)%x)
   ! end do
   ! write(1,*) envelope_bulk_cap%cps%T, envelope_bulk_cap%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Conf cap envelope
   write(3, "(A)") trim(header_cap)
   do i=1,size(envelope_conf_cap%points)
      write(3,fmt_cap) envelope_conf_cap%points(i)%T, envelope_conf_cap%points(i)%Py*conv_bar_kPa, &
      envelope_conf_cap%points(i)%Px*conv_bar_kPa, envelope_conf_cap%points(i)%Pcap*conv_bar_kPa,&
      log(envelope_conf_cap%points(i)%y/envelope_conf_cap%points(i)%x)
   end do
   write(3,*) envelope_conf_cap%cps%T, envelope_conf_cap%cps%P*conv_bar_kPa

   ! ===========================================================================
   ! DEW ENVELOPES
   ! ---------------------------------------------------------------------------
   !! initialize the phase envelope
   ! sat_point_bulk = saturation_temperature(model_bulk, z, P=0.5_pr, kind="dew", t0=200._pr)
   sat_point_conf = saturation_temperature(model_conf, z, P=0.5_pr, kind="dew", t0=200._pr)
   ! init_point_bulk_cap = capillary_init_point(model_bulk, z, Parachor, ang_cont, rp, "dew",&
   ! T0_in=200._pr, P_in=0.5_pr)
   init_point_conf_cap = capillary_init_point(model_conf, z, Parachor, ang_cont, rp, "dew",&
   T0_in=200._pr, P_in=0.5_pr)

   !! Calculate the phase envelope
   ! envelope_bulk = pt_envelope_2ph(model_bulk, z, sat_point_bulk)
   envelope_conf = pt_envelope_2ph(model_conf, z, sat_point_conf, 5)
   ! envelope_bulk_cap = nano_pt_envelope_2ph(model_bulk, z, rp, ang_cont, Parachor, init_point_bulk_cap, 5000)
   envelope_conf_cap = nano_pt_envelope_2ph(model_conf, z, rp, ang_cont, Parachor, init_point_conf_cap, 5)

   !! writing envelope
   !! --------------------------------------------------------------------------
   !! Bulk envelope
   ! write(2, "(A)") trim(header)
   ! do i=1,size(envelope_bulk%points)
   !    write(2,fmt) envelope_bulk%points(i)%T, envelope_bulk%points(i)%P*conv_bar_kPa, &
   !    log(envelope_bulk%points(i)%y/envelope_bulk%points(i)%x)
   ! end do
   ! write(2,*) envelope_bulk%cps%T, envelope_bulk%cps%P*conv_bar_kPa
   
   !! --------------------------------------------------------------------------
   !! Conf envelope
   write(2, "(A)") trim(header)
   do i=1,size(envelope_conf%points)
      write(2,fmt) envelope_conf%points(i)%T, envelope_conf%points(i)%P*conv_bar_kPa, &
      log(envelope_conf%points(i)%y/envelope_conf%points(i)%x)
   end do
   write(2,*) envelope_conf%cps%T, envelope_conf%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Cap envelopre
   ! write(2, "(A)") trim(header_cap)
   ! do i=1,size(envelope_bulk_cap%points)
   !    write(2,fmt_cap) envelope_bulk_cap%points(i)%T, envelope_bulk_cap%points(i)%Py*conv_bar_kPa, &
   !    envelope_bulk_cap%points(i)%Px*conv_bar_kPa, envelope_bulk_cap%points(i)%Pcap*conv_bar_kPa,&
   !    log(envelope_bulk_cap%points(i)%y/envelope_bulk_cap%points(i)%x)
   ! end do
   ! write(2,*) envelope_bulk_cap%cps%T, envelope_bulk_cap%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Conf cap envelope
   write(4, "(A)") trim(header_cap)
   do i=1,size(envelope_conf_cap%points)
      write(4,fmt_cap) envelope_conf_cap%points(i)%T, envelope_conf_cap%points(i)%Py*conv_bar_kPa, &
      envelope_conf_cap%points(i)%Px*conv_bar_kPa, envelope_conf_cap%points(i)%Pcap*conv_bar_kPa,&
      log(envelope_conf_cap%points(i)%y/envelope_conf_cap%points(i)%x)
   end do
   write(4,*) envelope_conf_cap%cps%T, envelope_conf_cap%cps%P*conv_bar_kPa

   contains
      ! Función para convertir números enteros a cadenas (para el índice i)
      character(len=10) function itoa(i)
         integer, intent(in) :: i
         write(itoa, '(I0)') i
      end function itoa
end program sandoval_fluid