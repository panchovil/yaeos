program paper_fluid
   use yaeos, only: pr, ArModel, PengRobinson78, PengRobinson78Nano,&
               saturation_pressure, saturation_temperature, k_wilson ,&
               EquilibriumState, PTEnvel2, pt_envelope_2ph,&
               NanoEquilibriumState, NanoPTEnvel2, nano_pt_envelope_2ph , find_hpl

   use Capillary_Initializer, only: capillary_init_point               
   implicit none

   ! ===========================================================================
   ! Variables definition
   ! ---------------------------------------------------------------------------
   integer, parameter :: nc=5  
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
   names(1) = "C1"
   names(2) = "C2-C4"
   names(3) = "C5-C7"
   names(4) = "C8-C9"
   names(5) = "nC10+"
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

   !! names="C1" "C2-C4" "C5-C7" "C8-C9" "C10+"
   
   
   !! composition vector                        
   z = (/0.2506_pr, 0.22_pr, 0.20_pr, 0.13_pr, 0.1994_pr/)
   z = z/sum(z)
   print*, sum(z), z

   !! Critical Pressure
   pc = (/45.40_pr, 42.54_pr, 33.76_pr, 30.91_pr, 21.58_pr/)

   !! Critical Temperature
   tc = (/190.60_pr, 363.30_pr, 511.56_pr, 579.34_pr, 788.74_pr/)

   !! Acentric Factor
   w = (/0.008_pr, 0.1432_pr, 0.2474_pr, 0.2861_pr, 0.6869_pr/)

   !! lij and kij matrix
   lij = 0.0_pr
   kij = 0.0_pr
   kij = reshape([ &
   0.0_pr, 0.0078_pr, 0.0242_pr, 0.0324_pr, 0.0779_pr, &  ! Fila 1
   0.0078_pr, 0.0_pr, 0.0046_pr, 0.0087_pr, 0.0384_pr, &  ! Fila 2
   0.0242_pr, 0.0046_pr, 0.0_pr, 0.0006_pr, 0.0169_pr, &  ! Fila 3
   0.0324_pr, 0.0087_pr, 0.0006_pr, 0.0_pr, 0.0111_pr, &  ! Fila 4
   0.0779_pr, 0.0384_pr, 0.0169_pr, 0.0111_pr, 0.0_pr  &  ! Fila 5
   ], shape=[nc, nc])
   
      
   !! Capillary pressure variables
   rp = 5._pr !! [nm]
   rp = rp/1E9 !! [m]
   ang_cont = 1.0472 !! angulo cualquiera de 60º en radianes
   LJ_par = 0.244_pr*((Tc/Pc)**(1.0/3.0))  !! [nm]
   Parachor = (/77.0_pr, 145.2_pr, 250.0_pr, 306.0_pr, 686.3_pr/) !! [cm^3/mol * (mN/m)^1/4]

   !Parachor = 40.1684*(0.151-0.0464*w)*(tc**(13._pr/12._pr))/(pc**(5._pr/6._pr)) !! [cm^3/mol * (mN/m)^1/4]

   
   !! Model definition
   model_bulk = PengRobinson78(tc, pc, w, kij, lij)
   model_conf = PengRobinson78Nano(LJ_par, rp*1E9, tc, pc, w, kij, lij)  ! pore radius in nanometers

   

   ! ===========================================================================
   ! BUBBLE ENVELOPES
   ! ---------------------------------------------------------------------------
   !! initialize the phase envelope
   ! sat_point_bulk = saturation_temperature(model_bulk, z, P=90._pr, kind="bubble", t0=310._pr)
   sat_point_conf = saturation_temperature(model_conf, z, P=90._pr, kind="bubble", t0=310._pr)
   init_point_bulk_cap = capillary_init_point(model_bulk, z, Parachor, ang_cont, rp, "bubble",&
   T0_in=180._pr, P_in=30._pr)
   ! init_point_conf_cap = capillary_init_point(model_conf, z, Parachor, ang_cont, rp, "bubble",&
   ! T0_in=310._pr, P_in=90._pr)

   !! Calculate the phase envelope
   ! envelope_bulk = pt_envelope_2ph(model_bulk, z, sat_point_bulk)
   envelope_conf = pt_envelope_2ph(model_conf, z, sat_point_conf, 500, 1000)!500, 500, 500, 500, 500, 500
   envelope_bulk_cap = nano_pt_envelope_2ph(model_bulk, z, rp, ang_cont, Parachor, init_point_bulk_cap, 1700, 1000) !500, 600, 700 900 P=70/ T=300/ 1300 P=30/ T=180/ 1700
   ! envelope_conf_cap = nano_pt_envelope_2ph(model_conf, z, rp, ang_cont, Parachor, init_point_conf_cap, 1450, 1000)!500, 600, 700, 1000, 1450, 1950

   !! writing envelope
   !! --------------------------------------------------------------------------
   !! Bulk envelope
   ! write(1, "(A)") trim(header)
   ! do i=1,size(envelope_bulk%points)
   !    write(1,fmt) envelope_bulk%points(i)%T, envelope_bulk%points(i)%P*conv_bar_kPa, &
   !    log(envelope_bulk%points(i)%y/envelope_bulk%points(i)%x)
   ! end do
   ! write(1,*) envelope_bulk%cps%T, envelope_bulk%cps%P*conv_bar_kPa
   
   !! --------------------------------------------------------------------------
   !! Conf envelope
   write(1, "(A)") trim(header)
   do i=1,size(envelope_conf%points)
      write(1,fmt) envelope_conf%points(i)%T, envelope_conf%points(i)%P*conv_bar_kPa, &
      log(envelope_conf%points(i)%y/envelope_conf%points(i)%x)
   end do
   write(1,*) envelope_conf%cps%T, envelope_conf%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Cap envelopre
   write(3, "(A)") trim(header_cap)
   do i=1,size(envelope_bulk_cap%points)
      write(3,fmt_cap) envelope_bulk_cap%points(i)%T, envelope_bulk_cap%points(i)%Py*conv_bar_kPa, &
      envelope_bulk_cap%points(i)%Px*conv_bar_kPa, envelope_bulk_cap%points(i)%Pcap*conv_bar_kPa,&
      log(envelope_bulk_cap%points(i)%y/envelope_bulk_cap%points(i)%x)
   end do
   write(3,*) envelope_bulk_cap%cps%T, envelope_bulk_cap%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Conf cap envelope
   ! write(3, "(A)") trim(header_cap)
   ! do i=1,size(envelope_conf_cap%points)
   !    write(3,fmt_cap) envelope_conf_cap%points(i)%T, envelope_conf_cap%points(i)%Py*conv_bar_kPa, &
   !    envelope_conf_cap%points(i)%Px*conv_bar_kPa, envelope_conf_cap%points(i)%Pcap*conv_bar_kPa,&
   !    log(envelope_conf_cap%points(i)%y/envelope_conf_cap%points(i)%x)
   ! end do
   ! write(3,*) envelope_conf_cap%cps%T, envelope_conf_cap%cps%P*conv_bar_kPa


   ! ===========================================================================
   ! DEW ENVELOPES
   ! ---------------------------------------------------------------------------
   !! initialize the phase envelope
   ! sat_point_bulk = saturation_temperature(model_bulk, z, P=0.5_pr, kind="dew", t0=500._pr)
   sat_point_conf = saturation_temperature(model_conf, z, P=0.5_pr, kind="dew", t0=500._pr)
   init_point_bulk_cap = capillary_init_point(model_bulk, z, Parachor, ang_cont, rp, "dew",&
   T0_in=400._pr, P_in=2.5_pr)
   ! init_point_conf_cap = capillary_init_point(model_conf, z, Parachor, ang_cont, rp, "dew",&
   ! T0_in=500._pr, P_in=1.5_pr)

   !! Calculate the phase envelope
   ! envelope_bulk = pt_envelope_2ph(model_bulk, z, sat_point_bulk)
   envelope_conf = pt_envelope_2ph(model_conf, z, sat_point_conf, points= 260, iterations=1000)!none, 450, 450, 400, 330, 260
   envelope_bulk_cap = nano_pt_envelope_2ph(model_bulk, z, rp, ang_cont, Parachor, init_point_bulk_cap, 1700, 1000) ! P=2,5/ 800, 850, 850, 850, 1400, 1700
   ! envelope_conf_cap = nano_pt_envelope_2ph(model_conf, z, rp, ang_cont, Parachor, init_point_conf_cap, 1800, 1000)!800, 850/P=0.5, 1000, 1200, 1800, 2500

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
   ! Cap envelopre
   write(4, "(A)") trim(header_cap)
   do i=1,size(envelope_bulk_cap%points)
      write(4,fmt_cap) envelope_bulk_cap%points(i)%T, envelope_bulk_cap%points(i)%Py*conv_bar_kPa, &
      envelope_bulk_cap%points(i)%Px*conv_bar_kPa, envelope_bulk_cap%points(i)%Pcap*conv_bar_kPa,&
      log(envelope_bulk_cap%points(i)%y/envelope_bulk_cap%points(i)%x)
   end do
   write(4,*) envelope_bulk_cap%cps%T, envelope_bulk_cap%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Conf cap envelope
   ! write(4, "(A)") trim(header_cap)
   ! do i=1,size(envelope_conf_cap%points)
   !    write(4,fmt_cap) envelope_conf_cap%points(i)%T, envelope_conf_cap%points(i)%Py*conv_bar_kPa, &
   !    envelope_conf_cap%points(i)%Px*conv_bar_kPa, envelope_conf_cap%points(i)%Pcap*conv_bar_kPa,&
   !    log(envelope_conf_cap%points(i)%y/envelope_conf_cap%points(i)%x)
   ! end do
   ! write(4,*) envelope_conf_cap%cps%T, envelope_conf_cap%cps%P*conv_bar_kPa


   print*, sum(z), z


   contains
      ! Función para convertir números enteros a cadenas (para el índice i)
      character(len=10) function itoa(i)
         integer, intent(in) :: i
         write(itoa, '(I0)') i
      end function itoa
end program paper_fluid