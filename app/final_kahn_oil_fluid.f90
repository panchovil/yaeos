program kahn_oil_fluid
   use yaeos, only: pr, ArModel, PengRobinson78, PengRobinson78Nano,&
               saturation_pressure, saturation_temperature, k_wilson ,&
               EquilibriumState, PTEnvel2, pt_envelope_2ph,&
               NanoEquilibriumState, NanoPTEnvel2, nano_pt_envelope_2ph , find_hpl

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
   real(pr) :: z1(nc), z2(nc), z3(nc), z4(nc)

   character(len=500) :: header, header_cap, names(nc)                          ! Format of envelopes
   character(len=1000) :: fmt, fmt_cap                                          ! Format of envelopes

   conv_bar_kPa = 100._pr
   ! ===========================================================================
   ! Format
   ! ---------------------------------------------------------------------------
   names(1) = "CO2"
   names(2) = "C1"
   names(3) = "C2-3"
   names(4) = "C4-6"
   names(5) = "C7-14"
   names(6) = "C15-25"
   names(7) = "C26+"
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

   !! names = 'CO2', 'C1', 'C2-3', 'C4-6', 'C7-14', 'C15-25', 'C26+'
   
   
   !! composition vector                        
   z = (/0.0169_pr, 0.1752_pr, 0.2244_pr, 0.1673_pr, 0.2422_pr, 0.1216_pr, 0.0524_pr/)



   z(1) = z(1)*100

   z = z/sum(z)
   z = z/sum(z)
   print*, sum(z), z
   !! Critical Pressure
   pc = (/73.76460200931301_pr, 46.001551764484894_pr, 44.68907955567195_pr, &
   34.17938032942354_pr, 21.868639497810967_pr, 16.036936568375182_pr, 15.212254648539314_pr/)
   !! Critical Temperature
   tc = (/304.2_pr, 174.44444444444446_pr, 347.26311111111113_pr, 459.73972222222227_pr, &
   595.1350555555555_pr, 729.9807777777778_pr, 910.1830555555556_pr/)
   !! Acentric Factor
   w = (/0.225_pr, 0.008_pr, 0.1331_pr, 0.2358_pr, 0.5977_pr, 0.9118_pr, 1.2444_pr/)
   !! lij and kij matrix
   lij = 0.0_pr
   kij = 0.0_pr
   kij(1, :) = (/ 0.0_pr, 0.085_pr, 0.085_pr, 0.085_pr, 0.104_pr, 0.104_pr, 0.104_pr /)
   kij(:, 1) = (/ 0.0_pr, 0.085_pr, 0.085_pr, 0.085_pr, 0.104_pr, 0.104_pr, 0.104_pr /)
   
      
   !! Capillary pressure variables
   rp = 10._pr !! [nm]
   rp = rp/1E9 !! [m]
   ang_cont = 1.0472 !! angulo cualquiera de 60º en radianes
   LJ_par = 0.244_pr*((Tc/Pc)**(1.0/3.0))  !! [nm]
   Parachor = 40.1684*(0.151-0.0464*w)*(tc**(13._pr/12._pr))/(pc**(5._pr/6._pr)) !! [cm^3/mol * (mN/m)^1/4]

   
   !! Model definition
   model_bulk = PengRobinson78(tc, pc, w, kij, lij)
   model_conf = PengRobinson78Nano(LJ_par, rp*1E9, tc, pc, w, kij, lij)  ! pore radius in nanometers

   

   ! ! ===========================================================================
   ! ! BUBBLE ENVELOPES
   ! ! ---------------------------------------------------------------------------
   ! !! initialize the phase envelope
   sat_point_bulk = saturation_temperature(model_bulk, z, P=1.5_pr, kind="bubble", t0=200._pr)
   ! sat_point_conf = saturation_temperature(model_conf, z, P=31.5_pr, kind="bubble", t0=300._pr)
   ! ! init_point_bulk_cap = capillary_init_point(model_bulk, z, Parachor, ang_cont, rp, "bubble",&
   ! ! T0_in=300._pr, P_in=31.5_pr)
   ! init_point_conf_cap = capillary_init_point(model_conf, z, Parachor, ang_cont, rp, "bubble",&
   ! T0_in=300._pr, P_in=31.5_pr)

   ! !! Calculate the phase envelope
   envelope_bulk = pt_envelope_2ph(model_bulk, z, sat_point_bulk)
   ! envelope_conf = pt_envelope_2ph(model_conf, z, sat_point_conf, 2000, 1000)
   ! ! envelope_bulk_cap = nano_pt_envelope_2ph(model_bulk, z, rp, ang_cont, Parachor, init_point_bulk_cap, 2000, 1000)
   ! envelope_conf_cap = nano_pt_envelope_2ph(model_conf, z, rp, ang_cont, Parachor, init_point_conf_cap, 1270, 1000)

   ! !! writing envelope
   ! !! --------------------------------------------------------------------------
   ! !! Bulk envelope
   write(1, "(A)") trim(header)
   do i=1,size(envelope_bulk%points)
      write(1,fmt) envelope_bulk%points(i)%T, envelope_bulk%points(i)%P*conv_bar_kPa, &
      log(envelope_bulk%points(i)%y/envelope_bulk%points(i)%x)
   end do
   write(1,*) envelope_bulk%cps%T, envelope_bulk%cps%P*conv_bar_kPa
   
   ! !! --------------------------------------------------------------------------
   ! !! Conf envelope
   ! write(1, "(A)") trim(header)
   ! do i=1,size(envelope_conf%points)
   !    write(1,fmt) envelope_conf%points(i)%T, envelope_conf%points(i)%P*conv_bar_kPa, &
   !    log(envelope_conf%points(i)%y/envelope_conf%points(i)%x)
   ! end do
   ! write(1,*) envelope_conf%cps%T, envelope_conf%cps%P*conv_bar_kPa

   ! !! --------------------------------------------------------------------------
   ! !! Cap envelopre
   ! ! write(3, "(A)") trim(header_cap)
   ! ! do i=1,size(envelope_bulk_cap%points)
   ! !    write(3,fmt_cap) envelope_bulk_cap%points(i)%T, envelope_bulk_cap%points(i)%Py*conv_bar_kPa, &
   ! !    envelope_bulk_cap%points(i)%Px*conv_bar_kPa, envelope_bulk_cap%points(i)%Pcap*conv_bar_kPa,&
   ! !    log(envelope_bulk_cap%points(i)%y/envelope_bulk_cap%points(i)%x)
   ! ! end do
   ! ! write(3,*) envelope_bulk_cap%cps%T, envelope_bulk_cap%cps%P*conv_bar_kPa

   ! !! --------------------------------------------------------------------------
   ! !! Conf cap envelope
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
   sat_point_bulk = saturation_temperature(model_bulk, z, P=0.5_pr, kind="dew", t0=600._pr)
   ! sat_point_conf = saturation_temperature(model_conf, z, P=2.5_pr, kind="dew", t0=620._pr)
   ! init_point_bulk_cap = capillary_init_point(model_bulk, z, Parachor, ang_cont, rp, "dew",&
   ! T0_in=620._pr, P_in=2.5_pr)
   ! init_point_conf_cap = capillary_init_point(model_conf, z, Parachor, ang_cont, rp, "dew",&
   ! T0_in=620._pr, P_in=2.5_pr)

   !! Calculate the phase envelope
   envelope_bulk = pt_envelope_2ph(model_bulk, z, sat_point_bulk)
   ! envelope_conf = pt_envelope_2ph(model_conf, z, sat_point_conf, 2000, 1000)
   ! envelope_bulk_cap = nano_pt_envelope_2ph(model_bulk, z, rp, ang_cont, Parachor, init_point_bulk_cap, 2000, 1000)
   ! envelope_conf_cap = nano_pt_envelope_2ph(model_conf, z, rp, ang_cont, Parachor, init_point_conf_cap, 20000, 1000)

   !! writing envelope
   !! --------------------------------------------------------------------------
   !! Bulk envelope
   write(2, "(A)") trim(header)
   do i=1,size(envelope_bulk%points)
      write(2,fmt) envelope_bulk%points(i)%T, envelope_bulk%points(i)%P*conv_bar_kPa, &
      log(envelope_bulk%points(i)%y/envelope_bulk%points(i)%x)
   end do
   write(2,*) envelope_bulk%cps%T, envelope_bulk%cps%P*conv_bar_kPa
   
   !! --------------------------------------------------------------------------
   !! Conf envelope
   ! write(2, "(A)") trim(header)
   ! do i=1,size(envelope_conf%points)
   !    write(2,fmt) envelope_conf%points(i)%T, envelope_conf%points(i)%P*conv_bar_kPa, &
   !    log(envelope_conf%points(i)%y/envelope_conf%points(i)%x)
   ! end do
   ! write(2,*) envelope_conf%cps%T, envelope_conf%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   ! Cap envelopre
   ! write(4, "(A)") trim(header_cap)
   ! do i=1,size(envelope_bulk_cap%points)
   !    write(4,fmt_cap) envelope_bulk_cap%points(i)%T, envelope_bulk_cap%points(i)%Py*conv_bar_kPa, &
   !    envelope_bulk_cap%points(i)%Px*conv_bar_kPa, envelope_bulk_cap%points(i)%Pcap*conv_bar_kPa,&
   !    log(envelope_bulk_cap%points(i)%y/envelope_bulk_cap%points(i)%x)
   ! end do
   ! write(4,*) envelope_bulk_cap%cps%T, envelope_bulk_cap%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Conf cap envelope
   ! write(4, "(A)") trim(header_cap)
   ! do i=1,size(envelope_conf_cap%points)
   !    write(4,fmt_cap) envelope_conf_cap%points(i)%T, envelope_conf_cap%points(i)%Py*conv_bar_kPa, &
   !    envelope_conf_cap%points(i)%Px*conv_bar_kPa, envelope_conf_cap%points(i)%Pcap*conv_bar_kPa,&
   !    log(envelope_conf_cap%points(i)%y/envelope_conf_cap%points(i)%x)
   ! end do
   ! write(4,*) envelope_conf_cap%cps%T, envelope_conf_cap%cps%P*conv_bar_kPa

   ! ===========================================================================
   ! HPL ENVELOPES
   ! ---------------------------------------------------------------------------
   !! initialize the phase envelope
   ! init_point_bulk_cap = capillary_init_point(model_bulk, z, Parachor, ang_cont, rp, "liquid-liquid",&
   ! T0_in=700._pr, P0_in=600._pr)
   ! init_point_conf_cap = capillary_init_point(model_conf, z, Parachor, ang_cont, rp, "liquid-liquid",&
   ! T0_in=700._pr, P0_in=600._pr)


   !! Calculate the phase envelope
   envelope_bulk = find_hpl(model_bulk, z, 700._pr, 1000._pr)
   ! envelope_conf = find_hpl(model_conf, z, 700._pr, 1000._pr)
   ! envelope_bulk_cap = nano_pt_envelope_2ph(model_bulk, z, rp, ang_cont, Parachor, init_point_bulk_cap, &
   ! points=20000, iterations=1000, specified_variable_0=nc+3, delta_0=-5.0_pr)
   ! envelope_conf_cap = nano_pt_envelope_2ph(model_conf, z, rp, ang_cont, Parachor, init_point_conf_cap, &
   ! points=20000, iterations=1000, specified_variable_0=nc+3, delta_0=-5.0_pr)

   !! writing envelope
   !! --------------------------------------------------------------------------
   !! Bulk envelope
   write(3, "(A)") trim(header)
   do i=1,size(envelope_bulk%points)
      write(3,fmt) envelope_bulk%points(i)%T, envelope_bulk%points(i)%P*conv_bar_kPa, &
      log(envelope_bulk%points(i)%y/envelope_bulk%points(i)%x)
   end do
   write(3,*) envelope_bulk%cps%T, envelope_bulk%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Conf envelope
   ! write(1, "(A)") trim(header)
   ! do i=1,size(envelope_conf%points)
   !    write(1,fmt) envelope_conf%points(i)%T, envelope_conf%points(i)%P*conv_bar_kPa, &
   !    log(envelope_conf%points(i)%y/envelope_conf%points(i)%x)
   ! end do
   ! write(1,*) envelope_conf%cps%T, envelope_conf%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Cap envelopre
   ! write(3, "(A)") trim(header_cap)
   ! do i=1,size(envelope_bulk_cap%points)
   !    write(3,fmt_cap) envelope_bulk_cap%points(i)%T, envelope_bulk_cap%points(i)%Py*conv_bar_kPa, &
   !    envelope_bulk_cap%points(i)%Px*conv_bar_kPa, envelope_bulk_cap%points(i)%Pcap*conv_bar_kPa,&
   !    log(envelope_bulk_cap%points(i)%y/envelope_bulk_cap%points(i)%x)
   ! end do
   ! write(3,*) envelope_bulk_cap%cps%T, envelope_bulk_cap%cps%P*conv_bar_kPa

   !! --------------------------------------------------------------------------
   !! Conf cap envelope
   ! write(3, "(A)") trim(header_cap)
   ! do i=1,size(envelope_conf_cap%points)
   !    write(3,fmt_cap) envelope_conf_cap%points(i)%T, envelope_conf_cap%points(i)%Py*conv_bar_kPa, &
   !    envelope_conf_cap%points(i)%Px*conv_bar_kPa, envelope_conf_cap%points(i)%Pcap*conv_bar_kPa,&
   !    log(envelope_conf_cap%points(i)%y/envelope_conf_cap%points(i)%x)
   ! end do
   ! write(3,*) envelope_conf_cap%cps%T, envelope_conf_cap%cps%P*conv_bar_kPa


   print*, sum(z), z


   contains
      ! Función para convertir números enteros a cadenas (para el índice i)
      character(len=10) function itoa(i)
         integer, intent(in) :: i
         write(itoa, '(I0)') i
      end function itoa
end program kahn_oil_fluid