program feos_fluid
   use yaeos, only: pr, ArModel, PengRobinson78, PengRobinson78Nano,&
               saturation_pressure, saturation_temperature, k_wilson ,&
               EquilibriumState, PTEnvel2, pt_envelope_2ph,&
               NanoEquilibriumState, NanoPTEnvel2, nano_pt_envelope_2ph, find_hpl 

   use Capillary_Initializer, only: capillary_init_point  
   implicit none

   integer, parameter :: nc=8  
   integer :: i
   class(ArModel), allocatable :: model_bulk, model_conf                        ! Thermodynamic model to be used
   type(EquilibriumState) :: sat_point_bulk, sat_point_conf                     ! Init without Capillary
   type(NanoEquilibriumState) :: init_point_bulk_cap, init_point_conf_cap       ! Init with Capillary
   type(PTEnvel2) :: envelope_bulk, envelope_conf                               ! PT Phase envelope without Capillary
   type(NanoPTEnvel2) :: envelope_bulk_cap ,envelope_conf_cap                   ! PT Phase envelope with Capillary

   real(pr) :: tc(nc), pc(nc), w(nc)                                            ! Component's critical constants
   real(pr) :: z(nc), kij(nc,nc), lij(nc,nc)                                    ! Termodynamic variables
   real(pr) :: rp, LJ_par(nc), ang_cont, Parachor(nc)                           ! Nano parameters
   real(pr) :: conv_bar_kPa                                                 

   character(len=500) :: header, header_cap, names(nc)                          ! Format of envelopes
   character(len=1000) :: fmt, fmt_cap                                          ! Format of envelopes

   conv_bar_kPa = 100._pr
   ! ===========================================================================
   ! Format
   ! ---------------------------------------------------------------------------
   names(1) = "CO2"
   names(2) = "C1-N2"
   names(3) = "C2-C3"
   names(4) = "C4"
   names(5) = "C5"
   names(6) = "C6"
   names(7) = "C7+n"
   names(8) = "Asph"

   header = "         T                     P"
   header_cap = "         T                     Py                     Px                     Pcap"

   ! Add lnK1", "lnK2", ..., "lnKn" dynamically
   do i = 1, nc
      header = trim(adjustl(header)) // "          lnK" // trim(adjustl(itoa(i))) // "(" // trim(adjustl(names(i))) // ")"
      header_cap = trim(adjustl(header_cap)) // "          lnK" // trim(adjustl(itoa(i))) // "(" // trim(adjustl(names(i))) // ")"
   end do

   ! Define fix format for the data 
   fmt = "(F25.8, F25.8"                      ! Fix space for "T" y "P"
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

   !! names="CO2" "C1-N2" "C2-C3" "C4" "C5" "C6" "C7+n" "Asph"
   
   !! composition vector                        
   z = (/0.0246_pr, 0.3694_pr, 0.0752_pr, 0.0193_pr, 0.0157_pr, 0.0162_pr, 0.47145_pr, 0.00815_pr/)
   !! Critical Pressure
   pc = (/73.7900_pr, 45.8300_pr, 45.4100_pr, 37.5400_pr, 33.8000_pr, 32.9000_pr, 12.4600_pr, 12.2900_pr/)
   !! Critical Temperature
   tc = (/304.0390_pr, 189.4280_pr, 339.8720_pr, 419.8170_pr, 465.0940_pr, 507.3170_pr, 860.3720_pr, 1424.8170_pr/)
   !! Acentric Factor
   w = (/0.225000_pr, 0.008500_pr, 0.127100_pr, 0.187800_pr, 0.239700_pr, 0.275000_pr, 1.022000_pr, 1.441000_pr/)
   !! lij and kij matrix
   lij = 0.0_pr
   kij = 0.0_pr
   kij(2,7) = 0.053_pr
   kij(2:5,8) = 0.135_pr
   kij(7,2) = 0.053_pr
   kij(8,2:5) = 0.135_pr
      
   !! Capillary pressure variables
   rp = 100 !! [nm]
   rp = rp/1E9 !! [m]
   ang_cont = 1.0472 !! angulo cualquiera de 60º en radianes
   LJ_par = 0.244_pr*((Tc/Pc)**(1.0/3.0))  !! [nm]
   Parachor = 40.1684*(0.151-0.0464*w)*(tc**(13._pr/12._pr))/(pc**(5._pr/6._pr)) !! [cm^3/mol * (mN/m)^1/4]

   !! Model definition
   model_bulk = PengRobinson78(tc, pc, w, kij, lij)
   model_conf = PengRobinson78Nano(LJ_par, rp*1E9, tc, pc, w, kij, lij)  ! pore radius in nanometers

   ! ===========================================================================
   ! BUBBLE ENVELOPES
   ! ---------------------------------------------------------------------------
   !! initialize the phase envelope
   sat_point_bulk = saturation_temperature(model_bulk, z, P=15._pr, kind="bubble", t0=175._pr)
   
   !! Calculate the phase envelope
   envelope_bulk = pt_envelope_2ph(model_bulk, z, sat_point_bulk)
   
   !! writing envelope
   !! --------------------------------------------------------------------------
   ! Bulk envelope
   write(1, "(A)") trim(header)
   do i=1,size(envelope_bulk%points)
      write(1,fmt) envelope_bulk%points(i)%T, envelope_bulk%points(i)%P*conv_bar_kPa, &
      log(envelope_bulk%points(i)%y/envelope_bulk%points(i)%x)
   end do
   !write(1,*) envelope_bulk%cps%T, envelope_bulk%cps%P*conv_bar_kPa
   !write(3,*) envelope_bulk

   ! ===========================================================================
   ! DEW ENVELOPES
   ! ---------------------------------------------------------------------------
   !! initialize the phase envelope
   sat_point_bulk = saturation_temperature(model_bulk, z, P=0.5_pr, kind="dew", t0=800._pr)

   !! Calculate the phase envelope
   envelope_bulk = pt_envelope_2ph(model_bulk, z, sat_point_bulk)
   
   !! writing envelope
   !! --------------------------------------------------------------------------
   ! Bulk envelope
   write(2, "(A)") trim(header)
   do i=1,size(envelope_bulk%points)
      write(2,fmt) envelope_bulk%points(i)%T, envelope_bulk%points(i)%P*conv_bar_kPa, &
      log(envelope_bulk%points(i)%y/envelope_bulk%points(i)%x)
   end do
   !write(1,*) envelope_bulk%cps%T, envelope_bulk%cps%P*conv_bar_kPa
   write(4,*) envelope_bulk

   envelope_bulk = find_hpl(model_bulk, z, t0=1100._pr, p0=1200._pr)
   write(3,*) envelope_bulk

   contains
      ! Función para convertir números enteros a cadenas (para el índice i)
      character(len=10) function itoa(i)
         integer, intent(in) :: i
         write(itoa, '(I0)') i
      end function itoa

end program feos_fluid