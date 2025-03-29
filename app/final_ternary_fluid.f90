program Ternary_fluid
   use yaeos, only: pr, ArModel, PTEnvel2, pt_envelope_2ph, EquilibriumState, &
   NanoEquilibriumState, nano_pt_envelope_2ph, NanoPTEnvel2, SoaveRedlichKwong, &
   saturation_pressure, saturation_temperature 
   use Capillary_Initializer, only: capillary_init_point
   implicit none

   integer, parameter :: nc=3
   class(ArModel), allocatable :: model ! Thermodynamic model to be used
   type(EquilibriumState) :: sat_point   ! Init bulk
   type(NanoEquilibriumState) :: init_point ! Init Nano

   type(PTEnvel2) :: envelope           ! PT Phase envelope
   type(NanoPTEnvel2) :: envelope_cap  ! PT Phase nano envelope
   real(pr) :: tc(nc), pc(nc), w(nc), z(nc), kij(nc,nc), lij(nc,nc) ! Component's values
   ! Capilar's values
   real(pr) :: rp, ang_cont
   real(pr), allocatable :: Parachor(:)
   integer :: i

   character(len=500) :: header, header_cap, names(nc)  ! Cadena para almacenar el encabezado
   character(len=1000) :: fmt, fmt_cap


   !! Formato de curvas

   !! names="CH4" "C50-C60" "Asph"

   names(1) = "CH4"
   names(2) = "C50-C60"
   names(3) = "Asph"
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



   !! Values loading
   !! composition vector                        
   z = (/0.621476711621686, 0.294896838294686, 0.0836264500836264/)
   !! Critical Temperature
   tc = (/190.56, 956.95, 1286.75/)
   !! Critical Pressure
   pc = (/45.99, 13.72, 18.11/)
   !! Acentric Factor
   w = (/0.0115, 1.313, 1.274/)
   lij = 0
   kij = 0
   Kij(1,3) = 0.017
   Kij(3,1) = 0.017
   
   !! Capillary pressure variables
   rp = 10 !! [nm]
   rp = rp/1E9 !! [m]
   ang_cont = 1.0472 !! angulo cualquiera de 60º en radianes
   Parachor = (/0.07405, 2.76, 3.2/)*1E3 !! [cm^3/mol * (mN/m)^1/4]


   !! model definition
   model = SoaveRedlichKwong(tc, pc, w, kij, lij)

   !! ------------------------- BUBBLE ENVELOPE ------------------
   !! initialize the phase envelope
   init_point = capillary_init_point(model, z, Parachor, ang_cont, rp, "bubble",&
   P0_in=900._pr, T_in=400._pr)
   
   sat_point = saturation_pressure(model, z, P0=900._pr, kind="bubble", t=400._pr)
   
   write(*,*) "----------------- Init Capillary Point Bubble------------------------"
   write(*,*) "iteraciones ",init_point%iters
   write(*,*) "kind ",init_point%kind
   write(*,*) "T inicial ",init_point%T
   write(*,*) "Vx inicial ",init_point%Vx
   write(*,*) "Vy inicial ",init_point%Vy
   write(*,*) "K inicial ",(init_point%y/init_point%x)
   write(*,*) "Final Point x:", init_point%x
   write(*,*) "Final Point y:", init_point%y
   write(*,*) "Pcap inicial ",init_point%Pcap
   write(*,*) "Px inicial ",init_point%Px
   write(*,*) "Py inicial ",init_point%Py
   !init_point%

   !! envelope calculation
   envelope_cap = nano_pt_envelope_2ph(model, z, rp, ang_cont,&
   Parachor, init_point, points=5000, iterations=1000)

   envelope = pt_envelope_2ph(model, z, sat_point, points=200)

   !!--------------writing-------------------

   write(1, "(A)") trim(header)
   do i=1,size(envelope%points)
      write(1,fmt) envelope%points(i)%T, envelope%points(i)%P, &
      log(envelope%points(i)%y/envelope%points(i)%x)
   end do
   write(*,*)"bulk bubble cps:", envelope%cps

   write(3, "(A)") trim(header_cap)
   do i=1,size(envelope_cap%points)
      write(3,fmt_cap) envelope_cap%points(i)%T, envelope_cap%points(i)%Py, &
      envelope_cap%points(i)%Px, envelope_cap%points(i)%Pcap,&
      log(envelope_cap%points(i)%y/envelope_cap%points(i)%x)
   end do
   write(*,*)"cap bubble cps:", envelope_cap%cps


   !! -------------------------- DEW ENVELOPE ---------------------------
   
   !! initialize the phase envelope
   init_point = capillary_init_point(model, z, Parachor, ang_cont, rp, "dew",&
   P_in=0.5_pr, T0_in=900._pr)

   sat_point = saturation_temperature(model, z, P=0.5_pr, kind="dew", t0=900._pr)


   write(*,*) "----------------- Init Capillary Point Dew------------------------"
   write(*,*) "iteraciones ",init_point%iters
   write(*,*) "kind ",init_point%kind
   write(*,*) "T inicial ",init_point%T
   write(*,*) "Vx inicial ",init_point%Vx
   write(*,*) "Vy inicial ",init_point%Vy
   write(*,*) "K inicial ",(init_point%y/init_point%x)
   write(*,*) "Final Point x:", init_point%x
   write(*,*) "Final Point y:", init_point%y
   write(*,*) "Pcap inicial ",init_point%Pcap
   write(*,*) "Px inicial ",init_point%Px
   write(*,*) "Py inicial ",init_point%Py
   
   !! envelope calculation
   envelope_cap = nano_pt_envelope_2ph(model, z, rp, ang_cont,&
   Parachor, init_point, points=5000, iterations=1000)
   
   envelope = pt_envelope_2ph(model, z, sat_point, points=100)

   !!--------------writing-------------------

   write(2, "(A)") trim(header)
   do i=1,size(envelope%points)
      write(2,fmt) envelope%points(i)%T, envelope%points(i)%P, &
      log(envelope%points(i)%y/envelope%points(i)%x)
   end do
   write(*,*)"bulk dew cps:", envelope%cps


   write(4, "(A)") trim(header_cap)
   do i=1,size(envelope_cap%points)
      write(4,fmt_cap) envelope_cap%points(i)%T, envelope_cap%points(i)%Py, &
      envelope_cap%points(i)%Px, envelope_cap%points(i)%Pcap,&
      log(envelope_cap%points(i)%y/envelope_cap%points(i)%x)
   end do
   write(*,*)"cap dew cps:", envelope_cap%cps

   contains
      character(len=10) function itoa(i)
         integer, intent(in) :: i
         write(itoa, '(I0)') i
      end function itoa
end program Ternary_fluid