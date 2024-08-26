module yaeos__equilibria_equilibrium_state
   use yaeos__constants, only: pr
   implicit none
      type :: EquilibriumState
      !! Description of a two-phase equilibria state.
      !!
      !! Contains the relevant information of an equilibrium point obtained
      !! from some kind of equilibria calculation.
      character(len=14) :: kind
      !! Kind of point ["bubble", "dew", "liquid-liquid", "split"]
      integer :: iters = 0
      !! Iterations needed to reach the state
      real(pr), allocatable :: y(:)
      !! Light-phase molar fractions
      real(pr), allocatable :: x(:)
      !! Heavy-phase molar fractions
      real(pr) :: Vx
      !! Heavy-phase volume [L/mol]
      real(pr) :: Vy
      !! Light-phase volume [L/mol]
      real(pr) :: T
      !! Temperature [K]
      real(pr) :: P
      !! Pressure [bar]
      real(pr) :: beta
      !! Mole fraction of light-phase
      integer :: ns
   contains
      private
      procedure, pass :: write => write_EquilibriumState
      generic, public :: write (FORMATTED) => write
   end type EquilibriumState
   type :: NanoEquilibriumState
      !! Description of a two-phase equilibria state.
      !!
      !! Contains the relevant information of an equilibrium point obtained
      !! from some kind of equilibria calculation.
      character(len=14) :: kind
      !! Kind of point ["bubble", "dew", "liquid-liquid", "split"]
      integer :: iters = 0
      !! Iterations needed to reach the state
      real(pr), allocatable :: y(:)
      !! Light-phase molar fractions
      real(pr), allocatable :: x(:)
      !! Heavy-phase molar fractions
      real(pr) :: Vx
      !! Heavy-phase volume [L/mol]
      real(pr) :: Vy
      !! Light-phase volume [L/mol]
      real(pr) :: T
      !! Temperature [K]
      real(pr) :: Pcap
      !! Capilar Pressure [bar]
      real(pr) :: Px
      !! Heavy-phase Pressure [bar]
      real(pr) :: Py
      !! Light-phase Pressure [bar]
      real(pr) :: beta
      !! Mole fraction of light-phase
      integer :: ns
      !! specification value
   contains
      private
      procedure, pass :: write => write_NanoEquilibriumState
      generic, public :: write (FORMATTED) => write
   end type NanoEquilibriumState
contains

   subroutine write_EquilibriumState(eq, unit, iotype, v_list, iostat, iomsg)
      class(EquilibriumState), intent(in) :: eq
      integer, intent(in) :: unit
      character(*), intent(in) :: iotype
      integer, intent(in)  :: v_list(:)
      integer, intent(out) :: iostat
      character(*), intent(inout) :: iomsg

      character(*), parameter :: nl = new_line("G")
      character(len=16) :: cadenaT, cadenaP

      !write(unit, *) eq%kind, eq%T, eq%P, eq%beta, eq%x, eq%y, eq%Vx, eq%Vy
      !write(unit,*)  eq%T, eq%P
      write(cadenaT,"(F16.8)")  eq%T
      write(cadenaP,"(F16.8)")  eq%P
      write(unit,*) cadenaT, cadenaP, "  iters", eq%iters, " ns :", eq%ns, log(eq%y/eq%x)

   end subroutine write_EquilibriumState
   subroutine write_NanoEquilibriumState(eq, unit, iotype, v_list, iostat, iomsg) 
      class(NanoEquilibriumState), intent(in) :: eq
      integer, intent(in) :: unit
      character(*), intent(in) :: iotype
      integer, intent(in)  :: v_list(:)
      integer, intent(out) :: iostat
      character(*), intent(inout) :: iomsg

      character(*), parameter :: nl = new_line("G")
      character(len=16) :: cadenaT, cadenaPy, cadenaPx, cadenaPcap

      !write(unit, *) eq%kind, eq%T, eq%Py, eq%Px, eq%Pcap, eq%beta, eq%x, eq%y, eq%Vx, eq%Vy
      !write(cadenaT,"(F8.4)")  eq%T, eq%Py, eq%Px, eq%Pcap
      write(cadenaT,"(F16.8)")  eq%T
      write(cadenaPy,"(F16.8)")  eq%Py
      write(cadenaPx,"(F16.8)")  eq%Px
      write(cadenaPcap,"(F16.8)")  eq%Pcap
      write(unit,*) cadenaT, ",", cadenaPy, ",", cadenaPx, ",", cadenaPcap, "  iters", eq%iters, " ns :", eq%ns, log(eq%y/eq%x)

   end subroutine write_NanoEquilibriumState
end module yaeos__equilibria_equilibrium_state
