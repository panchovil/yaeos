module Capillary_Initializer
   !! # Initializer
   !! The objective of this module is to initialize the capillary envelopes with a fixed
   !! pore radius, where confinement effects are minimal. Then, it reduces the pore radius
   !! step by step, calculating the converged point at each iteration until it reaches
   !! the required pore radius.
   
   
   ! use yaeos, only: pr, R, &
   !    SoaveRedlichKwong, PengRobinson76, PengRobinson78, RKPR, PengRobinson78Nano,&
   !    EquilibriumState, ArModel, PTEnvel2, &
   !    pt_envelope_2ph, saturation_pressure, saturation_temperature, k_wilson ,&
   !    NanoEquilibriumState, nano_pt_envelope_2ph, NanoPTEnvel2
   ! use yaeos, only: pr, R, &
   !     EquilibriumState, ArModel, PTEnvel2, &
   !     pt_envelope_2ph, saturation_pressure, saturation_temperature, k_wilson ,&
   !     NanoEquilibriumState, nano_pt_envelope_2ph, NanoPTEnvel2
   use yaeos__equilibria, only: EquilibriumState, saturation_pressure, saturation_temperature,&
        NanoEquilibriumState
   use yaeos__constants, only: pr, R
   use yaeos__equilibria_boundaries_nano_phase_envelopes_pt, only:&
      NanoPTEnvel2, nano_pt_envelope_2ph
   use yaeos__equilibria_boundaries_phase_envelopes_pt, only:&
      PTEnvel2, find_hpl
   use yaeos__models, only: ArModel
   
   implicit none
   
   
   ! class(ArModel), allocatable :: model, model_nano ! Thermodynamic model to be used
   ! type(EquilibriumState) :: sat_point, sat_point_nano!, sat_point_auto_eos            ! Init
   ! type(NanoEquilibriumState) :: init_point_cap, init_point_nano_cap
   ! !type(PTEnvel2) :: envelope, envelope_nano!, envelope_auto_eos
   ! type(NanoPTEnvel2) :: envelope_cap ,envelope_nano_cap   
   


   contains

      type(NanoEquilibriumState) function capillary_init_point&
      (model, z_in, Par_in, ang_cont, rp_in, kind_in, P_in, T_in, P0_in, T0_in)!, type_sat
         
         character(len=*), intent(in) :: kind_in!, type_sat
         real(pr), optional, intent(in) :: P_in, T_in, P0_in, T0_in

         class(ArModel), intent(in) :: model
         real(pr), intent(in) :: rp_in, ang_cont, z_in(:), Par_in(:)
         type(EquilibriumState) :: auxiliary_sat_point
         type(NanoEquilibriumState) :: auxiliary_point, final_point
         type(NanoPTEnvel2) :: auxiliary_env
         type(PTEnvel2) :: auxiliary_hpl_env

         real(pr) :: rp_init =  1E-6 !! [m] equivalente a 1000 nm
         real(pr) :: rp_step =  0.9  !! paso de reduccion
         real(pr) :: rp_aux          !! radio de poro auxiliar para iteracion
         integer :: final_iters, mid_iters !! final iters tiene las iteracciones totales de todos los puntos calculados, mid iters la cantidad de veces que se hizo un calculo de punto

         final_iters = 0
         mid_iters = 0

         print*, "Start Laplace"
         !! se calcula el punto de saturacion sin capilaridad
         if (present(P0_in) .and. present(T_in)) then
            auxiliary_sat_point = saturation_pressure(model, z_in, P0=P0_in, kind=kind_in, t=T_in)
         else if (present(P_in) .and. present(T0_in)) then
            auxiliary_sat_point = saturation_temperature(model, z_in, P=P_in, kind=kind_in, t0=T0_in)
         else if (present(P0_in) .and. present(T0_in)) then
            auxiliary_hpl_env = find_hpl(model, z_in, T0_in, P0_in, 1)
            auxiliary_sat_point = auxiliary_hpl_env%points(1) !!no es un punto de saturacion pero bueno
         else
            print*, "Error: Los parámetros ingresados no son válidos.&
             Asegúrese de ingresar P0 y T o P y T0 o kind correctamente."
            stop 1
         end if

         mid_iters = mid_iters + 1
         final_iters = final_iters + auxiliary_sat_point%iters

         !! se calcula el "punto de saturacion capilar" con el radio de poro con bajo efecto de confinamiento
         call Laplace_init(auxiliary_sat_point, Par_in, ang_cont, rp_init, auxiliary_point)
         
         !! se ajusta el primer punto de saturacion con full Newton
         auxiliary_env = nano_pt_envelope_2ph(model, z_in, rp_init, ang_cont,&
         Par_in, auxiliary_point, points=1, iterations=1000)
         auxiliary_point = auxiliary_env%points(1)

         final_iters = final_iters + auxiliary_point%iters
         mid_iters = mid_iters + 1

         !! Iteración reduciendo el radio de poro  
         rp_aux = rp_init
         do while (rp_aux > rp_in)  
            
            mid_iters = mid_iters + 1

            rp_aux = rp_aux * rp_step  ! se reduce el radio  
            print*, rp_aux*1E9

            auxiliary_env = nano_pt_envelope_2ph(model, z_in, rp_aux, ang_cont, &  
            Par_in, auxiliary_point, points=1, iterations=1000)  
            auxiliary_point = auxiliary_env%points(1)  

            final_iters = final_iters + auxiliary_point%iters
            !print*, rp_aux
         end do  
         print*, "radio de poro en m del ultimo valor antes de correcion con rp de entrada", rp_aux !! antes de correcion de rp
         
         !! Ultimo calculo para asegurar el punto con el radio de poro de entrada
         auxiliary_env = nano_pt_envelope_2ph(model, z_in, rp_in, ang_cont, &  
         Par_in, auxiliary_point, points=1, iterations=1000)  

         mid_iters = mid_iters + 1

         !! Corrección: Allocación de `final_point` antes de la asignación final
         if (.not. allocated(final_point%x)) allocate(final_point%x(size(auxiliary_env%points(1)%x)))
         if (.not. allocated(final_point%y)) allocate(final_point%y(size(auxiliary_env%points(1)%y))) 

         final_point = auxiliary_env%points(1)  
         final_iters = final_iters + final_point%iters


         ! Debugging: Ver valores antes de retornar
         ! print*, "Final Point T:", final_point%T
         ! print*, "Final Point Px:", final_point%Px
         ! print*, "Final Point Py:", final_point%Py
         ! print*, "Final Point x:", final_point%x
         ! print*, "Final Point y:", final_point%y

         print*, "iteraciones totales", final_iters
         print*, "cantidad de puntos calculados", mid_iters


         capillary_init_point = final_point

      end function capillary_init_point
      subroutine Laplace_init(sat_point_in, Par_in, ang_cont_in, rp_in, init_point_out)
         ! use yaeos, only: pr, EquilibriumState, NanoEquilibriumState
         ! implicit none
         type(EquilibriumState), intent(in) :: sat_point_in
         type(NanoEquilibriumState), intent(out) :: init_point_out
         real(pr), intent(in) :: rp_in, ang_cont_in, Par_in(:)
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
end module Capillary_Initializer