program phase_diagram_nano
   !! Program for calculation of phase diagrams. 
   !use forsus, only: Substance, forsus_dir, forsus_default_dir
   use yaeos, only: pr, &
      SoaveRedlichKwong, PengRobinson76, PengRobinson78, RKPR, PengRobinson78Nano,&
      EquilibriumState, ArModel, PTEnvel2, &
      pt_envelope_2ph, saturation_pressure, saturation_temperature, k_wilson 

   implicit none

   ! ===========================================================================
   ! Variables definition
   ! ---------------------------------------------------------------------------
   integer, parameter :: nc=4  
   class(ArModel), allocatable :: model, model_nano ! Thermodynamic model to be used
   type(EquilibriumState) :: sat_point, sat_point_nano              ! Init
   type(PTEnvel2) :: envelope, envelope_nano                       ! PT Phase envelope
   real(pr) :: tc(nc), pc(nc), w(nc)                ! Component's critical constants
   real(pr) :: z(nc), kij(nc,nc), lij(nc,nc)        ! Termodynamic variables
   real(pr) :: rp, LJ_par(nc)                       ! Nano parameters
   ! ===========================================================================
   ! Compound definition
   ! ---------------------------------------------------------------------------
   z = [0.1_pr, 0.4_pr, 0.3_pr, 0.2_pr] !CH4, CO2, C4H10(Butano), C6H14(Hexano)
   tc = [190.56_pr, 304.13_pr, 452.2_pr, 507.9_pr]
   pc = [45.99_pr, 73.8_pr, 37.96_pr, 30.25_pr]
   w = [0.011_pr, 0.225_pr, 0.199_pr, 0.297_pr]
   LJ_par = [0.3758_pr, 0.300_pr, 0.443_pr, 0.481_pr] !nm 
   rp = 100_pr !nm

   kij = reshape([0.000_pr, 0.110_pr, 0.025_pr, 0.030_pr, &
                  0.110_pr, 0.000_pr, 0.120_pr, 0.140_pr, &
                  0.025_pr, 0.120_pr, 0.000_pr, 0.015_pr, &
                  0.030_pr, 0.140_pr, 0.015_pr, 0.000_pr], [nc,nc]) 
   lij = kij / 2 
   
   ! Model definition
   model = PengRobinson78(tc, pc, w, kij, lij)
   model_nano = PengRobinson78Nano(LJ_par, rp, tc, pc, w, kij, lij)
   
   ! Calculate a dew point at low pressure to later 
   ! initialize the phase envelope
   ! sat_point = saturation_temperature(model, z, P=1._pr, kind="dew", t0=300._pr)
   ! sat_point_nano = saturation_temperature(model_nano, z, P=1._pr, kind="dew", t0=300._pr)
   sat_point = saturation_pressure(model, z, P0=500._pr, kind="bubble", t=300._pr)
   sat_point_nano = saturation_pressure(model_nano, z, P0=500._pr, kind="bubble", t=300._pr)
   ! Calculate phase envelope
   envelope = pt_envelope_2ph(model, z, sat_point)
   envelope_nano = pt_envelope_2ph(model_nano, z, sat_point_nano, 500)

   write(1, "(*(A,2x))") "kind","T", "P", "beta","x", "y", "Vx", "Vy"
   write(*, *) envelope%points(1)
   write(1,*) envelope
   
   write(2, "(*(A,2x))") "kind","T", "P", "beta","x", "y", "Vx", "Vy"
   write(*, *) envelope_nano%points(1)
   write(2,*) envelope_nano

end program phase_diagram_nano
