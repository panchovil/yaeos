program sandoval_phase_test
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
   integer, parameter :: nc=7  
   integer :: i
   class(ArModel), allocatable :: model, model_nano ! Thermodynamic model to be used
   type(EquilibriumState) :: sat_point, sat_point_nano              ! Init
   type(PTEnvel2) :: envelope, envelope_nano                       ! PT Phase envelope
   real(pr) :: tc(nc), pc(nc), w(nc)                ! Component's critical constants
   real(pr) :: z(nc), kij(nc,nc), lij(nc,nc)        ! Termodynamic variables
   real(pr) :: rp, LJ_par(nc)                       ! Nano parameters
   real(pr) :: pcc(nc), delP(nc), tcc(nc), delT(nc)
   ! ===========================================================================
   ! Compound definition
   ! ---------------------------------------------------------------------------
   !! names="N2" "CH4" "C2H6" "C3H8" "nC4" "nC5" "nC6"
   !! composition vector                        
   z = (/0.014, 0.943, 0.027, 0.0074, 0.0049, 0.0027, 0.001/)
   !! Critical Pressure
   pc = (/34.0, 45.99, 48.72, 42.48, 37.96, 33.70, 30.25/)
   !! Critical Temperature
   tc = (/126.20, 190.56, 305.32, 369.83, 425.12, 469.70, 507.60/)
   !! Acentric Factor
   w = (/0.0377, 0.0115, 0.0995, 0.1523, 0.2002, 0.2515, 0.3013/)
   lij = 0
   kij = 0
   Kij(1,2) = 0.0278
   Kij(1,3) = 0.4070
   Kij(1,4) = 0.0763
   Kij(1,5) = 0.0700
   Kij(1,6) = 0.0787
   Kij(1,7) = 0.1496

   Kij(2,1) = 0.0278
   Kij(2,3) = -0.0078
   Kij(2,4) = 0.0090
   Kij(2,5) = 0.0056
   Kij(2,6) = 0.0190
   Kij(2,7) = 0.0374
   
   rp = 150 !nm
   LJ_par = (/0.364_pr, 0.375_pr, 0.443_pr, 0.472_pr, 0.506_pr, 0.529_pr, 0.550_pr/) 
   
   !! Critical Pressure shift
   delP = 0.9793*((rp/LJ_par)**-0.6366)
   !! Critical Pressure modified
   pcc = -((delP*pc)-pc)
   !! Critical Temperature shift
   delT = 0.7597*((rp/LJ_par)**-0.7708)
   !! Critical Pressure modified
   tcc = -((delT*tc)-tc)

   ! Model definition
   model = PengRobinson78(tc, pc, w, kij, lij)
   model_nano = PengRobinson78Nano(LJ_par, rp, tc, pc, w, kij, lij)
   
   ! Calculate a dew point at low pressure to later 
   ! initialize the phase envelope
   ! sat_point = saturation_temperature(model, z, P=1._pr, kind="dew", t0=300._pr)
   ! sat_point_nano = saturation_temperature(model_nano, z, P=1._pr, kind="dew", t0=300._pr)
   print*, "pc", pc
   print*, "pcc", pcc
   print*, "tc", tc
   print*, "tcc", tcc

   ! ------------------------- BUBBLE ENVELOPE ------------------

   sat_point = saturation_temperature(model, z, P=1.5_pr, kind="bubble", t0=115._pr)
   sat_point_nano = saturation_temperature(model_nano, z, P=1.5_pr, kind="bubble", t0=115._pr)
   ! Calculate phase envelope
   envelope = pt_envelope_2ph(model, z, sat_point)
   envelope_nano = pt_envelope_2ph(model_nano, z, sat_point_nano, points=500, iterations=1000)

   !write(1, "(*(A,2x))") "kind","T", "P", "beta","x", "y", "Vx", "Vy"
   !write(*, *) envelope%points(1)
   do i=1,size(envelope%points)
      write(1,*) envelope%points(i)%T ,envelope%points(i)%P, log(envelope%points(i)%y/envelope%points(i)%x)
   end do
   write(1,*), envelope%cps
   !write(2, "(*(A,2x))") "kind","T", "P", "beta","x", "y", "Vx", "Vy"
   !write(*, *) envelope_nano%points(28)
   do i=1,size(envelope_nano%points)
      write(3,*) envelope_nano%points(i)%T ,envelope_nano%points(i)%P
   end do
   
   !write(2,*) envelope_nano
   !print*, size(envelope_nano%points)

   !! -------------------------- DEW ENVELOPE ---------------------------

   sat_point = saturation_temperature(model, z, P=0.5_pr, kind="dew", t0=200._pr)
   sat_point_nano = saturation_temperature(model_nano, z, P=0.5_pr, kind="dew", t0=200._pr)

   envelope = pt_envelope_2ph(model, z, sat_point)
   envelope_nano = pt_envelope_2ph(model_nano, z, sat_point_nano, points=500, iterations=500)

   ! write(3, "(*(A,2x))") "kind","T", "P", "beta","x", "y", "Vx", "Vy"
   ! write(*, *) envelope%points(1)
   ! write(3,*) envelope
   
   ! write(4, "(*(A,2x))") "kind","T", "P", "beta","x", "y", "Vx", "Vy"
   ! write(*, *) envelope_nano%points(1)
   ! write(4,*) envelope_nano

   do i=1,size(envelope%points)
      write(2,*) envelope%points(i)%T ,envelope%points(i)%P
   end do

   do i=1,size(envelope_nano%points)
      write(4,*) envelope_nano%points(i)%T ,envelope_nano%points(i)%P
   end do
   

end program sandoval_phase_test