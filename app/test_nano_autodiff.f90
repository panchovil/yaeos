module hyperdual_pr78_nano
   use yaeos__adiff_hyperdual_ar_api, only: ArModelAdiff, hyperdual
   use hyperdual_mod
   use yaeos__constants, only: pr, R
   use yaeos__substance, only: Substances
   implicit none

   type, extends(ArModelAdiff) :: PR78_nano_autodiff
      !! PengRobinson 78 EoS with adsortion and fluid/wall interaction

      ! Mixing rule Parameters
      real(pr), allocatable :: kij(:, :), lij(:, :)

      ! EoS parameters
      real(pr), allocatable :: ac(:), b(:), k(:)
      real(pr), allocatable :: tc(:), pc(:), w(:)
      real(pr), allocatable :: alpha_ads_i(:), rp, LJ_par(:)
   contains
      procedure :: Ar => arfun
      procedure :: get_v0 => v0
      procedure :: volume => volume
   end type PR78_nano_autodiff

   real(pr), parameter :: del1 = 1._pr + sqrt(2._pr)
   real(pr), parameter :: del2 = 1._pr - sqrt(2._pr)

contains

   type(PR78_nano_autodiff) function setup(LJ_par, rp, tc, pc, w, kij, lij) result(self)
        !! Function to obtain a defined PR78 model modified with setted up parameters
        !! as function of Tc, Pc, w, LJ_par and rp
      real(pr) :: tc(:)
      real(pr) :: pc(:)
      real(pr) :: w(:)
      real(pr) :: kij(: , :)
      real(pr) :: lij(: , :)
      real(pr) :: rp 
      real(pr) :: LJ_par(:)

      allocate(self%components%tc(size(tc)))
      allocate(self%components%pc(size(tc)))
      allocate(self%components%w(size(tc)))
      allocate(self%kij(size(tc), size(tc)))
      allocate(self%lij(size(tc), size(tc)))
      allocate(self%ac(size(tc)))
      allocate(self%b(size(tc)))
      allocate(self%k(size(tc)))
      allocate(self%alpha_ads_i(size(tc)))
      allocate(self%LJ_par(size(tc)))

      self%components%tc = tc
      self%components%pc = pc
      self%components%w = w
      
      self%LJ_par = LJ_par
      self%rp = rp
      self%alpha_ads_i = (1-(0.7597_pr*((rp/LJ_par)**-0.7708_pr)))/(1-(0.9793_pr*((rp/LJ_par)**(-0.6366_pr))))
      self%ac = 0.45723553_pr * R**2 * (self%components%tc**2 / self%components%pc) * &
      (1-(0.7597_pr*((rp/LJ_par)**-0.7708_pr)))**2/(1-(0.9793_pr*((rp/LJ_par)**(-0.6366_pr))))
      self%b = 0.07779607_pr * R * self%components%tc/(self%components%pc * self%alpha_ads_i)
        
      where (self%components%w <= 0.491)
         self%k = 0.37464 + 1.54226 * self%components%w - 0.26992 * self%components%w**2
      elsewhere
         self%k = 0.379642 + 1.48503 * self%components%w - 0.164423 * self%components%w**2&
          + 0.016666 * self%components%w**3
      end where

      self%kij = kij
      self%lij = lij
   end function

   function arfun(self, n, v, t) result(ar)
      !! Residual Helmholtz calculation for a generic cubic with
      !! quadratic mixing rules.
      class(PR78_nano_autodiff) :: self
      type(hyperdual), intent(in) :: n(:), v, t
      type(hyperdual) :: ar

      type(hyperdual) :: amix, a(size(n)), ai(size(n)), n2(size(n))
      type(hyperdual) :: bmix
      !type(hyperdual) :: b_v, nij
      type(hyperdual) :: nij
      type(hyperdual) :: alpha_ads_mix

      integer :: i, j

      ! Associate allows us to keep the later expressions simple.
      associate(&
         pc => self%components%pc, ac => self%ac, b => self%b, k => self%k,&
         kij => self%kij, lij => self%lij, tc => self%components%tc, & 
         alpha_ads_i => self%alpha_ads_i) 

         ! Soave alpha function
         a = 1.0_pr + k * (1.0_pr - sqrt(t/tc))
         a = ac * a ** 2
         ai = sqrt(a)

         ! Quadratic Mixing Rule
         amix = 0.0_pr
         bmix = 0.0_pr

         do i=1,size(n)-1
               do j=i+1,size(n)
                  nij = n(i) * n(j)
                  amix = amix + 2 * nij * (ai(i) * ai(j)) * (1 - kij(i, j))
                  bmix = bmix + nij * (b(i) + b(j)) * (1 - lij(i, j))
               end do
         end do
      

         amix = amix + sum(n**2*a)
         bmix = bmix + sum(n**2 * b)

         bmix = bmix/sum(n)

         ! Linear Mixing Rule
         alpha_ads_mix = 0.0_pr
         alpha_ads_mix = sum(n * alpha_ads_i)
         alpha_ads_mix = alpha_ads_mix/sum(n)
         !b_v = bmix/v

         ! Generic Nano Cubic Ar function
         ar = (- sum(n) * T * (R*(((1.0_pr/alpha_ads_mix)-1.0_pr)*log(V)+&
         (log(alpha_ads_mix*V-Bmix)-log(alpha_ads_mix*V))/alpha_ads_mix)) &
         - amix * (log((alpha_ads_mix*V + del1*Bmix)/(alpha_ads_mix*V + &
         del2*Bmix))/(Bmix*alpha_ads_mix*(del1 - del2))))
         
         
         ! ar = (&
         !    - sum(n) * log(1.0_pr - b_v) &
         !    - amix / (R*T*bmix)*1.0_pr / (del1 - del2) &
         !    * log((1.0_pr + del1 * b_v) / (1.0_pr + del2 * b_v)) &
         !    ) * (R * T)

      end associate
   end function arfun

   function v0(self, n, p, t)
      !! Initialization of liquid volume solving with covolume. This also
      !! helps the Michelsen volume solver
      class(PR78_nano_autodiff), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: p
      real(pr), intent(in) :: t
      real(pr) :: v0

      v0 = sum(n * self%b) / sum(n)
   end function v0

   subroutine volume(eos, n, P, T, V, root_type)
      !! In the case of models that have a "covolume" value, using the solver
      !! of Michelsen is a better option that the default.
      use yaeos__models_solvers, only: volume_michelsen
      class(PR78_nano_autodiff), intent(in) :: eos
      real(pr), intent(in)  :: n(:), P, T
      real(pr), intent(out) :: V
      character(len=*), intent(in) :: root_type

      call volume_michelsen(eos, n, P, T, V, root_type)
   end subroutine
end module hyperdual_pr78_nano


program test_nano_autodiff
   use hyperdual_pr78_nano
   implicit none
   class(PR78_nano_autodiff), allocatable :: auto_eos
   
   integer, parameter :: n=4
   ! real(pr) :: z(n)
   ! real(pr) :: v=1.0, t=150.0, p
   type(hyperdual) :: v, t, p , z(n), ar_result
   
   real(pr) :: tc(size(z)), pc(size(z)), w(size(z)), kij(size(z), size(z)), lij(size(z), size(z)), LJ_par(size(z)), rp



!CH4, CO2, C4H10(Butano), C6H14(Hexano)
   tc = [190.56_pr, 304.13_pr, 452.2_pr, 507.9_pr]
   pc = [45.99_pr, 73.8_pr, 37.96_pr, 30.25_pr]
   w = [0.011_pr, 0.225_pr, 0.199_pr, 0.297_pr]
   LJ_par = [0.3758_pr, 0.300_pr, 0.443_pr, 0.481_pr] !nm 
   rp = 100.0_pr !nm

   kij = reshape([0.000_pr, 0.110_pr, 0.025_pr, 0.030_pr, &
                  0.110_pr, 0.000_pr, 0.120_pr, 0.140_pr, &
                  0.025_pr, 0.120_pr, 0.000_pr, 0.015_pr, &
                  0.030_pr, 0.140_pr, 0.015_pr, 0.000_pr], [n,n]) 
   lij = kij / 2 

   auto_eos = setup(LJ_par=LJ_par, rp=rp, tc=tc, pc=pc, w=w, kij=kij, lij=lij)
   
   z = [0.1_pr, 0.4_pr, 0.3_pr, 0.2_pr] 
   t = 150._pr
   v = 1.0_pr
   ar_result = auto_eos%Ar(z, v, t)
   !arfun(auto_eos, z, v, t)

   print*, ar_result
end program test_nano_autodiff