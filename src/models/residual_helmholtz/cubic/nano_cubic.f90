module yaeos__models_ar_nanocubic
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   use yaeos__substance, only: Substances
   use yaeos__models_ar_genericcubic, only: AlphaFunction, CubicEoS, v0, CubicMixRule, volume
   use yaeos__models_ar_cubic_quadratic_mixing, only: QMR, Dmix, D1mix_constant, Bmix
   implicit none

   ! type, extends (QMR) :: CubicMixRuleNano
   ! !! blablablabla
   ! contains
   !    procedure :: alpha_ads_mix_linear
   ! end type CubicMixRuleNano
   type, extends(QMR) :: CubicMixRuleNano
   !! blablablabla
   contains
      procedure :: Dmix => Dmix_wrapper !! Attractive parameter mixing rule
      procedure :: Bmix => Bmix_wrapper !! Repulsive parameter mixing rule
      procedure :: D1mix => D1mix_constant_wrapper
      procedure :: alpha_ads_mix_linear
   end type CubicMixRuleNano


   type, extends(CubicEoS) :: CubicEoSNano
       
   !! blablablabla
      !class(CubicMixRuleNano),allocatable :: mixrule_nano
      real(pr), allocatable :: alpha_ads_i(:)
   contains
      procedure :: residual_helmholtz => GenericCubic_Ar_Nano
      procedure :: get_v0 => v0_wrapper
      !procedure :: volume => volume
      !procedure :: init_nano
   end type CubicEoSNano
      ! type, extends(CubicEoS) :: CubicEoSNano
   !    class(CubicMixRuleNano), allocatable :: mixrule
   !    real(pr), allocatable :: alpha_ads(:) !! 

   ! contains
   !    procedure :: residual_helmholtz => GenericCubic_Ar_Nano
   ! end type
contains

   ! subroutine init_nano(self)
   !    class(CubicEoSNano), intent(inout) :: self


   !    ! Inicializa `mixrule` como una instancia de CubicMixRuleNano
   !    allocate(self%mixrule_nano, source=CubicMixRuleNano())

   !    ! Inicializa otros atributos específicos de CubicEoSNano si es necesario
   !    ! Ejemplo:
   !    ! allocate(self%alpha_ads_i(size_variable))
   ! end subroutine init_nano

   subroutine alpha_ads_mix_linear(self, n, alpha_ads_i, alpha_ads_mix, dalpha_ads_mixi, dalpha_ads_mixij)
      class(CubicMixRuleNano), intent(in) :: self  
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: alpha_ads_i(:)
      real(pr), intent(out) :: alpha_ads_mix, dalpha_ads_mixi(:), dalpha_ads_mixij(:, :)
      !real(pr) :: bij(size(n), size(n))

      real(pr) :: totn

      integer :: i, j, nc

      nc = size(n)
      TOTN = sum(n)

      alpha_ads_mix = 0
      dalpha_ads_mixi = 0
      dalpha_ads_mixij = 0

      alpha_ads_mix = sum(n*alpha_ads_i)
      alpha_ads_mix = alpha_ads_mix/TOTN

      do i = 1,nc
         dalpha_ads_mixi(i) = (alpha_ads_i(i) - alpha_ads_mix)/TOTN
      end do
      
      do i = 1, nc
         do j = 1, nc
            dalpha_ads_mixij(i, j) = (-dalpha_ads_mixi(i)-dalpha_ads_mixi(j))/TOTN
         end do
      end do
   end subroutine alpha_ads_mix_linear

   subroutine GenericCubic_Ar_Nano(&
      self, n, V, T, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2&
      )
      !! Residual Helmholtz Energy for a generic Cubic Equation of State.
      !!
      !! Calculates the residual Helmholtz Energy for a generic Cubic EoS as
      !! defined by Michelsen and Møllerup:
      !!
      !! \[
      !!   P = \frac{RT}{V-b} 
      !!       - \frac{a_c\alpha(T_r)}{(V+b\delta_1)(V+b\delta_2)}
      !! \]
      !!
      !! This routine assumes that the \(\delta_1\) is not a constant parameter
      !! (as it uses to be in classical Cubic EoS) to be compatible with the
      !! three parameter EoS RKPR where \(delta_1\) is not a constant and
      !! has its own mixing rule.
      !!
      use yaeos__constants, only: R
      class(CubicEoSNano), intent(in) :: self
      real(pr), intent(in) :: n(:) !! Number of moles
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: t !! Temperature [K]

      real(pr), optional, intent(out) :: ar !! Residual Helmholtz
      real(pr), optional, intent(out) :: arv !! \(\frac{dAr}{dV}\)
      real(pr), optional, intent(out) :: ArT !! \(\frac{dAr}{dT}\)
      real(pr), optional, intent(out) :: artv !! \(\frac{d^2Ar}{dTdV}\)
      real(pr), optional, intent(out) :: arv2 !! \(\frac{d^2Ar}{dV^2}\)
      real(pr), optional, intent(out) :: ArT2 !! \(\frac{d^2Ar}{dT^2}\)
      real(pr), optional, intent(out) :: Arn(size(n)) !! \(\frac{dAr}{dn_i}\)
      real(pr), optional, intent(out) :: ArVn(size(n)) !! \(\frac{d^2Ar}{dVdn_i}\)
      real(pr), optional, intent(out) :: ArTn(size(n)) !! \(\frac{d^2Ar}{dTdn_i}\)
      real(pr), optional, intent(out) :: Arn2(size(n), size(n)) !! \(\frac{d^2Ar}{dn_{ij}}\)

      real(pr) :: Bmix, dBi(size(n)), dBij(size(n), size(n))
      real(pr) :: D, dDi(size(n)), dDij(size(n), size(n)), dDidT(size(n)), dDdT, dDdT2
      real(pr) :: totn
      ! estas ultimas puede por el momento no servirian
      real(pr) d1, dD1i(size(n)), dD1ij(size(n), size(n))
      real(pr) :: auxD2, fD1, fBD1, fVD1, fD1D1
      real(pr) d2
      
      real(pr) :: alpha_ads_mix, dalpha_ads_mixi(size(n)), dalpha_ads_mixij(size(n), size(n))

      real(pr) :: f_mod, g_mod, f_modv, g_modv, f_modv2, g_modv2

      real (pr), dimension(size(n)) :: f_modn, g_modn, g_modvn, f_modvn, f_modvn_AUX_1, f_modvn_AUX_2
      real (pr), dimension(size(n), size(n)) :: g_modn2, f_modn2, g_modn2_AUX, f_modn2_AUX_1, f_modn2_AUX_2

      real(pr) :: Tr(size(n)), a(size(n)), dadt(size(n)), dadt2(size(n))


      integer :: i, j, nc
      
      ! !! Prueba deri numerica
      ! real(pr) :: X(size(n)+2), del, F
      ! real(pr) :: X_back_back(size(X)), X_back_ford(size(X)), X_ford_back(size(X)), X_ford_ford(size(X)) 
      ! real(pr), dimension(2*size(n)) :: FdX_back_back, FdX_ford_ford, FdX_back_ford, FdX_ford_back
      ! real(pr) :: FdX(size(X)), dx(size(X))
      ! real(pr) :: FdX2(size(X)), X_back(size(X)), X_ford(size(X))

      ! real(pr) :: dF_mat(size(n),size(n))
      ! real(pr) :: FdX_mat_back_back(size(n),size(n)), FdX_mat_back_ford(size(n),size(n)),&
      ! FdX_mat_ford_back(size(n),size(n)), FdX_mat_ford_ford(size(n),size(n)), numdiff_vec_cruz(size(n)*2)
      
      ! X=(/n,v,t/)
      ! del=1.e-3_pr






      nc = size(n)
      TOTN = sum(n)

      Tr = T/self%components%Tc
     
      ! ========================================================================
      ! Attractive parameter and derivatives
      ! ------------------------------------------------------------------------
      call self%alpha%alpha(Tr, a, dadt, dadt2)
      a = self%ac * a
      dadt = self%ac * dadt / self%components%Tc
      dadt2 = self%ac * dadt2 / self%components%Tc**2
      
      ! ========================================================================
      ! Mixing rules
      ! ------------------------------------------------------------------------
      call self%mixrule%D1mix(n, self%del1, D1, dD1i, dD1ij)
      call self%mixrule%Bmix(n, self%b, Bmix, dBi, dBij)
      call self%mixrule%Dmix(&
         n, T, a, dadt, dadt2, D, dDdT, dDdT2, dDi, dDidT, dDij&
         )
      D2 = (1._pr - D1)/(1._pr + D1)
      !call self%mixrule%alpha_ads_mix_linear(n ,self%alpha_ads_i, alpha_ads_mix, dalpha_ads_mixi, dalpha_ads_mixij)
      
      select type(mixrule => self%mixrule)
         class is (CubicMixRuleNano)
        ! Ahora sabemos que mixrule es de tipo CubicMixRuleNano
            call mixrule%alpha_ads_mix_linear(n, self%alpha_ads_i, alpha_ads_mix, dalpha_ads_mixi, dalpha_ads_mixij)
         class default
        ! Opcional: Manejo de casos donde mixrule no es de tipo CubicMixRuleNano
            error stop "mixrule no es del tipo CubicMixRuleNano, no se puede llamar alpha_ads_mix_linear."
      end select

      ! ========================================================================
      ! Main functions defined by Møllerup and modified 
      ! The f's and g's used here are for Ar, not F (reduced Ar)
      ! This requires to multiply by R all g, f
      ! ------------------------------------------------------------------------
      f_mod = log((alpha_ads_mix*V + D1*Bmix)/(alpha_ads_mix*V + D2*Bmix))/(Bmix*alpha_ads_mix*(D1 - D2))
      g_mod = R*(((1/alpha_ads_mix)-1)*log(V)+(log(alpha_ads_mix*V-Bmix)-log(alpha_ads_mix*V))/alpha_ads_mix)
      g_modv = R*((1/(V*alpha_ads_mix-Bmix))-1/V)
      f_modv = -1/((V*alpha_ads_mix + D1*Bmix)*(V*alpha_ads_mix + D2*Bmix))
      g_modv2 = R*(1/V**2 - alpha_ads_mix/(V*alpha_ads_mix - Bmix)**2)
      f_modv2 = alpha_ads_mix*(Bmix*D1+Bmix*D2+(2*V*alpha_ads_mix))/&
      (((Bmix*D1+V*alpha_ads_mix)**2)*((Bmix*D2+V*alpha_ads_mix)**2))

      g_modn = R/(alpha_ads_mix**2)*(((dBi*alpha_ads_mix-Bmix*dalpha_ads_mixi)/(Bmix-V*alpha_ads_mix))-dalpha_ads_mixi*&
      log(V)+dalpha_ads_mixi*(log(alpha_ads_mix*V)-log((alpha_ads_mix*V-Bmix))))
      f_modvn_AUX_1 = -(dalpha_ads_mixi/alpha_ads_mix)+dBi/Bmix
      f_modvn_AUX_2 = dalpha_ads_mixi/alpha_ads_mix+dBi/Bmix
      f_modn = -f_modv*(V*(f_modvn_AUX_1))-f_mod*(f_modvn_AUX_2)
      g_modvn = R*((dBi-V*dalpha_ads_mixi)/((Bmix-V*alpha_ads_mix)**2))
      f_modvn = -(f_modvn_AUX_1*(f_modv2*V+f_modv))-f_modvn_AUX_2*f_modv

      if (present(Arn2)) then
         do i = 1, nc
            do j = 1, i
               g_modn2_AUX(i,j) = (dBi(i)*dalpha_ads_mixi(j)+dBij(i,j)*alpha_ads_mix-&
               Bmix*dalpha_ads_mixij(i,j)-dBi(j)*dalpha_ads_mixi(i))/(Bmix-V*alpha_ads_mix)
               
               g_modn2(i,j) = ((-2*dalpha_ads_mixi(j)*g_modn(i))/alpha_ads_mix)+(R/(alpha_ads_mix**2))*&
               (dalpha_ads_mixij(i,j)*(log(V*alpha_ads_mix)-log(V*alpha_ads_mix-Bmix))-&
               dalpha_ads_mixi(i)*((Bmix-V*alpha_ads_mix)*g_modvn(j)/R-dalpha_ads_mixi(j)/alpha_ads_mix)-&
               dalpha_ads_mixij(i,j)*log(V)+g_modn2_AUX(i,j)-(dBi(i)*alpha_ads_mix-Bmix*dalpha_ads_mixi(i))*g_modvn(j)/R)
               
               f_modn2_AUX_1(i,j) = dBij(i,j)/Bmix-dalpha_ads_mixij(i,j)/alpha_ads_mix+&
               (dalpha_ads_mixi(i)*dalpha_ads_mixi(j))/((alpha_ads_mix**2))-(dBi(i)*dBi(j))/((Bmix**2))

               f_modn2_AUX_2(i,j) = dBij(i,j)/Bmix+dalpha_ads_mixij(i,j)/alpha_ads_mix-&
               (dalpha_ads_mixi(i)*dalpha_ads_mixi(j))/((alpha_ads_mix**2))-(dBi(i)*dBi(j))/((Bmix**2))
               
               f_modn2(i,j) = -V*(f_modvn(j)*f_modvn_AUX_1(i)+f_modv*f_modn2_AUX_1(i,j))-&
               (f_modn(j)*f_modvn_AUX_2(i)+f_mod*f_modn2_AUX_2(i,j))
               
               g_modn2_AUX(j,i) = g_modn2_AUX(i,j)
               g_modn2(j,i) = g_modn2(i,j)
               f_modn2_AUX_1(j,i) = f_modn2_AUX_1(i,j)
               f_modn2_AUX_2(j,i) = f_modn2_AUX_2(i,j)
               f_modn2(j,i) = f_modn2(i,j)
            end do
         end do
      end if

      ! ========================================================================
      ! Reduced Helmholtz Energy and derivatives
      ! ------------------------------------------------------------------------
      if (present(Ar)) Ar = -TOTN*g_mod*T - D*f_mod
      if (present(ArV)) ArV = -TOTN*g_modv*T - D*f_modv
      if (present(ArV2)) ArV2 = -TOTN*g_modv2*T - D*f_modv2

      if (present(Arn))  Arn(:)  = -T*(g_mod + TOTN*g_modn(:)) - dDi(:)*f_mod - D*f_modn(:)
      if (present(ArVn)) ArVn(:) = -T*(g_modv + TOTN*g_modvn(:)) - dDi(:)*f_modv - D*f_modvn(:)
      if (present(ArTn)) ArTn(:) = -(g_mod + TOTN*g_modn(:)) - f_mod*dDidT(:) - f_modn(:)*dDdT

      if (present(Arn2)) then
         do i = 1, nc
            do j = 1, i
               Arn2(i, j) = -T*(g_modn(j)+g_modn(i)+TOTN*g_modn2(i,j))-dDij(i,j)*f_mod-&
               dDi(i)*f_modn(j)-dDi(j)*f_modn(i)-D*f_modn2(i,j)

               Arn2(j, i) = Arn2(i, j)
            end do
         end do
      end if

      ! !! Prueba der numerica Arvn y Artn

      ! if (present(ArVn)) then
      !    ! print*, "aaaaaaaaaaaaaaaaa", n, v, t
      !    ! print*, "aaaaaaaaaaaaaaaasda", X
      !    do i=1,nc
      !       dx = 0
      !       dx(i) = del * X(i)
      !       dx(nc+1) = del * X(nc+1)
      !       X_back_back = 0
      !       X_back_back = X-dx
      !       X_ford_ford = 0
      !       X_ford_ford = X+dx
      !       X_back_ford = X
      !       X_back_ford(i) = X(i)-dx(i)
      !       X_back_ford(nc+1) = X_back_ford(nc+1)+dx(nc+1)
      !       X_ford_back = X
      !       X_ford_back(i) = X(i)+dx(i)
      !       X_ford_back(nc+1) = X_ford_back(nc+1)-dx(nc+1)
            
      !       call self%residual_helmholtz(X_back_back(:nc), X_back_back(nc+1),&
      !        X_back_back(nc+2), Ar=FdX_back_back(i))
      !       call self%residual_helmholtz(X_ford_ford(:nc), X_ford_ford(nc+1),&
      !        X_ford_ford(nc+2), Ar=FdX_ford_ford(i))
      !       call self%residual_helmholtz(X_back_ford(:nc), X_back_ford(nc+1),&
      !        X_back_ford(nc+2), Ar=FdX_back_ford(i))
      !       call self%residual_helmholtz(X_ford_back(:nc), X_ford_back(nc+1),&
      !        X_ford_back(nc+2), Ar=FdX_ford_back(i))         
      !       ArVn(i)=(FdX_ford_ford(i)-FdX_ford_back(i)-&
      !       FdX_back_ford(i)+FdX_back_back(i))/(4*dx(i)*dx(nc+1))
      !    end do
      !    !ArVn(:) = numdiff_vec_cruz(:nc)
      ! end if
      ! if (present(ArTn)) then
      !    do i=1,nc
      !       dx = 0
      !       dx(i) = del * X(i)
      !       dx(nc+2) = del * X(nc+2)
      !       X_back_back = 0
      !       X_back_back = X-dx
      !       X_ford_ford = 0
      !       X_ford_ford = X+dx
      !       X_back_ford = X
      !       X_back_ford(i) = X(i)-dx(i)
      !       X_back_ford(nc+2) = X_back_ford(nc+2)+dx(nc+2)
      !       X_ford_back = X
      !       X_ford_back(i) = X(i)+dx(i)
      !       X_ford_back(nc+2) = X_ford_back(nc+2)-dx(nc+2)
      !       call self%residual_helmholtz(X_back_back(:nc), X_back_back(nc+1),&
      !        X_back_back(nc+2), Ar=FdX_back_back(nc+i))
      !       call self%residual_helmholtz(X_ford_ford(:nc), X_ford_ford(nc+1),&
      !        X_ford_ford(nc+2), Ar=FdX_ford_ford(nc+i))
      !       call self%residual_helmholtz(X_back_ford(:nc), X_back_ford(nc+1),&
      !        X_back_ford(nc+2), Ar=FdX_back_ford(nc+i))
      !       call self%residual_helmholtz(X_ford_back(:nc), X_ford_back(nc+1),&
      !        X_ford_back(nc+2), Ar=FdX_ford_back(nc+i))         
            
      !       ArTn(i)=(FdX_ford_ford(nc+i)-FdX_ford_back(nc+i)-&
      !       FdX_back_ford(nc+i)+FdX_back_back(nc+i))/(4*dx(i)*dx(nc+2))
      !    end do
      !    !ArTn(:) = numdiff_vec_cruz(nc+1:)
      ! end if


      ! !! prueba der numerica Arn2
      ! if (present(Arn2)) then
      !    call self%residual_helmholtz(X(:nc), X(nc+1), X(nc+2), Ar=F)

      !    do i = 1, nc
      !       do j = 1, nc
      !          if (i==j) then
      !             dx=0
      !             dx(i) = del * X(i)
      !             X_back = 0
      !             X_back = X-dx
      !             X_ford = 0
      !             X_ford = X+dx
      !             call self%residual_helmholtz(X_back(:nc), X_back(nc+1), &
      !             X_back(nc+2), Ar=FdX(i))
      !             call self%residual_helmholtz(X_ford(:nc), X_ford(nc+1), &
      !             X_ford(nc+2), Ar=FdX2(i))
      !             Arn2(i,j)=(FdX2(i) - 2*F + FdX(i))/((dx(i))**2)
      !          else
      !             dx=0
      !             dx(i) = del * X(i)
      !             dx(j) = del * X(j)
      !             X_back_back = 0
      !             X_back_back = X-dx
      !             X_ford_ford = 0
      !             X_ford_ford = X+dx
      !             X_back_ford = X
      !             X_back_ford(i) = X(i)-dx(i)
      !             X_back_ford(j) = X_back_ford(j)+dx(j)
      !             X_ford_back = X
      !             X_ford_back(i) = X(i)+dx(i)
      !             X_ford_back(j) = X_ford_back(j)-dx(j)
      !             call self%residual_helmholtz(X_back_back(:nc), &
      !             X_back_back(nc+1), X_back_back(nc+2), Ar=FdX_mat_back_back(i,j))
      !             call self%residual_helmholtz(X_ford_ford(:nc), &
      !             X_ford_ford(nc+1), X_ford_ford(nc+2), Ar=FdX_mat_ford_ford(i,j))
      !             call self%residual_helmholtz(X_back_ford(:nc), &
      !             X_back_ford(nc+1), X_back_ford(nc+2), Ar=FdX_mat_back_ford(i,j))
      !             call self%residual_helmholtz(X_ford_back(:nc), &
      !             X_ford_back(nc+1), X_ford_back(nc+2), Ar=FdX_mat_ford_back(i,j)) 
      !             Arn2(i,j) = (FdX_mat_ford_ford(i,j)-FdX_mat_ford_back(i,j)-&
      !             FdX_mat_back_ford(i,j)+FdX_mat_back_back(i,j))/(4*dx(i)*dx(j))
      !          end if
      !       end do
      !    end do
      ! end if

      ! TEMPERATURE DERIVATIVES
      if (present(ArT))  ArT = -TOTN*g_mod - dDdT*f_mod
      if (present(ArTV)) ArTV = -TOTN*g_modv - dDdT*f_modv
      if (present(ArT2)) ArT2 = -dDdT2*f_mod



   end subroutine GenericCubic_Ar_Nano

   function v0_wrapper(self, n, p, t)
      class(CubicEoSNano), intent(in) :: self
      real(pr), intent(in) :: n(:), p, t
      real(pr) :: v0_wrapper

      !real(pr) :: dbi(size(n)), dbij(size(n), size(n))
      !call self%mixrule%Bmix(n, self%b, v0, dbi, dbij)
      v0_wrapper = v0(self, n, p, t)
   end function

   subroutine D1mix_constant_wrapper(self, n, d1i, D1, dD1i, dD1ij)
      class(CubicMixRuleNano), intent(in) :: self !! Mixing rule
      real(pr), intent(in) :: n(:) !! Moles vector
      real(pr), intent(in) :: d1i(:) !! \(\delta_1\) parameter
      real(pr), intent(out) :: D1 !! Mixture's \(\Delta_1\)
      real(pr), intent(out) :: dD1i(:) !! \(\frac{dDelta_1}{dn_i} = 0\)
      real(pr), intent(out) :: dD1ij(:, :) !! \(\frac{d^2Delta_1}{dn_{ij}} = 0\)
      call D1mix_constant(self, n, d1i, D1, dD1i, dD1ij)
   end subroutine

   subroutine Bmix_wrapper(self, n, bi, B, dBi, dBij)
      class(CubicMixRuleNano), intent(in) :: self !! Mixing rule
      real(pr), intent(in) :: n(:) !! Moles vector.
      real(pr), intent(in) :: bi(:) !! Pure components repulsive parameters.
      real(pr), intent(out) :: B !! Mixture repulsive parameter.
      real(pr), intent(out) :: dBi(:) !! \(\frac{dB}{dn_i}\)
      real(pr), intent(out) :: dBij(:, :) !!\(\frac{d^2B}{dn_{ij}}\)
      call Bmix(self, n, bi, B, dBi, dBij)
   end subroutine

   subroutine Dmix_wrapper(self, n, T, &
      ai, daidt, daidt2, &
      D, dDdT, dDdT2, dDi, dDidT, dDij)
      class(CubicMixRuleNano), intent(in) :: self !! Mixing rule
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: n(:) !! Moles vector [mol]
      real(pr), intent(in) :: ai(:) !! Pure components attractive parameters \(a_i\)
      real(pr), intent(in) :: daidt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: daidt2(:) !! \(\frac{d^2a_i}{dT^2}\)

      real(pr), intent(out) :: D !! Mixture attractive parameter \(n^2a_{mix}\)
      real(pr), intent(out) :: dDdT !! \(\frac{dD}{dT}\)
      real(pr), intent(out) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)
      real(pr), intent(out) :: dDi(:) !! \(\frac{dD}{dn_i}\)
      real(pr), intent(out) :: dDidT(:) !! \(\frac{d^2D}{dTn_i}\)
      real(pr), intent(out) :: dDij(:, :)!! \(\frac{d^2D}{dn_{ij}}\)
      
      ! Variable temporal para manejar la conversión
      !class(QMR), pointer :: base_ptr

      ! Apuntar a la parte base del objeto
      !base_ptr => self

      ! Llamar a Dmix usando la referencia al tipo base
      call Dmix(self, n, T, &
      ai, daidt, daidt2, &
      D, dDdT, dDdT2, dDi, dDidT, dDij)
   end subroutine

end module
      !--------------------------- aca lo viejo ------------------------------------------
      ! real(pr) :: Bmix, dBi(size(n)), dBij(size(n), size(n))
      ! real(pr) :: D, dDi(size(n)), dDij(size(n), size(n)), dDidT(size(n)), dDdT, dDdT2

      ! real(pr) :: totn
      ! real(pr) d1, dD1i(size(n)), dD1ij(size(n), size(n))
      ! real(pr) :: auxD2, fD1, fBD1, fVD1, fD1D1
      ! real(pr) d2

      ! real(pr) :: f, g, fv, fB, gv, fv2, gv2, AUX, FFB, FFBV, FFBB

      ! real(pr) :: Tr(size(n)), a(size(n)), dadt(size(n)), dadt2(size(n))


      ! integer :: i, j, nc

      ! nc = size(n)
      ! TOTN = sum(n)

     
      ! Tr = T/self%components%Tc
     
      ! ! ========================================================================
      ! ! Attractive parameter and derivatives
      ! ! ------------------------------------------------------------------------
      ! call self%alpha%alpha(Tr, a, dadt, dadt2)
      ! a = self%ac * a
      ! dadt = self%ac * dadt / self%components%Tc
      ! dadt2 = self%ac * dadt2 / self%components%Tc**2
      
      ! ! ========================================================================
      ! ! Mixing rules
      ! ! ------------------------------------------------------------------------
      ! call self%mixrule%D1mix(n, self%del1, D1, dD1i, dD1ij)
      ! call self%mixrule%Bmix(n, self%b, Bmix, dBi, dBij)
      ! call self%mixrule%Dmix(&
      !    n, T, a, dadt, dadt2, D, dDdT, dDdT2, dDi, dDidT, dDij&
      !    )
      ! D2 = (1._pr - D1)/(1._pr + D1)

      ! ========================================================================
      ! Main functions defined by Møllerup
      ! The f's and g's used here are for Ar, not F (reduced Ar)
      ! This requires to multiply by R all g, f
      ! ------------------------------------------------------------------------
      ! f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
      ! g = R*log(1 - Bmix/V)
      ! fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
      ! fB = -(f + V*fv)/Bmix
      ! gv = R*Bmix/(V*(V - Bmix))
      ! fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
      ! gv2 = R*(1/V**2 - 1/(V - Bmix)**2)

      ! ! DERIVATIVES OF f WITH RESPECT TO DELTA1
      ! auxD2 = (1 + 2/(1 + D1)**2)
      ! fD1 = (1/(V + D1*Bmix) + 2/(V + D2*Bmix)/(1 + D1)**2) - f*auxD2
      ! fD1 = fD1/(D1 - D2)
      ! fBD1 = -(fB*auxD2 + D1/(V + D1*Bmix)**2 + 2*D2/(V + D2*Bmix)**2/(1 + D1)**2)
      ! fBD1 = fBD1/(D1 - D2)
      ! fVD1 = -(fV*auxD2 + 1/(V + D1*Bmix)**2 + 2/(V + D2*Bmix)**2/(1 + D1)**2)/(D1 - D2)
      ! fD1D1 = 4*(f - 1/(V + D2*Bmix))/(1 + D1)**3 + Bmix*(-1/(V + D1*Bmix)**2 &
      !       + 4/(V + D2*Bmix)**2/(1 + D1)**4) - 2*fD1*(1 + 2/(1 + D1)**2)
      !       fD1D1 = fD1D1/(D1 - D2)

      ! AUX = R*T/(V - Bmix)
      ! FFB = TOTN*AUX - D*fB
      ! FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
      ! FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2

      ! ========================================================================
      ! Reduced Helmholtz Energy and derivatives
      ! ------------------------------------------------------------------------
      ! if (present(Ar)) Ar = -TOTN*g*T - D*f
      ! if (present(ArV)) ArV = -TOTN*gv*T - D*fv
      ! if (present(ArV2)) ArV2 = -TOTN*gv2*T - D*fv2

      ! if (present(Arn))  Arn(:)  = -g*T + FFB*dBi(:) - f*dDi(:) - D*fD1 * dD1i(:)
      ! if (present(ArVn)) ArVn(:) = -gv*T + FFBV*dBi(:) - fv*dDi(:) - D*fVD1*dD1i(:)
      ! if (present(ArTn)) ArTn(:) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(:) - f*dDidT(:) - dDdT*fD1*dD1i(:)

      ! if (present(Arn2)) then
      !    do i = 1, nc
      !       do j = 1, i
      !          Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
      !             + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
      !          Arn2(i, j) = Arn2(i, j) - D*fBD1*(dBi(i)*dD1i(j) + dBi(j)*dD1i(i)) &
      !                   - fD1*(dDi(i)*dD1i(j) + dDi(j)*dD1i(i)) &
      !                   - D*fD1*dD1ij(i, j) - D*fD1D1*dD1i(i)*dD1i(j)
      !          Arn2(j, i) = Arn2(i, j)
      !       end do
      !    end do
      ! end if

      ! ! TEMPERATURE DERIVATIVES
      ! if (present(ArT))  ArT = -TOTN*g - dDdT*f
      ! if (present(ArTV)) ArTV = -TOTN*gv - dDdT*fV
      ! if (present(ArT2)) ArT2 = -dDdT2*f
