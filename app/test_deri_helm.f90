program test_deri
   !use yaeos, only: AlphaSoave, CubicEoSNano, GenericCubic_Ar, ArModel
   ! use yaeos__constants, only: pr, R
   ! use yaeos__models_ar_nanocubic, only: CubicEoSNano
   ! use yaeos__models_ar_cubic_implementations, only: PengRobinson78Nano
   !use yaeos 
   use yaeos, only: pr, R, Substances, AlphaSoave, CubicEoSNano, CubicEoS, GenericCubic_Ar_Nano, GenericCubic_Ar
   use yaeos, only: ArModel, PengRobinson78Nano
   implicit none
   class(CubicEoSNano),allocatable :: eos
   !type(CubicEoSNano), target :: eos

   integer, parameter :: n=4
   real(pr) :: z(n)
   real(pr) :: v=1.0, t=150.0, p

   real(pr) :: ar
   real(pr) :: art, arv, arv2, art2, artv
   real(pr) :: arn(size(z)), arvn(size(z)), artn(size(z)), arn2(size(z),size(z)) 
   real(pr) :: X(n+2) 
   !real(pr) :: lnfug(n), dlnphidp(n), dlnphidt(n), dlnphidn(n, n)


   real(pr) :: tc(size(z)), pc(size(z)), w(size(z)), kij(size(z), size(z)), lij(size(z), size(z)), LJ_par(size(z)), rp

   z = [0.1_pr, 0.4_pr, 0.3_pr, 0.2_pr] !CH4, CO2, C4H10(Butano), C6H14(Hexano)
   tc = [190.56_pr, 304.13_pr, 452.2_pr, 507.9_pr]
   pc = [45.99_pr, 73.8_pr, 37.96_pr, 30.25_pr]
   w = [0.011_pr, 0.225_pr, 0.199_pr, 0.297_pr]
   LJ_par = [0.3758_pr, 0.300_pr, 0.443_pr, 0.481_pr] !nm 
   rp = 100_pr !nm

   kij = reshape([0.000_pr, 0.110_pr, 0.025_pr, 0.030_pr, &
                  0.110_pr, 0.000_pr, 0.120_pr, 0.140_pr, &
                  0.025_pr, 0.120_pr, 0.000_pr, 0.015_pr, &
                  0.030_pr, 0.140_pr, 0.015_pr, 0.000_pr], [n,n]) 
   lij = kij / 2 

   !class(CubicEoSNano), allocatable :: eos
   !eos = PengRobinson78(LJ_par, rp, tc, pc, w, kij, lij)
   eos = PengRobinson78Nano(LJ_par, rp, tc, pc, w, kij, lij)

   v = 200
   t = 575

!    call eos%lnphi_vt(&
!         z, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
!    )
   ! Ar=0
   ! ArV=0
   ! arv2=0
   ! Art=0
   ! art2=0
   ! artv =0
   ! arn=0
   ! artn=0
   ! arvn=0
   ! arn2=0
   X=(/z,v,t/)
   print*, "X",X

   test_numdiff: block

         real(pr) :: F, numdiff_esc(size(X)), dF_esc(size(X)), dF_vec_2orden(2), numdiff_vec_2orden(2)
         real(pr) :: dF_vec_cruz(2*n+1), numdiff_vec_cruz(2*n+1)
         real(pr) :: X_back_back(size(X)), X_back_ford(size(X)), X_ford_back(size(X)), X_ford_ford(size(X)) 
         real(pr) :: FdX_back_back(2*n+1), FdX_ford_ford(2*n+1), FdX_back_ford(2*n+1), FdX_ford_back(2*n+1)
         real(pr) :: FdX(size(X)), dx(size(X)), dFdS(size(X))
         real(pr) :: FdX2(size(X)), X_back(size(X)), X_ford(size(X))
         real(pr) :: numdiff_mat(n,n), dF_mat(n,n)
         real(pr) :: FdX_mat_back_back(n,n), FdX_mat_back_ford(n,n),&
         FdX_mat_ford_back(n,n), FdX_mat_ford_ford(n,n)
         integer :: i, j
         integer :: loc(1)
         real(pr) :: maxerr

         do i=1,size(X)
            dx = 0
            dx(i) = 1.e-3_pr * X(i)
            X_back = 0
            X_back = X-dx
            X_ford = 0
            X_ford = X+dx
            call eos%residual_helmholtz(X_back(:n), X_back(n+1), X_back(n+2), Ar=FdX(i))
            call eos%residual_helmholtz(X_ford(:n), X_ford(n+1), X_ford(n+2), Ar=FdX2(i))


            numdiff_esc(i) = (FdX2(i) - FdX(i))/(2*dx(i))
         end do
         call eos%residual_helmholtz(&
         X(:n), X(n+1), X(n+2), Ar=F, ArV=dF_esc(n+1), ArT=dF_esc(n+2), Arn=dF_esc(:n),&
         Arv2=dF_vec_2orden(1), Art2=dF_vec_2orden(2), ArVn=dF_vec_cruz(:n),&
         ArTn=dF_vec_cruz(n+1:), ArTV=dF_vec_cruz(2*n+1), Arn2=dF_mat)
         
         numdiff_vec_2orden(1)=(FdX2(n+1) - 2*F + FdX(n+1))/((X(n+1)*1.e-3_pr)**2) !f(n,t,v+dx)
         numdiff_vec_2orden(2)=(FdX2(n+2) - 2*F + FdX(n+2))/((X(n+2)*1.e-3_pr)**2)

         do i=1,n
            dx = 0
            dx(i) = 1.e-3_pr * X(i)
            dx(n+1) = 1.e-3_pr * X(n+1)
            X_back_back = 0
            X_back_back = X-dx
            X_ford_ford = 0
            X_ford_ford = X+dx
            X_back_ford = X
            X_back_ford(i) = X(i)-dx(i)
            X_back_ford(n+1) = X_back_ford(n+1)+dx(n+1)
            X_ford_back = X
            X_ford_back(i) = X(i)+dx(i)
            X_ford_back(n+1) = X_ford_back(n+1)-dx(n+1)
            
            call eos%residual_helmholtz(X_back_back(:n), X_back_back(n+1), X_back_back(n+2), Ar=FdX_back_back(i))
            call eos%residual_helmholtz(X_ford_ford(:n), X_ford_ford(n+1), X_ford_ford(n+2), Ar=FdX_ford_ford(i))
            call eos%residual_helmholtz(X_back_ford(:n), X_back_ford(n+1), X_back_ford(n+2), Ar=FdX_back_ford(i))
            call eos%residual_helmholtz(X_ford_back(:n), X_ford_back(n+1), X_ford_back(n+2), Ar=FdX_ford_back(i))         
            numdiff_vec_cruz(i)=(FdX_ford_ford(i)-FdX_ford_back(i)-FdX_back_ford(i)+FdX_back_back(i))/(4*dx(i)*dx(n+1))
         end do
         do i=1,n
            dx = 0
            dx(i) = 1.e-3_pr * X(i)
            dx(n+2) = 1.e-3_pr * X(n+2)
            X_back_back = 0
            X_back_back = X-dx
            X_ford_ford = 0
            X_ford_ford = X+dx
            X_back_ford = X
            X_back_ford(i) = X(i)-dx(i)
            X_back_ford(n+2) = X_back_ford(n+2)+dx(n+2)
            X_ford_back = X
            X_ford_back(i) = X(i)+dx(i)
            X_ford_back(n+2) = X_ford_back(n+2)-dx(n+2)
            call eos%residual_helmholtz(X_back_back(:n), X_back_back(n+1), X_back_back(n+2), Ar=FdX_back_back(n+i))
            call eos%residual_helmholtz(X_ford_ford(:n), X_ford_ford(n+1), X_ford_ford(n+2), Ar=FdX_ford_ford(n+i))
            call eos%residual_helmholtz(X_back_ford(:n), X_back_ford(n+1), X_back_ford(n+2), Ar=FdX_back_ford(n+i))
            call eos%residual_helmholtz(X_ford_back(:n), X_ford_back(n+1), X_ford_back(n+2), Ar=FdX_ford_back(n+i))         
            numdiff_vec_cruz(n+i)=(FdX_ford_ford(n+i)-FdX_ford_back(n+i)-FdX_back_ford(n+i)+FdX_back_back(n+i))/(4*dx(i)*dx(n+2))
         end do
         dx = 0
         dx(n+1) = 1.e-3_pr * X(n+1)
         dx(n+2) = 1.e-3_pr * X(n+2)
         X_back_back = 0
         X_back_back = X-dx
         X_ford_ford = 0
         X_ford_ford = X+dx
         X_back_ford = X
         X_back_ford(n+1) = X(n+1)-dx(n+1)
         X_back_ford(n+2) = X_back_ford(n+2)+dx(n+2)
         X_ford_back = X
         X_ford_back(n+1) = X(n+1)+dx(n+1)
         X_ford_back(n+2) = X_ford_back(n+2)-dx(n+2)
         ! print*, X_back_back
         ! print*, X_ford_ford
         ! print*, X_ford_back
         ! print*, X_back_ford

         call eos%residual_helmholtz(X_back_back(:n), X_back_back(n+1), X_back_back(n+2), Ar=FdX_back_back(2*n+1))
         call eos%residual_helmholtz(X_ford_ford(:n), X_ford_ford(n+1), X_ford_ford(n+2), Ar=FdX_ford_ford(2*n+1))
         call eos%residual_helmholtz(X_back_ford(:n), X_back_ford(n+1), X_back_ford(n+2), Ar=FdX_back_ford(2*n+1))
         call eos%residual_helmholtz(X_ford_back(:n), X_ford_back(n+1), X_ford_back(n+2), Ar=FdX_ford_back(2*n+1)) 
         numdiff_vec_cruz(2*n+1)=(FdX_ford_ford(2*n+1)-FdX_ford_back(2*n+1)-FdX_back_ford(2*n+1)+&
         FdX_back_back(2*n+1))/(4*dx(n+1)*dx(n+2))


         do i=1,n
            do j=1,n
               if (i==j) then
                  dx=0
                  dx(i) = 1.e-3_pr * X(i)
                  X_back = 0
                  X_back = X-dx
                  X_ford = 0
                  X_ford = X+dx
                  call eos%residual_helmholtz(X_back(:n), X_back(n+1), X_back(n+2), Ar=FdX(i))
                  call eos%residual_helmholtz(X_ford(:n), X_ford(n+1), X_ford(n+2), Ar=FdX2(i))
                  numdiff_mat(i,j)=(FdX2(i) - 2*F + FdX(i))/((dx(i))**2)
               else
                  dx=0
                  dx(i) = 1.e-3_pr * X(i)
                  dx(j) = 1.e-3_pr * X(j)
                  X_back_back = 0
                  X_back_back = X-dx
                  X_ford_ford = 0
                  X_ford_ford = X+dx
                  X_back_ford = X
                  X_back_ford(i) = X(i)-dx(i)
                  X_back_ford(j) = X_back_ford(j)+dx(j)
                  X_ford_back = X
                  X_ford_back(i) = X(i)+dx(i)
                  X_ford_back(j) = X_ford_back(j)-dx(j)
                  call eos%residual_helmholtz(X_back_back(:n), X_back_back(n+1), X_back_back(n+2), Ar=FdX_mat_back_back(i,j))
                  call eos%residual_helmholtz(X_ford_ford(:n), X_ford_ford(n+1), X_ford_ford(n+2), Ar=FdX_mat_ford_ford(i,j))
                  call eos%residual_helmholtz(X_back_ford(:n), X_back_ford(n+1), X_back_ford(n+2), Ar=FdX_mat_back_ford(i,j))
                  call eos%residual_helmholtz(X_ford_back(:n), X_ford_back(n+1), X_ford_back(n+2), Ar=FdX_mat_ford_back(i,j)) 
                  numdiff_mat(i,j) = (FdX_mat_ford_ford(i,j)-FdX_mat_ford_back(i,j)-FdX_mat_back_ford(i,j)+&
                  FdX_mat_back_back(i,j))/(4*dx(i)*dx(j))
               end if
            end do
         end do
               


         loc = maxloc(abs(numdiff_esc - dF_esc))
         maxerr = abs((numdiff_esc(loc(1)) - dF_esc(loc(1))&
             )/numdiff_esc(loc(1)))
         if (maxerr > 0.01_pr) then
            print *, "ERROR: PXEnvel2 Numerical differentiation failed"
            loc = maxloc(abs(numdiff_esc - dF_esc))
            print *, loc
            print *, dF_esc(loc(1)), numdiff_esc(loc(1))
            ! error stop 1
         end if
         print*, "Arn, ArV, ArT analitica", dF_esc
         print*, "Arn, ArV, ArT numerica", numdiff_esc
         print*, "Arn, ArV, ArT dif", abs(dF_esc-numdiff_esc)
         print*, "Arv2, ArT2 analitica", dF_vec_2orden
         print*, "Arv2, ArT2 numerica", numdiff_vec_2orden
         print*, "Arv2, ArT2 dif", abs(dF_vec_2orden-numdiff_vec_2orden)
         print*, "ArVn, ArTn, ArTV analitica", dF_vec_cruz
         print*, "ArVn, ArTn, ArTV numerica", numdiff_vec_cruz
         print*, "ArVn, ArTn, ArTV dif", abs(dF_vec_cruz-numdiff_vec_cruz)
         print*, "Arn2 analitica"
         do i=1,n
            print*, dF_mat(i,:)
         end do
         print*, "Arn2 numerica"
         do i=1,n
            print*, numdiff_mat(i,:)
         end do
         print*, "Arn2 dif"
         do i=1,n
            print*, abs(dF_mat(i,:)-numdiff_mat(i,:))
         end do

   end block test_numdiff

         ! real(pr) :: F(size(X)), df(size(X), size(X)), numdiff(size(X), size(X))
         ! real(pr) :: FdX(size(X)), dx(size(X)), dFdS(size(X))
         ! real(pr) :: FdX2(size(X))
               ! call eos%residual_helmholtz(&
            ! X(:n), X(n+1), X(n+2), Ar=FdX(i), ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
            ! ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
            ! )
            ! call foo(X - dx, ns, S0, FdX, df, dFdS)
            ! call foo(X + dx, ns, S0, FdX2, df, dFdS)
            ! call foo(X, ns, S0, F, df, dFdS)


         ! loc = maxloc(abs(numdiff - df))
         ! maxerr = abs(&
         !    (numdiff(loc(1), loc(2)) - df(loc(1), loc(2))&
         !    )/numdiff(loc(1), loc(2)))
         ! if (maxerr > 0.01_pr) then
         !    print *, "ERROR: PXEnvel2 Numerical differentiation failed"
         !    loc = maxloc(abs(numdiff - df))
         !    print *, loc
         !    print *, df(loc(1), loc(2)), numdiff(loc(1), loc(2))
         !    ! error stop 1
         ! end if


   call eos%residual_helmholtz(&
           X(:n), X(n+1), X(n+2), Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
           ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
   )

   ! print *, "Ar: ", ar

   ! print *, "ArV: ", arV
   ! print *, "ArT: ", arT

   ! print *, "ArT2: ", arT2
   ! print *, "ArV2: ", ArV2
   
   ! print *, "ArTV: ", ArTV
   
   ! print *, "Arn: ", Arn

   ! print *, "ArVn: ", ArVn
   ! print *, "ArTn: ", ArTn

   ! print *, "Arn2: ", Arn2
   ! print *, size(z)

   !print*, "lnfug", lnfug
end program test_deri