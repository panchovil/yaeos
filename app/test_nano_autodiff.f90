


program test_nano_autodiff
   use hyperdual_pr78_nano
   use yaeos, only: pr, R, Substances, AlphaSoave, CubicEoSNano, CubicEoS, &
   GenericCubic_Ar_Nano, GenericCubic_Ar
   use yaeos, only: ArModel, PengRobinson78Nano
   implicit none
   class(PR78_nano_autodiff), allocatable :: auto_eos
   class(CubicEoSNano), allocatable :: eos
   
   integer, parameter :: n=4
   real(pr) :: z(n)
   real(pr) :: v, t, p
   ! type(hyperdual) :: v, t, p , z(n), ar_result
   real(pr) :: ar
   real(pr) :: art, arv, arv2, art2, artv
   real(pr) :: arn(size(z)), arvn(size(z)), artn(size(z)), arn2(size(z),size(z)) 
   real(pr) :: auto_ar
   real(pr) :: auto_art, auto_arv, auto_arv2, auto_art2, auto_artv
   real(pr) :: auto_arn(size(z)), auto_arvn(size(z)), auto_artn(size(z)), auto_arn2(size(z),size(z)) 


   real(pr) :: tc(size(z)), pc(size(z)), w(size(z)), kij(size(z), size(z)), lij(size(z), size(z)), LJ_par(size(z)), rp
   integer :: i


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
   eos = PengRobinson78Nano(LJ_par=LJ_par, rp=rp, tc=tc, pc=pc, w=w, kij=kij, lij=lij)


   z = [0.5_pr, 0.2_pr, 0.2_pr, 0.1_pr] 
   t = 500._pr
   v = 400.0_pr
   !ar_result = auto_eos%Ar(z, v, t)
   !arfun(auto_eos, z, v, t)
   !print*, ar_result
   Ar=0
   ArV=0
   arv2=0
   Art=0
   art2=0
   artv=0
   arn=0
   artn=0
   arvn=0
   arn2=0
   
   auto_ar=0
   auto_arv=0
   auto_arv2=0
   auto_art=0
   auto_art2=0
   auto_artv=0
   auto_arn=0
   auto_artn=0
   auto_arvn=0
   auto_arn2=0

   call auto_eos%residual_helmholtz(n=z, v=v, t=t, Ar=auto_ar, ArV=auto_arv, ArT=auto_art,&
   ArTV=auto_artv, ArV2=auto_arv2, art2=auto_art2, Arn=auto_Arn, ArVn=auto_arvn, &
   ArTn=auto_artn, Arn2=auto_arn2)
   call eos%residual_helmholtz(&
   n=z, v=v, t=t, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
   ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
   )



   print *, "Ar: ", ar
   print *, "auto_Ar: ", auto_ar
   print *, "Ar diff", abs(ar-auto_ar)
   print *, "-----------------------------------------------------------"


   print *, "ArV: ", arV
   print *, "auto_ArV: ", auto_arv
   print *, "ArV diff", abs(arv-auto_arv)
   print *, "-----------------------------------------------------------"

   print *, "ArT: ", arT
   print *, "auto_ArT: ", auto_art
   print *, "ArT diff", abs(art-auto_art)
   print *, "-----------------------------------------------------------"

   print *, "ArT2: ", arT2
   print *, "auto_ArT2: ", auto_art2
   print *, "ArT2 diff", abs(art2-auto_art2)
   print *, "-----------------------------------------------------------"

   print *, "ArV2: ", ArV2
   print *, "auto_ArV2: ", auto_arv2
   print *, "ArV2 diff", abs(arv2-auto_arv2)
   print *, "-----------------------------------------------------------"
   
   print *, "ArTV: ", ArTV
   print *, "auto_ArTV: ", auto_arTv
   print *, "ArTV diff", abs(artv-auto_artv)
   print *, "-----------------------------------------------------------"
   
   print *, "Arn: ", Arn
   print *, "auto_Arn: ", auto_arn
   print *, "Arn diff", abs(arn-auto_arn)
   print *, "-----------------------------------------------------------"

   print *, "ArVn: ", ArVn
   print *, "auto_ArVn: ", auto_arvn
   print *, "ArVn diff", abs(arvn-auto_arvn)
   print *, "-----------------------------------------------------------"

   print *, "ArTn: ", ArTn
   print *, "auto_ArTn: ", auto_arTn
   print *, "ArTn diff", abs(artn-auto_artn)
   print *, "-----------------------------------------------------------"

   print *, "Arn2: "
   do i=1,n
      print*, Arn2(i,:)
   end do
   print *, "auto_Arn2: "
   do i=1,n
      print*, auto_arn2(i,:)
   end do
   print *, "Arn2 diff "
   do i=1,n
      print*, abs(Arn2(i,:)-auto_arn2(i,:))
   end do
   
   

end program test_nano_autodiff