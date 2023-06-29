program main
    use constants, only: pr
    use pengrobinson76, only: setup_pr76
    use cubic_eos, only: ac, b, wmod => w, k
    use ar_models, only: residual_helmholtz
    use yaeos_thermo_properties, only: get_volume, pressure
    implicit none

    integer, parameter :: n = 2
    integer :: i

    real(pr) :: tc(n), pc(n), w(n), z(n), kij(n, n), lij(n, n)
    real(pr) :: t, v, v_vap, p , v_liq, vf=1500, dv

    z = [0.3, 0.7]
    tc = [190, 310]
    pc = [14, 30]
    w = [0.001, 0.03]

    kij = 0
    lij = 0

    call setup_pr76(n, tc, pc, w, kij, lij)

    t = 150._pr

    i = 0
    v = 7.0e-2
    vf = 3
    
    dv = 0.01
    do while (v <= vf)
        if (v < 5) then
            dv = vf/100000
        else
            dv = 0.5
        end if
        
        v = v + i * dv
        call pressure(z, v, t, p)
        call get_volume(z, p, t, v_vap, "vapor")
        call get_volume(z, p, t, v_liq, "liquid")
        print *, v, p, v_vap, v_liq
        i = i + 1
    end do
end program main

! module test_data
!     use constants, only: pr
!     
!     real(pr), parameter :: t_srk(100)= [
!    1246.3400827135777       
!    622.75421732772509   
!    414.89193604900760   
!    310.96055043246241   
!    248.60152278853064   
!    207.02867388688671   
!    177.33364119984182   
!    155.06224346319507   
!    137.73990511402485   
!    123.88195267413609   
!    112.54353583182481   
!    103.09477036185250   
!    95.099582567883672   
!    88.246491236432846   
!    82.307076809232512   
!    77.110024615962686   
!    72.524329318958635   
!    68.448097876348697   
!    64.800888411569744   
!    61.518347422027908   
!    58.548379213855860   
!    55.848360036281861   
!    53.383078963554446   
!    51.123193556385850   
!    49.044056164648886   
!    47.124811083297473   
!    45.347692340304803   
!    43.697471958614479   
!    42.161022370791478   
!    40.726966350742693   
!    39.385394700588307   
!    38.127636871244626   
!   0.10192619015766519   
!   0.10192570043676165   
!   0.10192521074734594   
!   0.10192472108941439   
!   0.10192423146296334   
!   0.10192374186798907   
!   0.10192325230448797   
!   0.10192276277245634   
!   0.10192227327189052   
!   0.10192178380278685   
!   0.10192129436514162   
!   0.10192080495895120   
!   0.10192031558421193   
!   0.10191982624092012   
!   0.10191933692907212   
!   0.10191884764866424   
!   0.10191835839969286   
!   0.10191786918215429   
!   0.10191737999604485   
!   0.10191689084136092   
!   0.10191640171809881   
!   0.10191591262625485   
!   0.10191542356582541   
!   0.10191493453680681   
!   0.10191444553919538   
!   0.10191395657298749   
!   0.10191346763817948   
!   0.10191297873476766   
!   0.10191248986274840   
!   0.10191200102211806   
!   0.10191151221287295   
!   0.10191102343500943   
!   0.10191053468852386   
!   0.10191004597341255   
!   0.10190955728967191   
!   0.10190906863729821   
!   0.10190858001628784   
!   0.10190809142663716   
!   0.10190760286834250   
!   0.10190711434140022   
!   0.10190662584580666   
!   0.10190613738155820   
!   0.10190564894865115   
!   0.10190516054708190   
!   0.10190467217684678   
!   0.10190418383794216   
!   0.10190369553036438   
!   0.10190320725410981   
!   0.10190271900917482   
!   0.10190223079555571   
!   0.10190174261324891   
!   0.10190125446225071   
!   0.10190076634255751   
!   0.10190027825416567   
!   0.10189979019707156   
!   0.10189930217127149   
!   0.10189881417676189   
!   0.10189832621353906   
!   0.10189783828159939   
!   0.10189735038093922   
!   0.10189686251155497   
!   0.10189637467344295   
!   0.10189588686659957   
!   0.10189539909102115   
!   0.10189491134670409   
!   0.10189442363364472   
!   0.10189393595183946   
!   0.10189344830128465   
! 
!     ]
!  
!     1.0000000000000000E-002       
!       2.0000000000000000E-002
!     2.9999999999999999E-002
!       4.0000000000000001E-002
!       5.0000000000000003E-002
!       5.9999999999999998E-002
!       7.0000000000000007E-002
!       8.0000000000000002E-002
!       8.9999999999999997E-002
!      0.10000000000000001     
!      0.11000000000000000     
!      0.12000000000000000     
!      0.13000000000000000     
!      0.14000000000000001     
!      0.14999999999999999     
!      0.16000000000000000     
!      0.17000000000000001     
!      0.17999999999999999     
!      0.19000000000000000     
!      0.20000000000000001     
!      0.20999999999999999     
!      0.22000000000000000     
!      0.23000000000000001     
!      0.23999999999999999     
!      0.25000000000000000     
!      0.26000000000000001     
!      0.27000000000000002     
!      0.28000000000000003     
!      0.28999999999999998     
!      0.29999999999999999     
!      0.31000000000000000     
!      0.32000000000000001     
!      0.33000000000000002     
!      0.34000000000000002     
!      0.34999999999999998     
!      0.35999999999999999     
!      0.37000000000000000     
!      0.38000000000000000     
!      0.39000000000000001     
!      0.40000000000000002     
!      0.40999999999999998     
!      0.41999999999999998     
!      0.42999999999999999     
!      0.44000000000000000     
!      0.45000000000000001     
!      0.46000000000000002     
!      0.46999999999999997     
!      0.47999999999999998     
!      0.48999999999999999     
!      0.50000000000000000     
!      0.51000000000000001     
!      0.52000000000000002     
!      0.53000000000000003     
!      0.54000000000000004     
!      0.55000000000000004     
!      0.56000000000000005     
!      0.56999999999999995     
!      0.57999999999999996     
!      0.58999999999999997     
!      0.59999999999999998     
!      0.60999999999999999     
!      0.62000000000000000     
!      0.63000000000000000     
!      0.64000000000000001     
!      0.65000000000000002     
!      0.66000000000000003     
!      0.67000000000000004     
!      0.68000000000000005     
!      0.68999999999999995     
!      0.69999999999999996     
!      0.70999999999999996     
!      0.71999999999999997     
!      0.72999999999999998     
!      0.73999999999999999     
!      0.75000000000000000     
!      0.76000000000000001     
!      0.77000000000000002     
!      0.78000000000000003     
!      0.79000000000000004     
!      0.80000000000000004     
!      0.81000000000000005     
!      0.81999999999999995     
!      0.82999999999999996     
!      0.83999999999999997     
!      0.84999999999999998     
!      0.85999999999999999     
!      0.87000000000000000     
!      0.88000000000000000     
!      0.89000000000000001     
!      0.90000000000000002     
!      0.91000000000000003     
!      0.92000000000000004     
!      0.93000000000000005     
!      0.93999999999999995     
!      0.94999999999999996     
!      0.95999999999999996     
!      0.96999999999999997     
!      0.97999999999999998     
!      0.98999999999999999     
!       1.0000000000000000 