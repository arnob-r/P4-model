! This code generates the parameter space of $gamma_0$-$xi_0$ by keeping $f$ fixed at 0.1 and calculate critical $x^*_u$.

program pppp
  implicit none

  integer :: mm, tn, r1, n1, n2, n3, i, j, k, n, i1,i2,i3
  real*8  :: t0, t, h, tt, alpha, kappa, gamma0, xi0, mA, mS, f, uh, ul, xu, ubar
  real*8  :: beta1(100,100), beta2(100,100), beta3(100,100)
  ! Adj1 = pollinators x plants (n1 x n2); Adj2 = plants x pests (n2 x n3)
  real*8  :: Adj1(100,100), Adj2(100,100)
  real*8  :: A0(100), P0(100), S0(100), A(100), P(100), S(100),A10(100), P10(100), S10(100)
  real*8  :: sum1(100), sum2(100), sum3(100), sum4(100), sum5(100), sum6(100), sum7(100)
  ! Degrees kept separate to avoid overwrite
  real*8  :: D1A(100), D1P(100), D2P(100), D2S(100)
  real*8  :: l, mu, r2, eps, den, den2, gain, gain2, hh
  ! RK storage (separate per guild to avoid overwrites)
  real*8  :: krkA(4,100), krkP(4,100), krkS(4,100)
  real*8  :: Anew(100), Pnew(100), Snew(100), number
  integer :: ISEED
  real*8  :: biiA, bijA, biiP, bijP, biiS, bijS, A00, P00, S00, flag(200), xstar

  ISEED = 1
  call srand(ISEED)

  mm = 20000
  tn = 90000
  n  = 100
  
  do r1 = 1, n
     r2 = rand()
  end do

  open(201, file='Pollinator_transition_g_xi_f_0.1.txt')                                       !!! Change

  n1 = 23   ! pollinators
  n2 = 36   ! plants
  n3 = 51   ! pests

  tt     = 0.5d0
  alpha  = 0.1d0
  kappa  = 0.0001d0     ! decay for A & S only
    
  mA     = 0.5d0  
  mS     = 0.9d0  
  
  uh     = 1.0d0  
  ul     = 0.5d0  
  
  mu     = 0.0000d0        
  hh     = 0.7d0
  
  
  f = 0.1d0          !0.2d0    ! 0.3d0                                                      !!! Change
  

 do i3 =1,30
 do i2 =1,30
    
    gamma0 = 0.1d0*i3
    xi0    = 0.1d0*i2
  
  do i1=1,100 
    
  xu = 0.01d0*i1                              ! BIFURCATION Parameter                         !!! Change
  ubar = ( (1.0d0 - xu)*uh + xu*ul )
  
!  write(*,*)i1
  
  

!  number = 1.0d0

  flag(1)=0.0D0
  
! Initial conditions (IC1: 0.2, 0.2, 0.2) Change
  do i = 1, n1
     A0(i) = 0.2d0
  end do
  do i = 1, n2
     P0(i) = 0.2d0
  end do
  do i = 1, n3
     S0(i) = 0.2d0
  end do

  biiA = 1.0d0; bijA = 0.1d0      
  biiP = 1.0d0; bijP = 0.1d0
  biiS = 1.0d0; bijS = 0.1d0


  ! Build competition matrices: diag=bii, off-diag=bij
  do i = 1, n1
     do j = 1, n1
        if (i .eq. j) then
           beta1(i,j) = biiA
        else
           beta1(i,j) = bijA
        end if
     end do
  end do
  do i = 1, n2
     do j = 1, n2
        if (i .eq. j) then
           beta2(i,j) = biiP
        else
           beta2(i,j) = bijP
        end if
     end do
  end do
  do i = 1, n3
     do j = 1, n3
        if (i .eq. j) then
           beta3(i,j) = biiS
        else
           beta3(i,j) = bijS
        end if
     end do
  end do

  ! --------- read adjacencies (PPo: n1 x n2, PPe: n2 x n3) ----------
  Adj1 = 0.0d0
  open(1001, file='PPo_adj_A.txt', status='old', action='read')
  do i = 1, n1
     read(1001, *) Adj1(i,1:n2)
  end do
  close(1001)

  Adj2 = 0.0d0
  open(1002, file='PPe_adj_A.txt', status='old', action='read')
  do i = 1, n2
     read(1002, *) Adj2(i,1:n3)
  end do
  close(1002)
  ! ------------------------------------------------------------------

  h  = 0.01d0
  t0 = 0.0d0

  ! Precompute degrees that depend only on adjacency (not on states)
  ! (We’ll reuse them every stage; max(.,1) guards avoid div-by-zero)
  do j = 1, n1
     D1A(j) = 0.0d0
     do k = 1, n2
        D1A(j) = D1A(j) + Adj1(j,k)
     end do
  end do
  do j = 1, n2
     D1P(j) = 0.0d0
     do k = 1, n1
        D1P(j) = D1P(j) + Adj1(k,j)
     end do
     D2P(j) = 0.0d0
     do k = 1, n3
        D2P(j) = D2P(j) + Adj2(j,k)
     end do
  end do
  do j = 1, n3
     D2S(j) = 0.0d0
     do k = 1, n2
        D2S(j) = D2S(j) + Adj2(k,j)
     end do
  end do

do j = 1,n1
!write(*,*)j, D1A(j)
end do

do j = 1,n2
!write(*,*)j, D1P(j)
end do

do j = 1,n2
!write(*,*)j, D2P(j)
end do

do j = 1,n3
!write(*,*)j, D2S(j)
end do
  !================== time stepping ==================
  do 100 i = 1, mm

    ! ---------- Stage 1 (k1) using A0, P0, S0 ----------
    ! A
    do j = 1, n1
       sum1(j) = 0.0d0
       do k = 1, n1
          sum1(j) = sum1(j) + beta1(j,k) * A0(k)
       end do
       sum2(j) = 0.0d0
       do k = 1, n2
          sum2(j) = sum2(j) + Adj1(j,k) * P0(k)
       end do
       den  = D1A(j)**(1.0d0-tt)       
       
       if(D1A(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (gamma0/den) * sum2(j)
       end if
       
       l = A0(j) * ( alpha - kappa - sum1(j)  &
           + (gain/(1.0d0 + hh*gain))         &
           - mA*(1.0d0 - f)*ubar ) + mu
       krkA(1,j) = l*h
    end do

    ! P
    do j = 1, n2
       sum3(j) = 0.0d0
       do k = 1, n2
          sum3(j) = sum3(j) + beta2(j,k) * P0(k)
       end do
       sum4(j) = 0.0d0
       do k = 1, n1
          sum4(j) = sum4(j) + Adj1(k,j) * A0(k)
       end do
       sum5(j) = 0.0d0
       do k = 1, n3
          sum5(j) = sum5(j) + Adj2(j,k) * S0(k)
       end do
       
       
       den  = D1P(j)**(1.0d0-tt)
       if(D1P(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (gamma0/den) * sum4(j)
       end if
       
       
       den2 = D2P(j)**(1.0d0-tt)
       if(D2P(j).eq.0)then
         gain2 = 0.0d0
       else   
         gain2 = (xi0/den2) * sum5(j)
       end if
       
       
       l = P0(j) * ( alpha - sum3(j)                 &
           + (gain/(1.0d0 + hh*gain))  - (gain2/(1.0d0 + hh*gain2))  ) + mu
       krkP(1,j) = l*h
    end do

    ! S
    do j = 1, n3
       sum6(j) = 0.0d0
       do k = 1, n3
          sum6(j) = sum6(j) + beta3(j,k) * S0(k)
       end do
       sum7(j) = 0.0d0
       do k = 1, n2
          sum7(j) = sum7(j) + Adj2(k,j) * P0(k)
       end do
       
       den  = D2S(j)**(1.0d0-tt)
       if(D2S(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (xi0/den) * sum7(j)
       end if
       
       
       l = S0(j) * ( alpha - kappa - sum6(j)  &
           + (gain/(1.0d0 + hh*gain))         &
           - mS*(1.0d0 - f)*ubar ) + mu
       krkS(1,j) = l*h
    end do

    ! ---------- Stage 2 (k2) at half-step ----------
    do j = 1, n1
       A(j) = A0(j) + krkA(1,j)/2.0d0
    end do
    do j = 1, n2
       P(j) = P0(j) + krkP(1,j)/2.0d0
    end do
    do j = 1, n3
       S(j) = S0(j) + krkS(1,j)/2.0d0
    end do
    t = t0 + h/2.0d0

    ! A
    do j = 1, n1
       sum1(j) = 0.0d0
       do k = 1, n1
          sum1(j) = sum1(j) + beta1(j,k) * A(k)
       end do
       sum2(j) = 0.0d0
       do k = 1, n2
          sum2(j) = sum2(j) + Adj1(j,k) * P(k)
       end do
       den  = D1A(j)**(1.0d0-tt)       
       
       if(D1A(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (gamma0/den) * sum2(j)
       end if
       
       l = A(j) * ( alpha - kappa - sum1(j)  &
           + (gain/(1.0d0 + hh*gain))        &
           - mA*(1.0d0 - f)*ubar ) + mu
       krkA(2,j) = l*h
    end do

    ! P
    do j = 1, n2
       sum3(j) = 0.0d0
       do k = 1, n2
          sum3(j) = sum3(j) + beta2(j,k) * P(k)
       end do
       sum4(j) = 0.0d0
       do k = 1, n1
          sum4(j) = sum4(j) + Adj1(k,j) * A(k)
       end do
       sum5(j) = 0.0d0
       do k = 1, n3
          sum5(j) = sum5(j) + Adj2(j,k) * S(k) 
       end do
       
       
       den  = D1P(j)**(1.0d0-tt)
       if(D1P(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (gamma0/den) * sum4(j)
       end if
       
       
       den2 = D2P(j)**(1.0d0-tt)
       if(D2P(j).eq.0)then
         gain2 = 0.0d0
       else   
         gain2 = (xi0/den2) * sum5(j)
       end if
       
       
       l = P(j) * ( alpha - sum3(j)                 &
           + (gain/(1.0d0 + hh*gain))  - (gain2/(1.0d0 + hh*gain2))  ) + mu
       krkP(2,j) = l*h
    end do

    ! S
    do j = 1, n3
       sum6(j) = 0.0d0
       do k = 1, n3
          sum6(j) = sum6(j) + beta3(j,k) * S(k)
       end do
       sum7(j) = 0.0d0
       do k = 1, n2
          sum7(j) = sum7(j) + Adj2(k,j) * P(k)
       end do
       
       den  = D2S(j)**(1.0d0-tt)
       if(D2S(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (xi0/den) * sum7(j)
       end if
       
       
       l = S(j) * ( alpha - kappa - sum6(j)  &
           + (gain/(1.0d0 + hh*gain))         &
           - mS*(1.0d0 - f)*ubar ) + mu
       krkS(2,j) = l*h
    end do

    ! ---------- Stage 3 (k3) at half-step ----------
    do j = 1, n1
       A(j) = A0(j) + krkA(2,j)/2.0d0
    end do
    do j = 1, n2
      P(j) = P0(j) + krkP(2,j)/2.0d0
    end do
    do j = 1, n3
      S(j) = S0(j) + krkS(2,j)/2.0d0
    end do
    t = t0 + h/2.0d0

    ! A
    do j = 1, n1
       sum1(j) = 0.0d0
       do k = 1, n1
          sum1(j) = sum1(j) + beta1(j,k) * A(k)
       end do
       sum2(j) = 0.0d0
       do k = 1, n2
          sum2(j) = sum2(j) + Adj1(j,k) * P(k)
       end do
       den  = D1A(j)**(1.0d0-tt)       
       
       if(D1A(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (gamma0/den) * sum2(j)
       end if
       
       l = A(j) * ( alpha - kappa - sum1(j)  &
           + (gain/(1.0d0 + hh*gain))        &
           - mA*(1.0d0 - f)*ubar ) + mu
       krkA(3,j) = l*h
    end do

    ! P
    do j = 1, n2
       sum3(j) = 0.0d0
       do k = 1, n2
          sum3(j) = sum3(j) + beta2(j,k) * P(k)
       end do
       sum4(j) = 0.0d0
       do k = 1, n1
          sum4(j) = sum4(j) + Adj1(k,j) * A(k)
       end do
       sum5(j) = 0.0d0
       do k = 1, n3
          sum5(j) = sum5(j) + Adj2(j,k) * S(k) 
       end do
       
       
       den  = D1P(j)**(1.0d0-tt)
       if(D1P(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (gamma0/den) * sum4(j)
       end if
       
       
       den2 = D2P(j)**(1.0d0-tt)
       if(D2P(j).eq.0)then
         gain2 = 0.0d0
       else   
         gain2 = (xi0/den2) * sum5(j)
       end if
       
       
       l = P(j) * ( alpha - sum3(j)                 &
           + (gain/(1.0d0 + hh*gain))  - (gain2/(1.0d0 + hh*gain2))  ) + mu
       krkP(3,j) = l*h
    end do

    ! S
    do j = 1, n3
       sum6(j) = 0.0d0
       do k = 1, n3
          sum6(j) = sum6(j) + beta3(j,k) * S(k)
       end do
       sum7(j) = 0.0d0
       do k = 1, n2
          sum7(j) = sum7(j) + Adj2(k,j) * P(k)
       end do
       
       den  = D2S(j)**(1.0d0-tt)
       if(D2S(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (xi0/den) * sum7(j)
       end if
       
       
       l = S(j) * ( alpha - kappa - sum6(j)  &
           + (gain/(1.0d0 + hh*gain))         &
           - mS*(1.0d0 - f)*ubar ) + mu
       krkS(3,j) = l*h
    end do

    ! ---------- Stage 4 (k4) at full step ----------
    do j = 1, n1
       A(j) = A0(j) + krkA(3,j)
    end do
    do j = 1, n2
       P(j) = P0(j) + krkP(3,j)
    end do
    do j = 1, n3
       S(j) = S0(j) + krkS(3,j)
    end do
    t = t0 + h

     ! A
    do j = 1, n1
       sum1(j) = 0.0d0
       do k = 1, n1
          sum1(j) = sum1(j) + beta1(j,k) * A(k)
       end do
       sum2(j) = 0.0d0
       do k = 1, n2
          sum2(j) = sum2(j) + Adj1(j,k) * P(k)
       end do
       den  = D1A(j)**(1.0d0-tt)       
       
       if(D1A(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (gamma0/den) * sum2(j)
       end if
       
       l = A(j) * ( alpha - kappa - sum1(j)  &
           + (gain/(1.0d0 + hh*gain))        &
           - mA*(1.0d0 - f)*ubar ) + mu
       krkA(4,j) = l*h
    end do

    ! P
    do j = 1, n2
       sum3(j) = 0.0d0
       do k = 1, n2
          sum3(j) = sum3(j) + beta2(j,k) * P(k)
       end do
       sum4(j) = 0.0d0
       do k = 1, n1
          sum4(j) = sum4(j) + Adj1(k,j) * A(k)
       end do
       sum5(j) = 0.0d0
       do k = 1, n3
          sum5(j) = sum5(j) + Adj2(j,k) * S(k) 
       end do
       
       
       den  = D1P(j)**(1.0d0-tt)
       if(D1P(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (gamma0/den) * sum4(j)
       end if
       
       
       den2 = D2P(j)**(1.0d0-tt)
       if(D2P(j).eq.0)then
         gain2 = 0.0d0
       else   
         gain2 = (xi0/den2) * sum5(j)
       end if
       
       
       l = P(j) * ( alpha - sum3(j)                 &
           + (gain/(1.0d0 + hh*gain))  - (gain2/(1.0d0 + hh*gain2))  ) + mu
       krkP(4,j) = l*h
    end do

    ! S
    do j = 1, n3
       sum6(j) = 0.0d0
       do k = 1, n3
          sum6(j) = sum6(j) + beta3(j,k) * S(k)
       end do
       sum7(j) = 0.0d0
       do k = 1, n2
          sum7(j) = sum7(j) + Adj2(k,j) * P(k)
       end do
       
       den  = D2S(j)**(1.0d0-tt)
       if(D2S(j).eq.0)then
         gain = 0.0d0
       else   
         gain = (xi0/den) * sum7(j)
       end if
       
       
       l = S(j) * ( alpha - kappa - sum6(j)  &
           + (gain/(1.0d0 + hh*gain))         &
           - mS*(1.0d0 - f)*ubar ) + mu
       krkS(4,j) = l*h
    end do
    
    ! ---------- Combine (single global step) ----------
    do j = 1, n1
       Anew(j) = A0(j) + ( krkA(1,j) + 2.0d0*(krkA(2,j)+krkA(3,j)) + krkA(4,j) ) / 6.0d0
    end do
    do j = 1, n2
       Pnew(j) = P0(j) + ( krkP(1,j) + 2.0d0*(krkP(2,j)+krkP(3,j)) + krkP(4,j) ) / 6.0d0
    end do
    do j = 1, n3
       Snew(j) = S0(j) + ( krkS(1,j) + 2.0d0*(krkS(2,j)+krkS(3,j)) + krkS(4,j) ) / 6.0d0
    end do

    ! advance one step and write
    do j = 1, n1
       A0(j) = Anew(j)
    end do
    do j = 1, n2
       P0(j) = Pnew(j)
    end do
    do j = 1, n3
       S0(j) = Snew(j)
    end do
    t0 = t0 + h
    
100 continue

     flag(i1+1)=A0(18)
     
!     write(*,*)xu,A0(18),  flag(i1+1),  flag(i1), flag(i1+1)-flag(i1)   
     if(abs( flag(i1+1) - flag(i1) ).gt.0.01d0)then
        xstar = xu -0.01d0                                                                                      !!! Change
        exit
     end if
     
     
end do

    write(201,*)gamma0, xi0, xstar
!    write(*,*)gamma0, xi0, xstar
    
end do
end do
end program pppp
