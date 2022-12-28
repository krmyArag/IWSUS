MODUlE module_sf_ahs_perp

CONTAINS

SUBROUTINE ahs_perp( HGT_ROOF, HWR, z, UA, N,                &
                     u1, u2, u3, u4, u5, u6, u7, u8, u9,       &
                     v1, v2, v3, v4, v5, v6, v7, v8, v9,       &
                     w1, w2, w3, w4, w5, w6, w7, w8, w9        &
                   )
implicit none

!integer, parameter :: N=100
integer :: N
real, intent(in ) :: HGT_ROOF, HWR, UA
real, dimension(N), intent(inout) :: z

real, dimension(N), intent(inout) :: u1, u2, u3, u4, u5, u6, u7, u8, u9
real, dimension(N), intent(inout) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
real, dimension(N), intent(inout) :: w1, w2, w3, w4, w5, w6, w7, w8, w9

print *, '+++++++++++perp++++++++++++++'

    call u7prof(HGT_ROOF, HWR, z, UA, u7, N)
    call u9prof(HGT_ROOF, HWR, z, UA, u9, N)
    u1=u7;u3=u9
    call u4prof(HGT_ROOF, HWR, z, UA, u4, N)
    call u6prof(HGT_ROOF, HWR, z, UA, u6, N)

    call v7prof(HGT_ROOF,HWR,z,UA,v7,N)
    call v9prof(HGT_ROOF,HWR,z,UA,v9,N)
    v4=0;v6=0;v1=-v7;v3=-v9

    call w7prof(HGT_ROOF, HWR, z, UA, w7, N)
    call w9prof(HGT_ROOF, HWR, z, UA, w9, N)
    w1=w7;w3=w9
    call w4prof(HGT_ROOF, HWR, z, UA, w4, N)
    call w6prof(HGT_ROOF, HWR, z, UA, w6, N)

    call u2prof(HGT_ROOF, HWR, z, UA, u2, N)
    u8=u2
    call v2prof(HGT_ROOF, HWR, z, UA, v2, N)
    v8=v2
    call w2prof(HGT_ROOF, HWR, z, UA, w2, N)
    w8=w2
    call u5prof(HGT_ROOF, HWR, z, UA, u5, N)
    call v5prof(HGT_ROOF, HWR, z, UA, v5, N)
    call w5prof(HGT_ROOF, HWR, z, UA, w5, N)


PRINT *, '******* AHS PROFILE ***********'
print *, 'MAX'
PRINT *, 'out_b', maxval(u1), maxval(v1), maxval(w1)
PRINT *, 'out_c', maxval(u2), maxval(v2), maxval(w2)
PRINT *, 'out_y', maxval(u3), maxval(v3), maxval(w3)
PRINT *, 'mid_b', maxval(u4), maxval(v4), maxval(w4)
PRINT *, 'mid_c', maxval(u5), maxval(v5), maxval(w5)
PRINT *, 'mid_y', maxval(u6), maxval(v6), maxval(w6)
PRINT *, 'en_b', maxval(u7), maxval(v7), maxval(w7)
PRINT *, 'en_c', maxval(u8), maxval(v8), maxval(w8)
PRINT *, 'en_y', maxval(u9), maxval(v9), maxval(w9)

print *, 'MIN'
PRINT *, 'out_b', minval(u1), minval(v1), minval(w1)
PRINT *, 'out_c', minval(u2), minval(v2), minval(w2)
PRINT *, 'out_y', minval(u3), minval(v3), minval(w3)
PRINT *, 'mid_b', minval(u4), minval(v4), minval(w4)
PRINT *, 'mid_c', minval(u5), minval(v5), minval(w5)
PRINT *, 'mid_y', minval(u6), minval(v6), minval(w6)
PRINT *, 'en_b', minval(u7), minval(v7), minval(w7)
PRINT *, 'en_c', minval(u8), minval(v8), minval(w8)
PRINT *, 'en_y', minval(u9), minval(v9), minval(w9)


END SUBROUTINE ahs_perp

subroutine u7prof(H, HWR, z, U0, uprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, uprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: u1, u2, u3, u4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4
    
!===================================

!                I      

    !------- ENDPOINT1--------------
    pos1=H*(0.2781*exp(-0.8123*HWR))

    z1=minloc(abs(z-pos1))
    
    !------- coeff of I ------------
    if (HWR .lt. 0.5) then
        u1=-0.5384*exp(0.5*HWR)+0.5304
    elseif (HWR .lt. 1) then
        u1=0.1944*HWR-0.2586
    else
        u1=0.09337*exp(-1.087*HWR)-0.09564
    end if

    !------- wind of I --------------
    u1=U0*u1
    a1=u1/log(pos1+1)

    uprof(1:z1(1))=a1*log(z(1:z1(1))+1)

!              II
    !------- ENDPOINT2--------------
    pos2=H
    z2=minloc(abs(z-pos2))

    !------- coeff of II ------------
    u2=U0*(0.2887*exp(-1.124*HWR)-0.1473)
    b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)**2
    uprof(z1(1)+1:z2(1))=a2*(z(z1(1)+1:z2(1))+b2)**2+c2

!              III
    !------- ENDPOINT2--------------
    pos3=H*(1.712*exp(-5.94*HWR)+1.612)
    z3=minloc(abs(z-pos3))
    u3=1.153*U0
    !------- coeff of III ------------
    b3=-pos3;c3=u3;a3=(u2-c3)/(pos2+b3)**2
    uprof(z2(1):z3(1))=a3*(z(z2(1):z3(1))+b3)**2+c3

!              IV
    !------- coeff of IV ------------
    c4=u0;b4=-0.1;a4=(u3-c4)/exp(b4*pos3)
    uprof(z3(1):N)=a4*exp(b4*z(z3(1):N))+c4

print *, '-----u7prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'u1', u1
print *, 'u2', u2
print *, 'u3', u3
print *, '----------------'
end subroutine u7prof

subroutine u9prof(H, HWR, z, U0, uprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, uprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3, pos4
    integer :: z1(1), z2(1), z3(1), z4(1) ! Endpoints
    real :: u1, u2, u3, u4, u5
    real :: a1, a2, a3, a4, a5,               &
    &       b1, b2, b3, b4, b5,               &
    &       c1, c2, c3, c4, c5
    
!===================================


!                I      

    !------- ENDPOINT1--------------
    pos1=0.1*H
    z1=minloc(abs(z-pos1))
    
    !------- coeff of I ------------
    u1=U0*(-0.3515*exp(-0.9082*HWR)+0.09882)
    a1=u1/log(pos1+1)

    !------- wind of I --------------
    uprof(1:z1(1))=a1*log(z(1:z1(1))+1)

!              II
    !------- ENDPOINT2--------------
    pos2=H*(-1.278*exp(-2.643*H)+0.8242)
    if (HWR .lt. 0.5) pos2=H*(-0.8667*HWR+0.9167)
    z2=minloc(abs(z-pos2))

    !------- coeff of II ------------
    u2=U0*(-0.3537*exp(-0.2449*HWR)+0.2514)
    if (HWR .lt. 0.5) u2=U0*(-0.4925*HWR+0.1848)
    if (HWR .gt.   1) u2=U0*(-0.3537*exp(-0.2449)+0.2514)
    a2=(u2-u1)/(pos2-pos1)
    b2=u2-a2*pos2

    !------- wind of II --------------
    uprof(z1(1):z2(1))=a2*z(z1(1):z2(1))+b2

!              III
    !------- ENDPOINT3--------------
    pos3=H;
    z3=minloc(abs(z-pos3))

    !------- coeff of III ------------
    u3=U0*(-0.09157*exp(1.713*HWR)+0.3768)
    if (HWR .gt. 1) u3=U0*(-0.09157*exp(1.173)+0.3768)
    b3=-pos2;c3=u2;a3=(u3-u2)/(pos3-pos2)**2

    !------- wind of III --------------
    uprof(z2(1):z3(1))=a3*(z(z2(1):z3(1))+b3)**2+c3

!              IV
    !------- ENDPOINT3--------------
    pos4=H*(0.4612*exp(-4.901*HWR)+1.364)
    z4=minloc(abs(z-pos4))

    !------- coeff of IV ------------
    u4=0.5688*U0
    b4=-pos3;c4=u3;a4=(u4-c4)/(pos4+b4)**2

    !------- wind of IV --------------
    uprof(z3(1):z4(1))=a4*(z(z3(1):z4(1))+b4)**2+c4

!              V
    !------- coeff of V ------------
    c5=U0;b5=-0.1;a5=(u4-c5)/exp(b5*pos4)

    !------- wind of IV --------------
    uprof(z4(1):N)=a5*exp(b5*z(z4(1):N))+c5

print *, '-----u9prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'u1', u1
print *, 'u2', u2
print *, 'u3', u3
print *, 'u4', u4
print *, '----------------'
end subroutine u9prof

subroutine u4prof(H, HWR, z, U0, uprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, uprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: u1, u2, u3, u4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4
    
!===================================


!                I      

    !------- ENDPOINT1--------------
    pos1=H*(-0.01956*exp(2.563*HWR)+0.3038);
    if (HWR .gt. 1)  pos1=H*(-0.01956*exp(2.563)+0.3038);
    z1=minloc(abs(z1-pos1))
    
    !------- coeff of I ------------
    u1=u0*(-0.6492*exp(-2.172*HWR)+0.1069)
    a1=u1/log(pos1+1)

    !------- wind of I --------------
    uprof(1:z1(1))=a1*log(z(1:z1(1))+1)

!              II
    !------- ENDPOINT2--------------
    pos2=H
    z2=minloc(abs(z-pos2))

    !------- coeff of II ------------
    u2=U0*(0.4947*HWR-0.327)
    if (HWR .ge. 0.5) u2=U0*(0.2591*exp(-3.461*HWR)-0.1188)
    b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)**2;

    uprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2;

!              III
    !------- ENDPOINT2--------------
    pos3=H*(1.069*exp(-2.746*HWR)+1.254)
    z3=minloc(abs(z-pos3))

    !------- coeff of III ------------
    u3=u0*(0.4118*exp(-1.276*HWR)+0.719)
    if (HWR .lt. 0.5) u3=u0*(0.118*HWR+0.8775)
    b3=-pos2;c3=u2;a3=(u3-c3)/(pos3+b3)**2

    uprof(z2(1):z3(1))=a3*(z(z2(1):z3(1))+b3)**2+c3

!              IV
    !------- coeff of IV ------------
    c4=U0;b4=-0.2;a4=(u3-c4)/exp(b4*pos3)

    uprof(z3(1):N)=a4*exp(b4*z(z3(1):N))+c4

!    print *, uprof

print *, '-----u4prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'u1', u1
print *, 'u2', u2
print *, 'u3', u3
print *, '----------------'

end subroutine u4prof

subroutine u6prof(H, HWR, z, U0, uprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, uprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3, pos4
    integer :: z1(1), z2(1), z3(1), z4(1) ! Endpoints
    real :: u1, u2, u3, u4, u5
    real :: a1, a2, a3, a4, a5,               &
    &       b1, b2, b3, b4, b5,              &
    &       c1, c2, c3, c4, c5
    
!===================================


!                I      

    !------- ENDPOINT1--------------
    pos1=H*(0.4667*HWR+0.08333)
    if (HWR .gt. 0.5) pos1=H*(4.551*exp(-5.545*HWR)+0.03222)
    z1=minloc(abs(z-pos1))
    
    !------- coeff of I ------------
    u1=U0*(0.01338*(exp(2.661*HWR))-0.1092)
    if (HWR .ge. 1) u1=U0*(0.01338*(exp(2.661))-0.1092)
    a1=u1/log(pos1+1)

    uprof(1:z1(1))=a1*log(z(1:z1(1))+1)

!              II
    !------- ENDPOINT2--------------
    pos2=0.8*H;z2=minloc(abs(z-pos2))

    !------- coeff of II ------------
    u2=U0*(-0.2438*exp(-3.828*HWR)-0.01533)
    if (HWR .lt. 0.5) u2=U0*(-0.04493*HWR-0.02887)
    b2=-pos1;c2=u1;a2=(u2-c2)/((pos2+b2)**2)

    uprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2
!              III
    !------- ENDPOINT2--------------
    pos3=H*(1.286*exp(-6.621*HWR)+1.054)
    if (HWR .le. 0.25) pos3=H*(1.286*exp(-6.621*0.25)+1.054)
    z3=minloc(abs(z-pos3))

    !------- coeff of III ------------
    u3=U0*(-0.54007*exp(0.9084*HWR)+0.5225)
    if (HWR .gt. 0.5) u3=U0*(-0.4343*exp(-1.968*HWR)-0.1667)
    b3=-pos3;c3=u3;a3=(u2-c3)/((pos2+b3)**2)

    uprof(z2(1):z3(1))=a3*(z(z2(1):z3(1))+b3)**2+c3

!              IV
    !------- coeff of IV ------------
    pos4=H*(3.289*exp(-5.56*HWR)+1.512);z4=minloc(abs(z-pos4))

    u4=0.6038*U0;b4=-pos3;c4=u3;a4=(u4-c4)/(pos4+b4)**2
    uprof(z3(1):z4(1))=a4*(z(z3(1):z4(1))+b4)**2+c4

!              V
    c5=U0;b5=-0.3;a5=(u4-c5)/exp(b5*pos4)
    uprof(z4(1):N)=a5*exp(b5*z(z4(1):N))+c5

!    print *, uprof

print *, '-----u6prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'u1', u1
print *, 'u2', u2
print *, 'u3', u3
print *, 'u4', u4
print *, '----------------'

end subroutine u6prof

subroutine v7prof(H, HWR, z, U0, vprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, vprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: v1, v2, v3, v4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4
    
!===================================

!                I      

    !------- ENDPOINT1--------------
    pos1=H*(0.4051*exp(-3.592*HWR)+0.03452)
    z1=minloc(abs(z-pos1))
    
    !------- coeff of I ------------
    v1=u0*(-0.5917*exp(-2.86*HWR)+0.09746)
    if (HWR .lt. 0.5) v1=u0*(-0.5917*exp(-2.86*0.5)+0.09746)
    a1=v1/log(pos1+1)

    vprof(1:z1(1))=a1*log(z(1:z1(1))+1)

!              II
    !------- ENDPOINT2--------------
    pos2=H*(1.079*exp(0.2143*HWR)-0.003025)
    if (HWR .gt. 1) pos2=H*(0.7993*exp(-0.3397*HWR)+0.7643)
    z2=minloc(abs(z-pos2))

    !------- coeff of II ------------
    v2=u0*(-0.06966*exp(-1.58*HWR)-0.1733)
    b2=-pos2;c2=v2;a2=(v1-c2)/(pos1+b2)**2

    vprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2

!              III
    b3=-0.1;a3=v2/exp(b3*pos2)

    vprof(z2(1):N)=A3*exp(b3*z(z2(1):N))
print *, '-----v7prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'v1', v1
print *, 'v2', v2
print *, '----------------'

end subroutine v7prof

subroutine v9prof(H, HWR, z, U0, vprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, vprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: v1, v2, v3, v4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4
    
!===================================


!                I      

    !------- ENDPOINT1--------------
    pos1=H*(1.519*exp(-4.394*HWR)+0.03125)
    if (HWR .lt. 0.5) pos1=H*(0.4*HWR)
    z1=minloc(abs(z-pos1))
    
    !------- coeff of I ------------
    v1=u0*(-0.02976*exp(2.984*HWR)+0.5034)
    if (HWR .gt. 1) v1=u0*(-0.02976*exp(2.984)+0.5034)
    a1=v1/log(pos1+1)

    vprof(1:z1(1))=a1*log(z(1:z1(1))+1)

!              II
    !------- ENDPOINT2--------------
    pos2=H*(-2.55*exp(0.2608*HWR)+3.722)
    if (HWR .gt. 0.5) pos2=H*(-2.55*exp(0.2608*0.5)+3.722)
    if (HWR .lt. 0.25) pos2=H
    z2=minloc(abs(z-pos2))

    !------- coeff of II ------------
    v2=0.1498*u0
    b2=-pos1;c2=v1;a2=(v2-c2)/(pos2+b2)**2

    vprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2

!              III
    !------- ENDPOINT2--------------
    pos3=H*(-2.391*exp(0.545*HWR)+4.673)
    if (HWR .gt. 0.5) pos3=H*(-2.391*exp(0.545*0.5)+4.673)
    z3=minloc(abs(z-pos3))

    !------- coeff of III ------------
    v3=u0*(0.1326*exp(-2.765*HWR)-0.1007)
    b3=-pos3;c3=v3;a3=(v2-c3)/(pos2+b3)**2

    vprof(z2(1):z3(1))=a3*(z(z2(1):z3(1))+b3)**2+c3

!              IV
    !------- coeff of IV ------------
    c4=0;b4=-0.1;a4=v3/exp(b4*pos3)

    vprof(z3(1):N)=a4*exp(b4*z(z3(1):N))

!    print *, uprof
print *, '-----v9prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'v1', v1
print *, 'v2', v2
print *, 'v3', v3
print *, '----------------'


end subroutine v9prof

subroutine w7prof(H, HWR, z, U0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: w1, w2, w3, w4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4
    
!===================================


!                I      

    !------- ENDPOINT1--------------
    pos1=H*(-5.284*exp(-6.056*HWR)+0.9939)
    if (HWR .lt. 0.5) pos1=H*(-5.284*exp(-6.056*0.5)+0.9939)
    z1=minloc(abs(z-pos1))
    
    !------- coeff of I ------------
    w1=u0*(0.04559*exp(2.746*HWR)+0.002081)
    if (HWR .ge. 0.5) w1=u0*(0.3718*exp(-5.31*HWR)+0.1556)
    b1=-pos1;c1=w1;a1=-c1/(b1)**2

    wprof(1:z1(1))=a1*(z(1:z1(1))+b1)**2+c1

!              II
    !------- ENDPOINT2--------------
    pos2=H*(0.3585*exp(-1.898*HWR)+1.057)
    z2=minloc(abs(z-pos2))

    !------- coeff of II ------------
    w2=u0*(-0.307*exp(-5.844*HWR)+0.1185)
    a2=(w2-w1)/(pos2-pos1);b2=w2-a2*pos2

    wprof(z1(1):z2(1))=a2*z(z1(1):z2(1))+b2

!              III
    !------- ENDPOINT2--------------
    pos3=H*(1.475*exp(-3.011*HWR)+1.318)
    z3=minloc(abs(z-pos3))

    !------- coeff of III ------------
    w3=u0*(-0.426*exp(-1.736*HWR)+0.4536)
    a3=(w3-w2)/(pos3-pos2);b3=w3-a3*pos3

    wprof(z2(1):z3(1))=a3*z(z2(1):z3(1))+b3

!              IV
    !------- coeff of IV ------------
    c4=u0*0.05;b4=-0.1;a4=(w3-c4)/exp(b4*pos3)
    wPROF(z3(1):N)=a4*exp(b4*z(z3(1):N))+c4

print *, '-----w7prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'w1', w1
print *, 'w2', w2
print *, 'w3', w3
print *, '----------------'
end subroutine w7prof

subroutine w9prof(H, HWR, z, U0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: w1, w2, w3, w4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4
    
!===================================

!                I      
    !------- ENDPOINT1--------------
    pos1=H*(-3.712*exp(-5.28*HWR)+0.9836)
    if (HWR .lt. 0.5) pos1=H*(-3.712*exp(-5.28*0.5)+0.9836)
    z1=minloc(abs(z-pos1))
    
    !------- coeff of I ------------
    w1=u0*(-0.04245*exp(4.028*HWR)-0.002003)
    if (HWR .gt. 0.5) w1=u0*(-2.098*exp(-5.081*HWR)-0.1547)
    b1=0;c1=0;a1=w1/(pos1)**2

    wprof(1:z1(1))=a1*(z(1:z1(1))+b1)**2+c1

!              II
    !------- ENDPOINT2--------------
    pos2=H*(8.236*exp(-5.263*HWR)+1.857)
    z2=minloc(abs(z-pos2))

    !------- coeff of II ------------
    w2=u0*(-0.7792*exp(-2.273*HWR)+0.2835)
    if (HWR .lt. 0.5) w2=u0*(-0.7792*exp(-2.273*0.5)+0.2835)
    b2=-pos2;c2=w2;a2=(w1-c2)/(pos1+b2)**2
    wprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2

!              III
    !------- ENDPOINT2--------------
    wprof(z2(1):N)=w2

print *, '-----w9prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'w1', w1
print *, 'w2', w2
print *, '----------------'
end subroutine w9prof

subroutine w4prof(H, HWR, z, U0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: w1, w2, w3, w4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4
    
!===================================

!                I      
    !------- ENDPOINT1--------------
    pos1=H*0.9;z1=minloc(abs(z-pos1))
    
    !------- coeff of I ------------
    w1=u0*(0.5417*HWR+0.03767)
    if (HWR .ge. 0.5) w1=u0*(9.341*exp(-8.185*HWR)+0.1524)
    b1=-pos1;c1=w1;a1=-c1/(b1)**2

    wprof(1:z1(1))=a1*(z(1:z1(1))+b1)**2+c1

!              II
    !------- ENDPOINT2--------------
    pos2=H*(0.5214*exp(-1.5*HWR)+1.055)
    z2=minloc(abs(z-pos2))

    !------- coeff of II ------------
    w2=u0*(0.225*HWR+0.05438)
    if (HWR .gt. 0.5) w2=u0*(0.904*exp(-5.519*HWR)+0.1201)
    a2=(w2-w1)/(pos2-pos1);b2=w2-a2*pos2

    wprof(z1(1):z2(1))=a2*z(z1(1):z2(1))+b2

!              III
    !------- ENDPOINT2--------------
    pos3=(1.534*exp(-2.733*HWR)+1.337)*H
    z3=minloc(abs(z-pos3))

    !------- coeff of III ------------
    w3=u0*(-0.4408*exp(-2.875*HWR)+0.5083)
    a3=(w3-w2)/(pos3-pos2);b3=w3-a3*pos3

    wprof(z2(1):z3(1))=a3*z(z2(1):z3(1))+b3

!              IV
    !------- coeff of IV ------------
    c4=u0/10;b4=-0.1;a4=(w3-c4)/exp(b4*pos3)
    wprof(z3(1):n)=a4*exp(b4*z(z3(1):n))+c4

!    print *, uprof
print *, '-----w4prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'w1', w1
print *, 'w2', w2
print *, 'w3', w3
print *, '----------------'


end subroutine w4prof

subroutine w6prof(H, HWR, z, U0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: w1, w2, w3, w4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4
    
!===================================

!                I      
    !------- ENDPOINT1--------------
    pos1=H*(-0.4067*HWR+0.7757)
    if (HWR .ge. 0.5) pos1=H*(0.0219*exp(1.688*HWR)+0.5214)
    z1=minloc(abs(z-pos1))
    
    !------- coeff of I ------------
    w1=0.2472*HWR+0.2858
    if (HWR .ge. 0.5) w1=4.77*exp(-5.175*HWR)+0.05055
    w1=u0*w1;b1=-pos1;c1=w1;a1=-c1/b1**2

    wprof(1:z1(1))=a1*(z(1:z1(1))+b1)**2+c1

!              II
    !------- ENDPOINT2--------------
    pos2=H*(4.424*exp(-4.943*HWR)+1.005)
    z2=minloc(abs(z-pos2))

    !------- coeff of II ------------
    w2=u0*(0.4383*exp(-4.059*HWR)-0.1127)
    a2=(w2-w1)/(pos2-pos1);b2=w2-a2*pos2

    wprof(z1(1):z2(1))=a2*z(z1(1):z2(1))+b2

!              III
    !------- ENDPOINT2--------------
    pos3=H*(2.654*exp(-2.872*HWR)+1.73)
    z3=minloc(abs(z-pos3))

    !------- coeff of III ------------
    w3=u0*(-1.275*exp(-0.6783*HWR)+0.9835)
    b3=-pos2;c3=w2;a3=(w3-c3)/(pos3+b3)**2

    wprof(z2(1):z3(1))=a3*(z(z2(1):z3(1))+b3)**2+c3

!              IV
    wprof(z3(1):N)=w3

print *, '-----w6prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'w1', w1
print *, 'w2', w2
print *, '----------------'
end subroutine w6prof

subroutine u2prof(H, HWR, z, U0, uprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, uprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3, pos4
    integer :: z1(1), z2(1), z3(1), z4(1) ! Endpoints
    real :: u1, u2, u3, u4, u5
    real :: a1, a2, a3, a4, a5,               &
    &       b1, b2, b3, b4, b5,              &
    &       c1, c2, c3, c4, c5
    
!===================================


pos1=H*(0.1789*exp(-1.181*HWR));
if (HWR .lt. 0.5) pos1=H*(0.1789*exp(-1.181*0.5));
z1=minloc(abs(z-pos1));

u1=u0*(-1.26*exp(-2.359*HWR));
if (HWR .lt. 0.5) u1=u0*(-0.7688*exp(0.7596*HWR)+0.7303);
a1=u1/log(pos1+1);
uprof(1:z1(1))=a1*log(z(1:z1(1))+1);

pos2=H;
z2=minloc(abs(z-pos2));

u2=u0*(0.9269*exp(-2.71*HWR)-0.1778);
if (HWR .lt. 0.5) u2=u0*(0.9269*exp(-2.71*0.5)-0.1778);
b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)**2;
uprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2;

pos3=H*(1.196*exp(-2.244*HWR)+1.69);
z3=minloc(abs(z-pos3));
u3=u0*1.15;
b3=-pos3;c3=u3;a3=(u2-c3)/(pos2+b3)**2;
uprof(z2(1):z3(1))=a3*(z(z2(1):z3(1))+b3)**2+c3;

c4=u0;
b4=-0.1;
a4=(u3-c4)/exp(b4*pos3);
uprof(z3(1):N)=a4*exp(b4*z(z3(1):N))+c4;

print *, '-----u2prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'u1', u1
print *, 'u2', u2
print *, 'u3', u3
print *, '----------------'

end subroutine u2prof

subroutine u5prof(H, HWR, z, U0, uprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, uprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3, pos4
    integer :: z1(1), z2(1), z3(1), z4(1) ! Endpoints
    real :: u1, u2, u3, u4, u5
    real :: a1, a2, a3, a4, a5,               &
    &       b1, b2, b3, b4, b5,              &
    &       c1, c2, c3, c4, c5
    
!===================================

pos1=H*(1.519*exp(-4.394*HWR)+0.03125);
if (HWR .lt. 0.5) pos1=H*(0.6628*exp(0.5*HWR)-0.651);
z1=minloc(abs(z-pos1));

u1=u0*(-1.035*exp(-1.36*HWR)+0.3217);
a1=u1/log(pos1+1);
uprof(1:z1(1))=a1*log(z(1:z1(1))+1);

pos2=H;
z2=minloc(abs(z-pos2));

u2=u0*(0.4586*exp(0.6027*HWR)-0.7929);
if (HWR .gt. 0.5) u2=u0*(0.4586*exp(0.6027*0.5)-0.7929);
b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)**2;
uprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2;

pos3=H*(1.712*exp(-5.94*HWR)+1.612);
z3=minloc(abs(z-pos3));
u3=1.153*u0;
b3=-pos3;c3=u3;a3=(u2-c3)/(pos2+b3)**2;
uprof(z2(1):z3(1))=a3*(z(z2(1):z3(1))+b3)**2+c3;

c4=u0;
b4=-0.1;
a4=(u3-c4)/exp(b4*pos3);
uprof(z3(1):N)=a4*exp(b4*z(z3(1):N))+c4;

print *, '-----u5prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'u1', u1
print *, 'u2', u2
print *, 'u3', u3
print *, '----------------'
end subroutine u5prof

subroutine v2prof(H, HWR, z, U0, vprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, vprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3, pos4
    integer :: z1(1), z2(1), z3(1), z4(1) ! Endpoints
    real :: u1, u2, u3, u4, u5
    real :: a1, a2, a3, a4, a5,               &
    &       b1, b2, b3, b4, b5,              &
    &       c1, c2, c3, c4, c5
    
!===================================

pos1=H*(0.5864*exp(-0.1093*HWR)-0.4646);
z1=minloc(abs(z-pos1));
u1=u0*(0.4113*exp(-1.663*HWR)+0.002);
a1=u1/log(pos1+1);
vprof(1:z1(1))=a1*log(z(1:z1(1))+1);

pos2=H*1.379;
z2=minloc(abs(z-pos2));

u2=u0*(-1.642*exp(-6.573*HWR)-0.129);
if (HWR .lt. 0.5) u2=u0*(0.3207*exp(-0.8731*HWR)-0.3977);
b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)**2;
vprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2;

c4=0;
b4=-0.1;
a4=(u2-c4)/exp(b4*pos2);
vprof(z2(1):N)=a4*exp(b4*z(z2(1):N))+c4;
print *, '-----v2prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'v1', u1
print *, 'v2', u2
print *, '----------------'
end subroutine v2prof

subroutine v5prof(H, HWR, z, U0, vprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, vprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3, pos4
    integer :: z1(1), z2(1), z3(1), z4(1) ! Endpoints
    real :: u1, u2, u3, u4, u5
    real :: a1, a2, a3, a4, a5,               &
    &       b1, b2, b3, b4, b5,              &
    &       c1, c2, c3, c4, c5
    
!===================================

pos1=H*(0.5864*exp(-0.1093*HWR)-0.4646);
z1=minloc(abs(z-pos1));

u1=u0*(0.4113*exp(-1.663*HWR)+0.002);
a1=u1/log(pos1+1);
vprof(1:z1(1))=a1*log(z(1:z1(1))+1);

pos2=H*1.379;
z2=minloc(abs(z-pos2));

u2=u0*(-1.642*exp(-6.573*HWR)-0.129);
if (HWR .lt. 0.5) u2=u0*(.3207*exp(-0.8731*HWR)-0.3977);
b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)**2;
vprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2;

c4=0;
b4=-0.1;
a4=(u2-c4)/exp(b4*pos2);
vprof(z2(1):N)=a4*exp(b4*z(z2(1):N))+c4;

print *, '-----v5prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'v1', u1
print *, 'v2', u2
print *, '----------------'
end subroutine v5prof

subroutine w2prof(H, HWR, z, U0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3, pos4
    integer :: z1(1), z2(1), z3(1), z4(1) ! Endpoints
    real :: u1, u2, u3, u4, u5
    real :: a1, a2, a3, a4, a5,               &
    &       b1, b2, b3, b4, b5,              &
    &       c1, c2, c3, c4, c5
    
!===================================

pos1=H*(1.44*exp(-0.3646*HWR));
if (HWR .lt. 0.5) pos1=H*(0.6936*exp(0.8402*HWR)+0.1443);
z1=minloc(abs(z-pos1));

u1=u0*(-0.1257*exp(-2.869*HWR)-0.0051);
a1=u1/log(pos1+1);
wprof(1:z1(1))=a1*log(z(1:z1(1))+1);

pos2=H*(4.052*exp(-3.593*HWR)+1.345);
z2=minloc(abs(z-pos2));

u2=u0*(-0.5287*exp(-0.6914*HWR)+0.5054);

b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)**2;
wprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2;

c4=0;
b4=-0.01;
a4=(u2-c4)/exp(b4*pos2);
wprof(z2(1):N)=a4*exp(b4*z(z2(1):N))+c4;

print *, '-----w2prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'w1', u1
print *, 'w2', u2
print *, '----------------'
end subroutine w2prof

subroutine w5prof(H, HWR, z, U0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, U0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i
    
!===================================
!           inner vars
!===================================
    real :: pos1, pos2, pos3, pos4
    integer :: z1(1), z2(1), z3(1), z4(1) ! Endpoints
    real :: u1, u2, u3, u4, u5
    real :: a1, a2, a3, a4, a5,               &
    &       b1, b2, b3, b4, b5,              &
    &       c1, c2, c3, c4, c5
    
!===================================

pos1=H*(0.8847*exp(-1.694*HWR)+0.1875);
if (HWR .lt. 1) pos1=H*(0.8847*exp(-1.694*1)+0.1875);
z1=minloc(abs(z-pos1));

u1=u0*(7.437*exp(-6.384*HWR)+0.02);
if (HWR .lt. 0.5) u1=u0*(7.437*exp(-6.384*0.5)+0.02);
a1=u1/log(pos1+1);
wprof(1:z1(1))=a1*log(z(1:z1(1))+1);

pos2=H*(2.259*exp(-2.836*HWR)+0.8526);
z2=minloc(abs(z-pos2));

u2=u0*(-0.422*exp(0.5*HWR)+0.6952);
if (HWR .gt. 1) u2=u0*(-0.422*exp(0.5*1)+0.6952);
if (HWR .lt. 0.5) u2=u0*(-0.422*exp(0.5*0.5)+0.6952);

b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)**2;
wprof(z1(1):z2(1))=a2*(z(z1(1):z2(1))+b2)**2+c2;

pos3=H*(2.326*exp(-2.575*HWR)+1.496);
z3=minloc(abs(z-pos3));

u3=u0*(-0.3842*exp(-0.9443*HWR)+0.4714);

b3=-pos2;c3=u2;a3=(u3-c3)/(pos3+b3)**2;
wprof(z2(1):z3(1))=a3*(z(z2(1):z3(1))+b3)**2+c3;

c4=0;
b4=-0.01;
a4=(u3-c4)/exp(b4*pos3);
wprof(z3(1):N)=a4*exp(b4*z(z3(1):N))+c4;

print *, '-----w5prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', u0
print *, 'w1', u1
print *, 'w2', u2
print *, 'w3', u3
print *, '----------------'
end subroutine w5prof

END MODULE module_sf_ahs_perp
