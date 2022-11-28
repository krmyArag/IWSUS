MODUlE module_sf_ahs_para

CONTAINS

SUBROUTINE ahs_para( HGT_ROOF, HWR, z, UA, N,                &
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

print *, '++++++++++   para   +++++++++++++'

    call u1prof(HGT_ROOF, HWR, z, UA, u1, N)
    u2=0;u3=-u1
    call u4prof(HGT_ROOF, HWR, z, UA, u4, N)
    u5=0;u6=-u4
    call u7prof(HGT_ROOF, HWR, z, UA, u7, N)
    u8=0;u9=-u7

    call v1prof(HGT_ROOF,HWR,z,UA,v1,N)
    call v2prof(HGT_ROOF,HWR,z,UA,v2,N)
    v3=v1
    call v4prof(HGT_ROOF,HWR,z,UA,v4,N)
    call v5prof(HGT_ROOF,HWR,z,UA,v5,N)
    v6=v4
    call v7prof(HGT_ROOF,HWR,z,UA,v7,N)
    call v8prof(HGT_ROOF,HWR,z,UA,v8,N)
    v9=v7

    call w1prof(HGT_ROOF,HWR,z,UA,w1,N)
    call w2prof(HGT_ROOF,HWR,z,UA,w2,N)
    w3=w1
    call w4prof(HGT_ROOF,HWR,z,UA,w4,N)
    call w5prof(HGT_ROOF,HWR,z,UA,w5,N)
    w6=w4
    call w7prof(HGT_ROOF,HWR,z,UA,v7,N)
    call w8prof(HGT_ROOF,HWR,z,UA,v8,N)
    w9=w7
    
!IF (.true.) then
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

!ENDIF

END SUBROUTINE ahs_para

subroutine u1prof(H, HWR, z, V0, uprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, uprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: u1, u2, u3, u4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H
        z1=minloc(abs(z-p1))
        u1=V0*(-0.04175*exp(-9.374*HWR)+0.02567)
    
        b1=-H;c1=u1;a1=-c1/b1**2
    
        uprof(1:z1(1))=0.66*(a1*(z(1:z1(1))+b1)**2+c1)
    
        c3=0;b3=-0.05;a3=(u1-c3)/exp(b3*p1);
        uprof(z1(1)+1:N)=0.66*(a3*exp(b3*z(z1(1)+1:N))+c3);

print *, '-----u1prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'u1', u1
print *, '----------------'

end subroutine u1prof

subroutine u4prof(H, HWR, z, V0, uprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, uprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: u1, u2, u3, u4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H;
        z1=minloc(abs(z-p1));
        u1=V0*(-0.04175*exp(-9.374*HWR)+0.02567);

        b1=-H;c1=u1;a1=-c1/b1**2;

        uprof(1:z1(1))=a1*(z(1:z1(1))+b1)**2+c1;

        c3=0;b3=-0.05;
        a3=(u1-c3)/exp(b3*p1);
        uprof(z1(1)+1:N)=a3*exp(b3*z(z1(1):N))+c3;

print *, '-----u4prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'u1', u1
print *, '----------------'
end subroutine u4prof

subroutine u7prof(H, HWR, z, V0, uprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, uprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: u1, u2, u3, u4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H*(0.7662*exp(-3.718*HWR)+0.0975);
        z1=minloc(abs(z-p1));
        u1=V0*(0.7217*exp(-8.564*HWR)-0.1092);
        if (HWR>=1) u1=V0*(0.7217*exp(-8.564*1)-0.1092);
        a1=u1/p1;
        uprof(1:z1(1))=a1*z(1:z1(1));

        p2=H;
        z2=minloc(abs(z-p2));
        u2=V0*(0.02603*exp(-3.252*HWR)+0.009276);

        b2=-p2;c2=u2;a2=(u1-c2)/(p1+b2)**2;
        uprof(z1(1)+1:z2(1))=a2*(z(z1(1)+1:z2(1))+b2)**2+c2;

        c3=0;
        b3=-0.05;
        a3=(u2-c3)/exp(b3*p2);
        uprof(z2(1)+1:N)=a3*exp(b3*z(z2(1):N))+c3;

print *, '-----u7prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'u1', u1
print *, 'u2', u2
print *, '----------------'
end subroutine u7prof

subroutine v1prof(H, HWR, z, V0, vprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, vprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: v1, v2, v3, v4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H*(0.5948*exp(-1.288*HWR)+0.1745);
        z1=minloc(abs(z-p1));
        v1=V0*(0.9);

        a1=v1/log(p1+1);
        vprof(1:z1(1))=0.9*a1*log(z(1:z1(1))+1);

        p2=H*(2.206*exp(-3.593*HWR)+1.173);
        z2=minloc(abs(z-p2));
        v2=V0*(-0.4812*exp(-0.529*HWR)+1.505);
        a2=(v2-v1)/(p2-p1);b2=v2-a2*p2;
        vprof(z1(1)+1:z2(1))=0.9*(a2*z(z1(1)+1:z2(1))+b2);

        c3=V0;
        b3=-0.05;
        a3=(v2-c3)/exp(b3*p2);
        vprof(z2(1):N)=0.9*(a3*exp(b3*z(z2(1):N))+c3);

print *, '-----v1prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'v1', v1
print *, 'v2', v2
print *, '----------------'
end subroutine v1prof

subroutine v2prof(H, HWR, z, V0, vprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, vprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: v1, v2, v3, v4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H*(0.901*exp(-4.313*HWR)+0.194);
        z1=minloc(abs(z-p1));
        v1=V0*(-0.258*exp(-1.673*HWR)+1.2)
        
        a1=v1/log(p1+1);
        vprof(1:z1(1))=0.9*a1*log(z(1:z1(1))+1);
        
        c2=V0;
        b2=-0.05;
        a2=(v1-c2)/exp(b2*p1);
        vprof(z1(1)+1:N)=0.9*(a2*exp(b2*z(z1(1)+1:N))+c2);

print *, '-----v2prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'v1', v1
print *, '----------------'
end subroutine v2prof

subroutine v4prof(H, HWR, z, V0, vprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, vprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: v1, v2, v3, v4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H*(0.5948*exp(-1.288*HWR)+0.1745);
        z1=minloc(abs(z-p1));
        v1=V0*(0.9);
        
        a1=v1/log(p1+1);
        vprof(1:z1(1))=a1*log(z(1:z1(1))+1);
        
        p2=H*(2.206*exp(-3.593*HWR)+1.173);
        z2=minloc(abs(z-p2));
        v2=V0*(-0.4812*exp(-0529*HWR)+1.505);
        a2=(v2-v1)/(p2-p1);b2=v2-a2*p2;
        vprof(z1(1)+1:z2(1))=a2*z(z1(1)+1:z2(1))+b2;

        c3=V0;
        b3=-0.05;
        a3=(v2-c3)/exp(b3*p2);
        vprof(z2(1)+1:N)=a3*exp(b3*z(z2(1)+1:N))+c3;

print *, '-----v4prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'v1', v1
print *, 'v2', v2
print *, '----------------'
end subroutine v4prof 

subroutine v5prof(H, HWR, z, V0, vprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, vprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: v1, v2, v3, v4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H*(0.901*exp(-4.313*HWR)+0.194);
        z1=minloc(abs(z-p1));
        v1=V0*(-0.258*exp(-1.673*HWR)+1.2);
        
        a1=v1/log(p1+1);
        vprof(1:z1(1))=a1*log(z(1:z1(1))+1);
        
        c2=V0;
        b2=-0.05;
        a2=(v1-c2)/exp(b2*p1);
        vprof(z1(1)+1:N)=a2*exp(b2*z(z1(1):N))+c2;

print *, '-----v5prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'v1', v1
print *, '----------------'
end subroutine v5prof 

subroutine v7prof(H, HWR, z, V0, vprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, vprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: v1, v2, v3, v4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H*(0.5948*exp(-1.288*HWR)+0.1745);
        z1=minloc(abs(z-p1));
        v1=V0*(1.05);
        
        a1=v1/log(p1+1);
        vprof(1:z1(1))=a1*log(z(1:z1(1))+1);
        
        p2=H;
        z2=minloc(abs(z-p2));
        v2=V0*(-0.0534*exp(-3.028*HWR)+1.102);
        a2=(v2-v1)/(p2-p1);b2=v2-a2*p2;
        vprof(z1(1)+1:z2(1))=a2*z(z1(1)+1:z2(1))+b2;
        
        c3=V0;
        b3=-0.05;
        a3=(v2-c3)/exp(b3*p2);
        vprof(z2(1)+1:N)=a3*exp(b3*z(z2(1)+1:N))+c3;

print *, '-----v7prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'v1', v1
print *, 'v2', v2
print *, '----------------'
end subroutine v7prof

subroutine v8prof(H, HWR, z, V0, vprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, vprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: v1, v2, v3, v4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H*(0.5697*exp(-1.487*HWR)+0.049);
        z1=minloc(abs(z-p1));
        v1=V0*(-0.205*exp(-1.619*HWR)+1.192);
        
        a1=v1/log(p1+1);
        vprof(1:z1(1))=a1*log(z(1:z1(1))+1);

        c2=V0;
        b2=-0.05;
        a2=(v1-c2)/exp(b2*p1);
        vprof(z1(1)+1:N)=a2*exp(b2*z(z1(1)+1:N))+c2;

print *, '-----v8prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'v1', v1
print *, '----------------'

end subroutine v8prof

subroutine w1prof(H, HWR, z, V0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: w1, w2, w3, w4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H*(1.222*exp(-3.644*HWR)+0.02068);
        z1=minloc(abs(z-p1));
        w1=V0*(-0.0484*exp(-0.8592*HWR)+0.00731);
        if (HWR<0.5) w1=V0*(-0.05433*HWR+0.003);
        b1=-p1;c1=w1;a1=-c1/b1**2;
        wprof(1:z1(1))=0.7*(a1*(z(1:z1(1))+b1)**2+c1);
        
        p2=H;
        z2=minloc(abs(z-p2));
        w2=V0*(-1.034*exp(-0.0328*HWR)+1.025);
        if (HWR<0.5) w2=V0*(-1.034*exp(-0.0328*0.5)+1.025);
        b2=-p2;c2=w2;a2=(w1-c2)/(p1+b2)**2;
        wprof(z1(1)+1:z2(1))=0.7*(a2*(z(z1(1)+1:z2(1))+b2)**2+c2);
        
        wprof(z2(1)+1:N)=w2*0.7;

print *, '-----w1prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'w1', w1
print *, 'w2', w2
print *, '----------------'
end subroutine w1prof

subroutine w2prof(H, HWR, z, V0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: w1, w2, w3, w4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=0.8*H;
        z1=minloc(abs(z-p1));
        w1=V0*(-0.06537*exp(-2.274*HWR)+0.04656);
        b1=-p1;c1=w1;a1=-c1/b1**2;
        wprof(1:z1(1))=0.7*(a1*(z(1:z1(1))+b1)**2+c1);
        
        c2=0;
        b2=-0.05;
        a2=(w1-c2)/exp(b2*p1);
        wprof(z1(1)+1:N)=0.7*(a2*exp(b2*z(z1(1)+1:N))+c2);

print *, '-----w2prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'w1', w1
print *, '----------------'
end subroutine w2prof

subroutine w4prof(H, HWR, z, V0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: w1, w2, w3, w4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H*(1.222*exp(-3.644*HWR)+0.02068);
        z1=minloc(abs(z-p1));
        w1=V0*(-0.0484*exp(-0.8592*HWR)+0.00731);
        if (HWR<0.5) w1=V0*(-0.05433*HWR+0.003);
        b1=-p1;c1=w1;a1=-c1/b1**2;
        wprof(1:z1(1))=a1*(z(1:z1(1))+b1)**2+c1;
        
        p2=H;
        z2=minloc(abs(z-p2));
        w2=V0*(-1.034*exp(-0.0328*HWR)+1.025);
        if (HWR<0.5) w2=V0*(-1.034*exp(-0.0328*0.5)+1.025);
        b2=-p2;c2=w2;a2=(w1-c2)/(p1+b2)**2;
        wprof(z1(1)+1:z2(1))=a2*(z(z1(1)+1:z2(1))+b2)**2+c2;
        
        wprof(z2(1)+1:N)=w2;

print *, '-----w4prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'w1', w1
print *, 'w2', w2
print *, '----------------'
end subroutine w4prof

subroutine w5prof(H, HWR, z, V0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: w1, w2, w3, w4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=0.8*H;
        z1=minloc(abs(z-p1));
        w1=V0*(-0.06537*exp(-2.274*HWR)+0.04656);
        b1=-p1;c1=w1;a1=-c1/b1**2;
        wprof(1:z1(1))=a1*(z(1:z1(1))+b1)**2+c1;
        
        c2=0;
        b2=-0.05;
        a2=(w1-c2)/exp(b2*p1);
        wprof(z1(1)+1:n)=a2*exp(b2*z(z1(1)+1:n))+c2;

print *, '-----w5prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'w1', w1
print *, '----------------'

end subroutine w5prof

subroutine w7prof(H, HWR, z, V0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: w1, w2, w3, w4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H;
        z1=minloc(abs(z-p1));
        w1=V0*(-0.1621*exp(-1.648*HWR)+0.156);
        b1=-p1;c1=w1;a1=-c1/b1**2;
        wprof(1:z1(1))=a1*(z(1:z1(1))+b1)**2+c1;

        c2=0;
        b2=-0.05;
        a2=(w1-c2)/exp(b2*p1);
        wprof(z1(1)+1:N)=a2*exp(b2*z(z1(1)+1:N))+c2;

print *, '-----w7prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'w1', w1
print *, '----------------'

end subroutine w7prof

subroutine w8prof(H, HWR, z, V0, wprof, N)
    implicit none

    real, intent(in ) :: H, HWR, V0
    integer, intent(in ) :: N
    real, intent(inout), dimension(N) :: z, wprof

    integer :: i

!===================================
!           inner vars
!===================================
    real :: p1, p2, p3
    integer :: z1(1), z2(1), z3(1) ! Endpoints
    real :: w1, w2, w3, w4
    real :: a1, a2, a3, a4,               &
    &       b1, b2, b3, b4,               &
    &       c1, c2, c3, c4

        p1=H;
        z1=minloc(abs(z-p1));
        w1=V0*(-0.156*exp(-3.475*HWR)+0.09242);
        
        b1=-p1;c1=w1;a1=-c1/b1**2;
        wprof(1:z1(1))=a1*(z(1:z1(1))+b1)**2+c1;
        
        c2=0;
        b2=-0.05;
        a2=(w1-c2)/exp(b2*p1);
        wprof(z1(1)+1:N)=a2*exp(b2*z(z1(1)+1:N))+c2;

print *, '-----w8prof-----'
print *, 'HGT', H, 'HWR', HWR
print *, 'u0', v0
print *, 'w1', w1
print *, '----------------'

end subroutine w8prof

END MODULE module_sf_ahs_para
