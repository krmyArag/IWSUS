MODULE module_sf_ahs

USE module_sf_ahs_perp, ONLY: ahs_perp ! perpendicular
USE module_sf_ahs_para, ONLY: ahs_para ! parallel
USE module_sf_ahs_obl, ONLY: ahs_obl ! oblique

CONTAINS

!=====================================================================================
!           AHS wind profile of urban canopy for WRF
!           Original Version: 2022/05/11 by ARAG
!=====================================================================================
!    Perpendicular
!=====================================================================================
!
! Input Data from WRF:
!
! Output Data to WRF:
! u_{POS}_AHS, v_{POS}_AHS, w_{POS}_AHS: u, v, w component at typical positions
! POS refers to the absolute positions as following:
!
!                                       ^
!                                     / | \  positive street-orientation
!           -----         -----         |  
!           |   | A  B  C |   |         | 
!           |   |         |   |         |
!           |   | D  O  E |   |         |
!           |   |         |   |         |
!           |   | F  G  H |   |         |
!           -----         -----         |
!                                       |
!                                       
!
!
!======================================================================================
! u_{num}_AHS, v_{num}_AHS, w_{num}_AHS: u, v, w component at relative
! positions, 
! num refers to the relative positions as following:
!
!                                        EXIT 
!           -----         -----
!           |   | 1  2  3 |   |
! Inflow    |   |         |   |
! ----->    |   | 4  5  6 |   |         MIDDLE
!   UA      |   |         |   |
!           |   | 7  8  9 |   |
!           -----         -----
!          leeward       windward
!                                      ENTRANCE
!                    ^
!                  / | \
!                    |
!                    |
!                 Inflow
!                   VA   
!
! The relative positions are determined by wind direction and street
! orientation.
!
!
!======================================================================================
SUBROUTINE ahs(                            &
                 ZA, U0, V0, HWR, HGT_ROOF, Orit, N                            &
                ,uA, uB, uC, uD, uE, uF, uG, uH, uO &
                ,vA, vB, vC, vD, vE, vF, vG, vH, vO &
                ,wA, wB, wC, wD, wE, wF, wG, wH, wO &
                ,U_Ctop, U_canopy, wall, road, roof          &

) 

        implicit none

        !! input
        integer, intent(in) :: N
        real, intent(in)  :: ZA, U0, V0
        real, intent(in)  :: HWR, HGT_ROOF, Orit

        !! output to pkuurban
    real, intent(inout) :: uA, uB, uC, uD, uE, uF, uG, uH, uO
    real, intent(inout) :: vA, vB, vC, vD, vE, vF, vG, vH, vO
    real, intent(inout) :: wA, wB, wC, wD, wE, wF, wG, wH, wO
    real :: U_Ctop, U_canopy, wall, road, roof
        
        !!!!!!!!!!!
    integer :: i
!    
!   wind profile at positions relative to street orientation
!
    real, dimension(N) :: z 
    real, dimension(N) :: u1_AHS, u2_AHS, u3_AHS, u4_AHS, u5_AHS, u6_AHS, u7_AHS, u8_AHS, u9_AHS
    real, dimension(N) :: v1_AHS, v2_AHS, v3_AHS, v4_AHS, v5_AHS, v6_AHS, v7_AHS, v8_AHS, v9_AHS
    real, dimension(N) :: w1_AHS, w2_AHS, w3_AHS, w4_AHS, w5_AHS, w6_AHS, w7_AHS, w8_AHS, w9_AHS
    
    real, dimension(N) :: res_z
!    
!   wind profile at absolute positions
!
    real, dimension(N) :: uA_AHS, uB_AHS, uC_AHS, uD_AHS, uE_AHS, uF_AHS, uG_AHS, uH_AHS, uO_AHS
    real, dimension(N) :: vA_AHS, vB_AHS, vC_AHS, vD_AHS, vE_AHS, vF_AHS, vG_AHS, vH_AHS, vO_AHS
    real, dimension(N) :: wA_AHS, wB_AHS, wC_AHS, wD_AHS, wE_AHS, wF_AHS, wG_AHS, wH_AHS, wO_AHS
    real :: UVDIR, dtheta, UV
    real, parameter :: PI=3.1415926536
    integer :: FLAG_DIR, FLAG_POS

!
!   wind in the RSL ( urban canopy upper than rooftop)
!
    real :: uA_top, uB_top, uC_top, uD_top, uE_top, uF_top, uG_top, uH_top, uO_top
    real :: vA_top, vB_top, vC_top, vD_top, vE_top, vF_top, vG_top, vH_top, vO_top
    real :: wA_top, wB_top, wC_top, wD_top, wE_top, wF_top, wG_top, wH_top, wO_top

!
!   wind speed within street canyon
!
    real :: uA_C, uB_C, uC_C, uD_C, uE_C, uF_C, uG_C, uH_C, uO_C



! initial z in urban canopy
    print *, 'top of urban canopy: ', ZA
    do i=1,N
!        print *, 'gerating z'
        z(i)=real(i)*real(ZA)/real(N)
    end do

print *, 'z1', z(1), 'z50', z(50), 'zMAX', z(N)
! 1. Determine relative wind direction to street orientation
print *, 'U0 =', U0, 'V0 =',V0
UVDIR=atan2(U0,V0)*180./PI
print *, 'UVDIR:', UVDIR
UV=sqrt(U0**2+V0**2)
print *, 'street orietation: ', Orit
print *, 'UV:', UV
dtheta=UVDIR-Orit
print *, 'UVDIR-Orit:', dtheta

FLAG_DIR=3 ! Oblique
if ((dtheta .gt. -22.5) .and. (dtheta .le. 22.5)) FLAG_DIR=1 !Parallel
if ((dtheta .gt. 180-22.5) .and. (dtheta .le. 180+22.5)) FLAG_DIR=1 !Parallel
if ((dtheta .gt. 90-22.5) .and. (dtheta .le. 90+22.5)) FLAG_DIR=2 !Perpendicular
if ((dtheta .gt. -90-22.5) .and. (dtheta .le. -90+22.5)) FLAG_DIR=2 !Perpendicular
print *, 'FLAG_DIR', FLAG_DIR

! 2. Caculate wind profile at each positions
select case (FLAG_DIR)
case (1)
 call ahs_para( HGT_ROOF, HWR, z, UV*cos(dtheta), N,                &   
                     u1_AHS, u2_AHS, u3_AHS, u4_AHS, u5_AHS, u6_AHS, u7_AHS, u8_AHS, u9_AHS,       &
                     v1_AHS, v2_AHS, v3_AHS, v4_AHS, v5_AHS, v6_AHS, v7_AHS, v8_AHS, v9_AHS,       &
                     w1_AHS, w2_AHS, w3_AHS, w4_AHS, w5_AHS, w6_AHS, w7_AHS, w8_AHS, w9_AHS       &
                     )
 !call AHS_para()
print *, 'parallel'
case (2)
 call ahs_perp( HGT_ROOF, HWR, z, UV*sin(dtheta), N,                &   
                     u1_AHS, u2_AHS, u3_AHS, u4_AHS, u5_AHS, u6_AHS, u7_AHS, u8_AHS, u9_AHS,       &
                     v1_AHS, v2_AHS, v3_AHS, v4_AHS, v5_AHS, v6_AHS, v7_AHS, v8_AHS, v9_AHS,       &
                     w1_AHS, w2_AHS, w3_AHS, w4_AHS, w5_AHS, w6_AHS, w7_AHS, w8_AHS, w9_AHS       &
                     )
print *, 'perpendicular'
case default
 call ahs_obl( HGT_ROOF, HWR, z, UV*sin(dtheta), UV*cos(dtheta), N,                &   
                     u1_AHS, u2_AHS, u3_AHS, u4_AHS, u5_AHS, u6_AHS, u7_AHS, u8_AHS, u9_AHS,       &
                     v1_AHS, v2_AHS, v3_AHS, v4_AHS, v5_AHS, v6_AHS, v7_AHS, v8_AHS, v9_AHS,       &
                     w1_AHS, w2_AHS, w3_AHS, w4_AHS, w5_AHS, w6_AHS, w7_AHS, w8_AHS, w9_AHS       &
                     )
! call AHS_oblique()
print *, 'oblique'
end select 

! 3. Value assignment at absolute positions (num -> POS)
FLAG_POS=1 !
if ((dtheta .gt. -22.5) .and. (dtheta .le. 90+22.5)) FLAG_POS=1 ! Positive
if ((dtheta .gt. 90+22.5) .and. (dtheta .le. 180+22.5)) FLAG_POS=2 ! the literal street-orientation symmetry
if ((dtheta .gt. 180+22.5) .and. (dtheta .le. 270+22.5)) FLAG_POS=3 ! central symmetry
if ((dtheta .gt. -180+22.5) .and. (dtheta .le. -90+22.5)) FLAG_POS=3 ! central symmetry
if ((dtheta .gt. -90+22.5) .and. (dtheta .le. -22.5)) FLAG_POS=4 ! street-orientation symmetry
print *, 'FLAG_POS', FLAG_POS
    
select case(FLAG_POS)
case (1)
    uA_AHS = u1_AHS; uB_AHS = u2_AHS; uC_AHS = u3_AHS
    uD_AHS = u4_AHS; uO_AHS = u5_AHS; uE_AHS = u6_AHS
    uF_AHS = u7_AHS; uG_AHS = u8_AHS; uH_AHS = u9_AHS
    vA_AHS = v1_AHS; vB_AHS = v2_AHS; vC_AHS = v3_AHS
    vD_AHS = v4_AHS; vO_AHS = v5_AHS; vE_AHS = v6_AHS
    vF_AHS = v7_AHS; vG_AHS = v8_AHS; vH_AHS = v9_AHS
    wA_AHS = w1_AHS; wB_AHS = w2_AHS; wC_AHS = w3_AHS
    wD_AHS = w4_AHS; wO_AHS = w5_AHS; wE_AHS = w6_AHS
    wF_AHS = w7_AHS; wG_AHS = w8_AHS; wH_AHS = w9_AHS
case (2)
    uA_AHS = u7_AHS; uB_AHS = u8_AHS; uC_AHS = u9_AHS
    uD_AHS = u4_AHS; uO_AHS = u5_AHS; uE_AHS = u6_AHS
    uF_AHS = u1_AHS; uG_AHS = u2_AHS; uH_AHS = u3_AHS
    vA_AHS = v7_AHS; vB_AHS = v8_AHS; vC_AHS = v9_AHS
    vD_AHS = v4_AHS; vO_AHS = v5_AHS; vE_AHS = v6_AHS
    vF_AHS = v1_AHS; vG_AHS = v2_AHS; vH_AHS = v3_AHS
    wA_AHS = w7_AHS; wB_AHS = w8_AHS; wC_AHS = w9_AHS
    wD_AHS = w4_AHS; wO_AHS = w5_AHS; wE_AHS = w6_AHS
    wF_AHS = w1_AHS; wG_AHS = w2_AHS; wH_AHS = w3_AHS
case (3)
    uA_AHS = u9_AHS; uB_AHS = u8_AHS; uC_AHS = u7_AHS
    uD_AHS = u6_AHS; uO_AHS = u5_AHS; uE_AHS = u4_AHS
    uF_AHS = u3_AHS; uG_AHS = u2_AHS; uH_AHS = u1_AHS
    vA_AHS = v9_AHS; vB_AHS = v8_AHS; vC_AHS = v7_AHS
    vD_AHS = v6_AHS; vO_AHS = v5_AHS; vE_AHS = v4_AHS
    vF_AHS = v3_AHS; vG_AHS = v2_AHS; vH_AHS = v1_AHS
    wA_AHS = w9_AHS; wB_AHS = w8_AHS; wC_AHS = w7_AHS
    wD_AHS = w6_AHS; wO_AHS = w5_AHS; wE_AHS = w4_AHS
    wF_AHS = w3_AHS; wG_AHS = w2_AHS; wH_AHS = w1_AHS
case (4)
    uA_AHS = u3_AHS; uB_AHS = u2_AHS; uC_AHS = u1_AHS
    uD_AHS = u6_AHS; uO_AHS = u5_AHS; uE_AHS = u4_AHS
    uF_AHS = u9_AHS; uG_AHS = u8_AHS; uH_AHS = u7_AHS
    vA_AHS = v3_AHS; vB_AHS = v2_AHS; vC_AHS = v1_AHS
    vD_AHS = v6_AHS; vO_AHS = v5_AHS; vE_AHS = v4_AHS
    vF_AHS = v9_AHS; vG_AHS = v8_AHS; vH_AHS = v7_AHS
    wA_AHS = w3_AHS; wB_AHS = w2_AHS; wC_AHS = w1_AHS
    wD_AHS = w6_AHS; wO_AHS = w5_AHS; wE_AHS = w4_AHS
    wF_AHS = w9_AHS; wG_AHS = w8_AHS; wH_AHS = w7_AHS
case default
    uA_AHS = u1_AHS; uB_AHS = u2_AHS; uC_AHS = u3_AHS
    uD_AHS = u4_AHS; uO_AHS = u5_AHS; uE_AHS = u6_AHS
    uF_AHS = u7_AHS; uG_AHS = u8_AHS; uH_AHS = u9_AHS
    vA_AHS = v1_AHS; vB_AHS = v2_AHS; vC_AHS = v3_AHS
    vD_AHS = v4_AHS; vO_AHS = v5_AHS; vE_AHS = v6_AHS
    vF_AHS = v7_AHS; vG_AHS = v8_AHS; vH_AHS = v9_AHS
    wA_AHS = w1_AHS; wB_AHS = w2_AHS; wC_AHS = w3_AHS
    wD_AHS = w4_AHS; wO_AHS = w5_AHS; wE_AHS = w6_AHS
    wF_AHS = w7_AHS; wG_AHS = w8_AHS; wH_AHS = w9_AHS
    print *, 'Warning: FLAG_POS'
end select

print *, 'Utop'
!print *, 'Atop'
    call UtopUcan(uA_AHS, vA_AHS, wA_AHS, z, N, uA_top, uA_C, HGT_ROOF, HWR, ZA &
                   ,uA, vA, wA )
print *, 'uA', uA_top 
!print *, 'Btop'
    call UtopUcan(uB_AHS, vB_AHS, wB_AHS, z, N, uB_top, uB_C, HGT_ROOF, HWR, ZA &
                   ,uB, vB, wB )
print *, 'uB', uB_top 
!print *, 'Ctop'
    call UtopUcan(uC_AHS, vC_AHS, wC_AHS, z, N, uC_top, uC_C, HGT_ROOF, HWR, ZA &
                   ,uC, vC, wC )
print *, 'uC', uC_top 
!print *, 'Dtop'
    call UtopUcan(uD_AHS, vD_AHS, wD_AHS, z, N, uD_top, uD_C, HGT_ROOF, HWR, ZA &
                   ,uD, vD, wD )
print *, 'uD', uD_top 
!print *, 'Etop'
    call UtopUcan(uE_AHS, vE_AHS, wE_AHS, z, N, uE_top, uE_C, HGT_ROOF, HWR, ZA &
                   ,uE, vE, wE )
print *, 'uE', uE_top 
!print *, 'Ftop'
    call UtopUcan(uF_AHS, vF_AHS, wF_AHS, z, N, uF_top, uF_C, HGT_ROOF, HWR, ZA &
                   ,uF,  vF, wF )
print *, 'uF', uF_top 
!print *, 'Gtop'
    call UtopUcan(uG_AHS, vG_AHS, wG_AHS, z, N, uG_top, uG_C, HGT_ROOF, HWR, ZA &
                   ,uG, vG, wG )
print *, 'uG', uG_top 
!print *, 'Htop'
    call UtopUcan(uH_AHS, vH_AHS, wH_AHS, z, N, uH_top, uH_C, HGT_ROOF, HWR, ZA &
                   ,uH, vH, wH )
print *, 'uH', uH_top 
!print *, 'Otop'
    call UtopUcan(uO_AHS, vO_AHS, wO_AHS, z, N, uO_top, uO_C, HGT_ROOF, HWR, ZA &
                   ,uO, vO, wO )
print *, 'uO', uO_top 
    U_Ctop=(uA_top + uB_top + uC_top + uD_top + uE_top + uF_top + uG_top + uH_top + uO_top) / 9
    U_canopy=(uA_C + uB_C + uC_C + uD_C + uE_C + uF_C + uG_C + uH_C + uO_C) / 9
    print *, "Utop", U_Ctop
    print *, "U_canopy", U_canopy
    
    call resistance( N, HGT_ROOF                      &
         , u1_AHS, u2_AHS, u3_AHS, u4_AHS, u5_AHS, u6_AHS, u7_AHS, u8_AHS, u9_AHS     &
         , v1_AHS, v2_AHS, v3_AHS, v4_AHS, v5_AHS, v6_AHS, v7_AHS, v8_AHS, v9_AHS     &
         , w1_AHS, w2_AHS, w3_AHS, w4_AHS, w5_AHS, w6_AHS, w7_AHS, w8_AHS, w9_AHS     &
         , z, wall, road, roof                &
)
    print *, 'RESISTANCE'
    print *, 'wall', wall, 'road', road, 'roof', roof


END SUBROUTINE ahs

subroutine UtopUcan(                           &
                    uprof, vprof, wprof, z, N, &
                    U_Ctop, U_canopy, H, HWR, ZA,   &
                    Ucan, Vcan, Wcan  &  ! average wind components within street
                    )
implicit none

integer, intent(in) :: N
real, intent(in), dimension(N) :: z, uprof, vprof, wprof
real, intent(in) :: H, HWR, ZA
real, intent(inout) :: U_Ctop, U_canopy ! average wind speed uvw of the street
real, intent(inout) :: Ucan, Vcan, Wcan

! inner vars
integer :: ztop(1) 
integer :: i
real :: dz
real :: Utop, Vtop, Wtop

! top
Utop=0
Vtop=0
Wtop=0
ztop=minloc(abs(z-H))
dz=z(2)-z(1)
print *, 'dz=',dz,'H=',H,'ztop=',ztop(1)
!print *, 'utop'
!print *, uprof
!print *, 'vtop'
!print *, vprof
!print *, 'vtop'
!print *, uprof
if (ZA .le. H) print *, '!!!!WARNING!!!! building height out of urban canopy top!'
do i=ztop(1),N-1
    Utop=Utop+(uprof(i)+uprof(i+1))*0.5*dz/abs(ZA-H)
    Vtop=Vtop+(vprof(i)+vprof(i+1))*0.5*dz/abs(ZA-H)
    Wtop=Wtop+(wprof(i)+wprof(i+1))*0.5*dz/abs(ZA-H)
enddo
U_Ctop=sqrt(Utop**2.0 + Vtop**2.0 + Wtop**2.0)

! Ucan
Ucan=0
Vcan=0
Wcan=0
do i=1,ztop(1)
    Ucan=Ucan+(uprof(i)+uprof(i+1))*0.5*dz/H
    Vcan=Vcan+(vprof(i)+vprof(i+1))*0.5*dz/H
    Wcan=Wcan+(wprof(i)+wprof(i+1))*0.5*dz/H
enddo
U_canopy=sqrt(Ucan**2.0 + Vcan ** 2.0 + Wcan**2.0)

end subroutine UtopUcan

SUBROUTINE resistance( N, HGT                      &
         , u1, u2, u3, u4, u5, u6, u7, u8, u9     &
         , v1, v2, v3, v4, v5, v6, v7, v8, v9     &
         , w1, w2, w3, w4, w5, w6, w7, w8, w9     &
         , z, wall, road, roof                &        
)

implicit none

integer, intent(in) :: N
real, intent(in) :: HGT

real, intent(in), dimension(N) :: z
real, intent(in), dimension(N) :: u1, u2, u3, u4, u5, u6, u7, u8, u9
real, intent(in), dimension(N) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
real, intent(in), dimension(N) :: w1, w2, w3, w4, w5, w6, w7, w8, w9

real, intent(inout) :: wall, roof, road

integer :: i, ztop(1), j
real, dimension(N,9) :: resWall, resRoad, resRoof, uu, vv, ww
real :: dz

uu(:,1)=u1
uu(:,2)=u2
uu(:,3)=u3
uu(:,4)=u4
uu(:,5)=u5
uu(:,6)=u6
uu(:,7)=u7
uu(:,8)=u8
uu(:,9)=u9
vv(:,1)=v1
vv(:,2)=v2
vv(:,3)=u3
vv(:,4)=u4
vv(:,5)=v5
vv(:,6)=v6
vv(:,7)=v7
vv(:,8)=v8
vv(:,9)=v9
ww(:,1)=w1
ww(:,2)=w2
ww(:,3)=w3
ww(:,4)=w4
ww(:,5)=w5
ww(:,6)=w6
ww(:,7)=w7
ww(:,8)=w8
ww(:,9)=w9
print *, 're-allocated uu vv ww'
resWall=1./(11.8+4.2*sqrt(uu**2+vv**2))
resRoad=1./(11.8+4.2*sqrt(uu**2+vv**2+ww**2))
resRoof=1./(11.8+4.2*sqrt(uu**2+vv**2+ww**2))
print *, 'resistance matrix finished'

ztop = minloc(abs(z-HGT))
dz=z(2)-z(1)
! RES for wall
roof=0;road=0;wall=0
do j=1,9
    DO i=ztop(1),N-1
        roof=roof+(resRoof(i,j)+resRoof(i+1,j))*dz*0.5/HGT/9
    ENDDO
    DO i=1,ztop(1)-1
        road=road+(resRoad(i,j)+resRoad(i+1,j))*dz*0.5/HGT/9
        wall=roof+(resWall(i,j)+resWall(i+1,j))*dz*0.5/HGT/9
        
    ENDDO
enddo
print *, 'roof', roof, 'road', road, 'wall', wall


END SUBROUTINE resistance

END MODULE module_sf_ahs
