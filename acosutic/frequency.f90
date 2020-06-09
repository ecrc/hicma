!*************************************************************************
!program Acoustic_code_hardCase

subroutine Acoustic_code_hardCase(nip, ntrain, zzsparse_temp, crhs) bind(C, name="AL4SAN_Matrix_Generation")
!====================================================================
  use, intrinsic :: iso_c_binding
  use global_com 
  use global_dim
  use gaussintegraldata
  !implicit none
  real(kind=dp)::mem_est
  integer(c_int) ::ntrain
  integer(c_int) ::nip
  complex(c_double_complex), dimension(ntrain*nip, ntrain*nip):: zzsparse_temp
  complex(c_double_complex), dimension(ntrain*nip):: crhs
!---------------------------------------------------------------------

if(nip.eq.3) then
open(unit=12,file='../hicma-dev/acosutic/mom_1.inp',status='old') ! main input file
elseif(nip.eq.6) then
open(unit=12,file='../hicma-dev/acosutic/mom_2.inp',status='old') ! main input file
elseif(nip.eq.12) then
open(unit=12,file='../hicma-dev/acosutic/mom_3.inp',status='old') ! main input file
else
     print*,"nippp should only be 1,2 or 3!!"
     stop
  endif
  call inmom_new   ! read in the control parameters for the MOM algorithm
!---------------------------------------------------------------------
 IF ((defile .eq. 0) .AND.  ( ntrain.gt. 1).AND.( ntrain.le. 120)) THEN !read ANSYS file
     open(unit=111,file='../hicma-dev/acosutic/geo_curve_tri_120.inp')  ! mesh file for curve surface triangles
 elseif ((defile .eq. 0) .AND.  ( ntrain.gt.120 ).AND.( ntrain.le. 396)) THEN !read ANSYS file  
     open(unit=111,file='../hicma-dev/acosutic/geo_curve_tri_396.inp')  ! mesh file for curve surface triangles
elseif ((defile .eq. 0) .AND.  ( ntrain.gt.396 ).AND.( ntrain.le. 622)) THEN !read ANSYS file
     open(unit=111,file='../hicma-dev/acosutic/geo_curve_tri_622.inp')  ! mesh file for curve surface triangles
elseif ((defile .eq. 0) .AND.  ( ntrain.gt.622 ).AND.( ntrain.le. 1126))THEN !read ANSYS file
     open(unit=111,file='../hicma-dev/acosutic/geo_curve_tri_1126.inp')  ! mesh file for curve surface triangles
elseif ((defile .eq. 0) .AND.  ( ntrain.gt.1126 ).AND.( ntrain.le. 2560))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic/geo_curve_tri_2560.inp')  ! mesh file for curve surface triangles
elseif ((defile .eq. 0) .AND.  ( ntrain.gt.2560 ).AND.( ntrain.le. 4032))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic/geo_curve_tri_4032.inp')  ! mesh file for curve surface triangles
elseif ((defile .eq. 0) .AND.  ( ntrain.gt.4032 ).AND.( ntrain.le. 7320))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic/geo_curve_tri_7320.inp')  ! mesh file for curve surface triangles
elseif ((defile .eq. 0) .AND.  ( ntrain.gt.7320 ).AND.( ntrain.le. 12288))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic/geo_curve_tri_12288.inp')  ! mesh file for curve surface triangles
else
     print*,"in case there is a need for a mesh larger than 12288 !!"
!open(unit=111,file='../hicma-dev/acosutic/geo_curve_tri_120.inp')  ! mesh file for curve surface triangles
  !   stop
  endif
call incurvefile ! read curve surface discretization info.
!---------------------------------------------------------------------
    call setgaussintpara(nipp,nipp)  !set gauss quadrature points for interpolation points
    call setgaussintpara11(nhdgqp_far) !set higher degree gauss quadrature points for source triangle_far
    call setgaussintpara22(nhdgqp_near) !set higher degree gauss quadrature points for source triangle_near
    call setgaussintline_new(nlqp,nclgp) !set line quadrature points for duffy line integral
!---------------------------------------------------------------------

!---------------------------------------------------------------------
  ntriangle=ntrain
  print*,"number of triangles:",ntriangle
  print*, 'NUMBER OF UNKNOWNS:', ntriangle*nipp
  allocate(zzsparse(1:(ntriangle*nipp)**2))   !z matrix
  allocate(JVector(1:ntriangle*nipp))   !1/J matrix
  zzsparse=(0.0_dp,0.0_dp)
  JVector=0.0_dp
  
  mem_est=(ntriangle*nipp)**2*complex_mem/1024.0_dp/1024.0_dp + (ntriangle*nipp)*complex_mem/1024.0_dp/1024.0_dp
  print*,'Near field matrices require (MB): ',mem_est
!---------------------------------------------------------------------

!---------------------------------------------------------------------

  call hard_nys_oldBasis_SST   !matrix filling
  
  call solve_nys_hard(zzsparse_temp, ntrain, nip, crhs)   !solve the equation and get solution
  
  !call near_sca_hard_near_scheme   !get near field at arbitrary position use near scheme
  
!---------------------------------------------------------------------

  print*,"the whole program is done!"
  print*, 'NUMBER OF UNKNOWNS:', ntriangle*nipp
  
end subroutine Acoustic_code_hardCase
!*************************************************************************

!*************************************************************************
subroutine inmom_new                                           
  use global_com 
  implicit none
  real(kind=dp)::the,phi,ethetad,ephid
  real(kind=dp)::epsilon_r,mu_r
  !integer::nhd_near


  pid=4.0_dp*atan(1.0_dp)  !pi=3.1415926

! read excitation information
  read(12,*) no_freqs
  read(12,*) fre_single  !use different mesh file at different frequency
  read(12,*) the,phi,ethetad,ephid  
  !incident direction of wave and polarization direction of electric field 
  ! direction that the wave is propagating towards
  ! magnitude of theta and phi component of electric field 
  the=the*pid/180.0_dp       
  phi=phi*pid/180.0_dp                     

  !************************************************
  vkk(1)=sin(the)*cos(phi)
  vkk(2)=sin(the)*sin(phi)
  vkk(3)=cos(the)

  vee(1)=cos(the)*cos(phi)*ethetad-sin(phi)*ephid
  vee(2)=cos(the)*sin(phi)*ethetad+cos(phi)*ephid
  vee(3)=-sin(the)*ethetad

  call real_cross_product(vkk,vee,vhh)  !vhh = H field * eta

  print*,'travelling to',vkk
  print*,'Efield components:',vee
  print*,'Efield magnitude:',sqrt(dot_product(vee,vee))
  print*,'Hfield components * eta:',vhh
  print*,'frequency (Hz):',fre_single
  !************************************************

  read(12,*) the_sta_degree,the_end_degree      !numbers of scanned theta and phi degrees
  read(12,*) phi_sta_degree,phi_end_degree      ! read frequency (Hz), factor2
  read(12,*) nthe_rcs,nphi_rcs                  ! read rcs and frequency domain current frequencies 
     
 

  the_sta=the_sta_degree*pid/180.0_dp
  the_end=the_end_degree*pid/180.0_dp
  phi_sta=phi_sta_degree*pid/180.0_dp
  phi_end=phi_end_degree*pid/180.0_dp

  read(12,*) epsilon_r,mu_r,loss_sigma  !relative parameter of background medium

  if (loss_sigma/=0.0_dp) then
     print*,'lossy part is incomplete... correct the fielde_dbl and aim.f90'
     stop
  end if

  cd=340.29_dp                        !sound speed of free space
  eps0d=epsilon_r/(pid*4.d-7*cd**2)   !permittivity of background medium
  rmu0d=pid*4.d-7*mu_r                !permeability of background medium
!  cd=cd/sqrt(epsilon_r*mu_r)          !light speed of background medium
  vwavelength=cd/fre_single           !wavelength of the single frequency
  wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
  cjvk=cjj*wavenumkk                  !j*k
  veta=sqrt(rmu0d/eps0d)              !wave impedance of background medium

  print*,"sound speed (m/s):",cd
!  print*,"permeability of background medium:",rmu0d
!  print*,"permittivity of background medium:",eps0d
!  print*,"wave impedance of background medium:",veta
  print*,"wavelength:",vwavelength
  print*,"wave number k:",wavenumkk



  read(12,*) order_tri
  if(order_tri.eq.1) then
     nnpt=3
     print*,"-------------------- 3 nodes for mesh plat triangle -------------------"
  elseif(order_tri.eq.2) then
     nnpt=6
     print*,"-------------------- 6 nodes for mesh curve triangle ------------------"
  else
     print*,"can not set higher order mesh now"
     stop
  endif

  
    read(12,*) order_old
    
    
    if(order_old.eq.0) then
        nipp=1
        print*,"---------------------- 1 points for interpolation ---------------------"
        
    elseif(order_old.eq.1) then
        nipp=3
        print*,"---------------------- 3 points for interpolation ---------------------"
        
    elseif(order_old.eq.2) then
        nipp=6
        print*,"---------------------- 6 points for interpolation ---------------------"        
        
    elseif(order_old.eq.3) then
        nipp=12
        print*,"---------------------- 12 points for interpolation ---------------------"
        
      
    else
        print*,"can not set higher order basis function now"
        stop
    endif
  
    print*,"nipp = ",nipp
    
  read(12,*) nhdgqp_far

  read(12,*) nhdgqp_near

  read(12,*) nlqp

  read(12,*) nclgp

  read(12,*) precis

  read(12,*) nitermax

  read(12,*) nearpd

  read(12,*) ngprcs

  read(12,*) denmrcs

  read(12,*) defile

  read(12,*) depconvert

  read(12,*) dejjmie

  read(12,*) deonep

  read(12,*) detotalp

  read(12,*) rnfobp

  read(12,*) the_sta_degree_nf,the_end_degree_nf

  read(12,*) phi_sta_degree_nf,phi_end_degree_nf

  read(12,*) nthe_nf,nphi_nf

  read(12,*) nusolver



  print*,"read inmom is done!"

  return
end subroutine inmom_new
!*************************************************************************

!*************************************************************************
subroutine incurvefile
!
!  read mesh file for curve surfaces 
!
  use global_geom,only:sunod,nsupan
  use global_com,only:dp,nnpt,aminx,aminy,aminz,amaxx,amaxy,amaxz,ntriangle,nnode

  implicit none
  integer::nnod,npat,i,j,pp,qq,tri,n,nod
  
  if(nnpt.eq.6) then

     do i=1,3
        read(111,*)         !ignore 3 irrelevant lines
     enddo
     read(111,*) npat,nnod      !npat=120; nnod=242
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!change the number of triangles 
     ntriangle=npat
    ! ntriangle=6
     nnode=nnod

     allocate(sunod(3,nnod))     !point coordinates of nodes
     allocate(nsupan(6,npat))    !serial numbers of triangle points. 
     sunod(:,:)=0.0_dp
     nsupan(:,:)=0

     pp=mod(npat,20)
     qq=(npat-pp)/20

     do i = 1 , qq
        read(111,*)             !ignore this and next 2 irrelevant lines
        read(111,*)
        read(111,*)
        do j=1,20
           read(111,*) tri,n,n,n,n,n,nsupan(1:3,tri),nsupan(4:5,tri),nsupan(6,tri)
        enddo
     enddo

     if (pp.ne.0)then 
        read(111,*)
        read(111,*)
        read(111,*)
        do i=1,pp
           read(111,*) tri,n,n,n,n,n,nsupan(1:3,tri),nsupan(4:5,tri),nsupan(6,tri)
        enddo
     endif

     pp=mod(nnod,20)
     qq=(nnod-pp)/20

     do  i=1,3
        read(111,*)
     enddo

     do i = 1 , qq
        read(111,*)	
        read(111,*)	
        do j=1,20
           read(111,*) nod,sunod(1:3,nod)
        enddo
     enddo

     if(pp.ne.0)then
        read(111,*)
        read(111,*)	
        do i=1,pp
           read(111,*) nod,sunod(1:3,nod)
        enddo
     endif

     close(111)


     aminx  = minval(sunod(1,1:nnod))
     aminy  = minval(sunod(2,1:nnod))
     aminz  = minval(sunod(3,1:nnod))

     amaxx  = maxval(sunod(1,1:nnod))
     amaxy  = maxval(sunod(2,1:nnod))
     amaxz  = maxval(sunod(3,1:nnod))

     print *,"the size of geometry"
     print *,'x: ',	 aminx,  amaxx
     print *,'y: ',	 aminy,  amaxy
     print *,'z: ',	 aminz,  amaxz

  else
     print*,"can not read higher order mesh now"
     stop
  endif

  print*,"read curve file is done!"

  return
end subroutine incurvefile
!*************************************************************************

!*************************************************************************
subroutine incurvefile_ideas
!
!  read mesh file for curve surfaces 
!
  use global_geom,only:sunod,nsupan
  use global_com,only:dp,nnpt,aminx,aminy,aminz,amaxx,amaxy,amaxz,ntriangle,nnode

  implicit none
  integer::nnod,npat,i,j,pp,qq,tri,n,nod
  
  if(nnpt.eq.6) then

     read(187,*) npat
     read(186,*) nnod

     ntriangle=npat
     nnode=nnod

     allocate(sunod(3,nnod))     !point coordinate
     allocate(nsupan(6,npat))    !serial numbers of triangle points
     sunod(:,:)=0.0_dp
     nsupan(:,:)=0

     do i = 1 , npat
        read(187,*) tri,nsupan(1,tri),nsupan(6,tri),nsupan(3,tri),nsupan(5,tri),nsupan(2,tri),nsupan(4,tri)
     enddo

     do i = 1 , nnod
        read(186,*) nod,sunod(1:3,nod)
     enddo

     close(186)
     close(187)

     aminx  = minval(sunod(1,1:nnod))
     aminy  = minval(sunod(2,1:nnod))
     aminz  = minval(sunod(3,1:nnod))

     amaxx  = maxval(sunod(1,1:nnod))
     amaxy  = maxval(sunod(2,1:nnod))
     amaxz  = maxval(sunod(3,1:nnod))

     print *,"the size of geometry"
     print *,'x: ',	 aminx,  amaxx
     print *,'y: ',	 aminy,  amaxy
     print *,'z: ',	 aminz,  amaxz

  else
     print*,"can not read higher order mesh now"
     stop
  endif

  print*,"read curve file is done!"

  return
end subroutine incurvefile_ideas
!*************************************************************************


!*************************************************************************
subroutine solve_nys_hard (zzsparse_temp, ntrain, nip, crhs)
  use global_com 
  use global_dim!,only:rj,crhs,zzsparse,JVector!,zzsparse_temp
  use global_geom,only:nsupan
  use gaussintegraldata
  implicit none
  integer::i,kk,q,j,k
  integer::ntrain, nip
  complex(c_double_complex), dimension(ntrain*nip, ntrain*nip):: zzsparse_temp
  complex(c_double_complex), dimension(ntrain*nip):: crhs
  integer::nodeo_sn(nnpt)
  real(kind=dp)::vos1(2),vos2(2),vos3(2)
  real(kind=dp)::vcsio(2),vippo(3),vuo(3),vvo(3),jacobo
  !complex(kind=dp)::zzsparse_temp(ntriangle*nipp,ntriangle*nipp)
  !complex(kind=dp)::crhs_temp(ntriangle*nipp)
  real(kind=dp)::ii_temp(2,2)
  integer::IPIV(ntriangle*nipp)
  integer::INFO

  vos1=(/1.0_dp,0.0_dp/)  !right angle triangle vertex points in (u,v) space	
  vos2=(/0.0_dp,1.0_dp/)  !the first number is u value, the second number is v value	
  vos3=(/0.0_dp,0.0_dp/)

  call rhs_nys_OldBasis(crhs, ntrain, nip)   !incident wave
  allocate(rj(ntriangle*nipp))    !solution rj, J_i^u, J_i^v


  !allocate(zzsparse_temp(ntriangle*nipp,ntriangle*nipp))  !z matrix


rj=(0.0_dp,0.0_dp)
  
!---------------------------------------------------------------------
!  if(ntriangle*nipp.gt.nusolver) then
     
     !Update leading diagonal entries according to the Main-text formulation  
 !    DO k=1,ntriangle*nipp
 !       zzsparse((k-1)*ntriangle*nipp+k)= zzsparse((k-1)*ntriangle*nipp+k)-0.5_dp*JVector(k)
  !   ENDDO 
     
   !  call fan_gmres(ntriangle*nipp)     !Iterative solver
     
    ! DEALLOCATE(zzsparse)   !zzsparse no longer needed. So release corresponding resources.
    ! DEALLOCATE(crhs)       !crhs no longer needed. So release corresponding resources.
    ! DEALLOCATE(JVector)    !JVector no longer needed. So release corresponding resources.
     
    ! print*,"Solution via iterative solver is done!"
     
  !else
      
     rj(1:(ntriangle*nipp)) = crhs(1:(ntriangle*nipp))

     !Convert zzsparse to a (ntriangle*nipp)-by-(ntriangle*nipp) square matrix and update leading diagonal entries according to the Main text formulation 
     DO k=1,ntriangle*nipp
        zzsparse_temp(k,1:(ntriangle*nipp))=zzsparse(((k-1)*ntriangle*nipp+1):(k*ntriangle*nipp))
  !      zzsparse_temp(k,K) = zzsparse_temp(k,K)-0.5_dp*JVector(k)
     ENDDO
   
  !PRINT*,zzsparse(1)
  !PRINT*,zzsparse(2)  
     DEALLOCATE(zzsparse)   !zzsparse no longer needed. So release corresponding resources.
     !DEALLOCATE(crhs)       !crhs no longer needed. So release corresponding resources.
     DEALLOCATE(JVector)    !JVector no longer needed. So release corresponding resources.

  !PRINT*,zzsparse(1)
  !PRINT*,zzsparse(2)
   
  !   call ZGESV(ntriangle*nipp,1,zzsparse_temp,ntriangle*nipp,&
  !              IPIV,rj,ntriangle*nipp,INFO)
  !PRINT*, zzsparse(1)
!  PRINT*,zzsparse_temp(1, 1)
!  PRINT*,zzsparse_temp(2, 1)
!  PRINT*,zzsparse_temp(3, 1)
 ! PRINT*,zzsparse_temp(4, 1)
 ! PRINT*,zzsparse_temp(360, 360)  



  PRINT*,zzsparse_temp(1, 1)
  PRINT*,zzsparse_temp(1, 2)
  PRINT*,zzsparse_temp(1, 3) 
  PRINT*,zzsparse_temp(1, 4)
     if(INFO>0) print*, 'Error in LU inversion',INFO
     if(INFO<0) print*, 'Error in LU inputs', -INFO
     
     print*,"direct solve is done!"
  !endif

   !   DEALLOCATE(zzsparse_temp)

  return
end subroutine solve_nys_hard                                                                 
!*************************************************************************







