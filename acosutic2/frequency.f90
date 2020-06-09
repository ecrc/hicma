!*************************************************************************
!program Acoustic_code_hardCase

subroutine Acoustic_code_hardCase(nip, ntrain) bind(C, name="Acoustic_Init")
!====================================================================
  use, intrinsic :: iso_c_binding
  use global_com 
  use global_dim
  use gaussintegraldata
  !implicit none
  real(kind=dp)::mem_est
  integer(c_int) ::ntrain
  integer(c_int) ::nip
  integer::i,j,p,q,k
  !complex(c_double_complex), dimension(ntrain*nip, ntrain*nip):: zzsparse_temp
  !complex(c_double_complex), dimension(ntrain*nip*ntrain*nip):: zzsparse
  !complex(c_double_complex), dimension(ntrain*nip):: JVector
  !complex(c_double_complex), dimension(ntrain*nip):: crhs

if(nip.eq.3) then
open(unit=12,file='../hicma-dev/acosutic2/mom_1.inp',status='old') ! main input file
elseif(nip.eq.6) then
open(unit=12,file='../hicma-dev/acosutic2/mom_2.inp',status='old') ! main input file
elseif(nip.eq.12) then
open(unit=12,file='../hicma-dev/acosutic2/mom_3.inp',status='old') ! main input file
else
     print*,"nippp should only be 1,2 or 3!!"
     stop
  endif
  call inmom_new   ! read in the control parameters for the MOM algorithm
!---------------------------------------------------------------------
!New Files from IDEAS
!---------------------------------------------------------------------
!print*, "Init ideas ntrain:", ntrain
if ((defile .eq. 0).AND. ( ntrain .eq. 4452)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_4452.inp')
    call incurvefile_ideas
    fre_single=1063_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 4800)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_4800.inp')
    call incurvefile_ideas    
    fre_single=1090_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 5104)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_5104.inp')
    call incurvefile_ideas
    fre_single=1134_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 5936)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_5936.inp')
    call incurvefile_ideas
    fre_single=1215_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 7632)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_7632.inp')
    call incurvefile_ideas
    fre_single=1372_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 7992)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_7992.inp')
    call incurvefile_ideas
    fre_single=1417_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 9758)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_9758.inp')
    call incurvefile_ideas
    fre_single=1563_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 10790)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_10790.inp')
    call incurvefile_ideas
    fre_single=1636_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 11460)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_11460.inp')
    call incurvefile_ideas
    fre_single=1701_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 12366)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_12366.inp')
    call incurvefile_ideas
    fre_single=1764_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 14212)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_14212.inp')
    call incurvefile_ideas
    fre_single=1890_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 15228)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_15228.inp')
    call incurvefile_ideas
    fre_single=1978_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 19020)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_19020.inp')
    call incurvefile_ideas
    fre_single=2181_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 19890)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_19890.inp')
    call incurvefile_ideas
    fre_single=2227_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 20232)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_20232.inp')
    call incurvefile_ideas
    fre_single=2268_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 21696)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_21696.inp')
    call incurvefile_ideas
    fre_single=2350_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 23616)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_23616.inp')
    call incurvefile_ideas
    fre_single=2430_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 27642)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_27642.inp')
    call incurvefile_ideas
    fre_single=2633_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 28956)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_28956.inp')
    call incurvefile_ideas
    fre_single=2699_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 29680)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_29680.inp')
    call incurvefile_ideas
    fre_single=2735_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 30702)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_30702.inp')
    call incurvefile_ideas
    fre_single=2771_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 31200)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_31200.inp')
    call incurvefile_ideas
    fre_single=2798_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 31780)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_31780.inp')
    call incurvefile_ideas
    fre_single=2827_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 36240)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_36240.inp')
    call incurvefile_ideas
    fre_single=3038_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 46208)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_46208.inp')
    call incurvefile_ideas
    fre_single=3402_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 51700)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_51700.inp')
    call incurvefile_ideas
    fre_single=3620_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 61920)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_61920.inp')
    call incurvefile_ideas
    fre_single=3953_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .eq. 78916)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_78916.inp')
    call incurvefile_ideas
    fre_single=4477_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 82926)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_82926.inp')
    call incurvefile_ideas
    fre_single=4586_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 98910)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_98910.inp')
    call incurvefile_ideas
    fre_single=5004_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .eq. 115640)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_115640.inp')
    call incurvefile_ideas
    fre_single=5435_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .le. 135548)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_135548.inp')
    call incurvefile_ideas
    fre_single=3858_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .le. 149448)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_149448.inp')
    call incurvefile_ideas
    fre_single=4046_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND. ( ntrain .le. 167706)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_167706.inp')
    call incurvefile_ideas
    fre_single=5044_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND. ( ntrain .le. 283864)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_283864.inp')
    call incurvefile_ideas
    fre_single=5507_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

!---------------------------------------------------------------------
elseif ((defile .eq. 0) .AND. ( ntrain .eq. 3170)) THEN !read ANSYS file
open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_3170.inp')
    call incurvefile
    fre_single=945_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk

elseif ((defile .eq. 0) .AND.( ntrain .eq. 5332)) THEN !read ANSYS file
open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_5332.lis')
    call incurvefile
    fre_single=1215_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 6224)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_6224.lis')
    call incurvefile
    fre_single=1308_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 8738)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_8738.lis')
    call incurvefile
    fre_single=1546_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 9420)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_9420.lis')
    call incurvefile
    fre_single=1605_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 9798)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_9798.lis')
    call incurvefile
    fre_single=1636_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 10616)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_10616.lis')
    call incurvefile
    fre_single=1701_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 11074)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_11074.lis')
    call incurvefile
    fre_single=1736_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 11548)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_11548.lis')
    call incurvefile
    fre_single=1772_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 12046)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_12046.lis')
    call incurvefile
    fre_single=1810_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
    print*,"I am here"   


elseif ((defile .eq. 0) .AND.( ntrain .eq. 13140)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_13140.lis')
    call incurvefile
    fre_single=1890_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 14410)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_14410.lis')
    call incurvefile
    fre_single=1978_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 15108)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_15108.lis')
    call incurvefile
    fre_single=2025_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 15914)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_15914.lis')
    call incurvefile
    fre_single=2074_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 16718)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_16718.lis')
    call incurvefile
    fre_single=2126_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 17598)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_17598.lis')
    call incurvefile
    fre_single=2181_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 18530)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_18530.lis')
    call incurvefile
    fre_single=2238_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 19546)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_19546.lis')
    call incurvefile
    fre_single=2299_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 20726)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_20726.lis')
    call incurvefile
    fre_single=2363_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 21920)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_21920.lis')
    call incurvefile
    fre_single=2430_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 23312)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_23312.lis')
    call incurvefile
    fre_single=2502_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 24742)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_24742.lis')
    call incurvefile
    fre_single=2577_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 26310)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_26310.lis')
    call incurvefile
    fre_single=2658_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 28076)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_28076.lis')
    call incurvefile
    fre_single=2744_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 29998)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_29998.lis')
    call incurvefile
    fre_single=2835_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 32096)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_32096.lis')
    call incurvefile
    fre_single=2933_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 33246)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_33246.lis')
    call incurvefile
    fre_single=2985_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


elseif ((defile .eq. 0) .AND.( ntrain .eq. 34474)) THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_34474.lis')
    call incurvefile
    fre_single=3038_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


!---------------------------------------------------------------------
elseif((defile .eq. 0) .AND. ( ntrain .eq. 120)) THEN !read ANSYS file
     open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_120.inp')  ! mesh file for curve surface triangles
     call incurvefile

elseif ((defile .eq. 0) .AND.( ntrain .eq. 396)) THEN !read ANSYS file
     open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_396.inp')  ! mesh file for curve surface triangles
     call incurvefile
elseif ((defile .eq. 0) .AND.( ntrain .eq. 622)) THEN !read ANSYS file
     open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_622.inp')  ! mesh file for curve surface triangles
    call incurvefile
elseif ((defile .eq. 0) .AND.( ntrain .eq. 1126))THEN !read ANSYS file
     open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_1126.inp')  ! mesh file for curve surface triangles
     call incurvefile
elseif ((defile .eq. 0) .AND. ( ntrain .eq. 1008)) THEN !read ANSYS fil ......chenage
     open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_1008.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
    fre_single=300_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
elseif ((defile .eq. 0) .AND. ( ntrain .eq. 3156)) THEN !read ANSYS fil ......chenage
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_3156.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
elseif ((defile .eq. 0) .AND. ( ntrain .eq. 252)) THEN !read ANSYS fil ......chenage
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_252.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
elseif ((defile .eq. 0) .AND. ( ntrain .eq.2560))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_2560.inp')  ! mesh file for curve surface triangles
    call incurvefile
elseif ((defile .eq. 0) .AND. ( ntrain .eq. 4032))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_4032.inp')  ! mesh file for curve surface triangles
   call incurvefile
elseif ((defile .eq. 0) .AND. ( ntrain .eq. 7320))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_7320.inp')  ! mesh file for curve surface triangles
    call incurvefile
elseif ((defile .eq. 0) .AND. ( ntrain .eq. 12288))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_12288.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
    fre_single=1810_dp
    !fre_single=1797_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
    print*,"I am here 12288"
elseif ((defile .eq. 0) .AND. ( ntrain .eq. 14338))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_14338.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
    fre_single=1978_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
elseif ((defile .eq. 0) .AND. ( ntrain .eq. 16180))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_16180.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
    fre_single=2074_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
elseif ((defile .eq. 0) .AND.  ( ntrain .eq. 22370))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_22370.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
    fre_single=2502_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
elseif ((defile .eq. 0) .AND.  ( ntrain .eq. 41258))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_41258.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
    fre_single=3007_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
elseif ((defile .eq. 0) .AND.  ( ntrain .eq. 49152))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_49152.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
    fre_single=3807_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
elseif ((defile .eq. 0) .AND.  ( ntrain .eq. 60204))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_60204.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
    fre_single=3680_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
elseif ((defile .eq. 0) .AND. ( ntrain .eq. 63380))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_63380.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
elseif ((defile .eq. 0) .AND.  ( ntrain .eq. 93590))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_93590.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
    fre_single=3720_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
elseif ((defile .eq. 0) .AND.  ( ntrain .eq. 196608))THEN !read ANSYS file
    open(unit=111,file='../hicma-dev/acosutic2/geo_curve_tri_196608.inp')  ! mesh file for curve surface triangles
    call incurvefile_ideas
    fre_single=4400_dp
    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk
!elseif (defile .eq. 1) THEN !read submarine file
!    open(unit=186,file='../hicma-dev/acosutic2/geo_node_submarine_nowing_h0.5m.lis')
!    open(unit=187,file='../hicma-dev/acosutic2/geo_surf_submarine_nowing_h0.5m.lis')
!    call incurvefile_ideas_sub
!    print*,"I am on submarine:"
!    open(unit=111,file='../hicma-dev/acosutic2/geo_mesh_Submarine_10818.inp')
!    call incurvefile_ideas
else
     print*," make sure the input file is avavilble !!"
  endif

!call incurvefile ! read curve surface discretization info.

!---------------------------------------------------------------------
    call setgaussintpara(nipp,nipp)  !set gauss quadrature points for interpolation points
    call setgaussintpara11(nhdgqp_far) !set higher degree gauss quadrature points for source triangle_far
    call setgaussintpara22(nhdgqp_near) !set higher degree gauss quadrature points for source triangle_near
    call setgaussintline_new(nlqp,nclgp) !set line quadrature points for duffy line integral
!---------------------------------------------------------------------

!---------------------------------------------------------------------
  ntriangle=ntrain
  !print*,"number of triangles:",ntriangle
  !print*, 'NUMBER OF UNKNOWNS:', ntriangle*nipp
  !allocate(zzsparse(1:(ntriangle*nipp)**2))   !z matrix
  !allocate(JVector(1:ntriangle*nipp))   !1/J matrix
  !zzsparse=(0.0_dp,0.0_dp)
  !JVector=0.0_dp
  
  mem_est=(ntriangle*nipp)**2*complex_mem/1024.0_dp/1024.0_dp + (ntriangle*nipp)*complex_mem/1024.0_dp/1024.0_dp
  !print*,'Near field matrices require (MB): ',mem_est
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!DO q=1,ntriangle
!      DO p=1,ntriangle
!  call hard_nys_oldBasis_SST(nip, ntrain, zzsparse,  p, q, JVector)  !matrix filling
!     ENDDO
! ENDDO 
  !call solve_nys_hard(zzsparse_temp, ntrain, nip, crhs)   !solve the equation and get solution
  
  !call near_sca_hard_near_scheme   !get near field at arbitrary position use near scheme
  
!---------------------------------------------------------------------

  !print*,"the whole program is done!"
  !print*, 'NUMBER OF UNKNOWNS:', ntriangle*nipp
  
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

  !print*,'travelling to',vkk
  !print*,'Efield components:',vee
  !print*,'Efield magnitude:',sqrt(dot_product(vee,vee))
  !print*,'Hfield components * eta:',vhh
  !print*,'frequency (Hz):',fre_single
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

    cd=340.29_dp
    vwavelength=cd/fre_single           !wavelength of the single frequency
    wavenumkk=2.0_dp*pid*fre_single/cd  !wave number k
    cjvk=cjj*wavenumkk


  read(12,*) order_tri
  if(order_tri.eq.1) then
     nnpt=3
   !  print*,"-------------------- 3 nodes for mesh plat triangle -------------------"
  elseif(order_tri.eq.2) then
     nnpt=6
    ! print*,"-------------------- 6 nodes for mesh curve triangle ------------------"
  else
     print*,"can not set higher order mesh now"
     stop
  endif

  
    read(12,*) order_old
    
    
    if(order_old.eq.0) then
        nipp=1
     !   print*,"---------------------- 1 points for interpolation ---------------------"
        
    elseif(order_old.eq.1) then
        nipp=3
      !  print*,"---------------------- 3 points for interpolation ---------------------"
        
    elseif(order_old.eq.2) then
        nipp=6
       ! print*,"---------------------- 6 points for interpolation ---------------------"        
        
    elseif(order_old.eq.3) then
        nipp=12
        !print*,"---------------------- 12 points for interpolation ---------------------"
        
      
    else
        print*,"can not set higher order basis function now"
        stop
    endif
  
    !print*,"nipp = ",nipp
    
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



  !print*,"read inmom is done!"

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

     ntriangle=npat
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

     !print *,"the size of geometry"
     !print *,'x: ',	 aminx,  amaxx
     !print *,'y: ',	 aminy,  amaxy
     !print *,'z: ',	 aminz,  amaxz

  else
     print*,"can not read higher order mesh now"
     stop
  endif

  !print*,"read curve file is done!"

  return
end subroutine incurvefile
!*************************************************************************

!*************************************************************************
subroutine incurvefile_ideas
use global_geom,only:sunod,nsupan
use global_com,only:dp,nnpt,aminx,aminy,aminz,amaxx,amaxy,amaxz,ntriangle,nnode,centerx,centery,centerz

implicit none
integer::nnod,npat,i,j,pp,qq,tri,n,nod
real(kind=dp),dimension(3)::rcent

if(nnpt.eq.6) then
read(111,*) nnod,npat
ntriangle=npat
nnode=nnod

allocate(sunod(3,nnod))     !point coordinate (x,y,z) for each node
allocate(nsupan(6,npat))    !serial numbers of triangle points (the six nodes for each triangle)
sunod(:,:)=0.0_dp
nsupan(:,:)=0


do j = 1 , nnod
read(111,*)sunod(1,j),sunod(2,j),sunod(3,j)
end do


do i = 1 , npat
read(111,*) nsupan(1,i),nsupan(6,i),nsupan(3,i),nsupan(5,i),nsupan(2,i),nsupan(4,i)
!ead(1872,*) nsupan(1,i),nsupan(3,i),nsupan(5,i),nsupan(2,i),nsupan(4,i),nsupan(6,i)
!read(111,*) nsupan(1,i),nsupan(5,i),nsupan(3,i),nsupan(6,i),nsupan(4,i),nsupan(2,i)
enddo

close(111)
amaxx  = maxval(sunod(1,1:nnod))
aminx  = minval(sunod(1,1:nnod))
centerx=(amaxx+aminx)/2

amaxy  = maxval(sunod(2,1:nnod))
aminy  = minval(sunod(2,1:nnod))
centery=(amaxy+aminy)/2

amaxz  = maxval(sunod(3,1:nnod))
aminz  = minval(sunod(3,1:nnod))
centerz=(amaxz+aminz)/2

do j = 1 , nnod
sunod(1,j)=sunod(1,j)-centerx
sunod(2,j)=sunod(2,j)-centery
sunod(3,j)=sunod(3,j)-centerz
end do



!print *,"the size of geometry"
!print *,'x: ',	 aminx,  amaxx
!print *,'y: ',	 aminy,  amaxy
!print *,'z: ',	 aminz,  amaxz

else
print*,"can not read higher order mesh now"
stop
endif

!print *,'Hello1',   aminz,  amaxz
!open(101,file='input_geo')

!write(101,*) nnod,npat

!do j = 1 , nnod
!write(101,*) sunod(1,j),sunod(2,j),sunod(3,j)
!end do

!print *,'Hello2',   aminz,  amaxz

!do i = 1 , npat
!write(101,*) nsupan(1,i),nsupan(6,i),nsupan(3,i),nsupan(5,i),nsupan(2,i),nsupan(4,i)
!write(101,*) nsupan(1,i),nsupan(3,i),nsupan(5,i),nsupan(2,i),nsupan(4,i),nsupan(6,i)
! write(101,*) nsupan(1,i),nsupan(5,i),nsupan(3,i),nsupan(6,i),nsupan(4,i),nsupan(2,i)
!enddo

!close(101)

!print*,"read curve file is done!"

return
end subroutine incurvefile_ideas
!*************************************************************************
subroutine incurvefile_ideas_sub
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

!print *,"the size of geometry"
!print *,'x: ',	 aminx,  amaxx
!print *,'y: ',	 aminy,  amaxy
!print *,'z: ',	 aminz,  amaxz

else
print*,"can not read higher order mesh now"
stop
endif

!print*,"read curve file is done!"

return
end subroutine incurvefile_ideas_sub
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

!  call rhs_nys_OldBasis(crhs, ntrain, nip)   !incident wave
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
  !   DO k=1,ntriangle*nipp
  !      zzsparse_temp(k,1:(ntriangle*nipp))=zzsparse(((k-1)*ntriangle*nipp+1):(k*ntriangle*nipp))
  !      zzsparse_temp(k,K) = zzsparse_temp(k,K)-0.5_dp*JVector(k)
   !  ENDDO
     
    ! DEALLOCATE(zzsparse)   !zzsparse no longer needed. So release corresponding resources.
     !DEALLOCATE(crhs)       !crhs no longer needed. So release corresponding resources.
    ! DEALLOCATE(JVector)    !JVector no longer needed. So release corresponding resources.

   
     !call ZGESV(ntriangle*nipp,1,zzsparse_temp,ntriangle*nipp,&
     !           IPIV,rj,ntriangle*nipp,INFO)
     
     if(INFO>0) print*, 'Error in LU inversion',INFO
     if(INFO<0) print*, 'Error in LU inputs', -INFO
     
     print*,"direct solve is done!"
  !endif

   !   DEALLOCATE(zzsparse_temp)

  return
end subroutine solve_nys_hard                                                                 
!*************************************************************************
