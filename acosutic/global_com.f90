module global_com   
  !
  !  This module contains all the constants which are common to the whole code
  !
USE, INTRINSIC :: ISO_C_BINDING 
  implicit none
  save
  integer,parameter::dp=kind(0.0d0),sp=kind(0.0)
  complex(kind=dp),parameter::cjj=(0.0_dp,1.0_dp)
  real(kind=dp)::the_sta,the_end,phi_sta,phi_end
  real(kind=dp)::the_sta_degree,the_end_degree,phi_sta_degree,phi_end_degree
  real(kind=dp)::the_sta_degree_nf,the_end_degree_nf,phi_sta_degree_nf,phi_end_degree_nf
  integer::nthe_rcs,nphi_rcs,nthe_nf,nphi_nf
  real(kind=dp)::alpha,veta 
  real(kind=dp)::pid,cd,eps0d,rmu0d,loss_sigma
  real(kind=dp)::fre_single
  real(kind=dp)::real_mem=8.0_dp,complex_mem=16.0_dp,int_mem=4.0_dp
  integer::order_tri,nnpt
  integer::order_basis,nipp,order_old
  integer::nhdgqp_far,nhdgqp_near,nlqp,nitermax,ngprcs,denmrcs,nclgp,defile,denscaf,detoj,deonep,nusolver,depconvert
  integer::dejjmie,detotalp, nhd_near
  real(kind=dp)::aminx,aminy,aminz,amaxx,amaxy,amaxz
  integer(kind = c_int)::ntriangle,nnode
  real(kind=dp)::vwavelength,wavenumkk,rnfobp
  complex(kind=dp)::cjvk
  real(kind=dp)::vkk(3),vee(3),vhh(3)
  integer::it_solv,restart
  real(kind=dp)::precis,nearpd
  integer::no_freqs


  ! TIMING
  real(kind=dp)::time_far_allinc=0.0_dp,time_near_allinc=0.0_dp,&
       time_it_allinc=0.0_dp,time_frcs_allinc=0.0_dp,&
       time_outcur_allinc=0.0_dp,time_solve_allinc=0.0_dp,&
       time_rhs_allinc=0.0_dp,time_pre_allinc=0.0_dp,&
       time_rhs_CKT_allinc=0.0_dp,time_pre_CKT_allinc=0.0_dp,&
       time_CKT_allinc=0.0_dp
  real(kind=dp)::tot_far_time=0.0_dp,tot_near_time=0.0_dp,&
       solve_time_it=0.0_dp,tot_frcs_time=0.0_dp,&
       tot_outcur_time=0.0_dp,tot_solve_time=0.0_dp,&
       rhs_time=0.0_dp,tot_pre_time=0.0_dp,&
       rhs_CKT_time=0.0_dp,tot_pre_CKT_time=0.0_dp,&
       tot_CKT_time=0.0_dp
end module global_com
