module global_geom 
  use global_com,only:dp
  implicit none
  save
  integer::nglunk,nsuunk,nwiunk,nwiseg,nwinodes,njs,nwi,ndielunk,npecunk 
  integer::nsup_unk
! pec body parameters
  real(kind=dp),dimension(:,:),allocatable::sunod
!  real(kind=dp),dimension(:,:),allocatable::unormal
!  real(kind=dp),dimension(:),allocatable::patch_area

  integer,dimension(3)::nsuinf
! nsuinf(3) is not equal to nsuunk!
  integer,dimension(:,:),allocatable::nsupan,nsupae
  integer,dimension(:,:),allocatable::nsuedn,nsuedp,basis_nodes
  integer,dimension(:,:),allocatable::super_nodes
! juction edge paramaters:
  integer::njunction_edges
  integer,dimension(:),allocatable::junction_store
  integer,dimension(:,:),allocatable::junction_patches
  integer::max_pae
! ******* WARNING: JUNCTION EDGES ARE INVISIBLE TO NSUPAE
!         the second patch of a junction edge does not list the junction edge
!         as one of its edges
  integer::disconnected_body_count
  type body_patch
     integer,pointer::id(:)
  end type body_patch
  type(body_patch)::patches_of_body(256) !at most 100 disconnected bodies

! ------------------- wiring --------------------
! radius of wires and node coordinates
  real(kind=dp),dimension(:),allocatable::rwiseg,rjbasis
  real(kind=dp),dimension(:,:),allocatable::winod
! nodes of segment,segments of basis, nodes of basis
  integer,dimension(:,:),allocatable::node_of_wseg,wbasis_seg,wbasis_nodes,&
       jbasis_nodes 
! ------------------- wiring --------------------
! ------------------- wiring and TL --------------------
 integer,dimension(:),allocatable::nofwseg_on_node
  type segment_of_node
     integer,pointer::p(:)
  end type segment_of_node
  type(segment_of_node),dimension(:),allocatable::wseg_on_node
! ------------------- wiring and TL --------------------
! Resistor/Impedance loading commons!
  integer::nnloaded
  integer,allocatable,dimension(:):: nload
  complex(kind=dp),allocatable,dimension(:)::rrload
! dielectric body parameters
  real(kind=dp),allocatable,dimension(:,:)::r_n    !location of nodes
  real(kind=dp),allocatable,dimension(:)::a_s    !triangle area
  real(kind=dp),allocatable,dimension(:)::v_v    !tetrahedron volume
!  real(kind=dp),allocatable,dimension(:,:)::rc_s   !triangle centroid
  ! the centroid might not be necessary
!  real(kind=dp),allocatable,dimension(:,:)::rc_v   !tetrahedron centroid
  real(kind=dp),allocatable,dimension(:)::epsr_v !tetrahedron relative permittivity
  real(kind=dp),allocatable,dimension(:)::kap_v  ! tetrahedron kappa
  real(kind=dp),allocatable,dimension(:)::sigma_v  ! tetrahedron kappa
  real(kind=dp)::loss_factor

  integer::N_nodes,N_tri,N_tetra
  integer,allocatable,dimension(:,:)::iv_n !pointer to tetrahedron's nodes
  integer,allocatable,dimension(:,:)::is_n !pointer to triangle's nodes
  integer,allocatable,dimension(:,:)::is_v !pointer to triangle's tetrahedrons
  integer,allocatable,dimension(:,:)::is_fn!pointer to triangle's free nodes
  integer,allocatable::spatch_of_vface(:),tetras_of_spatch(:,:),faces_of_tetra(:,:)
  real(kind=dp)::diel_const

! frill parameters (currently not operational)
contains
    function sunod1(k)
    implicit none
    integer,intent(in)::k
    real(kind=dp)::sunod1(1:3)
    sunod1(1:3)=sunod(1:3,k)
    return
  end function sunod1

  function unormal1(k)
    implicit none
    integer,intent(in)::k
    real(kind=dp)::unormal1(1:3)
    read(unit=95,rec=k) unormal1(1:3)
    return
  end function unormal1

  function patch_area1(k)
    implicit none
    integer,intent(in)::k
    real(kind=dp)::patch_area1
    read(unit=94,rec=k) patch_area1
    return
  end function patch_area1
end module global_geom 
