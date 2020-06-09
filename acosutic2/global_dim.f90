module global_dim
 use,intrinsic:: iso_c_binding
use global_com
!use global_dim_test
  implicit none
  save
  complex(kind=dp),dimension(:),allocatable::rj,jj_mie
  !complex(kind=dp),dimension(:),allocatable::zzsparse
  !complex(c_double_complex), dimension(:):: zzsparse
  !REAL(kind=dp),dimension(:),allocatable::JVector
end module global_dim
