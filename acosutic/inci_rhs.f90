!*************************************************************************
  subroutine rhs_nys_OldBasis(crhs, ntrain, nip)
    use global_com 
    use global_geom,only:nsupan,sunod
    use global_dim
    use gaussintegraldata
    implicit none
    integer::j,q,ielement1,ielement2,ielement0
    real(kind=dp)::vos1(2),vos2(2),vos3(2)
    integer::nodeo_sn(nnpt)
    real(kind=dp)::vcsio(2)
    real(kind=dp)::vippo(3),vhh1(3)
    real(kind=dp)::vuo(3),vvo(3),jacobo,vnormalo(3),vcpuvo(3),zeta(2,nipp),rNodes_field(3,nnpt)
    integer::ntrain, nip
    complex(c_double_complex), dimension(ntrain*nip):: crhs 
    vos1=(/1.0_dp,0.0_dp/)  !right angle triangle vertex points in (u,v) space	
    vos2=(/0.0_dp,1.0_dp/)  !the first number is u value, the second number is v value	
    vos3=(/0.0_dp,0.0_dp/)

  !  allocate(crhs(ntriangle*nipp))
    crhs=(0.0_dp,0.0_dp)
    
    do q=1,ntriangle
       nodeo_sn(:)=nsupan(:,q)  !serial numbers of curve triangle vertex points
       rNodes_field(:,:) = sunod(:,nodeo_sn)
       do j=1,nipp
                vcsio = alphao(j)*vos1+betao(j)*vos2+gammao(j)*vos3
           
                CALL mapGaussPtTo3Ddomain(vcsio(1),vcsio(2),nnpt,rNodes_field,vippo)
                ielement0=(q-1)*nipp+j
                crhs(ielement0)=(-1.0_dp)*exp(cjvk*dot_product(vkk,vippo))

       enddo
    enddo

    return
  end subroutine rhs_nys_OldBasis
!*************************************************************************
