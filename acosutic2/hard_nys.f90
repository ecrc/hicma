!*************************************************************************
!int *nip, int *ntrian, double _Complex *zzsparse_tempt, double _Complex *crhs, int *q, int *p
SUBROUTINE hard_nys_oldBasis_SST(nip, ntrain, zzsparse,  m, n, local_nt, nb) bind(C, name="AL4SAN_Matrix_Generation")
! This subroutine fills zpast
  use, intrinsic :: iso_c_binding
  use global_com
  !use global_dim,only:JVector
  use global_geom,only:nsupan,sunod
  use gaussintegraldata
!  implicit none
  integer::i,j,p,q,k, counter, local_nt
  integer(c_int) ::ntrain
  integer(c_int) ::nip
  real(kind=dp)::vos1(2),vos2(2),vos3(2)
  real(kind=dp)::vcsio(2),vcsis(2),vcsik(2), quadPt(nlqp,2),quadWgt(nlqp,2)
  integer::nodeo_sn(nnpt),nodes_sn(nnpt)
  real(kind=dp)::vippo(3),vipps(3),vippsk(3),unitVec(3),unitVec0(3),unitVeck(3)
  real(kind=dp)::vros_ipp(3),ros_ipp,rpq_min
  real(kind=dp)::jacobo,jacobk,zeta10,zeta20,up0,vp0,nCapDotR,nCapDotRk,nCapDotRf,nCapDotRfk 
  real(kind=dp)::vrrr(3),rrr ,vrrrk(3),rrrk ,vrrrf(3),rrrf,rpf(3),vrrrfk(3),rrrfk,rpfk(3)
  integer::ielement11,ielement12,ielement21,ielement22,ielement00
  complex(kind=dp)::z11,z12,z21,z22,z00, intVal,taylorSum
  real(kind=dp)::ipolator,ipolator0,ipolatort,ipolatork,quadPt_v(nclgp),quadWgt_v(nclgp)                
  real(kind=dp)::vrr0(3),rr0,rNodes_src(3,nnpt),rNodes_field(3,nnpt),quadPt_u(nlqp),quadWgt_u(nlqp)
  integer::mode,nums(2)
  real(kind=dp)::coeffs(nipp,nipp),coeffo(nipp,nipp),rp0(3),jacobo0,jacobot
  complex(c_double_complex), dimension(nb*nb):: zzsparse
!  complex(c_double_complex), dimension(ntrain*nip):: JVector
  REAL(kind=dp),dimension(nipp)::JVector
!************************************************
 !allocate(JVector(1:nipp)) 
  vos1=(/1.0_dp,0.0_dp/)  !right angle triangle vertex points in (u,v) space	
  vos2=(/0.0_dp,1.0_dp/)  !the first number is u value, the second number is v value	
  vos3=(/0.0_dp,0.0_dp/)

  CALL findPolyCoeffs(alphas,betas,gammas,nipp,coeffs) 
  CALL findPolyCoeffs(alphao,betao,gammao,nipp,coeffo) 
  
  quadPt_u = xxo(:)
  quadPt_v = xxs(:)
  quadWgt_u = wwo(:)
  quadWgt_v = wws(:)
  nums = (/nlqp,nclgp/)
  
  quadPt(:,1) = xxo(:)
  quadPt(:,2) = xxs(:)
  quadWgt(:,1) = wwo(:)
  quadWgt(:,2) = wws(:)
!  counter =1
  counter =0
  gcount=0
  !PRINT*,m, n
!loop patches
  DO q=m,(m+local_nt-1)
    !here nsupan (nodeo_sn) is the six (observation) nodes for each patch
     nodeo_sn(:)=nsupan(:,q)  !serial numbers of field triangular patches
     rNodes_field(:,:) = sunod(:,nodeo_sn)      
!     PRINT*,q,"of",ntriangle

    DO p=n,(n+local_nt-1)
        !here nsupan (nodes_sn) is the six (source) nodes for each patch
        nodes_sn(:)=nsupan(:,p)
        rNodes_src(:,:) = sunod(:,nodes_sn)
                
        IF(p.ne.q) THEN
            rpq_min=88888888.0_dp
            DO j=1,nipp
                vcsio = alphao(j)*vos1+betao(j)*vos2+gammao(j)*vos3
                CALL mapGaussPtTo3Ddomain(vcsio(1),vcsio(2),nnpt,rNodes_field,vippo)
                
                DO i=1,nipp
                    vcsis = alphas(i)*vos1+betas(i)*vos2+gammas(i)*vos3 
                    CALL mapGaussPtTo3Ddomain(vcsis(1),vcsis(2),nnpt,rNodes_src,vipps)
                    vros_ipp(:)=vippo(:)-vipps(:)
                    ros_ipp=SQRT(DOT_PRODUCT(vros_ipp,vros_ipp))
                    IF(ros_ipp.lt.rpq_min) rpq_min=ros_ipp
                ENDDO
            ENDDO

!*************************************************************************
 !far patch interaction
!*************************************************************************
           IF(rpq_min.gt.nearpd)  THEN   !far patch interaction
              !PRINT*,"computing for far-patch = ",p
              !counter=1
              DO j=1,nipp
                 ! ielement11=(q-1)*nipp+j
                  ielement11=j
                  vcsio = alphao(j)*vos1+betao(j)*vos2+gammao(j)*vos3
                  CALL mapGaussPtTo3Ddomain(vcsio(1),vcsio(2),nnpt,rNodes_field,vippo)
                  CALL calcJacobian(vcsio(1),vcsio(2),nnpt,rNodes_field,jacobot)
                  CALL basis_function(vcsio,nipp,coeffo,j,ipolatort)
                  JVector(ielement11)=1.0_dp/jacobot*ipolatort
                  
                  DO i=1,nipp
                      !PRINT*,"i = ",i
                      !ielement00=((q-1)*nipp+j-1)*ntriangle*nipp+(p-1)*nipp+i
                      !ielement00=offset+counter
                      !ielement00=counter
                      !counter=counter+1
                      !ielement00=((j-1)*nb+i+counter)
                      ielement00=((i-1)*nb+j+counter)
                      vcsis = alphas(i)*vos1+betas(i)*vos2+gammas(i)*vos3  
                      CALL mapGaussPtTo3Ddomain(vcsis(1),vcsis(2),nnpt,rNodes_src,vipps)
                      CALL calcUnitNormal(vcsis(1),vcsis(2),nnpt,rNodes_src,unitVec)
                      
                      vrrr(:)=vippo(:)-vipps(:)
                      rrr=SQRT(DOT_PRODUCT(vrrr,vrrr))        !|r-rp|
                      nCapDotR = DOT_PRODUCT(unitVec,vrrr)
                      
!                      zzsparse(ielement00) = ws(i)*0.5_dp/(4.0_dp*pid)*nCapDotR*(exp(cjvk*rrr)*(1.0_dp-cjvk*rrr)/(rrr**3)) !matrix fill for far intraction
 IF ((p .eq. q) .AND.  (j .eq. i)) THEN
                      zzsparse(ielement00) = (ws(i)*0.5_dp/(4.0_dp*pid)*nCapDotR*(exp(cjvk*rrr)*(1.0_dp-cjvk*rrr)/(rrr**3)))&
-0.5_dp*JVector(ielement11) !matrix fill for far intraction
ELSE
                      zzsparse(ielement00) = ws(i)*0.5_dp/(4.0_dp*pid)*nCapDotR*(exp(cjvk*rrr)*(1.0_dp-cjvk*rrr)/(rrr**3))                 
   ENDIF
                  ENDDO
                  
              ENDDO
              
!*************************************************************************
!near patch interaction...by setting higher degree gauss quadrature points for source triangle
!*************************************************************************
           ELSE
              !PRINT*,"computing for near-patch = ",p
              mode = 2
              !counter=1
              DO j=1,nipp
                  !PRINT*,"j = ",j
                  !ielement11=(q-1)*nipp+j
                  ielement11=j
                  vcsio = alphao(j)*vos1+betao(j)*vos2+gammao(j)*vos3
                  CALL mapGaussPtTo3Ddomain(vcsio(1),vcsio(2),nnpt,rNodes_field,vippo)
                  CALL calcJacobian(vcsio(1),vcsio(2),nnpt,rNodes_field,jacobot)
                  CALL basis_function(vcsio,nipp,coeffo,j,ipolatort)
                  JVector(ielement11)=1.0_dp/jacobot*ipolatort
                  
                  DO i=1,nipp
                      !PRINT*,"i = ",i
                      !ielement00=((q-1)*nipp+j-1)*ntriangle*nipp+(p-1)*nipp+i
                      !ielement00=offset+counter
                      !ielement00=counter
                      !counter=counter+1
                      !ielement00=((j-1)*nb+i+counter)
                      ielement00=((i-1)*nb+j+counter)
                        z00 = (0.0_dp,0.0_dp)
                        DO k=1,nhdgqp_near
                              vcsik = alphas22(k)*vos1+betas22(k)*vos2+gammas22(k)*vos3
                              CALL mapGaussPtTo3Ddomain(vcsik(1),vcsik(2),nnpt,rNodes_src,vippsk)
                              CALL calcUnitNormal(vcsik(1),vcsik(2),nnpt,rNodes_src,unitVec)
                              CALL basis_function(vcsik,nipp,coeffs,i,ipolator)
                           
                              vrrrk(:)=vippo(:)-vippsk(:)
                              rrrk=SQRT(DOT_PRODUCT(vrrrk,vrrrk))        !|r-rp| 
                              nCapDotR = DOT_PRODUCT(unitVec,vrrrk)
                                    
                              z00 = z00 + ws22(k)*0.5_dp/(4.0_dp*pid)*nCapDotR*(exp(cjvk*rrrk)/&
                                               (rrrk**3))*(1.0_dp-cjvk*rrrk)*ipolator
                           
                        ENDDO
     !                   zzsparse(ielement00) = z00
 IF ((p .eq. q) .AND.  (j .eq. i)) THEN
                        zzsparse(ielement00) = z00-0.5_dp*JVector(ielement11)
 Else
    zzsparse(ielement00) = z00
 ENDIF 
                  ENDDO
               ENDDO

           ENDIF
           
!*************************************************************************
!self patch interaction...by applying Improved Guiggiani's method.
!*************************************************************************
        ELSE     !self patch interaction
              !PRINT*,"computing for self-patch = ",p
              !counter=1
              DO j=1,nipp
                  !PRINT*,"j = ",j
                  !ielement11=(q-1)*nipp+j
                  ielement11=j
                  vcsio = alphao(j)*vos1+betao(j)*vos2+gammao(j)*vos3
                  CALL mapGaussPtTo3Ddomain(vcsio(1),vcsio(2),nnpt,rNodes_field,vippo)
                  CALL calcJacobian(vcsio(1),vcsio(2),nnpt,rNodes_field,jacobot)
                  CALL basis_function(vcsio,nipp,coeffo,j,ipolatort)
                  JVector(ielement11)=1.0_dp/jacobot*ipolatort
                  
                  DO i=1,nipp
                      !PRINT*,"i = ",i
                      !ielement00=((q-1)*nipp+j-1)*ntriangle*nipp+(p-1)*nipp+i
                      !ielement00=offset+counter
                      !ielement00=counter
                      !counter=counter+1
                      !ielement00=((j-1)*nb+i+counter)
                      ielement00=((i-1)*nb+j+counter)
                      CALL OriginalGuiggiani_hard(vcsio,nipp,nums,quadPt_u,quadWgt_u,&
                                                quadPt_v,quadWgt_v,coeffs,i,rNodes_field,nnpt,cjvk,z00)
     !                 zzsparse(ielement00) = z00

 IF ((p .eq. q) .AND.  (j .eq. i)) THEN
                      zzsparse(ielement00) = z00-0.5_dp*JVector(ielement11)
 ELSE
           zzsparse(ielement00) = z00
ENDIF 
                  ENDDO
              ENDDO
              
        ENDIF
        counter=counter+local_nt*nipp*nipp
    ENDDO
     gcount=gcount+nipp
     counter=gcount
     !fortran but double check
     !counter=gcount+local_nt*nipp*nipp
     !gcount=gcount+local_nt*nipp*nipp

     !counter=gcount+nb*nipp
     !gcount=gcount+nb*nipp
  ENDDO

  !PRINT*,"Calculations for hard-boundary condition are done!"
  !  DEALLOCATE(JVector)
  RETURN
  
END SUBROUTINE hard_nys_oldBasis_SST
!*************************************************************************
