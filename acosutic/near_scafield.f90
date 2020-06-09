
!*************************************************************************
 SUBROUTINE near_sca_hard_near_scheme
  use global_com
  use global_dim,only:rj
  use global_geom,only:nsupan,sunod
  use gaussintegraldata
  implicit none  
  integer::p,i,k,kkk,nnn,mmm
  integer::nodes_sn(nnpt)
  REAL(kind=dp)::vos1(2),vos2(2),vos3(2),vcsik(2),vcsis(2),jacobk      
  REAL(kind=dp)::vus(3),vvs(3),vipps(3),vobp(3),vippsk(3)         
  REAL(kind=dp)::dthe,dphi,phi_degree,the_degree
  REAL(kind=dp),dimension(3)::vekk,vethe,vephi
  integer::kph,kth
  !complex(kind=dp)::jdotp_the,jdotp_phi,esca_the1,esca_phi1
  !complex(kind=dp)::esca_the,esca_phi
  complex(kind=dp)::nearsca_hard,temp,temp_N1,temp_N2
  REAL(kind=dp)::ipolatork
  REAL(kind=dp)::vrrr(3),rrr, vrrrk(3),rrrk ,rpf(3)
  !integer::IPIV(nipp)
  REAL(kind=dp)::rNodes_src(3,nnpt),unitVec(3),unitVec0(3),nCapDotR,nCapDotRf
 
  REAL(kind=dp)::coeffs(nipp,nipp),rp0(3),zeta10,zeta20,up0,vp0,&
                jacobo0,jacobo
                   
  
  vos1=(/1.0_dp,0.0_dp/)  !right angle triangle vertex points in (u,v) space	
  vos2=(/0.0_dp,1.0_dp/)  !the first number is u value, the second number is v value	
  vos3=(/0.0_dp,0.0_dp/)
  
  CALL findPolyCoeffs(alphas,betas,gammas,nipp,coeffs)
  
  OPEN(unit=5702,file='/home/omairyrm/hicma-dev/acosutic/near_sca_hard_near_scheme_OldBasis_1_120.out')

  IF (nthe_nf==1) THEN 
     dthe=0.0_dp
  ELSE
     dthe=(the_end_degree_nf-the_sta_degree_nf)/REAL(nthe_nf-1,dp)
  ENDIF

  IF (nphi_nf==1) THEN 
     dphi=0.0_dp
  ELSE
     dphi=(phi_end_degree_nf-phi_sta_degree_nf)/REAL(nphi_nf-1,dp) 
  ENDIF
  
  DO kph=1,nphi_nf ! Loop over phi
     phi_degree=phi_sta_degree_nf+(kph-1)*dphi ! Compute phi
     vephi(1)=-sin(phi_degree*pid/180.0_dp)                                            
     vephi(2)=cos(phi_degree*pid/180.0_dp)                                              
     vephi(3)=0.0_dp                                                    

     DO kth=1,nthe_nf ! Loop over theta
        the_degree=the_sta_degree_nf+(kth-1)*dthe ! Compute theta
        ! Compute r,phi,and theta directed unit vectors. 
        vekk(1)=sin(the_degree*pid/180.0_dp)*cos(phi_degree*pid/180.0_dp)                                
        vekk(2)=sin(the_degree*pid/180.0_dp)*sin(phi_degree*pid/180.0_dp)                                
        vekk(3)=cos(the_degree*pid/180.0_dp)
        vethe(1)=cos(the_degree*pid/180.0_dp)*cos(phi_degree*pid/180.0_dp)                                  
        vethe(2)=cos(the_degree*pid/180.0_dp)*sin(phi_degree*pid/180.0_dp)                                  
        vethe(3)=-sin(the_degree*pid/180.0_dp)

        vobp(1)=rnfobp*sin(the_degree*pid/180.0_dp)*cos(phi_degree*pid/180.0_dp)  !observation point position
        vobp(2)=rnfobp*sin(the_degree*pid/180.0_dp)*sin(phi_degree*pid/180.0_dp)
        vobp(3)=rnfobp*cos(the_degree*pid/180.0_dp)

        nearsca_hard=(0.0_dp,0.0_dp) 
        
        DO p=1,ntriangle
           nodes_sn(:)=nsupan(:,p)
           rNodes_src(:,:) = sunod(:,nodes_sn)
           
           DO i=1,nipp
               
               DO k=1,nhdgqp_far
                  vcsik = alphas11(k)*vos1+betas11(k)*vos2+gammas11(k)*vos3
                  CALL mapGaussPtTo3Ddomain(vcsik(1),vcsik(2),nnpt,rNodes_src,vippsk)
                  CALL calcUnitNormal(vcsik(1),vcsik(2),nnpt,rNodes_src,unitVec)
                  CALL basis_function(vcsik,nipp,coeffs,i,ipolatork)
                  !CALL calcJacobian(vcsik(1),vcsik(2),nnpt,rNodes_src,jacobk)
              
                  vrrrk(:)=vobp(:)-vippsk(:)
                  rrrk=SQRT(DOT_PRODUCT(vrrrk,vrrrk))        !|r-rp| 
                  nCapDotR = DOT_PRODUCT(unitVec,vrrrk)
                  
                  nearsca_hard = nearsca_hard + ws11(k)*0.5_dp/(4.0_dp*pid)*nCapDotR*(exp(cjvk*rrrk)/&
                                   (rrrk**3))*(1.0_dp-cjvk*rrrk)*rj((p-1)*nipp+i)*ipolatork !*jacobk
                  
               ENDDO
              
           ENDDO
           
        ENDDO
   
        WRITE(5702,*) the_degree, phi_degree, REAL(nearsca_hard), AIMAG(nearsca_hard), ABS(nearsca_hard)

     ENDDO
  ENDDO

  CLOSE(5702)

  PRINT*,"Calculations for hard near scattered field near scheme are done!"

  RETURN
END SUBROUTINE near_sca_hard_near_scheme
!*************************************************************************





