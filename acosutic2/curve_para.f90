!===========================================
subroutine MatVecProd_real(Mat,m,vec)     !real product of a matrix and vector

  use global_com,only:dp

  implicit none
  INTEGER,INTENT(IN):: m
  REAL(kind=dp), INTENT(INOUT):: vec(m)
  REAL(kind=dp), INTENT(IN):: Mat(m,m)
  
  INTEGER:: i
  REAL(kind=dp):: temp(m),tempM(m)
  
  DO i = 1,m
        tempM(:) = Mat(i,:)
        temp(i) = DOT_PRODUCT(tempM,vec)
  ENDDO
  vec(:) = temp(:) 
  
END SUBROUTINE MatVecProd_real
!!******************************************************************************************END        

    
    
    
SUBROUTINE real_cross_product(va,vb,vc)     !cross product of two real numbers

  use global_com,only:dp

  implicit none
  real(kind=dp):: va(3),vb(3),vc(3)

  vc = 0.0_dp
  vc(1)=va(2)*vb(3)-va(3)*vb(2) 
  vc(2)=va(3)*vb(1)-va(1)*vb(3)
  vc(3)=va(1)*vb(2)-va(2)*vb(1) 

END subroutine real_cross_product
!******************************************************************************************END        
    
    
    
    
    
!===========================================      
subroutine getinverse(a,c,n)
  !
  ! Inverse matrix
  ! Method: Based on DOolittle LU factorization for Ax=b
  ! Alex G. December 2009
  !-----------------------------------------------------------
  ! input ...
  ! a(n,n) - array of coefficients for matrix A
  ! n      - dimension
  ! output ...
  ! c(n,n) - inverse matrix of A
  ! comments ...
  ! the original matrix a(n,n) will be destroyed 
  ! during the calculation
  !
  use global_com,only:dp

  implicit none 
  integer n
  real(kind=dp):: a(n,n), c(n,n)
  real(kind=dp):: L(n,n), U(n,n), b(n), d(n), x(n)
  real(kind=dp):: coeff
  integer i, j, k

  ! step 0: initialization for matrices L and U and b
  ! Fortran 90/95 aloows such operations on matrices
  L=0.0_dp
  U=0.0_dp
  b=0.0_dp

  ! step 1: forward elimination
  DO k=1, n-1
     DO i=k+1,n
        coeff=a(i,k)/a(k,k)
        L(i,k) = coeff
        DO j=k+1,n
           a(i,j) = a(i,j)-coeff*a(k,j)
        END DO
     END DO
  END DO

  ! Step 2: prepare L and U matrices 
  ! L matrix is a matrix of the elimination coefficient
  ! + the diagonal elements are 1.0
  DO i=1,n
     L(i,i) = 1.0_dp
  END DO
  ! U matrix is the upper triangular part of A
  DO j=1,n
     DO i=1,j
        U(i,j) = a(i,j)
     END DO
  END DO

  ! Step 3: compute columns of the inverse matrix C
  DO k=1,n
     b(k)=1.0_dp
     d(1) = b(1)
     ! Step 3a: Solve Ld=b using the forward substitution
     DO i=2,n
        d(i)=b(i)
        DO j=1,i-1
           d(i) = d(i) - L(i,j)*d(j)
        END DO
     END DO
     ! Step 3b: Solve Ux=d using the back substitution
     x(n)=d(n)/U(n,n)
     DO i = n-1,1,-1
        x(i) = d(i)
        DO j=n,i+1,-1
           x(i)=x(i)-U(i,j)*x(j)
        END DO
        x(i) = x(i)/u(i,i)
     END DO
     ! Step 3c: fill the solutions x(n) into column k of C
     DO i=1,n
        c(i,k) = x(i)
     END DO
     b(k)=0.0_dp
  END DO
END subroutine getinverse
!******************************************************************************************END    

    
 !***********************************************************************************************    
SUBROUTINE OriginalGuiggiani_hard(vcsio,nipp,nums,quadPt_u,quadWgt_u,&
                    quadPt_v,quadWgt_v,coeffs,i,rNodes_field,nnpt,cjvk,z00)

    USE global_com,ONLY:dp,pid
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: nnpt,nipp,i,nums(2)
    COMPLEX(kind=dp), INTENT(OUT) :: z00
    REAL(kind=dp), INTENT(IN)::vcsio(2),coeffs(nipp,nipp)
    REAL(kind=dp), INTENT(IN)::rNodes_field(3,nnpt),quadPt_u(nums(1))
    REAL(kind=dp), INTENT(IN)::quadPt_v(nums(2)),quadWgt_u(nums(1))
    REAL(kind=dp), INTENT(IN)::quadWgt_v(nums(2))
    COMPLEX(kind=dp), INTENT(IN) :: cjvk
    
    INTEGER :: aa,ff,uu,vv,N_uu,N_vv
    REAL(kind=dp)::U1(3),U2(3),nCap(3),beta(3),nCapDotR
    REAL(kind=dp)::etaVecS(2),eta(2),zeta(2),temp_T(2),Lij,rho
    REAL(kind=dp):: R_bar(3),rrr,z,thetaBarPrime,th_bar,th
    REAL(kind=dp)::xxu(nums(1)),xxv(nums(2)),rhoCap_th,vipps(3),vippo(3)
    REAL(kind=dp):: thetaBar(2,3),theta(2,3),h_alpha(3)
    REAL(kind=dp):: wwv(nums(2)),wwu(nums(2)),rhoCap(3)
    COMPLEX(kind=dp) :: z11,F_rho_th    
    
    !OPEN(unit=5001,file='polar_Gauss_points.out')
    
   !parameters for 1-D gauss quadrature
    xxu(:) = quadPt_u
    xxv(:)= quadPt_v
    wwu(:) = quadWgt_u
    wwv(:)= quadWgt_v
    N_uu = nums(1)
    N_vv = nums(2)
   
    CALL calcParams(vcsio,thetaBar,h_alpha,beta,rhoCap)     !Relevant quantities needed

  !Computing singular integral using Original Guiggiani formula
    z00 = 0.0_dp 
    
    DO aa = 1,3     !per sub-triangle

        z11 = 0.0_dp
        
        !Along x direction
        DO uu = 1,N_uu
 
            th_bar = (thetaBar(1,aa)-thetaBar(2,aa))*xxu(uu) + thetaBar(2,aa) 
            
            th = th_bar + beta(aa)
            
            rhoCap_th = h_alpha(aa)/cos(th_bar)
            thetaBarPrime = thetaBar(1,aa)-thetaBar(2,aa)
            
            !Along rho direction
            DO vv = 1,N_vv
                
               !map (rho,theta) to eta coordinates
                rho = rhoCap_th*xxv(vv)
                zeta(1) = vcsio(1) + rho*cos(th)
                zeta(2) = vcsio(2) + rho*sin(th)
                
                !Error detection mechanism
                IF ((zeta(1)<0.0_dp).OR.(zeta(2)<0.0_dp).OR.(zeta(1)>1.0_dp).OR.(zeta(2)>1.0_dp)) THEN
                    PRINT*,"zeta = ", zeta
                    PRINT*,"Error detected in zeta variable"
                    STOP
                ENDIF                
                
               !unit normal evaluated at source point in the zeta plane
                CALL calcUnitNormal(zeta(1),zeta(2),nnpt,rNodes_field,nCap)
                
               !corresponding local distance
                CALL mapGaussPtTo3Ddomain(zeta(1),zeta(2),nnpt,rNodes_field,vipps)
                CALL mapGaussPtTo3Ddomain(vcsio(1),vcsio(2),nnpt,rNodes_field,vippo)
                R_bar = vippo-vipps
                rrr = SQRT(DOT_PRODUCT(R_bar,R_bar))
                nCapDotR = DOT_PRODUCT(nCap,R_bar)
                
                !basis function, Lij
                CALL basis_function(zeta,nipp,coeffs,i,Lij)         
                
               !singular integral kernel evaluated at source point
                F_rho_th = nCapDotR*exp(cjvk*rrr)*(1.0_dp-cjvk*rrr)&
                                        /(4.0_dp*pid*rrr**3)*Lij*rho                
               ! I(1) summation
                z11 = z11 + wwu(uu)*wwv(vv)*F_rho_th*rhoCap_th*thetaBarPrime
                            
                
            ENDDO
        ENDDO
        
        z00 = z00 + z11        !Total singular integral value
    ENDDO
    
    RETURN
  
    
    CONTAINS
    
        SUBROUTINE calcParams(vcsio,thetaBar,h_alpha,beta,rhoCap)
            !Subroutine to find the point on the closest edge to the projected field point
            
            IMPLICIT NONE
    
            REAL(kind=dp), INTENT(IN) :: vcsio(2)   
            REAL(kind=dp), INTENT(OUT)::beta(3),h_alpha(3),thetaBar(2,3),rhoCap(3)
    
            INTEGER :: aa
            REAL(kind=dp)::Rs1(3),Rs2(3),Rs3(3),R31(3),R23(3),R12(3)
            REAL(kind=dp)::zetaCap(3,2),nCap(3),zetaVec1(2),zetaVec2(2),zetaVec3(2),zetaVecS(2)
            REAL(kind=dp)::mCap(3,3),temp(3)
            
           !compute eta_2
            zetaVec1 = (/1.0_dp, 0.0_dp/)
            zetaVec2 = (/0.0_dp, 1.0_dp/)
            zetaVec3 = (/0.0_dp, 0.0_dp/)
            zetaVecS = vcsio
    
           !local unit vectors in the zeta domain
            zetaCap(:,1) = (/1.0_dp,0.0_dp,0.0_dp/)
            zetaCap(:,2) = (/0.0_dp,1.0_dp,0.0_dp/)
            nCap = (/0.0_dp,0.0_dp,1.0_dp/)
            
            Rs1 = (zetaVec1(1)-zetaVecS(1))*zetaCap(:,1) + (zetaVec1(2)-zetaVecS(2))*zetaCap(:,2)
            Rs2 = (zetaVec2(1)-zetaVecS(1))*zetaCap(:,1) + (zetaVec2(2)-zetaVecS(2))*zetaCap(:,2)
            Rs3 = (zetaVec3(1)-zetaVecS(1))*zetaCap(:,1) + (zetaVec3(2)-zetaVecS(2))*zetaCap(:,2)
            R12 = (zetaVec2(1)-zetaVec1(1))*zetaCap(:,1) + (zetaVec2(2)-zetaVec1(2))*zetaCap(:,2)
            R23 = (zetaVec3(1)-zetaVec2(1))*zetaCap(:,1) + (zetaVec3(2)-zetaVec2(2))*zetaCap(:,2)
            R31 = (zetaVec1(1)-zetaVec3(1))*zetaCap(:,1) + (zetaVec1(2)-zetaVec3(2))*zetaCap(:,2)
    
            !
            rhoCap=(/SQRT(DOT_PRODUCT(Rs3,Rs3)),SQRT(DOT_PRODUCT(Rs1,Rs1)),&
                        SQRT(DOT_PRODUCT(Rs2,Rs2))/)
    
            CALL real_cross_product(R23,nCap,temp)
            mCap(:,1) = temp/SQRT(DOT_PRODUCT(temp,temp))
    
            CALL real_cross_product(R31,nCap,temp)
            mCap(:,2) = temp/SQRT(DOT_PRODUCT(temp,temp))
    
            CALL real_cross_product(R12,nCap,temp)
            mCap(:,3) = temp/SQRT(DOT_PRODUCT(temp,temp))
            
            h_alpha(1) = DOT_PRODUCT(Rs3,mCap(:,1))
            h_alpha(2) = DOT_PRODUCT(Rs1,mCap(:,2))
            h_alpha(3) = DOT_PRODUCT(Rs2,mCap(:,3))
    
            thetaBar(:,1) = (/ACOS(h_alpha(1)/rhoCap(1)),&
                             -ACOS(h_alpha(1)/rhoCap(3))/)
            thetaBar(:,2) = (/ACOS(h_alpha(2)/rhoCap(2)),&
                             -ACOS(h_alpha(2)/rhoCap(1))/)
            thetaBar(:,3) = (/ACOS(h_alpha(3)/rhoCap(3)),&
                             -ACOS(h_alpha(3)/rhoCap(2))/)
    
            beta(3) = ACOS(DOT_PRODUCT((R31/SQRT(DOT_PRODUCT(R31,R31))),&
                            mCap(:,3)))
            
            beta(2) = 3.0_dp/2.0_dp*pid
            beta(1) =  beta(3) + thetaBar(1,3) + abs(thetaBar(2,1)) 
            
        END SUBROUTINE calcParams        
                
END SUBROUTINE OriginalGuiggiani_hard   
!*****************************************************************************************END      
  
    
    
 !***********************************************************************************************    
SUBROUTINE ImprovedGuiggianiMethod(vcsio,nipp,nums,quadPt_u,quadWgt_u,&
                    quadPt_v,quadWgt_v,coeffs,i,rNodes_field,nnpt,cjvk,z00)

    USE global_com,ONLY:dp,pid
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: nnpt,nipp,i,nums(2)
    COMPLEX(kind=dp), INTENT(OUT) :: z00
    REAL(kind=dp), INTENT(IN)::vcsio(2),coeffs(nipp,nipp)
    REAL(kind=dp), INTENT(IN)::rNodes_field(3,nnpt),quadPt_u(nums(1))
    REAL(kind=dp), INTENT(IN)::quadPt_v(nums(2)),quadWgt_u(nums(1))
    REAL(kind=dp), INTENT(IN)::quadWgt_v(nums(2))
    COMPLEX(kind=dp), INTENT(IN) :: cjvk
    
    INTEGER :: aa,ff,uu,vv,N_uu,N_vv, wrt1=1, wrt2=2, derivativeDeg=1
    REAL(kind=dp)::U1(3),U2(3),nCap(3),beta(3),etaVec2(2),nCapDotR
    REAL(kind=dp)::etaVecS(2),eta(2),zeta(2),temp_T(2),Lij0,Lij,f_1,rho
    REAL(kind=dp):: R_bar(3),rrr,z,sigma_z,thetaBarPrime,th_bar,th,rr0
    REAL(kind=dp)::xxu(nums(1)),xxv(nums(2)),rhoCap_th,vipps(3),vippo(3)
    REAL(kind=dp):: Z_alpha(2,3),theta(2,3),h_alpha(3),R_cap(3)
    REAL(kind=dp)::Det_T_1,wwu(nums(1)),wwv(nums(2)),tVec_1(2,2),rhoCap(3)
    COMPLEX(kind=dp) :: z11,z22,F_rho_th    
    
   !parameters for 1-D gauss quadrature
    xxu(:) = quadPt_u
    xxv(:)= quadPt_v
    wwu(:) = quadWgt_u
    wwv(:)= quadWgt_v
    N_uu = nums(1)
    N_vv = nums(2)
    
   ! Column vectors of Jacobian matrix from physical to zeta plane
    CALL calcPosVectDerivative(derivativeDeg,wrt1,vcsio(1),vcsio(2),nnpt,&
                                rNodes_field,U1)   !dr/dZeta1
    CALL calcPosVectDerivative(derivativeDeg,wrt2,vcsio(1),vcsio(2),&
                                nnpt,rNodes_field,U2)   !dr/dZeta2
    
   !Relevant quantities needed
    CALL computeParams_IGM(U1,U2,Det_T_1,Z_alpha,etaVecS,h_alpha,tVec_1,beta,rhoCap,etaVec2)
    
  !Computing singular integral using improved Guiggiani formula
   
    z00 = 0.0_dp        
    DO aa = 1,3             !per subtriangle
        z11 = 0.0_dp
        
        
        DO uu = 1,N_uu              !Along theta direction
            z = (Z_alpha(1,aa)-Z_alpha(2,aa))*xxu(uu) + Z_alpha(2,aa)     ! z(x)
            
           !m = 3 
            sigma_z = (z**3)/(z**3+(1.0_dp-z)**3)       ! sigma(z) for m = 3
            thetaBarPrime = 3.0_dp*pid*(Z_alpha(1,aa)-Z_alpha(2,aa))*&
                            sigma_z*((1.0_dp-z)**2)/(z*(z**3 + (1.0_dp-z)**3))        
            
           !!m = 2  
            !sigma_z = (z**2)/(z**2+(1.0_dp-z)**2)       ! sigma(z) for m = 2
            !thetaBarPrime = 2.0_dp*pid*(Z_alpha(1,aa)-Z_alpha(2,aa))*sigma_z*&
            !                (1.0_dp-z)/(z*(z**2 + (1.0_dp-z)**2))               
            
            th_bar = pid*sigma_z - 0.5_dp*pid                   
            th = th_bar + beta(aa)
            
            rhoCap_th = h_alpha(aa)/cos(th_bar)
            
            
            DO vv = 1,N_vv          !Along rho direction
                
               !map (rho,theta) to eta coordinates
                rho = rhoCap_th*xxv(vv)
                eta(1) = etaVecS(1) + rho*cos(th)
                eta(2) = etaVecS(2) + rho*sin(th)
                temp_T = eta
                
               !map eta to zeta plane
                CALL MatVecProd_real(tVec_1,2,temp_T)
                zeta = temp_T
                
                !Error detection mechanism
                IF ((zeta(1)<0.0_dp).OR.(zeta(2)<0.0_dp).OR.(zeta(1)>1.0_dp).OR.(zeta(2)>1.0_dp)) THEN
                    PRINT*,"zeta = ", zeta
                    PRINT*,"vcsio = ", vcsio
                    PRINT*,"Error detected in zeta variable"
                    STOP
                ENDIF                
                
               !unit normal evaluated at source point in the zeta plane
                CALL calcUnitNormal(zeta(1),zeta(2),nnpt,rNodes_field,nCap)
                
               !corresponding local distance
                CALL mapGaussPtTo3Ddomain(zeta(1),zeta(2),nnpt,rNodes_field,vipps)
                CALL mapGaussPtTo3Ddomain(vcsio(1),vcsio(2),nnpt,rNodes_field,vippo)
                R_bar = vippo-vipps
                rrr = SQRT(DOT_PRODUCT(R_bar,R_bar))
                nCapDotR = DOT_PRODUCT(nCap,R_bar)
                
                !basis function, Lij
                CALL basis_function(zeta,nipp,coeffs,i,Lij)         
                
               !singular integral kernel evaluated at source point
                F_rho_th = nCapDotR*exp(cjvk*rrr)*(1.0_dp-cjvk*rrr)&
                                        /(4.0_dp*pid*rrr**3)*Lij*rho
                
               ! I(1) summation
                z11 = z11 + wwu(uu)*wwv(vv)*F_rho_th*rhoCap_th*thetaBarPrime*Det_T_1 
                
            ENDDO
        ENDDO
        
        z00 = z00 + z11        !Total singular integral value
    ENDDO
    
    RETURN
  
    
  CONTAINS
 
        SUBROUTINE computeParams_IGM(U1,U2,Det_T_1,Z_alpha,etaVecS,&
                                 h_alpha,tVec_1,beta,rhoCap,etaVec2)
            !Subroutine to find the point on the closest edge to the projected field point
            
            IMPLICIT NONE
    
            REAL(kind=dp), INTENT(IN) :: U1(3),U2(3)    
            REAL(kind=dp), INTENT(OUT)::Z_alpha(2,3)
            REAL(kind=dp), INTENT(OUT)::etaVecS(2),tVec_1(2,2)
            REAL(kind=dp), INTENT(OUT)::beta(3),h_alpha(3),Det_T_1
            REAL(kind=dp), INTENT(OUT)::rhoCap(3),etaVec2(2)
    
            INTEGER :: aa
            REAL(kind=dp)::Tvec(2,2),temp_T(2)
            REAL(kind=dp)::Rs1(3),Rs2(3),Rs3(3),R31(3),R23(3),R12(3)
            REAL(kind=dp)::etaCap(3,2),nCap(3),sigma(2)
            REAL(kind=dp)::mCap(3,3),thetaBar(2,3),temp(3),lambda,gamma
            
           !compute eta_2
            lambda = SQRT(DOT_PRODUCT(U1,U1))/SQRT(DOT_PRODUCT(U2,U2))
            gamma = ACOS(DOT_PRODUCT(U1,U2)/(SQRT(DOT_PRODUCT(U1,U1))&
                         *SQRT(DOT_PRODUCT(U2,U2))))
            etaVec2 = (/(COS(gamma)/lambda), (SIN(gamma)/lambda)/)
    
           !Define transformation jacobian from zeta to eta (i.e., Tvec) and vice versa, (i.e., tvec_1)
            temp_T = (/-etaVec2(1)/etaVec2(2), 1.0_dp/etaVec2(2)/)
            tVec_1 = RESHAPE((/1.0_dp, 0.0_dp,temp_T(1),temp_T(2)/),(/2,2/))
            Tvec = RESHAPE((/1.0_dp,0.0_dp,etaVec2(1),etaVec2(2)/),(/2,2/))
    
           !image of singular point in the eta plane.
            temp_T = vcsio      !singular point in the zeta plane
            CALL MatVecProd_real(Tvec,2,temp_T)
            etaVecS = temp_T
    
           !local unit vectors in the eta domain
            etaCap(:,1) = (/1.0_dp,0.0_dp,0.0_dp/)
            etaCap(:,2) = (/0.0_dp,1.0_dp,0.0_dp/)
            nCap = (/0.0_dp,0.0_dp,1.0_dp/)
            
           !Distance vectors
            Rs1 = (1.0_dp-etaVecS(1))*etaCap(:,1) - etaVecS(2)*etaCap(:,2)
            Rs2 = (etaVec2(1)-etaVecS(1))*etaCap(:,1) +&
                    (etaVec2(2)-etaVecS(2))*etaCap(:,2)
            Rs3 = -etaVecS(1)*etaCap(:,1) - etaVecS(2)*etaCap(:,2)
            R12 = (etaVec2(1)-1.0_dp)*etaCap(:,1) + etaVec2(2)*etaCap(:,2)
            R23 = -etaVec2(1)*etaCap(:,1) - etaVec2(2)*etaCap(:,2)
            R31 = etaCap(:,1)
    
            rhoCap=(/SQRT(DOT_PRODUCT(Rs3,Rs3)),SQRT(DOT_PRODUCT(Rs1,Rs1)),&
                        SQRT(DOT_PRODUCT(Rs2,Rs2))/)
    
            CALL real_cross_product(R23,nCap,temp)
            mCap(:,1) = temp/SQRT(DOT_PRODUCT(temp,temp))
    
            CALL real_cross_product(R31,nCap,temp)
            mCap(:,2) = temp/SQRT(DOT_PRODUCT(temp,temp))
    
            CALL real_cross_product(R12,nCap,temp)
            mCap(:,3) = temp/SQRT(DOT_PRODUCT(temp,temp))
            
            h_alpha(1) = DOT_PRODUCT(Rs3,mCap(:,1))
            h_alpha(2) = DOT_PRODUCT(Rs1,mCap(:,2))
            h_alpha(3) = DOT_PRODUCT(Rs2,mCap(:,3))
    
            CALL calcAngles(mCap,rhoCap,h_alpha,etaVecS,etaVec2,Rs1,Rs2,Rs3,R31,R23,R12,beta,thetaBar)
            
            DO aa = 1,3
                sigma(:) = (thetaBar(:,aa)+0.5_dp*pid)/pid
                
               !m = 2 
                !DO ff = 1,2
                !    IF (sigma(ff).EQ.0.5_dp) THEN
                !        Z_alpha(ff,aa) = 0.5_dp
                !    ELSE
                !        Z_alpha(ff,aa) = (sigma(ff)-sqrt(sigma(ff)-&
                !               sigma(ff)**2))/(2.0_dp*sigma(ff)-1.0_dp)                       
                !    ENDIF
                !ENDDO
                
               !m = 3
                Z_alpha(:,aa) = (sigma**3-sigma**2)**(1.0_dp/3.0_dp) + &
                    (sigma**2-sigma)/((sigma**3-sigma**2)**(1.0_dp/3.0_dp))&
                    + sigma    
                
            ENDDO
            Det_T_1 = tVec_1(1,1)*tVec_1(2,2)-tVec_1(1,2)*tVec_1(2,1)
            
        END SUBROUTINE computeParams_IGM 
                                 
        
                         
      SUBROUTINE calcAngles(mCap,rhoCap,h_alpha,etaVecS,etaVec2,Rs1,Rs2,Rs3,R31,R23,R12,beta,thetaBar)
       
            !Subroutine to find the point on the closest edge to the projected field point
            
            IMPLICIT NONE
    
            REAL(kind=dp), INTENT(IN) :: h_alpha(3),mCap(3,3),rhoCap(3),etaVec2(2),etaVecS(2)
            REAL(kind=dp), INTENT(IN) :: Rs1(3),Rs2(3),Rs3(3),R31(3),R23(3),R12(3)
            REAL(kind=dp), INTENT(OUT):: beta(3),thetaBar(2,3)
    
            REAL(kind=dp)::theta_tot(3),psi(3),area(3),ss,area_tot   
            REAL(kind=dp)::len_Rs1,len_Rs2,len_Rs3,len_R31,len_R23,len_R12
            
           !Lengths of the distance vectors in eta plane.
            len_Rs1 = SQRT(DOT_PRODUCT(Rs1,Rs1))
            len_Rs2 = SQRT(DOT_PRODUCT(Rs2,Rs2))
            len_Rs3 = SQRT(DOT_PRODUCT(Rs3,Rs3))
            len_R31 = SQRT(DOT_PRODUCT(R31,R31))
            len_R23 = SQRT(DOT_PRODUCT(R23,R23))
            len_R12 = SQRT(DOT_PRODUCT(R12,R12))
            
           !Angles subtended at the singular image point by each sub-triangle
            theta_tot(1) = ACOS((len_Rs2**2+len_Rs3**2-len_R23**2)/(2.0_dp*len_Rs2*len_Rs3))
            theta_tot(2) = ACOS((len_Rs1**2+len_Rs3**2-len_R31**2)/(2.0_dp*len_Rs1*len_Rs3))
            theta_tot(3) = ACOS((len_Rs1**2+len_Rs2**2-len_R12**2)/(2.0_dp*len_Rs1*len_Rs2))
            
           !
            psi(1) = ACOS((len_R31**2+len_Rs1**2-len_Rs3**2)/(2.0_dp*len_R31*len_Rs1))
            psi(2) = ACOS((len_Rs2**2+len_R23**2-len_Rs3**2)/(2.0_dp*len_Rs2*len_R23))
            psi(3) = ACOS((len_R12**2+len_Rs2**2-len_Rs1**2)/(2.0_dp*len_R12*len_Rs2))
            
            thetaBar(:,1) = (/ACOS(h_alpha(1)/rhoCap(1)),&
                             -ACOS(h_alpha(1)/rhoCap(3))/)
            thetaBar(:,2) = (/ACOS(h_alpha(2)/rhoCap(2)),&
                             -ACOS(h_alpha(2)/rhoCap(1))/)
            thetaBar(:,3) = (/ACOS(h_alpha(3)/rhoCap(3)),&
                             -ACOS(h_alpha(3)/rhoCap(2))/)
            beta(2) = 3.0_dp/2.0_dp*pid
            
            !Calculation of beta and thetaBar based on eta-plane triangle type.
            IF (etaVec2(1)>1) THEN
                
                beta(3) = -ACOS(DOT_PRODUCT((R31/SQRT(DOT_PRODUCT(R31,R31))),&
                            mCap(:,3)))
                
                 IF (theta_tot(3).LE.abs(thetaBar(1,3))) THEN
                    thetaBar(2,3) = thetaBar(1,3) -  theta_tot(3)
                    thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                    thetaBar(2,2) = -(theta_tot(2) - thetaBar(1,2))
                    beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                    !PRINT*,"triangle-type 1a  "
                    
                ELSEIF (theta_tot(2).LE.abs(thetaBar(2,2))) THEN
                    thetaBar(1,2) = thetaBar(2,2) + theta_tot(2)
                    beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                    !PRINT*,"triangle-type 1b  "
                    
                ELSE
                    thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3)) 
                    thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                    thetaBar(2,2) = -(theta_tot(2) - thetaBar(1,2))
                    beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                    !PRINT*,"triangle-type 1c  "
                ENDIF
                
            ELSEIF (etaVec2(1)<0) THEN  
                thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3))
                thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                beta(3) = ACOS(DOT_PRODUCT((R31/SQRT(DOT_PRODUCT(R31,R31))),&
                            mCap(:,3)))
                
                 IF (theta_tot(1).LE.abs(thetaBar(2,1))) THEN
                    thetaBar(1,1) = theta_tot(1) - abs(thetaBar(2,1)) 
                    thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3))
                    thetaBar(2,2) = -(theta_tot(2) - thetaBar(1,2))
                    beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                    !PRINT*,"triangle-type 2a  "
                    
                ELSEIF (theta_tot(2).LE.abs(thetaBar(1,2))) THEN
                    thetaBar(2,2) = thetaBar(1,2) - theta_tot(2)
                    thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                    thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3))                    
                    beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                    !PRINT*,"triangle-type 2b  " 
                    
                ELSE
                    thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3)) 
                    thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                    thetaBar(2,2) = -(theta_tot(2) - thetaBar(1,2))
                    beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                    !PRINT*,"triangle-type 2c  "
                ENDIF


            ELSE
                IF (len_R23.GT.len_R12) THEN  
                    beta(3) = ACOS(DOT_PRODUCT((R31/SQRT(DOT_PRODUCT(R31,R31))),&
                                mCap(:,3)))
                
                    thetaBar(2,2) = -(theta_tot(2) - thetaBar(1,2))
                    
                    IF (etaVec2(2)<1) THEN
                        
                        IF (theta_tot(3).LE.abs(thetaBar(2,3))) THEN
                           thetaBar(1,3) = thetaBar(2,3) + theta_tot(3)
                           thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                           beta(1) = beta(3) + ACOS(DOT_PRODUCT(mCap(:,1),mCap(:,3))) 
                           !PRINT*,"triangle-type 3a  "
                           
                        ELSEIF (theta_tot(1).LE.abs(thetaBar(1,1))) THEN
                           thetaBar(2,1) = thetaBar(1,1) - theta_tot(1)
                           beta(1) = beta(3) + thetaBar(1,3) - abs(thetaBar(2,1))
                           !PRINT*,"triangle-type 3b  " 
                           
                        ELSE
                           thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3)) 
                           thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                           beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                           !PRINT*,"triangle-type 3c  "
                           
                        ENDIF
                    ELSE
                        thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3)) 
                        thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                        beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                        !PRINT*,"triangle-type 3d  "
                    ENDIF
                    
                        
                ELSEIF (len_R23.LT.len_R12) THEN
                    
                    beta(3) = ACOS(DOT_PRODUCT((R31/SQRT(DOT_PRODUCT(R31,R31))),&
                                mCap(:,3)))
                
                    thetaBar(2,2) = -(theta_tot(2) - thetaBar(1,2))
                    
                    IF (etaVec2(2)<1) THEN
                        
                        IF (theta_tot(3).LE.abs(thetaBar(2,3))) THEN
                           thetaBar(1,3) = thetaBar(2,3) + theta_tot(3)
                           thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                           beta(1) = beta(3) + ACOS(DOT_PRODUCT(mCap(:,1),mCap(:,3))) 
                           !PRINT*,"triangle-type 4a  "
                           
                        ELSEIF (theta_tot(1).LE.abs(thetaBar(1,1))) THEN
                           thetaBar(2,1) = thetaBar(1,1) - theta_tot(1)
                           beta(1) = beta(3) + thetaBar(1,3) - abs(thetaBar(2,1))
                           !PRINT*,"triangle-type 4b  " 
                           
                        ELSE
                           thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3)) 
                           thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                           beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                           !PRINT*,"triangle-type 4c  "
                        ENDIF
                        
                    ELSE
                        thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3)) 
                        thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                        beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                        !PRINT*,"triangle-type 4d  "
                        
                    ENDIF
                ELSE
                    beta(3) = ACOS(DOT_PRODUCT((R31/SQRT(DOT_PRODUCT(R31,R31))),&
                                mCap(:,3)))
                
                    thetaBar(2,2) = -(theta_tot(2) - thetaBar(1,2))
                    
                    IF (etaVec2(2)<1) THEN
                        
                        IF (theta_tot(3).LE.abs(thetaBar(2,3))) THEN
                           thetaBar(1,3) = thetaBar(2,3) + theta_tot(3)
                           thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                           beta(1) = beta(3) + ACOS(DOT_PRODUCT(mCap(:,1),mCap(:,3))) 
                           !PRINT*,"triangle-type 5a  "
                        ELSEIF (theta_tot(1).LE.abs(thetaBar(1,1))) THEN
                           thetaBar(2,1) = thetaBar(1,1) - theta_tot(1)
                           beta(1) = beta(3) + thetaBar(1,3) - abs(thetaBar(2,1))
                           !PRINT*,"triangle-type 5b  " 
                        ELSE
                           thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3)) 
                           thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                           beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                           !PRINT*,"triangle-type 5c  "
                        ENDIF
                    ELSE
                        thetaBar(2,3) = -(theta_tot(3) - thetaBar(1,3)) 
                        thetaBar(2,1) = -(theta_tot(1) - thetaBar(1,1))
                        beta(1) = beta(3) + thetaBar(1,3) + abs(thetaBar(2,1))
                        !PRINT*,"triangle-type 5d  "
                    ENDIF
                    
                ENDIF
            ENDIF
            
           !Areas of each sub-triangles and area of eta-plane triangle
            area(1) = 0.5_dp*len_Rs2*len_Rs3*sin(theta_tot(1))
            area(2) = 0.5_dp*len_Rs1*len_Rs3*sin(theta_tot(2))
            area(3) = 0.5_dp*len_Rs1*len_Rs2*sin(theta_tot(3))
            area_tot = 0.5_dp*etaVec2(2)*len_R31
            
            !Error-detection mechanism
            IF (abs(area_tot-sum(area)).GT.1.0E-10_dp) THEN
                PRINT*,"An overlap of subtiangles detected."
                PRINT*,"Total Area = ",area_tot
                PRINT*,"sub_area1+sub_area2+sub_area3  = ",sum(area)
                PRINT*,"Error detected in 'calcAngles'!!."
                STOP
            ENDIF                
            
    END SUBROUTINE calcAngles 
      
END SUBROUTINE ImprovedGuiggianiMethod   
!*****************************************************************************************END                
    
    

!***********************************************************************************************    
SUBROUTINE calcJacobian(zeta,eta,numNodes,rNodes,J) 
    !Subroutine for calculating the Jacobian and unit normal vector, respectively, at area coordinates:(zeta,eta) 
    USE global_com,ONLY:dp
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: numNodes
    REAL(kind=dp), INTENT(OUT) ::   J
    REAL(kind=dp),INTENT(IN) :: zeta, eta, rNodes(3,numNodes)
    REAL(kind=dp) :: dr_bar1(3), dr_bar2(3), temp_J(3)
    INTEGER :: i, wrt1=1, wrt2=2, derivativeDeg=1

    CALL calcPosVectDerivative(derivativeDeg,wrt1,zeta,eta,numNodes,rNodes,dr_bar1)
    CALL calcPosVectDerivative(derivativeDeg,wrt2,zeta,eta,numNodes,rNodes,dr_bar2)
    
    CALL real_cross_product(dr_bar1,dr_bar2,temp_J)  !cross-product of edge vectors
    J = SQRT(DOT_PRODUCT(temp_J,temp_J))             ! Jacobian
    
END SUBROUTINE calcJacobian   
!******************************************************************************************END         
    
 
 !***********************************************************************************************    
SUBROUTINE calcUnitVec(zeta,eta,numNodes,rNodes,dr_bar1,dr_bar2) 
    !Subroutine for calculating the Jacobian and unit normal vector, respectively, at area coordinates:(zeta,eta) 
    USE global_com,ONLY:dp
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: numNodes
    REAL(kind=dp), INTENT(OUT) ::   dr_bar1(3), dr_bar2(3)
    REAL(kind=dp),INTENT(IN) :: zeta, eta, rNodes(3,numNodes)
    REAL(kind=dp) ::  temp_J(3)
    INTEGER :: i, wrt1=1, wrt2=2, derivativeDeg=1

    CALL calcPosVectDerivative(derivativeDeg,wrt1,zeta,eta,numNodes,rNodes,dr_bar1)
    CALL calcPosVectDerivative(derivativeDeg,wrt2,zeta,eta,numNodes,rNodes,dr_bar2)
    
END SUBROUTINE calcUnitVec   
!******************************************************************************************END     
    
    

 !***********************************************************************************************    
SUBROUTINE calcUnitNormal(zeta,eta,numNodes,rNodes,N) 
    !Subroutine for calculating the Jacobian and unit normal vector, respectively, at area coordinates:(zeta,eta) 
    USE global_com,ONLY:dp
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: numNodes
    REAL(kind=dp), INTENT(OUT) ::   N(3)
    REAL(kind=dp),INTENT(IN) :: zeta, eta, rNodes(3,numNodes)
    REAL(kind=dp) :: dr_bar1(3), dr_bar2(3), temp_J(3),J
    INTEGER :: i, wrt1=1, wrt2=2, derivativeDeg=1

    CALL calcPosVectDerivative(derivativeDeg,wrt1,zeta,eta,numNodes,rNodes,dr_bar1)
    CALL calcPosVectDerivative(derivativeDeg,wrt2,zeta,eta,numNodes,rNodes,dr_bar2)
    
    CALL real_cross_product(dr_bar1,dr_bar2,temp_J)  !cross-product of edge vectors
    J = SQRT(DOT_PRODUCT(temp_J,temp_J))             ! Jacobian
    
    N = temp_J/J
    
END SUBROUTINE calcUnitNormal   
!******************************************************************************************END     

    
    
!***********************************************************************************************    
SUBROUTINE calcPosVectDerivative(derivativeDeg,wrt,zeta,eta,numNodes,rNodes,dr_bar)
    !Subroutine for interpolating derivatives of r in terms of given nodes rNodes
    ! triangular patch given the area coordinates:(zeta,eta) 
    USE global_com,ONLY:dp
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: numNodes, wrt, derivativeDeg
    REAL(kind=dp), INTENT(OUT) :: dr_bar(3)  
    REAL(kind=dp),INTENT(IN) :: zeta, eta, rNodes(3,numNodes)
    REAL(kind=dp) :: dN(numNodes) 
    INTEGER :: i
    
!-------------Determine the degree and type of derivative --------------------------
    IF (wrt == 1) THEN          !derivative w.r.t zeta
        
        !First derivatives w.r.t zeta of shape functions 
         INNER1:IF (derivativeDeg==1) THEN
            INNER2:IF (numNodes==6)  THEN
                        dN(1) = 4.0_dp*zeta-1.0_dp
                        dN(2) = 0.0_dp
                        dN(3) = -3.0_dp+4.0_dp*zeta+4.0_dp*eta
                        dN(4) = 4.0_dp*eta
                        dN(5) = -4.0_dp*eta
                        dN(6) = 4.0_dp*(1.0_dp-2.0_dp*zeta-eta)
                    ELSEIF (numNodes==3)  THEN
                        dN(1) = 1.0_dp
                        dN(2) = 0.0_dp
                        dN(3) = -1.0_dp
                    ELSE 
                        PRINT*,"1st derivatives for order > 3 or < 1, not supported for now."
                        STOP
                    ENDIF INNER2
            
               !Second derivatives w.r.t zeta of shape functions   
         ELSEIF (derivativeDeg==2) THEN
                INNER3:IF (numNodes==6)  THEN             
                            dN(1) = 4.0_dp
                            dN(2) = 0.0_dp
                            dN(3) = 4.0_dp
                            dN(4) = 0.0_dp
                            dN(5) = 0.0_dp
                            dN(6) = -8.0_dp
                        ELSEIF (numNodes==3)  THEN
                            dN(1) = 0.0_dp
                            dN(2) = 0.0_dp
                            dN(3) = 0.0_dp
                        ELSE 
                            PRINT*,"2nd derivatives for order > 3 or < 1, not supported for now."
                            STOP
                        ENDIF INNER3
                ELSE 
                    PRINT*,"Derivative degrees must be <= 3."
                    STOP
                ENDIF INNER1

    ELSEIF (wrt == 2) THEN      !derivative w.r.t eta
        
        !First derivatives w.r.t eta of shape functions 
         INNERa:IF (derivativeDeg==1) THEN
             INNERb:IF (numNodes==6)  THEN             
                        dN(1) = 0.0_dp
                        dN(2) = 4.0_dp*eta-1.0_dp
                        dN(3) = -3.0_dp+4.0_dp*zeta+4.0_dp*eta
                        dN(4) = 4.0_dp*zeta
                        dN(5) = 4.0_dp*(1.0_dp-zeta-2.0_dp*eta)
                        dN(6) = -4.0_dp*zeta
                    ELSEIF (numNodes==3)  THEN
                        dN(1) = 0.0_dp
                        dN(2) = 1.0_dp
                        dN(3) = -1.0_dp
                    ELSE 
                        PRINT*,"1st derivatives for order > 3 or < 1, not supported for now."
                        STOP
                    ENDIF INNERb
            
               !Second derivatives w.r.t eta of shape functions   
                ELSEIF (derivativeDeg==2) THEN
                 INNERc:IF (numNodes==6)  THEN
                            dN(1) = 0.0_dp
                            dN(2) = 4.0_dp
                            dN(3) = 4.0_dp
                            dN(4) = 0.0_dp
                            dN(5) = -8.0_dp
                            dN(6) = 0.0_dp
                        ELSEIF (numNodes==3)  THEN
                            dN(1) = 0.0_dp
                            dN(2) = 0.0_dp
                            dN(3) = 0.0_dp
                        ELSE 
                            PRINT*,"2nd derivatives for order > 3 or < 1, not supported for now."
                            STOP
                        ENDIF INNERc
                ELSE 

                    PRINT*,"Derivative degrees must be <= 3."
                    STOP
                ENDIF INNERa
                
    ELSE                !2nd derivative w.r.t zeta and eta
           IF (numNodes==6)  THEN
                dN(1) = 0.0_dp
                dN(2) = 0.0_dp
                dN(3) = 4.0_dp
                dN(4) = 4.0_dp
                dN(5) = -4.0_dp
                dN(6) = -4.0_dp
            ELSEIF (numNodes==3)  THEN
                dN(1) = 0.0_dp
                dN(2) = 0.0_dp
                dN(3) = 0.0_dp
            ELSE 
                PRINT*,"2nd derivatives for order > 3 or < 1, not supported for now."
                STOP
            ENDIF
    ENDIF
   
!-----------------Compute position-vector derivative ---------------------------------
    dr_bar(:) = 0.0_dp         !initialize variable
    
    DO i = 1,numNodes 
        dr_bar(:) = dr_bar(:) + dN(i)*rNodes(:,i)
    ENDDO

END SUBROUTINE calcPosVectDerivative 
!******************************************************************************************END    
    
    

!***********************************************************************************************    
SUBROUTINE mapGaussPtTo3DDOmain(zetaSrc,etaSrc,numNodes,rNodes,r)
    !Subroutine for interpolating r in terms of given nodes rj
    ! triangular patch given the area coordinates:(zetaSrc,etaSrc) 
    USE global_com,ONLY:dp
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: numNodes
    REAL(kind=dp), INTENT(OUT) :: r(3)
    REAL(kind=dp),INTENT(IN) :: zetaSrc, etaSrc, rNodes(3,numNodes)
    REAL(kind=dp) :: shape_func(numNodes) 
    INTEGER :: i
    
    CALL shape_function(zetaSrc,etaSrc,numNodes,shape_func)
    
    r(:) = 0.0_dp
    DO i=1,numNodes
        r(:) = r(:) + shape_func(i)*rNodes(:,i)
    ENDDO
    
END SUBROUTINE mapGaussPtTo3DDOmain    
!******************************************************************************************END    


!**********************************************************************************************    
SUBROUTINE shape_function(zeta,eta,numNodes,shape_func)
    !Subroutine for implementing silvester-Lagrange shape function for a 1,3, & 6-node 
    ! triangular patch given an area coordinate:(zetaSrc,eta) 
    USE global_com,ONLY:dp
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: numNodes
    REAL(kind=dp), INTENT(OUT) :: shape_func(numNodes)
    REAL(kind=dp),INTENT(IN) :: zeta, eta
    REAL(kind=dp) :: AreaCoords(3)
    
    AreaCoords(1) = zeta
    AreaCoords(2) = eta
    AreaCoords(3) = 1.0_dp-zeta-eta
    
    IF (numNodes==1) THEN
        shape_func(1) = 1.0_dp
        
    ELSEIF (numNodes==3) THEN
        shape_func(1) = AreaCoords(1)
        shape_func(2) = AreaCoords(2)
        shape_func(3) = AreaCoords(3)
        
    ELSEIF (numNodes==6) THEN
        shape_func(1) = AreaCoords(1)*(2.0_dp*AreaCoords(1)-1.0_dp)
        shape_func(2) = AreaCoords(2)*(2.0_dp*AreaCoords(2)-1.0_dp)
        shape_func(3) = AreaCoords(3)*(2.0_dp*AreaCoords(3)-1.0_dp)
        shape_func(4) = 4.0_dp*AreaCoords(1)*AreaCoords(2)
        shape_func(5) = 4.0_dp*AreaCoords(2)*AreaCoords(3)
        shape_func(6) = 4.0_dp*AreaCoords(1)*AreaCoords(3)
        
    ELSEIF (numNodes==10) THEN
        shape_func(1) = 0.5_dp*AreaCoords(1)*(3.0_dp*AreaCoords(1)-1.0_dp)*(3.0*AreaCoords(1)-2.0_dp)
        shape_func(2) = 0.5_dp*AreaCoords(2)*(3.0_dp*AreaCoords(2)-1.0_dp)*(3.0*AreaCoords(2)-2.0_dp)
        shape_func(3) = 0.5_dp*AreaCoords(3)*(3.0_dp*AreaCoords(3)-1.0_dp)*(3.0*AreaCoords(3)-2.0_dp)
        shape_func(4) = 9.0_dp/2.0_dp*AreaCoords(2)*AreaCoords(1)*(3.0_dp*AreaCoords(1)-1.0_dp)
        shape_func(5) = 9.0_dp/2.0_dp*AreaCoords(1)*AreaCoords(2)*(3.0_dp*AreaCoords(2)-1.0_dp)
        shape_func(6) = 9.0_dp/2.0_dp*AreaCoords(3)*AreaCoords(2)*(3.0_dp*AreaCoords(2)-1.0_dp)
        shape_func(7) = 9.0_dp/2.0_dp*AreaCoords(2)*AreaCoords(3)*(3.0_dp*AreaCoords(3)-1.0_dp)
        shape_func(8) = 9.0_dp/2.0_dp*AreaCoords(1)*AreaCoords(3)*(3.0_dp*AreaCoords(3)-1.0_dp)
        shape_func(9) = 9.0_dp/2.0_dp*AreaCoords(3)*AreaCoords(1)*(3.0_dp*AreaCoords(1)-1.0_dp)
        shape_func(10) = 27.0_dp*AreaCoords(1)*AreaCoords(2)*AreaCoords(3)

    
    ELSE
        PRINT*,"Shape function undefined for order > 3 for now"
        STOP
    ENDIF
    
    END SUBROUTINE shape_function
!******************************************************************************************END    

 
!***********************************************************************************************    
SUBROUTINE basis_function(zeta,numNodes,coeffs,i,basisFunc)
    !Subroutine for interpolating the grid-robust Lagrange basis functions at any desired (interior) point
    ! on a patch with area coordinates: (zeta,eta). (See G. Kang et al, 2001).
 
    USE global_com,ONLY:dp
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: numNodes,i
    REAL(kind=dp), INTENT(OUT) :: basisFunc
    REAL(kind=dp), INTENT(IN) :: zeta(2),coeffs(numNodes,numNodes)
    INTEGER :: k
    REAL(kind=dp) :: polyvector(numNodes)


    basisFunc=0.0_dp
    IF(numNodes.eq.6) THEN
        polyvector(1)=1.0_dp
        polyvector(2)=zeta(1)
        polyvector(3)=zeta(2)
        polyvector(4)=zeta(1)*zeta(2)
        polyvector(5)=zeta(1)*zeta(1)
        polyvector(6)=zeta(2)*zeta(2)
        DO k=1,6
            basisFunc=basisFunc+coeffs(i,k)*polyvector(k)
        ENDDO
    ELSEIF(numNodes.eq.12) THEN
        polyvector(1)=1.0_dp
        polyvector(2)=zeta(1)
        polyvector(3)=zeta(2)
        polyvector(4)=zeta(1)*zeta(1)
        polyvector(5)=zeta(1)*zeta(2)
        polyvector(6)=zeta(2)*zeta(2)
        polyvector(7)=zeta(1)*zeta(1)*zeta(1)
        polyvector(8)=zeta(1)*zeta(1)*zeta(2)
        polyvector(9)=zeta(1)*zeta(2)*zeta(2)
        polyvector(10)=zeta(2)*zeta(2)*zeta(2)
        polyvector(11)=zeta(1)*zeta(1)*zeta(1)*zeta(2)
        polyvector(12)=zeta(1)*zeta(2)*zeta(2)*zeta(2)
        DO k=1,12
            basisFunc=basisFunc+coeffs(i,k)*polyvector(k)
        ENDDO
     ELSEIF(numNodes.eq.3) THEN
        polyvector(1)=1.0_dp
        polyvector(2)=zeta(1)
        polyvector(3)=zeta(2)
        DO k=1,3
           basisFunc=basisFunc+coeffs(i,k)*polyvector(k)
        ENDDO
      ELSEIF(numNodes.eq.1) THEN
         polyvector(1)=1.0_dp
         DO k=1,1
            basisFunc=basisFunc+coeffs(i,k)*polyvector(k)
         ENDDO
      ELSE
         PRINT*,"can not set polynomial vector for higher order now"
         STOP
      ENDIF
                       
END SUBROUTINE basis_function   
!******************************************************************************************END    
    
    
!***********************************************************************************************    
SUBROUTINE findPolyCoeffs(alphas,betas,gammas,numNodes,coeffs) 
    USE global_com,ONLY:dp
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: numNodes
    REAL(kind=dp), INTENT(OUT) :: coeffs(numNodes,numNodes)
    REAL(kind=dp), INTENT(IN) :: alphas(numNodes),betas(numNodes),gammas(numNodes)
    REAL(kind=dp) :: vcsis(2),vos1(2),vos2(2),vos3(2)
    REAL(kind=dp) :: polymatrix(numNodes,numNodes),idMat(numNodes,numNodes)
    INTEGER :: IPIV(numNodes),INFO,k,i
    
    vos1=(/1.0_dp,0.0_dp/)  !right angle triangle vertex points in (u,v) space	
    vos2=(/0.0_dp,1.0_dp/)  !the first number is u value, the second number is v value	
    vos3=(/0.0_dp,0.0_dp/)

    
    IF(numNodes.eq.6) THEN
         DO i=1,6
            vcsis= alphas(i)*vos1+betas(i)*vos2+gammas(i)*vos3 
            polymatrix(1,i)=1.0_dp
            polymatrix(2,i)=vcsis(1)
            polymatrix(3,i)=vcsis(2)
            polymatrix(4,i)=vcsis(1)*vcsis(2)
            polymatrix(5,i)=vcsis(1)*vcsis(1)
            polymatrix(6,i)=vcsis(2)*vcsis(2)
         ENDDO

         idMat(1:6,1:6)=0.0_dp
         DO k=1,6
            idMat(k,k)=1.0_dp
         ENDDO

         CALL dGESV(6,6,polymatrix,6,IPIV,idMat,6,INFO)

      ELSEIF(numNodes.eq.12) THEN
         DO i=1,12
            vcsis= alphas(i)*vos1+betas(i)*vos2+gammas(i)*vos3 
            polymatrix(1,i)=1.0_dp
            polymatrix(2,i)=vcsis(1)
            polymatrix(3,i)=vcsis(2)
            polymatrix(4,i)=vcsis(1)*vcsis(1)
            polymatrix(5,i)=vcsis(1)*vcsis(2)
            polymatrix(6,i)=vcsis(2)*vcsis(2)
            polymatrix(7,i)=vcsis(1)*vcsis(1)*vcsis(1)
            polymatrix(8,i)=vcsis(1)*vcsis(1)*vcsis(2)
            polymatrix(9,i)=vcsis(1)*vcsis(2)*vcsis(2)
            polymatrix(10,i)=vcsis(2)*vcsis(2)*vcsis(2)
            polymatrix(11,i)=vcsis(1)*vcsis(1)*vcsis(1)*vcsis(2)
            polymatrix(12,i)=vcsis(1)*vcsis(2)*vcsis(2)*vcsis(2)
         ENDDO

         idMat(1:12,1:12)=0.0_dp
         DO k=1,12
            idMat(k,k)=1.0_dp
         ENDDO

         CALL dGESV(12,12,polymatrix,12,IPIV,idMat,12,INFO)

      ELSEIF(numNodes.eq.3) THEN
         DO i=1,3
            vcsis= alphas(i)*vos1+betas(i)*vos2+gammas(i)*vos3 
            polymatrix(1,i)=1.0_dp
            polymatrix(2,i)=vcsis(1)
            polymatrix(3,i)=vcsis(2)
         ENDDO

         idMat(1:3,1:3)=0.0_dp
         DO k=1,3
            idMat(k,k)=1.0_dp
         ENDDO

         CALL dGESV(3,3,polymatrix,3,IPIV,idMat,3,INFO)

      ELSEIF(numNodes.eq.1) THEN
         DO i=1,1
            vcsis= alphas(i)*vos1+betas(i)*vos2+gammas(i)*vos3 
            polymatrix(1,i)=1.0_dp
         ENDDO

         idMat(1:1,1:1)=0.0_dp
         DO k=1,1
            idMat(k,k)=1.0_dp
         ENDDO

         CALL dGESV(1,1,polymatrix,1,IPIV,idMat,1,INFO)

      ELSE
         PRINT*,"can not set polynomial matrix for higher order now"
         STOP
      ENDIF
      
      coeffs(:,:) = idMat(:,:)

END SUBROUTINE findPolyCoeffs    
!******************************************************************************************END    
    
    
 
    
    

    