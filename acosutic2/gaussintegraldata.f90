
 module gaussintegraldata

   use global_com,only: dp

   implicit none      
   real(kind=dp),allocatable:: gammao(:),alphao(:),betao(:),wo(:)
   real(kind=dp),allocatable:: gammas(:),alphas(:),betas(:),ws(:)
   integer ngo,ngs

   !real(kind=dp),allocatable:: gammao11(:),alphao11(:),betao11(:),wo11(:)
   real(kind=dp),allocatable:: gammas11(:),alphas11(:),betas11(:),ws11(:)
   !integer ngo11,ngs11
   integer ngs11 

    !real(kind=dp),allocatable:: gammao22(:),alphao22(:),betao22(:),wo22(:)
    real(kind=dp),allocatable:: gammas22(:),alphas22(:),betas22(:),ws22(:)
    !integer ngo22,ngs22
    integer ngs22

    real(kind=dp),allocatable:: gammao2(:),alphao2(:),betao2(:),wo2(:)
   real(kind=dp),allocatable:: gammas2(:),alphas2(:),betas2(:),ws2(:)
   integer ngo2,ngs2

   real(kind=dp),allocatable:: gammao3(:),alphao3(:),betao3(:),wo3(:)
   real(kind=dp),allocatable:: gammas3(:),alphas3(:),betas3(:),ws3(:)
   integer ngo3,ngs3

   integer nlg,nlgo,nlgs
   real(kind=dp),allocatable,dimension(:):: wwwl,xxxl,wwo,xxo,wws,xxs 

 contains

   !****************************************************
   subroutine setgaussintpara(ng1,  ng2) 
   !****************************************************
     implicit none
     integer::ng1,ng2

     if(allocated(wo).and.ngo.eq.ng1.and.ngs.eq.ng2) return;

     if(allocated(wo))	 then 
        deallocate(wo,gammao,alphao,betao)
        deallocate(ws,gammas,alphas,betas)
     endif

     ngo = ng1
     ngs = ng2

     allocate(wo(ngo),gammao(ngo) , alphao(ngo),betao(ngo))
     allocate(ws(ngs),gammas(ngs) , alphas(ngs),betas(ngs))
		
     wo = 0.0_dp;  gammao = 0.0_dp;  alphao = 0.0_dp;  betao  = 0.0_dp
     ws = 0.0_dp;  gammas = 0.0_dp;  alphas = 0.0_dp;  betas  = 0.0_dp
		
     call ginttrg(ngo, alphao, betao, gammao, wo)
     call ginttrg(ngs, alphas, betas, gammas, ws)
		
   endsubroutine setgaussintpara
       !****************************************************
       !****************************************************
   
   subroutine setgaussintpara11(ng2)
       !****************************************************
     implicit none
     ! integer::ng1,ng2
     integer::ng2
    ! if(allocated(wo11).and.ngo11.eq.ng1.and.ngs11.eq.ng2) return;

     if(allocated(ws11))	 then
       ! deallocate(wo11,gammao11,alphao11,betao11)
        deallocate(ws11,gammas11,alphas11,betas11)
     endif

     !ngo11 = ng1
     ngs11 = ng2

     !allocate(wo11(ngo11),gammao11(ngo11) , alphao11(ngo11),betao11(ngo11))
     allocate(ws11(ngs11),gammas11(ngs11) , alphas11(ngs11),betas11(ngs11))
		
     !wo11 = 0.0_dp;  gammao11 = 0.0_dp;  alphao11 = 0.0_dp;  betao11  = 0.0_dp
     ws11 = 0.0_dp;  gammas11 = 0.0_dp;  alphas11 = 0.0_dp;  betas11  = 0.0_dp
		
     !call ginttrg(ngo11, alphao11, betao11, gammao11, wo11)
     call ginttrg(ngs11, alphas11, betas11, gammas11, ws11)
		
   endsubroutine setgaussintpara11
       !****************************************************
       !****************************************************

!****************************************************

    subroutine setgaussintpara22(  ng2)
!****************************************************
    implicit none
    !integer::ng1,ng2
    integer::ng2
    !if(allocated(wo22).and.ngo22.eq.ng1.and.ngs22.eq.ng2) return;

    if(allocated(ws22))	 then
      !  deallocate(wo22,gammao22,alphao22,betao22)
        deallocate(ws22,gammas22,alphas22,betas22)
    endif

    !ngo22 = ng1
    ngs22 = ng2

    !allocate(wo22(ngo22),gammao22(ngo22) , alphao22(ngo22),betao22(ngo22))
    allocate(ws22(ngs22),gammas22(ngs22) , alphas22(ngs22),betas22(ngs22))

    !wo22 = 0.0_dp;  gammao22 = 0.0_dp;  alphao22 = 0.0_dp;  betao22  = 0.0_dp
    ws22 = 0.0_dp;  gammas22 = 0.0_dp;  alphas22 = 0.0_dp;  betas22  = 0.0_dp

    !call ginttrg(ngo22, alphao22, betao22, gammao22, wo22)
    call ginttrg(ngs22, alphas22, betas22, gammas22, ws22)

    endsubroutine setgaussintpara22
    !****************************************************
    !****************************************************

   subroutine setgaussintpara2(ng1,  ng2)
       !****************************************************
     implicit none
     integer::ng1,ng2

     if(allocated(wo2).and.ngo2.eq.ng1.and.ngs2.eq.ng2) return;

     if(allocated(wo2))	 then 
        deallocate(wo2,gammao2,alphao2,betao2)
        deallocate(ws2,gammas2,alphas2,betas2)
     endif

     ngo2 = ng1
     ngs2 = ng2

     allocate(wo2(ngo2),gammao2(ngo2) , alphao2(ngo2),betao2(ngo2))
     allocate(ws2(ngs2),gammas2(ngs2) , alphas2(ngs2),betas2(ngs2))
		
     wo2 = 0.0_dp;  gammao2 = 0.0_dp;  alphao2 = 0.0_dp;  betao2  = 0.0_dp
     ws2 = 0.0_dp;  gammas2 = 0.0_dp;  alphas2 = 0.0_dp;  betas2  = 0.0_dp
		
     call ginttrg(ngo2, alphao2, betao2, gammao2, wo2)
     call ginttrg(ngs2, alphas2, betas2, gammas2, ws2)
		
   endsubroutine setgaussintpara2
       !****************************************************


   !****************************************************
   subroutine setgaussintpara3(ng1,  ng2) 
   !****************************************************
     implicit none
     integer::ng1,ng2

     if(allocated(wo3).and.ngo3.eq.ng1.and.ngs3.eq.ng2) return;

     if(allocated(wo3))	 then 
        deallocate(wo3,gammao3,alphao3,betao3)
        deallocate(ws3,gammas3,alphas3,betas3)
     endif

     ngo3 = ng1
     ngs3 = ng2

     allocate(wo3(ngo3),gammao3(ngo3) , alphao3(ngo3),betao3(ngo3))
     allocate(ws3(ngs3),gammas3(ngs3) , alphas3(ngs3),betas3(ngs3))
		
     wo3 = 0.0_dp;  gammao3 = 0.0_dp;  alphao3 = 0.0_dp;  betao3  = 0.0_dp
     ws3 = 0.0_dp;  gammas3 = 0.0_dp;  alphas3 = 0.0_dp;  betas3  = 0.0_dp
		
     call ginttrg(ngo3, alphao3, betao3, gammao3, wo3)
     call ginttrg(ngs3, alphas3, betas3, gammas3, ws3)
		
   endsubroutine setgaussintpara3
       !****************************************************



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine setgaussintline(nlg1)
      !******************************************
     implicit none
     integer::nlg1   
     logical b

     b = allocated(wwwl)
     if(b.and.nlg.eq.nlg1) return;
     if(b)	 deallocate(wwwl,xxxl)

     nlg = nlg1
     allocate(wwwl(nlg),xxxl(nlg))
     wwwl = 0.0_dp
     xxxl = 0.0_dp
!       call gauleg(-1.,1.,xxxl,wwwl,nlg)
     call gauleg(0.0_dp,1.0_dp,xxxl,wwwl,nlg)
   endsubroutine setgaussintline


      !******************************************
   subroutine setgaussintline_new(nlg1,nlg2) 	
      !******************************************
     implicit none
     integer::nlg1,nlg2 
     logical b

     b = allocated(wwo)
     if(b.and.nlgo.eq.nlg1.and.nlgs.eq.nlg2) return;
     if(b)	 deallocate(wwo,xxo,wws,xxs)

     nlgo = nlg1
     nlgs = nlg2
     allocate(wwo(nlgo),xxo(nlgo))
     allocate(wws(nlgs),xxs(nlgs))
     wwo = 0.0_dp ; wws = 0.0_dp
     xxo = 0.0_dp ; xxs = 0.0_dp
     call gauleg(0.0_dp,1.0_dp,xxo,wwo,nlgo)
     call gauleg(0.0_dp,1.0_dp,xxs,wws,nlgs)
   endsubroutine setgaussintline_new




!****************************************************************************
   subroutine ginttrg(n,alpha,beta,gamma,ww)
!****************************************************************************
    implicit none
    integer::n
    real(kind=dp):: alpha(n), beta(n),gamma(n),ww(n)
    real(kind=dp)::v1,v2,v3,v4,v5,v6,v7,v8,v9

    alpha = 0.0_dp;	beta  = 0.0_dp
    gamma = 0.0_dp;	ww    = 0.0_dp 

    if(n.eq.1) then ! p =1
       ww(1) = 1.d0
       v1 = 1.d0/3.d0
       alpha(1) = v1; beta(1) = v1;	gamma(1) = v1
    endif
    
    if(n.eq.3) then ! p =2
       ww(1:3) = 1.d0/3.d0
       v1 = 2.d0/3.d0;
       v2 = 1.d0/6.d0;

       alpha(1) = v1; alpha(2) = v2;  alpha(3) = v2;  
       beta(1) = v2;  beta(2) = v1;   beta(3) = v2;    
       gamma(1) = v2; gamma(2) = v2;  gamma(3) = v1;   
    endif
    
    if(n.eq.4) then ! p = 3
       ww(1) = -0.5625d0
       v1 = 1.d0/3.d0
       alpha(1) = v1; 
       beta(1) = v1;	
       gamma(1) = v1
       ww(2:4) = 0.52083333333333333d0

       alpha(2) = 0.6d0; beta(2) = 0.2d0;   gamma(2) = 0.2d0;  
       alpha(3) = 0.2d0; beta(3) = 0.6d0;   gamma(3) = 0.2d0;   
       alpha(4) = 0.2d0; beta(4) = 0.2d0;   gamma(4) = 0.6d0;
    endif


    if(n.eq.6) then ! p =4
       ww(1:3)  = 0.223381589678011d0
       ww(4:6)  = 0.109951743655322d0
       v1= 0.108103018168070d0;
       v3= 0.816847572980459d0;
       v2= 0.445948490915965d0;
       v4= 0.091576213509771d0;
       alpha( 1)  =  v1;      beta( 1)  =  v2;   gamma(1)  = v2
       alpha( 2)  =  v2;      beta( 2)  =  v1;   gamma(2)  = v2      
       alpha( 3)  =  v2;      beta( 3)  =  v2;   gamma(3)  = v1  

       alpha( 4)  =  v3;      beta( 4)  =  v4;   gamma(4)  = v4
       alpha( 5)  =  v4;      beta( 5)  =  v3;   gamma(5)  = v4        
       alpha( 6)  =  v4;      beta( 6)  =  v4;   gamma(6)  = v3      

    endif

    if(n.eq.7) then ! p =5  !see in graglia paper

       ww(1) = 0.225d0
       alpha(1) = 1.d0/3.d0;beta(1) = 1.d0/3.d0;gamma(1) = 1.d0/3.d0;

       ww(2:4) = (155.d0+sqrt(15.0d0))/1200.0d0 !0.132394152788506 !
       ww(5:7) = (155.d0-sqrt(15.0d0))/1200.0d0 !0.125939180544827 !

       v1 = (9.0d0-2.0d0*sqrt(15.0d0))/21.0d0
       v2 = (6.0d0-   sqrt(15.0d0))/21.0d0
       v3 = (6.0d0+   sqrt(15.0d0))/21.0d0
       v4 = (9.0d0+2.0d0*sqrt(15.0d0))/21.0d0

       
       
       alpha(2) = v1;  alpha(3) = v3;  alpha(4) = v3; 
       beta(2) = v3;   beta(3) = v1;   beta(4) = v3; 
       gamma(2) = v3;  gamma(3) = v3;  gamma(4) = v1;

       alpha(5) = v4;  alpha(6) = v2;  alpha(7) = v2;
       beta(5) = v2;   beta(6) = v4;   beta(7) = v2;
       gamma(5) = v2;  gamma(6) = v2;  gamma(7) = v4;
    endif

    !      if(n.eq.9) then ! p =3  !see in graglia paper
    !	  
    !	   ww(1:9) = 1./9.
    !      v1 =  1./18.
    !	v2 =  2./9.
    !	v3 =  7./18.
    !	v4 = 13./18.
    !
    !	alpha(1) = v1;  beta(1) = v2;
    !      alpha(2) = v2;  beta(2) = v1;  
    !	alpha(3) = v1;  beta(3) = v4;
    !	alpha(4) = v4;  beta(4) = v1; 
    !      alpha(5) = v2;  beta(5) = v3;
    !	alpha(6) = v3;  beta(6) = v2;
    !	alpha(7) = v3;  beta(7) = v3;
    !	alpha(8) = v2;  beta(8) = v4;
    !	alpha(9) = v4;  beta(9) = v2;
    !	gamma = 1. - beta -alpha
    !
    !	endif


    if(n.eq.9) then ! p =3  !×ô¼ºíæµ¼£¬±íê¾½«ò»¸öèý½çðî·ö³éµèãæ»ýµä9¸öèý½çðî

       ww(1:9) = 1.0d0/9.0d0
       v1 =  1.0d0/9.0d0
       v2 =  2.0d0/9.0d0
       v3 =  4.0d0/9.0d0
       v4 =  5.0d0/9.0d0
       v5 =  7.0d0/9.0d0

       
       alpha(1) = v1;  beta(1) = v5;
       alpha(2) = v5;  beta(2) = v1;  
       alpha(3) = v1;  beta(3) = v1;

       alpha(4) = v3;  beta(4) = v1; 
       alpha(5) = v1;  beta(5) = v3;
       alpha(6) = v3;  beta(6) = v3;

       alpha(7) = v2;  beta(7) = v4;
       alpha(8) = v2;  beta(8) = v2;
       alpha(9) = v4;  beta(9) = v2;
       gamma = 1. - beta -alpha

    endif
    
    if(n.eq.12) then ! p = 6
       ww(1:3)   =    0.116786275726379d0
       ww(4:6)   =    0.050844906370207d0
       ww(7:12)  =    0.082851075618374d0   !¿éòô¿¼âçóãn0 n1 n2à´¿øöæ

       v1  = 0.501426509658179d0;  v2 = 0.249286745170910d0
       v4  = 0.873821971016996d0;  v5 = 0.063089014491502d0
       v7  = 0.053145049844817d0;  v8 = 0.310352451033784d0
       v9  = 1. - v7 - v8

       alpha( 1)  =  v1;  beta( 1)  =  v2;    gamma( 1)  =  v2;
       alpha( 2)  =  v2;  beta( 2)  =  v1;    gamma( 2)  =  v2;       
       alpha( 3)  =  v2;  beta( 3)  =  v2;    gamma( 3)  =  v1;     

       alpha( 4)  =  v4;  beta( 4)  =  v5;    gamma( 4)  =  v5;
       alpha( 5)  =  v5;  beta( 5)  =  v4;    gamma( 5)  =  v5;         
       alpha( 6)  =  v5;  beta( 6)  =  v5;    gamma( 6)  =  v4;        

       alpha( 7)  =  v7;  beta( 7)  =  v8;    gamma( 7)  =  v9;
       alpha( 8)  =  v7;  beta( 8)  =  v9;    gamma( 8)  =  v8;    
       alpha( 9)  =  v8;  beta( 9)  =  v7;    gamma( 9)  =  v9;      
       alpha(10)  =  v8;  beta(10)  =  v9;    gamma(10)  =  v7;      
       alpha(11)  =  v9;  beta(11)  =  v8;    gamma(11)  =  v7;     
       alpha(12)  =  v9;  beta(12)  =  v7;    gamma(12)  =  v8;    

    endif

    if(n.eq.13) then ! p =7
       ww(1)     =   -0.149570044467682d0
       ww(2:4)   =    0.175615257433208d0
       ww(5:7)   =    0.053347235608838d0
       ww(8:13)  =    0.077113760890257d0   !¿éòô¿¼âçóãn0 n1 n2à´¿øöæ

       alpha(1)  = 0.333333333333333d0;  beta(1) = 0.333333333333333d0
       alpha(2)  = 0.479308067841920d0;	beta(2) = 0.260345966079040d0
       alpha(5)  = 0.869739794195568d0;  beta(5) = 0.065130102902216d0
       alpha(8)  = 0.048690315425316d0;	beta(8) = 0.312865496004874d0

       gamma(1)  = 1. - alpha(1) - beta(1)
       gamma(2)  = 1. - alpha(2) - beta(2)
       gamma(5)  = 1. - alpha(5) - beta(5)
       gamma(8)  = 1. - alpha(8) - beta(8)

       alpha( 3)  =  beta(2);         beta( 3)  = alpha(2)           
       alpha( 4)  =  beta(2);         beta( 4)  =  beta(2)           
       alpha( 6)  =  beta(5);         beta( 6)  = alpha(5)           
       alpha( 7)  =  beta(5);         beta( 7)  =  beta(5)           

       alpha( 9)  = alpha(8);         beta( 9)  = gamma(8)           
       alpha(10)  =  beta(8);         beta(10)  = alpha(8)           
       alpha(11)  =  beta(8);         beta(11)  = gamma(8)           
       alpha(12)  = gamma(8);         beta(12)  =  beta(8)           
       alpha(13)  = gamma(8);         beta(13)  = alpha(8)           

       gamma(3)     = 1. - alpha(3) - beta(3)
       gamma(4)     = 1. - alpha(4) - beta(4)
       gamma(6)     = 1. - alpha(6) - beta(6)
       gamma(7)     = 1. - alpha(7) - beta(7)

       gamma(9 )    = 1. - alpha( 9) - beta( 9)
       gamma(10)    = 1. - alpha(10) - beta(10)
       gamma(11)    = 1. - alpha(11) - beta(11)
       gamma(12)    = 1. - alpha(12) - beta(12)
       gamma(13)    = 1. - alpha(13) - beta(13)
    endif

    if(n.eq.16) then ! p = 8
       ww(1:3)   =  0.095091634267285d0
       ww(4:6)   =  0.103217370534718d0
       ww(7:12)  =  0.027230314174435d0   

       ww(13:15) =  0.032458497623198d0
       ww(16)    =  0.144315607677787d0



       
       v1  = 0.081414823414554d0;  v2 = 0.459292588292723d0
       v4  = 0.658861384496480d0;  v5 = 0.170569307751760d0
       v7  = 0.008394777409958d0;  v8 = 0.263112829634638d0;
       v9  = 1. - v7 - v8

       alpha( 1)  =  v1;  beta( 1)  =  v2;    gamma( 1)  =  v2;
       alpha( 2)  =  v2;  beta( 2)  =  v1;    gamma( 2)  =  v2;       
       alpha( 3)  =  v2;  beta( 3)  =  v2;    gamma( 3)  =  v1;     

       alpha( 4)  =  v4;  beta( 4)  =  v5;    gamma( 4)  =  v5;
       alpha( 5)  =  v5;  beta( 5)  =  v4;    gamma( 5)  =  v5;         
       alpha( 6)  =  v5;  beta( 6)  =  v5;    gamma( 6)  =  v4;        

       alpha( 7)  =  v7;  beta( 7)  =  v8;    gamma( 7)  =  v9;
       alpha( 8)  =  v7;  beta( 8)  =  v9;    gamma( 8)  =  v8;    
       alpha( 9)  =  v8;  beta( 9)  =  v7;    gamma( 9)  =  v9;      
       alpha(10)  =  v8;  beta(10)  =  v9;    gamma(10)  =  v7;      
       alpha(11)  =  v9;  beta(11)  =  v8;    gamma(11)  =  v7;     
       alpha(12)  =  v9;  beta(12)  =  v7;    gamma(12)  =  v8;    

       v1  = 0.898905543365938d0;  v2 = 0.050547228317031d0

       alpha(13)  =  v1;  beta(13)  =  v2;    gamma(13)  =  v2;
       alpha(14)  =  v2;  beta(14)  =  v1;    gamma(14)  =  v2;       
       alpha(15)  =  v2;  beta(15)  =  v2;    gamma(15)  =  v1;   


       v1 = 1.0d0/3.0d0
       alpha(16)  =  v1;  beta(16)  =  v1;    gamma(16)  =  v1;   


    endif
    
    if(n.eq.25) then ! p = 10
       ww(1:3)   =  0.036725957756467d0
       ww(4:6)   =  0.045321059435528d0
       ww(7:12)  =  0.072757916845420d0   
       ww(13:18) =  0.028327242531057d0   
       ww(19:24) =  0.009421666963733d0   

       ww(25)    =  0.090817990382754d0

       v1  = 0.028844733232685d0;  v2 = 0.485577633383657d0
       v4  = 0.781036849029926d0;  v5 = 0.109481575485037d0
       v7  = 0.141707219414880d0;  v8 = 0.307939838764121d0;
       v9  = 1. - v7 - v8

       alpha( 1)  =  v1;  beta( 1)  =  v2;    gamma( 1)  =  v2;
       alpha( 2)  =  v2;  beta( 2)  =  v1;    gamma( 2)  =  v2;       
       alpha( 3)  =  v2;  beta( 3)  =  v2;    gamma( 3)  =  v1;     

       alpha( 4)  =  v4;  beta( 4)  =  v5;    gamma( 4)  =  v5;
       alpha( 5)  =  v5;  beta( 5)  =  v4;    gamma( 5)  =  v5;         
       alpha( 6)  =  v5;  beta( 6)  =  v5;    gamma( 6)  =  v4;        

       alpha( 7)  =  v7;  beta( 7)  =  v8;    gamma( 7)  =  v9;
       alpha( 8)  =  v7;  beta( 8)  =  v9;    gamma( 8)  =  v8;    
       alpha( 9)  =  v8;  beta( 9)  =  v7;    gamma( 9)  =  v9;      
       alpha(10)  =  v8;  beta(10)  =  v9;    gamma(10)  =  v7;      
       alpha(11)  =  v9;  beta(11)  =  v8;    gamma(11)  =  v7;     
       alpha(12)  =  v9;  beta(12)  =  v7;    gamma(12)  =  v8;    

       v7  = 0.025003534762686d0;  v8 = 0.246672560639903d0;
       v9  = 1. - v7 - v8
       alpha(13)  =  v7;  beta(13)  =  v8;    gamma(13)  =  v9;
       alpha(14)  =  v7;  beta(14)  =  v9;    gamma(14)  =  v8;    
       alpha(15)  =  v8;  beta(15)  =  v7;    gamma(15)  =  v9;      
       alpha(16)  =  v8;  beta(16)  =  v9;    gamma(16)  =  v7;      
       alpha(17)  =  v9;  beta(17)  =  v8;    gamma(17)  =  v7;     
       alpha(18)  =  v9;  beta(18)  =  v7;    gamma(18)  =  v8;    	 

       v7  = 0.009540815400299d0;  v8 = 0.066803251012200d0;
       v9  = 1. - v7 - v8
       alpha(19)  =  v7;  beta(19)  =  v8;    gamma(19)  =  v9;
       alpha(20)  =  v7;  beta(20)  =  v9;    gamma(20)  =  v8;    
       alpha(21)  =  v8;  beta(21)  =  v7;    gamma(21)  =  v9;      
       alpha(22)  =  v8;  beta(22)  =  v9;    gamma(22)  =  v7;      
       alpha(23)  =  v9;  beta(23)  =  v8;    gamma(23)  =  v7;     
       alpha(24)  =  v9;  beta(24)  =  v7;    gamma(24)  =  v8;    

       v1 = 1.0d0/3.0d0
       alpha(25)  =  v1;  beta(25)  =  v1;    gamma(25)  =  v1;   


    endif
    ! 0000000000000000000000000000000000000000000
    if(n.eq.42) then ! p =14 

       ww( 1: 3) =    0.021883581369429d0
       ww( 4: 6) =    0.032788353544125d0
       ww( 7: 9) =    0.051774104507292d0   !¿éòô¿¼âçóãn0 n1 n2à´¿øöæ
       ww(10:12) =    0.042162588736993d0
       ww(13:15) =    0.014433699669777d0   
       ww(16:18) =    0.004923403602400d0

       ww(19:24) =    0.024665753212564d0
       ww(25:30) =    0.038571510787061d0
       ww(31:36) =    0.014436308113534d0
       ww(37:42) =    0.005010228838501d0


       alpha( 1) = 0.022072179275643d0;  beta( 1) = (1.d0 - alpha( 1))/2.d0
       alpha( 4) = 0.164710561319092d0;  beta( 4) = (1.d0 - alpha( 4))/2.d0
       alpha( 7) = 0.453044943382323d0;  beta( 7) = (1.d0 - alpha( 7))/2.d0
       alpha(10) = 0.645588935174913d0;  beta(10) = (1.d0 - alpha(10))/2.d0
       alpha(13) = 0.876400233818255d0;  beta(13) = (1.d0 - alpha(13))/2.d0
       alpha(16) = 0.961218077502598d0;	beta(16) = (1.d0 - alpha(16))/2.d0


       alpha(19) = 0.057124757403648d0;  beta(19) = 0.172266687821356d0
       alpha(25) = 0.092916249356972d0;  beta(25) = 0.336861459796345d0
       alpha(31) = 0.014646950055654d0;  beta(31) = 0.298372882136258d0
       alpha(37) = 0.001268330932872d0;  beta(37) = 0.118974497696957d0



       gamma  = 1.d0 - alpha - beta


       alpha( 2)  =  beta( 1);         beta( 2)  = alpha( 1)           
       alpha( 3)  =  beta( 1);         beta( 3)  = gamma( 1)

       alpha( 5)  =  beta( 4);         beta( 5)  = alpha( 4)           
       alpha( 6)  =  beta( 4);         beta( 6)  = gamma( 4)           

       alpha( 8)  =  beta( 7);         beta( 8)  = alpha( 7)           
       alpha( 9)  =  beta( 7);         beta( 9)  = gamma( 7) 

       alpha(11)  =  beta(10);         beta(11)  = alpha(10)           
       alpha(12)  =  beta(10);         beta(12)  = gamma(10)   

       alpha(14)  =  beta(13);         beta(14)  = alpha(13)           
       alpha(15)  =  beta(13);         beta(15)  = gamma(13)   

       alpha(17)  =  beta(16);         beta(17)  = alpha(16)           
       alpha(18)  =  beta(16);         beta(18)  = gamma(16)   



       alpha(20)  = alpha(19);         beta(20)  = gamma(19)           
       alpha(21)  =  beta(19);         beta(21)  = alpha(19)           
       alpha(22)  =  beta(19);         beta(22)  = gamma(19)           
       alpha(23)  = gamma(19);         beta(23)  =  beta(19)           
       alpha(24)  = gamma(19);         beta(24)  = alpha(19)           

       alpha(26)  = alpha(25);         beta(26)  = gamma(25)           
       alpha(27)  =  beta(25);         beta(27)  = alpha(25)           
       alpha(28)  =  beta(25);         beta(28)  = gamma(25)           
       alpha(29)  = gamma(25);         beta(29)  =  beta(25)           
       alpha(30)  = gamma(25);         beta(30)  = alpha(25)   

       alpha(32)  = alpha(31);         beta(32)  = gamma(31)           
       alpha(33)  =  beta(31);         beta(33)  = alpha(31)           
       alpha(34)  =  beta(31);         beta(34)  = gamma(31)           
       alpha(35)  = gamma(31);         beta(35)  =  beta(31)           
       alpha(36)  = gamma(31);         beta(36)  = alpha(31)    

       alpha(38)  = alpha(37);         beta(38)  = gamma(37)           
       alpha(39)  =  beta(37);         beta(39)  = alpha(37)           
       alpha(40)  =  beta(37);         beta(40)  = gamma(37)           
       alpha(41)  = gamma(37);         beta(41)  =  beta(37)           
       alpha(42)  = gamma(37);         beta(42)  = alpha(37)    

       gamma  = 1.d0 - alpha - beta

    endif
    ! 0000000000000000000000000000000000000000000000000
    if(ww(1).eq.0.0d0)  then
       print*,"error in ng and please use the ng numbers in this subroutine"
       stop
    endif

  end subroutine ginttrg


! ************************************************************	
      subroutine gauleg(x1,x2,x,w,n)
! ************************************************************	
      implicit none
      integer n
      real(kind=dp):: x1,x2,x(n),w(n)
      real(kind=dp),parameter:: eps=3.d-14
      integer i,j,m
      real(kind=dp):: p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=dcos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z

          z=z1-p1/pp
        if(abs(z-z1).gt.eps)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      end subroutine gauleg
! ************************************************************	


 endmodule gaussintegraldata
