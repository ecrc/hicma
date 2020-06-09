!*************************************************************************
subroutine fan_gmres(nd)
  use global_com,only: dp,nitermax,precis
  use global_dim
  implicit none
  integer,parameter :: m = 30
  integer::itertotal,iterout,n,nd,k
  complex(kind=dp)::cx(nd)
  complex(kind=dp)::cb(m+1)
  complex(kind=dp)::cy(nd,m+1),ch(m+1,m) 
  complex(kind=dp)::cs(m)
  real(kind=dp)::rc(m)
  complex(kind=dp)::ctemp
  real(kind=dp)::ay0, be, bea
  real(kind=dp)::xtime,xtim
  real(kind=dp),external:: scnrm22

  print *, '.......in fan_gmres algorithm ......'
  open(unit=161,file='con_fangmres.out')
  PRINT*, "Starting iteration of solution..."
  call cpu_time(xtime)

  cx    = (0.0_dp,0.0_dp)
  ay0   = scnrm22(nd,crhs)
  itertotal = 0

!print*,"nd",nd

!************************************************************************
  do iterout = 1, nitermax  
!************************************************************************
     itertotal = itertotal +1
!====================================================================
     call            convol(cx(1:nd),cy(1:nd,m+1),nd)
     cy(1:nd,m+1)  = crhs(1:nd) -cy(1:nd,m+1)
     cb(1)         = scnrm22(nd,cy(1:nd,m+1))  
     be            = cb(1)/ay0
!====================================================================
     !write(*,*) 'restart', iterout,itertotal,'  truebe',be
     !write(161,*) itertotal,be
     if(be.le.precis.or.itertotal>nitermax) exit
     cy(1:nd,1) =  cy(1:nd,m+1)/cb(1)
!************************************************************************
     do n = 1, m
!************************************************************************
        itertotal = itertotal +1    
!========================step1=======================================
        call convol(cy(1:nd,n),cy(1:nd,n+1),nd) 
!========================step2=======================================
  	do k=1,n
           ch(k,n) = dot_product(cy(1:nd,k),cy(1:nd,n+1))
           cy(1:nd,n+1) = cy(1:nd,n+1) - cy(1:nd,k)*ch(k,n)  
        enddo
        ch(n+1,n)  =  scnrm22(nd,cy(1:nd,n+1)) 
        if(n.lt.m)  cy(1:nd,n+1) = cy(1:nd,n+1)/ch(n+1,n)
        if(n.ne.1) then
           do k=1,n-1
              ctemp       = rc(k)*ch(k+1,n) -        cs(k) *ch(k,n)
              ch(k,n)     = rc(k)*ch(k,n)   + conjg (cs(k))*ch(k+1,n)
              ch(k+1,n)   = ctemp 
           enddo
        endif
!========================step3=======================================
        call givens_loc(ch(n,n),ch(n+1,n),rc(n),cs(n))
!========================step4=======================================
        cb(n+1) =  -cs(n)*cb(n)
        cb(n)   =   rc(n)*cb(n)   
        ch(n,n)   = rc(n)*ch(n,n)+ conjg (cs(n))*ch(n+1,n)  
        ch(n+1,n) = (0.0_dp,0.0_dp)
!========================step5=======================================
        bea = abs(cb(n+1))/ay0
        !write(*,*) itertotal,bea
        !write(161,*) itertotal,bea
        if(n.eq.m.or.bea.le.precis) then
           call cutrisub(n, ch(1:n,1:n), cb(1:n), cb(1:n))
           cy(1:nd,m+1) = matmul(cy(1:nd,1:n),cb(1:n))
           cx = cx + cy(1:nd,m+1)	     
           exit
        endif
!===================================================================
     enddo   !n }
  enddo !iterout}
!===================================================================
  if(be>precis) then
     write( *,*) 'convergence is not achieved'
     write(161,*) 'convergence is not achieved'
  else
     write( *,*) 'convergence achieved'
     write(161,*) 'convergence achieved'
  endif
  call cpu_time(xtim)
  write(*, *) "gmres itr time...", xtim-xtime
  write(161,*) "gmres itr time...", xtim-xtime 

!------------------------solution------------------------
  rj = cx

  close(161)

  return
end subroutine fan_gmres
!*************************************************************************

!*************************************************************************
subroutine givens_loc(z1,z2,c,s)
  use global_com,only: dp
  implicit none
  complex(kind=dp)::  z1,z2,s
  real(kind=dp)::      c,vnormz

  vnormz = sqrt ( abs(z1)**2 + abs(z2)**2)
  if (abs(z1).ne.0.) then
     c =     abs(z1)/vnormz
     s =     z2/z1* c
  else if (abs(z2).ne.0.) then
     c = 0.
     s = z2/ abs(z2)
  else
     c = 1.
     s = 0.
  endif

end subroutine givens_loc
!*************************************************************************

!*************************************************************************
function  scnrm22(nd,cy9) 
  use global_com,only: dp
  implicit none
  integer nd,i
  complex(kind=dp)::  cy9(nd)
  real(kind=dp):: scnrm22

  scnrm22=0.0_dp
  do i=1,nd
     scnrm22= scnrm22+( conjg( cy9(i) )*cy9(i))
  end do
  scnrm22=sqrt(scnrm22)
  !	write(*,*) scnrm22
  !	pause

  return
end function scnrm22
!*************************************************************************

!*************************************************************************
subroutine cutrisub(n, cma, cvx, cvy)
  use global_com,only: dp
  implicit none
  integer n
  complex(kind=dp):: cma(n, n), cvx(n), cvy(n)
  integer i, j

  do i = n, 1, - 1
     cvy(i) = cvx(i)
     do j = i + 1, n
        cvy(i) = cvy(i) - cma(i, j) * cvy(j)
     enddo
     cvy(i) = cvy(i) / cma(i, i)
  enddo

end subroutine cutrisub
!*************************************************************************

!*************************************************************************
subroutine convol(cin,cout,nd)  !!!convolnear
  use global_com,only: dp,ntriangle,nipp
  use global_dim,only: zzsparse
  implicit none
  integer::nd,i,kk1,kk2
  complex(kind=dp):: cin(nd),cout(nd)

  cout  = (0.0_dp,0.0_dp)

  do i = 1, ntriangle*nipp
     kk1 = (i-1)*ntriangle*nipp+1
     kk2 = i*ntriangle*nipp
     cout(i) = dot_product(conjg(zzsparse(kk1:kk2)),cin(1:ntriangle*nipp))
  enddo

end subroutine convol
!*************************************************************************
