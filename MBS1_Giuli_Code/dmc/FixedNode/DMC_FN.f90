MODULE modello
  real(8) :: a, LAMB
  real(8) :: t2,V
  integer :: L, Nstep
  
  public :: wfcG, init, L, LAMB, GreenF, Nstep, WFCandH
contains
  
  subroutine init(x)
    integer, intent(out) :: x
    write(*,*) 'Enter a,t2/t,V/t,L,X_0,Nstep'
    read(*,*) a,t2,V,L,x,Nstep
    LAMB=V*L
    
    if( (x<1).or.(x>L) ) stop 'invalid x_0'
  end subroutine init
  
  
  real(8)  function wfcG(x)
    integer, intent(in) :: x
    if( (x.ge.1).and.(x.le.L) ) then
       wfcG = dexp(-a*x)
    else
       wfcG = 0.0d0
    endif
  end function wfcG

  subroutine WFCandH(psi,H,S,x,eLoc)
    real(8), intent(out) :: psi(-2:2), H(-2:2), S(-2:2), eLoc
    integer, intent(in) :: x
    integer ::i
    
    do i=-2,2
       psi(i)=wfcG(x+i)
    enddo
    
    H(0) = 0
    H(1) = -1.0d0
    H(-1)= -1.0d0
    H(2) = t2
    H(-2)= t2
    
    S=psi*H/psi(0)
    H(0)=+V*x
    eLoc = sum(psi*H)/psi(0)
    where(S>0) H = 0.0d0
    H(0)=H(0)+0.5d0*(sum(S)+sum(abs(S)))!-eLog per eff

    
    
  end subroutine WFCandH

  subroutine GreenF(x,b,psi,H)
    integer :: k
    integer, intent(out) :: x
    real(8), intent(out) :: psi(-2:2),H(-2:2),b
    real(8) :: z,zsup,G(-2:2)

    G=-H*psi/psi(0)
    G(0)=LAMB+G(0)
    !supposing for the moment G>0
    b=sum(G)
    call random_number(z)
    z=z*b
    zsup=G(-2)
    do k=-2,2
       if(z.le.zsup) then
          x=x+k
          exit
       endif
       zsup=zsup+G(k+1)
    enddo
    
    
  end subroutine GreenF

END MODULE modello

PROGRAM DMC

  use modello
  implicit none
  integer :: x=10, xtrial, N, ifreq
  integer :: i,j,k
  real(8) :: H(-2:2),psi(-2:2),S(-2:2), b,  eLoc

  call init(x)
  call WFCandH(psi,H,S,x,eLoc)
  
  open(14,file='Process.dat',status='new')
  
  do i=1,Nstep
     call GreenF(x,b,psi,H)
     call WFCandH(psi,H,S,x,eloc)

     write(14,*) b/LAMB, eLoc,  x 
  enddo
  close(14)
END PROGRAM DMC
