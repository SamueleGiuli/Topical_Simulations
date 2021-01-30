module modello

  implicit none
  
  real(8) :: a,b !variational parameters
  real(8) :: t2,V !model parameters
  integer :: L
  
  public :: wfc, init, L
contains
  
  subroutine init()
!    write(*,*) 'Enter a,b,t2/t,V/t,L'
    read(*,*) a,b,t2,V,L
  end subroutine init
  
real(8)  function wfc(x)
    integer, intent(in) :: x
    if( (x.ge.1).and.(x.le.L) ) then
       wfc = dexp(-a*x)*x**b
    else
       wfc = 0.0d0
    endif
  end function wfc

end module modello

program vmc
  
  use modello
  implicit none
  integer :: x=10, xtrial, N, ifreq
  integer :: i,j,k
  real(8) :: r, z(2), e, psi, psiTrial, psiL, psiR, psiLL, psiRR
  real(8) :: O(2)
  call init()

  !write(*,*) 'Enter N, ifreq'
  read(*,*)  N, ifreq

  psi = wfc(x)
  psiR=wfc(x+1)
  psiL=wfc(x-1)
  psiRR=wfc(x+2)
  psiLL=wfc(x-2)

  open(12,file='data.dat',status='new',action='write')
  write(12,*) N/ifreq, a ,b
  
  do i=1,N

     call random_number(z)
     
     if( (x.ge.1).and.(x.le.L) ) then
        
        if(z(1)<0.5) then
           psiTrial = psiL
           xTrial = x-1
        else
           psiTrial = psiR
           xTrial = x+1
        endif
        
        r=(psiTrial/psi)**2
        
        if(z(2)<r) then
           x=xTrial
           psi=psiTrial
           psiR=wfc(x+1)
           psiL=wfc(x-1)
        endif
     else
        stop 'Error in walking'
     endif
     
     if (mod(i,ifreq).eq.0) then
        !Stoc.Ren.
        O(1) = -x
        O(2) = log(dble(x))
        psiRR=wfc(x+2)
        psiLL=wfc(x-2)
        e = V*x + ( -(psiR+psiL) +t2*(psiRR+psiLL) )/psi
        write(12,*) e,x,O(1),O(2)
     endif
     
  enddo

  close(12)
     
     

  
end program vmc
