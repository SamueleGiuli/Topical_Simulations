!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Fermi-Pasta-Ulam-Tsingou Problem
! Samuele Giuli
! h(i) = 0.5*v(i)**2  + 0.5*(u(i)-u(i+-1))^2 + 0.25*alpha*(u(i)-u(i+-1))^4
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module subr
  implicit none

  integer :: N, m, nstep, nstepPS, Nomega=200
  real(8), dimension(:), allocatable :: u, v, f
  real(8), dimension(:), allocatable :: q, Qk, Pk, Ek, Ekfin, OmegaK
  complex(8), dimension(:,:), allocatable :: C
  real(8) :: A, alpha, dt, T, ekin, epot, etot, Ektot, EpN, omeg0
  real(8), parameter :: pi=3.14159265358979d0
  logical :: displacement, pos
  
  public :: u, v,  f, N, nstep, dt, C, OmegaK, nstepps
  contains
  
  subroutine init()
    integer :: i,j
    
    print*, 'Enter the number of sites (N)>'
    read*, N
    allocate(u(0:N+1),v(1:N),f(1:N),q(1:N))
    allocate(Qk(1:N),Pk(1:N),Ek(1:N),Ekfin(1:N),OmegaK(1:N))
    allocate(C(1:Nomega,1:N))


    C=0.d0
    Ekfin=0.d0
    Ektot=0.d0
    
    !compute modes' frequencies
    do j=1,N
       q(j) = j*pi/dble(N+1)
       OmegaK(j) = 2.d0*sin(q(j)/2.d0) 
    enddo
    
    print*, 'Enter the anharmonic coefficient (alpha)>'
    read*, alpha
    print*, 'Enter the time step (dt) and the simulation time (T)>'
    read*, dt, T
    nstep = int(T/dt)  !num of simulation steps that must be even
    if(mod(nstep,2)==1) nstep=nstep+1
    
    nstepPS= nstep/2 !num os PowerSpectrum integration steps
    omeg0=2*pi/dble(nstepPS)
    print*, 'For the intial dispacemente enter .t. for a normal mode of the linear system  and .f. for a random distribution>'
    read*, displacement
    if(displacement) then
       print*, 'Enter the mode (m) and Energy per site (EpN)'
       read*, m, EpN

       if( (m>N).or.(m<1) ) stop 'm must be in [1:N]'

       A=sqrt(EpN*N/dble(N+1))/sin(pi*m/dble(2*N+2))
       !Initial configuration of m-th mode
       do j=1,N
          
          v(j) = 0.d0
          u(j) = A*sin(q(m)*dble(j))
       enddo
    else
       print*,'Use the previous configuration?'
       read*, pos
       if(pos) then
          open(40,file='previous.dat',status='old',action='read')
          do j=1,N
             read(40,*) u(j)
          enddo
          close(40)
       else
          print*, 'Enter the amplitude'
          read*, A       
          !Random configuration
          call random_number(u)
          u = (u-0.5d0)*2.0d0
          open(40,file='previous.dat')
          do j=1,N
             write(40,*) u(j)
          enddo
       endif
       v = 0.d0
    endif
    !
    !Fixed Boundary Conditions
    u(0) = 0.d0
    u(N+1)=0.d0
    
    call force()
    
  end subroutine init


  !This subroutines evolves the u,v parameters
  !using the Velocity Verlet Algorithm
  subroutine move()
    integer :: i

    do i=1,N
       u(i) = u(i) + v(i)*dt +0.5d0*f(i)*dt**2
       v(i) = v(i) +0.5d0*f(i)*dt
    enddo
    call force()
    do i=1,N
       v(i) = v(i) +0.5d0*f(i)*dt
    enddo
  end subroutine move

  
  !This subroutine compute the force
  !at the present configuration
  subroutine force()
    integer :: j
    do j=1,N
       f(j) = (u(j+1)+u(j-1)-2.0d0*u(j)) +alpha*( (u(j+1)-u(j))**3.d0 +(u(j-1)-u(j))**3.d0 )
       ! this is the force acting on j-th oscillator      
    enddo
  end subroutine force

  !This subroutine compute the energy
  !at the present configuration
  subroutine energy()
    integer :: i
    
    etot = 0.d0
    ekin = 0.d0
    epot = 0.5d0*(u(0)-u(1))**2 + 0.25d0*alpha*(u(0)-u(1))**4

    do i=1,N
       ekin = ekin + 0.5d0*v(i)**2
       epot = epot + 0.5d0*(u(i)-u(i+1))**2 + 0.25d0*alpha*(u(i)-u(i+1))**4
    enddo
    etot = ekin + epot
  end subroutine energy

  !This subroutines compute the harmonic modes'
  !contributions at the present configuration
  subroutine modes(step)
    integer :: i,j
    integer, intent(in) :: step

    !Power Spectrum r the last steps
    if(step > nstepPS) then
       do i=1, Nomega
          C(i,:)=C(i,:)+u(1:N)*exp(omeg0*(step-nstepPS)*i*(0.d0,1.d0))
       enddo
    endif

    !Normal Modes
    Qk=0.d0
    Pk=0.d0
    Ek=0.d0
    do i=1,N
       do j=1,N
          Qk(i)=Qk(i)+u(j)*sin(pi*i*dble(j)/dble(N+1))
          Pk(i)=Pk(i)+v(j)*sin(pi*i*dble(j)/dble(N+1))
       enddo
       Qk(i) = Qk(i)*sqrt(2.d0/dble(N+1))
       Pk(i) = Pk(i)*sqrt(2.d0/dble(N+1))
       Ek(i) = 0.5d0*(Pk(i)**2 + (OmegaK(i)*Qk(i))**2 )
    enddo
    Ektot = sum(Ek)
    Ek(:) = Ek(:)/sum(Ek)
    Ekfin=Ekfin+Ek
  end subroutine modes

  subroutine mean_var(Ek,dim,sigma,mean)
    integer, intent(in) :: dim
    real(8), intent(in) :: Ek(1:dim)
    real(8), intent(inout):: sigma, mean
    integer :: i
    real(8) :: var
    mean = sum(Ek)/dim
    var=0.d0
    do i=1,dim
       var=var+(Ek(i)-mean)**2
    enddo
    var=var/dble(N-1)
    sigma=sqrt(var)    
  end subroutine mean_var
end module subr


program FPUT

  use subr
  implicit none
  integer :: i,j
  real(8) :: sigma, mean

  
  call init()
  !open(12,file='position.dat')     !positions animation
  open(13,file='energy.dat')        !Etot
  !open(14,file='anim_modes.dat')   !modes animation
  open(15,file='modes_t.dat')       !E_k (t)
  open(16,file='power_spectrum.dat')!Energy power spectrum
  open(17,file='avg_modes.dat')
  open(18,file='freq.dat')
  
  write(13,*) '# TimeStep, ETOT, EKIN, EPOT'
  
  !do j=0,N+1
  !   write(12,*) j, u(j)
  !enddo
  !write(12,*) ''
  !write(12,*) ''

  do i=1, nstep
     call move()

     !write positions
!     do j=0,N+1
!        write(12,*) j, u(j)
!     enddo
!     write(12,*) ''
!     write(12,*) ''

     !compute modes and partial energies
     call energy()
     call modes(i)
     if( mod(i,10)==1 ) then

        !CODE FOR ANIMATION
!        do j=1,N
!           write(14,*) j, Ek(j)
!        enddo  
!        write(14,*) ''
!        write(14,*) ''

        
        !Write energy in two different ways
        write(13,*) i*dt, Ektot, etot

        !Write the first 10 modes
        write(15,fmt='(11E16.8)') i*dt, Ek(1:10)
        write(17,fmt='(11E16.8)') i*dt, Ekfin(1:10)/dble(i)
        
        !call mean_var(Ek,N,sigma,mean)
        !if( (sigma/mean<0.05) ) then
        !   write(17,*) m, alpha, EpN, i*dt
        !   print*, 'm, a, epn, i*dt'
        !   print*, m, alpha, EpN, i*dt
        !   stop
        !endif
                 
     endif     
  enddo

  !Write the modes energy averaged in time
  do i=1,Nomega
     write(16,fmt='(2E15.8)') omeg0*i/dt , real(dot_product(C(i,1:n),C(i,1:n)))/dble(N)/(nstepPS**2)
  enddo
  print*, dot_product(C(1,1:n),C(1,1:n))
  do i=1,N
     write(18, '(E15.8,I3)') OmegaK(i), 1
  enddo
  
  
  
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  
end program FPUT
