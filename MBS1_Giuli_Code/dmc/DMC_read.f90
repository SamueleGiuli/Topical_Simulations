

program DMC_read
  implicit none

  integer :: index,mindex,j,k,m
  integer,dimension(:),allocatable :: x
  integer :: Nstep,p,Nbin,BinDim,BinSkip,Nused
  real(8) :: Eavg, SigEavg, Xavg, SigXavg, g, eLoc
  real(8),dimension(:),allocatable :: E2bin,EGbin, Gbin, Xbin, b

!  write(*,*) ' Nstep,BinDim,BinSkip,p'
  read(*,*) Nstep,BinDim,BinSkip,p
  m=p/2
  Nbin=Nstep/BinDim
  Nused=Nbin-BinSkip
  
  allocate(EGbin(BinSkip+1:NBin),E2bin(Binskip+1:Nbin),Gbin(BinSkip+1:NBin),Xbin(BinSkip+1:NBin))
  allocate(b(1:p+1))
  allocate(x(1:m))

  open(14,file='Process.dat',status='old')
!  open(15,file='ResultDMC.dat',status='old')
  
  b=0.0d0
  eLoc=0.0d0
  EGbin=0.0d0
  E2bin=0.0d0
  Gbin=0.0d0
  Xbin=0.0d0
  index=1
  mindex=1
  
  do j=1,Nbin
     do k=1,BinDim
        read(14,*) b(index), eLoc, x(mindex) 
        index = mod(index,p+1)+1
        mindex= mod(mindex,m)+1
        if(j.le.BinSkip) cycle
        g=product(b)/(b(index))
        Xbin(j) = Xbin(j) + g*dble(x(mindex))
        EGbin(j) = EGbin(j) + g*eLoc
        Gbin(j) = Gbin(j) + g
!        print*,EGbin(j),Gbin(j), j, k,g
     end do         
  end do
  close(14)

  
  Eavg=sum(EGbin)/sum(Gbin) 
  SigEavg=dsqrt( sum(Gbin*(EGbin/Gbin - Eavg)**2 ) )/dsqrt(sum(Gbin)*dble(Nused-1))
  
  Xavg=sum(Xbin)/sum(Gbin) 
  SigXavg=dsqrt( sum(Gbin*(Xbin/Gbin - Xavg)**2 ) )/dsqrt(sum(Gbin)*dble(Nused-1))

  
  write(*,'(I8,4E15.8)') p, Eavg, SigEavg, Xavg, SigXavg
  
end program DMC_read
