

program DMC_read
  implicit none

  integer :: index,indexprev,mindex,j,k,m
  integer,dimension(:),allocatable :: x
  integer :: Nstep,p,Nbin,BinDim,BinSkip,Nused
  real(8) :: Eavg, SigEavg, Xavg, SigXavg, Savg, SigSavg, g, eLoc, ES, sES, gsign
  real(8),dimension(:),allocatable :: E2bin,EGbin, Gbin, Xbin, b, Sbin

!  write(*,*) ' Nstep,BinDim,BinSkip,p'
  read(*,*) Nstep,BinDim,BinSkip,p
  m=p/2
  Nbin=Nstep/BinDim
  Nused=Nbin-BinSkip
  
  allocate(EGbin(BinSkip+1:NBin),E2bin(Binskip+1:Nbin),Gbin(BinSkip+1:NBin),Xbin(BinSkip+1:NBin),Sbin(BinSkip+1:NBin))
  allocate(b(1:p+1))
  allocate(x(1:m))

  open(14,file='Process.dat',status='old')
  
  b=0.0d0
  eLoc=0.0d0
  EGbin=0.0d0
  E2bin=0.0d0
  Gbin=0.0d0
  Xbin=0.0d0
  Sbin=0.0d0
  index=1
  mindex=1
  
  do j=1,Nbin
     do k=1,BinDim
        read(14,*) b(index), eLoc, x(mindex)
        mindex= mod(mindex,m)+1
        indexprev=index
        index = mod(index,p+1)+1
        
        if(j.le.BinSkip) cycle
        
        g=(product(b)/(b(indexprev)))
        Gbin(j) = Gbin(j) + abs(g)
        
        Xbin(j) = Xbin(j) + g*dble(x(mindex))
        EGbin(j) = EGbin(j) + g*eLoc

        Sbin(j) = Sbin(j) + g
        
     end do         
  end do
  close(14)

  

  
  Eavg=sum(EGbin)/sum(Gbin) 
  SigEavg=dsqrt( sum(Gbin*(EGbin/Gbin - Eavg)**2 ) )/dsqrt(sum(Gbin)*dble(Nused-1))
  !SigEavg=dsqrt( sum((EGbin/Gbin-Eavg)**2))/dsqrt(dble(Nused-1))
  
  Xavg=sum(Xbin)/sum(Gbin) 
  SigXavg=dsqrt( sum(Gbin*(Xbin/Gbin - Xavg)**2 ) )/dsqrt(sum(Gbin)*dble(Nused-1))
  !SigXavg=dsqrt( sum(Xbin**2/Gbin)/sum(Gbin)-sum(Xbin/Gbin)**2  )/dsqrt(dble(Nused-1))
  !SigXavg=dsqrt(sum((Xbin/Gbin-Xavg)**2))/dsqrt(dble(Nused-1))
  
  Savg=sum(Sbin)/sum(Gbin) 
  SigSavg=dsqrt( sum(Gbin*(Sbin/Gbin - Savg)**2 ) )/dsqrt(sum(Gbin)*dble(Nused-1))
  !SigEavg=dsqrt( sum(Sbin**2/Gbin)/sum(Gbin)-sum(Sbin/Gbin)**2  )/dsqrt(dble(Nused-1))
  !SigSavg=dsqrt(sum((Sbin/Gbin-Savg)**2))/dsqrt(dble(Nused-1))
  
  ES=Eavg/Savg
  sES=ES*dsqrt( (SigEavg/Eavg)**2 + (SigSavg/Savg)**2  )

  open(15,file='p_sim.dat',status='new')
  write(15,fmt='(1I6,8E15.8)') p, ES, sES,Eavg,SigEavg, Xavg, SigXavg, Savg, SigSavg
  close(15)

end program DMC_read

