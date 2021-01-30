program VMC_read

  implicit none
  integer :: N, i, ibin, BinDim, BinSkip, Nbin, Nused
  real(8) :: Eavg,SEavg, Oavg(2), eOavg(2), a, b, empt
  real(8) :: dEda(2),Delta,epsil
  real(8), dimension(:), allocatable :: e, Etot
  real(8), dimension(:,:),allocatable:: O, Otot, eOtot
  integer, dimension(:), allocatable :: x

!  write(*,*) 'Enter BinDim, BinSkip'
  read(*,*) BinDim, BinSkip, Delta
  
  open(12,file='data.dat',status='old',action='read')
  read(12,*) N, a, b
!  print*, 'N:', N
  Nbin = (N/BinDim -BinSkip)
  Nused=Nbin-BinSkip
!  write(*,*) 'Nbin: ',Nbin
  
  allocate(e(BinDim),x(BinDim),O(2,BinDim))
  allocate(Etot(BinSkip+1:Nbin),Otot(2,BinSkip+1:Nbin),eOtot(2,BinSkip+1:Nbin))

  do ibin=1,Nbin
     
     do i=1,Bindim
        read(12,*) e(i),x(i), O(1,i), O(2,i)
     enddo

     if(ibin.le.BinSkip) cycle
     
     Etot(ibin)=sum(e)/dble(BinDim)
     Otot(1,ibin)=sum(O(1,:))/dble(BinDim)
     Otot(2,ibin)=sum(O(2,:))/dble(BinDim)
     eOtot(1,ibin)=sum(O(1,:)*e)/dble(BinDim)
     eOtot(2,ibin)=sum(O(2,:)*e)/dble(BinDim)
     
     
     
        
  enddo
  close(12)

  Eavg = sum(Etot)/dble(Nused)
  SEavg= sqrt( sum(Etot**2)/dble(Nused) - Eavg**2)/sqrt(dble(Nused-1))
!  write(*,*) 'Eavg: ', Eavg, '+-', SEavg

  do i=1,2
     Oavg(i) = sum(Otot(i,:))/dble(Nused)
     eOavg(i)= sum(eOtot(i,:))/dble(Nused)
     dEda(i) =2.0d0*( eOavg(i)-Eavg*Oavg(i)  )
  enddo

  
  write(*,*) 'a:', a, 'b:', b
  write(*,*) 'Energy:', Eavg, ' +- ', SEavg
  
  open(13,file='graph.dat',status='old',access='append')
  write(13,*) a, b, Eavg, SEavg
  close(13)

  !Updating a and b
  a= a-Delta*dEda(1)
  !b= b-Delta*dEda(2)

  open(14,file='variable.dat',status='new')
  write(14,*) a, b
  close(14)
  
  
end program VMC_read
