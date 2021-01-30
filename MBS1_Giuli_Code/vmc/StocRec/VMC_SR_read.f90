program VMC_read

  implicit none
  integer :: N,Nvar, i,j,k, ibin, BinDim, BinSkip, Nbin, Nused,INFO,LWORK
  real(8) :: Eavg,SEavg, Oavg(2), eOavg(2), empt, SdE(2)
  real(8) :: dEda(2), Delta, EpsReg,detS
  real(8), dimension(:), allocatable :: e,Avar, WORK
  real(8), dimension(:,:),allocatable:: O, S,checkS
  real(8), dimension(:), allocatable  :: Etot
  integer, dimension(:), allocatable :: x, IPIV

  Nvar=2!Number of variables
  LWORK=4*Nvar
  allocate(Avar(Nvar))
  EpsReg=0.001d0
  
!  write(*,*) 'Enter BinDim, BinSkip'
  read(*,*) BinDim, BinSkip, Delta
  
  open(12,file='data.dat',status='old',action='read')
  read(12,*) N, Avar(1), Avar(2)
  Nbin = (N/BinDim -BinSkip)
  Nused=Nbin-BinSkip !Num of Bin used, some discarded for equilib
  
  allocate(e(BinDim),x(BinDim),O(2,BinDim),S(2,2),checkS(2,2))
  allocate(Etot(BinSkip+1:Nbin))
  allocate(IPIV(Nvar),WORK(LWORK))
  S=0.0d0
  Oavg=0.0d0
  eOavg=0.0d0
  
  do ibin=1,Nbin!cycle over bins
     
     do i=1,Bindim
        read(12,*) e(i),x(i), O(1,i), O(2,i)
     enddo

     if(ibin.le.BinSkip) cycle
     
     Etot(ibin)=sum(e)/BinDim

     do j=1,Nvar
        Oavg(j)=Oavg(j) + sum(O(j,:))/BinDim
        eOavg(j)=eOavg(j) + sum(O(j,:)*e)/BinDim
     enddo
     
     !For Stoc. Rec.
     do j=1,Nvar
        do k=j,Nvar
           S(j,k) = S(j,k) + sum(O(j,:)*O(k,:))/BinDim
           S(k,j) = S(j,k)
        enddo
     enddo
          
  enddo
  close(12)
  
  Eavg = sum(Etot)/dble(Nused)
  SEavg= dsqrt( sum(Etot**2)/dble(Nused) - Eavg**2)/sqrt(dble(Nused))
  do i=1,Nvar
     Oavg(i) =  Oavg(i)/dble(Nused)
     eOavg(i)= eOavg(i)/dble(Nused)
     dEda(i) = 2.0d0*( eOavg(i)-Eavg*Oavg(i)  )
  enddo

  !For Stoc Rec
  do j=1,Nvar
     do k=j,Nvar
        S(j,k) = S(j,k)/dble(Nused) -Oavg(j)*Oavg(k)
        S(k,j) = S(j,k)
     enddo
  enddo
  
!  print*,'S'
!  print*, '|',S(1,1),'  ',S(1,2),'|'
!  print*, '|',S(2,1),'  ',S(2,2),'|'
  
  do j=1,Nvar
     S(j,j)=S(j,j)+EpsReg
  enddo

  !Inversion for 2x2 S:
  detS=S(1,1)*S(2,2)-S(2,1)*S(1,2)
  SdE(1)= ( S(2,2)*dEda(1)-S(1,2)*dEda(2) )/detS
  SdE(2)= (-S(2,1)*dEda(1)+S(1,1)*dEda(2) )/detS

  
  
!Generical inversion for S
!  call dgetrf(Nvar,Nvar,S,Nvar,IPIV,INFO)
!  if(INFO>0) stop 'error during dgetrf'
!  call dgetri(Nvar,S,Nvar,IPIV,WORK,LWORK,INFO)
!  if(INFO>0) stop 'error during inversion, info+'
!  if(INFO<0) stop 'error during inversion, info-'
! SdE=matmul(S,dEda)
  
!  checkS=matmul(checkS,S)
!  print*, '|',checkS(1,1),'  ',checkS(1,2),'|'
!  print*, '|',checkS(2,1),'  ',checkS(2,2),'|'

  
  
  

  Avar = Avar-Delta*SdE
  open(14,file='variable.dat',status='new')
  do i=1,Nvar
     write(14,*) Avar(i)
  enddo
  close(14)
  
  write(*,*) 'variables:', Avar
  write(*,*) 'Energy:', Eavg, ' +- ', SEavg
  
  open(13,file='graph.dat',status='old',access='append')
  write(13,*) Avar(1), Avar(2), Eavg, SEavg, Oavg(1),dEda(1),dEda(2),SdE(1),SdE(2)
  close(13)
  
end program VMC_read
