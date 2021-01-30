module sub_mod

  implicit none

  public :: S_i, S_z, Get_Norms, Get_Ham_Diag, Get_Ham_NoDiag, Lanczos
  contains


    
!********** S_z at position "i" state "N" **********
  integer function S_i(N,i)
    integer(16), intent(in) :: N
    integer, intent(in) :: i
    integer :: j,jp

    if( btest(N,i) ) then
       S_i = 1
    else
       S_i = -1
    end if
    
  end function S_i
    
!********** Total S_z of state "N" **********
  integer function S_z(N,L)
    integer(16), intent(in) :: N
    integer, intent(in) :: L
    integer :: i, Stot, s_trial
    Stot=0
    do i=0,L-1
       
       Stot = Stot + S_i(N,i)
    enddo
    S_z = Stot
  end function S_z

  !********** Set S_i to Val in status_swapped **********
  

  subroutine SetS_i(status_swapped,i,Val)
    integer(16), intent(out) :: status_swapped
    integer,intent(in) :: i,Val
    integer :: j,jp
    j=2*i
    jp=j+1
    if (Val==0)then
       status_swapped=ibclr(status_swapped,j)
       status_swapped=ibclr(status_swapped,jp)
    else if(Val==1)then
       status_swapped=ibset(status_swapped,j)
       status_swapped=ibclr(status_swapped,jp)
    else if(Val==-1)then
       status_swapped=ibclr(status_swapped,j)
       status_swapped=ibset(status_swapped,jp)
    else
       stop 'Error in Setting S_i'
    end if
    
  end subroutine SetS_i   

  
  !********** Computing Norms at fixed Sz Qx **********
  subroutine Get_Norms(NormsSzQx,SzSpace,N_Sz,Qx,L)
    integer, intent(in) :: L
    integer(16),intent(in) :: N_Sz
    integer(16), dimension(1:N_Sz), intent(in) :: SzSpace
    real(8), dimension(1:N_Sz), intent(out) :: NormsSzQx
    real(8), intent(in) :: Qx
    integer :: i,j,j_period,k
    integer(16) :: status_i, status_shift
    complex(8), dimension(:), allocatable :: status
    complex(8), parameter :: img = dcmplx(0.0d0,1.0d0)
    
    do i = 1,N_Sz
       status_i = SzSpace(i)
       status_shift = status_i

       !get periodicity
       do j = 1,L
          status_shift = mod(int(2,16)*status_shift,int(2,16)**L-1)
          if( status_shift == status_i) then
             j_period = j
             exit
          endif
       enddo

       !get norm
       if(j_period == L) then
          NormsSzQx(i) = sqrt(dble(L))
       else
          allocate(status(0:j_period))
          status = dcmplx(0.0d0,0.0d0)
          do j=1,L
             k=mod(j,j_period)
             status(k) = status(k)+exp( Qx*j*img )
          enddo
          NormsSzQx(i) = dsqrt(sum(abs(status)**2)  )
          deallocate(status)
       endif
    end do
  end subroutine Get_Norms

!********** Applying the S+ S- operator at site "i" to "status_j" **********
  
  subroutine SpSm(status_j,status_swapped,i,info,L)
    integer(16), intent(in) :: status_j
    integer(16), intent(out) :: status_swapped
    integer, intent(in) :: i,L
    logical, intent(out) :: info
    integer :: Si, Sip, ip
    ip=mod(i+1,L)
    Si=S_i(status_j,i)
    Sip=S_i(status_j,ip)
    info=.false.
    status_swapped=status_j
    if(Si==-1 .and. Sip==+1)then
       status_swapped=ibset(status_swapped,i)
       status_swapped=ibclr(status_swapped,ip)
    else
       info=.true.
    end if   
  
  end subroutine SpSm


!********** Building the non-diagonal elements of Hamiltonian **********  
  subroutine Get_Ham_NoDiag(H_NoDiag,H_NoDiag_index,SpaceSzQx,NormsSzQx,N_SzQx,L,Qx)
    integer, intent(in) :: L
    integer(16),intent(in) :: N_SzQx
    real(8), intent(in) :: Qx
    integer(16), dimension(1:N_SzQx), intent(in) :: SpaceSzQx
    real(8), dimension(1:N_SzQx), intent(in) :: NormsSzQx
    integer(8), dimension(1:N_SzQx,1:L), intent(out) :: H_NoDiag_index
    complex(8), dimension(1:N_SzQx,1:L), intent(out) :: H_NoDiag
    integer :: i,j,k,j_x, low, high
    integer(16) :: status_j, status_swap, status_rep
    logical :: info, goal_not_found
    
    H_NoDiag_index=0
    H_NoDiag = cmplx(0.0d0,0.0d0)
    
    do j=1,N_SzQx
       status_j = SpaceSzQx(j)
       do i=0,L-1
          call SpSm(status_j,status_swap,i,info,L)
          if(info) cycle
          
          status_rep=status_swap
          j_x=0
          k=0
          do while(k<L)
             k=k+1
             status_swap= mod(int(2,16)*status_swap,int(2,16)**L-1)
             if(status_swap.le.status_rep)then
                if(status_swap==status_rep) exit
                status_rep = status_swap
                j_x=k
             endif
          end do

          low = 1
          high = N_SzQx
          goal_not_found = .true.
          
          do while( (low.le.high).and.(goal_not_found) )
             k = (low + high)/2
             
             if ( status_rep == SpaceSzQx(k) ) then
                goal_not_found = .false.
                !S+S-
                H_NoDiag(j,i+1) = 0.5d0*exp(+cmplx(0.0d0,1.0d0)*Qx*j_x)*NormsSzQx(k)/NormsSzQx(j)
                H_NoDiag_index(j,i+1) = k
             else if ( status_rep < SpaceSzQx(k) ) then
                high = k - 1
             else
                low = k + 1
             end if
          end do
          
       end do
    end do
    



  end subroutine Get_Ham_NoDiag
!********** Building the diagonal elements of Hamiltonian **********  
  subroutine  Get_Ham_Diag(H_Diag,SpaceSzQx,N_SzQx,L)
    integer, intent(in) :: L
    integer(16),intent(in) :: N_SzQx
    complex(8), dimension(1:N_SzQx), intent(out) :: H_Diag
    integer(16), dimension(1:N_SzQx), intent(in)  :: SpaceSzQx
    integer :: i,j,j_x,k, Si,Sip, j_swapped
    integer(16) :: status_j, status_swapped, status_rep
    logical :: info
    
    H_Diag = cmplx(0.0d0)
    do j=1,N_SzQx
       status_j=SpaceSzQx(j)

       !Diagonal Term Sz Sz
       Si = S_i(status_j,0)
       do i=1,L-1
          Sip = S_i(status_j,i)
          H_Diag(j) = H_Diag(j) +Sip*Si*0.25d0
          Si=Sip
       enddo
       
       !PBC:
       H_Diag(j) = H_Diag(j) +Si*S_i(status_j,0)*0.25d0

    end do
       
  end subroutine Get_Ham_Diag



  !********** LANCZOC CHAIN **********
  subroutine Lanczos(alpha,beta,N_lanc,N_SzQx,H_Diag,H_NoDiag,H_NoDiag_index,L)
    integer, intent(in) :: N_lanc,L
    integer(16), intent(in) :: N_SzQx
    real(8), dimension(1:N_lanc), intent(out) :: alpha,beta
    complex(8),dimension( 1:N_SzQx ) :: H_Diag
    complex(8),dimension(1:N_SzQx,1:L) :: H_NoDiag
    integer(8),dimension(1:N_SzQx,1:L) :: H_NoDiag_index
    complex(8),dimension(1:N_SzQx) :: Psi, PsiPlus, PsiMinus
    integer :: i,j,k,index,iseed(1:4)
    complex(8) :: alpha_i
    
    
    !Gives random values to Psi
    iseed(1) = mod(irand(),4095)
    iseed(2)=  mod(irand(),4095)
    iseed(3)= mod(irand(),4095)
    iseed(4)= mod(irand(),4095)
    if( mod(iseed(4),2)==0) iseed(4)=iseed(4)+1
    call clarnv(2,iseed,2*N_SzQx,Psi)
    Psi=Psi/dsqrt(sum(abs(Psi)**2))
    PsiMinus=cmplx(0.0d0)
    PsiPlus=cmplx(0.0d0)
    alpha=cmplx(0.0d0)
    beta =cmplx(0.0d0)
    !print*,'Psi:',Psi
    !print*,'Hdiag',H_Diag
    
    do i=1,N_lanc-1

       
       !Computing alpha(i)
       alpha_i = sum(H_Diag*Psi*conjg(Psi))
       do j=1,N_SzQx
          do k=1,L
             !S+S-
             index = H_NoDiag_index(j,k)
             if(index.ne.0) alpha_i = alpha_i + 2.0d0*real(conjg(Psi(index))*H_NoDiag(j,k)*Psi(j))
          end do
       end do
       alpha(i) = real(alpha_i)
       
       
       !Computing (H-alpha_i)Psi-beta
       PsiPlus= (H_Diag-alpha(i))*Psi-beta(i)*PsiMinus
       
       do j=1,N_SzQx
          do k=1,L
             !S+S-
             index = H_NoDiag_index(j,k)
             if(index .ne. 0) PsiPlus(index) = PsiPlus(index) + H_NoDiag(j,k)*Psi(j)
             if(index .ne. 0) PsiPlus(j) = PsiPlus(j) + conjg(H_NoDiag(j,k))*Psi(index)
          enddo
       enddo
       beta(i+1) = dsqrt(sum(abs(PsiPlus)**2))
       PsiPlus = PsiPlus/beta(i+1)

       PsiMinus = Psi
       Psi = PsiPlus

    end do
    
    
    alpha_i = sum(H_Diag*Psi*conjg(Psi))
    do j=1,N_SzQx
       do k=1,L
          !S+S- e S-S+
          index = H_NoDiag_index(j,k)
          if(index .ne. 0)  alpha_i = alpha_i + 2.0d0*real(conjg(Psi(index))*H_NoDiag(j,k)*Psi(j))
       end do
    end do
    
    alpha(N_Lanc) = real(alpha_i)
    
  end subroutine Lanczos
  
end module sub_mod





!##########   PROGRAM   ##########

program HEIS_S1


  use sub_mod
  implicit none

  integer :: i,j,k,lwork,l_nodiag, info
  integer :: L, Sz, N_lanc
  integer(16) :: status, status_shift, N_L, N_Sz,N_Sz_all, S_status, N_SzQx,i_build
  integer(16), dimension(:), allocatable :: SpaceSz, SpaceSzQx,States_Builds
  integer(8), dimension(:,:), allocatable :: H_NoDiag_index
  real(8), dimension(:), allocatable :: momenta, NormsSzQx, NormsSz, eng, work, alpha, beta
  logical, dimension(:), allocatable :: InStates
  complex(8),dimension(:) , allocatable  ::  H_Diag
  complex(8),dimension(:,:) , allocatable  :: H_NoDiag
  real(8) :: Qx, E_gs, E_1, t1, t2
  real(8), parameter :: TwoPi = 8.d0*datan(1.d0)
  
  write(*,*) 'Insert L, 2*Sz,N_lanc'
  read(*,*) L,Sz,N_lanc
  N_L = int(2,16)**L-1
  print*, 'N_L', N_L
  allocate( momenta(1:L) )
  do i=1,L
     momenta(i) = (i-1)*TwoPi/L
  enddo


  !Creating Representatives For Sz=0
  t1 = time()
  N_Sz = 0
  open(13,file='Sz.dat')
  do status=0,N_L
     if( S_z(status,L) == Sz ) then
        status_shift = status
        do i=1, L
           status_shift = mod( int(2,16)*status_shift, int(2,16)**L-1)
           if (status_shift < status) then
              exit
           else if(status_shift == status) then
              write(13,'(I14)') status
              N_Sz = N_Sz+1
              exit
           endif 
        end do   
     end if
  enddo
  close(13)
  t2=time()
  
  print*, 'states written in file Sz.dat in ',t2-t1,' seconds'
  print*, 'Found N_Sz =',N_Sz

  !Reading Representatives for Sz=0
  allocate(SpaceSz(1:N_Sz))
  t1=time()
  open(13,file='Sz.dat',status='old',action='read')
  do i_build=1,N_Sz
     read(13,'(I14)') SpaceSz(i_build)
  end do
  close(13)
  t2=time()
  print*, 'states read from Sz.dat in',t2-t1,'seconds'



  !Array for Representative 
  allocate(NormsSz(1:N_Sz))
  allocate(InStates(1:N_Sz))

  !LOOP for momenta
  call cpu_time(t1)
  open(13,file='energies.dat')
  do i=1,L
     Qx=momenta(i)
     print*, '... Doing Qx:',Qx
     NormsSz = 0.0d0
     call Get_Norms(NormsSz,SpaceSz,N_Sz,Qx,L)
     InStates = (NormsSz > 1.0e-12)
     N_SzQx = count( InStates )
     print*, 'Qx:',Qx,'  - N:',N_SzQx
     allocate( SpaceSzQx(1:N_SzQx), NormsSzQx(1:N_SzQx) )

     !Writing the representative states
     !with non-zero norm
     j=0
     do status=1,N_Sz
        if( InStates(status) ) then
           j=j+1
           SpaceSzQx(j) = SpaceSz(status)
           NormsSzQx(j) = NormsSz(status)
        endif
     enddo


     l_nodiag = 2*L
     allocate( H_Diag( 1:N_SzQx) )
     allocate(H_NoDiag(1:N_SzQx,1:L),H_NoDiag_index(1:N_SzQx,1:L))

     print*,'step 1'
     call Get_Ham_Diag(H_Diag,SpaceSzQx,N_SzQx,L)
     print*,'step 2'
     call Get_Ham_NoDiag(H_NoDiag,H_NoDiag_index,SpaceSzQx,NormsSzQx,N_SzQx,L,Qx)
     lwork=2*N_lanc-2
     allocate(alpha(1:N_lanc),beta(1:N_lanc),eng(1:N_lanc),work(1:lwork))
     print*,'step 3'
     call Lanczos(alpha,beta,N_lanc,N_SzQx,H_Diag,H_NoDiag,H_NoDiag_index,L)
     info=0
     print*,'step 4'
     call dstev('N',N_lanc,alpha,beta(2:),eng,N_lanc,work,info)
     if(info.ne.0) then
        print*,info
        stop 'error in diag'
     end if

     
     print*, 'For Qx:',Qx, 'E min lanc:',alpha(1),alpha(2)
     do j=1,N_lanc
        write(13,*) alpha(j)
     enddo
     if(i==1) then
        E_GS = alpha(1)
        E_1 = alpha(2)
     else if(alpha(1)<E_1)then
        if( alpha(1)<E_GS ) then
           E_1=E_GS
           E_GS=alpha(1)
        else
           E_1=alpha(1)
        end if
     end if
        
     
     deallocate(eng)
     deallocate(work)
     deallocate(alpha,beta)
     deallocate(H_NoDiag)
     deallocate(H_NoDiag_index)
     deallocate(H_Diag)
     deallocate(SpaceSzQx)
     deallocate(NormsSzQx)
  end do
  call cpu_time(t2)

  open(14,file='outputs.dat')
  write(14,*) L, N_lanc, E_gs, E_1
  close(14)

  print*, 'E_gs:',E_gs
  print*, 'E_1: ',E_1
  print*, 'E_1 - E_gs :', E_1-E_gs
  print*, ' Lanczon in:',t2-t1,'seconds'


  
  
end program HEIS_S1
