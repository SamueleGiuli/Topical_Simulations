module sub_mod

  
  implicit none
  integer :: L,Sz,N_lanc,l_nodiag
  real(8) :: b
  public :: S_i, S_z, Get_Norms, Get_Ham_Diag, Get_Ham_NoDiag, Lanczos, init,L,Sz,N_lanc, l_nodiag
  contains

    subroutine init()
      
      write(*,*) 'Insert L, Sz,N_lanc,b  '
      read(*,*) L,Sz,N_lanc,b

    end subroutine init
      
!********** S_z at position "i" state "N" **********
  integer function S_i(N,i)
    integer(16), intent(in) :: N
    integer, intent(in) :: i
    integer :: j,jp
    j=2*L-1-2*i
    jp=j-1
    if ( xor( btest(N,j) , btest(N,jp)  ) ) then
       if(btest(N,j)) then
          S_i = -1
       else
          S_i = +1
       endif
    else if( btest(N,j) ) then
       S_i = 1000
    else
       S_i = 0
    endif
  end function S_i
    
!********** Total S_z of state "N" **********
  integer function S_z(N,L)
    integer(16), intent(in) :: N
    integer, intent(in) :: L
    integer :: i, Stot, s_trial
    Stot=0
    do i=0,L-1
       s_trial = S_i(N,i)
       
       !finding a 11 state means not counting it
       if (s_trial > 2) then
          Stot = 2*L
          exit
       endif
       
       Stot = Stot + S_i(N,i)
    enddo
    S_z = Stot
  end function S_z

  !********** Set S_i to Val in status_swapped **********
  
  subroutine SetS_i(status_swapped,i,Val)
    integer(16), intent(out) :: status_swapped
    integer,intent(in) :: i,Val
    integer :: j,jp
    j=2*L-1-2*i
    jp=j-1
    if (Val==0)then
       status_swapped=ibclr(status_swapped,j)
       status_swapped=ibclr(status_swapped,jp)
    else if(Val==1)then
       status_swapped=ibclr(status_swapped,j)
       status_swapped=ibset(status_swapped,jp)
    else if(Val==-1)then
       status_swapped=ibset(status_swapped,j)
       status_swapped=ibclr(status_swapped,jp)
    else
       stop 'Error in Setting S_i'
    end if
    
  end subroutine SetS_i   
  
  !********** Computing Norms at fixed Sz Qx **********
  subroutine Get_Norms(NormsSzQx,SzSpace,N_Sz,Qx)
    integer(16),intent(in) :: N_Sz
    integer(16), dimension(1:N_Sz), intent(in) :: SzSpace
    real(8), dimension(1:N_Sz), intent(out) :: NormsSzQx
    real(8), intent(in) :: Qx
    integer :: j
    integer(16) :: status_i, status_shift, i,j_period,k
    complex(8), dimension(:), allocatable :: status
    
    do i = 1,N_Sz
       status_i = SzSpace(i)
       status_shift = status_i

       !get periodicity
       do j = 1,L
          status_shift = mod(int(4,16)*status_shift,int(4,16)**L-1)
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
             status(k) = status(k)+exp( Qx*j*dcmplx(0.d0,1.d0) )
          enddo
          NormsSzQx(i) = dsqrt(sum(abs(status)**2)  )
          deallocate(status)
       endif
    end do
  end subroutine Get_Norms

!********** Applying the S+ S- operator at site "i" to "status_j" **********
  
  subroutine SpSm(status_j,status_swapped,i,info)
    integer(16), intent(in) :: status_j
    integer(16), intent(out) :: status_swapped
    integer, intent(in) :: i
    logical, intent(out) :: info
    integer :: Si, Sip, ip
    ip=mod(i+1,L)
    Si=S_i(status_j,i)
    Sip=S_i(status_j,ip)
    info=.false.
    status_swapped=status_j
    
    select case (Si)
    case (1)
       info=.true.
    case (0)
       
       select case(Sip)
       case(1) ! 0 | 1
          !print*,'0|1'
          call SetS_i(status_swapped,i,1)
          call SetS_i(status_swapped,ip,0)
       case(0) ! 0 | 0
          !print*,'0|0'
          call SetS_i(status_swapped,i,1)
          call SetS_i(status_swapped,ip,-1)
       case(-1)
          info=.true.
       case DEFAULT
          stop 'error'
       end select
       
    case (-1)
       
       select case(Sip)
       case(1) !-1 | 1
          !print*,'-1|1'
          call SetS_i(status_swapped,i,0)
          call SetS_i(status_swapped,ip,0)
       case(0) !-1 | 0
          !print*,'-1|0'
          call SetS_i(status_swapped,i,0)
          call SetS_i(status_swapped,ip,-1)
       case(-1)
          info=.true.
       case DEFAULT
          stop 'error'
       end select
       
    case DEFAULT
       stop 'error'
    end select 
    
  
  end subroutine SpSm


!********** Building the non-diagonal elements of Hamiltonian **********  
  subroutine Get_Ham_NoDiag(H_NoDiag,H_NoDiag_index,SpaceSzQx,NormsSzQx,N_SzQx,Qx,b)
    integer(16),intent(in) :: N_SzQx
    real(8), intent(in) :: Qx,b
    integer(16), dimension(1:N_SzQx), intent(in) :: SpaceSzQx
    real(8), dimension(1:N_SzQx), intent(in) :: NormsSzQx
    integer(8), dimension(1:N_SzQx,1:l_nodiag), intent(out) :: H_NoDiag_index
    complex(8), dimension(1:N_SzQx,1:l_nodiag), intent(out) :: H_NoDiag
    integer(16) :: j,k_rep, low, high
    integer :: i,ip,j_x,k
    integer(16) :: status_j, status_swap, status_rep,status_swap_2
    logical :: info, goal_not_found
    
    H_NoDiag_index = 0
    H_NoDiag = dcmplx(0.0d0,0.0d0)
    
    do j=1,N_SzQx
       status_j = SpaceSzQx(j)
       do i=0,L-1

!********* Trying S+S-
          call SpSm(status_j,status_swap,i,info)
          if(info) cycle
          status_swap_2 = status_swap
          ip=mod(i+1,L)
          status_rep=status_swap
          j_x=0
          k=0
          do while(k<L)
             k=k+1
             status_swap= mod(int(4,16)*status_swap,int(4,16)**L -1)
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
             k_rep = (low + high)/2
             if ( status_rep == SpaceSzQx(k_rep) ) then
                goal_not_found = .false.
                ! ( S+S- ) + (SzSz S+S-) + (S+S-SzSz)
                H_NoDiag(j,i+1) = exp(-dcmplx(0.0d0,1.0d0)*Qx*j_x)*NormsSzQx(k_rep)/NormsSzQx(j)*&
                     dcmplx(1.0 - b*(S_i(status_swap_2,i)*S_i(status_swap_2,ip) + S_i(status_j,i)*S_i(status_j,ip) ) )                
                H_NoDiag_index(j,i+1) = k_rep
                
             else if ( status_rep < SpaceSzQx(k_rep) ) then
                high = k_rep - 1
             else
                low = k_rep + 1
             end if
             
          end do
         ! if(goal_not_found) print*, 'ERROR IN FINDING THE REPR STATE',j,i,status_rep,status_j


!********* Trying (S+S-)**2

          !ERRORE DI PRIMA
          !call SpSm(status_j,status_swap,i,info)
          !call SpSm(status_j,status_swap,i,info)
          call SpSm(status_swap_2,status_swap,i,info)
          if(info) cycle
          
          status_rep=status_swap
          j_x = 0
          k = 0
          do while(k<L)
             k=k+1
             status_swap= mod(int(4,16)*status_swap,int(4,16)**L-1)
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
             k_rep = (low + high)/2
             
             if ( status_rep == SpaceSzQx(k_rep) ) then
                goal_not_found = .false.
                ! ( S+S- )**2
                H_NoDiag(j,i+1+L) = -b*exp(-dcmplx(0.0d0,1.0d0)*Qx*j_x)*NormsSzQx(k_rep)/NormsSzQx(j)
                H_NoDiag_index(j,i+1+L) = k_rep
                
             else if ( status_rep < SpaceSzQx(k_rep) ) then
                high = k_rep - 1
             else
                low = k_rep + 1
             end if
          end do

          !if(goal_not_found) print*, 'ERROR IN FINDING THE REPR STATE 2',j,i

       end do
       
    end do
    



  end subroutine Get_Ham_NoDiag
!********** Building the diagonal elements of Hamiltonian **********  
  subroutine  Get_Ham_Diag(H_Diag,SpaceSzQx,N_SzQx,b)
    real(8), intent(in) :: b
    integer(16),intent(in) :: N_SzQx
    complex(8), dimension(1:N_SzQx), intent(out) :: H_Diag
    integer(16), dimension(1:N_SzQx), intent(in)  :: SpaceSzQx
    integer(16) :: j, Si,Sip,status_j
    integer :: i
    
    H_Diag = dcmplx(0.0d0)
    do j=1,N_SzQx
       status_j=SpaceSzQx(j)

       Si = S_i(status_j,0)
       do i=1,L-1
          Sip = S_i(status_j,i)
          !Diagonal Term Sz Sz and (S_z S_z)^2
          H_Diag(j) = H_Diag(j) +Sip*Si - b*(Sip*Si)**2

          !Diagonal Term S-S+S+S-
          if( Si<1 .and. Sip>-1) then
             H_Diag(j) = H_Diag(j) - b
             !Diagonal Term S+S-S-S+
          end if
          if( Si>-1 .and. Sip<1) then
             H_Diag(j) = H_Diag(j) - b
          end if
             
          Si=Sip
       enddo
       !PBC:
       Sip = S_i(status_j,0)
       H_Diag(j) = H_Diag(j) +Si*Sip -b*(Si*Sip)**2
       
       !Diagonal Term S-S+S+S-
       if( Si<1 .and. Sip>-1) then
          H_Diag(j) = H_Diag(j) - b
       end if
       !Diagonal Term S+S-S-S+
       if( Si>-1 .and. Sip<1) then
          H_Diag(j) = H_Diag(j) - b
       end if
      
    end do
       
  end subroutine Get_Ham_Diag



  !********** LANCZOC CHAIN **********
  subroutine Lanczos(alpha,beta,N_SzQx,H_Diag,H_NoDiag,H_NoDiag_index)
    integer(16), intent(in) :: N_SzQx
    real(8), dimension(1:N_lanc), intent(out) :: alpha,beta
    complex(8),dimension( 1:N_SzQx ) :: H_Diag
    complex(8),dimension(1:N_SzQx,1:l_nodiag) :: H_NoDiag
    integer(8),dimension(1:N_SzQx,1:l_nodiag) :: H_NoDiag_index
    complex(8),dimension(1:N_SzQx) :: Psi, PsiPlus, PsiMinus
    integer(16) :: i,j,k,index
    integer :: iseed(1:4)
    complex(8) :: alpha_i
    
    
    !Gives random values to Psi
    iseed(1) = mod(irand(),4095)
    iseed(2) = mod(irand(),4095)
    iseed(3) = mod(irand(),4095)
    iseed(4) = mod(irand(),4095)
    if( mod(iseed(4),2)==0) iseed(4)=iseed(4)+1
    call clarnv(2,iseed,2*N_SzQx,Psi)
    Psi=Psi/dsqrt(sum(abs(Psi)**2))
    PsiMinus= dcmplx(0.0d0)
    PsiPlus = dcmplx(0.0d0)
    alpha   = dcmplx(0.0d0)
    beta    = dcmplx(0.0d0)
    !print*,'Psi:',Psi
    !print*,'Hdiag',H_Diag
    
    do i=1,N_lanc-1

       
       !Computing alpha(i)
       alpha_i = sum(H_Diag*Psi*conjg(Psi))
       do j=1,N_SzQx
          do k=1,l_nodiag
             !S+S-
             index = H_NoDiag_index(j,k)
             if(index.ne.0) alpha_i = alpha_i + 2.0d0*real(conjg(Psi(index))*H_NoDiag(j,k)*Psi(j))
          end do
       end do
       alpha(i) = real(alpha_i)
       
       
       !Computing (H-alpha_i)Psi-beta
       PsiPlus= (H_Diag-alpha(i))*Psi-beta(i)*PsiMinus
       
       do j=1,N_SzQx
          do k=1,l_nodiag
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
       do k=1,l_nodiag
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

  integer :: i,j,lwork, info
  integer(16) :: status, status_shift, N_L, N_Sz, N_SzQx,i_build
  integer(16), dimension(:), allocatable :: SpaceSz, SpaceSzQx
  integer(8), dimension(:,:), allocatable :: H_NoDiag_index
  real(8), dimension(:), allocatable :: momenta, NormsSzQx, NormsSz, eng, work, alpha, beta
  logical, dimension(:), allocatable :: InStates
  complex(8),dimension(:) , allocatable  ::  H_Diag
  complex(8),dimension(:,:) , allocatable  :: H_NoDiag
  real(8) :: Qx, E_gs, E_1, t1, t2
  real(8), parameter :: TwoPi = 8.d0*datan(1.d0)

  call init()
  N_L = int(4,16)**L-1
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
           status_shift = mod( int(4,16)*status_shift, int(4,16)**L-1)
           if (status_shift < status) then
              exit
           else if(status_shift == status) then
              write(13,'(I16)') status
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
     read(13,'(I16)') SpaceSz(i_build)
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
     call Get_Norms(NormsSz,SpaceSz,N_Sz,Qx)
     InStates = (NormsSz > 1.0e-12)
     N_SzQx = count( InStates )
     print*, 'Qx:',Qx,'  - N:',N_SzQx
     allocate( SpaceSzQx(1:N_SzQx), NormsSzQx(1:N_SzQx) )

     !Writing the representative states
     !with non-zero norm
     j=0
     do status=1, N_Sz
        
        if( InStates(status) ) then
           j=j+1
           SpaceSzQx(j) = SpaceSz(status)
           NormsSzQx(j) = NormsSz(status)
        endif
        
     enddo


     l_nodiag = 2*L
     allocate( H_Diag( 1:N_SzQx) )
     allocate(H_NoDiag(1:N_SzQx,1:l_nodiag),H_NoDiag_index(1:N_SzQx,1:l_nodiag))

     print*,'step 1'
     call Get_Ham_Diag(H_Diag,SpaceSzQx,N_SzQx,b)
     print*,'step 2'
     call Get_Ham_NoDiag(H_NoDiag,H_NoDiag_index,SpaceSzQx,NormsSzQx,N_SzQx,Qx,b)
     lwork=2*N_lanc-2
     allocate(alpha(1:N_lanc),beta(1:N_lanc),eng(1:N_lanc),work(1:lwork))
     print*,'step 3'
     print*,'diag:'
     !do status=1,N_SzQx
     ! print*, SpaceSzQx(status),H_Diag(status)
     !enddo
     !print*,'No_Diag:'
     ! do status=1,N_SzQx
     !    do i_build=1,l_nodiag
     !       if(H_NoDiag_index(status,i_build) .ne. 0)  then
     !          print*, SpaceSzQx(status),SpaceSzQx(H_NoDiag_index(status,i_build)),H_NoDiag(status,i_build) 
     !       end if
     !    end do
     ! end do
     call Lanczos(alpha,beta,N_SzQx,H_Diag,H_NoDiag,H_NoDiag_index)
     info=0
     print*,'step 4'
     call dstev('N',N_lanc,alpha,beta(2:),eng,N_lanc,work,info)
     if(info.ne.0) then
        print*,info
        stop 'error in diag'
     end if

     
     print*, 'For Qx:',Qx, 'E min lanc:',alpha(1)
     do j=1,N_lanc
        write(13,*) alpha(j)
     enddo
     if(i==1) then
        E_GS = alpha(1)
        E_1 = alpha(2)
     else if(alpha(1)<E_1)then
        if(alpha(1)<E_GS) then
           E_1 = E_GS
           E_GS = alpha(1)
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
  write(14,*) L, b, N_lanc, E_gs, E_1
  close(14)

  print*, 'E_gs:',E_gs
  print*, 'E_1: ',E_1
  print*, 'E_1 - E_gs :', E_1-E_gs
  print*, ' Lanczon in:',t2-t1,'seconds'


  
  
end program HEIS_S1
