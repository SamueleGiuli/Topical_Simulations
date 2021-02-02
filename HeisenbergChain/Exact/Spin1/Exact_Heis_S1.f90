module sub_mod

  implicit none

  public :: S_i, S_z, cancel_shifts, Get_Norms, Get_Ham
  contains


    
!********** S_z at position "i" state "N" **********
  integer function S_i(N,i)
    integer(8), intent(in) :: N
    integer, intent(in) :: i
    integer :: j,jp
    j=2*i
    jp=2*i+1
    if ( xor( btest(N,j) , btest(N,jp)  ) ) then
       if(btest(N,j)) then
          S_i = 1
       else
          S_i = -1
       endif
    else if( btest(N,j) ) then
       S_i = 1000
    else
       S_i = 0
    endif
  end function S_i
    
!********** Total S_z of state "N" **********
  integer function S_z(N,L)
    integer(8), intent(in) :: N
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
    integer(8), intent(out) :: status_swapped
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
  
  !********** Cancel periodic copies of the repres **********
  
  subroutine cancel_shifts(status,StatesSz,L,N_L)
    integer, intent(in) :: L,N_L
    logical, dimension(0:N_L), intent(out) :: StatesSz
    integer(8), intent(in) :: status
    integer(8) :: stat_shift
    integer :: i
    stat_shift = status
    do i=1,L
       stat_shift = mod(4*stat_shift,4**L)+(4*stat_shift/int(4**L))
       if(stat_shift == status) then
          exit
       endif
       StatesSz(stat_shift) = .false.
    enddo
    
  end subroutine cancel_shifts
    

  
  !********** Computing Norms at fixed Sz Qx **********
  subroutine Get_Norms(NormsSzQx,SzSpace,N_Sz,Qx,L)
    integer, intent(in) :: N_Sz,L
    integer(8), dimension(1:N_Sz), intent(in) :: SzSpace
    real(8), dimension(1:N_Sz), intent(out) :: NormsSzQx
    real(8), intent(in) :: Qx
    integer :: i,j,j_period,k
    integer(8) :: status_i, status_shift
    complex(8), dimension(:), allocatable :: status
    complex(8), parameter :: img = dcmplx(0.0d0,1.0d0)
    
    do i = 1,N_Sz
       status_i = SzSpace(i)
       status_shift = status_i

       !get periodicity
       do j = 1,L
          status_shift = mod(4*status_shift,4**L)+(4*status_shift/int(4**L))
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
    integer(8), intent(in) :: status_j
    integer(8), intent(out) :: status_swapped
    integer, intent(in) :: i,L
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


!********** Building the Hamiltonian **********  
  subroutine Get_Ham(H,SpaceSzQx,NormsSzQx,N_SzQx,L,Qx)
    integer, intent(in) :: N_SzQx,L
    real(8), intent(in) :: Qx
    complex(8), dimension(1:N_SzQx,1:N_SzQx), intent(out) :: H
    real(8), dimension(1:N_SzQx), intent(in) :: NormsSzQx
    integer(8), dimension(1:N_SzQx), intent(in) :: SpaceSzQx
    integer :: i,j,j_x,k, Si,Sip, j_swapped
    integer(8) :: status_j, status_swapped, status_rep
    logical :: info
    
    H = cmplx(0.0d0)
    do j=1,N_SzQx
       status_j=SpaceSzQx(j)

       !Diagonal Term Sz Sz
       Si = S_i(status_j,0)
       do i=1,L-1
          Sip = S_i(status_j,i)
          H(j,j) = H(j,j) +Sip*Si
          Si=Sip
       enddo
       
       !PBC:
       H(j,j) = H(j,j) +Si*S_i(status_j,0)
          
       
       !Non-diagonal S+ S- and H.c.
       !PBC ends at L-1
       !OBC ends at L-2
       do i=0,L-1
          call SpSm(status_j,status_swapped,i,info,L)
          if (info) cycle
          !print*, 'swap:', i, status_j, status_swapped
          status_rep=status_swapped
          info=.true.
          j_x=0
          k=0
          do while(k<L)
             k=k+1
             status_swapped= mod(4*status_swapped,4**L)+(4*status_swapped/int(4**L))
             if(status_swapped.le.status_rep)then
                if(status_swapped==status_rep) exit
                status_rep = status_swapped
                !print*, k, status_rep
                j_x=k
             endif
          end do

          

          do k=1,N_SzQx
             if(status_rep == SpaceSzQx(k) ) then
                !print*, 'in h:',status_j, status_rep, exp(+cmplx(0.0d0,1.0d0)*Qx*j_x)*NormsSzQx(j)/NormsSzQx(k)
                !here 1/2 for S+S- and 2 from S+- acting on spin 1
                H(k,j) = H(k,j)+1.0d0*exp(+cmplx(0.0d0,1.0d0)*Qx*j_x)*NormsSzQx(k)/NormsSzQx(j)
                H(j,k) = H(j,k)+1.0d0*exp(-cmplx(0.0d0,1.0d0)*Qx*j_x)*NormsSzQx(k)/NormsSzQx(j)
                exit
             end if
          end do
       end do
    end do
    
       
    

  end subroutine Get_Ham
    
    
    
          
          
  
end module sub_mod





!##########   PROGRAM   ##########

program HEIS_S1


  use sub_mod
  implicit none

  integer :: i,j,k,lwork,l_rwork
  integer :: L, Sz, N_Sz,S_status, N_L, N_SzQx
  integer(8) :: status
  logical :: info
  integer(8), dimension(:), allocatable :: SpaceSz, SpaceSzQx
  real(8), dimension(:), allocatable :: momenta, NormsSzQx, NormsSz, eng, rwork
  logical, dimension(:), allocatable :: InStates
  complex(8),dimension(:,:), allocatable :: H
  complex(8),dimension(:) , allocatable  :: work
  real(8) :: Qx, E_gs, E_1
  real(8), parameter :: TwoPi = 8.d0*datan(1.d0)
  
  write(*,*) 'Insert L, Sz'
  read(*,*) L,Sz
  N_L = 4**L-1
  allocate(InStates(0:N_L) , momenta(1:L) )
  do i=1,L
     momenta(i) = (i-1)*TwoPi/L
  enddo
!  print*, 'mom:',momenta
  
  InStates = .true.
  
  do status=0,N_L
     if( InStates(status) ) then
        S_status = S_z(status,L)
        if( S_status .ne. Sz )  then
           InStates(status)=.false.
        else
           call cancel_shifts(status,InStates,L,N_L)
        endif
     endif
     
     !print*, status, S_status, StatesSz(status)
  enddo

  N_Sz = count(InStates(0:N_L))
  print*,'Found N_Sz =',N_Sz
  allocate(SpaceSz(1:N_Sz))
  SpaceSz=0
  
  i=0
  do status=0,N_L
     if( InStates(status) ) then
        i=i+1
        SpaceSz(i) = status
        !print*, status, ' è dentro SzSpace all indice ',i
     endif
  enddo

  !Now the Sz states have been found
  deallocate(InStates)

  !Finding StatesSzQx
  allocate(NormsSz(1:N_Sz))
  allocate(InStates(1:N_Sz))
  
  open(13,file='energies.dat')
  do i=1,L
     Qx=momenta(i)
     print*, '... Doing Qx:',Qx
     NormsSz = 0.0d0
     call Get_Norms(NormsSz,SpaceSz,N_Sz,Qx,L)
     InStates = (NormsSz > 1.0e-14)
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

     
     allocate( H( 1:N_SzQx , 1:N_SzQx ) )
     H = cmplx(0.0d0)
     call Get_Ham(H,SpaceSzQx,NormsSzQx,N_SzQx,L,Qx)
     if(Qx==0) then
        print*, 'H'
        print*, H(1,1), H(1,2)
        print*, H(2,1), H(2,2)
     end if
        
     lwork=2*N_SzQx
     l_rwork = 3*N_SzQx-1
     allocate(eng(1:N_SzQx),work(1:lwork),rwork(1:l_rwork))
     call zheev ('N','U',N_SzQx,H,N_SzQx,eng,work,lwork,rwork,info)
     if (info) stop 'error diagonaliz'

     print*, 'For Qx:',Qx, 'E min:',eng(1)
     do j=1,N_SzQx
        write(13,*) eng(j)
     enddo
     if(i==1) then
        E_GS = eng(1)
        E_1 = eng(2)
     else if(eng(1)<E_1)then
        if(eng(1)<E_GS) then
           E_1=E_GS
           E_GS=eng(1)
        else
           E_1=eng(1)
        end if
     end if
        
     
     deallocate(eng)
     deallocate(work)
     deallocate(rwork)
     deallocate(H)
     deallocate(SpaceSzQx)
     deallocate(NormsSzQx)
  end do
  



  print*, 'E_gs:',E_gs
  print*, 'E_1: ',E_1
  print*, 'E_1 - E_gs per site:',(E_1-E_gs)/L
  open(13,file='outputs.dat')
  write(13,*) L, E_GS, E_1
  close(13)

  
  
end program HEIS_S1
