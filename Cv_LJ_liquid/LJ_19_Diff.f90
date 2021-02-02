module modulo1
use ConstAndKind
  implicit none
private
public ::  Termostat

contains

subroutine Termostat(Ttarget,nbody,vel,mekin,termotype,freq,it)
real(kind=kr), intent(out), dimension(:,:) :: vel
real(kind=kr), intent(in) :: mekin, Ttarget, freq
real(kind=kr) :: Scaling, Sigma,Random
logical, intent(in) ::termotype
integer, intent(in) :: nbody, it
integer :: i,j,t,modulo

if(termotype) then
   t=int(freq**(-1))
   modulo=mod(it,t)
!Qui scaling delle velocità
   if(modulo .eq. 0) then
      Scaling=sqrt((3*nbody*Ttarget)/(2*mekin))
      vel=vel*Scaling
   end if

else
   
!Qui termostato di andersen
   do i=1, nbody
      CALL RANDOM_NUMBER(Random)
      if(Random .lt. freq) then
         Sigma= sqrt(Ttarget)
         do j=1,3
            vel(j,i)=gaus(Sigma,0.d0)
         end do
      end if
   end do
   
end if
 
end subroutine Termostat  
  
end module modulo1

program corpi3d
use ConstAndKind, rk=>kr
use Modulo1  !stesso kind (signif.)
use Potenziale
use Funzioni
implicit none
integer :: nstep,k,it,nbody,nsave,ios,i,j,jmol,a,b,nh,nvel !nsave=numdatiout
real(kind=rk), dimension(:),allocatable :: gr
integer, dimension(:),allocatable :: distr
real(kind=rk) :: dt,mepot,mekin,massa=1.,alfa=1.,vmax,L,Lcella,del,rad,dvol,Ttarget, Pit, vdel, vmaxdist, freq
real(kind=rk),dimension(:,:),allocatable :: pos, poszero, pospbc
real(kind=rk),dimension(:),allocatable :: velcm
real(kind=rk),dimension(:,:),allocatable :: ekin,vel,f
logical :: dyn,anneal,cont,ret,stat,Tcontrol,vstat,pres,TermoType
character (len=3) :: celltype, cellnum
character (len=8) :: fname

!PI=3.14159265358979323846
del=0

write(unit=*,fmt="(a)",advance="no") "Continuation run?"
read*, cont

write(unit=*,fmt="(a)",advance="no")"dynamics?"
read*, dyn
if (dyn) then
  write(unit=*,fmt="(a)",advance="no")"annealing?"
  read*, anneal
  if (anneal) then
    write(unit=*,fmt="(a)",advance="no")"scaling parameter?"
    read*, alfa
 else
    
    write(unit=*,fmt="(a)",advance="no")"Do you want to obtain a fixed temperature?"
    read*, Tcontrol
    
    if(Tcontrol) then
       write(unit=*,fmt="(a)",advance="no")"Which temperature do you want?"
       read*, Ttarget
       write(unit=*,fmt="(a)",advance="no")"Enter .t. for Temperature rescaling or .f. for Anderson Termostat"
       read*, TermoType
       write(unit=*,fmt="(a)",advance="no")"Enter the frequency of the Termostat"
       read*,freq
    end if
 end if
end if


write(unit=*,fmt="(a)",advance="no")"Lenght of periodic cubic side?"
read*, L

write(unit=*,fmt="(a)",advance="no")"n of bodies?"
read*, nbody

if(.not.cont) then
write(unit=*,fmt="(a)",advance="no")"Which type of cell?"
read*, celltype

write(unit=*,fmt="(a)",advance="no")"Lenght of the primitive cell"
read*, Lcella
end if




!Scelgo se avere g(r)
write(unit=*,fmt="(a)",advance="no") "Compute g(r)?"
read*, stat
if (stat.eqv..true.) then
write(unit=*,fmt="(a)",advance="no")"How many hits for the g(r)?"
read*, nh
del=L/(2*nh)
end if


write(unit=*,fmt="(a)",advance="no")"how many temporal step in output?"
read*, nsave

allocate(vel(3,nbody))
allocate(ekin(3,nbody))
allocate(pos(3,nbody))
allocate(poszero(3,nbody))
allocate(pospbc(3,nbody))
allocate(f(3,nbody))
allocate(velcm(3))
allocate(gr(nh))
gr=0

write(unit=*,fmt="(a)",advance="no")"time step: "
read*,dt
write(unit=*,fmt="(a)",advance="no")"n.step: "
read*,nstep
!print*,"mass of the bodies: "
!read*,massa
write(unit=*,fmt="(a)",advance="no") "Return the atoms inside the periodic cell?"
read*, ret

write(unit=*,fmt="(a)",advance="no") "Compute velocity distribution?"
read*, vstat

if(vstat) then
   
   write(unit=*,fmt="(a)",advance="no") "How many hit for the distribution?"
   read*, nvel

   write(unit=*,fmt="(a)",advance="no") "Max Velocity expected for the distribution?"
   read*, vmaxdist
   
   vdel=vmaxdist/nvel
   allocate(distr(nvel))
   distr=0
end if


write(unit=*,fmt="(a)",advance="no") "Do you want to compute the Pressure?"
read*, pres
Pit=0

vel=0.


if (cont) then
  read(unit=11) pos,vel
else

   write(cellnum,'(I3)') nbody
   fname = celltype//cellnum//".d"

   if(nbody<100) then
      fname(4:4) = "0"
      if (nbody<10) fname(5:5) = "0"
   end if
   if(celltype.eq."SC") fname(3:3) = "0"
   write(unit=*,fmt="(8a)",advance="no") fname
   
   open(10,file=fname)
   
  do
    read(unit=10,fmt=*,iostat=ios) pos 
    if (ios<=0) exit
 end do
 
 pos=pos*Lcella
 
  if (dyn) then
    call random_number(vel)
    write(unit=*,fmt="(a)",advance="no")"vmax? "
    read*, vmax 
    vel=2*vmax*(vel-0.5)
    do i=1,3
      velcm(i) = sum(vel(i,:))/nbody
      vel(i,:) = vel(i,:) - velcm(i)
    end do
  end if
end if

poszero=pos
gr=0

write(unit=15,fmt=*) 0, 0
call interazione(pos,nbody,f,mepot,L,del,.false.,gr,nh,pres,Pit)


do it = 1,nstep
   if (dyn) then
      Pit=0 !Pressione a zero prima di ogni step
    do i = 1,nbody
      pos(:,i) = pos(:,i) + vel(:,i) * dt + 0.5* f(:,i)/massa * dt**2
      vel(:,i) = vel(:,i) + 0.5 * dt * f(:,i)/massa !veloc verlain
    end do
    call interazione(pos,nbody,f,mepot,L,del,stat,gr,nh,pres,Pit) !serve per finire di calcolare v+
    
    do i=1,nbody
      vel(:,i) = vel(:,i) + 0.5 * dt * f(:,i)/massa
      ekin(:,i) = 0.5 * massa * (vel(:,i))**2 !vettore perchè contributo coordinate
   end do

   
   mekin = sum(ekin)
   

   
     !Calcolo coefficiente di diffusione
   if(stat) write(unit=15,fmt=*) it, (sum((pos-poszero)**2)/nbody), (sum((pos-poszero)**2/nbody))/(6*it)
   

   
   !Riscalo le temperature
   if(Tcontrol)  then
      call Termostat(Ttarget,nbody,vel,mekin,termotype,freq,it)
   end if
   

    !Applico le PBC
    if(ret) then
       do i=1,nbody
          do b=1,3
             pospbc(b,i)=mod(pos(b,i),L)
             if(pospbc(b,i)<0) pospbc(b,i)=pospbc(b,i)+L
          end do
       end do
    end if

! VECCHIO CALCOLO PRESSIONE    
!!!! call Pressione(nbody,mekin,L,f,pospbc,Pit)
! NUOVO CALCOLO PRESSIONE
Pit=(Pit+2*mekin)/(3*L**3)

    
    if (mod(it,nstep/nsave).eq.0) then !non scrivo ad ogni passo ma solo in modo da averne nsave a distanze uguali
      write(unit=1,fmt=*)it,it*dt,pos,vel
      write(unit=2,fmt=*)it,mekin,mepot,mekin+mepot, Pit
      write(unit=7,fmt=*) nbody
      write(unit=7,fmt=*)
      do jmol =1,nbody
         if(ret) then
            write(unit=7,fmt=*) "Ar",pospbc(:,jmol)
         else
            write(unit=7,fmt=*) "Ar",pos(:,jmol)

         end if
      end do

    end if
    if (anneal) vel=alfa*vel
 else
    do i = 1,nbody
      pos(:,i) = pos(:,i) + f(:,i)/massa * dt**2
    end do
    call interazione(pos,nbody,f,mepot,L,del,.false.,gr,nh,.false.,Pit)
    if (mod(it,nstep/nsave).eq.0) write(unit=1,fmt=*)it,pos
    if (mod(it,nstep/nsave).eq.0) write(unit=2,fmt=*)it,mepot
 end if
 

end do !Fine evoluzione temporale

dvol=0
rad=0

if(vstat) then
   call Distribuzione(vel,distr,nbody,vdel,nvel)

   do k=1,nvel
      write(unit=13,fmt=*) vdel*(k-0.5), Distr(k)
   end do
end if

!Scrittura G(r)
if(stat.eqv..true.) then
   do k=1,nh
      rad=del*(k+0.5)
      dvol=4*PI*rad**2
      gr(k)=gr(k)*(L**3)/(dvol*nstep*nbody**2)   
      write(unit=12,fmt=*) rad,gr(k) 
   end do

end if





do i=1,nbody
   do j=i+1,nbody
      write(unit=3,fmt=*) i,j,sqrt(dot_product(pos(:,i)-pos(:,j),pos(:,i)-pos(:,j)))
   end do
end do

write(unit=4) pos,vel

end program corpi3d

!ACLUNE INFORMAZIONI UTILI

!unità 1  : it,it*dt,pos,vel
!unità 2  : it,mekin,mepot,mekin+mepot, Pit
!unità 3  : i,j,Rij
!unità 4  : pos, vel per run di continuazione da linkare a unità 11
!unità 7  : nbody   E POI
!           "Ar", pos OPPURE pospbc
!unità 12 : rad, gr(rad)
!unità 13 : vel, Distr(vel)
!unità 15 : it, D(it)


