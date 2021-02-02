module Potenziale
  use ConstAndKind
  implicit none
  private
  public :: interazione
contains

  subroutine interazione(pos,nbody,f,upot,L,del,stat,gr,nh,pres,Pit)
    
real(kind=kr), intent(in), dimension(:,:) :: pos
integer, intent(in) :: nbody
integer, intent(in) :: nh
real(kind=kr), intent(in) :: del,L
real(kind=kr), intent(out), dimension(:) :: gr
real(kind=kr), intent(out) :: upot, Pit
real(kind=kr), intent(out), dimension(:,:) :: f
real(kind=kr), dimension(size(pos,1)) :: posij,posI,Fagg
real(kind=kr) :: rij,ig,dr
integer :: i,j,a,k
logical, intent(in) :: stat, pres
upot = 0
f = 0
k=0
dr=0
do i=1,nbody
do j=1,nbody
  if( i==j ) cycle !passa al prossimo, si potrebbe ottimizzare considerando Fij=Fji
  posij = pos(:,i)-pos(:,j)
  posI=mod(posij,L)


  
  do a=1,3
     if(posI(a).gt.(L/2) ) posI(a)=posI(a)-L
     if(posI(a).lt.(-L/2)) posI(a)=posI(a)+L
  end do



  rij=sqrt( dot_product(posI,posI) ) !modulo
  upot = upot + 4*(rij**(-12)-rij**(-6))
  Fagg= 24*(2.*rij**(-14)-rij**(-8))*posI
  f(:,i) = f(:,i) +Fagg ! forza atomo i

  
  !Calcolo della G(r)
  
 if( (i.lt.j).and.(stat))  then               
    if(rij.lt.(L/2)) then
       k=int(rij/del)+1
       gr(k)=gr(k)+2
    end if             
 end if

  !Calcolo della Pressione dovuta al viriale
 if(pres.and.(i.lt.j)) then
    Pit=Pit+dot_product(posI,Fagg)
 end if
 
  
end do
end do
upot = upot/2
end subroutine interazione
  
end module Potenziale
