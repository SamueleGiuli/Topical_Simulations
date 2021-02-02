module Funzioni
  use ConstAndKind
  implicit none
  private
  public :: Pressione, Distribuzione
contains

subroutine Pressione(nbody,mekin,L,f,pos,Pit)
  integer, intent(in) :: nbody
  real(kind=kr), intent(in) :: mekin, L
  real(kind=kr), intent(in), dimension(:,:) :: f, pos
  real(kind=kr), intent(out) :: Pit
  Pit = (2*mekin + sum(pos*f))/(3*L**3)
end subroutine Pressione

subroutine Distribuzione(vel,Distr,nbody,vdel,nvel)
  real(kind=kr), intent(in), dimension(:,:) :: vel
  integer, intent(out), dimension(:)  :: Distr
  real(kind=kr), intent(in) :: vdel
  real(kind=kr) :: V
  integer, intent(in) :: nbody,nvel
  integer :: a, k
  Do a=1,nbody
     k=0
     V= dot_product(vel(:,a),vel(:,a))
     k= int(sqrt(V)/vdel)+1
     Distr(k)=Distr(k)+1           
  end do
end subroutine Distribuzione



end module Funzioni
