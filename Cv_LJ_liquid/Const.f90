module ConstAndKind
  implicit none
  public gaus
  
  integer, parameter :: kr=selected_real_kind(12)
  real(kind=kr) ::  PI=3.14159265358979323846
  real(kind=kr) ::  twoPI=6.28318530717958647692

  contains
function gaus(s,m) result(g)  !Metodo di Box-Muller
  implicit none
    real(kind=kr), intent(in) :: s,m ! input
    real(kind=kr) :: g, a, b ! output
    CALL RANDOM_NUMBER(a)
    CALL RANDOM_NUMBER(b)
    g = s*(sqrt(-2.0 * log(a)) * cos(twoPI*b)) +m
end function gaus
  
  
end module ConstAndKind
