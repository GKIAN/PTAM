program main
  use ptam
  implicit none

  integer, parameter :: MK = 4

  real(MK) :: dk = 0.010_MK
  real(MK) :: kc = 17.0_MK, kl = 28.5_MK
  integer ::  npt = 10, nkc, nkl
  real(MK) :: rptam, r

  nkc = int(kc / dk) + 1
  nkl = int(kl / dk) + 1

  call ptamInit(npt, nkc, nkl)
  rptam = ptamRun(exfunc, dk)
  call ptamFinal()

  r = exfunc()
  write(*, '(A, G0)') 'The PTAM result is: ', rptam
  write(*, '(A, G0)') 'The analytical result is: ', r
  write(*, '(A, G0)') 'The relative error is: ', (rptam - r)/r

  contains

    real(MK) function exfunc(k) result(res)
      real(MK), intent(in), optional :: k
      real(MK) :: a = 0.001_MK, b = 3.0_MK
      if(present(k)) then
        res = exp( - a * k) * bessel_j1(b * k)
      else
        res = 1.0_MK / b * (1.0_MK - a / sqrt(a * a + b * b))
      end if
    end function exfunc

end program main
