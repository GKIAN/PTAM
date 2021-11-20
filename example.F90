!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Thu 30 Sep 2021 10:33:09 AM CST
!-------------------------------------------------------------------------------

program main

  implicit none
  integer, parameter :: MK = 4

  real(kind = MK), parameter :: dk = 0.010_MK, kc = 17.0_MK, kl = 28.5_MK
  integer, parameter :: npt = 10

  real(kind = MK) :: k, S(3), dS, vpt(npt)
  integer :: ipt
  logical :: ga = .false.

  real(kind = MK) :: rptam, ranal

  ranal = exfunc()

  k = dk
  S = 0.0_MK
  do while(.true.)
    dS = exfunc(k) * dk
    S(2) = S(2) + dS
    if(k > kc) exit
    S(1) = S(2)
    k = k + dk
  end do

  ipt = 0
  do while(k < kl)
    k = k + dk
    dS = exfunc(k) * dk
    S(3) = S(2) + dS

    ga = ptamGather(S, ipt, vpt)
    if(ga) exit

    S(1) = S(2)
    S(2) = S(3)
  end do

  rptam = ptamReduce(ipt, vpt)

  write(*, '(A, G0)') 'The PTAM result is: ', rptam
  write(*, '(A, G0)') 'The analytical result is: ', ranal
  write(*, '(A, G0)') 'The relative error is: ', (rptam - ranal)/ranal

  contains

    real(kind = MK) function exfunc(k) result(res)
      real(kind = MK), intent(in), optional :: k
      real(kind = MK) :: a = 0.001_MK, b = 3.0_MK
      if(present(k)) then
        res = exp( - a * k) * bessel_j1(b * k)
      else
        res = 1.0_MK / b * (1.0_MK - a / sqrt(a * a + b * b))
      end if
    end function exfunc

    logical function ptamGather(S, ipt, vpt) result(gotAll)
      real(kind = MK), intent(in) :: S(:)
      integer, intent(inout) :: ipt
      real(kind = MK), intent(inout) :: vpt(:)
      real(kind = MK) :: toput, a(3)
      logical :: ifput
      ifput = .false.
      if(S(1) == S(2) .and. S(2) == S(3)) then
        ifput = .true.
        toput = S(1)
      else if((S(2) - S(1)) * (S(2) - S(3)) > 0.0_MK .or. S(2) == S(3)) then
        ifput = .true.
        a(1) = 2.0_MK * S(3) - 4.0_MK * S(2) + 2.0_MK * S(1)
        if(a(1) == 0.0_MK) then
          ! if there is only a tiny difference between S(1) with S(2) when S(2)
          !&  is equal to S(3)
          toput = S(2)
        else
          a(2) = 4.0_MK * S(2) - S(3) - 3.0_MK * S(1)
          a(3) = S(1)
          toput = a(3) - a(2) * a(2) / (4.0_MK * a(1))
        end if
      end if
      if(ifput) then
        ipt = ipt + 1
        if(ipt > npt) then
          vpt(1:npt - 1) = vpt(2:npt)
          ipt = npt
        end if
        vpt(ipt) = toput
      end if
      gotAll = (ipt >= npt)
    end function ptamGather

    real(kind = MK) function ptamReduce(ipt, vpt) result(r)
      integer, intent(in) :: ipt
      real(kind = MK), intent(inout) :: vpt(:)
      integer :: i, j
      do i = ipt - 1, 1, - 1
        do j = 1, i, 1
          vpt(j) = (vpt(j) + vpt(j + 1)) / 2.0_MK
        end do
      end do
      r = vpt(1)
    end function ptamReduce

end program main

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
