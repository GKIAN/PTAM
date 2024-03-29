!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Wed 02 Jun 2021 10:36:07 AM CST
!-------------------------------------------------------------------------------

module ptam

  implicit none
  private

  integer, parameter :: MK = 4

  integer :: numPTs = 10
  logical :: isFlat = .false., isNotSo = .true.
  integer :: nkCriti = 1, nkLimit = 9999
  real(MK) :: dk = 0.0_MK

  real(MK), allocatable :: valuePTs(:)

  public ptamInit, ptamFinal, ptamRun

  contains

    subroutine ptamInit(dk_, npt, nkc, nkl)
      real(MK), intent(in) :: dk_
      integer, intent(in), optional :: npt, nkc, nkl

      dk = dk_
      if(present(npt)) numPTs = npt
      if(present(nkc)) nkCriti = nkc
      if(present(nkl)) nkLimit = nkl

      isFlat = .false.
      isNotSo = .true.
      allocate(valuePTs(numPTs))
      valuePTs = 0.0_MK
    end subroutine ptamInit

    subroutine ptamFinal()
      if(allocated(valuePTs)) deallocate(valuePTs)
    end subroutine ptamFinal

    real(MK) function ptamRun(func) result(res)
      real(MK), external :: func
      real(MK) :: k(3), S(3), kExtre
      integer :: iPT = 0
      integer i, j

      S(1) = 0.0_MK
      do i = 1, nkCriti
        k(1) = (i - 1) * dk
        S(1) = S(1) + func(k(1)) * dk
      end do

      k(2) = k(1) + dk
      k(3) = k(2) + dk

      S(2) = S(1) + func(k(2)) * dk
      S(3) = S(2) + func(k(3)) * dk

      i = nkCriti + 2
      do while(iPT < numPTs)
        if(i > nkLimit) then
          if(isFlat) then
            exit
          else
            write(*, '(A)') '>> Number of peaks and troughs NOT enough! '
            res = - huge(res)
            return
          end if
        end if

        if(S(1) == S(2) .and. S(2) == S(3)) then
          ! almost not-oscillatory integrand, a flat peak/trough
          if(.not. isFlat) then
            iPT = 0 ! if nearly stationary, only take the flat parts into account
            isFlat = .true.
          end if
          if(isNotSo) then
            ! for a flat peak/trough, only take one value at the first point
            !           (*) - * - *     *
            ! that is,  /           and  \
            !          *                 (*) - * - *
            iPT = iPT + 1
            valuePTs(iPT) = S(1)
#ifdef DEBUG
            write(*, '(A, I0, 6(2X, G0))') 'iPT = ', iPT, k, S
#endif
          end if
          isNotSo = .false.
        else
          if(isFlat) then
            ! the integrand is not invariable at all time
            isNotSo = .true.
          else if((S(2) - S(1)) * (S(2) - S(3)) > 0.0_MK .or. S(2) == S(3)) then
            ! only if the integrand is absolutely oscillatory,
            !             * - *     *             * - *            *
            ! only take  /      and  \      , not      \  or      /
            !           *             * - *             *    * - *
            iPT = iPT + 1
            valuePTs(iPT) = getExtreme_(k, S, ex = kExtre)
#ifdef DEBUG
            write(*, '(A, I0, 6(2X, G0))') 'iPT = ', iPT, k, S
            write(*, '(A, 3(2X, G0))') '  ', S(2) - S(1), S(2) - S(3), kExtre
#endif
          end if
        end if

        k(1) = k(2); k(2) = k(3)
        S(1) = S(2); S(2) = S(3)

        k(3) = k(3) + dk ! i.e., (i - 1) * dk
        S(3) = S(3) + func(k(3)) * dk
        i = i + 1
      end do

      do i = iPT - 1, 1, - 1
        do j = 1, i, 1
          valuePTs(j) = (valuePTs(j) + valuePTs(j + 1)) / 2.0_MK
        end do
      end do
      res = valuePTs(1)

    end function ptamRun

    real(MK) function getExtreme_(k, S, ex) result(ey)
      real(MK), intent(in) :: k(3), S(3)
      real(MK), intent(out), optional :: ex
      real(MK) :: a(3)
      a(1) = 2.0_MK * S(3) - 4.0_MK * S(2) + 2.0_MK * S(1)
      if(a(1) == 0.0_MK) then
        ! if there is only a tiny difference between S(1) with S(2) when S(2)
        !&  is equal to S(3)
        ex = k(2)
        ey = S(2)
      else
        a(2) = 4.0_MK * S(2) - S(3) - 3.0_MK * S(1)
        a(3) = S(1)
        ex = k(1) - a(2) / (2.0_MK * a(1)) * (k(3) - k(1))
        ey = a(3) - a(2) * a(2) / (4.0_MK * a(1))
      end if
    end function getExtreme_

  end module ptam

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
