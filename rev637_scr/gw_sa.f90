subroutine gw2sw(j)
  use parm
!  integer                              ::   j              !hru index

!  integer                              ::   isub           !subbasin index


  real                                 ::   qbs_max        !maximum baseflow

!! SWAT variables:
!!    dep_chan(:)      |m              |average daily water depth in channel
!!    sa_q(:)          |mm             |average daily water depth in channel
!  isub = inum1
  qbs_max = ( (shallst(j)-gwqmn(j)) /gw_spyld(j) - dep_chan(inum1)*1000. ) * sa_cond(j)
  shallst(j) = shallst(j) - qbs_max
  sa_q(j) = qbs_max * alpha_bf(j) + sa_q0(j) * (1. - alpha_bf(j))
  sa_q0(j) = sa_q(j)

  if (sa_q(j)>0.) then
    ! if it's positive, it will be added as baseflow
    gw_q(j) = sa_q(j)
  else
    ! if it's negative, it will be accounted for as transmission loss
    gw_q(j) = 0.
  endif

end subroutine


subroutine sw2gw(jrch)
  use parm
  integer                              ::   jrch           !stream reach index

  integer                              ::   j              !hru index
  real                                 :: rttlc1, rttlc2, rttlc0   !

!! SWAT variables:                       ::   isub           !subbasin index
!!    rttlc       |m^3 H2O       |transmission losses from reach on day





!! variables:
  rttlc = 0.
  do j = hru1(isub), hru1(isub) + hrutot(isub) - 1
    if (sa_q(j) < 0.) then
      rttlc = rttlc - sa_q(j) * hru_km(j) * 1.e-3
    end if
  end do

  rttlc2 = rttlc * rchstor(jrch) / (rtwtr + rchstor(jrch))

  if (rchstor(jrch) <= rttlc2) then
    rttlc2 = min(rttlc2, rchstor(jrch))
    rchstor(jrch) = rchstor(jrch) - rttlc2
    rttlc1 = rttlc - rttlc2
    if (rtwtr <= rttlc1) then
      rttlc1 = min(rttlc1, rtwtr)
      rtwtr = rtwtr - rttlc1
    else
      rtwtr = rtwtr - rttlc1
    end if
  else
    rchstor(jrch) = rchstor(jrch) - rttlc2
    rttlc1 = rttlc - rttlc2
    if (rtwtr <= rttlc1) then
      rttlc1 = min(rttlc1, rtwtr)
      rtwtr = rtwtr - rttlc1
    else
      rtwtr = rtwtr - rttlc1
    end if
  end if
  rttlc0 = rttlc1 + rttlc2  !actual loss

  if (rttlc0 < rttlc) then

    ! inflow and storage is smaller than the potential streambed leakage
    ! it will be redistributed to the shallow aquifer
    rttlc0 = rttlc - rttlc0
    do j = hru1(isub), hru1(isub) + hrutot(isub) - 1
      shallst(j) = shallst(j) - rttlc0 * hru_fr(j) / hru_km(j) * 1.e3
    end do
  end if
  rttlc = rttlc1 + rttlc2



end subroutine
