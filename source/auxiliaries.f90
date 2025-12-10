!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 29/04/2019
!  For storing variables used by multiple physics schemes.
module auxiliaries
    use types, only: p
    use params

    implicit none

    private
    public precnv, precls, snowcv, snowls, cbmf, tsr, ssrd, ssr, slrd, slr, olr, slru
    public ustr, vstr, shf, evap, hfluxn
    public qdif, denvvs_out, iptop, iptop_out, psa_out, se_out, qa_out, qsat_out
    public mss_out, mse0_out, mse1_out, mss0_out, mss2_out, msthr_out, qthr0out, qthr1out, ktop1out, ktop2out
    public tt_cnv, qt_cnv
    public phis0_out, fmask_out, tsea_out, ssrd_out, slrd_out, shf_out, slru_out, hfluxn_out, tsfc_out, tskin_out, t0_out, stl_am_out, soilw_am_out, tskin_out_2, hfluxn_debug, denvvs_debug, qsat0_debug, dtskin_debug, rlus_debug, shf_debug, evap_debug, t0_debug, dthl_debug, q_debug, soilw_am_debug_final, qsat0_debug_final, q1_debug_final, tskin_debug_final, psa_debug_final, tt_rsw_out

    ! Physical variables shared among all physics schemes
    real(p), dimension(ix,il)   :: precnv !! Convective precipitation  [g/(m^2 s)], total
    real(p), dimension(ix,il)   :: precls !! Large-scale precipitation [g/(m^2 s)], total
    real(p), dimension(ix,il)   :: snowcv !! Convective precipitation  [g/(m^2 s)], snow only
    real(p), dimension(ix,il)   :: snowls !! Large-scale precipitation [g/(m^2 s)], snow only
    real(p), dimension(ix,il)   :: cbmf   !! Cloud-base mass flux
    real(p), dimension(ix,il)   :: tsr    !! Top-of-atmosphere shortwave radiation (downward)
    real(p), dimension(ix,il)   :: ssrd   !! Surface shortwave radiation (downward-only)
    real(p), dimension(ix,il)   :: ssr    !! Surface shortwave radiation (net downward)
    real(p), dimension(ix,il)   :: slrd   !! Surface longwave radiation (downward-only)
    real(p), dimension(ix,il)   :: slr    !! Surface longwave radiation (net upward)
    real(p), dimension(ix,il)   :: olr    !! Outgoing longwave radiation (upward)
    real(p), dimension(ix,il,3) :: slru   !! Surface longwave emission (upward)
    real(p), dimension(ix,il)   :: qdif
    real(p), dimension(ix,il,0:2)   :: denvvs_out
    integer, dimension(ix,il)   :: iptop
    integer, dimension(ix,il)   :: iptop_out
    real(p), dimension(ix,il)   :: psa_out
    real(p), dimension(ix,il,kx) :: se_out
    real(p), dimension(ix,il,kx) :: qa_out
    real(p), dimension(ix,il,kx) :: qsat_out
    real(p), dimension(ix,il,2:kx) :: mss_out
    real(p), dimension(ix,il) :: mse0_out
    real(p), dimension(ix,il) :: mse1_out
    real(p), dimension(ix,il) :: mss0_out
    real(p), dimension(ix,il) :: mss2_out
    real(p), dimension(ix,il) :: msthr_out
    real(p), dimension(ix,il) :: qthr0out
    real(p), dimension(ix,il) :: qthr1out
    integer, dimension(ix,il) :: ktop1out
    integer, dimension(ix,il) :: ktop2out
    real(p), dimension(ix,il,kx) :: tt_cnv
    real(p), dimension(ix,il,kx) :: qt_cnv
    
    real(p), dimension(ix,il) :: phis0_out
    real(p), dimension(ix,il) :: fmask_out
    real(p), dimension(ix,il) :: tsea_out
    real(p), dimension(ix,il) :: ssrd_out
    real(p), dimension(ix,il) :: slrd_out
    real(p), dimension(ix,il,3) :: shf_out
    real(p), dimension(ix,il,3) :: slru_out
    real(p), dimension(ix,il,2) :: hfluxn_out
    real(p), dimension(ix,il) :: tsfc_out
    real(p), dimension(ix,il) :: tskin_out
    real(p), dimension(ix,il) :: t0_out
    real(p), dimension(ix,il) :: stl_am_out
    real(p), dimension(ix,il) :: soilw_am_out

    real(p), dimension(ix,il) :: tskin_out_2
    real(p), dimension(ix,il) :: hfluxn_debug
    real(p), dimension(ix,il) :: denvvs_debug
    real(p), dimension(ix,il) :: qsat0_debug
    real(p), dimension(ix,il) :: dtskin_debug
    real(p), dimension(ix,il) :: rlus_debug
    real(p), dimension(ix,il) :: shf_debug
    real(p), dimension(ix,il) :: evap_debug
    real(p), dimension(ix,il) :: t0_debug
    real(p), dimension(ix,il) :: dthl_debug
    real(p), dimension(ix,il) :: q_debug

    real(p), dimension(ix,il) :: soilw_am_debug_final
    real(p), dimension(ix,il) :: qsat0_debug_final
    real(p), dimension(ix,il) :: q1_debug_final
    real(p), dimension(ix,il) :: tskin_debug_final
    real(p), dimension(ix,il) :: psa_debug_final

    real(p), dimension(ix,il,kx) :: tt_rsw_out

    ! Third dimension -> 1:land, 2:sea, 3: weighted average
    real(p), dimension(ix,il,3) :: ustr   !! U-stress
    real(p), dimension(ix,il,3) :: vstr   !! V-stress
    real(p), dimension(ix,il,3) :: shf    !! Sensible heat flux
    real(p), dimension(ix,il,3) :: evap   !! Evaporation [g/(m^2 s)]
    real(p), dimension(ix,il,2) :: hfluxn !! Net heat flux into surface
end module
