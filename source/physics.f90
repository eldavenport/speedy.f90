module physics
    use types, only: p
    use params

    implicit none

    private
    public initialize_physics, get_physical_tendencies
    public rh_out, tt_lsc, tt_rsw, tt_rlw, tt_pbl_out, tt_pbl, stl_am_a, coa_factor_a, ssrd_a, alb_l_factor_a, psa_a, flux_out_a, qcloud_out, icltop_out, icltop_2, cloudc_out, clstr_out

    real(p), dimension(ix,il,kx) :: rh ! making rh available to input_output
    real(p), dimension(ix,il,kx) :: rh_out
    real(p), dimension(ix,il,kx) :: tt_lsc
    real(p), dimension(ix,il,kx) :: tt_rsw
    real(p), dimension(ix,il,kx) :: tt_rlw
    real(p), dimension(ix,il,kx) :: tt_pbl_out
    real(p), dimension(ix,il,kx) :: tt_pbl
    real(p), dimension(ix,il) :: stl_am_a
    real(p), dimension(il) :: coa_factor_a
    real(p), dimension(ix,il) :: ssrd_a
    real(p), dimension(ix,il) :: alb_l_factor_a
    real(p), dimension(ix,il) :: psa_a
    real(p), dimension(ix,il,kx,2) :: flux_out_a
    real(p), dimension(ix,il) :: qcloud_out
    integer, dimension(ix,il) :: icltop_out
    integer, dimension(ix,il) :: icltop_2
    real(p), dimension(ix,il) :: cloudc_out
    real(p), dimension(ix,il) :: clstr_out

contains
    ! Initialize physical parametrization routines
    subroutine initialize_physics
        use physical_constants, only: grav, cp, p0, sigl, sigh, grdsig, grdscp, wvi
        use geometry, only: hsg, fsg, dhs

        integer :: k

        ! 1.2 Functions of sigma and latitude
        sigh(0) = hsg(1)

        do k = 1, kx
            sigl(k) = log(fsg(k))
            sigh(k) = hsg(k+1)
            grdsig(k) = grav/(dhs(k)*p0)
            grdscp(k) = grdsig(k)/cp
        end do

        ! Weights for vertical interpolation at half-levels(1,kx) and surface
        ! Note that for phys.par. half-lev(k) is between full-lev k and k+1
        ! Fhalf(k) = Ffull(k)+WVI(K,2)*(Ffull(k+1)-Ffull(k))
        ! Fsurf = Ffull(kx)+WVI(kx,2)*(Ffull(kx)-Ffull(kx-1))
        do k = 1, kx-1
            wvi(k,1) = 1./(sigl(k+1)-sigl(k))
            wvi(k,2) = (log(sigh(k))-sigl(k))*wvi(k,1)
        end do

        wvi(kx,1) = 0.
        wvi(kx,2) = (log(0.99)-sigl(kx))*wvi(kx-1,1)
    end

    !> Compute physical parametrization tendencies for u, v, t, q and add them
    !  to the dynamical grid-point tendencies
    subroutine get_physical_tendencies(vor, div, t, q, phi, psl, utend, vtend, ttend, qtend)
        use auxiliaries, only: precnv, precls, cbmf, tsr, ssrd, ssr, slrd, slr, olr, slru, ustr, vstr, shf, evap, hfluxn, qdif, denvvs_out, iptop, iptop_out, psa_out, se_out, qa_out, qsat_out, mss_out, mse0_out, mse1_out, mss0_out, mss2_out, msthr_out, qthr0out, qthr1out, ktop1out, ktop2out, tt_cnv, qt_cnv, phis0_out, fmask_out, tsea_out, ssrd_out, slrd_out, shf_out, slru_out, hfluxn_out, tsfc_out, tskin_out, t0_out, stl_am_out, soilw_am_out, tskin_out_2, hfluxn_debug, denvvs_debug, qsat0_debug, dtskin_debug, rlus_debug, shf_debug, evap_debug, t0_debug, dthl_debug, q_debug, soilw_am_debug_final, qsat0_debug_final, q1_debug_final, tskin_debug_final, psa_debug_final
        use physical_constants, only: sigh, grdsig, grdscp, cp
        use geometry, only: fsg
        use boundaries, only: phis0
        use land_model, only: fmask_l
        use sea_model, only: sst_am, ssti_om, sea_coupling_flag
        use sppt, only: mu, gen_sppt
        use convection, only: get_convection_tendencies
        use large_scale_condensation, only: get_large_scale_condensation_tendencies
        use shortwave_radiation, only: get_shortwave_rad_fluxes, clouds, compute_shortwave
        use longwave_radiation, only: &
                get_downward_longwave_rad_fluxes, get_upward_longwave_rad_fluxes
        use surface_fluxes, only: get_surface_fluxes
        use vertical_diffusion, only: get_vertical_diffusion_tend
        use humidity, only: spec_hum_to_rel_hum
        use spectral, only: spec_to_grid, uvspec

        complex(p), intent(in) :: vor(mx,nx,kx) !! Vorticity
        complex(p), intent(in) :: div(mx,nx,kx) !! Divergence
        complex(p), intent(in) :: t(mx,nx,kx)   !! Temperature
        complex(p), intent(in) :: q(mx,nx,kx)   !! Specific Humidity
        complex(p), intent(in) :: phi(mx,nx,kx) !! Geopotential
        complex(p), intent(in) :: psl(mx,nx)    !! ln(Surface pressure)

        real(p), intent(inout) :: utend(ix,il,kx) !! Zonal velocity tendency
        real(p), intent(inout) :: vtend(ix,il,kx) !! Meridional velocity tendency
        real(p), intent(inout) :: ttend(ix,il,kx) !! Temperature tendency
        real(p), intent(inout) :: qtend(ix,il,kx) !! Specific humidity tendency

        complex(p), dimension(mx,nx) :: ucos, vcos
        real(p), dimension(ix,il) :: pslg, rps, gse
        real(p), dimension(ix,il,kx) :: ug, vg, tg, qg, phig, utend_dyn, vtend_dyn, ttend_dyn, qtend_dyn
        real(p), dimension(ix,il,kx) :: se, qsat
        real(p), dimension(ix,il) :: psg, ts, tskin, u0, v0, t0, cloudc, clstr, cltop, prtop
        real(p), dimension(ix,il,kx) :: qt_lsc, ut_pbl, vt_pbl,&
            & qt_pbl
        integer :: icltop(ix,il,2), icnv(ix,il), i, j, k
        real(p) :: sppt_pattern(ix,il,kx)

        ! Keep a copy of the original (dynamics only) tendencies
        utend_dyn = utend
        vtend_dyn = vtend
        ttend_dyn = ttend
        qtend_dyn = qtend

        ! =========================================================================
        ! Compute grid-point fields
        ! =========================================================================

        ! Convert model spectral variables to grid-point variables
        do k = 1, kx
            call uvspec(vor(:,:,k), div(:,:,k), ucos, vcos)
            ug(:,:,k)   = spec_to_grid(ucos, 2)
            vg(:,:,k)   = spec_to_grid(vcos, 2)
            tg(:,:,k)   = spec_to_grid(t(:,:,k), 1)
            qg(:,:,k)   = spec_to_grid(q(:,:,k), 1)
            phig(:,:,k) = spec_to_grid(phi(:,:,k), 1)
        end do

        pslg = spec_to_grid(psl, 1)

        ! =========================================================================
        ! Compute thermodynamic variables
        ! =========================================================================

        psg = exp(pslg)
        rps = 1.0/psg

        qg = max(qg, 0.0)
        se = cp*tg + phig

        do k = 1, kx
            call spec_hum_to_rel_hum(tg(:,:,k), psg, fsg(k), qg(:,:,k), rh(:,:,k), qsat(:,:,k))
        end do

        ! =========================================================================
        ! Precipitation
        ! =========================================================================

        ! Deep convection
        call get_convection_tendencies(psg, se, qg, qsat, iptop, cbmf, precnv, tt_cnv, qt_cnv, qdif, psa_out, se_out, qa_out, qsat_out, mss_out, mse0_out, mse1_out, mss0_out, mss2_out, msthr_out, qthr0out, qthr1out, ktop1out, ktop2out)

        iptop_out = iptop

        do k = 2, kx
            tt_cnv(:,:,k) = tt_cnv(:,:,k)*rps*grdscp(k)
            qt_cnv(:,:,k) = qt_cnv(:,:,k)*rps*grdsig(k)
        end do

        icnv = kx - iptop

        ! Large-scale condensation
        call get_large_scale_condensation_tendencies(psg, qg, qsat, iptop, precls, tt_lsc, qt_lsc)

        ttend = ttend + tt_cnv + tt_lsc
        qtend = qtend + qt_cnv + qt_lsc

        ! =========================================================================
        ! Radiation (shortwave and longwave) and surface fluxes
        ! =========================================================================

        ! Compute shortwave tendencies and initialize lw transmissivity
        ! The shortwave radiation may be called at selected time steps
        fmask_out = 0.0*fmask_l
        if (compute_shortwave) then
            gse = (se(:,:,kx-1) - se(:,:,kx))/(phig(:,:,kx-1) - phig(:,:,kx))

            call clouds(qg, rh, precnv, precls, iptop, gse, fmask_l, icltop, icltop_2, cloudc, clstr, qcloud_out)

            do i = 1, ix
                do j = 1, il
                    cltop(i,j) = sigh(icltop(i,j,1) - 1)*psg(i,j)
                    prtop(i,j) = float(iptop(i,j))
                end do
            end do

            call get_shortwave_rad_fluxes(psg, qg, icltop, cloudc, clstr, ssrd, ssr, tsr, tt_rsw, flux_out_a)

            do k = 1, kx
                tt_rsw(:,:,k) = tt_rsw(:,:,k)*rps*grdscp(k)
            end do
            fmask_out = fmask_l
        end if

        icltop_out = icltop(:,:,1)
        cloudc_out = cloudc
        clstr_out = clstr

        ! Compute downward longwave fluxes
        call get_downward_longwave_rad_fluxes(tg, slrd, tt_rlw)

        rh_out = rh
        phis0_out = phis0 ! FIXME
        tsea_out = sst_am
        ssrd_out = ssrd
        slrd_out = slrd

        ! Compute surface fluxes and land skin temperature
        call get_surface_fluxes(psg, ug, vg, tg, qg, rh, phig, phis0, fmask_l, sst_am, ssrd, slrd, ustr, vstr, shf, evap, slru, hfluxn, ts, tskin, u0, v0, t0, .true., denvvs_out, stl_am_out, soilw_am_out, tskin_out_2, hfluxn_debug, denvvs_debug, qsat0_debug, dtskin_debug, rlus_debug, shf_debug, evap_debug, t0_debug, dthl_debug, q_debug, soilw_am_debug_final, qsat0_debug_final, q1_debug_final, tskin_debug_final, psa_debug_final, stl_am_a, coa_factor_a, ssrd_a, alb_l_factor_a, psa_a)

        shf_out = shf
        slru_out = slru
        hfluxn_out = hfluxn
        tsfc_out = ts
        tskin_out = tskin
        t0_out = t0

        sea_coupling_flag = 0

        ! Recompute sea fluxes in case of anomaly coupling
        if (sea_coupling_flag > 0) then
           call get_surface_fluxes(psg, ug, vg, tg, qg, rh, phig, phis0, fmask_l, ssti_om, ssrd, slrd, ustr, vstr, shf, evap, slru, hfluxn, ts, tskin, u0, v0, t0, .false., denvvs_out, stl_am_out, soilw_am_out, tskin_out_2, hfluxn_debug, denvvs_debug, qsat0_debug, dtskin_debug, rlus_debug, shf_debug, evap_debug, t0_debug, dthl_debug, q_debug, soilw_am_debug_final, qsat0_debug_final, q1_debug_final, tskin_debug_final, psa_debug_final, stl_am_a, coa_factor_a, ssrd_a, alb_l_factor_a, psa_a)
        end if

        ! Compute upward longwave fluxes, convert them to tendencies and add
        ! shortwave tendencies
        call get_upward_longwave_rad_fluxes(tg, ts, slrd, slru(:,:,3), slr, olr, tt_rlw)

        do k = 1, kx
            tt_rlw(:,:,k) = tt_rlw(:,:,k)*rps*grdscp(k)
        end do

        ttend = ttend + tt_rsw + tt_rlw

        ! =========================================================================
        ! Planetary boundary later interactions with lower troposphere
        ! =========================================================================

        ! Vertical diffusion and shallow convection
        call get_vertical_diffusion_tend(se, rh, qg, qsat, phig, icnv, ut_pbl, vt_pbl, &
            & tt_pbl, qt_pbl)

        tt_pbl_out = tt_pbl

        ! Add tendencies due to surface fluxes
        ut_pbl(:,:,kx) = ut_pbl(:,:,kx) + ustr(:,:,3)*rps*grdsig(kx)
        vt_pbl(:,:,kx) = vt_pbl(:,:,kx) + vstr(:,:,3)*rps*grdsig(kx)
        tt_pbl(:,:,kx) = tt_pbl(:,:,kx)  + shf(:,:,3)*rps*grdscp(kx)
        qt_pbl(:,:,kx) = qt_pbl(:,:,kx) + evap(:,:,3)*rps*grdsig(kx)

        utend = utend + ut_pbl
        vtend = vtend + vt_pbl
        ttend = ttend + tt_pbl
        qtend = qtend + qt_pbl

        ! Add SPPT noise
        if (sppt_on) then
            sppt_pattern = gen_sppt()

            ! The physical contribution to the tendency is *tend - *tend_dyn, where * is u, v, t, q
            do k = 1,kx
                utend(:,:,k) = (1 + sppt_pattern(:,:,k)*mu(k))*(utend(:,:,k) - utend_dyn(:,:,k)) &
                        & + utend_dyn(:,:,k)
                vtend(:,:,k) = (1 + sppt_pattern(:,:,k)*mu(k))*(vtend(:,:,k) - vtend_dyn(:,:,k)) &
                        & + vtend_dyn(:,:,k)
                ttend(:,:,k) = (1 + sppt_pattern(:,:,k)*mu(k))*(ttend(:,:,k) - ttend_dyn(:,:,k)) &
                        & + ttend_dyn(:,:,k)
                qtend(:,:,k) = (1 + sppt_pattern(:,:,k)*mu(k))*(qtend(:,:,k) - qtend_dyn(:,:,k)) &
                        & + qtend_dyn(:,:,k)
            end do
        end if
    end
end module
