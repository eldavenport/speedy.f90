!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 08/05/2019
!  For performing input and output.
module input_output
    use types, only: p, sp
    use netcdf
    use params
    use physics, only: &
        ! SWRadiationData
        qcloud_out, fsol_out, rsds_out, rsns_out, ozone_out, ozupp_out, zenit_out, stratz_out, gse_out, &
        icltop_out, cloudc_out, cloudstr_out, ftop_sw_out, dfabs_sw_out, &
        ! LWRadiationData
        dfabs_lw_out, ftop_lw_out, &
        ! ConvectionData
        se_out, iptop_out, cbmf_out, qdif_out, precnv_out, &
        ! ModRadConData
        ablco2_out, alb_l_out, alb_s_out, albsfc_out, snowc_out, &
        tau2_out, st4a_out, stratc_out, flux_out, &
        ! HumidityData
        rh_out, qsat_out, &
        ! CondensationData
        precls_out, dtlsc_out, dqlsc_out, &
        ! SurfaceFluxData
        ustr_out, vstr_out, shf_out, evap_out, rlus_out, rlds_out, rlns_out, hfluxn_out, &
        tsfc_out, tskin_out, u0_out, v0_out, t0_out

    implicit none

    private
    public output

contains
    !> Writes a snapshot of all prognostic (and other) variables to a NetCDF file.
    subroutine output(timestep, vor, div, t, ps, tr, phi)
        use geometry, only: radang, fsg
        use physical_constants, only: p0, grav
        use date, only: model_datetime, start_datetime
        use spectral, only: spec_to_grid, uvspec

        integer, intent(in) :: timestep           !! The time step that is being written
        complex(p), intent(in) :: vor(mx,nx,kx,2)    !! Vorticity
        complex(p), intent(in) :: div(mx,nx,kx,2)    !! Divergence
        complex(p), intent(in) :: t(mx,nx,kx,2)      !! Temperature
        complex(p), intent(in) :: ps(mx,nx,2)        !! log(normalized surface pressure)
        complex(p), intent(in) :: tr(mx,nx,kx,2,ntr) !! Tracers
        complex(p), intent(in) :: phi(mx,nx,kx)      !! Geopotential

        complex(p), dimension(mx,nx)     :: ucos, vcos
        real(p), dimension(ix,il,kx)  :: u_grid, v_grid, t_grid, q_grid, phi_grid
        real(p), dimension(ix,il)     :: ps_grid
        real(sp), dimension(ix,il,kx) :: u_out, v_out, t_out, q_out, phi_out
        real(sp), dimension(ix,il)    :: ps_out
        character(len=15) :: file_name = 'yyyymmddhhmm.nc'
        character(len=32) :: time_template = 'hours since yyyy-mm-dd hh:mm:0.0'
        integer :: k, ncid
        integer :: timedim, latdim, londim, levdim
        integer :: timevar, latvar, lonvar, levvar, uvar, vvar, tvar, qvar, phivar, psvar
        ! SWRadiationData vars
        integer :: qcloud_var, fsol_var, rsds_var, rsns_var, ozone_var, ozupp_var, zenit_var, stratz_var, gse_var
        integer :: icltop_var, cloudc_var, cloudstr_var, ftop_sw_var, dfabs_sw_var
        ! LWRadiationData vars
        integer :: dfabs_lw_var, ftop_lw_var
        ! ConvectionData vars
        integer :: se_var, iptop_var, cbmf_var, qdif_var, precnv_var
        ! ModRadConData vars
        integer :: ablco2_var, alb_l_var, alb_s_var, albsfc_var, snowc_var
        integer :: tau2_var_1, tau2_var_2, tau2_var_3, tau2_var_4
        integer :: st4a_var_1, st4a_var_2
        integer :: stratc_var_1, stratc_var_2
        integer :: flux_var_1, flux_var_2, flux_var_3, flux_var_4
        ! HumidityData vars
        integer :: rh_var, qsat_var
        ! CondensationData vars
        integer :: precls_var, dtlsc_var, dqlsc_var
        ! SurfaceFluxData vars
        integer :: ustr_var_1, ustr_var_2, ustr_var_3
        integer :: vstr_var_1, vstr_var_2, vstr_var_3
        integer :: shf_var_1, shf_var_2, shf_var_3
        integer :: evap_var_1, evap_var_2, evap_var_3
        integer :: rlus_var_1, rlus_var_2, rlus_var_3
        integer :: rlds_var, rlns_var
        integer :: hfluxn_var_1, hfluxn_var_2
        integer :: tsfc_var, tskin_var, u0_var, v0_var, t0_var

        ! Construct file_name
        write (file_name(1:4),'(i4.4)') model_datetime%year
        write (file_name(5:6),'(i2.2)') model_datetime%month
        write (file_name(7:8),'(i2.2)') model_datetime%day
        write (file_name(9:10),'(i2.2)') model_datetime%hour
        write (file_name(11:12),'(i2.2)') model_datetime%minute

        ! Construct time string
        write (time_template(13:16),'(i4.4)') start_datetime%year
        write (time_template(18:19),'(i2.2)') start_datetime%month
        write (time_template(21:22),'(i2.2)') start_datetime%day
        write (time_template(24:25),'(i2.2)') start_datetime%hour
        write (time_template(27:28),'(i2.2)') start_datetime%minute

        ! Create NetCDF output file
        call check(nf90_create(file_name, nf90_clobber, ncid))

        ! Define time
        call check(nf90_def_dim(ncid, "time", nf90_unlimited, timedim))
        call check(nf90_def_var(ncid, "time", nf90_real4, timedim, timevar))
        call check(nf90_put_att(ncid, timevar, "units", time_template))

        ! Define space
        call check(nf90_def_dim(ncid, "lon", ix, londim))
        call check(nf90_def_dim(ncid, "lat", il, latdim))
        call check(nf90_def_dim(ncid, "lev", kx, levdim))
        call check(nf90_def_var(ncid, "lon", nf90_real4, londim, lonvar))
        call check(nf90_put_att(ncid, lonvar, "long_name", "longitude"))
        call check(nf90_def_var(ncid, "lat", nf90_real4, latdim, latvar))
        call check(nf90_put_att(ncid, latvar, "long_name", "latitude"))
        call check(nf90_def_var(ncid, "lev", nf90_real4, levdim, levvar))
        call check(nf90_put_att(ncid, levvar, "long_name", "atmosphere_sigma_coordinate"))

        ! Define prognostic fields
        call check(nf90_def_var(ncid, "u", nf90_real4, (/ londim, latdim, levdim, timedim /), uvar))
        call check(nf90_put_att(ncid, uvar, "long_name", "eastward_wind"))
        call check(nf90_put_att(ncid, uvar, "units", "m/s"))
        call check(nf90_def_var(ncid, "v", nf90_real4, (/ londim, latdim, levdim, timedim /), vvar))
        call check(nf90_put_att(ncid, vvar, "long_name", "northward_wind"))
        call check(nf90_put_att(ncid, vvar, "units", "m/s"))
        call check(nf90_def_var(ncid, "t", nf90_real4, (/ londim, latdim, levdim, timedim /), tvar))
        call check(nf90_put_att(ncid, tvar, "long_name", "air_temperature"))
        call check(nf90_put_att(ncid, tvar, "units", "K"))

        call check(nf90_def_var(ncid, "q", nf90_real4, (/ londim, latdim, levdim, timedim /), qvar))
        call check(nf90_put_att(ncid, qvar, "long_name", "specific_humidity"))
        call check(nf90_put_att(ncid, qvar, "units", "1"))

        call check(nf90_def_var(ncid, "phi", nf90_real4, (/ londim, latdim, levdim, timedim /), &
            & phivar))
        call check(nf90_put_att(ncid, phivar, "long_name", "geopotential_height"))
        call check(nf90_put_att(ncid, phivar, "units", "m"))

        call check(nf90_def_var(ncid, "ps", nf90_real4, (/ londim, latdim, timedim /), psvar))
        call check(nf90_put_att(ncid, psvar, "long_name", "surface_air_pressure"))
        call check(nf90_put_att(ncid, psvar, "units", "Pa"))

        ! Define diagnostic fields
        ! SWRadiationData - 2D fields
        call check(nf90_def_var(ncid, "qcloud", nf90_real4, (/ londim, latdim, timedim /), qcloud_var))
        call check(nf90_def_var(ncid, "fsol", nf90_real4, (/ londim, latdim, timedim /), fsol_var))
        call check(nf90_def_var(ncid, "rsds", nf90_real4, (/ londim, latdim, timedim /), rsds_var))
        call check(nf90_def_var(ncid, "rsns", nf90_real4, (/ londim, latdim, timedim /), rsns_var))
        call check(nf90_def_var(ncid, "ozone", nf90_real4, (/ londim, latdim, timedim /), ozone_var))
        call check(nf90_def_var(ncid, "ozupp", nf90_real4, (/ londim, latdim, timedim /), ozupp_var))
        call check(nf90_def_var(ncid, "zenit", nf90_real4, (/ londim, latdim, timedim /), zenit_var))
        call check(nf90_def_var(ncid, "stratz", nf90_real4, (/ londim, latdim, timedim /), stratz_var))
        call check(nf90_def_var(ncid, "gse", nf90_real4, (/ londim, latdim, timedim /), gse_var))
        call check(nf90_def_var(ncid, "icltop", nf90_int, (/ londim, latdim, timedim /), icltop_var))
        call check(nf90_def_var(ncid, "cloudc", nf90_real4, (/ londim, latdim, timedim /), cloudc_var))
        call check(nf90_def_var(ncid, "cloudstr", nf90_real4, (/ londim, latdim, timedim /), cloudstr_var))
        call check(nf90_def_var(ncid, "ftop_sw", nf90_real4, (/ londim, latdim, timedim /), ftop_sw_var))
        ! SWRadiationData - 3D fields
        call check(nf90_def_var(ncid, "dfabs_sw", nf90_real4, (/ londim, latdim, levdim, timedim /), dfabs_sw_var))

        ! LWRadiationData
        call check(nf90_def_var(ncid, "dfabs_lw", nf90_real4, (/ londim, latdim, levdim, timedim /), dfabs_lw_var))
        call check(nf90_def_var(ncid, "ftop_lw", nf90_real4, (/ londim, latdim, timedim /), ftop_lw_var))

        ! ConvectionData
        call check(nf90_def_var(ncid, "se", nf90_real4, (/ londim, latdim, levdim, timedim /), se_var))
        call check(nf90_def_var(ncid, "iptop", nf90_int, (/ londim, latdim, timedim /), iptop_var))
        call check(nf90_def_var(ncid, "cbmf", nf90_real4, (/ londim, latdim, timedim /), cbmf_var))
        call check(nf90_def_var(ncid, "qdif", nf90_real4, (/ londim, latdim, timedim /), qdif_var))
        call check(nf90_def_var(ncid, "precnv", nf90_real4, (/ londim, latdim, timedim /), precnv_var))

        ! ModRadConData - 2D fields
        call check(nf90_def_var(ncid, "alb_l", nf90_real4, (/ londim, latdim, timedim /), alb_l_var))
        call check(nf90_def_var(ncid, "alb_s", nf90_real4, (/ londim, latdim, timedim /), alb_s_var))
        call check(nf90_def_var(ncid, "albsfc", nf90_real4, (/ londim, latdim, timedim /), albsfc_var))
        call check(nf90_def_var(ncid, "snowc", nf90_real4, (/ londim, latdim, timedim /), snowc_var))
        ! ModRadConData - 3D fields with extra dimensions (split by channel)
        call check(nf90_def_var(ncid, "tau2_1", nf90_real4, (/ londim, latdim, levdim, timedim /), tau2_var_1))
        call check(nf90_def_var(ncid, "tau2_2", nf90_real4, (/ londim, latdim, levdim, timedim /), tau2_var_2))
        call check(nf90_def_var(ncid, "tau2_3", nf90_real4, (/ londim, latdim, levdim, timedim /), tau2_var_3))
        call check(nf90_def_var(ncid, "tau2_4", nf90_real4, (/ londim, latdim, levdim, timedim /), tau2_var_4))
        call check(nf90_def_var(ncid, "st4a_1", nf90_real4, (/ londim, latdim, levdim, timedim /), st4a_var_1))
        call check(nf90_def_var(ncid, "st4a_2", nf90_real4, (/ londim, latdim, levdim, timedim /), st4a_var_2))
        ! ModRadConData - 2D fields with extra dimensions
        call check(nf90_def_var(ncid, "stratc_1", nf90_real4, (/ londim, latdim, timedim /), stratc_var_1))
        call check(nf90_def_var(ncid, "stratc_2", nf90_real4, (/ londim, latdim, timedim /), stratc_var_2))
        call check(nf90_def_var(ncid, "flux_1", nf90_real4, (/ londim, latdim, timedim /), flux_var_1))
        call check(nf90_def_var(ncid, "flux_2", nf90_real4, (/ londim, latdim, timedim /), flux_var_2))
        call check(nf90_def_var(ncid, "flux_3", nf90_real4, (/ londim, latdim, timedim /), flux_var_3))
        call check(nf90_def_var(ncid, "flux_4", nf90_real4, (/ londim, latdim, timedim /), flux_var_4))

        ! HumidityData
        call check(nf90_def_var(ncid, "rh", nf90_real4, (/ londim, latdim, levdim, timedim /), rh_var))
        call check(nf90_def_var(ncid, "qsat", nf90_real4, (/ londim, latdim, levdim, timedim /), qsat_var))

        ! CondensationData
        call check(nf90_def_var(ncid, "precls", nf90_real4, (/ londim, latdim, timedim /), precls_var))
        call check(nf90_def_var(ncid, "dtlsc", nf90_real4, (/ londim, latdim, levdim, timedim /), dtlsc_var))
        call check(nf90_def_var(ncid, "dqlsc", nf90_real4, (/ londim, latdim, levdim, timedim /), dqlsc_var))

        ! SurfaceFluxData - 2D fields with channels (split)
        call check(nf90_def_var(ncid, "ustr_1", nf90_real4, (/ londim, latdim, timedim /), ustr_var_1))
        call check(nf90_def_var(ncid, "ustr_2", nf90_real4, (/ londim, latdim, timedim /), ustr_var_2))
        call check(nf90_def_var(ncid, "ustr_3", nf90_real4, (/ londim, latdim, timedim /), ustr_var_3))
        call check(nf90_def_var(ncid, "vstr_1", nf90_real4, (/ londim, latdim, timedim /), vstr_var_1))
        call check(nf90_def_var(ncid, "vstr_2", nf90_real4, (/ londim, latdim, timedim /), vstr_var_2))
        call check(nf90_def_var(ncid, "vstr_3", nf90_real4, (/ londim, latdim, timedim /), vstr_var_3))
        call check(nf90_def_var(ncid, "shf_1", nf90_real4, (/ londim, latdim, timedim /), shf_var_1))
        call check(nf90_def_var(ncid, "shf_2", nf90_real4, (/ londim, latdim, timedim /), shf_var_2))
        call check(nf90_def_var(ncid, "shf_3", nf90_real4, (/ londim, latdim, timedim /), shf_var_3))
        call check(nf90_def_var(ncid, "evap_1", nf90_real4, (/ londim, latdim, timedim /), evap_var_1))
        call check(nf90_def_var(ncid, "evap_2", nf90_real4, (/ londim, latdim, timedim /), evap_var_2))
        call check(nf90_def_var(ncid, "evap_3", nf90_real4, (/ londim, latdim, timedim /), evap_var_3))
        call check(nf90_def_var(ncid, "rlus_1", nf90_real4, (/ londim, latdim, timedim /), rlus_var_1))
        call check(nf90_def_var(ncid, "rlus_2", nf90_real4, (/ londim, latdim, timedim /), rlus_var_2))
        call check(nf90_def_var(ncid, "rlus_3", nf90_real4, (/ londim, latdim, timedim /), rlus_var_3))
        call check(nf90_def_var(ncid, "rlds", nf90_real4, (/ londim, latdim, timedim /), rlds_var))
        call check(nf90_def_var(ncid, "rlns", nf90_real4, (/ londim, latdim, timedim /), rlns_var))
        call check(nf90_def_var(ncid, "hfluxn_1", nf90_real4, (/ londim, latdim, timedim /), hfluxn_var_1))
        call check(nf90_def_var(ncid, "hfluxn_2", nf90_real4, (/ londim, latdim, timedim /), hfluxn_var_2))
        ! SurfaceFluxData - simple 2D fields
        call check(nf90_def_var(ncid, "tsfc", nf90_real4, (/ londim, latdim, timedim /), tsfc_var))
        call check(nf90_def_var(ncid, "tskin", nf90_real4, (/ londim, latdim, timedim /), tskin_var))
        call check(nf90_def_var(ncid, "u0", nf90_real4, (/ londim, latdim, timedim /), u0_var))
        call check(nf90_def_var(ncid, "v0", nf90_real4, (/ londim, latdim, timedim /), v0_var))
        call check(nf90_def_var(ncid, "t0", nf90_real4, (/ londim, latdim, timedim /), t0_var))

        call check(nf90_enddef(ncid))

        ! Write dimensions to file
        call check(nf90_put_var(ncid, timevar, timestep*24.0/real(nsteps,sp),               (/ 1 /)))
        call check(nf90_put_var(ncid, lonvar, (/ (3.75*k, k = 0, ix-1) /),                 (/ 1 /)))
        call check(nf90_put_var(ncid, latvar, (/ (radang(k)*90.0/asin(1.0), k = 1, il) /), (/ 1 /)))
        call check(nf90_put_var(ncid, levvar, (/ (fsg(k), k = 1, 8) /),                    (/ 1 /)))

        ! Convert prognostic fields from spectral space to grid point space
        do k = 1, kx
           call uvspec(vor(:,:,k,1), div(:,:,k,1), ucos, vcos)
           u_grid(:,:,k)   = spec_to_grid(ucos, 2)
           v_grid(:,:,k)   = spec_to_grid(vcos, 2)
           t_grid(:,:,k)   = spec_to_grid(t(:,:,k,1), 1)
           q_grid(:,:,k)   = spec_to_grid(tr(:,:,k,1,1), 1)
           phi_grid(:,:,k) = spec_to_grid(phi(:,:,k), 1)
        end do
        ps_grid = spec_to_grid(ps(:,:,1), 1)

        ! Output date
        print '(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)',&
            & 'Write gridded dataset for year/month/date/hour/minute: ', &
            & model_datetime%year,'/',model_datetime%month,'/',model_datetime%day,'/', &
            & model_datetime%hour,'/',model_datetime%minute

        ! Preprocess output variables
        u_out = real(u_grid, sp)
        v_out = real(v_grid, sp)
        t_out = real(t_grid, sp)
        q_out = real(q_grid*1.0e-3, sp) ! kg/kg
        phi_out = real(phi_grid/grav, sp)   ! m
        ps_out = real(p0*exp(ps_grid), sp)! Pa

        ! Write prognostic variables to file
        call check(nf90_put_var(ncid, uvar, u_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, vvar, v_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tvar, t_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qvar, q_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, phivar, phi_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, psvar, ps_out, (/ 1, 1, 1 /)))

        ! Write diagnostic variables to file
        ! SWRadiationData - 2D fields
        call check(nf90_put_var(ncid, qcloud_var, real(qcloud_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, fsol_var, real(fsol_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, rsds_var, real(rsds_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, rsns_var, real(rsns_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, ozone_var, real(ozone_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, ozupp_var, real(ozupp_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, zenit_var, real(zenit_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, stratz_var, real(stratz_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, gse_var, real(gse_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, icltop_var, icltop_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, cloudc_var, real(cloudc_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, cloudstr_var, real(cloudstr_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, ftop_sw_var, real(ftop_sw_out, sp), (/ 1, 1, 1 /)))
        ! SWRadiationData - 3D fields
        call check(nf90_put_var(ncid, dfabs_sw_var, real(dfabs_sw_out, sp), (/ 1, 1, 1, 1 /)))

        ! LWRadiationData
        call check(nf90_put_var(ncid, dfabs_lw_var, real(dfabs_lw_out, sp), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, ftop_lw_var, real(ftop_lw_out, sp), (/ 1, 1, 1 /)))

        ! ConvectionData
        call check(nf90_put_var(ncid, se_var, real(se_out, sp), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, iptop_var, iptop_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, cbmf_var, real(cbmf_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qdif_var, real(qdif_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, precnv_var, real(precnv_out, sp), (/ 1, 1, 1 /)))

        ! ModRadConData - 2D fields
        call check(nf90_put_var(ncid, alb_l_var, real(alb_l_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, alb_s_var, real(alb_s_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, albsfc_var, real(albsfc_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, snowc_var, real(snowc_out, sp), (/ 1, 1, 1 /)))
        ! ModRadConData - 3D fields with extra dimensions (split by channel)
        call check(nf90_put_var(ncid, tau2_var_1, real(tau2_out(:,:,:,1), sp), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tau2_var_2, real(tau2_out(:,:,:,2), sp), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tau2_var_3, real(tau2_out(:,:,:,3), sp), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tau2_var_4, real(tau2_out(:,:,:,4), sp), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, st4a_var_1, real(st4a_out(:,:,:,1), sp), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, st4a_var_2, real(st4a_out(:,:,:,2), sp), (/ 1, 1, 1, 1 /)))
        ! ModRadConData - 2D fields with extra dimensions
        call check(nf90_put_var(ncid, stratc_var_1, real(stratc_out(:,:,1), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, stratc_var_2, real(stratc_out(:,:,2), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, flux_var_1, real(flux_out(:,:,1), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, flux_var_2, real(flux_out(:,:,2), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, flux_var_3, real(flux_out(:,:,3), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, flux_var_4, real(flux_out(:,:,4), sp), (/ 1, 1, 1 /)))

        ! HumidityData
        call check(nf90_put_var(ncid, rh_var, real(rh_out, sp), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qsat_var, real(qsat_out, sp), (/ 1, 1, 1, 1 /)))

        ! CondensationData
        call check(nf90_put_var(ncid, precls_var, real(precls_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, dtlsc_var, real(dtlsc_out, sp), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, dqlsc_var, real(dqlsc_out, sp), (/ 1, 1, 1, 1 /)))

        ! SurfaceFluxData - 2D fields with channels (split)
        call check(nf90_put_var(ncid, ustr_var_1, real(ustr_out(:,:,1), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, ustr_var_2, real(ustr_out(:,:,2), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, ustr_var_3, real(ustr_out(:,:,3), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, vstr_var_1, real(vstr_out(:,:,1), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, vstr_var_2, real(vstr_out(:,:,2), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, vstr_var_3, real(vstr_out(:,:,3), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, shf_var_1, real(shf_out(:,:,1), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, shf_var_2, real(shf_out(:,:,2), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, shf_var_3, real(shf_out(:,:,3), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, evap_var_1, real(evap_out(:,:,1), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, evap_var_2, real(evap_out(:,:,2), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, evap_var_3, real(evap_out(:,:,3), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, rlus_var_1, real(rlus_out(:,:,1), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, rlus_var_2, real(rlus_out(:,:,2), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, rlus_var_3, real(rlus_out(:,:,3), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, rlds_var, real(rlds_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, rlns_var, real(rlns_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, hfluxn_var_1, real(hfluxn_out(:,:,1), sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, hfluxn_var_2, real(hfluxn_out(:,:,2), sp), (/ 1, 1, 1 /)))
        ! SurfaceFluxData - simple 2D fields
        call check(nf90_put_var(ncid, tsfc_var, real(tsfc_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tskin_var, real(tskin_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, u0_var, real(u0_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, v0_var, real(v0_out, sp), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, t0_var, real(t0_out, sp), (/ 1, 1, 1 /)))

        call check(nf90_close(ncid))
    end subroutine

    !> Handles any errors from the NetCDF API.
    subroutine check(ierr)
        integer, intent(in) :: ierr

        if (ierr /= nf90_noerr) then
            print *, trim(adjustl(nf90_strerror(ierr)))
            stop
        end if
    end subroutine
end module
