!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 08/05/2019
!  For performing input and output.
module input_output
    use types, only: p, sp
    use netcdf
    use params
    use physics, only: rh_out, tt_lsc, tt_rsw, tt_rlw, tt_pbl_out, tt_pbl, stl_am_a, coa_factor_a, ssrd_a, alb_l_factor_a, psa_a, flux_out_a, qcloud_out, icltop_out, icltop_2, cloudc_out, clstr_out
    use auxiliaries, only: precnv, qdif, evap, denvvs_out, iptop_out, psa_out, se_out, qa_out, qsat_out, mss_out, mse0_out, mse1_out, mss0_out, mss2_out, msthr_out, qthr0out, qthr1out, ktop1out, ktop2out, tt_cnv, qt_cnv, phis0_out, fmask_out, tsea_out, ssrd_out, slrd_out, shf_out, slru_out, hfluxn_out, tsfc_out, tskin_out, t0_out, stl_am_out, soilw_am_out, tskin_out_2, hfluxn_debug, denvvs_debug, qsat0_debug, dtskin_debug, rlus_debug, shf_debug, evap_debug, t0_debug, dthl_debug, q_debug, soilw_am_debug_final, qsat0_debug_final, q1_debug_final, tskin_debug_final, psa_debug_final

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
        integer :: timevar, latvar, lonvar, levvar, uvar, vvar, tvar, qvar, phivar, psvar, rhvar, precnvvar, qdifvar, evapvar_0, evapvar_1, evapvar_2, denvvsvar_0, denvvsvar_1, denvvsvar_2, iptopvar, psa_outvar, se_outvar, qa_outvar, qsat_outvar, mss_outvar, mse0_outvar, mse1_outvar, mss0_outvar, mss2_outvar, msthr_outvar, qthr0outvar, qthr1outvar, ktop1outvar, ktop2outvar, tt_cnvvar, qt_cnvvar, phis0var, fmaskvar, tseavar, ssrdvar, slrdvar, shfvar0, shfvar1, shfvar2, slruvar0, slruvar1, slruvar2, hfluxnvar0, hfluxnvar1, tsfcvar, tskinvar, t0var, stl_amvar, soilw_amvar, tskin_out_2var, hfluxn_debugvar, denvvs_debugvar, qsat0_debugvar, dtskin_debugvar, rlus_debugvar, shf_debugvar, evap_debugvar, t0_debugvar, dthl_debugvar, q_debugvar, soilw_am_debug_finalvar, qsat0_debug_finalvar, q1_debug_finalvar, tskin_debug_finalvar, psa_debug_finalvar, tt_lscvar, tt_rswvar, tt_rlwvar, tt_pbl_outvar, tt_pblvar, stl_am_avar, coa_factor_avar, ssrd_avar, alb_l_factor_avar, psa_avar, flux_out_a0var, flux_out_a1var, qcloud_outvar, icltop_outvar, icltop_2var, cloudc_outvar, clstr_outvar

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

        call check(nf90_def_var(ncid, "precnv", nf90_real4, (/ londim, latdim, timedim /), precnvvar))
        call check(nf90_put_att(ncid, precnvvar, "long_name", "convective_precipitation"))
        call check(nf90_put_att(ncid, precnvvar, "units", "g/(m^2 s)"))

        call check(nf90_def_var(ncid, "qdif", nf90_real4, (/ londim, latdim, timedim /), qdifvar))
        call check(nf90_put_att(ncid, qdifvar, "long_name", "qdif"))
        call check(nf90_put_att(ncid, qdifvar, "units", "idk"))

        call check(nf90_def_var(ncid, "evap_0", nf90_real4, (/ londim, latdim, timedim /), evapvar_0))
        call check(nf90_put_att(ncid, evapvar_0, "long_name", "evaporation_0"))
        call check(nf90_put_att(ncid, evapvar_0, "units", "idk"))
        
        call check(nf90_def_var(ncid, "evap_1", nf90_real4, (/ londim, latdim, timedim /), evapvar_1))
        call check(nf90_put_att(ncid, evapvar_1, "long_name", "evaporation_1"))
        call check(nf90_put_att(ncid, evapvar_1, "units", "idk"))
        
        call check(nf90_def_var(ncid, "evap_2", nf90_real4, (/ londim, latdim, timedim /), evapvar_2))
        call check(nf90_put_att(ncid, evapvar_2, "long_name", "evaporation_2"))
        call check(nf90_put_att(ncid, evapvar_2, "units", "idk"))

        call check(nf90_def_var(ncid, "denvvs_0", nf90_real4, (/ londim, latdim, timedim /), denvvsvar_0))
        call check(nf90_put_att(ncid, denvvsvar_0, "long_name", "denvvs_0"))
        call check(nf90_put_att(ncid, denvvsvar_0, "units", "idk"))
        
        call check(nf90_def_var(ncid, "denvvs_1", nf90_real4, (/ londim, latdim, timedim /), denvvsvar_1))
        call check(nf90_put_att(ncid, denvvsvar_1, "long_name", "denvvs_1"))
        call check(nf90_put_att(ncid, denvvsvar_1, "units", "idk"))
        
        call check(nf90_def_var(ncid, "denvvs_2", nf90_real4, (/ londim, latdim, timedim /), denvvsvar_2))
        call check(nf90_put_att(ncid, denvvsvar_2, "long_name", "denvvs_2"))
        call check(nf90_put_att(ncid, denvvsvar_2, "units", "idk"))

        call check(nf90_def_var(ncid, "iptop", nf90_real4, (/ londim, latdim, timedim /), iptopvar))
        call check(nf90_put_att(ncid, iptopvar, "long_name", "cloud top index"))
        call check(nf90_put_att(ncid, iptopvar, "units", "1"))

        call check(nf90_def_var(ncid, "psa_out", nf90_real4, (/ londim, latdim, timedim /), psa_outvar))
        call check(nf90_put_att(ncid, psa_outvar, "long_name", "psa_out"))
        call check(nf90_put_att(ncid, psa_outvar, "units", "idk"))

        call check(nf90_def_var(ncid, "se_out", nf90_real4, (/ londim, latdim, levdim, timedim /), se_outvar))
        call check(nf90_put_att(ncid, se_outvar, "long_name", "se_out"))
        call check(nf90_put_att(ncid, se_outvar, "units", "idk"))

        call check(nf90_def_var(ncid, "qa_out", nf90_real4, (/ londim, latdim, levdim, timedim /), qa_outvar))
        call check(nf90_put_att(ncid, qa_outvar, "long_name", "qa_out"))
        call check(nf90_put_att(ncid, qa_outvar, "units", "idk"))

        call check(nf90_def_var(ncid, "qsat_out", nf90_real4, (/ londim, latdim, levdim, timedim /), qsat_outvar))
        call check(nf90_put_att(ncid, qsat_outvar, "long_name", "qsat_out"))
        call check(nf90_put_att(ncid, qsat_outvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "mss_out", nf90_real4, (/ londim, latdim, levdim-1, timedim /), mss_outvar))
        ! call check(nf90_put_att(ncid, mss_outvar, "long_name", "mss_out"))
        ! call check(nf90_put_att(ncid, mss_outvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "mse0_out", nf90_real4, (/ londim, latdim, timedim /), mse0_outvar))
        ! call check(nf90_put_att(ncid, mse0_outvar, "long_name", "mse0_out"))
        ! call check(nf90_put_att(ncid, mse0_outvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "mse1_out", nf90_real4, (/ londim, latdim, timedim /), mse1_outvar))
        ! call check(nf90_put_att(ncid, mse1_outvar, "long_name", "mse1_out"))
        ! call check(nf90_put_att(ncid, mse1_outvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "mss0_out", nf90_real4, (/ londim, latdim, timedim /), mss0_outvar))
        ! call check(nf90_put_att(ncid, mss0_outvar, "long_name", "mss0_out"))
        ! call check(nf90_put_att(ncid, mss0_outvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "mss2_out", nf90_real4, (/ londim, latdim, timedim /), mss2_outvar))
        ! call check(nf90_put_att(ncid, mss2_outvar, "long_name", "mss2_out"))
        ! call check(nf90_put_att(ncid, mss2_outvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "msthr_out", nf90_real4, (/ londim, latdim, timedim /), msthr_outvar))
        ! call check(nf90_put_att(ncid, msthr_outvar, "long_name1", "msthr_out"))
        ! call check(nf90_put_att(ncid, msthr_outvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "qthr0out", nf90_real4, (/ londim, latdim, timedim /), qthr0outvar))
        ! call check(nf90_put_att(ncid, qthr0outvar, "long_name", "qthr0out"))
        ! call check(nf90_put_att(ncid, qthr0outvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "qthr1out", nf90_real4, (/ londim, latdim, timedim /), qthr1outvar))
        ! call check(nf90_put_att(ncid, qthr1outvar, "long_name", "qthr1out"))
        ! call check(nf90_put_att(ncid, qthr1outvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "ktop1out", nf90_real4, (/ londim, latdim, timedim /), ktop1outvar))
        ! call check(nf90_put_att(ncid, ktop1outvar, "long_name", "ktop1"))
        ! call check(nf90_put_att(ncid, ktop1outvar, "units", "1"))

        ! call check(nf90_def_var(ncid, "ktop2out", nf90_real4, (/ londim, latdim, timedim /), ktop2outvar))
        ! call check(nf90_put_att(ncid, ktop2outvar, "long_name", "ktop2"))
        ! call check(nf90_put_att(ncid, ktop2outvar, "units", "1"))

        call check(nf90_def_var(ncid, "tt_cnv", nf90_real4, (/ londim, latdim, levdim, timedim /), tt_cnvvar))
        call check(nf90_put_att(ncid, tt_cnvvar, "long_name", "tt_cnv"))
        call check(nf90_put_att(ncid, tt_cnvvar, "units", "idk"))

        call check(nf90_def_var(ncid, "qt_cnv", nf90_real4, (/ londim, latdim, levdim, timedim /), qt_cnvvar))
        call check(nf90_put_att(ncid, qt_cnvvar, "long_name", "qt_cnv"))
        call check(nf90_put_att(ncid, qt_cnvvar, "units", "idk"))

        call check(nf90_def_var(ncid, "phis0", nf90_real4, (/ londim, latdim, timedim /), phis0var))
        call check(nf90_put_att(ncid, phis0var, "long_name", "phis0"))
        call check(nf90_put_att(ncid, phis0var, "units", "idk"))

        call check(nf90_def_var(ncid, "fmask*compute_shortwave", nf90_real4, (/ londim, latdim, timedim /), fmaskvar))
        call check(nf90_put_att(ncid, fmaskvar, "long_name", "fmask"))
        call check(nf90_put_att(ncid, fmaskvar, "units", "idk"))

        call check(nf90_def_var(ncid, "tsea", nf90_real4, (/ londim, latdim, timedim /), tseavar))
        call check(nf90_put_att(ncid, tseavar, "long_name", "tsea"))
        call check(nf90_put_att(ncid, tseavar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "ssrd", nf90_real4, (/ londim, latdim, timedim /), ssrdvar))
        ! call check(nf90_put_att(ncid, ssrdvar, "long_name", "ssrd"))
        ! call check(nf90_put_att(ncid, ssrdvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "slrd", nf90_real4, (/ londim, latdim, timedim /), slrdvar))
        ! call check(nf90_put_att(ncid, slrdvar, "long_name", "slrd"))
        ! call check(nf90_put_att(ncid, slrdvar, "units", "idk"))

        call check(nf90_def_var(ncid, "shf_0", nf90_real4, (/ londim, latdim, timedim /), shfvar0))
        call check(nf90_put_att(ncid, shfvar0, "long_name", "shf_0"))
        call check(nf90_put_att(ncid, shfvar0, "units", "idk"))

        call check(nf90_def_var(ncid, "shf_1", nf90_real4, (/ londim, latdim, timedim /), shfvar1))
        call check(nf90_put_att(ncid, shfvar1, "long_name", "shf_1"))
        call check(nf90_put_att(ncid, shfvar1, "units", "idk"))

        call check(nf90_def_var(ncid, "shf_2", nf90_real4, (/ londim, latdim, timedim /), shfvar2))
        call check(nf90_put_att(ncid, shfvar2, "long_name", "shf_2"))
        call check(nf90_put_att(ncid, shfvar2, "units", "idk"))

        ! call check(nf90_def_var(ncid, "slru_0", nf90_real4, (/ londim, latdim, timedim /), slruvar0))
        ! call check(nf90_put_att(ncid, slruvar0, "long_name", "slru_0"))
        ! call check(nf90_put_att(ncid, slruvar0, "units", "idk"))

        ! call check(nf90_def_var(ncid, "slru_1", nf90_real4, (/ londim, latdim, timedim /), slruvar1))
        ! call check(nf90_put_att(ncid, slruvar1, "long_name", "slru_1"))
        ! call check(nf90_put_att(ncid, slruvar1, "units", "idk"))

        ! call check(nf90_def_var(ncid, "slru_2", nf90_real4, (/ londim, latdim, timedim /), slruvar2))
        ! call check(nf90_put_att(ncid, slruvar2, "long_name", "slru_2"))
        ! call check(nf90_put_att(ncid, slruvar2, "units", "idk"))

        call check(nf90_def_var(ncid, "hfluxn_0", nf90_real4, (/ londim, latdim, timedim /), hfluxnvar0))
        call check(nf90_put_att(ncid, hfluxnvar0, "long_name", "hfluxn_0"))
        call check(nf90_put_att(ncid, hfluxnvar0, "units", "idk"))

        call check(nf90_def_var(ncid, "hfluxn_1", nf90_real4, (/ londim, latdim, timedim /), hfluxnvar1))
        call check(nf90_put_att(ncid, hfluxnvar1, "long_name", "hfluxn_1"))
        call check(nf90_put_att(ncid, hfluxnvar1, "units", "idk"))

        call check(nf90_def_var(ncid, "tsfc", nf90_real4, (/ londim, latdim, timedim /), tsfcvar))
        call check(nf90_put_att(ncid, tsfcvar, "long_name", "tsfc"))
        call check(nf90_put_att(ncid, tsfcvar, "units", "idk"))

        call check(nf90_def_var(ncid, "tskin", nf90_real4, (/ londim, latdim, timedim /), tskinvar))
        call check(nf90_put_att(ncid, tskinvar, "long_name", "tskin"))
        call check(nf90_put_att(ncid, tskinvar, "units", "idk"))

        call check(nf90_def_var(ncid, "t0", nf90_real4, (/ londim, latdim, timedim /), t0var))
        call check(nf90_put_att(ncid, t0var, "long_name", "t0"))
        call check(nf90_put_att(ncid, t0var, "units", "idk"))

        ! call check(nf90_def_var(ncid, "stl_am", nf90_real4, (/ londim, latdim, timedim /), stl_amvar))
        ! call check(nf90_put_att(ncid, stl_amvar, "long_name", "stl_am"))
        ! call check(nf90_put_att(ncid, stl_amvar, "units", "idk"))

        call check(nf90_def_var(ncid, "soilw_am", nf90_real4, (/ londim, latdim, timedim /), soilw_amvar))
        call check(nf90_put_att(ncid, soilw_amvar, "long_name", "soilw_am"))
        call check(nf90_put_att(ncid, soilw_amvar, "units", "idk"))

        call check(nf90_def_var(ncid, "tskin_out_2", nf90_real4, (/ londim, latdim, timedim /), tskin_out_2var))
        call check(nf90_put_att(ncid, tskin_out_2var, "long_name", "tskin_out_2"))
        call check(nf90_put_att(ncid, tskin_out_2var, "units", "idk"))

        call check(nf90_def_var(ncid, "hfluxn_debug", nf90_real4, (/ londim, latdim, timedim /), hfluxn_debugvar))
        call check(nf90_put_att(ncid, hfluxn_debugvar, "long_name", "hfluxn_debug"))
        call check(nf90_put_att(ncid, hfluxn_debugvar, "units", "idk"))

        call check(nf90_def_var(ncid, "denvvs_debug", nf90_real4, (/ londim, latdim, timedim /), denvvs_debugvar))
        call check(nf90_put_att(ncid, denvvs_debugvar, "long_name", "denvvs_debug"))
        call check(nf90_put_att(ncid, denvvs_debugvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "qsat0_debug", nf90_real4, (/ londim, latdim, timedim /), qsat0_debugvar))
        ! call check(nf90_put_att(ncid, qsat0_debugvar, "long_name", "qsat0_debug"))
        ! call check(nf90_put_att(ncid, qsat0_debugvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "dtskin_debug", nf90_real4, (/ londim, latdim, timedim /), dtskin_debugvar))
        ! call check(nf90_put_att(ncid, dtskin_debugvar, "long_name", "dtskin_debug"))
        ! call check(nf90_put_att(ncid, dtskin_debugvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "rlus_debug", nf90_real4, (/ londim, latdim, timedim /), rlus_debugvar))
        ! call check(nf90_put_att(ncid, rlus_debugvar, "long_name", "rlus_debug"))
        ! call check(nf90_put_att(ncid, rlus_debugvar, "units", "idk"))

        call check(nf90_def_var(ncid, "shf_debug", nf90_real4, (/ londim, latdim, timedim /), shf_debugvar))
        call check(nf90_put_att(ncid, shf_debugvar, "long_name", "shf_debug"))
        call check(nf90_put_att(ncid, shf_debugvar, "units", "idk"))

        call check(nf90_def_var(ncid, "evap_debug", nf90_real4, (/ londim, latdim, timedim /), evap_debugvar))
        call check(nf90_put_att(ncid, evap_debugvar, "long_name", "evap_debug"))
        call check(nf90_put_att(ncid, evap_debugvar, "units", "idk"))

        call check(nf90_def_var(ncid, "t0_debug", nf90_real4, (/ londim, latdim, timedim /), t0_debugvar))
        call check(nf90_put_att(ncid, t0_debugvar, "long_name", "t0_debug"))
        call check(nf90_put_att(ncid, t0_debugvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "dthl_debug", nf90_real4, (/ londim, latdim, timedim /), dthl_debugvar))
        ! call check(nf90_put_att(ncid, dthl_debugvar, "long_name", "dthl_debug"))
        ! call check(nf90_put_att(ncid, dthl_debugvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "q_debug", nf90_real4, (/ londim, latdim, timedim /), q_debugvar))
        ! call check(nf90_put_att(ncid, q_debugvar, "long_name", "q_debug"))
        ! call check(nf90_put_att(ncid, q_debugvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "soilw_am_debug_final", nf90_real4, (/ londim, latdim, timedim /), soilw_am_debug_finalvar))
        ! call check(nf90_put_att(ncid, soilw_am_debug_finalvar, "long_name", "soilw_am_debug_final"))
        ! call check(nf90_put_att(ncid, soilw_am_debug_finalvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "qsat0_debug_final", nf90_real4, (/ londim, latdim, timedim /), qsat0_debug_finalvar))
        ! call check(nf90_put_att(ncid, qsat0_debug_finalvar, "long_name", "qsat0_debug_final"))
        ! call check(nf90_put_att(ncid, qsat0_debug_finalvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "q1_debug_final", nf90_real4, (/ londim, latdim, timedim /), q1_debug_finalvar))
        ! call check(nf90_put_att(ncid, q1_debug_finalvar, "long_name", "q1_debug_final"))
        ! call check(nf90_put_att(ncid, q1_debug_finalvar, "units", "idk"))

        call check(nf90_def_var(ncid, "tskin_debug_final", nf90_real4, (/ londim, latdim, timedim /), tskin_debug_finalvar))
        call check(nf90_put_att(ncid, tskin_debug_finalvar, "long_name", "tskin_debug_final"))
        call check(nf90_put_att(ncid, tskin_debug_finalvar, "units", "idk"))

        ! call check(nf90_def_var(ncid, "psa_debug_final", nf90_real4, (/ londim, latdim, timedim /), psa_debug_finalvar))
        ! call check(nf90_put_att(ncid, psa_debug_finalvar, "long_name", "psa_debug_final"))
        ! call check(nf90_put_att(ncid, psa_debug_finalvar, "units", "idk"))

        call check(nf90_def_var(ncid, "tt_lsc", nf90_real4, (/ londim, latdim, levdim, timedim /), tt_lscvar))
        call check(nf90_put_att(ncid, tt_lscvar, "long_name", "tt_lsc"))
        call check(nf90_put_att(ncid, tt_lscvar, "units", "idk"))

        call check(nf90_def_var(ncid, "tt_rsw", nf90_real4, (/ londim, latdim, levdim, timedim /), tt_rswvar))
        call check(nf90_put_att(ncid, tt_rswvar, "long_name", "tt_rsw"))
        call check(nf90_put_att(ncid, tt_rswvar, "units", "idk"))

        call check(nf90_def_var(ncid, "tt_rlw", nf90_real4, (/ londim, latdim, levdim, timedim /), tt_rlwvar))
        call check(nf90_put_att(ncid, tt_rlwvar, "long_name", "tt_rlw"))
        call check(nf90_put_att(ncid, tt_rlwvar, "units", "idk"))

        call check(nf90_def_var(ncid, "tt_pbl_out", nf90_real4, (/ londim, latdim, levdim, timedim /), tt_pbl_outvar))
        call check(nf90_put_att(ncid, tt_pbl_outvar, "long_name", "tt_pbl_out"))
        call check(nf90_put_att(ncid, tt_pbl_outvar, "units", "idk"))

        call check(nf90_def_var(ncid, "tt_pbl", nf90_real4, (/ londim, latdim, levdim, timedim /), tt_pblvar))
        call check(nf90_put_att(ncid, tt_pblvar, "long_name", "tt_pbl"))
        call check(nf90_put_att(ncid, tt_pblvar, "units", "idk"))

        call check(nf90_def_var(ncid, "stl_am_a", nf90_real4, (/ londim, latdim, timedim /), stl_am_avar))
        call check(nf90_put_att(ncid, stl_am_avar, "long_name", "stl_am_a"))
        call check(nf90_put_att(ncid, stl_am_avar, "units", "idk"))

        call check(nf90_def_var(ncid, "coa_factor_a", nf90_real4, (/ latdim, timedim /), coa_factor_avar))
        call check(nf90_put_att(ncid, coa_factor_avar, "long_name", "coa_factor_a"))
        call check(nf90_put_att(ncid, coa_factor_avar, "units", "idk"))
        
        call check(nf90_def_var(ncid, "ssrd_a", nf90_real4, (/ londim, latdim, timedim /), ssrd_avar))
        call check(nf90_put_att(ncid, ssrd_avar, "long_name", "ssrd_a"))
        call check(nf90_put_att(ncid, ssrd_avar, "units", "idk"))
        
        call check(nf90_def_var(ncid, "alb_l_factor_a", nf90_real4, (/ londim, latdim, timedim /), alb_l_factor_avar))
        call check(nf90_put_att(ncid, alb_l_factor_avar, "long_name", "alb_l_factor_a"))
        call check(nf90_put_att(ncid, alb_l_factor_avar, "units", "idk"))
        
        call check(nf90_def_var(ncid, "psa_a", nf90_real4, (/ londim, latdim, timedim /), psa_avar))
        call check(nf90_put_att(ncid, psa_avar, "long_name", "psa_a"))
        call check(nf90_put_att(ncid, psa_avar, "units", "idk"))

        call check(nf90_def_var(ncid, "flux_out_a0", nf90_real4, (/ londim, latdim, levdim, timedim /), flux_out_a0var))
        call check(nf90_put_att(ncid, flux_out_a0var, "long_name", "flux_out_a0"))
        call check(nf90_put_att(ncid, flux_out_a0var, "units", "idk"))

        call check(nf90_def_var(ncid, "flux_out_a1", nf90_real4, (/ londim, latdim, levdim, timedim /), flux_out_a1var))
        call check(nf90_put_att(ncid, flux_out_a1var, "long_name", "flux_out_a1"))
        call check(nf90_put_att(ncid, flux_out_a1var, "units", "idk"))

        call check(nf90_def_var(ncid, "qcloud", nf90_real4, (/ londim, latdim, timedim /), qcloud_outvar))
        call check(nf90_put_att(ncid, qcloud_outvar, "long_name", "qcloud"))
        call check(nf90_put_att(ncid, qcloud_outvar, "units", "idk"))

        call check(nf90_def_var(ncid, "icltop", nf90_real4, (/ londim, latdim, timedim /), icltop_outvar))
        call check(nf90_put_att(ncid, icltop_outvar, "long_name", "icltop"))
        call check(nf90_put_att(ncid, icltop_outvar, "units", "idk"))

        call check(nf90_def_var(ncid, "icltop_out", nf90_real4, (/ londim, latdim, timedim /), icltop_2var))
        call check(nf90_put_att(ncid, icltop_2var, "long_name", "icltop_out"))
        call check(nf90_put_att(ncid, icltop_2var, "units", "idk"))

        call check(nf90_def_var(ncid, "cloudc", nf90_real4, (/ londim, latdim, timedim /), cloudc_outvar))
        call check(nf90_put_att(ncid, cloudc_outvar, "long_name", "cloudc"))
        call check(nf90_put_att(ncid, cloudc_outvar, "units", "idk"))

        call check(nf90_def_var(ncid, "clstr", nf90_real4, (/ londim, latdim, timedim /), clstr_outvar))
        call check(nf90_put_att(ncid, clstr_outvar, "long_name", "clstr"))
        call check(nf90_put_att(ncid, clstr_outvar, "units", "idk"))

        call check(nf90_def_var(ncid, "rh", nf90_real4, (/ londim, latdim, levdim, timedim /), rhvar))
        call check(nf90_put_att(ncid, rhvar, "long_name", "relative_humidity"))
        call check(nf90_put_att(ncid, rhvar, "units", "1"))

        call check(nf90_def_var(ncid, "phi", nf90_real4, (/ londim, latdim, levdim, timedim /), &
            & phivar))
        call check(nf90_put_att(ncid, phivar, "long_name", "geopotential_height"))
        call check(nf90_put_att(ncid, phivar, "units", "m"))

        call check(nf90_def_var(ncid, "ps", nf90_real4, (/ londim, latdim, timedim /), psvar))
        call check(nf90_put_att(ncid, psvar, "long_name", "surface_air_pressure"))
        call check(nf90_put_att(ncid, psvar, "units", "Pa"))

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
        call check(nf90_put_var(ncid, rhvar, rh_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, precnvvar, precnv, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qdifvar, qdif, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, evapvar_0, evap(:,:,1), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, evapvar_1, evap(:,:,2), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, evapvar_2, evap(:,:,3), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, denvvsvar_0, denvvs_out(:,:,0), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, denvvsvar_1, denvvs_out(:,:,1), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, denvvsvar_2, denvvs_out(:,:,2), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, iptopvar, iptop_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, psa_outvar, psa_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, se_outvar, se_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qa_outvar, qa_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qsat_outvar, qsat_out, (/ 1, 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, mss_outvar, mss_out, (/ 1, 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, mse0_outvar, mse0_out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, mse1_outvar, mse1_out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, mss0_outvar, mss0_out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, mss2_outvar, mss2_out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, msthr_outvar, msthr_out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, qthr0outvar, qthr0out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, qthr1outvar, qthr1out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, ktop1outvar, ktop1out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, ktop2outvar, ktop2out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tt_cnvvar, tt_cnv, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qt_cnvvar, qt_cnv, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, phis0var, phis0_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, fmaskvar, fmask_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tseavar, tsea_out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, ssrdvar, ssrd_out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, slrdvar, slrd_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, shfvar0, shf_out(:,:,1), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, shfvar1, shf_out(:,:,2), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, shfvar2, shf_out(:,:,3), (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, slruvar0, slru_out(:,:,1), (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, slruvar1, slru_out(:,:,2), (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, slruvar2, slru_out(:,:,3), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, hfluxnvar0, hfluxn_out(:,:,1), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, hfluxnvar1, hfluxn_out(:,:,2), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tsfcvar, tsfc_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tskinvar, tskin_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, t0var, t0_out, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, stl_amvar, stl_am_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, soilw_amvar, soilw_am_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tskin_out_2var, tskin_out_2, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, hfluxn_debugvar, hfluxn_debug, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, denvvs_debugvar, denvvs_debug, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, qsat0_debugvar, qsat0_debug, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, dtskin_debugvar, dtskin_debug, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, rlus_debugvar, rlus_debug, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, shf_debugvar, shf_debug, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, evap_debugvar, evap_debug, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, t0_debugvar, t0_debug, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, dthl_debugvar, dthl_debug, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, q_debugvar, q_debug, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, soilw_am_debug_finalvar, soilw_am_debug_final, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, qsat0_debug_finalvar, qsat0_debug_final, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, q1_debug_finalvar, q1_debug_final, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tskin_debug_finalvar, tskin_debug_final, (/ 1, 1, 1 /)))
        ! call check(nf90_put_var(ncid, psa_debug_finalvar, psa_debug_final, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tt_lscvar, tt_lsc, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tt_rswvar, tt_rsw, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tt_rlwvar, tt_rlw, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tt_pbl_outvar, tt_pbl_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, tt_pblvar, tt_pbl, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, stl_am_avar, stl_am_a, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, coa_factor_avar, coa_factor_a, (/ 1, 1 /)))
        call check(nf90_put_var(ncid, ssrd_avar, ssrd_a, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, alb_l_factor_avar, alb_l_factor_a, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, psa_avar, psa_a, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, flux_out_a0var, flux_out_a(:,:,:,1), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, flux_out_a1var, flux_out_a(:,:,:,2), (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qcloud_outvar, qcloud_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, icltop_outvar, icltop_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, icltop_2var, icltop_2, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, cloudc_outvar, cloudc_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, clstr_outvar, clstr_out, (/ 1, 1, 1 /)))

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
