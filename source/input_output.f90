!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 08/05/2019
!  For performing input and output.
module input_output
    use types, only: p, sp
    use netcdf
    use params
    use physics, only: rh
    use auxiliaries, only: precnv, qdif, evap, denvvs_out, iptop, psa_out, se_out, qa_out, qsat_out

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
        integer :: timevar, latvar, lonvar, levvar, uvar, vvar, tvar, qvar, phivar, psvar, rhvar, precnvvar, qdifvar, evapvar_0, evapvar_1, evapvar_2, denvvsvar_0, denvvsvar_1, denvvsvar_2, iptopvar, psa_outvar, se_outvar, qa_outvar, qsat_outvar

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
        call check(nf90_put_var(ncid, rhvar, rh, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, precnvvar, precnv, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qdifvar, qdif, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, evapvar_0, evap(:,:,1), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, evapvar_1, evap(:,:,2), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, evapvar_2, evap(:,:,3), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, denvvsvar_0, denvvs_out(:,:,0), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, denvvsvar_1, denvvs_out(:,:,1), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, denvvsvar_2, denvvs_out(:,:,2), (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, iptopvar, iptop, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, psa_outvar, psa_out, (/ 1, 1, 1 /)))
        call check(nf90_put_var(ncid, se_outvar, se_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qa_outvar, qa_out, (/ 1, 1, 1, 1 /)))
        call check(nf90_put_var(ncid, qsat_outvar, qsat_out, (/ 1, 1, 1, 1 /)))

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
