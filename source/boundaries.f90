!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 29/04/2019
!  For reading and storing boundary conditions.
!
!  This module handles the reading and storing of the land-sea mask, the
!  surface geopotential (i.e. the orography), the filtered surface geopotential
!  (i.e. the smoothed orography) and the bare-land annual-mean albedo.
module boundaries
    use types, only: p, sp
    use params
    use netcdf

    implicit none

    private
    public initialize_boundaries, fillsf, forchk
    public fmask, phi0, phis0, alb0
    public load_boundary_file

    real(p) :: fmask(ix,il) !! Original (fractional) land-sea mask

    ! Time invariant surface fields
    real(p) :: phi0(ix,il)  !! Unfiltered surface geopotential
    real(p) :: phis0(ix,il) !! Spectrally-filtered surface geopotential
    real(p) :: alb0(ix,il)  !! Bare-land annual-mean albedo

    !> Interface for reading boundary files.
    interface load_boundary_file
    module procedure load_boundary_file_2d
    module procedure load_boundary_file_one_month_from_year
    module procedure load_boundary_file_one_month_from_long
    end interface

contains
    !> Loads the given 2D field from the given boundary file.
    function load_boundary_file_2d(file_name, field_name) result(field)
        character(len=*), intent(in) :: file_name  !! The NetCDF file to read from
        character(len=*), intent(in) :: field_name !! The field to read

        integer :: ncid, varid
        real(sp), dimension(ix,il) :: raw_input
        real(p), dimension(ix,il)  :: field

        ! Open boundary file, read variable and then close
        call check(nf90_open(file_name, nf90_nowrite, ncid))
        call check(nf90_inq_varid(ncid, field_name, varid))
        call check(nf90_get_var(ncid, varid, raw_input, start = (/ 1, 1 /), count =  (/ ix, il /)))
        call check(nf90_close(ncid))
        field = raw_input(:,il:1:-1)

        ! Fix undefined values
        where (field <= -999) field = 0.0
    end function

    !> Loads the given 2D field at the given month from the given monthly
    !  boundary file.
    function load_boundary_file_one_month_from_year(file_name, field_name, month) result(field)
        character(len=*), intent(in) :: file_name  !! The NetCDF file to read from
        character(len=*), intent(in) :: field_name !! The field to read
        integer, intent(in)          :: month      !! The month to read

        integer :: ncid, varid
        real(sp), dimension(ix,il,12) :: raw_input
        real(p), dimension(ix,il)     :: field

        ! Open boundary file, read variable and then close
        call check(nf90_open(file_name, nf90_nowrite, ncid))
        call check(nf90_inq_varid(ncid, field_name, varid))
        call check(nf90_get_var(ncid, varid, raw_input))
        call check(nf90_close(ncid))
        field = raw_input(:,il:1:-1,month)

        ! Fix undefined values
        where (field <= -999) field = 0.0
    end

    !> Loads the given 2D field at the given month from the given boundary file
    !  of a given length.
    !
    !  This is used for reading the SST anomalies from a particular month of a
    !  particular year. The SST anomalies are stored in a long multidecadal
    !  file and the total number of months in this file must be passed as an
    !  argument (`length`).
    function load_boundary_file_one_month_from_long(file_name, field_name, month, length) &
        & result(field)
        character(len=*), intent(in) :: file_name  !! The NetCDF file to read from
        character(len=*), intent(in) :: field_name !! The field to read
        integer, intent(in)          :: month      !! The month to read
        integer, intent(in)          :: length     !! The total length of the file in number of
                                                !! months

        integer :: ncid, varid
        real(sp), dimension(ix,il,length) :: raw_input
        real(p), dimension(ix,il)         :: field

        ! Open boundary file, read variable and then close
        call check(nf90_open(file_name, nf90_nowrite, ncid))
        call check(nf90_inq_varid(ncid, field_name, varid))
        call check(nf90_get_var(ncid, varid, raw_input))
        call check(nf90_close(ncid))
        field = raw_input(:,il:1:-1,month)

        ! Fix undefined values
        where (field <= -999) field = 0.0
    end


    !> Initialize boundary conditions (land-sea mask, surface geopotential
    !  and surface albedo).
    subroutine initialize_boundaries
        use physical_constants, only: grav

        ! Read surface geopotential (i.e. orography)
        phi0 = grav*load_boundary_file("surface.nc", "orog")

        ! Also store spectrally truncated surface geopotential
        call spectral_truncation(phi0, phis0)

        ! Read land-sea mask
        fmask = load_boundary_file("surface.nc", "lsm")

        ! Annual-mean surface albedo
        alb0 = load_boundary_file("surface.nc", "alb")
    end subroutine

    !> Check consistency of surface fields with land-sea mask and set undefined
    !  values to a constant (to avoid over/underflow).
    subroutine forchk(fmask, nf, fmin, fmax, fset, field)
        real(p), intent(in)    :: fmask(ix,il)    !! The fractional land-sea mask
        integer, intent(in)    :: nf              !! The number of input 2D fields
        real(p), intent(in)    :: fmin            !! The minimum allowable value
        real(p), intent(in)    :: fmax            !! The maximum allowable value
        real(p), intent(in)    :: fset            !! Replacement for undefined values
        real(p), intent(inout) :: field(ix,il,nf) !! The output field

        integer :: i, j, jf, nfault

        do jf = 1, nf
            nfault = 0

            do i = 1, ix
                do j = 1, il
                    if (fmask(i,j) > 0.0) then
                        if (field(i,j,jf) < fmin .or. field(i,j,jf) > fmax) then
                            nfault = nfault + 1
                        end if
                    else
                        field(i,j,jf) = fset
                    end if
                end do
            end do
        end do
    end subroutine

    !> Compute a spectrally-filtered grid-point field.
    subroutine spectral_truncation(fg1, fg2)
        use spectral, only: grid_to_spec, spec_to_grid

        real(p), intent(inout) :: fg1(ix,il) !! Original grid-point field
        real(p), intent(inout) :: fg2(ix,il) !! Filtered grid-point field

        complex(p) :: fsp(mx,nx)
        integer :: n, m, total_wavenumber

        fsp = grid_to_spec(fg1)

        do n = 1, nx
            do m = 1, mx
                total_wavenumber = m + n - 2
                if (total_wavenumber > trunc) fsp(m,n) = (0.0, 0.0)
            end do
        end do

        fg2 = spec_to_grid(fsp, 1)
    end subroutine

    !> Replace missing values in surface fields.
    !  @note It is assumed that non-missing values exist near the Equator.
    subroutine fillsf(sf, fmis)
        real(p), intent(inout) :: sf(ix,il) !! Field to replace missing values in
        real(p), intent(in)    :: fmis      !! Replacement for missing values

        real(p) :: sf2(0:ix+1)

        integer :: hemisphere, j, j1, j2, j3, i, nmis
        real(p) :: fmean

        do hemisphere = 1, 2
           if (hemisphere == 1) then
                j1 = il/2
                j2 = 1
                j3 = -1
            else
                j1 = j1+1
                j2 = il
                j3 = 1
            end if

            do j = j1, j2, j3
                sf2(1:ix) = sf(:,j)

                nmis = 0
                do i = 1, ix
                    if (sf(i,j) < fmis) then
                        nmis = nmis + 1
                        sf2(i) = 0.0
                    end if
                end do

                if (nmis < ix) fmean = sum(sf2(1:ix))/float(ix - nmis)

                do i = 1, ix
                    if (sf(i,j) < fmis) sf2(i) = fmean
                end do

                sf2(0)      = sf2(ix)
                sf2(ix+1) = sf2(1)
                do i = 1, ix
                    if (sf(i,j) < fmis) sf(i,j) = 0.5*(sf2(i-1) + sf2(i+1))
                end do
            end do
        end do
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
