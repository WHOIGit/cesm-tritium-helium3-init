program tritium_dep_test

    !---------------------------------------------------------------------------
    ! Read GNIP tritium deposition field
    ! Created by Ivan Lima on Mon, 16 Oct 2017 11:06:30 -0400
    !---------------------------------------------------------------------------

    use netcdf
    implicit none

    integer, parameter ::                            &
            char_len       = 160,                    &
            int_kind       = kind(1),                &
            log_kind       = kind(.true.),           &
            long_int_kind  = selected_int_kind(10),  &
            dbl_kind       = selected_real_kind(12), &
            i4             = selected_int_kind(6),   &
            i8             = selected_int_kind(13),  &
            r4             = selected_real_kind(6),  &
            r8             = selected_real_kind(13), &
            POP_CharLength = 256,                    &
            POP_Logical    = kind(.true.),           &
            POP_i4         = selected_int_kind(6),   &
            POP_i8         = selected_int_kind(13),  &
            POP_r4         = selected_real_kind(6),  &
            POP_r8         = selected_real_kind(13), &
            POP_r16        = selected_real_kind(26), &
            nlat           = 17,                     &
            ntime          = 733,                    &
            nx_block       = 320,                    & ! nlon
            ny_block       = 384                       ! ntlat

    real (r8), parameter ::        &
            c0     =     0.0_r8   ,&
            c1     =     1.0_r8   ,&
            c2     =     2.0_r8   ,&
            c3     =     3.0_r8   ,&
            c4     =     4.0_r8   ,&
            c5     =     5.0_r8   ,&
            c8     =     8.0_r8   ,&
            c10    =    10.0_r8   ,&
            c16    =    16.0_r8   ,&
            c1000  =  1000.0_r8   ,&
            c10000 = 10000.0_r8   ,&
            c1p5   =     1.5_r8   ,&
            p33    =  c1/c3       ,&
            p5     =  0.500_r8    ,&
            p25    =  0.250_r8    ,&
            p125   =  0.125_r8    ,&
            p001   =  0.001_r8    ,&
            eps    =  1.0e-10_r8  ,&
            eps2   =  1.0e-20_r8  ,&
            bignum =  1.0e+30_r8

    real, parameter :: rad2deg = 57.29577951308232_r8

    character(char_len) ::                                                                                        &
            tritium_dep_file ,                                                                                    &
            horiz_grid_file  = '/bali/data/ilima/cesm/cesm1_input/ocn/pop/gx1v6/grid/horiz_grid_20010402.ieeer8', &
            nc_file_name     = '/home/ivan/Python/data/gx1v6.nc'

    integer (int_kind) :: &
            model_year,   & ! arbitrary model year
            data_year,    & ! year in data that corresponds to model_year
            gnip_data_len   ! length of atmospheric tritium deposition record

    real (r8), parameter :: &
            max_gnip_extension = 2.0_r8
    ! maximum number of years that gnip record will be extrapolated

    real (r8), dimension(nlat)                 :: gnip_lat     ! lower bound of latitude bins
    real (r8), dimension(ntime)                :: gnip_date    ! date/time in decimal years
    real (r8), dimension(nlat,ntime)           :: gnip_tritium ! tritium concentration in precipitation (TU)
    real (POP_r8), dimension(:,:), allocatable :: ULAT, ULON   ! latitude,longitude of U points
    real (r8), dimension(nx_block, ny_block)   :: KMT, TLAT

    logical (log_kind), dimension(nx_block, ny_block) :: LAND_MASK

    integer (int_kind) ::      &
            iyear        = 1,  &
            iday_of_year = 1,  &
            days_in_year = 365

    real (r8) :: frac_day = 0.0_r8

    ! integer (int_kind), dimension(:), allocatable :: &
    !         data_ind
    integer (int_kind) :: data_ind ! index into data for current timestep

    ! call read_horiz_grid(horiz_grid_file)
    ! print *, minval(ULAT), maxval(ULAT), shape(ULAT)
    ! print *, minval(ULON), maxval(ULON), shape(ULON)

    call read_netcdf(nc_file_name)
    ! print *, minval(KMT), maxval(KMT), shape(KMT)
    ! print *, minval(TLAT), maxval(TLAT), shape(TLAT)

    LAND_MASK = (KMT.gt.0)

    call tr3he_init
    call tr3he_init_sflux
    do iyear = 1, 10
        call tr3he_set_sflux
    enddo

contains

    !***************************************************************************

    subroutine tr3he_init

        tritium_dep_file = 'tritium_dep_gnip.txt' ! file name for tritium deposition data
        model_year       = 1
        data_year        = 1950

    end subroutine tr3he_init

    !***************************************************************************

    subroutine tr3he_init_sflux

        call read_tritium_dep_data

    end subroutine tr3he_init_sflux

    !***************************************************************************

    subroutine read_tritium_dep_data

        integer :: ios, ix, nu
        character(char_len) :: header
        character(9) :: label

        nu = 10

        print *, 'Reading ', tritium_dep_file
        open(nu, file=trim(tritium_dep_file), access='sequential', &
                form="formatted", iostat=ios)
        read(nu,'(a)') header
        read(nu,*) label, gnip_lat(:)
        read(nu,'(a)') header
        ! print '(a9,17f6.1)', label, gnip_lat(:)
        do ix = 1, ntime
            read(nu,*) gnip_date(ix), gnip_tritium(:,ix)
        end do
        close(nu)
        gnip_data_len = size(gnip_date)
        ! print '(a,733f8.2)', 'gnip_date: ', gnip_date
        ! print *, gnip_date
        ! print *, gnip_date(165), gnip_tritium(:,165)

    end subroutine read_tritium_dep_data

    !***************************************************************************

    subroutine tr3he_set_sflux

        !-----------------------------------------------------------------------
        !  local variables
        !-----------------------------------------------------------------------

        real (r8), dimension(nx_block,ny_block) :: &
                tritium_dep  ! tritium concentration in precipitation

        logical (log_kind), save :: first = .true.

        if (first) then
            ! allocate( data_ind(max_blocks_clinic) )
            data_ind = -1
            first = .false.
        endif

            ! call comp_tritium_dep(iblock, LAND_MASK(:,:,iblock), data_ind(iblock), &
            !         tritium_dep)
        call comp_tritium_dep(LAND_MASK, data_ind, tritium_dep)

    end subroutine tr3he_set_sflux

    !***************************************************************************

    ! subroutine comp_tritium_dep(iblock, LAND_MASK, data_ind, tritium_dep)
    subroutine comp_tritium_dep(LAND_MASK, data_ind, tritium_dep)

        ! Set tritium concentration in precipitation according to latitude and date

        ! use time_management, only : iyear, iday_of_year, frac_day, days_in_year

        !-----------------------------------------------------------------------
        ! INPUT PARAMETERS:
        !-----------------------------------------------------------------------

        logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
                LAND_MASK          ! land mask for this block

        integer (int_kind) :: &
                iblock          ! block index

        !-----------------------------------------------------------------------
        ! INPUT/OUTPUT PARAMETERS:
        !-----------------------------------------------------------------------

        integer (int_kind) :: &
                data_ind ! data_ind is the index into data for current timestep,
                         ! i.e data_ind is largest integer less than pcfc_data_len s.t.
                         !  pcfc_date(i) <= iyear + (iday_of_year-1+frac_day)/days_in_year
                         !                  - model_year + data_year
                         !  note that data_ind is always strictly less than pcfc_data_len
                         !  and is initialized to -1 before the first call

        ! !OUTPUT PARAMETERS:
        real (r8), dimension(nx_block,ny_block), intent(out) :: &
                tritium_dep ! tritium concentration in precipitation (TU)

        !-----------------------------------------------------------------------
        !  local variables
        !-----------------------------------------------------------------------

        integer (int_kind) :: &
                i, j              ! loop indices

        real (r8) :: &
                mapped_date,    & ! date of current model timestep mapped to data timeline
                weight            ! weighting for temporal interpolation

        real (r8), dimension(nlat) :: gnip_tritium_curr ! tritium concentration for current time step

        !-----------------------------------------------------------------------
        !  Generate mapped_date and check to see if it is too large.
        !  The check for mapped_date being too small only needs to be done
        !  on the first time step.
        !-----------------------------------------------------------------------

        mapped_date = iyear + 0.45 + (iday_of_year-1+frac_day)/days_in_year &
                - model_year + data_year

        ! if (mapped_date >= gnip_date(gnip_data_len) + max_gnip_extension) &
        !         print *, 'exit_POP(sigAbort, model date maps too far beyond pcfc_date(end))'

        !-----------------------------------------------------------------------
        ! Set tritium concentrations before GNIP record
        !-----------------------------------------------------------------------

        if (mapped_date < gnip_date(1)) then
            data_ind = 1
            gnip_tritium_curr = gnip_tritium(:,1)
            return
        endif

        !-----------------------------------------------------------------------
        !  On first time step, perform linear search to find data_ind.
        !-----------------------------------------------------------------------

        if (data_ind == -1) then
            do data_ind = gnip_data_len-1,1,-1
                if (mapped_date >= gnip_date(data_ind)) exit
            end do
        endif

        !-----------------------------------------------------------------------
        !  See if data_ind need to be updated,
        !  but do not set it to gnip_data_len.
        !-----------------------------------------------------------------------

        if (data_ind < gnip_data_len) then
            ! if (mapped_date >= gnip_date(data_ind+1)) data_ind = data_ind + 1
            do data_ind = gnip_data_len,1,-1
                if (mapped_date >= gnip_date(data_ind)) exit
            end do
        endif

        !-----------------------------------------------------------------------
        !  Generate tritium concentrations for current time step.
        !-----------------------------------------------------------------------

        if (data_ind < gnip_data_len) then
            weight = (mapped_date - gnip_date(data_ind)) &
                    / (gnip_date(data_ind+1) - gnip_date(data_ind))

            gnip_tritium_curr = &
                    weight * gnip_tritium(:,data_ind+1) + (c1-weight) * gnip_tritium(:,data_ind)
        else
            gnip_tritium_curr = gnip_tritium(:,gnip_data_len)
        endif

        ! print *, gnip_data_len, data_ind
        print '(2f8.2)', mapped_date, gnip_date(data_ind)
        print '(17f9.3)', gnip_tritium(:,data_ind)
        print '(17f9.3)', gnip_tritium_curr

        do i=1, nlat
            ! print '(f6.1, f9.3)', gnip_lat(i), gnip_tritium_curr(i)
            where (TLAT > gnip_lat(i) .and. LAND_MASK) tritium_dep = gnip_tritium_curr(i)
        enddo

        ! open(20, file='dep.txt', access='sequential', form="formatted")
        ! write(20,*) tritium_dep
        ! close(20)

    end subroutine comp_tritium_dep

    !***************************************************************************

    subroutine read_netcdf(file_name)

        ! Read variables from netCDF file

        character (*), intent(in) :: file_name ! netCDF file name
        integer :: ncid, varid, x, y, ierr

        ierr = nf90_open(file_name, NF90_NOWRITE, ncid) ! open file
        ierr = nf90_inq_varid(ncid, "KMT", varid)       ! read KMT
        ierr = nf90_get_var(ncid, varid, KMT)
        ierr = nf90_inq_varid(ncid, "TLAT", varid)      ! read TLAT
        ierr = nf90_get_var(ncid, varid, TLAT)
        ierr = nf90_close(ncid)                         ! close file

    end subroutine read_netcdf

    !***************************************************************************

    subroutine read_horiz_grid(horiz_grid_file)

        ! Reads horizontal grid information from input grid file

        character (*), intent(in) :: horiz_grid_file ! filename of file containing grid data

        integer (i4) :: &
                nu,     & ! i/o unit number
                ioerr,  & ! i/o error flag
                reclength ! record length

        if (.not. allocated(ULAT)) then
            allocate (ULAT(nx_block,ny_block), ULON(nx_block,ny_block))
        endif

        inquire(iolength=reclength) ULAT
        nu    = 20
        ioerr = 0
        open(nu,file=trim(horiz_grid_file),status='old', form='unformatted', &
                access='direct', recl=reclength, iostat=ioerr)
        read(nu,rec=1,iostat=ioerr) ULAT
        read(nu,rec=2,iostat=ioerr) ULON
        close(nu)
        ULAT = ULAT * rad2deg ! convert from radians to degrees
        ULON = ULON * rad2deg

    end subroutine read_horiz_grid

    !***************************************************************************

end program tritium_dep_test
