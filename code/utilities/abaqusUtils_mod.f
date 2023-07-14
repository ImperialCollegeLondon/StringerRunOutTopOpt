!DEC$ FREEFORM

module abaqusUtils_mod
    use iso_fortran_env
    use abaqusXIT_mod
    use stringUtils_mod
    implicit none

    character(len=50), parameter :: keywords(25) =  &
        [ Character(len=50) ::                      &
            "NNODEL",                               &
            "COORDS",                               &
            "UELS",                                 &
            "UELSG",                                &
            "CON",                                  &
            "BOUNDARY",                             &
            "PROPS",                                &
            "ORIENTATIONS",                         &
            "THICKNESSES",                          &
            "TRANSFORM_SHELL",                      &
            "TRANSFORM_LS",                         &
            "NDSG",                                 &
            "NNODESUEL",                            &
            "NNODESLS",                             &
            "LSNDS",                                &
            "LSINIT",                               &
            "MESHSIZE",                             &
            "ALPHA",                                &
            "VFRAC",                                &
            "ADJCOORDS",                            &
            "ADJUELS",                              &
            "ADJUELSG",                             &
            "ADJCON",                               &
            "ADJNDSG",                              &
            "NNODESADJUEL"                          &
        ]

    type :: INPFile 
        integer                   :: nnodes         = 0
        integer                   :: nnodesls       = 0

        integer                   :: nnodesuel      = 0
        integer                   :: nnodesadjuel   = 0

        integer                   :: nuels          = 0
        integer                   :: nadjuels       = 0

        integer                   :: nndeluels      = 0
        integer                   :: nndeladjuels   = 0

        logical                   :: hasBoundary = .false.
        logical                   :: hasDisbond  = .false.

        real(real64), allocatable :: x(:,:)
        integer,      allocatable :: con(:,:)
        integer,      allocatable :: uels(:)
        integer,      allocatable :: uels_global(:)
        integer,      allocatable :: nds_global(:)

        real(real64), allocatable :: props(:)
        real(real64), allocatable :: orientations(:)
        real(real64), allocatable :: thicknesses(:)
        real(real64)              :: transform_shell(3,3)   = 0.0d0

        real(real64)              :: transform_ls(3,3)      = 0.0d0
        real(real64)              :: lsorigin(3)            = 0.0d0
        integer,      allocatable :: lsnodes(:)
        real(real64), allocatable :: lsinit(:)

        integer,      allocatable :: boundary(:)
        
        real(real64)              :: meshsize               = 0.0d0
        real(real64)              :: alpha                  = 0.0d0
        real(real64)              :: vfrac                  = 0.0d0
       
        real(real64), allocatable :: adjx(:,:)
        integer,      allocatable :: adjcon(:,:)
        integer,      allocatable :: adjuels(:)
        integer,      allocatable :: adjuels_global(:)
        integer,      allocatable :: adjnds_global(:)

        real(real64), allocatable :: k(:,:,:)
        real(real64), allocatable :: dkda(:,:,:)

    contains
        procedure :: finalize => inpfile_destroy

        ! procedure :: isInit => inpfile_isInit

        final :: inpfile_dtor
    end type

    interface INPFile 
        procedure :: inpfile_ctor
    end interface

contains

    function inpfile_ctor(filename, outdir) result(this)
        implicit none 

        type(INPFile) :: this 
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: outdir

        integer                         :: i, ifl, sf1, sf2
        integer                         :: ios, commaIdx 
        integer                         :: nnodel(2), auxId, auxN, dtyp(1)
        integer, allocatable            :: con(:), adjcon(:) 
        real(real64)                    :: coords(3), auxR
        character(len=5096)             :: lineRaw
        character(len=:), allocatable   :: line
        character(len=:), allocatable   :: keyword
        character(len=50)               :: kwd
        integer                         :: dataType
        integer                         :: eid, row, col, elid(1)
        real(real64)                    :: val

        ! vars init
        ios         = 0
        commaIdx    = 0
        nnodel      = 0
        auxId       = 0
        auxN        = 0
        dtyp        = 0
        coords      = 0.0d0
        auxR        = 0.0d0
        line        = ""
        keyword     = ""
        dataType    = -1
        eid         = 0
        row         = 0
        col         = 0
        val         = 0.0d0
        elid        = 0 

        call sroLog(" >> Sro::InpFile << | ctor")

        open(newunit=ifl, file=outdir // filename // ".info", status="old", action="read", iostat=ios)
        if (ios /= 0) then
            call sroXIT(" >>>>>> inpfile_ctor <<<<<< Error opening info file: " //&
                outdir // filename // ".info ")
        end if

        ! info file
        do
            read(ifl, '(a)', iostat=ios) lineRaw
            if (ios /= 0) exit
            line = adjustl(trim(lineRaw))

            if (line(1:1)  == "*") then

                auxId = 0
                commaIdx = index(line,",")
                if (commaIdx /= 0) then
                    keyword = trim(line(2:(commaIdx - 1)))
                else
                    keyword = trim(line(2:))
                end if
                write(kwd, '(a)') trim(upper(keyword))

                dtyp = findloc(keywords, kwd)
                dataType = dtyp(1)
            else 
                line = commaToSpace(line)
                select case(dataType)
                case(1) ! nnodel 

                    read(line, *) nnodel
                    allocate(con(nnodel(1)), source= 0)
                    allocate(adjcon(nnodel(2)), source= 0)
                    this % nndeluels = nnodel(1)
                    this % nndeladjuels = nnodel(2)
                    
                case(2) ! coords

                    if (auxId == 0) then
                        read(line, *) auxN 
                        allocate(this % x(auxN,3), source=0.0d0)
                        allocate(this % nds_global(auxN), source=0)
                        this % nnodes = auxN
                        auxId = auxId + 1
                    else
                        read(line, *) coords
                        this % x(auxId, :) = coords
                        auxId = auxId + 1
                    end if
                    
                case(3) ! uels

                    if (auxId == 0) then 
                        read(line, *) auxN
                        allocate(this % con(auxN, nnodel(1)), source = 0)
                        allocate(this % uels(auxN), source = 0)
                        allocate(this % uels_global(auxN), source = 0)
                        this % nuels = auxN
                        auxId = auxId + 1
                    else
                        read(line, *) this % uels(auxId)
                        auxId = auxId + 1
                    end if

                case(4) ! uelsg

                    if (auxId == 0) auxId = 1
                    read(line, *) this % uels_global(auxId)
                    auxId = auxId + 1

                case(5) ! con

                    read(line, *) con 
                    this % con(auxId+1, :) = con 
                    auxId = auxId + 1

                case(6) ! boundary

                    this % hasBoundary = .true.
                    if (auxId == 0) then 
                        read(line, *) auxN 
                        allocate(this % boundary(auxN), source=0)
                        auxId = auxId + 1
                    else 
                        read(line, *) this % boundary(auxId)
                        auxId = auxId + 1
                    end if

                case(7) ! props

                    if (auxId == 0) then 
                        read(line, *) auxN 
                        allocate(this % props(auxN), source=0.0d0)
                        auxId = auxId + 1
                    else 
                        read(line, *) this % props(auxId)
                        auxId = auxId + 1
                    end if

                case(8) ! orientations

                    if (auxId == 0) then 
                        read(line, *) auxN 
                        allocate(this % orientations(auxN), source=0.0d0)
                        auxId = auxId + 1
                    else 
                        read(line, *) this % orientations(auxId)
                        auxId = auxId + 1
                    end if

                case(9) ! thicknesses

                    if (auxId == 0) then 
                        read(line, *) auxN 
                        allocate(this % thicknesses(auxN), source=0.0d0)
                        auxId = auxId + 1
                    else 
                        read(line, *) this % thicknesses(auxId)
                        auxId = auxId + 1
                    end if

                case(10) ! transform_shell

                    read(line, *) coords 
                    this % transform_shell(auxId+1, :) = coords 
                    auxId = auxId + 1

                case(11) ! transform_ls

                    read(line, *) coords 
                    if (auxId+1 <= 3) then 
                        this % transform_ls(auxId+1, :) = coords 
                    else 
                        this % lsorigin = coords 
                    end if
                    auxId = auxId + 1

                case(12) ! ndsg

                    read(line, *) auxN 
                    this % nds_global(auxId+1) = auxN
                    auxId = auxId + 1

                case(13) ! nnodesuel

                    read(line, *) this % nnodesuel

                case(14) ! nnodesls

                    read(line, *) this % nnodesls
                
                case(15) ! lsnds

                    if (auxId == 0) then 
                        read(line, *) auxN
                        allocate(this % lsnodes(auxN), source=0)
                        allocate(this % lsinit(auxN), source=0.0d0)
                        auxId = auxId + 1
                    else 
                        read(line, *) this % lsnodes(auxId)
                        auxId = auxId + 1
                    end if
                
                case(16) ! lsinit

                    read(line, *) auxR 
                    this % lsinit(auxId + 1) = auxR 
                    auxId = auxId + 1

                case(17) ! meshsize

                    read(line,*) this % meshsize

                case(18) ! alpha

                    read(line,*) this % alpha

                case(19) ! vfrac

                    read(line,*) this % vfrac
                
                case(20) ! adjcoords

                    this % hasDisbond = .true.
                    if (auxId == 0) then
                        read(line, *) auxN 
                        allocate(this % adjx(auxN,3), source=0.0d0)
                        allocate(this % adjnds_global(auxN), source=0)
                        auxId = auxId + 1
                    else
                        read(line, *) coords
                        this % adjx(auxId, :) = coords
                        auxId = auxId + 1
                    end if

                case(21) ! adjuels

                    this % hasDisbond = .true.
                    if (auxId == 0) then 
                        read(line, *) auxN
                        allocate(this % adjcon(auxN, this % nndeladjuels), source = 0)
                        allocate(this % adjuels(auxN), source = 0)
                        allocate(this % adjuels_global(auxN), source = 0)
                        allocate(this % k(auxN, 24, 24), source=0.0d0)
                        allocate(this % dkda(auxN, 24, 24), source=0.0d0)
                        this % nadjuels = auxN
                        auxId = auxId + 1
                    else
                        read(line, *) this % adjuels(auxId)
                        auxId = auxId + 1
                    end if
                    
                case(22) ! adjuelsg

                    if (auxId == 0) auxId = 1
                    read(line, *) this % adjuels_global(auxId)
                    auxId = auxId + 1

                case(23) ! adjcon

                    read(line, *) adjcon 
                    this % adjcon(auxId+1, :) = adjcon 
                    auxId = auxId + 1

                case(24) ! adjndsg

                    read(line, *) auxN 
                    this % adjnds_global(auxId+1) = auxN
                    auxId = auxId + 1

                case(25) ! nnodesadjuel

                    read(line, *) this % nnodesadjuel

                end select
            end if
        end do

        if (this % hasDisbond) then
            open(newunit=sf1, file=outdir // filename // &
                "_stiff.mtx", status="old", action="read", iostat=ios)
            if (ios /= 0) then
                call sroXIT(" >>>>>> inpfile_ctor <<<<<< Error opening mtx file: " // &
                    outdir // filename(:(len_trim(filename) - 4)) // "_stiff.mtx ")
            end if

            open(newunit=sf2, file=outdir // filename // &
                "_stiff_da.mtx", status="old", action="read", iostat=ios)
            if (ios /= 0) then
                call sroXIT(" >>>>>> inpfile_ctor <<<<<< Error opening mtx file: " // &
                    outdir // filename(:(len_trim(filename) - 4)) // "_stiff_da.mtx ")
            end if

            ! sf1
            do
                read(sf1, '(a)', iostat=ios) lineRaw
                if (ios /= 0) exit
                line = adjustl(trim(lineRaw))

                read(line, *) eid, row, col, val 

                elid = findloc(this % adjuels_global, eid)
                if (elid(1) /= 0) then
                    this % k(elid(1), row, col) = val
                    this % k(elid(1), col, row) = val
                else 
                    call sroXIT(" >>>>>> inpfile_ctor <<<<<< " &
                    // "Unable to match element id while reading stiffness data")
                end if

            end do

            ! sf2
            do
                read(sf2, '(a)', iostat=ios) lineRaw
                if (ios /= 0) exit
                line = adjustl(trim(lineRaw))

                read(line, *) eid, row, col, val 

                elid = findloc(this % adjuels_global, eid)
                if (elid(1) /= 0) then
                    this % dkda(elid(1), row, col) = val
                    this % dkda(elid(1), col, row) = val
                else 
                    call sroXIT(" >>>>>> inpfile_ctor <<<<<< " &
                    // "Unable to match element id while reading stiffness (da) data")
                end if

            end do

            close(sf1)
            close(sf2)
        end if

        close(ifl)
        if (allocated(con)) deallocate(con)
        if (allocated(adjcon)) deallocate(adjcon)
    end function

    subroutine inpfile_dtor(this)
        implicit none 

        type(INPFile), intent(inout) :: this 

        call sroLog(" >> Sro::InpFile << | dtor")

        call this % finalize
    end subroutine

    subroutine inpfile_destroy(this)
        implicit none 

        class(INPFile), intent(inout) :: this 

        call sroLog(" >> Sro::InpFile << | destroy")

        if (allocated(this % x)) deallocate(this % x)
        if (allocated(this % con)) deallocate(this % con)
        if (allocated(this % nds_global)) deallocate(this % nds_global)
        if (allocated(this % uels)) deallocate(this % uels)
        if (allocated(this % uels_global)) deallocate(this % uels_global)

        if (allocated(this % boundary)) deallocate(this % boundary)
        
        if (allocated(this % orientations)) deallocate(this % orientations)
        if (allocated(this % thicknesses)) deallocate(this % thicknesses)

        if (allocated(this % lsnodes)) deallocate(this % lsnodes)
        if (allocated(this % lsinit)) deallocate(this % lsinit)

        if (allocated(this % adjx)) deallocate(this % adjx)
        if (allocated(this % adjcon)) deallocate(this % adjcon)
        if (allocated(this % adjnds_global)) deallocate(this % adjnds_global)
        if (allocated(this % adjuels)) deallocate(this % adjuels)
        if (allocated(this % adjuels_global)) deallocate(this % adjuels_global)

        if (allocated(this % k)) deallocate(this % k)
        if (allocated(this % dkda)) deallocate(this % dkda)

    end subroutine

end module