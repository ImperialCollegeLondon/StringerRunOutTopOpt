!DEC$ FREEFORM

module stringUtils_mod
    use iso_fortran_env
    use abaqusXIT_mod
    implicit none 

    public :: upper
    public :: commaToSpace
    character(len=26), parameter :: uca = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=26), parameter :: lca = 'abcdefghijklmnopqrstuvwxyz'

contains 

    pure character function upperCh(ch)
        implicit none
        character, intent(in) :: ch 
        integer :: ix 

        ! vars init
        ix = 0
        upperCh = 'a'

        upperCh = ch 
        ix = index(lca, ch)
        if (ix /= 0) upperCh = uca(ix:ix)
    end function

    pure function upper(str) result(ustr)
        implicit none 
        character(len=*), intent(in) :: str 
        character(len=len(str)) :: ustr 
        integer :: idx

        ! vars init
        idx = 0
        ustr = str

        do idx = 1, len(str)
            ustr(idx:idx) = upperCh(str(idx:idx))
        end do
    end function

    pure function commaToSpace(str) result(rstr)
        implicit none 
        character(len=*), intent(in) :: str 
        character(len=len(str)) :: rstr 
        integer :: ii 

        ! vars init
        ii = 0
        rstr = str

        do ii = 1, len(str)
            if (str(ii:ii) == ",") rstr(ii:ii) = " "
        end do
    end function

    pure function int2str(int) result(str)
        implicit none
        character(len=:), allocatable :: str 
        integer, intent(in) :: int 

        character(len=50) :: aux

        ! vars init
        str = ""

        write(aux, *) int 
        str = trim(adjustl(aux))
    end function

    pure function float2str(flt) result(str)
        implicit none
        character(len=:), allocatable :: str 
        real(real64), intent(in) :: flt 

        character(len=50) :: aux

        ! vars init
        str = ""

        ! write(aux, *) flt
        if (abs(flt) < 1.0E-9) then
            write(aux, '(a)') '0.0'
        else
            ! clamps values between -1E20 and 1E20
            write(aux, '(ES13.5)') max( min(flt, 1.0d20) , -1.0d20) 
        end if
        str = trim(adjustl(aux))
    end function

    pure function intArray2str(int) result(str)
        implicit none
        character(len=:), allocatable :: str 
        integer, intent(in) :: int(:)

        integer :: i 

        !vars init
        str = ""
        i = 0

        do i = 1, size(int)
            str = str // int2str(int(i)) // ' '
        end do
        str = trim(adjustl(str))
    end function

    pure function floatArray2str(flt) result(str)
        implicit none
        character(len=:), allocatable :: str 
        real(real64), intent(in) :: flt(:)

        integer :: i 

        ! vars init
        str = ""
        i = 0
        
        do i = 1, size(flt)
            str = str // float2str(flt(i)) // ' '
        end do
        str = trim(adjustl(str))
    end function

    pure function logical2str(bool) result(str)
        implicit none 
        character(len=:), allocatable :: str
        logical, intent(in) :: bool

        str = "FALSE"
        if (bool) str = "TRUE"
    end function

end module