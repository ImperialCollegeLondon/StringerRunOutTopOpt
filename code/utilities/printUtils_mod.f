!DEC$ FREEFORM

module printUtils_mod
    use iso_fortran_env
    use abaqusXIT_mod
    use stringUtils_mod
    implicit none 

    character(len=50), parameter :: sep1   = "--------------------------------------------------" 
    character(len=50), parameter :: sep2   = "==================================================" 
    character(len=50), parameter :: sep3   = "//////////////////////////////////////////////////" 
    
    character(len= 8), parameter :: fmt001 = "(1x, i6)" 
    character(len=11), parameter :: fmt002 = "(1x, e10.3)" 
    character(len= 8), parameter :: fmt003 = "(1x, l6)" 
    
    character(len= 3), parameter :: fmt100 = "(a)" 
    character(len=11), parameter :: fmt101 = "(a, 1x, i6)" 
    character(len=14), parameter :: fmt102 = "(a, 1x, e10.3)" 
    character(len=11), parameter :: fmt103 = "(a, 1x, l6)" 
    character(len=12), parameter :: fmt104 = "(a, 1x, i12)"
    
    character(len=19), parameter :: fmt201 = "(1x,  1(e10.3, 2x))"
    character(len=19), parameter :: fmt202 = "(1x,  2(e10.3, 2x))"
    character(len=19), parameter :: fmt203 = "(1x,  3(e10.3, 2x))"
    character(len=19), parameter :: fmt204 = "(1x,  4(e10.3, 2x))"
    character(len=19), parameter :: fmt205 = "(1x,  5(e10.3, 2x))"
    character(len=19), parameter :: fmt206 = "(1x,  6(e10.3, 2x))"
    character(len=19), parameter :: fmt207 = "(1x,  7(e10.3, 2x))"
    character(len=19), parameter :: fmt208 = "(1x,  8(e10.3, 2x))"
    character(len=19), parameter :: fmt209 = "(1x,  9(e10.3, 2x))"
    character(len=19), parameter :: fmt210 = "(1x, 10(e10.3, 2x))"
    character(len=19), parameter :: fmt211 = "(1x, 11(e10.3, 2x))"
    character(len=19), parameter :: fmt212 = "(1x, 12(e10.3, 2x))"
    character(len=19), parameter :: fmt213 = "(1x, 13(e10.3, 2x))"
    character(len=19), parameter :: fmt214 = "(1x, 14(e10.3, 2x))"
    character(len=19), parameter :: fmt215 = "(1x, 15(e10.3, 2x))"
    character(len=19), parameter :: fmt216 = "(1x, 16(e10.3, 2x))"
    character(len=19), parameter :: fmt217 = "(1x, 17(e10.3, 2x))"
    character(len=19), parameter :: fmt218 = "(1x, 18(e10.3, 2x))"
    character(len=19), parameter :: fmt219 = "(1x, 19(e10.3, 2x))"
    character(len=19), parameter :: fmt220 = "(1x, 20(e10.3, 2x))"
    character(len=19), parameter :: fmt221 = "(1x, 21(e10.3, 2x))"
    character(len=19), parameter :: fmt222 = "(1x, 22(e10.3, 2x))"
    character(len=19), parameter :: fmt223 = "(1x, 23(e10.3, 2x))"
    character(len=19), parameter :: fmt224 = "(1x, 24(e10.3, 2x))"
    character(len=19), parameter :: fmt225 = "(1x, 25(e10.3, 2x))"
    character(len=19), parameter :: fmt226 = "(1x, 26(e10.3, 2x))"
    character(len=19), parameter :: fmt227 = "(1x, 27(e10.3, 2x))"
    character(len=19), parameter :: fmt228 = "(1x, 28(e10.3, 2x))"
    character(len=19), parameter :: fmt229 = "(1x, 29(e10.3, 2x))"
    character(len=19), parameter :: fmt230 = "(1x, 30(e10.3, 2x))"
    
    character(len=25), parameter :: fmt301 = "(a, /, 1x,  1(e10.3, 2x))"
    character(len=25), parameter :: fmt302 = "(a, /, 1x,  2(e10.3, 2x))"
    character(len=25), parameter :: fmt303 = "(a, /, 1x,  3(e10.3, 2x))"
    character(len=25), parameter :: fmt304 = "(a, /, 1x,  4(e10.3, 2x))"
    character(len=25), parameter :: fmt305 = "(a, /, 1x,  5(e10.3, 2x))"
    character(len=25), parameter :: fmt306 = "(a, /, 1x,  6(e10.3, 2x))"
    character(len=25), parameter :: fmt307 = "(a, /, 1x,  7(e10.3, 2x))"
    character(len=25), parameter :: fmt308 = "(a, /, 1x,  8(e10.3, 2x))"
    character(len=25), parameter :: fmt309 = "(a, /, 1x,  9(e10.3, 2x))"
    character(len=25), parameter :: fmt310 = "(a, /, 1x, 10(e10.3, 2x))"
    character(len=25), parameter :: fmt311 = "(a, /, 1x, 11(e10.3, 2x))"
    character(len=25), parameter :: fmt312 = "(a, /, 1x, 12(e10.3, 2x))"
    character(len=25), parameter :: fmt313 = "(a, /, 1x, 13(e10.3, 2x))"
    character(len=25), parameter :: fmt314 = "(a, /, 1x, 14(e10.3, 2x))"
    character(len=25), parameter :: fmt315 = "(a, /, 1x, 15(e10.3, 2x))"
    character(len=25), parameter :: fmt316 = "(a, /, 1x, 16(e10.3, 2x))"
    character(len=25), parameter :: fmt317 = "(a, /, 1x, 17(e10.3, 2x))"
    character(len=25), parameter :: fmt318 = "(a, /, 1x, 18(e10.3, 2x))"
    character(len=25), parameter :: fmt319 = "(a, /, 1x, 19(e10.3, 2x))"
    character(len=25), parameter :: fmt320 = "(a, /, 1x, 20(e10.3, 2x))"
    character(len=25), parameter :: fmt321 = "(a, /, 1x, 21(e10.3, 2x))"
    character(len=25), parameter :: fmt322 = "(a, /, 1x, 22(e10.3, 2x))"
    character(len=25), parameter :: fmt323 = "(a, /, 1x, 23(e10.3, 2x))"
    character(len=25), parameter :: fmt324 = "(a, /, 1x, 24(e10.3, 2x))"
    character(len=25), parameter :: fmt325 = "(a, /, 1x, 25(e10.3, 2x))"
    character(len=25), parameter :: fmt326 = "(a, /, 1x, 26(e10.3, 2x))"
    character(len=25), parameter :: fmt327 = "(a, /, 1x, 27(e10.3, 2x))"
    character(len=25), parameter :: fmt328 = "(a, /, 1x, 28(e10.3, 2x))"
    character(len=25), parameter :: fmt329 = "(a, /, 1x, 29(e10.3, 2x))"
    character(len=25), parameter :: fmt330 = "(a, /, 1x, 30(e10.3, 2x))" 

    character(len=16), parameter :: fmt401 = "(1x,  1(i6, 2x))"
    character(len=16), parameter :: fmt402 = "(1x,  2(i6, 2x))"
    character(len=16), parameter :: fmt403 = "(1x,  3(i6, 2x))"
    character(len=16), parameter :: fmt404 = "(1x,  4(i6, 2x))"
    character(len=16), parameter :: fmt405 = "(1x,  5(i6, 2x))"
    character(len=16), parameter :: fmt406 = "(1x,  6(i6, 2x))"
    character(len=16), parameter :: fmt407 = "(1x,  7(i6, 2x))"
    character(len=16), parameter :: fmt408 = "(1x,  8(i6, 2x))"
    character(len=16), parameter :: fmt409 = "(1x,  9(i6, 2x))"
    character(len=16), parameter :: fmt410 = "(1x, 10(i6, 2x))"
    character(len=16), parameter :: fmt411 = "(1x, 11(i6, 2x))"
    character(len=16), parameter :: fmt412 = "(1x, 12(i6, 2x))"
    character(len=16), parameter :: fmt413 = "(1x, 13(i6, 2x))"
    character(len=16), parameter :: fmt414 = "(1x, 14(i6, 2x))"
    character(len=16), parameter :: fmt415 = "(1x, 15(i6, 2x))"
    character(len=16), parameter :: fmt416 = "(1x, 16(i6, 2x))"
    character(len=16), parameter :: fmt417 = "(1x, 17(i6, 2x))"
    character(len=16), parameter :: fmt418 = "(1x, 18(i6, 2x))"
    character(len=16), parameter :: fmt419 = "(1x, 19(i6, 2x))"
    character(len=16), parameter :: fmt420 = "(1x, 20(i6, 2x))"
    character(len=16), parameter :: fmt421 = "(1x, 21(i6, 2x))"
    character(len=16), parameter :: fmt422 = "(1x, 22(i6, 2x))"
    character(len=16), parameter :: fmt423 = "(1x, 23(i6, 2x))"
    character(len=16), parameter :: fmt424 = "(1x, 24(i6, 2x))"
    character(len=16), parameter :: fmt425 = "(1x, 25(i6, 2x))"
    character(len=16), parameter :: fmt426 = "(1x, 26(i6, 2x))"
    character(len=16), parameter :: fmt427 = "(1x, 27(i6, 2x))"
    character(len=16), parameter :: fmt428 = "(1x, 28(i6, 2x))"
    character(len=16), parameter :: fmt429 = "(1x, 29(i6, 2x))"
    character(len=16), parameter :: fmt430 = "(1x, 30(i6, 2x))"
    
    character(len=22), parameter :: fmt501 = "(a, /, 1x,  1(i6, 2x))"
    character(len=22), parameter :: fmt502 = "(a, /, 1x,  2(i6, 2x))"
    character(len=22), parameter :: fmt503 = "(a, /, 1x,  3(i6, 2x))"
    character(len=22), parameter :: fmt504 = "(a, /, 1x,  4(i6, 2x))"
    character(len=22), parameter :: fmt505 = "(a, /, 1x,  5(i6, 2x))"
    character(len=22), parameter :: fmt506 = "(a, /, 1x,  6(i6, 2x))"
    character(len=22), parameter :: fmt507 = "(a, /, 1x,  7(i6, 2x))"
    character(len=22), parameter :: fmt508 = "(a, /, 1x,  8(i6, 2x))"
    character(len=22), parameter :: fmt509 = "(a, /, 1x,  9(i6, 2x))"
    character(len=22), parameter :: fmt510 = "(a, /, 1x, 10(i6, 2x))"
    character(len=22), parameter :: fmt511 = "(a, /, 1x, 11(i6, 2x))"
    character(len=22), parameter :: fmt512 = "(a, /, 1x, 12(i6, 2x))"
    character(len=22), parameter :: fmt513 = "(a, /, 1x, 13(i6, 2x))"
    character(len=22), parameter :: fmt514 = "(a, /, 1x, 14(i6, 2x))"
    character(len=22), parameter :: fmt515 = "(a, /, 1x, 15(i6, 2x))"
    character(len=22), parameter :: fmt516 = "(a, /, 1x, 16(i6, 2x))"
    character(len=22), parameter :: fmt517 = "(a, /, 1x, 17(i6, 2x))"
    character(len=22), parameter :: fmt518 = "(a, /, 1x, 18(i6, 2x))"
    character(len=22), parameter :: fmt519 = "(a, /, 1x, 19(i6, 2x))"
    character(len=22), parameter :: fmt520 = "(a, /, 1x, 20(i6, 2x))"
    character(len=22), parameter :: fmt521 = "(a, /, 1x, 21(i6, 2x))"
    character(len=22), parameter :: fmt522 = "(a, /, 1x, 22(i6, 2x))"
    character(len=22), parameter :: fmt523 = "(a, /, 1x, 23(i6, 2x))"
    character(len=22), parameter :: fmt524 = "(a, /, 1x, 24(i6, 2x))"
    character(len=22), parameter :: fmt525 = "(a, /, 1x, 25(i6, 2x))"
    character(len=22), parameter :: fmt526 = "(a, /, 1x, 26(i6, 2x))"
    character(len=22), parameter :: fmt527 = "(a, /, 1x, 27(i6, 2x))"
    character(len=22), parameter :: fmt528 = "(a, /, 1x, 28(i6, 2x))"
    character(len=22), parameter :: fmt529 = "(a, /, 1x, 29(i6, 2x))"
    character(len=22), parameter :: fmt530 = "(a, /, 1x, 30(i6, 2x))" 

    interface printVar
        procedure :: printMatR, printMatI, printVecR, printVecI
    end interface

contains

    subroutine printMatR(v, name, unit)
        implicit none
        real(real64) :: v(:,:)
        character(len=*), intent(in) :: name
        integer, intent(in) :: unit

        integer :: s(2), ii, jj
        character(len=5) :: cs1, cs2
        character(len=:), allocatable :: fmt

        ! vars init
        s   = 0
        ii  = 0
        jj  = 0
        cs1 = "   "
        cs2 = "   "
        fmt = ""

        s = shape(v)
        write(cs1, "(i5)") s(1)
        write(cs2, "(i5)") s(2)
        write(unit, fmt100) name // " [" // cs1 // ", " //cs2// "] <Real64>"

        select case(s(2))
        case(1)
            fmt = fmt201
        case(2)
            fmt = fmt202
        case(3)
            fmt = fmt203
        case(4)
            fmt = fmt204
        case(5)
            fmt = fmt205
        case(6)
            fmt = fmt206
        case(7)
            fmt = fmt207
        case(8)
            fmt = fmt208
        case(9)
            fmt = fmt209
        case(10)
            fmt = fmt210
        case(11)
            fmt = fmt211
        case(12)
            fmt = fmt212
        case(13)
            fmt = fmt213
        case(14)
            fmt = fmt214
        case(15)
            fmt = fmt215
        case(16)
            fmt = fmt216
        case(17)
            fmt = fmt217
        case(18)
            fmt = fmt218
        case(19)
            fmt = fmt219
        case(20)
            fmt = fmt220
        case(21)
            fmt = fmt221
        case(22)
            fmt = fmt222
        case(23)
            fmt = fmt223
        case(24)
            fmt = fmt224
        case(25)
            fmt = fmt225
        case(26)
            fmt = fmt226
        case(27)
            fmt = fmt227
        case(28)
            fmt = fmt228
        case(29)
            fmt = fmt229
        case(30)
            fmt = fmt230
        end select
        if (s(2) <= 30) then
            do ii = 1, s(1)
                write(unit, fmt) v(ii,:)
            end do
        else
            fmt = "(1x, " // int2str(s(2)) // "(e10.3, 2x))"
            do ii = 1, s(1)
                write(unit, fmt) v(ii,:)
            end do
        end if
        write(unit, fmt100) ""
    end subroutine

    subroutine printMatI(v, name, unit)
        implicit none
        integer :: v(:,:)
        character(len=*), intent(in) :: name
        integer, intent(in) :: unit

        integer :: s(2), ii, jj
        character(len=5) :: cs1, cs2
        character(len=:), allocatable :: fmt

        ! vars init
        s   = 0
        ii  = 0
        jj  = 0
        cs1 = "   "
        cs2 = "   "
        fmt = ""

        s = shape(v)
        write(cs1, "(i5)") s(1)
        write(cs2, "(i5)") s(2)
        write(unit, fmt100) name // " [" // cs1 // ", " //cs2// "] <Integer>"

        select case(s(2))
        case(1)
            fmt = fmt401
        case(2)
            fmt = fmt402
        case(3)
            fmt = fmt403
        case(4)
            fmt = fmt404
        case(5)
            fmt = fmt405
        case(6)
            fmt = fmt406
        case(7)
            fmt = fmt407
        case(8)
            fmt = fmt408
        case(9)
            fmt = fmt409
        case(10)
            fmt = fmt410
        case(11)
            fmt = fmt411
        case(12)
            fmt = fmt412
        case(13)
            fmt = fmt413
        case(14)
            fmt = fmt414
        case(15)
            fmt = fmt415
        case(16)
            fmt = fmt416
        case(17)
            fmt = fmt417
        case(18)
            fmt = fmt418
        case(19)
            fmt = fmt419
        case(20)
            fmt = fmt420
        case(21)
            fmt = fmt421
        case(22)
            fmt = fmt422
        case(23)
            fmt = fmt423
        case(24)
            fmt = fmt424
        case(25)
            fmt = fmt425
        case(26)
            fmt = fmt426
        case(27)
            fmt = fmt427
        case(28)
            fmt = fmt428
        case(29)
            fmt = fmt429
        case(30)
            fmt = fmt430
        end select
        if (s(2) <= 30) then
            do ii = 1, s(1)
                write(unit, fmt) v(ii,:)
            end do
        else
            fmt = "(1x, " // int2str(s(2)) // "(i6, 2x))"
            do ii = 1, s(1)
                write(unit, fmt) v(ii,:)
            end do
        end if
        write(unit, fmt100) ""
    end subroutine

    subroutine printVecR(v, name, unit)
        implicit none
        real(real64) :: v(:)
        character(len=*), intent(in) :: name
        integer, intent(in) :: unit

        integer :: s
        character(len=5) :: cs1
        character(len=:), allocatable :: fmt

        ! vars init
        s   = 0
        cs1 = "   "
        fmt = ""

        s = size(v)
        write(cs1, "(i5)") s
        write(unit, fmt100) name // " [" // cs1 // "] <Real64>"

        select case(s)
        case(1)
            fmt = fmt201
        case(2)
            fmt = fmt202
        case(3)
            fmt = fmt203
        case(4)
            fmt = fmt204
        case(5)
            fmt = fmt205
        case(6)
            fmt = fmt206
        case(7)
            fmt = fmt207
        case(8)
            fmt = fmt208
        case(9)
            fmt = fmt209
        case(10)
            fmt = fmt210
        case(11)
            fmt = fmt211
        case(12)
            fmt = fmt212
        case(13)
            fmt = fmt213
        case(14)
            fmt = fmt214
        case(15)
            fmt = fmt215
        case(16)
            fmt = fmt216
        case(17)
            fmt = fmt217
        case(18)
            fmt = fmt218
        case(19)
            fmt = fmt219
        case(20)
            fmt = fmt220
        case(21)
            fmt = fmt221
        case(22)
            fmt = fmt222
        case(23)
            fmt = fmt223
        case(24)
            fmt = fmt224
        case(25)
            fmt = fmt225
        case(26)
            fmt = fmt226
        case(27)
            fmt = fmt227
        case(28)
            fmt = fmt228
        case(29)
            fmt = fmt229
        case(30)
            fmt = fmt230
        end select
        if (s <= 30) then
            write(unit, fmt) v
        else
            fmt = "(1x, " // int2str(s) // "(e10.3, 2x))"
            write(unit, fmt) v
        end if
        write(unit, fmt100) ""
    end subroutine

    subroutine printVecI(v, name, unit)
        implicit none
        integer :: v(:)
        character(len=*), intent(in) :: name
        integer, intent(in) :: unit

        integer :: s
        character(len=5) :: cs1
        character(len=:), allocatable :: fmt

        ! vars init
        s   = 0
        cs1 = "   "
        fmt = ""

        s = size(v)
        write(cs1, "(i5)") s
        write(unit, fmt100) name // " [" // cs1 // "] <Integer>"

        select case(s)
        case(1)
            fmt = fmt401
        case(2)
            fmt = fmt402
        case(3)
            fmt = fmt403
        case(4)
            fmt = fmt404
        case(5)
            fmt = fmt405
        case(6)
            fmt = fmt406
        case(7)
            fmt = fmt407
        case(8)
            fmt = fmt408
        case(9)
            fmt = fmt409
        case(10)
            fmt = fmt410
        case(11)
            fmt = fmt411
        case(12)
            fmt = fmt412
        case(13)
            fmt = fmt413
        case(14)
            fmt = fmt414
        case(15)
            fmt = fmt415
        case(16)
            fmt = fmt416
        case(17)
            fmt = fmt417
        case(18)
            fmt = fmt418
        case(19)
            fmt = fmt419
        case(20)
            fmt = fmt420
        case(21)
            fmt = fmt421
        case(22)
            fmt = fmt422
        case(23)
            fmt = fmt423
        case(24)
            fmt = fmt424
        case(25)
            fmt = fmt425
        case(26)
            fmt = fmt426
        case(27)
            fmt = fmt427
        case(28)
            fmt = fmt428
        case(29)
            fmt = fmt429
        case(30)
            fmt = fmt430
        end select
        if (s <= 30) then
            write(unit, fmt) v
        else
            fmt = "(1x, " // int2str(s) // "(i6, 2x))"
            write(unit, fmt) v
        end if
        write(unit, fmt100) ""
    end subroutine

end module