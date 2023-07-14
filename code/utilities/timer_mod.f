!DEC$ FREEFORM

module timer_mod
    use iso_fortran_env
    implicit none
    
    type :: Timer 
        integer :: c0
        integer :: rate

    contains
        procedure :: getTime    => timer_getTime
        procedure :: getTimeStr => timer_getTimeStr

    end type

    interface Timer
        procedure :: timer_ctor
    end interface
    
contains

    function timer_ctor() result(this)
        implicit none

        type(Timer) :: this

        call system_clock(this % c0, this % rate)
    end function

    function timer_getTime(this) result(t)
        implicit none

        class(Timer) :: this
        real(real64) :: t 

        integer :: c

        t = 0.0d0
        c = 0

        call system_clock(c)
        t = (real(c) - real(this % c0)) / this % rate 
    end function

    function timer_getTimeStr(this, tin) result(tstr)
        implicit none

        class(Timer) :: this
        character(len=:), allocatable :: tstr
        real(real64), optional, intent(in) :: tin

        real(real64) :: t
        integer :: h
        integer :: m
        real(real64) :: s

        character(len=50) :: aux
        character(len=20) :: ch
        character(len=20) :: cm
        character(len=20) :: cs

        tstr = ""
        h = 0
        m = 0
        s = 0.0d0
        
        t = this % getTime()
        if (present(tin)) t = tin

        s = t
        h = floor(s) / 3600
        s = mod(s, 3600.0d0)
        m = floor(s) / 60
        s = mod(s, 60.0d0)
        
        write(ch, "(i4)") h
        write(cm, "(i2)") m
        write(cs, "(f3.3)") s
        write(aux, "(i4, a, i3, a, f7.3, a)") h, "h", m, "m", s, "s"
        
        tstr = trim(aux)
    end function
   
end module