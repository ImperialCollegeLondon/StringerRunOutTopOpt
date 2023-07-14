!DEC$ FREEFORM

module abaqusXIT_mod
    use iso_fortran_env
    use timer_mod

    logical :: use_separate_file = .false.
    integer :: errunit = 6

    type(Timer) :: logtimer

contains 

#ifdef SRO_DEBUG
    subroutine XIT()
        error stop 
    end subroutine
#endif

    subroutine sroXIT(msg)
        character(len=*) :: msg
#ifdef SRO_DEBUG
        print *, msg
        print *, " >>>>>> sroXIT <<<<<< Terminating due to error"
#else
        dimension intv(1), realv(1)
        character*8 charv(1)
        
        if (use_separate_file) then 
            write(errunit, '(a)') logtimer % getTimeStr() // " --[ERROR]-- " // msg
        end if

        call stdb_abqerr(-3, logtimer % getTimeStr() // msg, intv, realv, charv)
        call stdb_abqerr(-3, logtimer % getTimeStr() // &
            " >>>>>> sroXIT <<<<<< Terminating due to error", intv, realv, charv)
#endif
        call XIT
    end subroutine

    subroutine sroLog(msg)
        character(len=*) :: msg

#ifdef SRO_DEBUG
        print *, msg 
#else 
        dimension intv(1), realv(1)
        character*8 charv(1)
        
        if (use_separate_file) then
            write(errunit, '(a)') logtimer % getTimeStr() // msg
        else
            call stdb_abqerr(1, logtimer % getTimeStr() // msg, intv, realv, charv)
        end if
#endif
    end subroutine
    
    subroutine sroSetLogFile(file_unit)
        integer :: file_unit 

        errunit = file_unit 
        use_separate_file = .true.
    end subroutine

    subroutine sroInitTimer()
        logtimer = Timer()
    end subroutine

end module