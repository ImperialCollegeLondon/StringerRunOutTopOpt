!DEC$ FREEFORM

module levelsetOptimization_mod
    use iso_fortran_env
    use component_mod
    implicit none
    
    type, extends(Component) :: LSOptimizer
        
        real(real64) :: lagrMult    = 0.0d0
        real(real64) :: penalty     = 0.0d0
        real(real64) :: lambda      = 0.0d0
        real(real64) :: lambdaInit  = 0.0d0
        real(real64) :: alpha       = 0.0d0
        real(real64) :: beta        = 0.0d0

        real(real64) :: comp        = 0.0d0
        real(real64) :: enrr        = 0.0d0
        real(real64) :: lagr        = 0.0d0

        real(real64) :: vol         = 0.0d0
        real(real64) :: volInit     = 0.0d0
        real(real64) :: volFrac     = 0.0d0
        real(real64) :: currentVolFrac  = 0.0d0

        integer      :: convergenceCounter = 0
        integer      :: iterationCounter = 0

        logical      :: initSet     = .false.

    contains
        
        procedure   :: update   => lsopt_update
        procedure   :: print    => lsopt_print

    end type
    
    interface LSOptimizer
        procedure :: lsopt_ctor
    end interface
    
contains

    function lsopt_ctor(alpha, vfrac, vol) result(this) 
        implicit none 

        type(LSOptimizer) :: this 
        real(real64), intent(in) :: alpha
        real(real64), intent(in) :: vfrac
        real(real64), intent(in) :: vol

        this % name         = "LSOptimizer"

        this % alpha        = alpha 
        this % volFrac      = vfrac

        this % lagrMult     = 1.0d0
        this % penalty      = 1.0d0
        this % lambda       = 0.0d0
        this % comp         = 0.0d0
        this % enrr         = 0.0d0
        this % lagr         = 0.0d0
        this % vol          = 0.0d0
        this % volInit      = vol
        this % currentVolFrac = 1.0d0
        this % beta         = 1.05d0
        this % lambdaInit   = 0.0d0
        this % convergenceCounter = 0
        this % iterationCounter = 1

        this % initSet      = .false.

    end function

    subroutine lsopt_update(this)
        implicit none
    
        class(LSOptimizer), intent(inout) :: this

        real(real64) :: vmax, dv, dvfrac, aux

        ! vars init
        vmax = 0.0d0
        dv   = 0.0d0
        dvfrac = 0.0d0
        aux = 0.0d0
    
        vmax = this % volFrac * this % volInit
        dv = this % vol - vmax
        dv = 2.0d0 * (dv / this % volInit)
        
#ifndef SRO_CONSTRAINT_STAB
        this % lagrMult = max(this % lagrMult + this % penalty * dv, 0.0d0)
        this % penalty = this % beta * this % penalty
        this % lambda = max(this % lagrMult + this % penalty * dv, 0.0d0)
#else
        this % lagrMult = min(max(this % lagrMult + this % penalty * dv, 0.0d0), 5.0d0)
        this % penalty = min(this % beta * this % penalty, 10.0d0)
        this % lambda = max(this % lagrMult + this % penalty * dv, 0.0d0)
#endif

        this % iterationCounter = this % iterationCounter + 1
        
        this % lagr =   (this % alpha * this % enrr)        + &
                        ((1 - this % alpha) * this % comp)  + &
                        (this % lagrMult * dv)              + &
                        (0.5d0 * this % penalty * dv * dv)
        
        this % currentVolFrac = this % vol / this % volInit

        if (.not. this % initSet) then
            this % lambdaInit = this % lambda
            this % initSet = .true.
        end if
        
    end subroutine

    subroutine lsopt_print(this, varName, unit)
        implicit none
        
        class(LSOptimizer), intent(inout) :: this
        character(len=*), intent(in) :: varName
        integer, optional, intent(in) :: unit
        
        integer :: uid
    
        ! vars init 
        uid = 0
        
        uid = 6
        if (present(unit)) uid = unit 
        
        write(uid, fmt100) ""
        write(uid, fmt100) varName // " <" // this % name // ">"
        write(uid, fmt100) ""
        write(uid, fmt102) "lagrMult ", this % lagrMult
        write(uid, fmt102) "penalty  ", this % lagrMult
        write(uid, fmt102) "lambda   ", this % lagrMult
        write(uid, fmt102) "alpha    ", this % lagrMult
        write(uid, fmt102) "comp ", this % lagrMult
        write(uid, fmt102) "enrr ", this % lagrMult
        write(uid, fmt102) "lagr ", this % lagrMult
        write(uid, fmt102) "vol     ", this % lagrMult
        write(uid, fmt102) "volInit ", this % lagrMult
        write(uid, fmt102) "volFrac ", this % lagrMult

    end subroutine
   
end module