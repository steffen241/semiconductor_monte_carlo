module integrals
  use physconst
  use fgsl
  use, intrinsic :: iso_c_binding
  implicit none

  contains

  function int_G1(x) result(res)
        real( kind = fgsl_double )   :: res, error
        real(kind=s)                 :: x
        integer( kind = fgsl_size_t) :: neval
        integer( kind = fgsl_int)    :: i
        type( fgsl_function )        :: func

        func = fgsl_function_init( G1, c_null_ptr )

        i = fgsl_integration_qng ( func,              &
                                 -0.0000000001_fgsl_double,       &
                                 x,       &
                                 1e-9_fgsl_double,      &
                                 1e-9_fgsl_double,      &
                                 res, error, neval )
    end function int_G1

    function G1( x, params ) bind(c)
        implicit none

        real( kind = c_double )        :: G1
        real( kind = c_double ), value :: x
        type( c_ptr ), value           :: params

        params = params
        G1 = ((1.0_s/(exp(x)-1.0_s))+1.0_s)*x**2.0_s
    end function G1

    function int_G2(x) result(res)
        real( kind = fgsl_double )   :: res, error
        real(kind=s)                 :: x
        integer( kind = fgsl_size_t) :: neval
        integer( kind = fgsl_int)    :: i
        type( fgsl_function )        :: func

        func = fgsl_function_init( G2, c_null_ptr )

        i = fgsl_integration_qng ( func,              &
                                 -0.0000000001_fgsl_double,       &
                                 x,       &
                                 1e-9_fgsl_double,      &
                                 1e-9_fgsl_double,      &
                                 res, error, neval )
    end function int_G2

    function G2( x, params ) bind(c)
        implicit none

        real( kind = c_double )        :: G2
        real( kind = c_double ), value :: x
        type( c_ptr ), value           :: params

        params = params
        G2 = ((1.0_s/(exp(x)-1.0_s))+1.0_s)*x**3.0_s
    end function G2

!    function int_F1(x) result(res)
!        real( kind = fgsl_double )   :: res, error
!        real(kind= fgsl_double)                 :: x
!        integer( kind = fgsl_size_t) :: nmax
!        integer( kind = fgsl_int)    :: i
!        type( fgsl_function )        :: func
!        type( fgsl_integration_workspace) :: wk
!
!        func = fgsl_function_init( F1, c_null_ptr )
!        wk = fgsl_integration_workspace_alloc(1000_fgsl_size_t)
!        nmax = 100
!
!        i = fgsl_integration_qags ( func,              &
!                                 -0.0_fgsl_double,       &
!                                 x,       &
!                                 1e-9_fgsl_double,      &
!                                 1e-9_fgsl_double,      &
!                                 nmax, wk, res, error )
!
!    end function int_F1

    function int_F1(x) result(res)
        real( kind = fgsl_double )   :: res, error
        real(kind=s)                 :: x
        integer( kind = fgsl_size_t) :: neval
        integer( kind = fgsl_int)    :: i
        type( fgsl_function )        :: func

        func = fgsl_function_init( F_1, c_null_ptr )

        i = fgsl_integration_qng ( func,              &
                                 -0.0000000001_fgsl_double,       &
                                 x,       &
                                 1e-5_fgsl_double,      &
                                 1e-5_fgsl_double,      &
                                 res, error, neval )
    end function int_F1

    function F_1( x, params ) bind(c)
        implicit none

        real( kind = c_double )        :: F_1
        real( kind = c_double ), value :: x
        type( c_ptr ), value           :: params

        params = params
        F_1 = (1.0_s/(exp(x)-1.0_s))*x**2.0_s
    end function F_1

    function int_F2(x) result(res)
        real( kind = fgsl_double )   :: res, error
        real(kind=s)                 :: x
        integer( kind = fgsl_size_t) :: neval
        integer( kind = fgsl_int)    :: i
        type( fgsl_function )        :: func

        func = fgsl_function_init( F_2, c_null_ptr )

        i = fgsl_integration_qng ( func,              &
                                 -0.0000000001_fgsl_double,       &
                                 x,       &
                                 1e-5_fgsl_double,      &
                                 1e-5_fgsl_double,      &
                                 res, error, neval )
    end function int_F2

    function F_2( x, params ) bind(c)
        implicit none

        real( kind = c_double )        :: F_2
        real( kind = c_double ), value :: x
        type( c_ptr ), value           :: params

        params = params
        F_2 = (1.0_s/(exp(x)-1.0_s))*x**3.0_s
    end function F_2

end module integrals
