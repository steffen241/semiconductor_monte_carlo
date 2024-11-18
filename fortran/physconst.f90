module physconst
    implicit none

    integer, parameter                  :: s   =   selected_real_kind(p=15)

    ! Material index
    integer,parameter                   :: NUM_MAT     = 6
    integer,parameter                   :: MATERIAL    = 2
    integer,parameter                   :: GAAS        = 1
    integer,parameter                   :: INGAAS      = 2
    integer,parameter                   :: INALAS      = 3
    integer,parameter                   :: INAS        = 4
    integer,parameter                   :: ALAS        = 5
    integer,parameter                   :: AIR         = 6

    integer,parameter                   :: G           = 1
    integer,parameter                   :: L           = 2
    integer,parameter                   :: X           = 3

    ! physical constants
    real(kind=s),parameter              :: m0      =   5.68562985e-6_s       ! electron rest mass in eVps^2/nm^2
    real(kind=s),parameter              :: hb      =   6.58211928e-4_s       ! Plank's wirkungsquantum in eVps
    real(kind=s),parameter              :: q       =   1.0_s                 ! electron charge in e
    real(kind=s),parameter              :: q_si    =   1.602176565e-19_s     ! electron charge in SI units
    real(kind=s),parameter              :: kB      =   8.6173323e-5_s        ! Boltzmann constant in eV/K
    real(kind=s),parameter              :: eps_vac =   5.52634959e-2_s;      ! Vacuum permittivity in e^2/eVnm
    real(kind=s)                        :: pi

    ! variables used for calculation of scattering rates
    integer                             :: size_E
    real(kind=s),dimension(300000)      :: E              !0              ! Energy values
    real(kind=s)              :: T_lattice = 300.0_s

    contains
    subroutine init_const
        integer                  :: i
        pi = 4.0_s*atan(1.0_s)
        size_E = size(E)
        ! Energy array
        do i=1,300000!0
            E(i) = i*1.0_s/100000_s !0_s !(i-1.0_s)*1.0_s/100000_s
        end do
    end subroutine init_const
end module physconst
