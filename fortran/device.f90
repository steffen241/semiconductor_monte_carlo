module device
    use physconst
    implicit none

#if DIM == 1
    type                                        :: device_1
    ! Device geometry
    real(kind=s),dimension(:),allocatable       :: nodes
    ! Variable which stores number of particles leaving the device through contacts
    integer,dimension(:,:),allocatable          :: particle_flow_out, particle_flow_in
    ! Place for output of poisson equation
    real(kind=s),dimension(:,:),allocatable     :: g_el_conc, g_el_pot, g_el_field
    ! Variables for macroscopic values
    real(kind=s),dimension(:,:),allocatable     :: g_vel
    real(kind=s),dimension(:,:),allocatable     :: g_energy
    integer,dimension(:,:),allocatable          :: g_num
    contains

    end type device_1

    type(device_1)              :: dev
#endif

#if DIM == 2
    type device_2
    ! Device geometry
!    real(kind=s),dimension(:),allocatable       :: nodes
!    real(kind=s),dimension(:),allocatable       :: time
    ! Variable which stores number of particles leaving the device through contacts
!    integer,dimension(:,:),allocatable          :: particle_flow_out, particle_flow_in
    ! Place for output of poisson equation
!    real(kind=s),dimension(:,:,:),allocatable     :: g_el_conc, g_el_pot
!    real(kind=s),dimension(:,:,:,:),allocatable   :: g_el_field
    ! Variables for macroscopic values
!    real(kind=s),dimension(:,:,:,:),allocatable     :: g_vel
!    real(kind=s),dimension(:,:,:),allocatable     :: g_energy
!    real(kind=s),dimension(:,:,:,:),allocatable     :: tun_energy, tun_prob
!    real(kind=s),dimension(:,:,:),allocatable   :: tun_dest
!    integer,dimension(:,:,:),allocatable          :: g_num
    contains

    end type device_2
    type(device_2)              :: dev
#endif

#if DIM == 3
    type device_3
    ! Device geometry
    real(kind=s),dimension(:),allocatable       :: nodes
    ! Variable which stores number of particles leaving the device through contacts
    !integer,dimension(:,:),allocatable          :: particle_flow_out
    ! Place for output of poisson equation
    real(kind=s),dimension(:,:,:,:),allocatable     :: g_el_conc, g_el_pot
    real(kind=s),dimension(:,:,:,:,:),allocatable   :: g_el_field
    ! Variables for macroscopic values
    real(kind=s),dimension(:,:),allocatable     :: g_vel
    real(kind=s),dimension(:,:),allocatable     :: g_energy
    contains

    end type device_3
    type(device_3)              :: dev
#endif

contains

end module device
