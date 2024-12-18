module config
    use physconst
    use device
    use rng_mod
    use omp_lib
    implicit none

    ! Set dimension
    !integer                         :: SDIM = 0             ! Bulk: 0, 1D: 1
    real(kind=s)            :: mole_x_tmp

    ! All the configuration parameter here
    type(rng_t), dimension(NUM_THREADS),save :: rng
    character(len=40)               :: bulk_filename
    character(len=40)               :: scatrate_filename, dev2_filename
    integer                         :: bulk_filename_len, scatrate_filename_len, dev2_filename_len

    ! Number of used carriers in the ensemble
    integer,parameter                          :: num_carriers = 600000!1000000
    integer,parameter                          :: max_carriers = 2800000 !1150000
    real(kind=s)                               :: p_charge                  ! Superparticle charge (1D: in e/nm^2)


    ! Electrical field
#if DIM == 1
    real(kind=s),dimension(:),allocatable      :: el_field
    real(kind=s),dimension(:),allocatable      :: el_pot
    real(kind=s),dimension(:),allocatable      :: el_charge
    real(kind=s),dimension(:),allocatable      :: el_conc
    integer,dimension(:,:),allocatable         :: particle_out
    integer,dimension(2)                       :: particle_in
    integer                                    :: c_pos
#endif
#if DIM == 2
    real(kind=s),dimension(:,:,:),allocatable    :: el_field
    real(kind=s),dimension(:),allocatable      :: el_pot
    real(kind=s),dimension(:,:),allocatable      :: el_pot2d
    real(kind=s),dimension(:,:),allocatable      :: el_charge
    real(kind=s),dimension(:,:),allocatable      :: el_conc
    integer,dimension(:,:),allocatable         :: particle_out
    integer,dimension(:),allocatable             :: particle_in
    integer                                      :: c_pos, n_out
#endif
#if DIM == 3
    real(kind=s),dimension(:,:,:,:),allocatable  :: el_field
    real(kind=s),dimension(:),allocatable        :: el_pot
    real(kind=s),dimension(:,:,:),allocatable      :: el_charge
    real(kind=s),dimension(:,:,:),allocatable      :: el_conc

    integer                                    :: c_pos
#endif


    !real(kind=s),dimension(3)                  :: e_field ! = (/0.0_s,2e5_s,0.0_s/)/1e9_s
!    real(kind=s),dimension(27)      :: k_array = [5.0_s,6.0_s,7.0_s,8.0_s,9.0_s, &
!                                               & 10.0_s,20.0_s,30.0_s,40.0_s,50.0_s,60.0_s,70.0_s,80.0_s,90.0_s, &
!                                               & 100.0_s,200.0_s,300.0_s,400.0_s,500.0_s,600.0_s,700.0_s,800.0_s,900.0_s, &
!                                               & 1000.0_s, 2000.0_s, 3000.0_s, 4000.0_s]*1e4_s
!    real(kind=s),dimension(43)   :: k_array = [1.0_s,2.0_s,3.0_s,4.0_s,5.0_s,6.0_s,7.0_s,8.0_s,9.0_s,10.0_s,11.0_s,12.0_s, &
!                                            & 13.0_s,14.0_s,15.0_s,17.5_s,20.0_s,22.5_s,25.0_s,27.5_s,30.0_s,35.0_s,40.0_s &
!                                            & ,45.0_s,50.0_s,55.0_s,60.0_s,65.0_s,70.0_s,75.0_s,80.0_s,85.0_s,90.0_s,95.0_s, &
!                                            & 100.0_s,110.0_s,120.0_s,130.0_s,140.0_s,150.0_s,160.0_s,180.0_s,200.0_s]*1e4_s/1e9_s
    real(kind=s),dimension(59)   :: k_array = [5.0_s,6.0_s,7.0_s,8.0_s,9.0_s,10.0_s,11.0_s,12.0_s,13.0_s,14.0_s,15.0_s,17.5_s, &
                                               & 20.0_s,22.5_s,25.0_s,27.5_s,30.0_s,35.0_s,40.0_s,45.0_s,50.0_s, &
                                               & 55.0_s,60.0_s,65.0_s,70.0_s,75.0_s,80.0_s,85.0_s,90.0_s,95.0_s,&
                                               & 100.0_s,110.0_s,120.0_s,130.0_s,140.0_s,150.0_s,160.0_s,180.0_s,200.0_s,220.0_s, &
                                               & 240.0_s,260.0_s,280.0_s, &
                                               & 300.0_s,320.0_s,340.0_s,360.0_s,380.0_s,400.0_s,420.0_s,440.0_s,460.0_s,480.0_s, &
                                               & 500.0_s,520.0_s,540.0_s,560.0_s,580.0_s,600.0_s]*1e4_s/1e9_s
!    real(kind=s),dimension(43)   :: k_array = [5.0_s,6.0_s,7.0_s,8.0_s,9.0_s,10.0_s,15.0_s,20.0_s,25.0_s,30.0_s,35.0_s,40.0_s, &
!                                            & 45.0_s,50.0_s,55.0_s,60.0_s,65.0_s,70.0_s,75.0_s,80.0_s,85.0_s,90.0_s,95.0_s, &
!                                            & 100.0_s,150.0_s,200.0_s,250.0_s,300.0_s,350.0_s,400.0_s,450.0_s,500.0_s,600.0_s, &
!                                            & 650.0_s,700.0_s,750.0_s,800.0_s,850.0_s,900.0_s,950.0_s,1000.0_s,1500.0_s, &
!                                            & 2000.0_s]*1e4_s/1e9_s


    integer,parameter               :: k_max                = 59 !59

    real(kind=s),parameter          :: bin_dE               = 0.03_s
    ! Variables for the simulation time
    integer(kind=8)                 :: n_step, step
    real(kind=s),parameter          :: t_sim                = 5.0_s !80.0_s !40.0_s
    real(kind=s),parameter          :: t_step               = 0.001_s !0.005_s!1.0
    integer                         :: num_synchronized
    real(kind=s)                    :: time_tmp

    ! Material parameter secton; choose models (1 for advanced temperature dependency)
    integer,parameter               :: GAAS_MODEL           = 1
    integer,parameter               :: INGAAS_MODEL         = 1
    integer,parameter               :: INALAS_MODEL         = 1

    !
    ! Device parameter section
    !
    ! Define contact types0V_100K_18
    integer,parameter               :: HETEROINTERFACE      = -1
    integer,parameter               :: INSULATOR            = 1
    integer,parameter               :: SCHOTTKY             = 2
    integer,parameter               :: OHMIC                = 3
    integer,parameter               :: OHMIC_OPEN           = 4

    ! Define orientation for 3D simulation
    integer,parameter               :: INSIDE               = -1
    integer,parameter               :: FRONT                = 1
    integer,parameter               :: BACK                 = 2
    integer,parameter               :: LEFT                 = 3
    integer,parameter               :: RIGHT                = 4
    integer,parameter               :: TOP                  = 5
    integer,parameter               :: BOTTOM               = 6


#if DIM == 0
    integer,parameter,dimension(1,1)                  :: m = 1
    integer,parameter                                 :: IMPURITY_SCAT = 1
    integer,parameter                                 :: pep_hdf5       = 1
    real(kind=s),dimension(1)               :: donor_density = 1e18_s/1e21_s
    integer,parameter,dimension(1)          :: m_material(1) = INGAAS !INALAS
    real(kind=s),dimension(1)               :: scat_m               = (/1/)     ! for the calculation of the scattering rate we have only two different materials
    real(kind=s),dimension(1)               :: scat_rate_comp       = (/1/)

#endif

#if DIM == 1
    real(kind=s),parameter                  :: dev_x                = 87.0_s           ! Device length: x-direction
    integer,dimension(:),allocatable        :: num_nodes
    real(kind=s),dimension(2)               :: mean_vel_ohmic

    ! Mesh configuration
    real(kind=s),parameter,dimension(2,6)             :: m  = reshape((/0.0_s,40.0_s,     & ! buffer
                                                                    &   40.0_s,60.0_s,   & ! channel
                                                                    &   60.0_s,63.0_s,   & ! spacer
                                                                    &   63.0_s,67.0_s,   & ! doping
                                                                    &   67.0_s,77.0_s,   & ! spacer
                                                                    &   77.0_s,87.0_s   &
                                                                    & /),(/2,6/))
    real(kind=s),parameter,dimension(6)               :: m_dx                 = 1.0_s         ! Size of mesh cell
    integer,dimension(:),allocatable                  :: m_region
    integer,parameter,dimension(6)                    :: m_material           = (/INALAS,INGAAS,INALAS,INALAS,INALAS,INALAS/)
    integer,dimension(6)                    :: scat_m               = (/1,2,1,1,1,1/)
    integer,dimension(6)                    :: scat_rate_comp       = (/1,1,0,0,0,0/)

    real(kind=s),dimension(:),allocatable   :: m_lx
    real(kind=s),dimension(:),allocatable   :: node_coord                         ! Give the coordinate of a mesh cell
    real(kind=s),dimension(:,:),allocatable :: node_dist                          ! Gives the distance for a given mesh cell to next and previous node
    real(kind=s),dimension(:),allocatable   :: prob_init

    ! Doping configuration (apply doping in the regions defined by the mesh)
    integer,parameter                                 :: IMPURITY_SCAT = 0
    real(kind=s),parameter,dimension(6)               :: donor_density = (/0.0_s,0.0_s,0.0_s,5e18_s,0.0_s,0.0_s/)/1e21_s
    real(kind=s),dimension(:),allocatable             :: donor

    ! Contact configuration
    integer,dimension(2)                    :: contact_type
    real(kind=s),dimension(2)               :: contact_pot
    real(kind=s),dimension(2)               :: contact_vel
#endif

#if DIM == 2
!    integer,parameter                       :: load_carriers        = 0
!    real(kind=s)                            :: dev_x                = 150.0_s           ! Devfice length: x-direction
!    real(kind=s),parameter                  :: buffer_h             = 0.0_s
!    real(kind=s),parameter                  :: dd_buffer            = 0.0_s
!    integer,parameter                       :: dd_b_int             = int(dd_buffer)
!    real(kind=s),parameter                  :: dev_y                = 350.0_s!+dd_buffer
!    integer,dimension(:),allocatable        :: num_nodes
!    ! Mesh configuration
!    real(kind=s),dimension(4,3)             :: m                    = reshape((/0.0_s,120.0_s,0.0_s,dev_y, &
!                                                                            &  120.0_s,130.0_s,0.0_s,dev_y, &
!                                                                            &  130.0_s,150.0_s,0.0_s,dev_y/),(/4,3/))
!    integer,dimension(3)                    :: scat_m               = (/2,1,2/)
!    integer,dimension(3)                    :: scat_rate_comp       = (/0,1,0/)
!
!    real(kind=s)                            :: m_dx                 = 1.0_s         ! Size of mesh cell
!    real(kind=s)                            :: m_dy                 = 2.0_s
!    real(kind=s)                            :: m_dx_half, m_dy_half
!
!    integer,dimension(:,:),allocatable        :: m_region
!    integer,dimension(3)                      :: m_material = (/INALAS,INGAAS,INALAS/)
!!    integer,dimension(9)                      :: m_material           = (/INALAS,INGAAS,INALAS,INGAAS,INALAS,INALAS,INALAS, &
!!                                                                        & INGAAS,INALAS/)!(/INALAS, INGAAS, INALAS/) (/INGAAS,INGAAS/)!
!    !integer,dimension(10)                      :: m_material    = &
!    !            & (/INALAS,INGAAS,INALAS,INGAAS,INALAS,INALAS,INALAS,INGAAS,AIR,INGAAS/)!(/INALAS, INGAAS, INALAS/) (/INGAAS,INGAAS/)!
!    real(kind=s),dimension(:),allocatable     :: m_lx
!    real(kind=s),dimension(:),allocatable     :: m_ly
!    integer                                   :: num_cells
!    real(kind=s),dimension(:,:),allocatable   :: node_coord, mesh_dist_lr, mesh_dist_ud     ! Give the coordinate of a mesh cell: node_coord(region, x idx, y idx)
!    real(kind=s),dimension(:,:,:),allocatable :: node_dist_lr, node_dist_ud           ! Gives the distance for a given mesh cell to next and previous node
!    real(kind=s),dimension(:,:),allocatable   :: prob_init
!
!    ! Doping configuration (apply doping in the regions defined by the mesh)
!    integer,parameter                         :: IMPURITY_SCAT = 1
!    real(kind=s),dimension(3)                 :: donor_density = (/0.0_s,1e18_s,0.0_s/)/1e21_s
!    real(kind=s),dimension(:,:),allocatable   :: donor

    integer,parameter                       :: load_carriers        = 0
    real(kind=s)                            :: dev_x                = 400.0_s           ! Devfice length: x-direction
    real(kind=s),parameter                  :: buffer_h             = 0.0_s
    real(kind=s),parameter                  :: dd_buffer            = 0.0_s
    integer,parameter                       :: dd_b_int             = int(dd_buffer)
    real(kind=s),parameter                  :: dev_y                = 700.0_s!+dd_buffer
    integer,dimension(:),allocatable        :: num_nodes
!     Mesh configuration
    real(kind=s),dimension(4,4)             :: m                    = reshape((/0.0_s,100.0_s,0.0_s,700.0_s,&
                                                                            &  100.0_s,300.0_s,0.0_s,250.0_s, &
                                                                            &  100.0_s,300.0_s,250.0_s,700.0_s, &
                                                                            &  300.0_s,400.0_s,0.0_s,700.0_s/),(/4,4/))
    integer,dimension(4)                    :: scat_m               = (/1,2,3,4/)
    integer,dimension(4)                    :: scat_rate_comp       = (/1,1,1,1/)

    real(kind=s)                            :: m_dx                 = 5.0_s         ! Size of mesh cell
    real(kind=s)                            :: m_dy                 = 5.0_s
    real(kind=s)                            :: m_dx_half, m_dy_half

    integer,dimension(:,:),allocatable        :: m_region
    integer,dimension(4)                      :: m_material = (/INALAS,INGAAS,INGAAS,INALAS/)
!    integer,dimension(9)                      :: m_material           = (/INALAS,INGAAS,INALAS,INGAAS,INALAS,INALAS,INALAS, &
!                                                                        & INGAAS,INALAS/)!(/INALAS, INGAAS, INALAS/) (/INGAAS,INGAAS/)!
!    integer,dimension(10)                      :: m_material    = &
!                & (/INALAS,INGAAS,INALAS,INGAAS,INALAS,INALAS,INALAS,INGAAS,AIR,INGAAS/)!(/INALAS, INGAAS, INALAS/) (/INGAAS,INGAAS/)!
    real(kind=s),dimension(:),allocatable     :: m_lx
    real(kind=s),dimension(:),allocatable     :: m_ly
    integer                                   :: num_cells
    real(kind=s),dimension(:,:),allocatable   :: node_coord, mesh_dist_lr, mesh_dist_ud     ! Give the coordinate of a mesh cell: node_coord(region, x idx, y idx)
    real(kind=s),dimension(:,:,:),allocatable :: node_dist_lr, node_dist_ud           ! Gives the distance for a given mesh cell to next and previous node
    real(kind=s),dimension(:,:),allocatable   :: prob_init

!     Doping configuration (apply doping in the regions defined by the mesh)
    integer,parameter                         :: IMPURITY_SCAT = 0
    real(kind=s),dimension(4)                 :: donor_density = (/0.0_s,1e17_s,2e16_s,0.0_s/)/1e21_s
    real(kind=s),dimension(:,:),allocatable   :: donor

    ! Contact configuration
    integer,dimension(:,:),allocatable       :: contact_type, contact_id, contact_open
    real(kind=s),dimension(:,:),allocatable  :: contact_pot

    real(kind=s),dimension(:,:),allocatable  :: ohmic_f
    !real(kind=s),dimension(2)               :: contact_vel

    ! Tunnel configuration ! (2,1)
    integer                                  :: n_tunnel       = 0 !5
    real(kind=s),parameter                   :: compute_tunnel = 0.05_s ! Compute tunneling rates every 50fs
    integer,dimension(2,4)                   :: tunnel_x       = reshape((/1,81, &
                                                                         & 561-80,561, &
                                                                         & 561-80,561, &
                                                                         & 1,81/),(/2,4/)) !,       &
                                                                        ! & 52,300/),(/2,5/))
    integer,dimension(2,4)                   :: tunnel_y       = reshape &
                   & ((/61-1,71+1+dd_b_int,61-1,71+1+dd_b_int,61-1,71+1+dd_b_int,61-1,71+1+dd_b_int/),(/2,4/))-5 !,61-1,70+1/),(/2,5/))-5.0_s
    integer,dimension(4)                     :: tunnel_dir     = (/2,2,1,1/)

    ! Pauli exclusion principle configuration
    real(kind=s),parameter                   :: compute_pep    = 0.1_s  ! Every 200fs the FD statistics should be solved
    real(kind=s),parameter                   :: save_pep       = 0.005_s ! How often should the values needed for computation extracted

    ! Output configuration
    ! Choose region to save
    integer,dimension(2)                     :: save_x         = (/1,81/)
    integer,dimension(2)                     :: save_y         = (/1,141/)
    integer                                  :: save_t         = 1! 200
    integer                                  :: l_save_x, l_save_y
    ! Save these computations?? 0: No; 1: Yes
    integer,parameter                        :: tunnel_hdf5    = 0

    ! Save electron temperature and quasi Fermi level?
    integer,parameter                        :: pep_hdf5       = 0
#endif

#if DIM == 3
    real(kind=s)                            :: dev_x                = 700.0_s           ! Device length: x-direction
    real(kind=s)                            :: dev_y                = 100.0_s
    real(kind=s)                            :: dev_z                = 100.0_s
    integer,dimension(:),allocatable        :: num_nodes

    real(kind=s),dimension(6,2)               :: m             = reshape((/0.0_s,250.0_s,0.0_s,100.0_s,0.0_s,100.0_s, &
                                                                         & 250.0_s,700.0_s,0.0_s,100.0_s,0.0_s,100.0_s/),(/6,2/)) !reshape((/0.0_s,20.0_s,80.0_s,20.0_s,80.0_s,1000.0_s/),(/3,2/))
    real(kind=s),dimension(1)                 :: m_dx                 = (/5.0_s/)         ! Size of mesh cell
    integer,dimension(:,:,:),allocatable      :: m_region
    integer,dimension(2)                      :: m_material           = (/INGAAS,INGAAS/)
    real(kind=s),dimension(:),allocatable     :: m_lx
    real(kind=s),dimension(:),allocatable     :: m_ly
    real(kind=s),dimension(:),allocatable     :: m_lz
    integer,dimension(:,:),allocatable      :: A_idx

    real(kind=s),dimension(:,:),allocatable     :: node_coord
    integer,dimension(:),allocatable            :: num_cells
    !real(kind=s),dimension(:,:,:),allocatable   :: node_coord, mesh_dist_lr, mesh_dist_ud     ! Give the coordinate of a mesh cell: node_coord(region, x idx, y idx)
    real(kind=s),dimension(:,:),allocatable     :: mesh_dist_lr, mesh_dist_ud, mesh_dist_bf       ! Gives the distance for a given mesh cell to next and previous node (x,y,z)
    real(kind=s),dimension(:,:,:),allocatable   :: prob_init

    ! Doping configuration (apply doping in the regions defined by the mesh)
    integer,parameter                       :: IMPURITY_SCAT = 0
    real(kind=s),dimension(2)               :: donor_density = (/1e17_s,2e16_s/)/1e21_s
    real(kind=s),dimension(:),allocatable   :: donor

    ! Contact configuration
    integer,dimension(:,:,:),allocatable       :: contact_type      !
    real(kind=s),dimension(:,:,:),allocatable       :: contact_pot      !
    real(kind=s),dimension(:,:,:),allocatable  :: ohmic_f
    !real(kind=s),dimension(2)               :: contact_vel
#endif

    contains

    subroutine init_config()
        integer,dimension(:),allocatable      :: seed
        integer                               :: n, i

        call omp_set_num_threads(NUM_THREADS)
        !print *, omp_get_dynamic()
        call omp_set_dynamic(.true.)
        !print *, omp_get_dynamic()
        do i=1,NUM_THREADS
            call rng_seed(rng(i), 932117+i)
        end do

        ! set seed for random numbers
        call random_seed(size=n)
        allocate(seed(n))
        seed = 51251421
        call random_seed(put=seed)

        n_step = int(t_sim/t_step+0.5_s)
        ! Read the cli arguments
#if DIM == 0
        ! Lattice temperature
        call get_command_argument(1, scatrate_filename)
        read (scatrate_filename, *) T_lattice !donor_density
        ! Energy bin width - 0.03 found in numerical experiments at various temperatures for GaAs
        !call get_command_argument(4, scatrate_filename)
        !read (scatrate_filename, *) bin_dE
        ! Output files
        call get_command_argument(3, scatrate_filename)
        scatrate_filename_len = len_trim(scatrate_filename)
        call get_command_argument(2, bulk_filename)
        bulk_filename_len = len_trim(bulk_filename)
#endif

#if DIM >=1
        ! Get arguments for investigating force problem
        call get_command_argument(1, scatrate_filename)     ! y-Position of second charge
        read (scatrate_filename, *) c_pos
        call get_command_argument(2, dev2_filename)
        dev2_filename_len = len_trim(dev2_filename)

#endif

#if DIM == 1
        ! Set contact type & potential
        contact_type(2) = INSULATOR
        contact_type(1) = SCHOTTKY
        contact_pot(1) = 0.0_s
#endif
    end subroutine init_config

#if DIM == 1
    subroutine init_mesh_1()
        integer                             :: diff_nodes, i, j, max_nodes, idx

        ! Get first the number of different mesh areas (spacing) and allocate arrays
        diff_nodes = size(m(1,:))
        !print *, 'diff', diff_nodes
        allocate(num_nodes(diff_nodes))

        do i=1,diff_nodes
            num_nodes(i) = int((m(2,i)-m(1,i))/m_dx(i))
            !if ((i .eq. 1) .or. (i .eq. diff_nodes)) then
            if (i .eq. diff_nodes) then
                num_nodes(i) = num_nodes(i)+1
            else
                num_nodes(i) = num_nodes(i)
            end if
        end do
        max_nodes = sum(num_nodes)
        !print *, max_nodes, num_nodes
        allocate(node_coord(max_nodes))
        allocate(node_dist(2,max_nodes))
        allocate(m_lx(max_nodes))
        node_dist = 0.0_s

        ! Fill now the coordinate of mesh cells (node_coord) and the cell size (divided by 2: distance to next mesh cell)
        do i=1,diff_nodes
        !print *, num_nodes(i)
            do j=1,num_nodes(i)
                idx = sum(num_nodes(1:i-1))+j
                node_coord(idx) = m(1,i)+(j-1)*m_dx(i) !-m_dx(i)/2.0_s
                ! Compute distance to left border
                if ((idx .eq. 1) .or. (idx .eq. max_nodes)) then
                    node_dist(1,idx) = m_dx(i)/2.0_s
                    node_dist(2,idx) = m_dx(i)/2.0_s
                    m_lx(idx) = node_dist(1,idx)
                else
                    node_dist(2,idx) = m_dx(i)/2.0_s                                                ! We don't care about neighbours, to right edge we have always full mesh distance
                    node_dist(1,idx) = node_coord(idx)-(node_coord(idx-1)+node_dist(2,idx-1))       ! Distance to left
                    m_lx(idx) = node_dist(2,idx)+node_dist(1,idx)
                end if
            end do
        end do

        ! With the known mesh spacings we can now fill unknown properties: donor
        allocate(donor(max_nodes))
        allocate(m_region(max_nodes))
        do i=1,diff_nodes
            do j=1,num_nodes(i)
                idx = sum(num_nodes(1:i-1))+j
                donor(idx) = donor_density(i)!*m_lx(idx)
                m_region(idx) = i
            end do
        end do

        call init_carrier_prob_1()

        allocate(dev%nodes(sum(num_nodes)))
        dev%nodes = node_coord

        !print *, m_lx
        !print *, m_region
        !print *, donor
        !print *, node_coord
        !print *, node_dist(1,:)
        !print *, node_dist(2,:)
        !print *, sum(num_nodes)
    end subroutine init_mesh_1

    subroutine init_carrier_prob_1()
        real(kind=s),dimension(sum(num_nodes))      :: n_doping
        real(kind=s)                                :: n_doping_tot
        integer                 :: i

        allocate(prob_init(sum(num_nodes)))
        n_doping_tot = sum(donor)
        n_doping = donor/n_doping_tot
        n_doping(1) = n_doping(1)/2
        n_doping(sum(num_nodes)) = n_doping(sum(num_nodes))/2
        prob_init = 0.0_s
        do i=1,sum(num_nodes)
            prob_init(i) = sum(n_doping(1:i))
        end do
        !print *, prob_init
    end subroutine init_carrier_prob_1
#endif

#if DIM == 2
    subroutine init_mesh_2()
        integer                             :: num_region, i, j, k
        integer     :: s_ohm, e_ohm

        ! Get first the number of different mesh areas (spacing) and allocate arrays
        num_region = size(m(1,:))
        allocate(num_nodes(2))

        !allocate(scat_rate_comp(size(m_scat)))
        !scat_rate_comp = 0
        !
        ! num_nodes(num_region, x/y direction)
        ! node_coord(num_region, x/y, idx)
        !
        m_dx_half = m_dx/2.0_s
        m_dy_half = m_dy/2.0_s
        ! X direction
        num_nodes(1) = int(dev_x/m_dx+1)!int((m(2,i)-m(1,i))/m_dx)
        ! Y direction
        num_nodes(2) = int(dev_y/m_dy+1)!int((m(4,i)-m(3,i))/m_dx)

        num_cells = num_nodes(1)*num_nodes(2)

        allocate(node_coord(2,maxval(num_nodes)))

        do i=1,num_nodes(1)
            node_coord(1,i) = (i-1)*m_dx
        end do
        do i=1,num_nodes(2)
            node_coord(2,i) = (i-1)*m_dy
        end do

        allocate(m_lx(num_nodes(1)))
        m_lx = m_dx
        m_lx(1) = m_dx/2.0_s
        m_lx(num_nodes(1)) = m_dx/2.0_s
        allocate(m_ly(num_nodes(2)))
        m_ly = m_dy
        m_ly(1) = m_dy/2.0_s
        m_ly(num_nodes(2)) = m_dy/2.0_s
        !print *, node_coord(1,1:num_nodes(1))
        !print *, node_coord(2,1:num_nodes(2))

        ! Allocate and fill arrays for handling contacts and heterointerfaces
        allocate(contact_type(num_nodes(1),num_nodes(2)))
        contact_type = -1
        allocate(contact_pot(num_nodes(1),num_nodes(2)))
        contact_pot = 0.0_s
        allocate(contact_id(num_nodes(1),num_nodes(2)))
        contact_id = 0
        allocate(contact_open(num_nodes(1),num_nodes(2)))
        contact_open = 0
        ! Set number of contacts
        n_out = 2

        s_ohm = 20
        e_ohm = 59

        ! Set contact type and potential
        ! down
        contact_type(1:num_nodes(1),1) = INSULATOR
        contact_type(s_ohm:e_ohm,1) = OHMIC
        contact_pot(s_ohm:e_ohm,1) = 0.0_s
        contact_id(s_ohm:e_ohm,1) = 1
        ! top
        contact_type(1:num_nodes(1),num_nodes(2)) = INSULATOR
        contact_type(s_ohm:e_ohm,num_nodes(2)) = SCHOTTKY
        contact_pot(s_ohm:e_ohm,num_nodes(2)) = -0.7_s
        !contact_open(s_ohm:e_ohm,num_nodes(2)) = 2
        contact_id(s_ohm:e_ohm,num_nodes(2)) = 2
         !contact_type(num_nodes(1)-20:num_nodes(1),num_nodes(2)) = SCHOTTKY
         !contact_pot(num_nodes(1)-20:num_nodes(1),num_nodes(2)) = 0.0_s
        !contact_type(81:121,num_nodes(2)) = SCHOTTKY !81:121
        !contact_pot(81:121,num_nodes(2)) = -0.7_s

        !contact_id(num_nodes(1)-21:num_nodes(1),num_nodes(2)) = 2
        ! left
        contact_type(1,1:num_nodes(2)) = INSULATOR
        !contact_type(1,1:num_nodes(2)) = SCHOTTKY ! 98
        !contact_pot(1,1:num_nodes(2)) = 1_s
        !contact_pot(1,40:50) = 0.0_s! -0.58_s
        !contact_id(1,39:51) = 2
        ! right
        contact_type(num_nodes(1),1:(num_nodes(2))) = INSULATOR
!    contact_type(num_nodes(1),50:num_nodes(2)) = SCHOTTKY
!    contact_pot(num_nodes(1),50:num_nodes(2)) = 0.0_s !-0.738_s
!    contact_type(num_nodes(1),75:100) = SCHOTTKY
!    contact_pot(num_nodes(1),75:100) = 0.0_s !-0.738_s
        !contact_id(num_nodes(1),39:51) = 3

        ! With the known mesh spacings we can now fill unknown properties: donor
        allocate(donor(num_nodes(1),num_nodes(2)))
        allocate(m_region(num_nodes(1),num_nodes(2)))

        !print *, 'm_region start'
        !!$omp parallel private(i,j,k) shared(m_region)
        !!$omp do
        do i=1,num_region
            do j=1,num_nodes(1)
                do k=1,num_nodes(2)
                    if ((node_coord(1,j) <= m(2,i)) .and. (node_coord(1,j) >= m(1,i))) then
                        if ((node_coord(2,k) <= m(4,i)) .and. (node_coord(2,k) >= m(3,i))) then
                            m_region(j,k) = i
                        end if
                    end if
                end do
            end do
        end do
        !!$omp end do
        !!$omp end parallel
        do i=1,num_nodes(1)
            do j=1,num_nodes(2)
                donor(i,j) = donor_density(m_region(i,j))
            end do
        end do
        !print *, donor(:,50)
        !print *, 'm_region filled'
        call init_carrier_prob_2()
        allocate(ohmic_f(num_nodes(1),num_nodes(2)))
        ohmic_f = 1.0_s
        !print *, node_coord(1,:)
        !print *, node_coord(2,:)
        !save_x = (/1,num_nodes(1)/)
        !save_x(2) = num_nodes(1)
        save_y = (/1,num_nodes(2)/)
        print *, save_x, save_y
        l_save_x = save_x(2)-save_x(1)+1
        l_save_y = save_y(2)-save_y(1)+1
        !print *, node_coord(2,:)
        !print *, prob_init(num_nodes(1),:)
        !print *, prob_init(num_nodes(1),:)
        print *, num_nodes(1), num_nodes(2), l_save_x, l_save_y
    end subroutine init_mesh_2

    subroutine init_carrier_prob_2()
        real(kind=s),dimension(:,:),allocatable      :: n_doping
        real(kind=s)                                 :: n_doping_tot
        integer                 :: i, j
        integer                 :: x_idx, y_idx, x, y

        allocate(prob_init(num_nodes(1),num_nodes(2)))
        allocate(n_doping(num_nodes(1),num_nodes(2)))
        !print *, donor(:,50)
        n_doping_tot = 0.0_s
        !do x_idx=1,num_nodes(1)
        !    do y_idx=1,num_nodes(2)
        !        if ((x_idx == num_nodes(1)) .or. (x_idx == 1)) then
        !            n_doping_tot = n_doping_tot+(donor(x_idx,y_idx))/2.0_s
        !        elseif ((y_idx == num_nodes(2)) .or. (y_idx == 1)) then
        !            n_doping_tot = n_doping_tot+(donor(x_idx,y_idx))/2.0_s
        !        else
        !            n_doping_tot = n_doping_tot+donor(x_idx,y_idx)
        !        end if
        !    end do
        !end do

        n_doping_tot = sum(donor)
        n_doping = donor/n_doping_tot
        !print *, sum(n_doping)
        !n_doping(1:num_nodes(1),1) = n_doping(1:num_nodes(1),1)/2.0_s
        !n_doping(1:num_nodes(1),num_nodes(2)) = n_doping(1:num_nodes(1),num_nodes(2))/2.0_s

        !n_doping(1,2:num_nodes(2)-1) = n_doping(1,2:num_nodes(2)-1)/2.0_s
        !n_doping(num_nodes(1),2:num_nodes(2)-1) = n_doping(num_nodes(1),2:num_nodes(2)-1)/2.0_s
        !print *, 'prob start'
        prob_init = 0.0_s
        !$omp parallel private(i,x_idx,y_idx,x,y) shared(prob_init)
        !$omp do
        do i=0,num_cells-1
            !print *, i
                y_idx = i/num_nodes(1)+1
                x_idx = i-((y_idx-1)*num_nodes(1))+1
                !print *, x_idx, y_idx
            do j=0,i
                y = j/num_nodes(1)+1
                x = j-((y-1)*num_nodes(1))+1
                !if (n_doping(x,y) > 0.0_s) then
                prob_init(x_idx,y_idx) = prob_init(x_idx,y_idx)+n_doping(x,y) !prob_init(x_idx,y_idx)+
                !end if
            end do
        end do
        !$omp end do
        !$omp end parallel
        !print *, 'prob end'
        !print *, prob_init
        !print *, sum(prob_init)
    end subroutine init_carrier_prob_2
#endif

#if DIM == 3
    subroutine init_mesh_3()
        integer                             :: num_region, i, j, k, max_nodes, idx, h
        real(kind=s)                        :: m_tmp

        allocate(num_nodes(3))
        num_nodes = 0

        num_region = size(m(1,:))

        ! num_nodes(num_region, x/y/z direction)
        ! node_coord(num_region, x/y/z, idx)

        ! X direction
        num_nodes(1) = int(dev_x/m_dx(1)+1) !int((m(2,1)-m(1,1))/m_dx(1))
        ! Y direction
        num_nodes(2) = int(dev_y/m_dx(1)+1) !int((m(4,1)-m(3,1))/m_dx(1))
        ! Z direction
        num_nodes(3) = int(dev_z/m_dx(1)+1) !int((m(6,1)-m(5,1))/m_dx(1))
        !if (m(2,1) >= dev_x-1e-5_s) then
        !    num_nodes(1) = num_nodes(1)+1
        !end if
        !if (m(4,1) >= dev_y-1e-5_s) then
        !    num_nodes(2) = num_nodes(2)+1
        !end if
        !if (m(6,1) >= dev_z-1e-5_s) then
        !    num_nodes(3) = num_nodes(3)+1
        !end if

        allocate(node_coord(3,maxval(num_nodes)))
        do i=1,3
            do j=1,num_nodes(i)
                node_coord(i,j) = dble(j-1)*m_dx(1)
            end do
        end do

        allocate(m_lx(num_nodes(1)))
        allocate(m_ly(num_nodes(2)))
        allocate(m_lz(num_nodes(3)))
        m_lx = m_dx(1)
        m_ly = m_dx(1)
        m_lz = m_dx(1)
        m_lx(1) = m_lx(1)/2
        m_ly(1) = m_ly(1)/2
        m_lz(1) = m_lz(1)/2
        m_lx(num_nodes(1)) = m_lx(num_nodes(1))/2
        m_ly(num_nodes(2)) = m_ly(num_nodes(2))/2
        m_lz(num_nodes(3)) = m_lz(num_nodes(3))/2

        allocate(num_cells(1))
        num_cells(1) = num_nodes(1)*num_nodes(2)*num_nodes(3)
        !print *, 'Num cells', num_cells(1)
        !print *, 'Node coord', node_coord(3,:)
        ! Allocate and fill arrays for handling contacts and heterointerfaces
        allocate(contact_type(6,maxval(num_nodes),maxval(num_nodes)))
        contact_type = -1
        allocate(contact_pot(6,maxval(num_nodes),maxval(num_nodes)))
        contact_pot = 0.0_s
        ! Set contact type and potential
        !print *, contact_type(BOTTOM,1:10,50:60)
        ! down (y = 1)
        contact_type(BOTTOM,1:num_nodes(1),1:num_nodes(3)) = INSULATOR
        !contact_pot(BOTTOM,1:num_nodes(1),1:num_nodes(3)) = 0.0_s
        ! top (y = max)
        contact_type(TOP,1:num_nodes(1),1:num_nodes(3)) = INSULATOR
        !contact_pot(TOP,1:num_nodes(1),1:num_nodes(3)) = 0.0_s
        ! left (x = 1)
        contact_type(LEFT,1:num_nodes(2),1:num_nodes(3)) = OHMIC
        contact_pot(LEFT,1:num_nodes(2),1:num_nodes(3)) = -0.58_s
        ! right (x = max)
        contact_type(RIGHT,1:num_nodes(2),1:num_nodes(3)) = SCHOTTKY
        contact_pot(RIGHT,1:num_nodes(2),1:num_nodes(3)) = -0.738_s
        ! front (z = 1)
        contact_type(FRONT,1:num_nodes(1),1:num_nodes(2)) = INSULATOR
        !contact_pot(FRONT,1:num_nodes(1),1:num_nodes(2)) = 1.0_s
        ! back (z = max)
        contact_type(BACK,1:num_nodes(1),1:num_nodes(2)) = INSULATOR
        !contact_pot(BACK,1:num_nodes(1),1:num_nodes(2)) = 0.0_s

        ! With the known mesh spacings we can now fill unknown properties: donor
        allocate(donor(num_region))
        allocate(m_region(num_nodes(1),num_nodes(2),num_nodes(3)))

        !print *, num_region
        do i=1,num_region
            do j=1,num_nodes(1)
                do k=1,num_nodes(2)
                    do h=1,num_nodes(3)
                        if ((node_coord(1,j) <= m(2,i)) .and. (node_coord(1,j) >= m(1,i))) then
                            if ((node_coord(2,k) <= m(4,i)) .and. (node_coord(2,k) >= m(3,i))) then
                                if ((node_coord(3,h) <= m(6,i)) .and. (node_coord(3,h) >= m(5,i))) then
                                    m_region(j,k,h) = i
                                end if
                            end if
                        end if
                    end do
                end do
            end do
        end do
        do i=1,num_region
            donor(i) = donor_density(i)
        end do

        call init_carrier_prob_3()
        allocate(ohmic_f(6,maxval(num_nodes),maxval(num_nodes)))
        ohmic_f = 1.0_s
        !print *, donor, m_region(50,10,10)
        !print *, m_lx
        !print *, node_coord(1,:)
    end subroutine init_mesh_3

    subroutine init_carrier_prob_3()
        real(kind=s)                                 :: n_doping_tot
        integer                                      :: i, j
        integer                                      :: x_idx, y_idx, z_idx, x, y

        allocate(A_idx(num_cells(1),3))
        allocate(prob_init(num_nodes(1),num_nodes(2),num_nodes(3)))
        !allocate(prob_init(num_cells(1)))
        n_doping_tot = 0.0_s

        i = 0
        do z_idx=1,num_nodes(3)
            do y_idx=1,num_nodes(2)
                do x_idx=1,num_nodes(1)
                    i = i+1
                    A_idx(i,1:3) = (/x_idx,y_idx,z_idx/)
                    n_doping_tot = n_doping_tot+donor_density(m_region(x_idx,y_idx,z_idx))
                end do
            end do
        end do

        prob_init = 0.0_s
        !do z_idx=1,num_nodes(3)
        !    do y_idx=1,num_nodes(2)
        !        do x_idx=1,num_nodes(1)
        !            prob_init(x_idx,y_idx,z_idx) = prob_init(x_idx,y_idx,z_idx)+ &
        !                                         & donor_density(m_region(x_idx,y_idx,z_idx))/n_doping_tot
        !        end do
        !    end do
        !end do

        do i=2,num_cells(1)
            x_idx = A_idx(i,1)
            y_idx = A_idx(i,2)
            z_idx = A_idx(i,3)
            prob_init(x_idx,y_idx,z_idx) = prob_init(A_idx(i-1,1),A_idx(i-1,2),A_idx(i-1,3))+ &
                        & donor_density(m_region(x_idx,y_idx,z_idx))/n_doping_tot
            !prob_init(i) = prob_init(i-1)+donor_density(m_region(x_idx,y_idx,z_idx))/n_doping_tot
        end do

        !print *, 'prob init', sum(prob_init), n_doping_tot
    end subroutine init_carrier_prob_3

    subroutine coord2idx(x_idx,y_idx,z_idx,orientation,coord1,coord2)
        integer,intent(IN)                      :: x_idx, y_idx, z_idx
        integer,intent(OUT)                     :: orientation, coord1, coord2

        if (x_idx == 1) then
            orientation = LEFT
            coord1 = y_idx
            coord2 = z_idx
        elseif (x_idx == num_nodes(1)) then
            orientation = RIGHT
            coord1 = y_idx
            coord2 = z_idx
        elseif (y_idx == 1) then
            orientation = BOTTOM
            coord1 = x_idx
            coord2 = z_idx
        elseif (y_idx == num_nodes(2)) then
            orientation = TOP
            coord1 = x_idx
            coord2 = z_idx
        elseif (z_idx == 1) then
            orientation = FRONT
            coord1 = x_idx
            coord2 = y_idx
        elseif (z_idx == num_nodes(3)) then
            orientation = BACK
            coord1 = x_idx
            coord2 = y_idx
        else
            orientation = INSIDE
        end if
    end subroutine coord2idx

#endif

#if DIM == 2
    subroutine gauss_input(t,f)
        real(kind=s)                :: t_width, A, b, c, f, t_tmp, t

        t_width = 10.0_s
        t_tmp = t-int(t/t_width)*t_width
        A = 8.0e-2_s
        b = 5.0_s
        c = 1.0_s
        f = A*exp(-(t_tmp-b)**2.0_s/(2.0_s*c**2.0_s))
        !if ((step*t_step >= 2.0_s) .and. (step*t_step <= 5.0_s)) then
        !    f = -1.0_s
        !else
        !    f = 0.0_s
        !end if
    end subroutine gauss_input


    subroutine schottky_diode_dc(t,f)
        real(kind=s)                :: t,f

        if (t < 20_s) then
            f = 0.40_s !0.37_s
        elseif ((t >= 20_s) .and. (t < 40.0_s)) then
            f = 0.40_s
        elseif ((t >= 40) .and. (t < 60.0_s)) then
            f = 0.43_s
        elseif ((t >= 60) .and. (t < 80.0_s)) then
            f = 0.46_s
        elseif ((t >= 80) .and. (t < 100.0_s)) then
            f = 0.49_s
        elseif ((t >= 100) .and. (t < 120.0_s)) then
            f = 0.52_s
        elseif ((t >= 120) .and. (t < 140.0_s)) then
            f = 0.55_s
        elseif ((t >= 140) .and. (t < 160.0_s)) then
            f = 0.58_s
        elseif ((t >= 160) .and. (t < 180.0_s)) then
            f = 0.61_s
        elseif ((t >= 180) .and. (t < 200.0_s)) then
            f = 0.64_s
        elseif ((t >= 200) .and. (t < 220.0_s)) then
            f = 0.67_s
        elseif ((t >= 220) .and. (t < 240.0_s)) then
            f = 0.7_s
        elseif ((t >= 240) .and. (t < 260.0_s)) then
            f = 0.73_s
        elseif ((t >= 260) .and. (t < 280.0_s)) then
            f = 0.76_s
        end if
    end subroutine

    subroutine hemt_dc(t,f)
        real(kind=s)                ::f, t

        if (t < 20.0_s) then
            f = 0.2_s
        elseif ((t >= 20) .and. (t < 30_s)) then
            f = 0.1_s
        elseif ((t >= 30) .and. (t < 40.0_s)) then
            f = 0.2_s
        elseif ((t >= 40) .and. (t < 50.0_s)) then
            f = 0.3_s
        elseif ((t >= 50) .and. (t < 60.0_s)) then
            f = 0.4_s
        elseif ((t >= 60) .and. (t < 70.0_s)) then
            f = 0.5_s
        elseif ((t >= 70) .and. (t < 80.0_s)) then
            f = 0.6_s
        elseif ((t >= 80) .and. (t < 85.0_s)) then
            f = 0.7_s
        elseif ((t >= 85) .and. (t < 90.0_s)) then
            f = 0.8_s
        elseif ((t >= 90) .and. (t < 95.0_s)) then
            f = 0.9_s
        elseif ((t >= 95) .and. (t < 100.0_s)) then
            f = 1.0_s
        elseif ((t >= 100) .and. (t < 105.0_s)) then
            f = 1.1_s
        elseif ((t >= 105) .and. (t < 110.0_s)) then
            f = 1.2_s
        elseif ((t >= 110) .and. (t < 115.0_s)) then
            f = 1.3_s
        elseif ((t >= 115) .and. (t < 120.0_s)) then
            f = 1.4_s
        elseif ((t >= 120) .and. (t < 125.0_s)) then
            f = 1.5_s
        elseif ((t >= 125) .and. (t < 130.0_s)) then
            f = 1.6_s
        elseif ((t >= 130) .and. (t < 135.0_s)) then
            f = 1.7_s
        elseif ((t >= 135) .and. (t < 140.0_s)) then
            f = 1.8_s
        elseif ((t >= 140) .and. (t < 145.0_s)) then
            f = 1.9_s
        elseif ((t >= 145) .and. (t < 150.0_s)) then
            f = 2.0_s
        end if
    end subroutine hemt_dc

    subroutine hemt_dc2(t,f)
        real(kind=s)                ::f, t

        if (t < 8.0_s) then
            f = 0.0_s
        elseif ((t >= 8.0) .and. (t < 10.0_s)) then
            f = 0.0_s
        elseif ((t >= 10.0) .and. (t < 12.0_s)) then
            f = 0.15_s
        elseif ((t >= 12.0) .and. (t < 14.0_s)) then
            f = 0.3_s
        elseif ((t >= 14.0) .and. (t < 16.0_s)) then
            f = 0.45_s
        elseif ((t >= 16.0) .and. (t < 18.0_s)) then
            f = 0.6_s
        elseif ((t >= 18.0) .and. (t < 20.0_s)) then
            f = 0.75_s
        elseif ((t >= 20.0) .and. (t < 22.0_s)) then
            f = 0.9_s
        elseif ((t >= 22.0) .and. (t < 24.0_s)) then
            f = 1.05_s
        elseif ((t >= 24.0) .and. (t < 26.0_s)) then
            f = 1.2_s
        elseif ((t >= 26.0) .and. (t < 28.0_s)) then
            f = 1.35_s
        elseif ((t >= 28.0) .and. (t < 30.0_s)) then
            f = 1.5_s
        elseif ((t >= 30.0) .and. (t < 32.0_s)) then
            f = 1.65_s
        elseif ((t >= 32.0) .and. (t < 34.0_s)) then
            f = 1.8_s
        elseif ((t >= 34.0) .and. (t < 36.0_s)) then
            f = 1.95_s
        end if
        f = dble(c_pos/1000.0_s)
    end subroutine hemt_dc2
#endif
end module config
