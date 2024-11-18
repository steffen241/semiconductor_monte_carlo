module save_hdf5
    use physconst
    use config
    use poisson
    use hdf5
    implicit none

#if DIM == 2
    !real(kind=s),dimension(:,:),allocatable      :: s_el_conc

    type                        :: save_2
        ! Variables for storing hdf5 information, mainly pointer
        integer(HID_T)          :: file_id       ! File identifier
        integer(HID_T)          :: dset_id       ! Dataset identifier
        integer(HID_T)          :: dspace_id!, field_dspace_id, part_dspace_id     ! Dataspace identifier
        integer(HID_T)          :: group_id      ! Group identifier
        integer(HID_T)          :: memspace, field_memspace, part_memspace, val3_memspace, tun_memspace, pauli_memspace
        CHARACTER(LEN=40)        :: filename = "device.h5" ! File name
        integer                 :: error, i, j
        integer,dimension(50,10) :: vdata
        ! Chunk parameters
        integer(HSIZE_T),dimension(2) :: part_voffs, part_vchunkcount
        integer(HSIZE_T),dimension(3) :: voffset, vchunkcount
        integer(HSIZE_T),dimension(4) :: field_voffs, field_vchunkcount, val3_voffs, val3_vchunkcount
        integer(HSIZE_T),dimension(4) :: tun_voffs, tun_vchunkcount, pauli_voffs!, pauli_vchunkcount

        ! Variables for storing simulation data
        ! Particles flowing in/out contacts
        integer,dimension(:),allocatable           :: part_in, part_out
        ! Number of particles in cell
        integer,dimension(:,:,:),allocatable       :: num_grid
        ! Velocity vector saved on grid
        real(kind=s),dimension(:,:,:),allocatable  :: vel_grid
        ! Mean kinetic energy
        real(kind=s),dimension(:,:),allocatable    :: energy_grid
        ! Tunneling probability
        real(kind=s),dimension(:,:,:),allocatable  :: tunnel_energy, tunnel_prob
        ! Fermi level and electron temperature
        integer                                    :: write_pep = 0
        real(kind=s),dimension(:,:,:),allocatable  :: fermi_lvl, el_temp

    contains
        procedure       :: save_init_2
        procedure       :: write_media_2
        procedure       :: calc_media_2
    end type save_2

    type(save_2)              :: sav
#endif

contains

#if DIM == 2
subroutine save_init_2(this)
    class(save_2)       :: this

    integer             :: error, i, j, num_step, cb_mat
    real(kind=s),dimension(:),allocatable   :: time
    real(kind=s),dimension(:,:,:),allocatable :: cb_offs
    integer,dimension(2)    :: cb_min_idx

    !allocate(s_el_conc(num_nodes(1),num_nodes(2)))
    !s_el_conc = 0.0_s

    allocate(this%part_in(n_out))
    allocate(this%part_out(n_out))

    allocate(this%num_grid(3,num_nodes(1),num_nodes(2)))
    allocate(this%vel_grid(2,num_nodes(1),num_nodes(2)))

    allocate(this%tunnel_energy(1000,3,n_tunnel))
    allocate(this%tunnel_prob(1000,3,n_tunnel))

    allocate(this%fermi_lvl(l_save_y,l_save_x,3))
    allocate(this%el_temp(l_save_y,l_save_x,3))

    allocate(cb_offs(num_nodes(2),num_nodes(1),3))
    cb_offs = 0.0_s
    do i=1,num_nodes(1)
        do j=1,num_nodes(2)
            cb_mat = m_material(m_region(i,j))
            cb_offs(j,i,1) = emin(cb_mat)
            cb_offs(j,i,2) = emin(cb_mat)+valley_offs(cb_mat,2)
            cb_offs(j,i,3) = emin(cb_mat)+valley_offs(cb_mat,3)
        end do
    end do
    cb_min_idx = minloc(cb_offs(:,:,1))
    cb_offs(:,:,2) = cb_offs(:,:,2)-cb_offs(cb_min_idx(1),cb_min_idx(2),2)
    cb_offs(:,:,3) = cb_offs(:,:,3)-cb_offs(cb_min_idx(1),cb_min_idx(2),3)

    num_step = int(n_step/save_t)
    allocate(time(num_step))
    do i=1,num_step
        time(i) = i*save_t*t_step
    end do

    this%filename = dev2_filename(1:dev2_filename_len)

    call h5open_f(error)
    call h5fcreate_f(this%filename, H5F_ACC_TRUNC_F, this%file_id, error)

    call h5gcreate_f(this%file_id, "setup", this%group_id, error)
    call h5gopen_f(this%file_id, "setup", this%group_id, error)

    call h5screate_simple_f(1, (/int(1,8)/), this%dspace_id, error)
    call h5dcreate_f(this%file_id, "setup/charge", H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(p_charge), (/int(1,8)/), error)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    call h5screate_simple_f(1, (/int(l_save_x,8)/), this%dspace_id, error)
    call h5dcreate_f(this%file_id, "setup/node_x", H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(node_coord(1,save_x(1):save_x(2))), (/int(l_save_x,8)/), error)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    call h5screate_simple_f(1, (/int(l_save_y,8)/), this%dspace_id, error)
    call h5dcreate_f(this%file_id, "setup/node_y", H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(node_coord(2,save_y(1):save_y(2))), (/int(l_save_y,8)/), error)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    call h5screate_simple_f(1, (/int(n_step/save_t,8)/), this%dspace_id, error)
    call h5dcreate_f(this%file_id, "setup/time", H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(time), (/int(n_step/save_t,8)/), error)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    call h5screate_simple_f(3, int((/num_nodes(2),num_nodes(1),3/),8), this%dspace_id, error)
    call h5dcreate_f(this%file_id, "setup/cb_offs", H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(cb_offs), int((/num_nodes(2),num_nodes(1),3/),8), error)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    call h5gcreate_f(this%file_id, "particles", this%group_id, error)
    call h5gopen_f(this%file_id, "particles", this%group_id, error)

    call h5screate_simple_f(2, int((/n_out,num_step/),8), this%dspace_id, error)
    call h5dcreate_f(this%file_id, "particles/in" , H5T_NATIVE_INTEGER, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(2, int((/n_out,1/),8), this%part_memspace, error)
    call h5dcreate_f(this%file_id, "particles/out" , H5T_NATIVE_INTEGER, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(2, int((/n_out,1/),8), this%part_memspace, error)
    call h5sclose_f(this%dspace_id, error)

    call h5screate_simple_f(3, int((/l_save_y,l_save_x,num_step/),8), this%dspace_id, error)
    call h5dcreate_f(this%file_id, "el_conc" , H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(3, int((/l_save_y,l_save_x,1/),8), this%memspace, error)
    call h5dcreate_f(this%file_id, "el_pot" , H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(3, int((/l_save_y,l_save_x,1/),8), this%memspace, error)
    call h5dcreate_f(this%file_id, "energy" , H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(3, int((/l_save_y,l_save_x,1/),8), this%memspace, error)
    call h5sclose_f(this%dspace_id, error)

    call h5screate_simple_f(4, int((/l_save_y,l_save_x,2,num_step/),8), this%dspace_id, error)
    call h5dcreate_f(this%file_id, "el_field" , H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(4, int((/l_save_y,l_save_x,2,1/),8), this%field_memspace, error)
    call h5dcreate_f(this%file_id, "velocity" , H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(4, int((/l_save_y,l_save_x,2,1/),8), this%field_memspace, error)
    call h5sclose_f(this%dspace_id, error)

    call h5screate_simple_f(4, int((/l_save_y,l_save_x,3,num_step/),8), this%dspace_id,  error)
    call h5dcreate_f(this%file_id, "particles/cell" , H5T_NATIVE_INTEGER, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(4, int((/l_save_x,l_save_y,3,1/),8), this%val3_memspace, error)
    call h5sclose_f(this%dspace_id, error)

    if (tunnel_hdf5 == 1) then
    call h5gcreate_f(this%file_id, "tunneling", this%group_id, error)
    call h5gopen_f(this%file_id, "tunneling", this%group_id, error)
    call h5screate_simple_f(4, int((/1000,3,n_tunnel,num_step/),8), this%dspace_id,  error)
    call h5dcreate_f(this%file_id, "tunneling/energy" , H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(4, int((/1000,3,n_tunnel,1/),8), this%tun_memspace, error)
    call h5dcreate_f(this%file_id, "tunneling/prob" , H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(4, int((/1000,3,n_tunnel,1/),8), this%tun_memspace, error)
    call h5sclose_f(this%dspace_id, error)
    end if

    if (pep_hdf5 == 1) then
    call h5gcreate_f(this%file_id, "pauli", this%group_id, error)
    call h5gopen_f(this%file_id, "pauli", this%group_id, error)
    call h5screate_simple_f(4, int((/l_save_y,l_save_x,3,int(t_sim/compute_pep)/),8), this%dspace_id,  error)
    call h5dcreate_f(this%file_id, "pauli/fermi_lvl" , H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(4, int((/l_save_y,l_save_x,3,1/),8), this%pauli_memspace, error)
    call h5dcreate_f(this%file_id, "pauli/el_temp" , H5T_NATIVE_REAL, this%dspace_id, this%dset_id, error)
    call h5screate_simple_f(4, int((/l_save_y,l_save_x,3,1/),8), this%pauli_memspace, error)
    call h5sclose_f(this%dspace_id, error)
    end if

    call h5fclose_f(this%file_id, error)
    !call h5close_f(error)
    this%part_voffs = 0
    this%part_vchunkcount = (/n_out,1/)
    this%voffset = 0
    this%vchunkcount = (/l_save_y,l_save_x,1/)
    this%field_voffs = 0
    this%field_vchunkcount = (/l_save_y,l_save_x,2,1/)
    this%val3_voffs = 0
    this%val3_vchunkcount = (/l_save_y,l_save_x,3,1/)
    this%tun_voffs = 0
    this%tun_vchunkcount = (/1000,3,n_tunnel,1/)
    this%pauli_voffs = 0
end subroutine save_init_2

subroutine write_media_2(this)
    class(save_2)   :: this
    integer         :: error
    real(kind=4),dimension(l_save_y,l_save_x,2)   :: val_t
    integer,dimension(l_save_y,l_save_x,3)        :: val3_t

    call h5fopen_f(this%filename, H5F_ACC_RDWR_F, this%file_id, error)

    call h5dopen_f(this%file_id, "particles/in", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%part_voffs, this%part_vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_INTEGER, this%part_in, int((/n_out,1/),8), error, this%part_memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    call h5dopen_f(this%file_id, "particles/out", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%part_voffs, this%part_vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_INTEGER, this%part_out, int((/n_out,1/),8), error, this%part_memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    val3_t(:,:,1) = transpose(this%num_grid(1,save_x(1):save_x(2),save_y(1):save_y(2)))
    val3_t(:,:,2) = transpose(this%num_grid(2,save_x(1):save_x(2),save_y(1):save_y(2)))
    val3_t(:,:,3) = transpose(this%num_grid(3,save_x(1):save_x(2),save_y(1):save_y(2)))
    call h5dopen_f(this%file_id, "particles/cell", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%val3_voffs, this%val3_vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_INTEGER, val3_t, &
    & int((/l_save_y,l_save_x,3/),8), error, this%val3_memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    call h5dopen_f(this%file_id, "el_conc", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%voffset, this%vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(transpose(el_conc(save_x(1):save_x(2),save_y(1):save_y(2)))), &
    & int((/l_save_y,l_save_x/),8), error, this%memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    call h5dopen_f(this%file_id, "el_pot", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%voffset, this%vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(transpose(el_pot2d(save_x(1):save_x(2),save_y(1):save_y(2)))), &
    & int((/l_save_y,l_save_x/),8), error, this%memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    call h5dopen_f(this%file_id, "energy", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%voffset, this%vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(transpose(this%energy_grid(save_x(1):save_x(2),save_y(1):save_y(2)))), &
    & int((/l_save_y,l_save_x/),8), error, this%memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    val_t(:,:,1) = sngl(transpose(el_field(1,save_x(1):save_x(2),save_y(1):save_y(2))))
    val_t(:,:,2) = sngl(transpose(el_field(2,save_x(1):save_x(2),save_y(1):save_y(2))))
    call h5dopen_f(this%file_id, "el_field", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%field_voffs, this%field_vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, val_t,int((/l_save_y,l_save_x,2/),8), error, this%field_memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    val_t(:,:,1) = sngl(transpose(this%vel_grid(1,save_x(1):save_x(2),save_y(1):save_y(2))))
    val_t(:,:,2) = sngl(transpose(this%vel_grid(2,save_x(1):save_x(2),save_y(1):save_y(2))))
    call h5dopen_f(this%file_id, "velocity", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%field_voffs, this%field_vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, val_t, &
    & int((/l_save_y,l_save_x,2/),8), error, this%field_memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)

    if (tunnel_hdf5 == 1) then
    call h5dopen_f(this%file_id, "tunneling/energy", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%tun_voffs, this%tun_vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(this%tunnel_energy), &
    & int((/1000,3,n_tunnel/),8), error, this%tun_memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)
    call h5dopen_f(this%file_id, "tunneling/prob", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%tun_voffs, this%tun_vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(this%tunnel_prob), &
    & int((/1000,3,n_tunnel/),8), error, this%tun_memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)
    end if

    if ((pep_hdf5 == 1) .and. (this%write_pep == 1)) then
    call h5dopen_f(this%file_id, "pauli/fermi_lvl", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%pauli_voffs, this%val3_vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(this%fermi_lvl), &
    & int((/l_save_y,l_save_x,3/),8), error, this%pauli_memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)
    call h5dopen_f(this%file_id, "pauli/el_temp", this%dset_id, error)
    call h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, this%pauli_voffs, this%val3_vchunkcount, error)
    call h5dwrite_f(this%dset_id, H5T_NATIVE_REAL, sngl(this%el_temp), &
    & int((/l_save_y,l_save_x,3/),8), error, this%pauli_memspace, this%dspace_id)
    call h5dclose_f(this%dset_id, error)
    call h5sclose_f(this%dspace_id, error)
    this%write_pep = 0
    this%pauli_voffs(4) = this%pauli_voffs(4)+1
    end if

    call h5fclose_f(this%file_id, error)
    this%part_voffs(2) = this%part_voffs(2)+1
    this%voffset(3) = this%voffset(3)+1
    this%field_voffs(4) = this%field_voffs(4)+1
    this%val3_voffs(4) = this%val3_voffs(4)+1
    this%tun_voffs(4) = this%tun_voffs(4)+1
!    this%pauli_voffs(4) = this%pauli_voffs(4)+1

end subroutine write_media_2

subroutine calc_media_2(this)
    class(save_2)           :: this
    integer                 :: i, j
    integer,dimension(:,:,:),allocatable :: num
    real(kind=s),dimension(:,:),allocatable   :: energy
    real(kind=s),dimension(:,:,:),allocatable :: vel
    real(kind=s)            :: cur_tmp

    allocate(num(3,num_nodes(1),num_nodes(2)))
    num = 0
    allocate(energy(num_nodes(1),num_nodes(2)))
    energy = 0
    allocate(vel(2,num_nodes(1),num_nodes(2)))
    vel = 0

    ! Collect in/out particles
    this%part_in = particle_in
    particle_in = 0
    this%part_out = sum(particle_out,1)
    particle_out = 0

    ! Number of particles in each cell, for each valley
    this%num_grid = 0

    ! Variable for saving tunneling probalities
    this%tunnel_energy = 0.0_s
    this%tunnel_prob = 0.0_s

    ! Calculate mean velocity on grid, get number of particles in each cells first
    !!$OMP PARALLEL DO REDUCTION(+:num)
    do j=1,max_carriers
        if ((c(j)%valley > 0) .and. (c(j)%op == 0)) then
            call get_r_idx_2(c(j)%r,c(j)%r_idx)
            if (c(j)%valley == G) then
                num(G,c(j)%r_idx(1),c(j)%r_idx(2)) = num(G,c(j)%r_idx(1),c(j)%r_idx(2))+1
            elseif (c(j)%valley == L) then
                num(L,c(j)%r_idx(1),c(j)%r_idx(2)) = num(L,c(j)%r_idx(1),c(j)%r_idx(2))+1
            else
                num(X,c(j)%r_idx(1),c(j)%r_idx(2)) = num(X,c(j)%r_idx(1),c(j)%r_idx(2))+1
            end if
            ! mean velocity (vector)
            vel(1,c(j)%r_idx(1),c(j)%r_idx(2)) = vel(1,c(j)%r_idx(1),c(j)%r_idx(2))+c(j)%vel(1)
            vel(2,c(j)%r_idx(1),c(j)%r_idx(2)) = vel(2,c(j)%r_idx(1),c(j)%r_idx(2))+c(j)%vel(2)
            ! compute mean energy in mesh cell
            energy(c(j)%r_idx(1),c(j)%r_idx(2)) = energy(c(j)%r_idx(1),c(j)%r_idx(2))+c(j)%energy
        end if
    end do
    !!$OMP END PARALLEL DO

    this%num_grid = num
    this%energy_grid = energy/sum(num,1)
    this%vel_grid(1,:,:) = vel(1,:,:)/sum(num,1)
    this%vel_grid(2,:,:) = vel(2,:,:)/sum(num,1)


    do i=1,n_tunnel
        this%tunnel_prob(:,:,i) = tun(i)%tunnel_prob
        this%tunnel_energy(:,:,i) = tun(i)%energy
    end do

    do i=1,3
        this%fermi_lvl(:,:,i) = transpose(p%fermi_lvl(i,save_x(1):save_x(2),save_y(1):save_y(2)))
        this%el_temp(:,:,i) = transpose(p%el_temp(i,save_x(1):save_x(2),save_y(1):save_y(2)))
    end do
    !this%fermi_lvl
    !this%el_temp
!    do i=1,num_nodes(1)
!        do j=1,num_nodes(2)
!                if (isnan(this%vel_grid(1,i,j))) this%vel_grid(1,i,j) = -1.0_s
!                if (isnan(this%vel_grid(2,i,j))) this%vel_grid(2,i,j) = -1.0_s
!        end do
!    end do

    do i=1,n_res
    if (res_2(i)%loc == BOTTOM) then
        cur_tmp = res_2(i)%current
        res_2(i)%current = sum(el_conc(res_2(i)%x(1)+2:res_2(i)%x(2)-2,1:5)) &
                     & /((res_2(i)%x(2)-res_2(i)%x(1)+1-2*2)*5.0_s)* &
                     & sum(this%vel_grid(2,res_2(i)%x(1)+2:res_2(i)%x(2)-2,1:5)) &
                     & /((res_2(i)%x(2)-res_2(i)%x(1)+1-2*2)*5.0_s)
        if (res_2(i)%current < 0.0_s) then
            !res_2(i)%current = -res_2(i)%current
        end if
        if ((res_2(i)%current)+1.0_s .eq. (res_2(i)%current))then
            print *, 'detected nan'
            res_2(i)%current = cur_tmp !0.0_s
        end if
        if (isnan(res_2(i)%current))then
            print *, 'detected nan'
            res_2(i)%current = cur_tmp !0.0_s
        end if
    end if
    if (res_2(i)%loc == TOP) then
        cur_tmp = res_2(i)%current
        res_2(i)%current = sum(el_conc(res_2(i)%x(1)+2:res_2(i)%x(2)-2,num_nodes(2)-4:num_nodes(2))) &
                     & /((res_2(i)%x(2)-res_2(i)%x(1)+1-2*2)*5.0_s)* &
                     & sum(this%vel_grid(2,res_2(i)%x(1)+2:res_2(i)%x(2)-2,num_nodes(2)-4:num_nodes(2))) &
                     & /((res_2(i)%x(2)-res_2(i)%x(1)+1-2*2)*5.0_s)
        !if (res_2(i)%current > 0.0_s) then
            res_2(i)%current = -res_2(i)%current
        !end if
        if ((res_2(i)%current)+1.0_s .eq. (res_2(i)%current))then
            print *, 'detected nan'
            res_2(i)%current = cur_tmp!0.0_s
        end if
        if (isnan(res_2(i)%current))then
            print *, 'detected nan'
            res_2(i)%current = cur_tmp!0.0_s
        end if
    end if
    !res_2(i)%current = 0.0_s
    !res_2(i)%current = res_2(i)%current)-15.0_s
    !if ((n_res == 2) .and. (el_pot2d(res_2(1)%x(1),1) == el_pot2d(res_2(2)%x(1),num_nodes(2)))) then
    !    res_2(i)%current = 0.0_s
    !end if
    end do
                                !sum(el_conc(res_2(i)%x(1)+2:res_2(i)%x(2)-2))
    !res_2(1)%current = (el_conc(45,1)*this%vel_grid(2,45,1)+el_conc(43,2)*this%vel_grid(2,43,2)+&
    !& el_conc(46,3)*this%vel_grid(2,46,3)+el_conc(48,4)*this%vel_grid(2,48,4)+el_conc(48,5)*this%vel_grid(2,48,5))/5.0_s
            print *, res_2(1)%current/res_2%doping
   ! res_2(2)%current = (el_conc(45,400)*this%vel_grid(2,45,400)+el_conc(43,399)*this%vel_grid(2,43,399)+&
!& el_conc(46,398)*this%vel_grid(2,46,398)+el_conc(48,397)*this%vel_grid(2,48,397)+el_conc(48,396)*this%vel_grid(2,48,396))/5.0_s
      !      print *, res_2(2)%current/res_2%doping
end subroutine calc_media_2
#endif

end module save_hdf5
