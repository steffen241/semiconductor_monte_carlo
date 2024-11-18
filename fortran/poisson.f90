module poisson
    use config
    use mc_core
    use device
    use reservoir
    use tunneling
    use pauli

    implicit none
    include "hips.inc"
    include "mpif.h"

    real(kind=s),dimension(:,:),allocatable    :: A
    integer,dimension(:,:,:),allocatable       :: A_sub2idx
    integer,dimension(:,:),allocatable         :: A_idx2sub
    real(kind=s),dimension(:),allocatable      :: u_pre
    integer                                    :: n, nnz
#if DIM == 1
    real(kind=s),dimension(:),allocatable      :: eps_dev
#endif
#if DIM == 2
    real(kind=s),dimension(:,:),allocatable    :: eps_dev
#endif
    contains
#if DIM == 1
    subroutine init_poisson_1()
        real(kind=s),dimension(:),allocatable  :: A_tmp
        integer                                :: i, j
        real(kind=s)                           :: dl, dr
        !integer                                :: max_nodes

        allocate(dev%particle_flow_out(n_step,2))
        dev%particle_flow_out = 0
        allocate(dev%particle_flow_in(n_step,2))
        dev%particle_flow_in = 0

        ! Allocate space for electric field, potential and charge
        allocate(el_field(sum(num_nodes)))
        allocate(el_pot(sum(num_nodes)))
        allocate(el_charge(sum(num_nodes)))
        allocate(el_conc(sum(num_nodes)))
        allocate(particle_out(NUM_THREADS,2))
        particle_out = 0
        particle_in = 0

        allocate(dev%g_el_conc(n_step,sum(num_nodes)))
        allocate(dev%g_el_pot(n_step,sum(num_nodes)))
        allocate(dev%g_el_field(n_step,sum(num_nodes)))

        dev%g_el_conc = 0.0_s
        dev%g_el_pot = 0.0_s
        dev%g_el_field = 0.0_s

        el_pot = 0.0_s
        el_field = 0.0_s !-1.0e4_s/1e9_s

        allocate(dev%g_vel(n_step,sum(num_nodes)))
        dev%g_vel = 0.0_s

        allocate(dev%g_energy(n_step,sum(num_nodes)))
        dev%g_energy = 0.0_s

        allocate(dev%g_num(n_step,sum(num_nodes)))
        dev%g_num = 0

        ! Calculate superparticle charge
        p_charge = 0.0_s
        j = size(m(1,:))
        print *, j
        do i=1,j
            p_charge = p_charge+donor_density(i)*(m(2,i)-m(1,i))/num_carriers
            !print *, p_charge, 'range:', (m(i,2)-m(i,1)+margin)
        end do
        print *, 'Superparticle charge: (e/nm^2):', p_charge

        allocate(A(sum(num_nodes),sum(num_nodes)))
        allocate(A_tmp(sum(num_nodes)))
        A = 0.0_s
        A_tmp = 0.0_s

        ! Fill matrix
        do i=2,sum(num_nodes)-1
            A_tmp = 0.0_s
            dl = node_coord(i)-node_coord(i-1)
            dr = node_coord(i+1)-node_coord(i)
            A_tmp(i-1) = 2.0_s/(dl*(dl+dr)) !-
            A_tmp(i) = -2/(dl*dr)
            A_tmp(i+1) = 2.0_s/(dr*(dr+dl)) !-
            A(i,:) = A_tmp !*4.0_s
            !print *, A(i,:)
        end do

        ! Set boundary conditions
        if ((contact_type(1) .eq. OHMIC) .or. (contact_type(1) .eq. SCHOTTKY)) then
            A(1,1) = 1.0_s
        else
            dr = node_coord(2)
            A(1,1) = 1.0_s/(dr*dr)
            A(1,2) = -1.0_s/(dr*dr)
        end if
        if ((contact_type(2) .eq. OHMIC) .or. (contact_type(2) .eq. SCHOTTKY)) then
            A(sum(num_nodes),sum(num_nodes)) = 1.0_s
        else
            dl = node_coord(sum(num_nodes))-node_coord(sum(num_nodes)-1)
            A(sum(num_nodes),sum(num_nodes)-1) = -1.0_s/(dl*dl)
            A(sum(num_nodes),sum(num_nodes)) = 1.0_s/(dl*dl)
        end if

        allocate(eps_dev(sum(num_nodes)))
        do i=1,sum(num_nodes)
            eps_dev(i) = eps_static(m_material(m_region(i)))
        end do
    end subroutine init_poisson_1

    subroutine solve_poisson_1()
        real(kind=8),dimension(sum(num_nodes))    :: ipiv, rho
        real(kind=8),dimension(sum(num_nodes))    :: pot_new
        integer                              :: info, i
        real(kind=s),dimension(:,:),allocatable  :: M
        real(kind=s)                        :: f

        allocate(M(sum(num_nodes),sum(num_nodes)))
        M = A

        ! Set up RHS
        rho = (-el_conc+donor)
        pot_new = -rho/eps_dev !(eps_static(INGAAS))

        ! Set boundary conditions
        if ((contact_type(1) .eq. OHMIC) .or. (contact_type(1) .eq. SCHOTTKY)) then
            call schottky_diode_dc(step*t_step,f)
            pot_new(1) = -f
            !pot_new(1) = contact_pot(1)
        else
            pot_new(1) = rho(1)/2.0_s
        end if
        if ((contact_type(2) .eq. OHMIC) .or. (contact_type(2) .eq. SCHOTTKY)) then
            pot_new(sum(num_nodes)) = contact_pot(2)
            !call schottky_diode_dc(step*t_step,f)
            !pot_new(sum(num_nodes)) = contact_pot(2)+f
        else
            pot_new(sum(num_nodes)) = rho(sum(num_nodes))/2.0_s
        end if

        !open (30, file='/tmp/rho.txt', access='APPEND', status='REPLACE', iostat=st2)
        !write (30,*) -rho/eps_static(INGAAS)
        !print *, el_conc
        !print *, donor
        !print *, -rho
        call dgesv(sum(num_nodes), 1, M, sum(num_nodes), ipiv, pot_new, sum(num_nodes), info)
        !print *, 'poisson info', info
        el_pot = pot_new

        do i=2,(sum(num_nodes)-1)
            el_field(i) = -(el_pot(i+1)-el_pot(i-1))/(node_coord(i+1)-node_coord(i-1))
            !print *, (node_coord(i+1)-node_coord(i-1))
        end do
        el_field(1) = el_field(2)
        el_field(sum(num_nodes)) = el_field(sum(num_nodes)-1)

        do i=1,sum(num_nodes)
            el_pot(i) = el_pot(i)-emin(m_material(m_region(i)))
        end do
        ! Save output in variables
        dev%g_el_conc(step,:) = el_conc
        dev%g_el_pot(step,:) = el_pot
        dev%g_el_field(step,:) = el_field
    end subroutine solve_poisson_1

    subroutine assign_charge_1()
        integer                         :: i, k, max_nodes
        real(kind=s)                    :: x

        max_nodes = sum(num_nodes)
        el_charge = 0.0_s
        el_conc = 0.0_s
        do k=1,max_carriers
            if ((c(k)%valley > 0) .and. (c(k)%op == 0)) then
                x = c(k)%r(1)/m_dx(m_region(c(k)%r_idx))
                i = int(x+1)
                el_charge(i) = el_charge(i)+(1-(x-(i-1)))
                el_charge(i+1) = el_charge(i+1)+(x-(i-1))
            end if
        end do

        !do i=1,max_carriers
        !    if (c(i)%valley .ne. 0) then
        !        cell_idx = minloc(abs(c(i)%r(1)-node_coord))
        !        el_charge(cell_idx(1)) = el_charge(cell_idx(1))+1.0_s
        !    end if
        !end do
        el_charge = el_charge*p_charge
        el_conc = el_charge/m_lx
    end subroutine assign_charge_1

    subroutine media()
        integer                                     :: j
        real(kind=s),dimension(:),allocatable       :: mean_vel, mean_vel_tot, mean_energy, mean_energy_tot
        integer,dimension(:),allocatable            :: num

        allocate(mean_vel(max_carriers))
        allocate(mean_energy(max_carriers))
        allocate(mean_vel_tot(sum(num_nodes)))
        allocate(mean_energy_tot(sum(num_nodes)))
        allocate(num(sum(num_nodes)))
        mean_vel_tot = 0.0_s
        mean_energy_tot = 0.0_s
        num = 0

        do j=1,max_carriers
            if ((c(j)%valley > 0) .and. (c(j)%op == 0)) then
                call get_r_idx(c(j)%r,c(j)%r_idx)
                mean_vel_tot(c(j)%r_idx) = mean_vel_tot(c(j)%r_idx)+c(j)%vel(1)
                mean_energy_tot(c(j)%r_idx) = mean_energy_tot(c(j)%r_idx)+c(j)%energy
                num(c(j)%r_idx) = num(c(j)%r_idx)+1
            end if
        end do
        !do x_idx=1,sum(num_nodes)
        !    if (num(x_idx)== 0) then
        !        num(x_idx) = 1e7
        !    end if
        !end do
        dev%particle_flow_out(step,:) = sum(particle_out,1)
        dev%particle_flow_in(step,:) = particle_in
        !dev%particle_flow_out(step,1) = sum(particle_out(:,1))
        !dev%particle_flow_out(step,2) = sum(particle_out(:,2))
        particle_out = 0
        particle_in = 0
        dev%g_vel(step,:) = mean_vel_tot/num
        dev%g_energy(step,:) = mean_energy_tot/num
        dev%g_num(step,:) = num
    end subroutine media
#endif

#if DIM == 2
    subroutine init_poisson_hips_2()
        integer                                :: x_idx, y_idx, nnz_idx, n_idx, A_idx, idnbr, ierr, id, ntasks, proc_id
        integer                                :: num_region, i
        real(kind=8),dimension(:),allocatable  :: aval, xx
        integer,dimension(:),allocatable       :: ia, ja
        real(kind=8)                           :: dl, dr, du, dd, prec

        call MPI_Comm_rank(MPI_COMM_WORLD, proc_id , ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, ntasks , ierr)

        !if (proc_id == 0) then
        ! Allocate space for electric field, potential and charge
        allocate(el_field(2,num_nodes(1),num_nodes(2)))
        allocate(el_pot(num_cells))
        allocate(el_charge(num_nodes(1),num_nodes(2)))
        allocate(el_conc(num_nodes(1),num_nodes(2)))
        allocate(el_pot2d(num_nodes(1),num_nodes(2)))
!
!        allocate(dev%g_el_conc(int(n_step/save_t),l_save_x,l_save_y)) !,num_nodes(1),num_nodes(2)))
!        allocate(dev%g_el_pot(int(n_step/save_t),l_save_x,l_save_y))
!        allocate(dev%g_el_field(int(n_step/save_t),2,l_save_x,l_save_y))
!        allocate(dev%time(int(n_step/save_t)))
!
!        dev%g_el_conc = 0.0_s
!        dev%g_el_pot = 0.0_s
!        dev%g_el_field = 0.0_s
!        dev%time = 0.0_s

        allocate(A_sub2idx(1,num_nodes(1),num_nodes(2)))
        allocate(A_idx2sub(num_cells,3))

!        allocate(dev%g_vel(int(n_step/save_t),2,l_save_x,l_save_y))
!        dev%g_vel = 0.0_s
!
!        allocate(dev%g_energy(int(n_step/save_t),l_save_x,l_save_y))
!        dev%g_energy = 0.0_s
!        allocate(dev%g_num(int(n_step/save_t),l_save_x,l_save_y))
!        dev%g_num = 0

        ! Calculate superparticle charge
        num_region = size(m(1,:))
        p_charge = 0.0_s
        do i=1,num_region
            p_charge = p_charge+donor_density(i)*((m(2,i)-m(1,i))*(m(4,i)-m(3,i)))/num_carriers
        end do
        print *, 'Superparticle charge: (e/nm):', p_charge

        call MPI_Comm_rank(MPI_COMM_WORLD, proc_id , ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, ntasks , ierr)

        allocate(eps_dev(num_nodes(1),num_nodes(2)))

        ! We use CSR format
        ! Compute number of non-zero elements and number of unknowns
        n = num_cells
        nnz = 0
        do y_idx=1,num_nodes(2)
            do x_idx=1,num_nodes(1)
                eps_dev(x_idx,y_idx) = eps_static(m_material(m_region(x_idx,y_idx)))
               ! print *, x_idx, y_idx, m_material(m_region(x_idx,y_idx)), eps_static(m_material(m_region(x_idx,y_idx)))
                if (contact_type(x_idx,y_idx) == -1) then
                    nnz = nnz+5
                end if
                if ((contact_type(x_idx,y_idx) == SCHOTTKY) .or. &
                & ((contact_type(x_idx,y_idx) == OHMIC) .and. contact_open(x_idx,y_idx) == 0))  then
                    nnz = nnz+1
                end if
                if ((contact_type(x_idx,y_idx) == INSULATOR) .or. &
                & ((contact_type(x_idx,y_idx) == OHMIC) .and. contact_open(x_idx,y_idx) == 1))  then
                    if ((x_idx == num_nodes(1) .and. (y_idx == num_nodes(2))) .or. &
                      & ((x_idx == 1) .and. (y_idx == num_nodes(2))) .or. &
                      & ((x_idx == num_nodes(1)) .and. (y_idx == 1)) .or. &
                      & ((x_idx == 1) .and. (y_idx == 1))) then
                        nnz = nnz+3
                    else
                        nnz = nnz+4
                    end if
                end if
            end do
        end do
        print *, 'n:', n, 'nnz:', nnz

        ! Fill sparse martrix, CSR format
        allocate(aval(nnz))
        allocate(ja(nnz))
        allocate(ia(n+1))
        allocate(xx(n))
        allocate(u_pre(num_cells))

        nnz_idx = 1
        n_idx = 1
        do y_idx=1,num_nodes(2)
            do x_idx=1,num_nodes(1)
                A_idx = (y_idx-1)*num_nodes(1)+x_idx
                A_sub2idx(1,x_idx,y_idx) = A_idx
                A_idx2sub(A_idx,1:3) = (/1,x_idx, y_idx/)
                ! Compute coefficients
                dl = m_dx !mesh_dist_lr(1,1,x_idx)
                dr = m_dx !mesh_dist_lr(1,2,x_idx)
                dd = m_dy !mesh_dist_ud(1,1,y_idx)
                du = m_dy !mesh_dist_ud(1,2,y_idx)
                if (contact_type(x_idx,y_idx) == -1) then
                    aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                    ja(nnz_idx) = A_idx-num_nodes(1)
                    ia(n_idx) = nnz_idx
                    nnz_idx = nnz_idx+1

                    aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                    ja(nnz_idx) = A_idx-1
                    nnz_idx = nnz_idx+1

                    aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                    ja(nnz_idx) = A_idx
                    nnz_idx = nnz_idx+1

                    aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                    ja(nnz_idx) = A_idx+1
                    nnz_idx = nnz_idx+1

                    aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                    ja(nnz_idx) = A_idx+num_nodes(1)
                    nnz_idx = nnz_idx+1
                end if

                ! Schottky or ohmic contacts, fixed potential
                if ((contact_type(x_idx,y_idx) == SCHOTTKY) .or. &
                & ((contact_type(x_idx,y_idx) == OHMIC) .and. contact_open(x_idx,y_idx) == 0))  then
                    aval(nnz_idx) = 1.0_s
                    ja(nnz_idx) = A_idx
                    ia(n_idx) = nnz_idx
                    u_pre(n_idx) = contact_pot(x_idx,y_idx)
                    nnz_idx = nnz_idx+1
                end if

                ! Neumann boundary condition, E_n = 0
                if ((contact_type(x_idx,y_idx) == INSULATOR) .or. &
                & ((contact_type(x_idx,y_idx) == OHMIC) .and. contact_open(x_idx,y_idx) == 1))  then
                            if ((x_idx == num_nodes(1)) .and. (y_idx == num_nodes(2))) then
                                ia(n_idx) = nnz_idx
                                aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 4.0/(dd*(du+dd))
                                ja(nnz_idx) = A_idx-num_nodes(1)
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 4.0_s/(dl*(dl+dr))
                                ja(nnz_idx) = A_idx-1
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !-2.0_s/(du*dd)!
                                ja(nnz_idx) = A_idx
                                nnz_idx = nnz_idx+1
                            elseif ((x_idx == 1) .and. (y_idx == num_nodes(2))) then
                                ia(n_idx) = nnz_idx
                                aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 4.0/(dd*(du+dd))
                                ja(nnz_idx) = A_idx-num_nodes(1)
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !-2.0_s/(du*dd)!
                                ja(nnz_idx) = A_idx
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 4.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 4.0_s/(dr*(dl+dr))
                                ja(nnz_idx) = A_idx+1
                                nnz_idx = nnz_idx+1
                            elseif ((x_idx == num_nodes(1)) .and. (y_idx == 1)) then
                                ia(n_idx) = nnz_idx
                                aval(nnz_idx) = 4.0_s/(dl*(dl+dr))
                                ja(nnz_idx) = A_idx-1
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !-2.0_s/(du*dd)!
                                ja(nnz_idx) = A_idx
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 4.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 4.0/(du*(du+dd))
                                ja(nnz_idx) = A_idx+num_nodes(1)
                                nnz_idx = nnz_idx+1
                            elseif ((x_idx == 1) .and. (y_idx == 1)) then
                                ia(n_idx) = nnz_idx
                                aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !-2.0_s/(du*dd)!
                                ja(nnz_idx) = A_idx
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 4.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 4.0_s/(dr*(dl+dr))
                                ja(nnz_idx) = A_idx+1
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 4.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 4.0/(du*(du+dd))
                                ja(nnz_idx) = A_idx+num_nodes(1)
                                nnz_idx = nnz_idx+1
                            elseif (y_idx == 1) then
                                ia(n_idx) = nnz_idx
                                aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                                ja(nnz_idx) = A_idx-1
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                                ja(nnz_idx) = A_idx
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                                ja(nnz_idx) = A_idx+1
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 4.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                                ja(nnz_idx) = A_idx+num_nodes(1)
                                nnz_idx = nnz_idx+1
                            elseif (y_idx == num_nodes(2)) then
                                ia(n_idx) = nnz_idx
                                aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                                ja(nnz_idx) = A_idx-num_nodes(1)
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                                ja(nnz_idx) = A_idx-1
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                                ja(nnz_idx) = A_idx
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                                ja(nnz_idx) = A_idx+1
                                nnz_idx = nnz_idx+1
                            elseif (x_idx == 1) then
                                ia(n_idx) = nnz_idx
                                aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                                ja(nnz_idx) = A_idx-num_nodes(1)
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                                ja(nnz_idx) = A_idx
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 4.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                                ja(nnz_idx) = A_idx+1
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                                ja(nnz_idx) = A_idx+num_nodes(1)
                                nnz_idx = nnz_idx+1
                            elseif (x_idx == num_nodes(1)) then
                                ia(n_idx) = nnz_idx
                                aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                                ja(nnz_idx) = A_idx-num_nodes(1)
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 4.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                                ja(nnz_idx) = A_idx-1
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                                ja(nnz_idx) = A_idx
                                nnz_idx = nnz_idx+1

                                aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                                ja(nnz_idx) = A_idx+num_nodes(1)
                                nnz_idx = nnz_idx+1
                            end if
                end if
                n_idx = n_idx+1
            end do
        end do
        ia(n+1) = nnz+1
        !print *, rhs
        !end if
print *, nnz_idx, n_idx, ia(n+1), n+1

        ! Init the HIPS package
        idnbr = 1
        call HIPS_INITIALIZE(idnbr, ierr)
        id = 0

        ! Set Options
        prec = 1e-5_s
        !call HIPS_SetDefaultOptions(id, HIPS_ITERATIVE, ierr)
        print *, n, ntasks
        if (ntasks > 1) then
            call HIPS_SetDefaultOptions(id, HIPS_HYBRID, ierr)
            !CALL HIPS_SetOptionINT(id, HIPS_PARTITION_TYPE, 0, ierr)
            CALL HIPS_SetOptionINT(id, HIPS_DOMSIZE, n/ntasks, ierr)
        else
            call HIPS_SetDefaultOptions(id, HIPS_ITERATIVE, ierr)
            call HIPS_SetOptionINT(id, HIPS_DOMNBR, ntasks, ierr)
        end if
        !call HIPS_SetOptionINT(id, HIPS_SYMMETRIC, 1, ierr)
        call HIPS_SetOptionREAL(id, HIPS_PREC, prec, ierr)
        call HIPS_SetOptionINT(id, HIPS_LOCALLY, 0, ierr)
        call HIPS_SetOptionINT(id, HIPS_ITMAX, 100, ierr)
        call HIPS_SetOptionINT(id, HIPS_KRYLOV_RESTART, 50, ierr)
        call HIPS_SetOptionINT(id, HIPS_VERBOSE, 0, ierr)
        call HIPS_SetOptionINT(id, HIPS_CHECK_GRAPH, 1, ierr)
        call HIPS_SetOptionINT(id, HIPS_CHECK_MATRIX, 1, ierr)

        call HIPS_SetOptionINT(id, HIPS_DISABLE_PRECOND, 0, ierr)

        call HIPS_GRAPHGLOBALCSR(id, n, int(ia,8), ja, 0, ierr)
        call HIPS_MATRIXGLOBALCSR(id, n, int(ia,8), ja, aval,  0, HIPS_ASSEMBLY_OVW, 0, ierr)

        !call HIPS_SETGLOBALRHS(id, u_pre, 0, 0, ierr)
        !call HIPS_EXITONERROR(ierr)
        !call HIPS_GETGLOBALSOLUTION(id, xx, 0, ierr)
        !call HIPS_EXITONERROR(ierr)

        !if (proc_id == 0) then
        !    dev%g_el_pot(step,:) = xx
        !end if
        !print *, rhs

        !call HIPS_CLEAN(id, ierr)
        !call HIPS_FINALIZE(ierr)
        !call HIPS_EXITONERROR(ierr)

    end subroutine init_poisson_hips_2

    subroutine solve_poisson_hips_2()
        integer                                :: proc_id, ierror, x_idx, y_idx, s_step
        integer,dimension(3)                   :: idx
        real(kind=8)                           :: rho, f


        call MPI_Comm_rank(MPI_COMM_WORLD, proc_id , ierror)
        if (proc_id == 0) then
        do y_idx=1,num_nodes(2)
                do x_idx=1,num_nodes(1)
                    rho = -el_conc(x_idx,y_idx)+donor(x_idx,y_idx)
                    !rho = 0.0_s
                    !print *, x_idx, y_idx, eps_dev(x_idx,y_idx)
                    u_pre(A_sub2idx(1,x_idx,y_idx)) = -rho/eps_dev(x_idx,y_idx)
                    if ((contact_type(x_idx,y_idx) == SCHOTTKY) .or. &
                    & ((contact_type(x_idx,y_idx) == OHMIC) .and. contact_open(x_idx,y_idx) == 0))  then
                        !if ((y_idx == num_nodes(2))) then !.and. (x_idx >= 200)) then
                        if ((y_idx == 1)) then
                            call hemt_dc2(step*t_step,f)
                            u_pre(A_sub2idx(1,x_idx,y_idx)) = -f
                        else
                            u_pre(A_sub2idx(1,x_idx,y_idx)) = contact_pot(x_idx,y_idx)
                        end if
                    end if
                    if ((contact_type(x_idx,y_idx) == INSULATOR) .or. &
                & ((contact_type(x_idx,y_idx) == OHMIC) .and. contact_open(x_idx,y_idx) == 1))  then
                        u_pre(A_sub2idx(1,x_idx,y_idx)) = rho/eps_dev(x_idx,y_idx)/1.0_s!rho/1.0_s !2.0_s !1.0_s

                    end if
                end do
            end do
        end if
        call HIPS_SETGLOBALRHS(0, u_pre, 0, 0, ierror)
        call HIPS_EXITONERROR(ierror)
        call HIPS_GETGLOBALSOLUTION(0, el_pot, 0, ierror)
        call HIPS_EXITONERROR(ierror)

        if (proc_id == 0) then
        ! E-Field in x-direction
        do x_idx=2,num_nodes(1)-1
            do y_idx=1,num_nodes(2)
                el_field(1,x_idx,y_idx) = -(el_pot(A_sub2idx(1,x_idx+1,y_idx))-el_pot(A_sub2idx(1,x_idx-1,y_idx)))/ &
                                                  &  (node_coord(1,x_idx+1)-node_coord(1,x_idx-1))
            end do
        end do
        ! E-Field in y-direction
        do y_idx=2,num_nodes(2)-1
        do x_idx=1,num_nodes(1)
             el_field(2,x_idx,y_idx) = -(el_pot(A_sub2idx(1,x_idx,y_idx+1))-el_pot(A_sub2idx(1,x_idx,y_idx-1)))/ &
                                                  &  (node_coord(2,y_idx+1)-node_coord(2,y_idx-1))
            end do
        end do
        do x_idx=1,num_nodes(1)
            el_field(2,x_idx,1) = el_field(2,x_idx,2)
            el_field(2,x_idx,num_nodes(2)) = el_field(2,x_idx,num_nodes(2)-1)
        end do
        do y_idx=1,num_nodes(2)
            el_field(1,1,y_idx) = el_field(1,2,y_idx)
            el_field(1,num_nodes(1),y_idx) = el_field(1,num_nodes(1)-1,y_idx)
        end do
        if (step > 0) then
        do x_idx=1,num_cells
            idx(1:2) = A_idx2sub(x_idx,2:3)
            el_pot2d(idx(1),idx(2)) = el_pot(x_idx) !u_pre(x_idx) !
            !el_conc(idx(1),idx(2)) = u_pre(x_idx) !
        end do

!!!!! Substituted with save-module, no global large variable needed - was too slow...
!        if (mod(step,save_t) == 0) then
!            s_step = int(step/save_t)
!            dev%time(s_step) = step*t_step
!            dev%g_el_pot(s_step,:,:) = el_pot2d(save_x(1):save_x(2),save_y(1):save_y(2)) !prob_init!
!            dev%g_el_field(s_step,1:2,:,:) = el_field(:,save_x(1):save_x(2),save_y(1):save_y(2))
!            dev%g_el_conc(s_step,:,:) = el_conc(save_x(1):save_x(2),save_y(1):save_y(2))
!        end if
        end if
        end if
    end subroutine solve_poisson_hips_2

!    subroutine media_2()
!        integer                                     :: j, s_step
!        real(kind=s),dimension(:,:),allocatable       :: mean_energy_tot ! mean_vel, mean_energy
!        real(kind=s),dimension(:,:,:),allocatable   :: mean_vel_tot
!        integer,dimension(:,:),allocatable            :: num
!        !integer                                     :: x_idx, y_idx
!
!        if (mod(step,save_t) == 0) then
!        !allocate(mean_vel(max_carriers))
!        !allocate(mean_energy(max_carriers))
!        allocate(mean_vel_tot(2,num_nodes(1),num_nodes(2)))
!        allocate(mean_energy_tot(num_nodes(1),num_nodes(2)))
!        allocate(num(num_nodes(1),num_nodes(2)))
!        mean_vel_tot = 0.0_s
!        mean_energy_tot = 0.0_s
!        num = 0
!
!        !$OMP PARALLEL DO REDUCTION(+:mean_energy_tot,num,mean_vel_tot)
!        do j=1,max_carriers
!            if ((c(j)%valley > 0) .and. (c(j)%op == 0)) then
!                call get_r_idx_2(c(j)%r,c(j)%r_idx)
!                ! compute mean energy in mesh cell
!                mean_energy_tot(c(j)%r_idx(1),c(j)%r_idx(2)) = mean_energy_tot(c(j)%r_idx(1),c(j)%r_idx(2))+c(j)%energy
!                num(c(j)%r_idx(1),c(j)%r_idx(2)) = num(c(j)%r_idx(1),c(j)%r_idx(2))+1
!                ! mean velocity (vector)
!                mean_vel_tot(1,c(j)%r_idx(1),c(j)%r_idx(2)) = mean_vel_tot(1,c(j)%r_idx(1),c(j)%r_idx(2))+c(j)%vel(1)
!                mean_vel_tot(2,c(j)%r_idx(1),c(j)%r_idx(2)) = mean_vel_tot(2,c(j)%r_idx(1),c(j)%r_idx(2))+c(j)%vel(2)
!            end if
!        end do
!        !$OMP END PARALLEL DO
!        !do x_idx=1,num_nodes(1)
!        !    do y_idx=1,num_nodes(2)
!        !        if (num(x_idx,y_idx)== 0) then
!        !            num(x_idx,y_idx) = 1e7
!        !        end if
!        !    end do
!        !end do
!        !dev%g_vel(step,:) = mean_vel_tot/num
!        !if (mod(step,save_t) == 0) then
!            s_step = int(step/save_t)
!            dev%g_energy(s_step,:,:) = mean_energy_tot(save_x(1):save_x(2),save_y(1):save_y(2))/&
!                                   & num(save_x(1):save_x(2),save_y(1):save_y(2))
!            dev%g_vel(s_step,1,:,:) = mean_vel_tot(1,save_x(1):save_x(2),save_y(1):save_y(2))/&
!                                   & num(save_x(1):save_x(2),save_y(1):save_y(2))
!            dev%g_vel(s_step,2,:,:) = mean_vel_tot(2,save_x(1):save_x(2),save_y(1):save_y(2))/&
!                                   & num(save_x(1):save_x(2),save_y(1):save_y(2))
!            dev%g_num(s_step,:,:) = num(save_x(1):save_x(2),save_y(1):save_y(2))
!                                   !print *, particle_out(:,1), particle_out(:,2)
!            dev%particle_flow_out(s_step,:) = sum(particle_out,1)
!            dev%particle_flow_in(s_step,:) = particle_in
!            particle_out = 0
!            particle_in = 0
!        end if
!
!        !print *, 'media_2 called', mean_energy_tot(65,15)
!        ! Computes the mean velocity, needed for ohmic contacts
!        !allocate(mean_vel(2,num_carriers+10000))
!        !mean_vel = 0.0_s
!        !where ((c(:)%valley .gt. 0) .and. (c(:)%r(1) .lt. 20.0_s)) mean_vel(1,:) = c(:)%vel(1)
!        !where ((c(:)%valley .gt. 0) .and. (c(:)%r(1) .lt. (dev_x-20.0_s))) mean_vel(2,:) = c(:)%vel(1)
!        !mean_vel_ohmic(1) = sum(mean_vel(1,:), mean_vel(1,:)>0)/(max(1,count(mean_vel(1,:)>0)))
!        !mean_vel_ohmic(2) = sum(mean_vel(2,:), mean_vel(2,:)>0)/(max(1,count(mean_vel(2,:)>0)))
!        !print *, mean_vel_ohmic
!    end subroutine media_2

    subroutine save_tun_pot()
        integer             :: i

        do i=1,n_tunnel
            !print *, i, size(el_pot2d(tun(i)%t_x_idx(1):tun(i)%t_x_idx(2), &
            !& tun(i)%t_y_idx(1):tun(i)%t_y_idx(2)),2)
            !print *, tun(i)%elpot_save_idx
            tun(i)%elpot_save(tun(i)%elpot_save_idx,:,:) = el_pot2d(tun(i)%t_x_idx(1):tun(i)%t_x_idx(2), &
            & tun(i)%t_y_idx(1):tun(i)%t_y_idx(2))
            if (tun(i)%elpot_save_idx >= tun(i)%steps_compute_tunnel) then !size(tun(i)%elpot_save(:,1,1))) then
                tun(i)%elpot_save_idx = 1
            else
                tun(i)%elpot_save_idx = tun(i)%elpot_save_idx+1
            end if
        end do
    end subroutine save_tun_pot

    subroutine save_pauli()
        integer         :: i, j, mat, x_idx, y_idx
        real(kind=s)    :: E_tmp, l_nonp, l_me
        real(kind=s),dimension(3)   :: k_tmp

        if (mod(step,p%steps_save_pep) == 0) then
            print *, 'Save pep', step, step*t_step
        do i=1,max_carriers
            if ((c(i)%valley > 0) .and. (c(i)%op == 0)) then
            p%num_particles(c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2)) = &
                    & p%num_particles(c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2))+1
            p%step_n_conc(p%step,c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2)) = &
                    & p%step_n_conc(p%step,c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2))+1.0_s
            p%step_k_drift(p%step,c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2),:) = &
                & p%step_k_drift(p%step,c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2),:)+c(i)%k
           !p%test_energy(p%step,c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2)) = &
           !     & p%test_energy(p%step,c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2))+c(i)%energy
            end if
        end do
        !print *, 'k 1', p%step_k_drift(p%step,1,1,10,:)
        do i=1,3
            do j=1,3
            do y_idx=1,num_nodes(2)
                do x_idx=1,num_nodes(1)
            p%step_k_drift(p%step,i,x_idx,y_idx,j) = p%step_k_drift(p%step,i,x_idx,y_idx,j)/dble(p%num_particles(i,x_idx,y_idx))
                end do
            end do
            end do
        end do
        !print *, p%step_k_drift(p%step,1,1,10,:)

        do i=1,max_carriers
            if ((c(i)%valley > 0) .and. (c(i)%op == 0)) then
            mat = m_material(m_region(c(i)%r_idx(1),c(i)%r_idx(2)))
            l_nonp = nonp(mat,c(i)%valley)
            l_me = me(mat,c(i)%valley)
            k_tmp = c(i)%k-p%step_k_drift(p%step,c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2),:)
            !k_tmp = c(i)%k-p%k_drift(c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2),:)
            if ((c(i)%valley == X) .or. (c(i)%valley == L)) then
            E_tmp = -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(k_tmp,k_tmp)/q)/ &
                            & (2.0_s*m0*l_nonp))
            else
            E_tmp = -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(k_tmp,k_tmp)/q)/ &
                            & (2.0_s*l_me*l_nonp))
            end if
            p%step_E_avg(p%step,c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2)) = &
                p%step_E_avg(p%step,c(i)%valley,c(i)%r_idx(1),c(i)%r_idx(2))+E_tmp
            end if
        end do
        !print *, p%step_E_avg(p%step,1,:,10)
        !print *, p%num_particles(1,:,10)
        do i=1,3
            p%step_E_avg(p%step,i,:,:) = p%step_E_avg(p%step,i,:,:)/dble(p%num_particles(i,:,:))
        end do
        !print *, p%step_E_avg(p%step,1,:,10)
        !print *, p%num_particles(1,30,77), p%step_k_drift(p%step,1,30,77,:), p%step_E_avg(p%step,1,30,77)
        p%num_particles = 0
        p%step = p%step+1
        end if
    end subroutine save_pauli

    subroutine assign_charge_2()
        integer                         :: k, i, j, x_idx, y_idx
        real(kind=s)                    :: x, y!, t_ohmic
        !integer                         :: countf, counti, count_rate

            !call system_clock(counti,count_rate)
        el_charge = 0.0_s
        el_conc = 0.0_s
        !t_ohmic = 0.0_s
        !$omp parallel do reduction(+:el_charge)
        do k=1,max_carriers
            if ((c(k)%valley > 0) .and. (c(k)%op == 0)) then
!                x = c(k)%r(1)/m_dx
!                y = c(k)%r(2)/m_dy
!                i = int(x+1)
!                j = int(y+1)
                i = int((c(k)%r(1)+m_dx_half)/m_dx)+1
                j = int((c(k)%r(2)+m_dy_half)/m_dy)+1
!                if (c(k)%r(2) > 399.0_s) then
!                    print *, i,j, c(k)%r
!                end if
                !if ((i > 351) .or. (i<1)) then
                !    print *, c(k)%loc_e_field, c(k)%r(1:2), int(m_region(c(k)%r_idx(1),c(k)%r_idx(2))), c(k)%valley, c(k)%mat
                !    print *, c(k)%cidx, c(k)%r(1:2),c(k)%r_idx,c(k)%energy, c(k)%k
                !    print *, c(k)%cidx, c(k)%last_scat_idx(1), c(k)%last_e(1), c(k)%last_k(1,:)
                !   print *, c(k)%cidx, c(k)%last_scat_idx(2), c(k)%last_e(2), c(k)%last_k(2,:)
                !end if
                if ((i < 1) .or. (i > num_nodes(1)) .or. (j < 1) .or. (j > num_nodes(2))) then
                    print *, 'particle left dev'
                    c(k)%valley = 0
                else
                el_charge(i,j) = el_charge(i,j)+1
!                el_charge(i,j) = el_charge(i,j)+(1-(x-(i-1)))*(1-(y-(j-1)))
!                el_charge(i,j+1) = el_charge(i,j+1)+(1-(x-(i-1)))*(y-(j-1))
!                el_charge(i+1,j+1) = el_charge(i+1,j+1)+(x-(i-1))*(y-(j-1))
!                el_charge(i+1,j) = el_charge(i+1,j)+(x-(i-1))*(1-(y-(j-1)))
                end if
                !print *, el_charge(i,j)+el_charge(i,j+1)+el_charge(i+1,j+1)+el_charge(i+1,j)
            end if
        end do
        !$omp end parallel do
            !    call system_clock(countf)
            !print *, 'assign charge int', real(countf-counti)/real(count_rate)
        do y_idx=1,num_nodes(2)
            do x_idx=1,num_nodes(1)
                el_charge(x_idx,y_idx) = el_charge(x_idx,y_idx)*p_charge
                el_conc(x_idx,y_idx) = el_charge(x_idx,y_idx)/(m_lx(x_idx)*m_ly(y_idx))
!
!                ! calculate factor for keeping ohmic contacts
!                if (step == 120) then
!                    !t_ohmic = 1! mod(step,100)
!                end if
 !               if ((contact_type(x_idx,y_idx) == OHMIC) .and. (t_ohmic >= 0.5_s) .and. (donor(x_idx,y_idx) > 0.0_s)) then
  !                  conc_tmp = sum(dev%g_el_conc(30:120,x_idx,y_idx))/90
   !                 ohmic_f(x_idx,y_idx) = 1.0_s-(conc_tmp-donor(x_idx,y_idx))/conc_tmp
!                    print *, conc_tmp, ohmic_f(x_idx,y_idx)
                    !if ((ohmic_f(x_idx,y_idx) <= 0.0_s) .or. (ohmic_f(x_idx,y_idx) > 1.0_s)) then
                    !    ohmic_f(x_idx,y_idx) = 1.0_s
                    !end if
                    !print *, conc_tmp, ohmic_f(x_idx,y_idx)
    !            end if
            end do
        end do

        if (step*t_step >= 0.05_s) then
            call test_conc_2()
        end if
       !call floating_conc()
    end subroutine assign_charge_2

    subroutine asdf()
!    subroutine init_poisson_2()
!        real(kind=s),dimension(:),allocatable  :: A_tmp!, num_cells
!        integer                                :: i, j, num_region, A_offs, A_idx, x_idx, y_idx, s_x, e_x, s_y, e_y
!        real(kind=s)                           :: dl, dr, dd, du
!        !integer                                :: max_nodes
!
!        num_region = size(m(1,:))
!        !allocate(num_cells(num_region))
!
!        !allocate(dev%particle_flow_out(n_step,2))
!        !dev%particle_flow_out = 0
!
!        ! Allocate space for electric field, potential and charge
!        allocate(el_field(2,sum(num_cells)))
!        allocate(el_pot(sum(num_cells)))
!        !allocate(el_charge(sum(num_nodes)))
!        !allocate(el_conc(sum(num_nodes)))
!
!        !allocate(dev%g_el_conc(n_step,sum(n    um_nodes)))
!        allocate(dev%g_el_pot(n_step,sum(num_cells)))
!        allocate(dev%g_el_field(n_step,2,sum(num_cells)))
!
!        !dev%g_el_conc = 0.0_s
!        dev%g_el_pot = 0.0_s
!        dev%g_el_field = 0.0_s
!
!        !el_field = 0.0_s !-1.0e4_s/1e9_s
!
!        !allocate(dev%g_vel(n_step,sum(num_nodes)))
!        !dev%g_vel = 0.0_s
!
!        !allocate(dev%g_energy(n_step,sum(num_nodes)))
!        !dev%g_energy = 0.0_s
!
!        ! Calculate superparticle charge
!        p_charge = 0.0_s
!        do i=1,num_region
!            p_charge = p_charge+donor_density(i)*((m(2,i)-m(1,i))*(m(4,i)-m(3,i)))/num_carriers
!        end do
!        print *, 'Superparticle charge: (e/nm):', p_charge
!
!    !print *, num_cells(1)
!        allocate(A(sum(num_cells),sum(num_cells)))
!        allocate(A_tmp(sum(num_cells)))
!        allocate(A_sub2idx(num_region,maxval(num_nodes(:,1)),maxval(num_nodes(:,2))))
!        allocate(A_idx2sub(sum(num_cells),3))
!        A_idx2sub = 0
!        A_sub2idx = 0
!        A = 0.0_s
!        A_tmp = 0.0_s
!
!        ! Fill matrix
!        A_offs = 0
!        i = 1
!            do y_idx=1,num_nodes(1,2)
!                do x_idx=1,num_nodes(1,1)
!                    A_idx = A_offs+((y_idx-1)*num_nodes(i,1)+x_idx)
!                    A_sub2idx(i,x_idx,y_idx) = A_idx
!                    A_idx2sub(A_idx,1:3) = (/i, x_idx, y_idx/)
!                    dl = mesh_dist_lr(i,1,x_idx)
!                    dr = mesh_dist_lr(i,2,x_idx)
!                    dd = mesh_dist_ud(i,1,y_idx)
!                    du = mesh_dist_ud(i,2,y_idx)
!
!                    ! only compute the matrix for inner nodes, rest has to be computed using boundary conditions
!                    if (((node_coord(i,1,x_idx) > 0.0_s+1e-5_s) .and. (node_coord(i,1,x_idx) < dev_x-1e-5_s)) .and. &
!                     &  ((node_coord(i,2,y_idx) > 0.0_s+1e-5_s) .and. (node_coord(i,2,y_idx) < dev_y-1e-5_s))) then
!                        !print *, x_idx, y_idx, (y_idx-1)*num_nodes(i,1)+x_idx
!                        A_tmp = 0.0_s
!                        A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
!                        A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
!                        A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
!                        A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
!                        A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
!                        A(A_idx,:) = A_tmp
!                    else
!                        if ((contact_type(i,x_idx,y_idx) == SCHOTTKY) .or. (contact_type(i,x_idx,y_idx) == OHMIC)) then
!                            A(A_idx,A_idx) = 1.0_s
!                        elseif (contact_type(i,x_idx,y_idx) == INSULATOR) then
!                            !print *, 'boundary:', i, x_idx, y_idx
!
!                            if ((x_idx == num_nodes(1,1)) .and. (y_idx == num_nodes(1,2))) then
!                            !    print *, 'neumann', y_idx, dl, dr, dd, du
!                                A_tmp = 0.0_s
!                                A_tmp(A_idx-1) = 4.0_s/(dl*(dl+dr))
!                                A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !-2.0_s/(du*dd)!
!                                !!A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
!                                A_tmp(A_idx-num_nodes(1,1)) = 4.0/(dd*(du+dd))
!                                !A_tmp(A_idx+num_nodes(1,1)) = 2.0/(du*(du+dd))
!                            elseif ((x_idx == 1) .and. (y_idx == num_nodes(1,2))) then
!                            !    print *, 'neumann', y_idx, dl, dr, dd, du
!                                A_tmp = 0.0_s
!                                !A_tmp(A_idx-1) = 4.0_s/(dl*(dl+dr))
!                                A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !-2.0_s/(du*dd)!
!                                A_tmp(A_idx+1) = 4.0_s/(dr*(dl+dr))
!                                A_tmp(A_idx-num_nodes(1,1)) = 4.0/(dd*(du+dd))
!                                !A_tmp(A_idx+num_nodes(1,1)) = 4.0/(du*(du+dd))
!                            elseif ((x_idx == num_nodes(1,1)) .and. (y_idx == 1)) then
!                            !    print *, 'neumann', y_idx, dl, dr, dd, du
!                                A_tmp = 0.0_s
!                                A_tmp(A_idx-1) = 4.0_s/(dl*(dl+dr))
!                                A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !-2.0_s/(du*dd)!
!                                !A_tmp(A_idx+1) = 4.0_s/(dr*(dl+dr))
!                                !A_tmp(A_idx-num_nodes(1,1)) = 4.0/(dd*(du+dd))
!                                A_tmp(A_idx+num_nodes(1,1)) = 4.0/(du*(du+dd))
!                            elseif ((x_idx == 1) .and. (y_idx == 1)) then
!                            !    print *, 'neumann', y_idx, dl, dr, dd, du
!                                A_tmp = 0.0_s
!                                !A_tmp(A_idx-1) = 4.0_s/(dl*(dl+dr))
!                                A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !-2.0_s/(du*dd)!
!                                A_tmp(A_idx+1) = 4.0_s/(dr*(dl+dr))
!                                !A_tmp(A_idx-num_nodes(1,1)) = 4.0/(dd*(du+dd))
!                                A_tmp(A_idx+num_nodes(1,1)) = 4.0/(du*(du+dd))
!                            elseif (y_idx == 1) then
!                                A_tmp = 0.0_s
!                                A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
!                                A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !1.0_s/200.0_s ! -2.0_s/(du*dd) !
!                                A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
!                                A_tmp(A_idx+num_nodes(1,1)) = 4.0/(du*(du+dd)) !-1.0_s/200.0_s
!                            elseif (y_idx == num_nodes(1,2)) then
!                            !    print *, 'neumann', y_idx, dl, dr, dd, du
!                                A_tmp = 0.0_s
!                                A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
!                                A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !-2.0_s/(du*dd)!
!                                A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
!                                A_tmp(A_idx-num_nodes(1,1)) = 4.0/(dd*(du+dd))
!                            elseif (x_idx == 1) then
!                                A_tmp = 0.0_s
!                                !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
!                                A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !1.0_s/200.0_s ! -2.0_s/(du*dd) !
!                                A_tmp(A_idx+1) = 4.0_s/(dr*(dl+dr))
!                                A_tmp(A_idx-num_nodes(1,1)) = 2.0/(dd*(du+dd))
!                                A_tmp(A_idx+num_nodes(1,1)) = 2.0/(du*(du+dd))
!                            elseif (x_idx == num_nodes(1,1)) then
!                            !    print *, 'neumann', y_idx, dl, dr, dd, du
!                                A_tmp = 0.0_s
!                                A_tmp(A_idx-1) = 4.0_s/(dl*(dl+dr))
!                                A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd) !-2.0_s/(du*dd)!
!                                !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
!                                A_tmp(A_idx-num_nodes(1,1)) = 2.0/(dd*(du+dd))
!                                A_tmp(A_idx+num_nodes(1,1)) = 2.0/(du*(du+dd))
!                            end if
!                            A(A_idx,:) = A_tmp
!                        else
!                            print *, 'something went wrong...'
!                        end if
!                    end if
!                end do
!            end do
!        !end do
!!        do i=1,num_region
!!            if (i > 1) then
!!                A_offs = A_offs+(num_nodes(i-1,1)*num_nodes(i-1,2))
!!            end if
!!            do x_idx=1,num_nodes(i,1)
!!                do y_idx=1,num_nodes(i,2)
!!                    A_idx = A_offs+((y_idx-1)*num_nodes(i,1)+x_idx)
!!                    A_sub2idx(i,x_idx,y_idx) = A_idx
!!                    A_idx2sub(A_idx,1:3) = (/i, x_idx, y_idx/)
!!                    if (i == 1) then
!!                        print *, A_idx, A_sub2idx(i,x_idx,y_idx), A_idx2sub(A_idx,1:3)
!!                        end if
!!                    ! only compute the matrix for inner nodes, rest has to be computed using boundary conditions
!!                    if (((node_coord(i,1,x_idx) > 0.0_s+1e-5_s) .and. (node_coord(i,1,x_idx) < dev_x-1e-5_s)) .and. &
!!                     &  ((node_coord(i,2,y_idx) > 0.0_s+1e-5_s) .and. (node_coord(i,2,y_idx) < dev_y-1e-5_s))) then
!!                        !print *, x_idx, y_idx, (y_idx-1)*num_nodes(i,1)+x_idx
!!                        A_tmp = 0.0_s
!!                        dl = mesh_dist_lr(i,1,x_idx)
!!                        dr = mesh_dist_lr(i,2,x_idx)
!!                        dd = mesh_dist_ud(i,1,y_idx)
!!                        du = mesh_dist_ud(i,2,y_idx)
!!                        A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
!!                        A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
!!                        A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
!!                        !A_tmp(A_idx-num_nodes(i,1)) = 2.0_s/(dd*(du+dd))
!!                        !A_tmp(A_idx+num_nodes(i,1)) = 2.0_s/(du*(du+dd))
!!                        A_tmp(A_idx-70) = 2.0_s/(dd*(du+dd))
!!                        A_tmp(A_idx+70) = 2.0_s/(du*(du+dd))
!!                        A(A_idx,:) = A_tmp
!!                        !print *, dl, dr, dd, du
!!                        !print *, x_idx, y_idx, A_idx, node_coord(i,1,x_idx), node_coord(i,2,y_idx)
!!                        !print *, x_idx, y_idx, A_idx, dl
!!                    else
!!                        if ((contact_type(i,x_idx,y_idx) == SCHOTTKY) .or. (contact_type(i,x_idx,y_idx) == OHMIC)) then
!!                            A(A_idx,A_idx) = 1.0_s
!!                            !print *, 'boundary:', i, x_idx, y_idx
!!                        end if
!!                    end if
!!                end do
!!            end do
!!        end do
!        !print *, 'A', A(58,57), A(58,58), A(58,59)
!
!        ! Precompute known vector u
!        allocate(u_pre(sum(num_cells)))
!        u_pre = 0.0_s
!        i = 1
!        do y_idx=1,num_nodes(1,2)
!            do x_idx=1,num_nodes(1,1)
!                if ((contact_type(i,x_idx,y_idx) == SCHOTTKY) .or. (contact_type(i,x_idx,y_idx) == OHMIC)) then
!                    u_pre(A_sub2idx(i,x_idx,y_idx)) = contact_pot(i,x_idx,y_idx)
!                    !print *, i, x_idx, y_idx, A_sub2idx(i,x_idx,y_idx), contact_pot(i,x_idx,y_idx)
!                end if
!            end do
!        end do
!        !do i=1,num_region
!        !    do x_idx=1,num_nodes(i,1)
!        !        do y_idx=1,num_nodes(i,2)
!        !            if ((contact_type(i,x_idx,y_idx) == SCHOTTKY) .or. (contact_type(i,x_idx,y_idx) == OHMIC)) then
!        !                    u_pre(A_sub2idx(i,x_idx,y_idx)) = contact_pot(i,x_idx,y_idx)
!        !                    !print *, i, x_idx, y_idx, A_sub2idx(i,x_idx,y_idx), contact_pot(i,x_idx,y_idx)
!        !            end if
!        !        end do
!        !    end do
!        !end do
!
!        !open (30, file='/tmp/A.txt', access='APPEND', status='REPLACE')
!        !write (30,*) A
!        print *, 'n cells', num_cells, sum(num_cells)
!        print *, A_sub2idx(1,25,6)
!        print *, 'idx2sub', A_idx2sub(150,1:3)
!    end subroutine init_poisson_2
!
!    subroutine solve_poisson_2()
!        real(kind=8),dimension(:),allocatable    :: ipiv, rho
!        real(kind=8),dimension(:),allocatable    :: pot_new
!        integer                                  :: info, i, l, x_idx, y_idx
!        real(kind=s),dimension(:,:),allocatable  :: M
!
!        l = sum(num_cells)
!        allocate(pot_new(l))
!        allocate(rho(l))
!        allocate(ipiv(l))
!        allocate(M(l,l))
!        M = (A)
!        pot_new = 0.0_s
!
!        ! Set up RHS
!        !rho = (-el_conc+donor)
!        !pot_new = -rho/(eps_static(INGAAS))
!
!        ! Set boundary conditions
!        pot_new = u_pre
!        pot_new(A_sub2idx(1,51,51)) = (-1.0_s/25.0_s)/(eps_static(INGAAS))
!        pot_new(A_sub2idx(1,51,c_pos)) = (-1.0_s/25.0_s)/(eps_static(INGAAS))
!
!        print *, 'charge 1', A_idx2sub(5570,:), 'charge 2', A_idx2sub(5575,:)
!        !open (30, file='/tmp/u.txt', access='APPEND', status='REPLACE')
!        !write (30,*) pot_new
!        !write (30,*) u_pre
!        !write (30,*) contact_type(2,:,:)
!        !open (30, file='/tmp/rho.txt', access='APPEND', status='REPLACE', iostat=st2)
!        !write (30,*) -rho/eps_static(INGAAS)
!
!        call dgesv(l, 1, M, l, ipiv, pot_new, l, info)
!        print *, 'poisson info', info
!
!        el_pot = pot_new
!        ! E-Field in x-direction
!        do x_idx=2,num_nodes(1,1)-1
!            do y_idx=1,num_nodes(1,2)
!                el_field(1,A_sub2idx(1,x_idx,y_idx)) = -(el_pot(A_sub2idx(1,x_idx+1,y_idx))-el_pot(A_sub2idx(1,x_idx-1,y_idx)))/ &
!                                                    &  (node_coord(1,1,x_idx+1)-node_coord(1,1,x_idx-1))
!            end do
!        end do
!        ! E-Field in y-direction
!        do y_idx=2,num_nodes(1,2)-1
!            do x_idx=1,num_nodes(1,1)
!                el_field(2,A_sub2idx(1,x_idx,y_idx)) = -(el_pot(A_sub2idx(1,x_idx,y_idx+1))-el_pot(A_sub2idx(1,x_idx,y_idx-1)))/ &
!                                                    &  (node_coord(1,2,y_idx+1)-node_coord(1,2,y_idx-1))
!            end do
!        end do
!        do x_idx=1,num_nodes(1,1)
!            el_field(2,A_sub2idx(1,x_idx,1)) = el_field(2,A_sub2idx(1,x_idx,2))
!            el_field(2,A_sub2idx(1,x_idx,num_nodes(1,2))) = el_field(2,A_sub2idx(1,x_idx,num_nodes(1,2)-1))
!        end do
!        do y_idx=1,num_nodes(1,2)
!            el_field(1,A_sub2idx(1,1,y_idx)) = el_field(1,A_sub2idx(1,2,y_idx))
!            el_field(1,A_sub2idx(1,num_nodes(1,1),y_idx)) = el_field(1,A_sub2idx(1,num_nodes(1,1)-1,y_idx))
!        end do
!
!        ! Save output in variables
!        !dev%g_el_conc(step,:) = el_conc
!        dev%g_el_pot(step,:) = el_pot
!        dev%g_el_field(step,1:2,:,:) = el_field
!    end subroutine solve_poisson_2
    end subroutine asdf
#endif

#if DIM == 3
    subroutine init_poisson_hips_3()
        integer                                :: x_idx, y_idx, z_idx, nnz_idx, n_idx, A_idx, idnbr, ierr, id, ntasks, proc_id
        integer                                :: num_region, i, orientation, coord1, coord2, num_xy
        real(kind=s),dimension(:),allocatable  :: aval
        integer,dimension(:),allocatable       :: ia, ja
        real(kind=s)                           :: dl, dr, du, dd, db, df, prec

        call MPI_Comm_rank(MPI_COMM_WORLD, proc_id , ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, ntasks , ierr)

        ! Allocate space for electric field, potential and charge
        allocate(el_field(3,num_nodes(1),num_nodes(2),num_nodes(3)))
        allocate(el_pot(num_cells(1)))
        allocate(el_charge(num_nodes(1),num_nodes(2),num_nodes(3)))
        allocate(el_conc(num_nodes(1),num_nodes(2),num_nodes(3)))

        allocate(dev%g_el_conc(n_step,num_nodes(1),num_nodes(2),num_nodes(3)))
        allocate(dev%g_el_pot(n_step,num_nodes(1),num_nodes(2),num_nodes(3)))
        allocate(dev%g_el_field(n_step,3,num_nodes(1),num_nodes(2),num_nodes(3)))

        allocate(A_idx2sub(sum(num_cells),3))
        allocate(A_sub2idx(num_nodes(1),num_nodes(2),num_nodes(3)))

        !dev%g_el_conc = 0.0_s
        !dev%g_el_pot = 0.0_s
        dev%g_el_field = 0.0_s

        !allocate(dev%g_vel(n_step,sum(num_nodes)))
        !dev%g_vel = 0.0_s

        !allocate(dev%g_energy(n_step,sum(num_nodes)))
        !dev%g_energy = 0.0_s

        ! Calculate superparticle charge
        num_region = size(m(1,:))
        p_charge = 0.0_s
        do i=1,num_region
            p_charge = p_charge+donor_density(i)*((m(2,i)-m(1,i))*(m(4,i)-m(3,i))*(m(6,i)-m(5,i)))/num_carriers
        end do
        print *, 'Superparticle charge: (e):', p_charge, num_region


        call MPI_Comm_rank(MPI_COMM_WORLD, proc_id , ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, ntasks , ierr)

        ! We use CSR format
        ! Compute number of non-zero elements and number of unknowns
        n = num_cells(1)
        nnz = 0
        print *, 'num_nodes', num_nodes(1:3)
        do z_idx=1,num_nodes(3)
            do y_idx=1,num_nodes(2)
                do x_idx=1,num_nodes(1)
                    call coord2idx(x_idx,y_idx,z_idx,orientation,coord1,coord2)
                    if (orientation == INSIDE) then
                        nnz = nnz+7
                    else if ((contact_type(orientation,coord1,coord2) == SCHOTTKY) .or. &
                          & (contact_type(orientation,coord1,coord2) == OHMIC)) then
                        nnz = nnz+1
                    else if (contact_type(orientation,coord1,coord2) == INSULATOR) then
                        if (((x_idx == 1) .and. (y_idx == 1) .and. (z_idx == 1)) .or. &
                         & ((x_idx == num_nodes(1)) .and. (y_idx == 1) .and. (z_idx == 1)) .or. &
                         & ((x_idx == 1) .and. (y_idx == num_nodes(2)) .and. (z_idx == 1)) .or. &
                         & ((x_idx == num_nodes(1)) .and. (y_idx == num_nodes(2)) .and. (z_idx == 1)) .or. &
                         & ((x_idx == 1) .and. (y_idx == 1) .and. (z_idx == num_nodes(3))) .or. &
                         & ((x_idx == num_nodes(1)) .and. (y_idx == 1) .and. (z_idx == num_nodes(3))) .or. &
                         & ((x_idx == 1) .and. (y_idx == num_nodes(2)) .and. (z_idx == num_nodes(3))) .or. &
                         & ((x_idx == num_nodes(1)) .and. (y_idx == num_nodes(2)) .and. (z_idx == num_nodes(3)))) then
                            nnz = nnz+4
                        elseif (((y_idx == 1) .and. (z_idx == 1)) .or. &
                         & ((y_idx == num_nodes(2)) .and. (z_idx == 1)) .or. &
                         & ((x_idx == 1) .and. (z_idx == 1)) .or. &
                         & ((x_idx == num_nodes(1)) .and. (z_idx == 1)) .or. &
                         & ((x_idx == 1) .and. (y_idx == 1)) .or. &
                         & ((x_idx == num_nodes(1)) .and. (y_idx == 1)) .or. &
                         & ((x_idx == 1) .and. (y_idx == num_nodes(2))) .or. &
                         & ((x_idx == num_nodes(1)) .and. (y_idx == num_nodes(2))) .or. &
                         & ((y_idx == 1) .and. (z_idx == num_nodes(3))) .or. &
                         & ((y_idx == num_nodes(2)) .and. (z_idx == num_nodes(3))) .or. &
                         & ((x_idx == 1) .and. (z_idx == num_nodes(3))) .or. &
                         & ((x_idx == num_nodes(1)) .and. (z_idx == num_nodes(3)))) then
                            nnz = nnz+5
                        else
                        nnz = nnz+6
                    end if
                    end if
                end do
            end do
        end do
        print *, 'n:', n, 'nnz:', nnz

        ! Fill sparse martrix, CSR format
        allocate(aval(nnz))
        allocate(ja(nnz))
        allocate(ia(n+1))
        allocate(u_pre(n))

        if (proc_id == 0) then

        nnz_idx = 1
        n_idx = 1
        num_xy = num_nodes(1)*num_nodes(2)
        do z_idx=1,num_nodes(3)
            do y_idx=1,num_nodes(2)
                do x_idx=1,num_nodes(1)
                    A_idx = (z_idx-1)*num_nodes(1)*num_nodes(2)+(y_idx-1)*num_nodes(1)+x_idx
                    A_sub2idx(x_idx,y_idx,z_idx) = A_idx
                    !print *, x_idx, y_idx, z_idx, A_idx
                    A_idx2sub(A_idx,1:3) = (/x_idx, y_idx, z_idx/)
                    ! Compute coefficients
                    dl = m_dx(1) !mesh_dist_lr(1,1,x_idx)
                    dr = m_dx(1) !mesh_dist_lr(1,2,x_idx)
                    dd = m_dx(1) !mesh_dist_ud(1,1,y_idx)
                    du = m_dx(1) !mesh_dist_ud(1,2,y_idx)
                    df = m_dx(1)
                    db = m_dx(1)
                    !print *, dr, dl, dd, du
                    call coord2idx(x_idx,y_idx,z_idx,orientation,coord1,coord2)

                    if (orientation == INSIDE) then
                        ia(n_idx) = nnz_idx

                        aval(nnz_idx) = 2.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                        ja(nnz_idx) = A_idx-num_xy
                        nnz_idx = nnz_idx+1

                        aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                        ja(nnz_idx) = A_idx-num_nodes(1)
                        nnz_idx = nnz_idx+1

                        aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                        ja(nnz_idx) = A_idx-1
                        nnz_idx = nnz_idx+1

                        aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                        ja(nnz_idx) = A_idx
                        nnz_idx = nnz_idx+1

                        aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                        ja(nnz_idx) = A_idx+1
                        nnz_idx = nnz_idx+1

                        aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                        ja(nnz_idx) = A_idx+num_nodes(1)
                        nnz_idx = nnz_idx+1

                        aval(nnz_idx) = 2.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                        ja(nnz_idx) = A_idx+num_xy
                        nnz_idx = nnz_idx+1

                    ! Schottky or ohmic contacts, fixed potential
                    else if ((contact_type(orientation,coord1,coord2) == SCHOTTKY) .or. &
                          & (contact_type(orientation,coord1,coord2) == OHMIC)) then
                        aval(nnz_idx) = 1.0_s
                        ja(nnz_idx) = A_idx
                        ia(n_idx) = nnz_idx
                        u_pre(n_idx) = contact_pot(orientation,coord1,coord2)
                        nnz_idx = nnz_idx+1

                    ! Neumann boundary condition, E_n = 0
                    else if (contact_type(orientation,coord1,coord2) == INSULATOR) then
                        if ((x_idx == 1) .and. (y_idx == 1) .and. (z_idx == 1)) then               ! Ecke vorne links unten
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(db*(df+db))
                            ja(nnz_idx) = A_idx+num_xy
                        elseif ((x_idx == num_nodes(1)) .and. (y_idx == 1) .and. (z_idx == 1)) then ! Ecke vorne rechts unten
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(db*(df+db))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == 1) .and. (y_idx == num_nodes(2)) .and. (z_idx == 1)) then ! Ecke vorne links oben
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(db*(df+db))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == num_nodes(1)) .and. (y_idx == num_nodes(2)) .and. (z_idx == 1)) then ! Ecke vorne rechts oben
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == 1) .and. (y_idx == 1) .and. (z_idx == num_nodes(3))) then                       ! Ecke hinten links unten
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(df*(db+df))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == num_nodes(1)) .and. (y_idx == 1) .and. (z_idx == num_nodes(3))) then            ! Ecke hinten rechts unten
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == num_nodes(1)) .and. (y_idx == 1) .and. (z_idx == num_nodes(3))) then            ! Ecke hinten links oben
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == num_nodes(1)) .and. (y_idx == num_nodes(2)) .and. (z_idx == num_nodes(3))) then ! Ecke hinten rechts oben
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                        elseif ((y_idx == 1) .and. (z_idx == 1)) then                       ! Kante, vorne unten
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((y_idx == num_nodes(2)) .and. (z_idx == 1)) then            ! Kante, vorne oben
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == 1) .and. (z_idx == 1)) then                       ! Kante, vorne links
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == num_nodes(1)) .and. (z_idx == 1)) then            ! Kante, vorne rechts
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == 1) .and. (y_idx == 1)) then                       ! Kante, unten links
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == num_nodes(1)) .and. (y_idx == 1)) then            ! Kante, unten rechts
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == 1) .and. (y_idx == num_nodes(2))) then            ! Kante, oben links
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == num_nodes(1)) .and. (y_idx == num_nodes(2))) then ! Kante, oben rechts
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif ((y_idx == 1) .and. (z_idx == num_nodes(3))) then            ! Kante, hinten unten
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                        elseif ((y_idx == num_nodes(2)) .and. (z_idx == num_nodes(3))) then ! Kante, hinten oben
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == 1) .and. (z_idx == num_nodes(3))) then            ! Kante, hinten links
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                        elseif ((x_idx == num_nodes(1)) .and. (z_idx == num_nodes(3))) then ! Kante, hinten rechts
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                        elseif (x_idx == 1) then
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif (x_idx == num_nodes(1)) then
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif (y_idx == 1) then
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif (y_idx == num_nodes(2)) then
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif (z_idx == 1) then
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 4.0_s/(db*(df+db))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_xy
                            nnz_idx = nnz_idx+1
                        elseif (z_idx == num_nodes(3)) then
                            ia(n_idx) = nnz_idx
                            aval(nnz_idx) = 4.0_s/(df*(db+df))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_xy
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dd*(du+dd))              !A_tmp(A_idx-num_nodes(1,1)) = 2.0_s/(dd*(du+dd))
                            ja(nnz_idx) = A_idx-num_nodes(1)
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dl*(dl+dr))              !A_tmp(A_idx-1) = 2.0_s/(dl*(dl+dr))
                            ja(nnz_idx) = A_idx-1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)-2.0_s/(df*db)    !A_tmp(A_idx) = -2.0_s/(dl*dr)-2.0_s/(du*dd)
                            ja(nnz_idx) = A_idx
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(dr*(dl+dr))              !A_tmp(A_idx+1) = 2.0_s/(dr*(dl+dr))
                            ja(nnz_idx) = A_idx+1
                            nnz_idx = nnz_idx+1
                            aval(nnz_idx) = 2.0_s/(du*(du+dd))              !A_tmp(A_idx+num_nodes(1,1)) = 2.0_s/(du*(du+dd))
                            ja(nnz_idx) = A_idx+num_nodes(1)
                            nnz_idx = nnz_idx+1
                        end if
                    end if
                    n_idx = n_idx+1
                end do
            end do
        end do
        ia(n+1) = nnz+1
        end if

        ! Init the HIPS package
        idnbr = 1
        call HIPS_INITIALIZE(idnbr, ierr)
        id = 0

        ! Set Options
        prec = 1e-6
        call HIPS_SetDefaultOptions(id, HIPS_ITERATIVE, ierr) !HIPS_ITERATIVE
        !call HIPS_SetOptionINT(id, HIPS_DOMSIZE, 100, ierr)
        call HIPS_SetOptionINT(id, HIPS_SYMMETRIC, 0, ierr)
        call HIPS_SetOptionREAL(id, HIPS_PREC, prec, ierr)
        call HIPS_SetOptionINT(id, HIPS_LOCALLY, 0, ierr)
        call HIPS_SetOptionINT(id, HIPS_ITMAX, 100, ierr)
        call HIPS_SetOptionINT(id, HIPS_KRYLOV_RESTART, 50, ierr)
        call HIPS_SetOptionINT(id, HIPS_VERBOSE, 2, ierr)
        call HIPS_SetOptionINT(id, HIPS_DOMNBR, ntasks, ierr)
        call HIPS_SetOptionINT(id, HIPS_CHECK_GRAPH, 1, ierr)
        call HIPS_SetOptionINT(id, HIPS_CHECK_MATRIX, 1, ierr)
         !call HIPS_SetOptionREAL(id, HIPS_DROPTOL0, 0.001_s, ierr)
         !call HIPS_SetOptionREAL(id, HIPS_DROPTOL1, 0.001_s, ierr)
         !call HIPS_SetOptionREAL(id, HIPS_DROPTOLE, 0.001_s, ierr)

        call HIPS_GRAPHGLOBALCSR(id, n, ia, ja, 0, ierr)
        call HIPS_MATRIXGLOBALCSR(id, n, ia, ja, aval,  0, HIPS_ASSEMBLY_OVW, 0, ierr)
    end subroutine init_poisson_hips_3

    subroutine solve_poisson_hips_3()
        integer                                :: proc_id, ierror, x_idx, y_idx, z_idx, orientation, coord1, coord2!, step
        integer,dimension(3)                   :: idx
        real(kind=s)                           :: rho

        !step = 1

        call MPI_Comm_rank(MPI_COMM_WORLD, proc_id , ierror)
        if (proc_id == 0) then
            do x_idx=1,num_nodes(1)
                do y_idx=1,num_nodes(2)
                    do z_idx=1,num_nodes(3)
                        rho = -el_conc(x_idx,y_idx,z_idx)+donor(m_region(x_idx,y_idx,z_idx))
                        call coord2idx(x_idx,y_idx,z_idx,orientation,coord1,coord2)
                        if (orientation == INSIDE) then
                            u_pre(A_sub2idx(x_idx,y_idx,z_idx)) = -rho/(eps_static(INGAAS))
                        elseif ((contact_type(orientation,coord1,coord2) == SCHOTTKY) .or. &
                              & (contact_type(orientation,coord1,coord2) == OHMIC)) then
                            u_pre(A_sub2idx(x_idx,y_idx,z_idx)) = contact_pot(orientation,coord1,coord2)
                        else
                            u_pre(A_sub2idx(x_idx,y_idx,z_idx)) = rho/2.0_s
                        end if
                    end do
                end do
            end do
            !call coord2idx(51,51,51,orientation,coord1,coord2)
            !print *, orientation, coord1, coord2
            !print *, u_pre(A_sub2idx(51,51,51)), el_conc(51,51,51)
        end if
        call HIPS_SETGLOBALRHS(0, u_pre, 0, 0, ierror)
        call HIPS_EXITONERROR(ierror)
        call HIPS_GETGLOBALSOLUTION(0, el_pot, 0, ierror)
        call HIPS_EXITONERROR(ierror)

        if (proc_id == 0) then
        do x_idx=1,num_cells(1)
            idx(1:3) = A_idx2sub(x_idx,1:3)
            dev%g_el_pot(step,idx(1),idx(2),idx(3)) = el_pot(x_idx)
        end do
        ! E-Field in x-direction
        do x_idx=2,num_nodes(1)-1
            do y_idx=1,num_nodes(2)
                do z_idx=1,num_nodes(3)
                    idx(1:3) = (/x_idx,y_idx,z_idx/)
                    el_field(1,x_idx,y_idx,z_idx) = -(dev%g_el_pot(step,idx(1)+1,idx(2),idx(3)) &
                                     & -dev%g_el_pot(step,idx(1)-1,idx(2),idx(3)))/&
                                     & (node_coord(1,x_idx+1)-node_coord(1,x_idx-1))
                end do
            end do
        end do
        ! E-Field in y-direction
        do x_idx=1,num_nodes(1)
            do y_idx=2,num_nodes(2)-1
                do z_idx=1,num_nodes(3)
                    idx(1:3) = (/x_idx,y_idx,z_idx/)
                    el_field(2,x_idx,y_idx,z_idx) = -(dev%g_el_pot(step,idx(1),idx(2)+1,idx(3)) &
                                     & -dev%g_el_pot(step,idx(1),idx(2)-1,idx(3)))/&
                                     & (node_coord(2,y_idx+1)-node_coord(2,y_idx-1))
                end do
            end do
        end do
        ! E-Field in z-direction
        do x_idx=1,num_nodes(1)
            do y_idx=1,num_nodes(2)
                do z_idx=2,num_nodes(3)-1
                    idx(1:3) = (/x_idx,y_idx,z_idx/)
                    el_field(3,x_idx,y_idx,z_idx) = -(dev%g_el_pot(step,idx(1),idx(2),idx(3)+1) &
                                     & -dev%g_el_pot(step,idx(1),idx(2),idx(3)-1))/&
                                     & (node_coord(3,z_idx+1)-node_coord(3,z_idx-1))
                end do
            end do
        end do

        el_field(1,1,1:num_nodes(2),1:num_nodes(3)) = el_field(1,2,1:num_nodes(2),1:num_nodes(3))
        el_field(1,num_nodes(1),1:num_nodes(2),1:num_nodes(3)) = el_field(1,num_nodes(1)-1, &
                            & 1:num_nodes(2),1:num_nodes(3))

        el_field(2,1:num_nodes(1),1,1:num_nodes(3)) = el_field(2,1:num_nodes(1),2,1:num_nodes(3))
        el_field(2,1:num_nodes(1),num_nodes(2),1:num_nodes(3)) = el_field(2,1:num_nodes(1), &
                            & num_nodes(2)-1,1:num_nodes(3))

        el_field(3,1:num_nodes(1),1:num_nodes(2),1) = el_field(3,1:num_nodes(1),1:num_nodes(2),2)
        el_field(3,1:num_nodes(1),1:num_nodes(2),num_nodes(3)) = el_field(2,1:num_nodes(1), &
                            & 1:num_nodes(2),num_nodes(3)-1)

        dev%g_el_conc(step,:,:,:) = el_conc
        dev%g_el_field(step,:,:,:,:) = el_field
        end if
    end subroutine solve_poisson_hips_3

    subroutine assign_charge_3()
        integer                         :: m, k, i, j, x_idx, y_idx, z_idx, orientation, coord1, coord2
        real(kind=s)                    :: x, y, z, conc_tmp, t_ohmic

        el_charge = 0.0_s
        el_conc = 0.0_s
        t_ohmic = 0.0_s
        do m=1,max_carriers
            if ((c(m)%valley > 0) .and. (c(m)%op == 0)) then
                !print *, c(m)%cidx, c(m)%r
                x = c(m)%r(1)/m_dx(1)
                y = c(m)%r(2)/m_dx(1)
                z = c(m)%r(3)/m_dx(1)
                i = int(x+1)
                j = int(y+1)
                k = int(z+1)
                el_charge(i,j,k) = el_charge(i,j,k)+(1-(x-(i-1)))*(1-(y-(j-1)))*(1-(z-(k-1)))
                el_charge(i,j+1,k) = el_charge(i,j+1,k)+(1-(x-(i-1)))*(y-(j-1))*(1-(z-(k-1)))
                el_charge(i,j+1,k+1) = el_charge(i,j+1,k+1)+(1-(x-(i-1)))*(y-(j-1))*(z-(k-1))
                el_charge(i,j,k+1) = el_charge(i,j,k+1)+(1-(x-(i-1)))*(1-(y-(j-1)))*(z-(k-1))
                el_charge(i+1,j+1,k+1) = el_charge(i+1,j+1,k+1)+(x-(i-1))*(y-(j-1))*(z-(k-1))
                el_charge(i+1,j,k) = el_charge(i+1,j,k)+(x-(i-1))*(1-(y-(j-1)))*(1-(z-(k-1)))
                el_charge(i+1,j+1,k) = el_charge(i+1,j+1,k)+(x-(i-1))*(y-(j-1))*(1-(z-(k-1)))
                el_charge(i+1,j,k+1) = el_charge(i+1,j,k+1)+(x-(i-1))*(1-(y-(j-1)))*(z-(k-1))
            end if
        end do

        do x_idx=1,num_nodes(1)
            do y_idx=1,num_nodes(2)
                do z_idx=1,num_nodes(3)
                    el_charge(x_idx,y_idx,z_idx) = el_charge(x_idx,y_idx,z_idx)*p_charge
                    el_conc(x_idx,y_idx,z_idx) = el_charge(x_idx,y_idx,z_idx)/(m_lx(x_idx)*m_ly(y_idx)*m_lz(z_idx))

                    ! calculate factor for keeping ohmic contacts
                    if (step == 120) then
                        t_ohmic = 1! mod(step,100)
                    end if
                    call coord2idx(x_idx,y_idx,z_idx,orientation,coord1,coord2)
                    if (orientation > 0) then
                    if ((contact_type(orientation,coord1,coord2) == OHMIC) .and. (t_ohmic >= 0.5_s)) then
                        conc_tmp = sum(dev%g_el_conc(40:step,x_idx,y_idx,z_idx))/80
                        ohmic_f(orientation,coord1,coord2) = 1.0_s-(conc_tmp-donor(m_region(x_idx,y_idx,z_idx)))/conc_tmp
                        !print *, conc_tmp, ohmic_f(x_idx,y_idx)
                        if ((ohmic_f(orientation,coord1,coord2) <= 0.0_s) .or. (ohmic_f(orientation,coord1,coord2) > 1.0_s)) then
                            ohmic_f(orientation,coord1,coord2) = 1.0_s
                        end if
                        !print *, conc_tmp, ohmic_f(x_idx,y_idx)
                    end if
                    end if
                    end do
            end do
        end do
    end subroutine assign_charge_3
#endif

#if DIM == 2
    subroutine inject_carrier_2()
        integer                             :: i, x_idx, y_idx
        real(kind=s)                        :: r, t_left
        integer                             :: mat!, loc

        do x_idx=1,num_nodes(1)
        do y_idx=1,num_nodes(2)
            if (contact_type(x_idx,y_idx) .eq. OHMIC) then
            i = 0
            !print *, 'inject start', y_idx, donor(x_idx,y_idx), contact_conc(x_idx,y_idx)
            !print *, ohmic_f(x_idx,y_idx)
            do while ((1*donor(x_idx,y_idx)-contact_conc(x_idx,y_idx)) .gt. (donor(x_idx,y_idx)*0.001)) !ohmic_f(x_idx,y_idx) ! 0.89
                    i = i+1
                    if (c(i)%valley .eq. 0) then
                        mat = m_material(m_region(x_idx,y_idx))
                        c(i)%t_id = 1
                        call c(i)%init(mat)
                        c(i)%cidx = i
                        !c(i)%vel(1) = abs(c(i)%vel(1))
                        !c(i)%k(1) = abs(c(i)%k(1))

                        !if (x_idx == 1) then
                        !    call random_number(r)
                        !    c(i)%vel(1) = sqrt((-1.0_s*log(r)*2.0_s*kB*T_lattice)/(me(mat,1)))
                        !    c(i)%k(1) = c(i)%vel(1)*(me(mat,1)*(1.0_s+2.0_s* &
                        !            & nonp(mat,1)*c(i)%energy))/hb
                        !    c(i)%r(1) = 0.0_s+2e-5_s
                        !end if
                        !if (x_idx == num_nodes(1)) then
                        !    call random_number(r)
                        !    c(i)%vel(1) = -sqrt((-1.0_s*log(r)*2.0_s*kB*T_lattice)/(me(mat,1)))
                        !    c(i)%k(1) = c(i)%vel(1)*(me(mat,1)*(1.0_s+2.0_s* &
                        !            & nonp(mat,1)*c(i)%energy))/hb
                        !    c(i)%r(1) = dev_x-2e-5_s
                        !   !print *, 'init', c(i)%cidx, c(i)%r(1:2), c(i)%vel(1:2), c(i)%k(1:2)
                        !end if
                        if (y_idx == num_nodes(2)) then
                            call random_number(r)
                            c(i)%vel(2) = -sqrt((-1.0_s*log(r)*2.0_s*kB*T_lattice)/(me(mat,1)))
                            c(i)%k(2) = c(i)%vel(2)*(me(mat,1)*(1.0_s+2.0_s* &
                                    & nonp(mat,1)*c(i)%energy))/hb
                            c(i)%r(2) = dev_y-2e-5_s
                            !print *, 'init', c(i)%cidx, c(i)%r(1:2), c(i)%vel(1:2), c(i)%k(1:2)
                        end if
                        call random_number(r)
                        !c(i)%r(2) = node_coord(2,y_idx) !(r-0.5_s)*m_ly(y_idx)+node_coord(2,y_idx)
                        !if (y_idx == 1) then
                        !        c(i)%r(2) = node_coord(2,y_idx) ! r*m_ly(y_idx)+
                        !end if
                        !if (y_idx == num_nodes(2)) then
                        !        c(i)%r(2) = node_coord(2,y_idx)-m_ly(y_idx)/2.0_s
                        !        !print *, c(i)%cidx, c(i)%r(1:2), c(i)%vel(1:2)
                        !end if
                        c(i)%r(1) = node_coord(1,x_idx) !(r-0.5_s)*m_ly(y_idx)+node_coord(2,y_idx)
                        if (x_idx == 1) then
                                c(i)%r(1) = node_coord(1,x_idx) ! r*m_ly(y_idx)+
                        end if
                        if (x_idx == num_nodes(1)) then
                                c(i)%r(1) = node_coord(1,x_idx)-m_lx(x_idx)/2.0_s
                                !print *, c(i)%cidx, c(i)%r(1:2), c(i)%vel(1:2)
                        end if
                        call random_number(r)
                        t_left = r*t_step
                        if (t_left .gt. t_rest(x_idx,y_idx)) then
                            t_left = t_rest(x_idx,y_idx)
                        end if
                        c(i)%e_time = step*t_step-t_left
                        do while ((c(i)%synchronized < 0.1) .and. (c(i)%valley > 0))
                            call c(i)%drift(1)
                        end do
                        !print *, c(i)%cidx, c(i)%r(1:2), c(i)%vel(1:2)
                        c(i)%synchronized = 0
                        call keep_ohmic_contacts_2()
                        !if (c(i)%cidx == 300907) then
                        !    print *, c(i)%cidx, c(i)%r(1:2)
                        !end if
                        !print *, 'after injection', donor(x_idx,y_idx), contact_conc(x_idx,y_idx)
                    end if
            end do
                    !print *, 'inject end', y_idx, donor(x_idx,y_idx), contact_conc(x_idx,y_idx)
            end if
        end do
        end do
    end subroutine inject_carrier_2
#endif

#if DIM == 3
    subroutine inject_carrier_3()
        integer                             :: i, x_idx, y_idx, z_idx, region, orientation, coord1, coord2
        real(kind=s)                        :: r, t_left

        !do x_idx=1,num_nodes(1)
        x_idx = 1
        !!$omp parallel private(y_idx, z_idx) shared(c,i)
        do y_idx=1,num_nodes(2)
        !!$omp do
        do z_idx=1,num_nodes(3)
            call coord2idx(x_idx,y_idx,z_idx,orientation,coord1,coord2)
            if (contact_type(orientation,coord1,coord2) .eq. OHMIC) then
            i = 0
            region = m_region(x_idx,y_idx,z_idx)
            !print *, 'inject start', z_idx, contact_conc(orientation,coord1,coord2)
            !print *, ohmic_f(x_idx,y_idx)
            do while ((ohmic_f(orientation,coord1,coord2)*donor(region)-contact_conc(orientation,coord1,coord2)) &
                            & .gt. (donor(region)*0.001)) !ohmic_f(x_idx,y_idx) ! 0.89
                    i = i+1
                    if (c(i)%valley .eq. 0) then
                        call c(i)%init(region)
                        !c(i)%t_id = 1
                        c(i)%t_id = omp_get_thread_num()+1
                        c(i)%cidx = i
                        !c(i)%vel(1) = abs(c(i)%vel(1))
                        !c(i)%k(1) = abs(c(i)%k(1))
                        if (x_idx == 1) then
                            !call random_number(r)
                            r = rng_uniform(rng(c(i)%t_id))
                            c(i)%vel(1) = sqrt((-1.0_s*log(r)*2.0_s*kB*T_lattice)/(me(region,1)))
                            c(i)%k(1) = c(i)%vel(1)*(me(region,1)*(1.0_s+2.0_s* &
                                    & nonp(region,1)*c(i)%energy))/hb
                            c(i)%r(1) = 0.0_s+2e-5_s
                        end if
                        if (x_idx == num_nodes(1)) then
                            !call random_number(r)
                            r = rng_uniform(rng(c(i)%t_id))
                            c(i)%vel(1) = -sqrt((-1.0_s*log(r)*2.0_s*kB*T_lattice)/(me(region,1)))
                            c(i)%k(1) = c(i)%vel(1)*(me(region,1)*(1.0_s+2.0_s* &
                                    & nonp(region,1)*c(i)%energy))/hb
                            c(i)%r(1) = dev_x-2e-5_s
                            !print *, 'init', c(i)%cidx, c(i)%r(1:2), c(i)%vel(1:2), c(i)%k(1:2)
                        end if
                        !call random_number(r)
                        c(i)%r(2) = node_coord(2,y_idx) !(r-0.5_s)*m_ly(y_idx)+node_coord(2,y_idx)
                        c(i)%r(3) = node_coord(3,z_idx)
                        if (y_idx == 1) then
                                c(i)%r(2) = node_coord(2,y_idx)+m_lx(1)/2.0_s ! r*m_ly(y_idx)+
                        end if
                        if (y_idx == num_nodes(2)) then
                                c(i)%r(2) = node_coord(2,y_idx)-m_lx(1)/2.0_s !-r*m_ly(y_idx)+
                                !print *, c(i)%cidx, c(i)%r(1:2), c(i)%vel(1:2)
                        end if
                        if (z_idx == 1) then
                                c(i)%r(3) = node_coord(3,z_idx)+m_lx(1)/2.0_s ! r*m_ly(y_idx)+
                        end if
                        if (z_idx == num_nodes(3)) then
                                c(i)%r(3) = node_coord(3,z_idx)-m_lx(1)/2.0_s !-r*m_ly(y_idx)+
                                !print *, c(i)%cidx, c(i)%r(1:2), c(i)%vel(1:2)
                        end if
                        !call random_number(r)
                        r = rng_uniform(rng(c(i)%t_id))
                        t_left = r*t_step
                        if (t_left .gt. t_rest(orientation,coord1,coord2)) then
                            t_left = t_rest(orientation,coord1,coord2)
                        end if
                        c(i)%e_time = step*t_step-t_left
                        do while ((c(i)%synchronized < 0.1) .and. (c(i)%valley > 0) .and. (c(i)%op == 0))
                            call c(i)%drift(1)
                        end do
                        c(i)%synchronized = 0
                        call keep_ohmic_contacts_3()
                        !print *, c(i)%cidx, c(i)%r, c(i)%e_time, c(i)%synchronized
                        !if (c(i)%cidx == 300907) then
                        !    print *, c(i)%cidx, c(i)%r(1:2)
                        !end if
                        !print *, 'after injection',  z_idx, contact_conc(orientation,coord1,coord2)
                    end if
            end do
                    !print *, 'inject end', z_idx,  z_idx, contact_conc(orientation,coord1,coord2)
            end if
        end do
        !!$omp end do
        end do
        !!$omp end parallel
        !print *, 'inject end', contact_conc(LEFT,1:num_nodes(2),1:num_nodes(3))
    end subroutine inject_carrier_3
#endif

end module poisson
