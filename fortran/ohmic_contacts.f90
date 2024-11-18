module ohmic_contacts
    use config
    use materialdef
    use device
    use time_2border

    implicit none
#if DIM == 1
    real(kind=s),dimension(2)               :: time_contacts, contact_conc, t_rest
#endif
#if DIM == 2
    real(kind=s),dimension(:,:),allocatable :: time_contacts, contact_conc, t_rest
    real(kind=s),dimension(:,:,:),allocatable :: time_contacts_tmp
#endif
#if DIM == 3
    real(kind=s),dimension(:,:,:),allocatable :: time_contacts, contact_conc, t_rest
    real(kind=s),dimension(:,:,:,:),allocatable :: time_contacts_tmp
#endif

    contains
#if DIM == 1
    subroutine init_contacts
        time_contacts = 0.0_s
        call keep_ohmic_contacts()
        !dev%particle_flow_out = 0
    end subroutine init_contacts

    subroutine time_in_contact(cidx, r_new, dr, t_delta)
        integer                             :: i, cidx
        real(kind=s)                        :: r_new, dr, t_delta, last_cell

        cidx = cidx
        last_cell = dev_x-node_dist(1,sum(num_nodes))-1e-5_s
        if (contact_type(1) .eq. OHMIC) then
            if ((r_new .lt. (m_lx(1)+1e-5_s)) .and. ((r_new-dr) .lt. (m_lx(1)+1e-5_s))) then
                time_contacts(1) = time_contacts(1)+t_delta
                !print *, cidx, r_new, r_new-dr, t_delta
            end if
        end if
!        if (contact_type(1) .eq. OHMIC) then
!            if ((r_new .lt. (m_lx(1)+1e-5_s))) then
!                !time_contacts(1) = time_contacts(1)+t_delta
!                !print *, 'only in first cell:', cidx, r_new, r_new-dr, t_delta
!            end if
!        end if
        if (contact_type(2) .eq. OHMIC) then
            if ((r_new .gt. last_cell) .and. &
             & ((r_new-dr) .gt. last_cell)) then
                time_contacts(2) = time_contacts(2)+t_delta
            end if
        end if
!        if (contact_type(2) .eq. OHMIC) then
!            if (r_new .gt. last_cell) then
!                !time_contacts(1) = time_contacts(1)+t_delta
!                !print *, 'only in last cell:', cidx, r_new, r_new-dr, t_delta
!            end if
!        end if
        !print *, time_contacts(1)
    end subroutine time_in_contact

    subroutine keep_ohmic_contacts()
        integer             :: max_nodes

        max_nodes = sum(num_nodes)
        contact_conc(1) = (time_contacts(1)/t_step)*p_charge/(m_lx(1))
        contact_conc(2) = (time_contacts(2)/t_step)*p_charge/m_lx(max_nodes)
        t_rest(1) = (donor(1)-contact_conc(1))*t_step*(m_lx(1))/p_charge
        t_rest(2) = (donor(max_nodes)-contact_conc(2))*t_step*m_lx(max_nodes)/p_charge

        !print *, contact_conc(1), contact_conc(2), time_contacts
    end subroutine keep_ohmic_contacts
#endif

#if DIM == 2
!    subroutine init_contacts_2()
!        integer                     :: n_out
!
!        allocate(time_contacts_tmp(NUM_THREADS,num_nodes(1),num_nodes(2)))
!        allocate(time_contacts(num_nodes(1),num_nodes(2)))
!        allocate(contact_conc(num_nodes(1),num_nodes(2)))
!        allocate(t_rest(num_nodes(1),num_nodes(2)))
!        time_contacts_tmp = 0.0_s
!        time_contacts = 0.0_s
!        contact_conc = 0.0_s
!        t_rest = 0.0_s
!
!        n_out = 2
!        allocate(particle_in(n_out))
!        particle_in = 0
!        allocate(particle_out(NUM_THREADS,n_out))
!        particle_out = 0
!        allocate(dev%particle_flow_in(int(n_step/save_t),n_out))
!        dev%particle_flow_in = 0
!        allocate(dev%particle_flow_out(int(n_step/save_t),n_out))
!        dev%particle_flow_out = 0
!    end subroutine

    subroutine time_in_contact_2(cidx, t_id, r_new, dr, t_delta)
        integer                             :: cidx, t_id, x_idx, y_idx
        integer,dimension(3)                :: idx_old_1, idx_old_2, idx_new_1, idx_new_2
        real(kind=s),dimension(3)           :: r_new, dr
        real(kind=s)                        :: t_delta

        call get_r_idx_2((r_new-dr)-1e-5_s,idx_old_1)
        call get_r_idx_2((r_new-dr)+1e-5_s,idx_old_2)
        call get_r_idx_2(r_new-1e-5_s,idx_new_1)
        call get_r_idx_2(r_new+1e-5_s,idx_new_2)

        x_idx = 0
        y_idx = 0
        if (idx_old_1(1) == idx_new_1(1)) then
            x_idx = idx_old_1(1)
        else if (idx_old_1(1) == idx_new_2(1)) then
            x_idx = idx_old_1(1)
        else if (idx_old_2(1) == idx_new_1(1)) then
            x_idx = idx_old_2(1)
        else if (idx_old_2(1) == idx_new_2(1)) then
            x_idx = idx_old_2(1)
        end if

        if (idx_old_1(2) == idx_new_1(2)) then
            y_idx = idx_old_1(2)
        else if (idx_old_1(2) == idx_new_2(2)) then
            y_idx = idx_old_1(2)
        else if (idx_old_2(2) == idx_new_1(2)) then
            y_idx = idx_old_2(2)
        else if (idx_old_2(2) == idx_new_2(2)) then
            y_idx = idx_old_2(2)
        end if

        if (x_idx == 0) then
            print *, cidx, t_id, r_new, r_new-dr
            print *, idx_new_1, idx_new_2, idx_old_1, idx_old_2
        end if

        if (contact_type(x_idx,y_idx) == OHMIC) then
            time_contacts_tmp(t_id,x_idx,y_idx) = time_contacts_tmp(t_id,x_idx,y_idx)+t_delta
            !print *, cidx, r_new(1:2),r_new(1:2)-dr(1:2), x_idx, y_idx, t_delta
        end if
    end subroutine time_in_contact_2

    subroutine keep_ohmic_contacts_2()
        integer             :: x_idx, y_idx

        do x_idx=1,num_nodes(1)
            do y_idx=1,num_nodes(2)
                if (contact_type(x_idx,y_idx) == OHMIC) then
                    time_contacts(x_idx,y_idx) = sum(time_contacts_tmp(1:NUM_THREADS,x_idx,y_idx))
                    contact_conc(x_idx,y_idx) = ((time_contacts(x_idx,y_idx)/t_step)*p_charge*2.0_s)/(m_lx(x_idx)*m_ly(y_idx))
                    t_rest(x_idx,y_idx) = (donor(x_idx,y_idx)-contact_conc(x_idx,y_idx))*t_step*(m_lx(x_idx)*m_ly(y_idx))/ &
                                    & (2.0_s*p_charge)
                    !print *, x_idx, y_idx, time_contacts(x_idx,y_idx), m_lx(x_idx)*m_ly(y_idx)
                end if
            end do
        end do
        !print *, contact_conc(1,:)
        !print *, t_rest(1,:)
    end subroutine keep_ohmic_contacts_2
#endif

#if DIM == 3
    subroutine init_contacts_3()
        allocate(time_contacts_tmp(NUM_THREADS,6,maxval(num_nodes),maxval(num_nodes)))
        allocate(time_contacts(6,maxval(num_nodes),maxval(num_nodes)))
        allocate(contact_conc(6,maxval(num_nodes),maxval(num_nodes)))
        allocate(t_rest(6,maxval(num_nodes),maxval(num_nodes)))
        time_contacts_tmp = 0.0_s
        time_contacts = 0.0_s
        contact_conc = 0.0_s
        t_rest = 0.0_s
    end subroutine

    subroutine time_in_contact_3(cidx, t_id, r_new, dr, t_delta)
        integer                             :: cidx, t_id, x_idx, y_idx, z_idx, orientation, coord1, coord2
        integer,dimension(3)                :: idx_old_1, idx_old_2, idx_new_1, idx_new_2
        real(kind=s),dimension(3)           :: r_new, dr
        real(kind=s)                        :: t_delta

        call get_r_idx_3((r_new-dr)-1e-5_s,idx_old_1)
        call get_r_idx_3((r_new-dr)+1e-5_s,idx_old_2)
        call get_r_idx_3(r_new-1e-5_s,idx_new_1)
        call get_r_idx_3(r_new+1e-5_s,idx_new_2)

        x_idx = 0
        y_idx = 0
        z_idx = 0
        if (idx_old_1(1) == idx_new_1(1)) then
            x_idx = idx_old_1(1)
        else if (idx_old_1(1) == idx_new_2(1)) then
            x_idx = idx_old_1(1)
        else if (idx_old_2(1) == idx_new_1(1)) then
            x_idx = idx_old_2(1)
        else if (idx_old_2(1) == idx_new_2(1)) then
            x_idx = idx_old_2(1)
        end if
        !print *, 'xidx', cidx, idx_old_1(1), idx_old_2(1), idx_new_1(1), idx_new_2(1), t_delta
        if (idx_old_1(2) == idx_new_1(2)) then
            y_idx = idx_old_1(2)
        else if (idx_old_1(2) == idx_new_2(2)) then
            y_idx = idx_old_1(2)
        else if (idx_old_2(2) == idx_new_1(2)) then
            y_idx = idx_old_2(2)
        else if (idx_old_2(2) == idx_new_2(2)) then
            y_idx = idx_old_2(2)
        end if

        if (idx_old_1(3) == idx_new_1(3)) then
            z_idx = idx_old_1(3)
        else if (idx_old_1(3) == idx_new_2(3)) then
            z_idx = idx_old_1(3)
        else if (idx_old_2(3) == idx_new_1(3)) then
            z_idx = idx_old_2(3)
        else if (idx_old_2(3) == idx_new_2(3)) then
            z_idx = idx_old_2(3)
        end if

        !print *, r_new(1:3), x_idx, y_idx, z_idx
        !print *, '              ', r_new-dr
        !call get_r_idx_3(r_new,idx_old_1)
        !print *, idx_old_1
        call coord2idx(x_idx,y_idx,z_idx,orientation,coord1,coord2)
        if (orientation > 0) then
        if (contact_type(orientation,coord1,coord2) == OHMIC) then
            time_contacts_tmp(t_id,orientation,coord1,coord2) = time_contacts_tmp(t_id,orientation,coord1,coord2)+t_delta
            !print *, cidx, r_new(1:2),r_new(1:2)-dr(1:2), x_idx, y_idx, t_delta
            !print *, 'ohmic', cidx, orientation,coord1,coord2,y_idx,z_idx, r_new(1:3), sum(time_contacts_tmp)
        end if
        end if
    end subroutine time_in_contact_3

    subroutine keep_ohmic_contacts_3()
        integer             :: x_idx, y_idx, z_idx, orientation, coord1, coord2

        do x_idx=1,num_nodes(1)
            do y_idx=1,num_nodes(2)
                do z_idx=1,num_nodes(3)
                    call coord2idx(x_idx,y_idx,z_idx,orientation,coord1,coord2)
                    if (orientation > 0) then
                    if (contact_type(orientation,coord1,coord2) == OHMIC) then
                time_contacts(orientation,coord1,coord2) = sum(time_contacts_tmp(1:NUM_THREADS,orientation,coord1,coord2))
                contact_conc(orientation,coord1,coord2) = ((time_contacts(orientation,coord1,coord2)/t_step)* &
                                                    & p_charge*1.0_s)/(m_lx(x_idx)*m_ly(y_idx)*m_lz(z_idx))
                t_rest(orientation,coord1,coord2) = (donor(m_region(x_idx,y_idx,z_idx))-contact_conc(orientation,coord1,coord2))* &
                                        & t_step*(m_lx(x_idx)*m_ly(y_idx)*m_lz(z_idx))/(1.0_s*p_charge)
                        !print *, x_idx, y_idx, time_contacts(x_idx,y_idx), m_lx(x_idx)*m_ly(y_idx)
                    end if
                    end if
                    end do
            end do
        end do
        !print *, contact_conc(1,:)
        !print *, t_rest(1,:)
    end subroutine keep_ohmic_contacts_3
#endif

end module ohmic_contacts
