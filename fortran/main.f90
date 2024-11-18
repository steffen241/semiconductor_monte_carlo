program main
    use omp_lib
    use physconst
    use materialdef
    use scattering_rates
    use mc_core
    use config
    !use gnufor2
    use show_output
    use poisson
    use time_2border
    use ohmic_contacts
    use device
    use integrals
    use reservoir
    use tunneling
    use save_hdf5
    use pauli

    implicit none
    integer                     :: counti, countf, count_rate
    integer                     :: counti2, countf2, count_rate2
    integer :: i, i_tun, j
    integer :: c_idx
    !integer                               :: i, j
    !real(kind=s),dimension()           :: tmp
    !real(kind=s),dimension(2,300000)        :: test, test2
    !real(kind=s)                          :: j2

    integer                         :: proc_id, nproc, ierror
    !** Init MPI environment **!
    call MPI_Init(ierror)

    call MPI_Comm_rank(MPI_COMM_WORLD, proc_id , ierror)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    !print *, proc_id, nproc

    ! initialize some stuff
    call init_const
    call init_config

!    call gauss_input(10.0_s,j2)
!    print *, j2

!#if TESTCC == 0

#if DIM < 4
        if (proc_id == 0) then
            print *, 'Start scat rates'
    call system_clock(counti,count_rate)
            call init_parameter
        print *, nonp(:,3)
    call system_clock(countf)
    print *, 'Time elapsed:', real(countf-counti)/real(count_rate)
    p%mat = INGAAS
    p%valley = 1
    !call test_fd()

            !print *, D(GAAS,1,1,:)
            !print *, D(GAAS,1,2,:)
            !print *, D(GAAS,1,3,:)
            print *, 'Start scat rates'
    call system_clock(counti,count_rate)
        call calc_scatrates
    call system_clock(countf)
    print *, 'Time elapsed:', real(countf-counti)/real(count_rate)
             call save_scatrates_hdf5
        call init_mc()
            !call save_materialdef
        end if
#endif

print *, eps_hf(INGAAS), eps_hf(INALAS)
    ! Bulk and 1D mode don't need MPI, only for solving poisson equation
    if (proc_id == 0) then
    !
    ! Bulk mode
    !
#if DIM == 0
    print *, 'Starting bulk simulation at lattice temperature: ', T_lattice, 'for material', m_material(1)
    !!$omp parallel private(i) shared(c)
    !!$omp do
    do i=1,num_carriers
        c(i)%t_id = omp_get_thread_num()+1
        call init(c(i),m_material(1))
        c(i)%cidx = i
    end do
    !!$omp end do
    !!$omp end parallel
    call system_clock(counti,count_rate)
    do j=1,1 !"25 !42 !49,49!1,49 !1,k_max
        !e_field = (/k_array(j),0.0_s,0.0_s/)
        do i=1,num_carriers
            !c(i)%loc_e_field = (/0.0_s,0.0_s,-k_array(j)/)
            c(i)%loc_e_field = (/0.0_s,0.0_s,0.0_s/)
            !c(i)%loc_e_field = (/k_array(j),k_array(j),k_array(j)/)/sqrt(3.0_s) ! Field in [111]-direction
        end do
        print *, 'Restarting simulation with E-Field:', k_array(j)*1e9, 'V/m', c(1)%loc_e_field
        c(:)%e_time = 0
        do step=1,n_step
            if (step > 20) then
                do i=1,num_carriers
                    c(i)%loc_e_field = (/0.0_s,0.0_s,-1e-2_s/)
                end do
            end if
            print *, 'Perform time step:', step*t_step, 'ps'
                    call system_clock(counti,count_rate)
            call emc()
                    call system_clock(countf)
                    print *, 'Time elapsed:', real(countf-counti)/real(count_rate)
            call calc_mean_val(int(j,4),int(step,4),c(1)%mat)
            call calc_valley_occ(int(j,4),int(step,4))
            call calc_pauli()
            call calc_fd()
            call save_pauli(int(j,4),int(step,4))
        end do
    end do
    call save_bulk_hdf5()
#endif

    !
    ! Device simulation
    !
#if DIM == 1
    print *, 'Starting device simulation in 1D'
    print *, me(INGAAS,1)/m0, me(INALAS,1)/m0
    print *, valley_offs(INGAAS,:), valley_offs(INALAS,:)
    call init_mesh_1
    call init_poisson_1
    call init_contacts
    ! initiliaze electron ensemble + initial solution of poisson equation
    !do i=1,num_carriers
    !    do while (c(i)%MAT == AIR)
    !        call c(i)%init_1()
    !    end do
    !    c(i)%cidx = i
    !end do
    !$omp parallel private(i) shared(c)
    !$omp do
    do i=1,num_carriers
        c(i)%t_id = omp_get_thread_num()+1
        do while (c(i)%MAT == AIR)
            call c(i)%init_1()
        end do
        c(i)%cidx = i
    end do
    !$omp end do
    !$omp end parallel
    call init_reservoir_1()
    call assign_charge_1()
    step = 1
    call solve_poisson_1()

    ! Do simulation, see config module for variables
    print *, 'Simulation at lattice temperature:', T_lattice
    do step=1,n_step
        print *, 'Doing time step:', step*t_step, 'ps'
        !call init_contacts()
        call emc()
        call res_drift_1()
        call assign_charge_1()
        call keep_ohmic_contacts()
        !call inject_carrier()
                call count_carriers()
        call assign_charge_1()
        call solve_poisson_1()
        call media()
    end do
    print *, 'Simulation end'
    call save_device()
    !end if
#endif
    end if




        !eps_static(AIR)      = eps_static(INALAS) !12.9_s*eps_vac !12.90_s*eps_vac
        !eps_hf(AIR)          = eps_hf(INALAS)
        print *, eps_static(AIR), eps_static(INALAS), eps_static(INGAAS), eps_vac
#if DIM == 2
        call init_mesh_2()
    if (proc_id == 0) then
        print *, 'Starting device simulation in 2D'
        print *, 'Simulation at lattice temperature:', T_lattice
        !do i=1,88
        !    print *, i, m_region(10,i)
        !end do
        print *, 'Mesh initialized'
        !call init_contacts_2()
        !print *, 'Contacts initialized'
                       call system_clock(counti,count_rate)
        if (load_carriers == 0) then
        !$omp parallel private(i) shared(c)
        !$omp do
        do i=1,(num_carriers)!-25000.0_s
            c(i)%t_id = omp_get_thread_num()+1
            do while (c(i)%MAT == AIR)
            call c(i)%init_2()
            end do
            c(i)%cidx = i
        end do
        !$omp end do
        !$omp end parallel
        else
            call load_carriers_2()
        end if
        print *, 'Carriers initialized'
                        call system_clock(countf)
                        print *, char(9), char(9), 'Init time elapsed:', real(countf-counti)/real(count_rate)
    end if


    call init_poisson_hips_2()
    if (proc_id == 0) then
        step = 1
    print *, 'Poisson solver initialized'
        call init_reservoir_2()
        call assign_charge_2()
        !print *, el_conc(43,num_nodes(2)-4:num_nodes(2))
        print *, 'Calculating tunneling probablity'
        call init_tunneling()
        call sav%save_init_2()
            call init_pauli()
            call test_fd()
            print *, 'Fermi-Dirac statistic initialized'
        !call tetest()
    end if
    call solve_poisson_hips_2()
    do step=1,n_step
        call system_clock(counti2,count_rate2)
        if (proc_id == 0) then

           !print *, donor(30,:)
                !call system_clock(counti,count_rate)
            print *, 'Doing time step:', step*t_step, 'ps'
        call system_clock(counti,count_rate)
            call emc()
        call system_clock(countf)
        print *, 'emc', real(countf-counti)/real(count_rate)
!            call keep_ohmic_contacts_2()
        call system_clock(counti,count_rate)
                call res_drift_2()
        call system_clock(countf)
        print *, 'res', real(countf-counti)/real(count_rate)
                call system_clock(counti,count_rate)
            call assign_charge_2()
                call system_clock(countf)
            print *, 'assign charge', real(countf-counti)/real(count_rate)
                call system_clock(counti,count_rate)
        end if
            call solve_poisson_hips_2()
        if (proc_id == 0) then
                call system_clock(countf)
            print *, 'poisson', real(countf-counti)/real(count_rate)
                !call system_clock(counti,count_rate)
            !call media_2()

            call save_tun_pot()
        if (pep_hdf5 == 1) then
            call save_pauli()
        end if
        if (n_tunnel > 0) then
        if (mod(step,tun(1)%steps_compute_tunnel) == 0) then
            do i_tun=1,n_tunnel
            print *, 'Computing tunneling rates: tun', i_tun
            call calc_tunnelrate(tun(i_tun),G)
            call calc_tunnelrate(tun(i_tun),L)
            call calc_tunnelrate(tun(i_tun),X)
            end do
        end if
        end if
        if ((mod(step,p%steps_compute_pep) == 0) .and. (pep_hdf5 == 1)) then
                call system_clock(counti,count_rate)
            print *, 'Solving self-consistently Fermi-Dirac statistics'
            call calc_fd
            p%step = 1
            sav%write_pep = 1
               call system_clock(countf)
                print *, char(9), char(9), 'Program time elapsed:', real(countf-counti)/real(count_rate)

        end if
            if (mod(step,save_t) == 0) then
                call count_carriers()
                call sav%calc_media_2()
                call sav%write_media_2()
            end if
!            if (mod(step,save_t*10) == 0) then
!                        open(12, file="test.txt", status="unknown", position="append", action="write")
!                do c_idx=1,max_carriers
!                    if ((c(c_idx)%valley == 1) .and. (c(c_idx)%op == 0) .and. (step*t_step >= 2.0_s) &
!                        & .and. (c(c_idx)%r(2) > 100.0_s) .and. (c(c_idx)%r(2) < 200.0_s)) then
!                        write (12,*) c(c_idx)%energy
!                    end if
!                end do
!            end if
        end if
        call system_clock(countf2)
        print *, 'Total', real(countf2-counti2)/real(count_rate2)
    end do
    if (proc_id == 0) then
        print *, 'Simulation end'
        ! Save particles states if not loaded at start
        if (load_carriers == 0) then
        print *, 'Save carrier'
        open(12, file="/scratch/schumans/13_carriers.txt", status="unknown", position="append", action="write")
        do i=1,max_carriers
            if (c(i)%valley > 0) then
            write (12,*) c(i)%r(1:2), c(i)%k, c(i)%energy, c(i)%valley, c(i)%v_subidx, c(i)%mat, c(i)%op, c(i)%res_id, &
                    & c(i)%res_scat_idx
            end if
        end do
        end if
    end if
#endif



#if DIM == 3
    call init_mesh_3()
    if (proc_id == 0) then
        print *, 'Starting device simulation in 3D'
        print *, 'Simulation at lattice temperature:', T_lattice
        step = 1
        call init_contacts_3()
        !$omp parallel private(i) shared(c)
        !$omp do
        !schedule(static,100)
        do i=1,num_carriers
            c(i)%t_id = omp_get_thread_num()+1
            call c(i)%init_3()
            c(i)%cidx = i
        end do
        !$omp end do
        !$omp end parallel
    end if
    call init_poisson_hips_3()
    if (proc_id == 0) then
        call assign_charge_3()
    end if
    call solve_poisson_hips_3()
!    if (proc_id == 0) then
!        call assign_force_3()
!        print *, sqrt(dot_product(c(2)%loc_e_field,c(2)%loc_e_field))
!    end if
    do step=1,n_step
        if (proc_id == 0) then
            print *, 'Doing time step:', step*t_step, 'ps'
                call system_clock(counti,count_rate)
            call emc()
            call keep_ohmic_contacts_3()
            call inject_carrier_3()
            time_contacts_tmp = 0.0_s
            call assign_charge_3()
            call count_carriers()
        end if
        call solve_poisson_hips_3()
        if (proc_id == 0) then
                call system_clock(countf)
                print *, char(9), char(9), 'Program time elapsed:', real(countf-counti)/real(count_rate)
        end if
    end do
    if (proc_id == 0) then
        print *, 'Simulation end'
        call save_device()
    end if
#endif
    call MPI_Finalize(ierror)
end program main
