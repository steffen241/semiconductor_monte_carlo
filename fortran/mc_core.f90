module mc_core
    use omp_lib
    use physconst
    use materialdef
    use scattering_rates
    use config
    use time_2border
    use ohmic_contacts
    use device
    use rng_mod
    use tunneling
    use pauli
    implicit none

   ! variables for carrier
    type                            :: electron
        integer                     :: t_id, visit                                          ! id if the thread, used for the random number generator
        real(kind=s),dimension(3)   :: r, r_old, r_old2 = 0.0_s                                       ! real space: x,y,z
        real(kind=s),dimension(3)   :: k, old_k, old_k2  = 0.0_s                                    ! reciprocal space: kx,ky,kz
        integer                     :: valley = 0                                           ! valley index
        integer                     :: v_subidx = 0                                         !
        real(kind=s)                :: energy = 0.0_s                                       ! electron energy
        integer                     :: mat    = AIR                                        ! Material index
        real(kind=s),dimension(3)   :: vel                                                  ! electron velocity vector
        real(kind=s)                :: e_time = 0                                           ! total propagation time of electron
        integer                     :: synchronized = 0                                     ! in asynchronous mode: 0, reached synchronous time step: 1
        integer                     :: e_idx = 0
        integer                     :: op = 0                                               ! operation: device = 0, reservoir = 1
        integer                     :: res_id = 0
        integer,dimension(2)        :: res_scat_idx = 0
        integer                     :: last_action
        real(kind=s),dimension(4)        :: tdelta

        integer,dimension(2)        :: last_scat_idx
        real(kind=s),dimension(2)   :: last_e
        real(kind=s),dimension(2,3) :: last_k
        real(kind=s)                :: ins_coord = 0

        ! Backup variables, in case the new state is rejected due to Pauli exclusion principle
        real(kind=s),dimension(3)   :: pep_r, pep_k, pep_vel
        real(kind=s)                :: pep_energy
        integer                     :: pep_valley, pep_vsubidx


#if DIM == 1
        integer                     :: r_idx = 0
#endif
#if DIM == 2
        integer,dimension(2)        :: r_idx = 0
#endif
#if DIM == 3
        integer,dimension(3)        :: r_idx = 0
#endif
        integer                     :: cidx = 0
        real(kind=s),dimension(3)   :: loc_e_field = 0                                      ! E-Field at carrier position
        !real                       :: energy_time = 0, vel_time = 0
    contains
        procedure                   :: init
#if DIM == 0
        procedure                   :: check_pep_0
#endif
#if DIM == 1
        procedure                   :: init_1
        procedure                   :: assign_force_1
#endif
#if DIM == 2
        procedure                   :: init_2
        procedure                   :: assign_force_2
        procedure                   :: check_pep
#endif
#if DIM == 3
        procedure                   :: init_3
        procedure                   :: assign_force_3
#endif
        procedure                   :: init_X
        procedure                   :: drift
        procedure                   :: res_drift
        procedure                   :: time_2e
        procedure                   :: compute_vel
        !procedure                   :: scatter_carrier
    end type electron

    type(electron),dimension(:),allocatable   :: c

    contains

    subroutine init_mc()
        allocate(c(max_carriers))
    end subroutine init_mc

    subroutine load_carriers_2()
        integer     :: lun = 10, i
        open(unit=lun, file="/scratch/schumans/10_carriers.txt", form="FORMATTED")
        do i=1,343670
            read(lun,*) c(i)%r(1:2), c(i)%k, c(i)%energy, c(i)%valley, c(i)%v_subidx, c(i)%mat, c(i)%op, c(i)%res_id, &
                    & c(i)%res_scat_idx
        end do
        print *, c(1)%r, c(1)%valley, c(1)%energy
    end subroutine load_carriers_2

    subroutine emc()
        integer(kind=8)             :: i
        !print *, 'schottky', SCHOTTKY

        !$omp parallel private(i) shared(c)
        !$omp do schedule(dynamic,80)
        !schedule(static,100)
        !schedule(dynamic,1000)
        !count_carriers = 0
            do i=1,max_carriers
                c(i)%t_id = omp_get_thread_num()+1
                !print *, c(i)%cidx, c(i)%t_id
                !!print *, 'start', c(i)%cidx, c(i)%r, c(i)%vel, c(i)%synchronized, c(i)%e_time
                do while ((c(i)%synchronized < 0.1) .and. (c(i)%valley > 0) .and. (c(i)%op == 0))
                    call c(i)%drift(0)
                end do
                c(i)%synchronized = 0
            end do
        !$omp end do
        !$omp end parallel
#if DIM == 1
        call keep_ohmic_contacts()
#endif
    end subroutine emc

    subroutine drift(this,mode)
        class(electron)             :: this
        integer                     :: idx, mode, ff_reject                    ! mode: 0 mc drift, mode: 1 injection (abort drift at border crossing
        integer,dimension(1)        :: choose_action, n_idx
        real(kind=s)                :: t_sync, t_scatter, t_delta, r
        real(kind=s)                :: l_me, l_nonp, e_old
        real(kind=s)                :: t_energy, t_real, el_fermi, el_temp, f_occ, E_tmp
        real(kind=s),dimension(3)   :: dk, dr, e_field, k_tmp, dr2
        integer                     :: addx1, addx2, subx1, subx2, addy1, addy2, suby1, suby2
        integer                     :: cross_hetero

#if DIM > 0
        integer,dimension(1)        :: nr_idx
#endif
#if DIM == 2
        integer,dimension(9)        :: n_mat
        integer                     :: did_tunnel
!        integer                     :: idx_verbose = 100020
#endif

#if DIM == 1
        if (mode == 0) then
            call get_r_idx(this%r,this%r_idx)
            call this%assign_force_1
        end if
#endif
#if DIM == 2
        call get_r_idx_2(this%r,this%r_idx)
        call get_e_idx(this)
        !if (mode == 0) then
            call this%assign_force_2
        !else
        !    this%loc_e_field(1) = el_field(1,this%r_idx(1),num_nodes(2)) !el_field(1,c(i)%ins_coord,num_nodes(2))
        !    this%loc_e_field(2) = (el_field(2,this%r_idx(1),num_nodes(2)-1)+el_field(2,this%r_idx(1),num_nodes(2)-2))/2.0_s
        !end if
#endif
#if DIM == 3
        call get_r_idx_3(this%r,this%r_idx)
        call this%assign_force_3
        !print *, this%loc_e_field
#endif

                if (this%valley < 1) then
                    print *, 3, this%energy, this%r(1:2), this%mat, this%valley
                end if

#if DIM > 0
        ! Calculate velocity
        call this%compute_vel()
#endif
       !print *, 'start', this%cidx
        ff_reject = 1
        cross_hetero = 0
        e_field = this%loc_e_field
        l_me = me(this%mat,this%valley)
        l_nonp = nonp(this%mat,this%valley)
        !e_old = this%energy
        !this%r_old2 = this%r_old
        !this%r_old = this%r
        !this%old_k2 = this%old_k
        !this%old_k = this%k

!        if (this%r_idx(1) < 30) then
            !print *, this%r, this%mat, this%energy, this%op
              !      this%valley = 0
!        end if

      !  if ((this%mat == INGAAS) .and. (this%r(2) > 59.49999999_s)) then !(this%r(2) > 39.5_s) .and. (this%r(2) < 59.5_s)) then !
      !      print *, this%k(1:2), this%old_k(1:2), this%old_k2(1:2)
       !     print *, 'you dont belong here', this%cidx, this%mat, this%valley, this%r(1:2), this%r_old(1:2), &
       !         & this%r_old2(1:2), this%vel(1:2), this%visit
       !     print *, this%last_action, this%energy, this%tdelta
       ! end if
!        if (this%energy .lt. 1e-5_s) then
!            this%energy = 1e-5_s
!        end if
        idx = E2idx(this%energy)

        !if (this%r(1) <= 2.5_s) then
        !   print *, 'start', this%cidx, this%r(1:2), time_contacts_tmp(1:NUM_THREADS,1,10:20)
        !end if
    !    if (this%cidx == 710) then
     !       print *, 'start', this%cidx, this%r(1:2), this%energy, this%valley
    !    end if

        if (idx <= 0) then
!          print *, 'e smaller 0',idx, this%energy, this%cidx, this%r_idx, mode, this%last_e, this%r(1:2), this%r_old(1:2), this%mat
            !this%energy = 1.0e-5_s
            !idx = E2idx(this%energy)
        else if (idx >= 300000) then
    !        print *, idx, this%energy, this%cidx, this%r
    !        print *, '----------------------------'
    !        print *, 'energy larger than 4 eV', mode, this%op, this%cidx, this%energy, this%valley, this%loc_e_field(1:2)
            this%energy = 1.99999
            idx = E2idx(this%energy)
    !        print *, '----------------------------'
        end if
        !print *, 'carrier start', this%r(1:2), this%energy

        ! Compute the time until next synchronous ensemble step
        t_sync = step*t_step-this%e_time
        !print *, t_sync, step, t_step, this%e_time
        !t_sync = t_step-this%e_time
        ! Calculate the time until a real space boundary is reached
#if DIM == 1
        call time2r_1(t_real,nr_idx,this%r,this%vel,this%r_idx,this%cidx)
#elif DIM == 2
        call time2r_2(t_real,nr_idx,this%r,this%r_idx,this%vel)
#elif DIM == 3
        call time2r_3(t_real,nr_idx,this%r,this%r_idx,this%vel)
#else
        t_real = 100000.0_s
#endif

        ! Get time until carrier hits energy boundary
        call time_2e(this,t_energy,n_idx)

        ! Compute the free flight time until next scattering event
        r = rng_uniform(rng(this%t_id))
#if DIM == 1
        t_scatter = -1/upper_bound(scat_m(m_region(this%r_idx)),this%valley,idx)*log(r)
#elif DIM == 2
        !t_scatter = -1/upper_bound(this%mat,this%valley,idx,m_region(this%r_idx(1),this%r_idx(2)))*log(r)
        !if (mode == 0) then
   !     if ((idx < 1) .or. (this%valley < 1) .or. (scat_m(m_region(this%r_idx(1),this%r_idx(2))) .ne. 1)) then
        !    print *, this%energy, this%mat, this%r, idx, m_region(this%r_idx(1),this%r_idx(2)), this%r_idx, this%valley
        !    print *, scat_m(m_region(this%r_idx(1),this%r_idx(2)))
        !    this%valley = 0
        !    t_scatter = 100.0_s
  !      else
            t_scatter = -1/upper_bound(scat_m(m_region(this%r_idx(1),this%r_idx(2))),this%valley,idx)*log(r)
 !       end if
        !else
        !    t_scatter = -1/upper_bound(scat_m(m_region(this%res_scat_idx(1),this%res_scat_idx(2))),this%valley,idx)*log(r)
        !end if
!        print *, 'region', this%cidx, this%energy, scat_m(m_region(this%r_idx(1),this%r_idx(2))), &
!            & upper_bound(scat_m(m_region(this%r_idx(1),this%r_idx(2))),this%valley,idx), t_scatter
#elif DIM == 3
        !t_scatter = -1/upper_bound(this%mat,this%valley,idx,m_region(this%r_idx(1),this%r_idx(2),this%r_idx(3)))*log(r)
        t_scatter = -1/upper_bound(scat_m(m_region(this%r_idx(1),this%r_idx(2),this%r_idx(3))),this%valley,idx)*log(r)
#else
        !t_scatter = -1/upper_bound(this%mat,this%valley,idx,1)*log(r)
        t_scatter = -1/upper_bound(1,this%valley,idx)*log(r)
#endif

        ! Select the minimum time and the event
        t_delta = min(t_sync,t_scatter,t_energy,t_real)
        choose_action = minloc([t_sync,t_scatter,t_energy,t_real])
        this%tdelta = (/t_sync,t_scatter,t_energy,t_real/)
        !this%last_action = choose_action(1)
        !print *, t_sync,t_scatter,t_energy,t_real!, this%loc_e_field

        ! Backup starting values, PEP
        !this%pep_r = this%r
        !this%pep_vel = this%vel
        !this%pep_k = this%k
        !this%pep_valley = this%valley
        !this%pep_vsubidx = this%v_subidx
        !this%pep_energy = this%energy

        ! Integrate now equations of motion
        ! Real space, neglect acceleration term
        dr = t_delta*this%vel
       ! dr2 = t_delta**2.0_s*q*e_field/(2.0_s*l_me*(1.0_s+2.0_s*l_nonp*this%energy))
       ! print *, dr, dr2
        this%r = this%r+dr

        ! Reciprocal space
        dk = -q/hb*e_field*t_delta
        if (this%valley == L) then
            if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
                 dk = -q/hb*t_delta*matmul(T_hv_L(this%MAT,1,:,:),e_field)
            end if
            if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
                 dk = -q/hb*t_delta*matmul(T_hv_L(this%MAT,2,:,:),e_field)
            end if
            if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
                 dk = -q/hb*t_delta*matmul(T_hv_L(this%MAT,3,:,:),e_field)
            end if
            if ((this%v_subidx == 7) .or. (this%v_subidx == 8)) then
                 dk = -q/hb*t_delta*matmul(T_hv_L(this%MAT,4,:,:),e_field)
            end if
        end if
        if (this%valley == X) then
            if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
                 dk = -q/hb*t_delta*matmul(T_hv(this%MAT,X,1,:,:),e_field)
            end if
            if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
                 dk = -q/hb*t_delta*matmul(T_hv(this%MAT,X,2,:,:),e_field)
            end if
            if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
                 dk = -q/hb*t_delta*matmul(T_hv(this%MAT,X,3,:,:),e_field)
            end if
        end if
        this%k = this%k+dk
        !this%old_k2 = dk

        ! Update new energy, gained in electric field
        if ((this%valley == X) .or. (this%valley == L)) then
            this%energy = -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(this%k,this%k)/q)/ &
                        & (2.0_s*m0*l_nonp))
        else
            this%energy = -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(this%k,this%k)/q)/ &
                        & (2.0_s*l_me*l_nonp))
        end if
!print *, this%energy
#if DIM == 0
        ! Calculate velocity
        call this%compute_vel()
#endif

        if (choose_action(1) == 1) then
            !call this%check_pep(0, ff_reject)
            !if (ff_reject == 0) then
                this%synchronized = 1
                call get_e_idx(this)
            !end if
        else if (choose_action(1) == 2) then
                call scatter_carrier(this, ff_reject)
                call get_e_idx(this)
                call compute_vel(this)
                !if (ff_reject == 1) then
                !    !print *, 'after scattering rejected', this%energy, this%pep_energy
                !    call this%check_pep(0, ff_reject)
                !    !print *, this%energy
                !end if
            !if (ff_reject == 0) then
            !end if
        else if (choose_action(1) == 3) then
            !call this%check_pep(0, ff_reject)
            !if (ff_reject == 0) then
                if ((n_idx(1) == 1) .or. (n_idx(1) == 2)) then
                    this%e_idx = this%e_idx+1
                else if (((n_idx(1) == 3) .or. (n_idx(1) == 4)) .and. (this%e_idx > 1)) then
                    this%e_idx = this%e_idx-1
                end if
            !end if
        end if

#if DIM == 1
        !if (mode == 0) then
        if (this%r(1) .lt. (0.0_s+1e-5_s)) then
            call action_border_1(this,1)
        elseif (this%r(1) .gt. (dev_x-1e-5_s)) then
            call action_border_1(this,2)
        elseif (choose_action(1) == 4) then
           ! print *, this%cidx, nr_idx, this%r(1), this%r_idx
            call compute_vel(this)
            addx1 = 1
            addx2 = 2
            addy1 = 1
            addy2 = 2
            subx1 = 1
            subx2 = 2
            suby1 = 1
            suby2 = 2
            if (this%r_idx == sum(num_nodes)) then
                addx1 = 0
                addx2 = 0
            end if
            if (this%r_idx == sum(num_nodes)-1) then
                addx2 = 0
            end if
            if (this%r_idx == 1) then
                subx1 = 0
                subx2 = 0
            end if
            if (this%r_idx == 2) then
                subx2 = 0
            end if
            if ((nr_idx(1) == 2) .and. (m_material(m_region(this%r_idx+addx1)) .ne. m_material(m_region(this%r_idx+addx2)))) then
                call handle_heterostructure(this,1,m_material(m_region(this%r_idx+addx1)),m_material(m_region(this%r_idx+addx2)))
            end if
            if ((nr_idx(1) == 1) .and. (m_material(m_region(this%r_idx)) .ne. m_material(m_region(this%r_idx+addx1)))) then
                    call handle_heterostructure(this,1,m_material(m_region(this%r_idx)),m_material(m_region(this%r_idx+addx1)))
            end if
            if ((nr_idx(1) == 3) .and. (m_material(m_region(this%r_idx-subx1)) .ne. m_material(m_region(this%r_idx)))) then
                    call handle_heterostructure(this,2,m_material(m_region(this%r_idx)),m_material(m_region(this%r_idx-subx1)))
            end if
            if ((nr_idx(1) == 4) .and. (m_material(m_region(this%r_idx-subx2)) .ne. m_material(m_region(this%r_idx-subx1)))) then
                call handle_heterostructure(this,2,m_material(m_region(this%r_idx-subx1)),m_material(m_region(this%r_idx-subx2)))
            end if
        end if
        !end if
#endif

#if DIM == 2
        if ((choose_action(1) == 4)) then
!            print *, 'border', this%cidx, nr_idx, this%r(1:2), this%r_idx
        !if (mode == 0) then
            call compute_vel(this)
            addx1 = 1
            addx2 = 2
            addy1 = 1
            addy2 = 2
            subx1 = 1
            subx2 = 2
            suby1 = 1
            suby2 = 2
            if (this%r_idx(1) == num_nodes(1)) then
                addx1 = 0
                addx2 = 0
            end if
            if (this%r_idx(1) == num_nodes(1)-1) then
                addx2 = 0
            end if
            if (this%r_idx(1) == 1) then
                subx1 = 0
                subx2 = 0
            end if
            if (this%r_idx(1) == 2) then
                subx2 = 0
            end if
            if (this%r_idx(2) == num_nodes(2)) then
                addy1 = 0
                addy2 = 0
            end if
            if (this%r_idx(2) == num_nodes(2)-1) then
                addy2 = 0
            end if
            if (this%r_idx(2) == 1) then
                suby1 = 0
                suby2 = 0
            end if
            if (this%r_idx(2) == 2) then
                suby2 = 0
            end if
            !print *, addx1, addx2, addy1, addy2, subx1, subx2, suby1, suby2
            n_mat(1) = m_material(m_region(this%r_idx(1)+addx1,this%r_idx(2)))
            n_mat(2) = m_material(m_region(this%r_idx(1)+addx2,this%r_idx(2)))
            n_mat(5) = m_material(m_region(this%r_idx(1)-subx1,this%r_idx(2)))
            n_mat(6) = m_material(m_region(this%r_idx(1)-subx2,this%r_idx(2)))
            n_mat(3) = m_material(m_region(this%r_idx(1),this%r_idx(2)+addy1))
            n_mat(4) = m_material(m_region(this%r_idx(1),this%r_idx(2)+addy2))
            n_mat(7) = m_material(m_region(this%r_idx(1),this%r_idx(2)-suby1))
            n_mat(8) = m_material(m_region(this%r_idx(1),this%r_idx(2)-suby2))
            n_mat(9) = this%mat !m_material(m_region(this%r_idx(1),this%r_idx(2)))
      !  if (this%cidx == idx_verbose) then
      !     print *, this%cidx, nr_idx, this%r, this%r_idx, n_mat
     !       end if

            ! Handle heterojunctions in x-direction
            if ((nr_idx(1) == 1) .and. (n_mat(9) .ne. n_mat(1))) then
            ! We check tunneling first
                !did_tunnel = par_tunnel(this)
                !if (did_tunnel == 0) then
                    call handle_heterostructure(this,1,n_mat(9),n_mat(1))
                !end if
                cross_hetero = 1
            !end if
            elseif ((nr_idx(1) == 2) .and. (n_mat(9) .ne. n_mat(2))) then
                call handle_heterostructure(this,1,n_mat(9),n_mat(2))
                cross_hetero = 1
            !end if
            elseif ((nr_idx(1) == 5) .and. (n_mat(9) .ne. n_mat(5))) then
                call handle_heterostructure(this,2,n_mat(9),n_mat(5))
                cross_hetero = 1
            !end if
            elseif ((nr_idx(1) == 6) .and. (n_mat(9) .ne. n_mat(6))) then
                call handle_heterostructure(this,2,n_mat(9),n_mat(6))
                cross_hetero = 1
            !end if
            ! Handle heterojunctions in y-direction
            elseif ((nr_idx(1) == 3) .and. (n_mat(9) .ne. n_mat(3))) then
                did_tunnel = par_tunnel(this,1,n_mat(3))
                if (did_tunnel == 0) then
                    call handle_heterostructure(this,3,n_mat(9),n_mat(3))
                end if
                cross_hetero = 1
            !end if
            elseif ((nr_idx(1) == 4) .and. (n_mat(9) .ne. n_mat(4))) then
                did_tunnel = par_tunnel(this,1,n_mat(4))
                if (did_tunnel == 0) then
                    call handle_heterostructure(this,3,n_mat(9),n_mat(4))
                end if
                cross_hetero = 1
            !end if
            elseif ((nr_idx(1) == 7) .and. (n_mat(9) .ne. n_mat(7))) then
                did_tunnel = par_tunnel(this,2,n_mat(7))
                if (did_tunnel == 0) then
                    call handle_heterostructure(this,4,n_mat(9),n_mat(7))
                end if
                cross_hetero = 1
            !end if
            elseif ((nr_idx(1) == 8) .and. (n_mat(9) .ne. n_mat(8))) then
                did_tunnel = par_tunnel(this,2,n_mat(8))
                if (did_tunnel == 0) then
                    call handle_heterostructure(this,4,n_mat(9),n_mat(8))
                end if
                cross_hetero = 1
            end if
        !end if
            ! If no heterostructure handling is required we move carrier in next cell
        if (cross_hetero == 0) then
            if (nr_idx(1) == 1) then
                this%r(1) = this%r(1)+1e-6_s
            elseif (nr_idx(1) == 5) then
                this%r(1) = this%r(1)-1e-6_s
            elseif (nr_idx(1) == 3) then
                this%r(2) = this%r(2)+1e-6_s
            elseif (nr_idx(1) == 7) then
                this%r(2) = this%r(2)-1e-6_s
            end if
        end if
        end if
        if ((mode == 0) .and. ((this%r(1) < 0.0_s+1e-5_s) .or. (this%r(2) < 0.0_s+1e-5_s) .or. &
          & (this%r(1) > dev_x-1e-5_s) .or. (this%r(2) > dev_y-1e-5_s))) then
            call action_border_2(this)
        end if
        !if (this%valley > 0) then
        !    call get_r_idx_2(this%r,this%r_idx)
        !    call this%check_pep(0, ff_reject)
        !end if
#endif

#if DIM == 3
        if ((this%r(1) < 0.0_s+1e-5_s) .or. (this%r(2) < 0.0_s+1e-5_s) .or. &
          & (this%r(1) > dev_x-1e-5_s) .or. (this%r(2) > dev_y-1e-5_s) .or. &
          & (this%r(3) > dev_z-1e-5_s) .or. (this%r(3) < 1e-5_s)) then
          !print *, this%cidx, this%r
            call action_border_3(this)
        end if
#endif
        !end if

        ! Update the propagated time
        !if (ff_reject == 0) then
            this%e_time = this%e_time+t_delta
        !end if
        !print *, 'carrier stop', choose_action, this%r(1:2), this%energy
    end subroutine drift

   subroutine res_drift(this,action_idx)
        class(electron)             :: this
        integer                     :: idx, mode, ff_reject                    ! mode: 0 mc drift, mode: 1 injection (abort drift at border crossing
        integer,dimension(1)        :: choose_action, n_idx
        real(kind=s)                :: t_sync, t_scatter, t_delta, r
        real(kind=s)                :: l_me, l_nonp, e_old
        real(kind=s)                :: t_energy, t_real, el_fermi, el_temp, f_occ, E_tmp
        real(kind=s),dimension(3)   :: dk, dr, e_field, k_tmp
        integer                     :: addx1, addx2, subx1, subx2, addy1, addy2, suby1, suby2
        integer                     :: cross_hetero
        integer,intent(out)         :: action_idx

#if DIM > 0
        integer,dimension(1)        :: nr_idx
#endif
#if DIM == 2
        integer,dimension(9)        :: n_mat
        integer                     :: did_tunnel
!        integer                     :: idx_verbose = 100020
#endif

!#if DIM == 1
!        if (mode == 0) then
!            call get_r_idx(this%r,this%r_idx)
!            call this%assign_force_1
!        end if
!#endif
#if DIM == 2
         call get_e_idx(this)
#endif
!#if DIM == 3
!        call get_r_idx_3(this%r,this%r_idx)
!        call this%assign_force_3
!        !print *, this%loc_e_field
!#endif



#if DIM > 0
        ! Calculate velocity
        call this%compute_vel()
#endif

        !ff_reject = 1
        !cross_hetero = 0
        e_field = this%loc_e_field
        l_me = me(this%mat,this%valley)
        l_nonp = nonp(this%mat,this%valley)
        !e_old = this%energy
        !this%r_old2 = this%r_old
        !this%r_old = this%r
        !this%old_k2 = this%old_k
        !this%old_k = this%k

        idx = E2idx(this%energy)
        if (idx < 1) then
            print *, this%energy, idx, this%r(1:2), this%mat, this%k, this%valley
            !this%valley = 0
            !return
        end if
        if (this%energy > 3.0_s) then
            this%energy = 2.9_s
            idx = E2idx(this%energy)
        end if

        ! Compute the time until next synchronous ensemble step
        t_sync = step*t_step-this%e_time
        !print *, t_sync, step, t_step, this%e_time
        !t_sync = t_step-this%e_time
        ! Calculate the time until a real space boundary is reached
#if DIM == 1
        call time2r_1(t_real,nr_idx,this%r,this%vel,this%r_idx,this%cidx)
#elif DIM == 2
        call time2r_2(t_real,nr_idx,this%r,this%r_idx,this%vel)
#elif DIM == 3
        call time2r_3(t_real,nr_idx,this%r,this%r_idx,this%vel)
#else
        t_real = 100000.0_s
#endif

        ! Get time until carrier hits energy boundary
        call time_2e(this,t_energy,n_idx)

        ! Compute the free flight time until next scattering event
        r = rng_uniform(rng(this%t_id))
#if DIM == 1
        t_scatter = -1/upper_bound(scat_m(m_region(this%r_idx)),this%valley,idx)*log(r)
#elif DIM == 2
        t_scatter = -1/upper_bound(scat_m(m_region(this%res_scat_idx(1),this%res_scat_idx(2))),this%valley,idx)*log(r)
        !if ((idx < 1) .or. (this%valley < 1) .or. (scat_m(m_region(this%r_idx(1),this%r_idx(2))) .ne. 1)) then
        !    print *, this%energy, this%mat, this%r, idx, m_region(this%r_idx(1),this%r_idx(2)), this%r_idx, this%valley
        !    print *, scat_m(m_region(this%r_idx(1),this%r_idx(2)))
        !    !this%valley = 0
        !    t_scatter = 100.0_s
        !end if
        !t_scatter = -1/upper_bound(scat_m(m_region(this%r_idx(1),this%r_idx(2))),this%valley,idx)*log(r)

!        print *, 'region', this%cidx, this%energy, scat_m(m_region(this%r_idx(1),this%r_idx(2))), &
!            & upper_bound(scat_m(m_region(this%r_idx(1),this%r_idx(2))),this%valley,idx), t_scatter
#elif DIM == 3
        !t_scatter = -1/upper_bound(this%mat,this%valley,idx,m_region(this%r_idx(1),this%r_idx(2),this%r_idx(3)))*log(r)
        t_scatter = -1/upper_bound(scat_m(m_region(this%r_idx(1),this%r_idx(2),this%r_idx(3))),this%valley,idx)*log(r)
#else
        !t_scatter = -1/upper_bound(this%mat,this%valley,idx,1)*log(r)
        t_scatter = -1/upper_bound(1,this%valley,idx)*log(r)
#endif

        !print *, e_field
        !print *, this%energy, this%res_scat_idx, m_region(this%res_scat_idx(1),this%res_scat_idx(2)), this%res_id
        !print *, t_sync,t_scatter,t_energy,t_real
        ! Select the minimum time and the event
        t_delta = min(t_sync,t_scatter,t_energy,t_real)
        choose_action = minloc([t_sync,t_scatter,t_energy,t_real])
        this%tdelta = (/t_sync,t_scatter,t_energy,t_real/)
        !this%last_action = choose_action(1)

        ! Integrate now equations of motion
        ! Real space, neglect acceleration term
        dr = t_delta*this%vel
        this%r = this%r+dr

        ! Reciprocal space
        dk = -q/hb*e_field*t_delta
        if (this%valley == L) then
            if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
                 dk = -q/hb*t_delta*matmul(T_hv_L(this%MAT,1,:,:),e_field)
            end if
            if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
                 dk = -q/hb*t_delta*matmul(T_hv_L(this%MAT,2,:,:),e_field)
            end if
            if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
                 dk = -q/hb*t_delta*matmul(T_hv_L(this%MAT,3,:,:),e_field)
            end if
            if ((this%v_subidx == 7) .or. (this%v_subidx == 8)) then
                 dk = -q/hb*t_delta*matmul(T_hv_L(this%MAT,4,:,:),e_field)
            end if
        end if
        if (this%valley == X) then
            if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
                 dk = -q/hb*t_delta*matmul(T_hv(this%MAT,X,1,:,:),e_field)
            end if
            if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
                 dk = -q/hb*t_delta*matmul(T_hv(this%MAT,X,2,:,:),e_field)
            end if
            if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
                 dk = -q/hb*t_delta*matmul(T_hv(this%MAT,X,3,:,:),e_field)
            end if
        end if
        this%k = this%k+dk

        ! Update new energy, gained in electric field
        if ((this%valley == X) .or. (this%valley == L)) then
            this%energy = -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(this%k,this%k)/q)/ &
                        & (2.0_s*m0*l_nonp))
        else
            this%energy = -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(this%k,this%k)/q)/ &
                        & (2.0_s*l_me*l_nonp))
        end if

        if (choose_action(1) == 1) then
                this%synchronized = 1
                call get_e_idx(this)
        else if (choose_action(1) == 2) then
                call scatter_carrier(this, ff_reject)
                call get_e_idx(this)
                if (this%valley < 1) then
                    print *, 3, this%energy, this%r(1:2), this%mat, this%valley
                end if
                call compute_vel(this)
        else if (choose_action(1) == 3) then
                if ((n_idx(1) == 1) .or. (n_idx(1) == 2)) then
                    this%e_idx = this%e_idx+1
                else if (((n_idx(1) == 3) .or. (n_idx(1) == 4)) .and. (this%e_idx > 1)) then
                    this%e_idx = this%e_idx-1
                end if
        end if

#if DIM == 2
        if ((choose_action(1) == 4)) then
            if (nr_idx(1) == 1) then
                this%r(1) = this%r(1)+1e-6_s
            elseif (nr_idx(1) == 5) then
                this%r(1) = this%r(1)-1e-6_s
            elseif (nr_idx(1) == 3) then
                this%r(2) = this%r(2)+1e-6_s
            elseif (nr_idx(1) == 7) then
                this%r(2) = this%r(2)-1e-6_s
            end if
        end if
#endif
        ! Update the propagated time
        !if (ff_reject == 0) then
            this%e_time = this%e_time+t_delta
        !end if
        action_idx = choose_action(1)
    end subroutine res_drift

#if DIM == 0
    subroutine check_pep_0(this, event, rejected)
        class(electron)         :: this
        real(kind=s)            :: el_fermi, el_temp, l_nonp, l_me, E_tmp, f_occ, r
        real(kind=s),dimension(3)   :: k_tmp
        integer                     :: event ! 0: after free flight; 1: after scattering
        integer,intent(out)         :: rejected ! 0: accepted; 1: new state rejected

        r = rng_uniform(rng(this%t_id))
        el_fermi = p%fermi_lvl(this%valley)
        el_temp = p%el_temp(this%valley)
        if ((el_fermi > -2.0_s) .and. (el_temp > -2.0_s)) then
            l_nonp = nonp(this%mat,this%valley)
            l_me = me(this%mat,this%valley)
            f_occ = 1.0_s/(exp((this%energy-el_fermi)/(el_temp*kB))+1.0_s)
        else
            f_occ = 0.0_s
        end if
        rejected = 0
        if ((r < f_occ) .and. (event == 0) .and. (this%op == 0)) then
            rejected = 1
            this%r = this%pep_r
            this%k = this%pep_k
            this%energy = this%pep_energy

        elseif ((r < f_occ) .and. (event == 1) .and. (this%op == 0)) then
            rejected = 1
        end if
    end subroutine check_pep_0
#endif
#if DIM == 2
    subroutine check_pep(this, event, rejected)
        class(electron)         :: this
        real(kind=s)            :: el_fermi, el_temp, l_nonp, l_me, E_tmp, f_occ, r
        real(kind=s),dimension(3)   :: k_tmp
        integer                     :: event ! 0: after free flight; 1: after scattering
        integer,intent(out)         :: rejected ! 0: accepted; 1: new state rejected

        r = rng_uniform(rng(this%t_id))
        el_fermi = p%fermi_lvl(this%valley,this%r_idx(1),this%r_idx(2))
        el_temp = p%el_temp(this%valley,this%r_idx(1),this%r_idx(2))
        if ((el_fermi > -2.0_s) .and. (el_temp > -2.0_s)) then
            l_nonp = nonp(this%mat,this%valley)
            l_me = me(this%mat,this%valley)
            !k_tmp = abs(this%k-p%step_k_drift(p%step,this%valley,this%r_idx(1),this%r_idx(2),:))!p%k_drift(this%valley,this%r_idx(1),this%r_idx(2),:)
            !if ((this%valley == X) .or. (this%valley == L)) then
            !E_tmp = -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(k_tmp,k_tmp)/q)/ &
            !                & (2.0_s*m0*l_nonp))
            !else
            !E_tmp = -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(k_tmp,k_tmp)/q)/ &
            !                & (2.0_s*l_me*l_nonp))
            !end if
            !f_occ = 1.0_s/(exp((E_tmp-el_fermi)/(el_temp*kB))+1.0_s)
            f_occ = 1.0_s/(exp((this%energy-el_fermi)/(el_temp*kB))+1.0_s)
        else
            f_occ = 0.0_s
        end if
!print *, 'check pep'
        rejected = 0
        if ((r < f_occ) .and. (event == 0) .and. (this%op == 0)) then
            rejected = 1
            !print *, f_occ, el_fermi, el_temp, this%energy, this%pep_energy
            this%r = this%pep_r
            this%k = this%pep_k
            this%energy = this%pep_energy
            !this%valley = this%pep_valley
            !this%v_subidx = this%pep_vsubidx
        elseif ((r < f_occ) .and. (event == 1) .and. (this%op == 0)) then
            rejected = 1
        !    print *, f_occ, el_fermi, el_temp, this%energy, this%pep_energy
        end if
    end subroutine check_pep
#endif

    subroutine compute_vel(this)
        class(electron)             :: this
        real(kind=s)                :: l_nonp, l_me

        !print *, this%cidx, this%op, this%r, this%k, this%valley, this%mat

        if (this%valley > 0) then

        l_me = me(this%mat,this%valley)
        l_nonp = nonp(this%mat,this%valley)

        if (this%valley == L) then
            if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
                this%vel = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,1,:,:),this%k)
            end if
            if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
                this%vel = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,2,:,:),this%k)
            end if
            if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
                this%vel = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,3,:,:),this%k)
            end if
            if ((this%v_subidx == 7) .or. (this%v_subidx == 8)) then
                this%vel = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,4,:,:),this%k)
            end if
        else if (this%valley == X) then
            if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
                this%vel = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv(this%MAT,X,1,:,:),this%k)
            end if
            if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
                this%vel = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv(this%MAT,X,2,:,:),this%k)
            end if
            if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
                this%vel = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv(this%MAT,X,3,:,:),this%k)
            end if
        else
           this%vel = hb*this%k/(l_me*(1.0_s+2.0_s*l_nonp*this%energy))
        end if

        else
            print *, this%energy, this%k, this%mat, this%valley, this%vel, this%op, this%r
        end if
    end subroutine compute_vel

#if DIM == 1
    subroutine action_border_1(this,border)
        class(electron)             :: this
        !real(kind=s),dimension(3)   :: vel_tmp, vel_n, k_n
        !real(kind=s)                :: l_nonp, f
        integer                     :: border !, subidx, border

        ! Delete carrier when hitting schottky or ohmic contact
        if ((contact_type(border) .eq. OHMIC) ) then
            this%valley = 0
            particle_out(this%t_id,border) = particle_out(this%t_id,border)+1
            !print *, particle_out(:,1)
            !print *, particle_out(:,2)
            !dev%particle_flow_out(step,border) = dev%particle_flow_out(step,border)+1
        elseif (contact_type(border) .eq. SCHOTTKY) then
            this%valley = 0
            particle_out(this%t_id,border) = particle_out(this%t_id,border)+1
 !           dev%particle_flow_out(step,border) = dev%particle_flow_out(step,border)+1
        else
            ! Left boundary
            if (this%r(1) < 0) then
                this%r(1) = this%r(1)*(-1.0_s)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(1) = this%k(1)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'left border, r:', this%r(1:2), 'v:', this%vel(1:2)
                    call border_L_yz(this)
                end if
            end if
            ! Right
            if (this%r(1) > dev_x) then
                !print *, 'right border, r:', this%valley, this%r(1:2), 'v:', this%vel(1:2)
                this%r(1) = dev_x-(this%r(1)-dev_x)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(1) = this%k(1)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'right border L valley, r:', this%cidx, this%v_subidx, this%r(1:2), 'v:', this%vel(1:2)
                    call border_L_yz(this)
                end if
            end if
            call this%compute_vel()
            !print *, 'after reflection', this%cidx,  this%v_subidx, this%r(1:2), this%vel(1:2)
        end if

!        subidx = 0
!        if (contact_type(border) .eq. INSULATOR) then
!        l_nonp = nonp(this%mat,this%valley)
!        ! Specular reflect momentum and thus velocity
!        if (border .eq. 1) then
!            this%r(1) = this%r(1)*(-1.0_s)
!            if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
!                this%k(1) = this%k(1)*(-1.0_s)
!            else if ((this%valley .eq. L) .and. (this%vel(1) .lt. 0)) then
!                if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
!                    subidx = 1
!                    vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,1,:,:),this%k)
!                end if
!                if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
!                    subidx = 2
!                    vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,2,:,:),this%k)
!                end if
!                if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
!                    subidx = 3
!                    vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,3,:,:),this%k)
!                end if
!                if ((this%v_subidx == 7) .or. (this%v_subidx == 8)) then
!                    subidx = 4
!                    vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,4,:,:),this%k)
!                end if
!                ! Invert x-component of velocity
!                vel_n = vel_tmp*(/-1.0_s,1.0_s,1.0_s/)
!                ! Get corresponding k-vector
!                k_n = matmul(T_hv_L_vel(this%MAT,subidx,:,:),vel_n)/(hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy)))
!                f = sqrt(dot_product(k_n,k_n))/sqrt(dot_product(this%k,this%k))
!                vel_n = vel_n/f
!                this%k = matmul(T_hv_L_vel(this%MAT,subidx,:,:),vel_n)/(hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy)))
!            end if
!        end if
!        if (border .eq. 2) then
!            this%r(1) = dev_x-(this%r(1)-dev_x)
!            if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
!                this%k(1) = this%k(1)*(-1.0_s)
!            else if ((this%valley .eq. L) .and. (this%vel(1) .gt. 0)) then
!                if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
!                    subidx = 1
!                    vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,1,:,:),this%k)
!                end if
!                if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
!                    subidx = 2
!                    vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,2,:,:),this%k)
!                end if
!                if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
!                    subidx = 3
!                    vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,3,:,:),this%k)
!                end if
!                if ((this%v_subidx == 7) .or. (this%v_subidx == 8)) then
!                    subidx = 4
!                    vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,4,:,:),this%k)
!                end if
!                ! Invert x-component of velocity
!                vel_n = vel_tmp*(/-1.0_s,1.0_s,1.0_s/)
!                ! Get corresponding k-vector
!                k_n = matmul(T_hv_L_vel(this%MAT,subidx,:,:),vel_n)/(hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy)))
!                f = sqrt(dot_product(k_n,k_n))/sqrt(dot_product(this%k,this%k))
!                vel_n = vel_n/f
!                this%k = matmul(T_hv_L_vel(this%MAT,subidx,:,:),vel_n)/(hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy)))
!            end if
!        end if

!        if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
!            if ((border .eq. 1) .and. (this%k(1) .lt. 0)) then
!                !print *, border, this%k(1)
!                this%k(1) = this%k(1)*(-1.0_s)
!                this%r(1) = this%r(1)*(-1.0_s)
!            else if (border .eq. 1) then
!                this%r(1) = this%r(1)*(-1.0_s)
!            end if
!            if ((border .eq. 2) .and. (this%k(1) .gt. 0)) then
!                this%k(1) = this%k(1)*(-1.0_s)
!                this%r(1) = dev_x-(this%r(1)-dev_x)
!            else if (border .eq. 2) then
!                this%r(1) = dev_x-(this%r(1)-dev_x)
!            end if
!        end if
!        if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
!            if (this%r_idx .eq. 0) then
!                this%k(1) = this%k(1)*(-1.0_s)
!                !this%vel(1) = this%vel(1)*(-1.0_s)
!                this%r_idx = 1
!                if ((this%r(1) .lt. 0) .and. (abs(this%r(1)) .lt. 1e-10)) then
!                    this%r(1) = 0.0_s
!                end if
!            end if
!            if (this%r_idx .eq. (sum(num_nodes)+1)) then
!                this%k(1) = this%k(1)*(-1.0_s)
!                !this%vel(1) = this%vel(1)*(-1.0_s)
!                this%r_idx = sum(num_nodes)
!                !print *, 'jo'
!            end if
!        ! We choose the new k-vector in the L-valleys according to a velocity which gives a specular
!        ! reflections at the plane
!        if (this%valley .eq. L) then
!            if (((border .eq. 1) .and. (this%vel(1) .lt. 0)) .or. &
!             & ((border .eq. 2) .and. (this%vel(1) .gt. 0))) then
!            this%r(1) = this%r(1)*(-1.0_s)
!            if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
!                subidx = 1
!                vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,1,:,:),this%k)
!            end if
!            if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
!                subidx = 2
!                vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,2,:,:),this%k)
!            end if
!            if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
!                subidx = 3
!                vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,3,:,:),this%k)
!            end if
!            if ((this%v_subidx == 7) .or. (this%v_subidx == 8)) then
!                subidx = 4
!                vel_tmp = hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy))*matmul(T_hv_L_t(this%MAT,4,:,:),this%k)
!            end if
!            ! Invert x-component of velocity
!            vel_n = vel_tmp*(/-1.0_s,1.0_s,1.0_s/)
!            ! Get corresponding k-vector
!            k_n = matmul(T_hv_L_vel(this%MAT,subidx,:,:),vel_n)/(hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy)))
!            f = sqrt(dot_product(k_n,k_n))/sqrt(dot_product(this%k,this%k))
!            vel_n = vel_n/f
!            this%k = matmul(T_hv_L_vel(this%MAT,subidx,:,:),vel_n)/(hb/(m0*(1.0_s+2.0_s*l_nonp*this%energy)))
!            end if
!        end if
            !if (abs(this%r(1)) .lt. 1e-8_s) then
            !    this%r(1) = 0.0_s
            !end if
            !if (abs(this%r(1)-dev_x) .lt. 1e-8_s) then
            !    this%r(1) = dev_x
            !end if
            !print *, 'wrong way', SCHOTTKY
!            end if
            ! Delete carrier when hitting schottky or ohmic contact
!            if (contact_type(border) .eq. OHMIC) then
!                dev%particle_flow_out(step,border) = dev%particle_flow_out(step,border)+1
!                this%valley = 0
!            end if
!            if (contact_type(border) .eq. SCHOTTKY) then
!                dev%particle_flow_out(step,border) = dev%particle_flow_out(step,border)+1
!                this%valley = 0
!                !print *, 'schottky', this%r(1), this%valley, this%cidx
!            end if
    end subroutine action_border_1
#endif

#if DIM == 2
    subroutine action_border_2(this)
        class(electron)             :: this
        integer                     :: c_id

        !s1 = int(mod(step,save_t)+1,4)
        !print *, s1, step, save_t

        call get_r_idx_2(this%r,this%r_idx)

        c_id = contact_id(this%r_idx(1),this%r_idx(2))

        ! Delete carrier when hitting schottky or ohmic contact
        if ((contact_type(this%r_idx(1),this%r_idx(2)) .eq. OHMIC)) then! .and. ((this%r(1) > 0.0_s) .and. (this%r(1) < dev_x))) then
            this%valley = 0
            particle_out(this%t_id,c_id) = particle_out(this%t_id,c_id)+1
        elseif (contact_type(this%r_idx(1),this%r_idx(2)) .eq. SCHOTTKY) then
            this%valley = 0
            !print *, this%r_idx, c_id
            particle_out(this%t_id,c_id) = particle_out(this%t_id,c_id)+1
        else
            ! Left boundary
            if (this%r(1) < 0.0_s) then
                this%r(1) = this%r(1)*(-1.0_s)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(1) = this%k(1)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'left border, r:', this%r(1:2), 'v:', this%vel(1:2)
                    call border_L_yz(this)
                end if
            end if
            ! Right
            if (this%r(1) > dev_x) then
                !print *, 'right border, r:', this%valley, this%r(1:2), 'v:', this%vel(1:2)
                this%r(1) = dev_x-(this%r(1)-dev_x)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(1) = this%k(1)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'right border L valley, r:', this%valley, this%r(1:2), 'v:', this%vel(1:2)
                    call border_L_yz(this)
                end if
            end if
            ! Bottom
            if (this%r(2) < 0.0_s) then
                this%r(2) = this%r(2)*(-1.0_s)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(2) = this%k(2)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'bottom border, r:', this%cidx, this%valley, this%v_subidx, this%r(1:2), &
                    !                & 'k:', this%k(1:3), 'v:', this%vel(1:2)
                    call border_L_xz(this)
                    !call this%compute_vel()
                    !print *, 'bottom ref', this%cidx, this%valley, this%v_subidx, this%r(1:2), 'k:', this%k(1:3), &
                    !                & 'v:', this%vel(1:2)
                end if
            end if
            ! Top
            if (this%r(2) > dev_y) then
                this%r(2) = dev_y-(this%r(2)-dev_y)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(2) = this%k(2)*(-1.0_s)
                end if
                if (this%valley == L) then
                   !print *, 'top border, r:', this%r(1:2), 'v:', this%vel(1:2)
                    call border_L_xz(this)
                end if
            end if
            call this%compute_vel()
            !print *, 'after reflection', this%cidx,  this%r(1:2), this%vel(1:2)
        end if
    end subroutine action_border_2
#endif

#if DIM == 3
    subroutine action_border_3(this)
        class(electron)             :: this
        integer                     :: orientation, coord1, coord2, num_xy

        call get_r_idx_3(this%r,this%r_idx)
        call coord2idx(this%r_idx(1),this%r_idx(2),this%r_idx(3),orientation,coord1,coord2)

        ! Delete carrier when hitting schottky or ohmic contact
        if (contact_type(orientation,coord1,coord2) .eq. OHMIC) then
            this%valley = 0
        elseif (contact_type(orientation,coord1,coord2) .eq. SCHOTTKY) then
            this%valley = 0
        else
            ! Left boundary
            if (this%r(1) < 0) then
                this%r(1) = this%r(1)*(-1.0_s)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(1) = this%k(1)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'left border, r:', this%r(1:2), 'v:', this%vel(1:2)
                    call border_L_yz(this)
                end if
            end if
            ! Right
            if (this%r(1) > dev_x) then
                !print *, 'right border, r:', this%valley, this%r(1:2), 'v:', this%vel(1:2)
                this%r(1) = dev_x-(this%r(1)-dev_x)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(1) = this%k(1)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'right border L valley, r:', this%valley, this%r(1:2), 'v:', this%vel(1:2)
                    call border_L_yz(this)
                end if
            end if
            ! Bottom
            if (this%r(2) < 0) then
                this%r(2) = this%r(2)*(-1.0_s)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(2) = this%k(2)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'bottom border, r:', this%cidx, this%valley, this%v_subidx, this%r(1:2), &
                    !                & 'k:', this%k(1:3), 'v:', this%vel(1:2)
                    call border_L_xz(this)
                    !call this%compute_vel()
                    !print *, 'bottom ref', this%cidx, this%valley, this%v_subidx, this%r(1:2), 'k:', this%k(1:3), &
                    !                & 'v:', this%vel(1:2)
                end if
            end if
            ! Top
            if (this%r(2) > dev_y) then
                this%r(2) = dev_y-(this%r(2)-dev_y)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(2) = this%k(2)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'top border, r:', this%r(1:2), 'v:', this%vel(1:2)
                    call border_L_xz(this)
                end if
            end if
            ! Front boundary
            if (this%r(3) < 0) then
                this%r(3) = this%r(3)*(-1.0_s)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(3) = this%k(3)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'front border, r:', this%r(1:2), 'v:', this%vel(1:2)
                    call border_L_xy(this)
                end if
            end if
            ! Back
            if (this%r(3) > dev_z) then
                !print *, 'right border, r:', this%valley, this%r(1:2), 'v:', this%vel(1:2)
                this%r(3) = dev_z-(this%r(3)-dev_z)
                if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                    this%k(3) = this%k(3)*(-1.0_s)
                end if
                if (this%valley == L) then
                    !print *, 'back border L valley, r:', this%valley, this%r(1:2), 'v:', this%vel(1:2)
                    call border_L_xy(this)
                end if
            end if
            call this%compute_vel()
            !print *, 'after reflection', this%cidx,  this%r(1:3), this%vel(1:3)
        end if
    end subroutine action_border_3
#endif

    subroutine border_L_yz(this)
        class(electron)             :: this
        real(kind=s),dimension(3)   :: k_org

        select case (this%v_subidx)
        case ( 1,2 )
            k_org = matmul(transpose(D(this%mat,1,:,:)),this%k)
            k_org(1) = k_org(1)*(-1.0_s)
            this%k = matmul(D(this%mat,2,:,:),k_org)
            this%v_subidx = 3
        case ( 3,4 )
            k_org = matmul(transpose(D(this%mat,2,:,:)),this%k)
            k_org(1) = k_org(1)*(-1.0_s)
            this%k = matmul(D(this%mat,1,:,:),k_org)
            this%v_subidx = 1
        case ( 5,6 )
            k_org = matmul(transpose(D(this%mat,3,:,:)),this%k)
            k_org(1) = k_org(1)*(-1.0_s)
            this%k = matmul(D(this%mat,4,:,:),k_org)
            this%v_subidx = 7
        case ( 7,8 )
            k_org = matmul(transpose(D(this%mat,4,:,:)),this%k)
            k_org(1) = k_org(1)*(-1.0_s)
            this%k = matmul(D(this%mat,3,:,:),k_org)
            this%v_subidx = 5
        end select
    end subroutine border_L_yz

    subroutine border_L_xz(this)
        class(electron)             :: this
        real(kind=s),dimension(3)   :: k_org

        select case (this%v_subidx)
        case ( 1,2 )
            k_org = matmul(transpose(D(this%mat,1,:,:)),this%k)
            k_org(2) = k_org(2)*(-1.0_s)
            this%k = matmul(D(this%mat,4,:,:),k_org)
            this%v_subidx = 7
        case ( 3,4 )
            k_org = matmul(transpose(D(this%mat,2,:,:)),this%k)
            k_org(2) = k_org(2)*(-1.0_s)
            this%k = matmul(D(this%mat,3,:,:),k_org)
            this%v_subidx = 5
        case ( 5,6 )
            k_org = matmul(transpose(D(this%mat,3,:,:)),this%k)
            k_org(2) = k_org(2)*(-1.0_s)
            this%k = matmul(D(this%mat,2,:,:),k_org)
            this%v_subidx = 3
        case ( 7,8 )
            k_org = matmul(transpose(D(this%mat,4,:,:)),this%k)
            k_org(2) = k_org(2)*(-1.0_s)
            this%k = matmul(D(this%mat,1,:,:),k_org)
            this%v_subidx = 1
        end select
    end subroutine

    subroutine border_L_xy(this)
        class(electron)             :: this
        real(kind=s),dimension(3)   :: k_org

        select case (this%v_subidx)
        case ( 1,2 )
            k_org = matmul(transpose(D(this%mat,1,:,:)),this%k)
            k_org(3) = k_org(3)*(-1.0_s)
            this%k = matmul(D(this%mat,3,:,:),k_org)
            this%v_subidx = 5
        case ( 3,4 )
            k_org = matmul(transpose(D(this%mat,2,:,:)),this%k)
            k_org(3) = k_org(3)*(-1.0_s)
            this%k = matmul(D(this%mat,4,:,:),k_org)
            this%v_subidx = 7
        case ( 5,6 )
            k_org = matmul(transpose(D(this%mat,3,:,:)),this%k)
            k_org(3) = k_org(3)*(-1.0_s)
            this%k = matmul(D(this%mat,1,:,:),k_org)
            this%v_subidx = 1
        case ( 7,8 )
            k_org = matmul(transpose(D(this%mat,4,:,:)),this%k)
            k_org(3) = k_org(3)*(-1.0_s)
            this%k = matmul(D(this%mat,2,:,:),k_org)
            this%v_subidx = 3
        end select
    end subroutine border_L_xy

    subroutine scatter_carrier(this, scat_rejected)
        class(electron)             :: this
        integer                     :: scat_idx, idx, i
        real(kind=s)                :: r, el_fermi, el_temp, f_occ, E_tmp, l_nonp, l_me
        integer                     :: m_idx
        integer,intent(out)         :: scat_rejected
        real(kind=s),dimension(3)   :: pep_k
        integer                     :: pep_valley, pep_vsubidx
        real(kind=s)                :: pep_energy

        ! Backup values before scattering event
        pep_k = this%k
        pep_valley = this%valley
        pep_vsubidx = this%v_subidx
        pep_energy = this%energy

        if (this%energy .gt. 3.0_s) then
        !    print *, '----------------------------'
        !    print *, 'energy larger than 3 eV', this%cidx, this%energy
            this%energy = 2.99999
        !    print *, '----------------------------'
        end if
        if (this%energy .lt. 0.0_s) then
        !    print *, '----------------------------'
        !    print *, 'energy larger than 3 eV', this%cidx, this%energy
            this%energy = 1e-5_s
        !    print *, '----------------------------'
        end if
#if DIM == 0
        m_idx = 1
#elif DIM == 1
        m_idx = int(m_region(this%r_idx))
#elif DIM == 2
        if (this%op == 0) then
            call get_r_idx_2(this%r,this%r_idx)
            m_idx = int(m_region(this%r_idx(1),this%r_idx(2)))
        else
            m_idx = int(m_region(this%res_scat_idx(1),this%res_scat_idx(2)))
            !print *, this%op, m_idx, this%res_scat_idx
        end if
#elif DIM == 3
        m_idx = int(m_region(this%r_idx(1),this%r_idx(2),this%r_idx(3)))
#endif
        idx = E2idx(this%energy)

        r = rng_uniform(rng(this%t_id))
        if (r == 1.0_s) then
            r = 0.999999999999_s
        end if

        if (r <= n_scat_rates(scat_m(m_idx),this%valley,idx,1)) then
            scat_idx = 1
        else
            do i=1,14
                if ((r > n_scat_rates(scat_m(m_idx),this%valley,idx,i)) .and. &
                  & (r <= n_scat_rates(scat_m(m_idx),this%valley,idx,i+1))) then
                    scat_idx = i+1
                    exit
                end if
            end do
        end if
        !this%last_scat_idx(1) = this%last_scat_idx(2)
        !this%last_scat_idx(2) = scat_idx
        !this%last_e(1) = this%last_e(2)
        !this%last_e(2) = this%energy
        !this%last_k(1,:) = this%last_k(2,:)
        !this%last_k(2,:) = this%k
        call perform_scattering(this,scat_idx)
#if DIM == 0
        if (pep_hdf5 == 1) then
            call this%check_pep_0(1, scat_rejected)
            if (scat_rejected == 1) then
                !print *, 'scattering rejected', this%k, pep_k, this%energy, pep_energy
                this%k = pep_k
                this%valley = pep_valley
                this%v_subidx = pep_vsubidx
                this%energy = pep_energy
            end if
        end if
#endif
#if DIM == 2
        if (pep_hdf5 == 1) then
            call this%check_pep(1, scat_rejected)
            if (scat_rejected == 1) then
                !print *, 'scattering rejected', this%k, pep_k, this%energy, pep_energy
                this%k = pep_k
                this%valley = pep_valley
                this%v_subidx = pep_vsubidx
                this%energy = pep_energy
            end if
            !else
            !    print *, 'scattering accepted', this%energy, pep_energy
            !end if
        end if
#endif
    end subroutine scatter_carrier

    subroutine perform_scattering(this,scat_idx)
        class(electron)             :: this
        integer                     :: scat_idx

        ! Handle the new k-vector and energy, according to the scattering type
        if (scat_idx == 1) then
            ! Self scattering: do nothing :-)
        ! Elastic acoustic phonon scattering
        !else if (scat_idx == 2) then
        !    call scatter_isotropic(this)
        ! Inelastic acoustic phonon scattering
        else if (scat_idx == 2) then
            call scatter_ac(this,2)
        else if (scat_idx == 3) then
            call scatter_ac(this,1)
        ! Nonpolar optical phonon emission
        else if (scat_idx == 4) then
            this%energy = this%energy-phonon_e(this%mat,1)
            call scatter_isotropic(this)
        ! Nonpolar optical phonon absorption
        else if (scat_idx == 5) then
            this%energy = this%energy+phonon_e(this%mat,1)
            call scatter_isotropic(this)
        ! Polar optical phonon emission/absorption
        else if ((scat_idx >= 6) .and. (scat_idx <= 7)) then
            !print *, 'Old:', this%energy, sqrt(dot_product(this%k,this%k))
            call scatter_pop(this,scat_idx)
        ! Intervalley scattering
        else if ((scat_idx >= 8) .and. (scat_idx <= 13)) then
            call scatter_iv(this,scat_idx)
        ! Alloy scattering
        else if (scat_idx == 14) then
            call scatter_isotropic(this)
        ! Impurity scattering
        else if (scat_idx == 15) then
            call scatter_imp(this)
        end if
    end subroutine perform_scattering

    ! Function which implements the IV scattering and calls scatter_isotropic, type: see index in scat_rates!
    subroutine scatter_iv(this,type)
        class(electron)            :: this
        integer                    :: old_val, type
        real(kind=s)               :: diff_LX

        diff_LX = valley_offs(this%mat,X)-valley_offs(this%mat,L)
        old_val = this%valley
        if (old_val == G) then
            !print *, "G-valley:", this%valley, this%energy
            if (type == 8) then
                this%energy = this%energy-phonon_e(this%mat,2)-valley_offs(this%mat,L)
                this%v_subidx = new_valley_subidx(this,L,this%v_subidx)
                this%valley = L
            else if (type == 9) then
                this%energy = this%energy+phonon_e(this%mat,2)-valley_offs(this%mat,L)
                this%v_subidx = new_valley_subidx(this,L,this%v_subidx)
                this%valley = L
            else if (type == 10) then
                this%energy = this%energy-phonon_e(this%mat,3)-valley_offs(this%mat,X)
                this%v_subidx = new_valley_subidx(this,X,this%v_subidx)
                this%valley = X
            else if (type == 11) then
                this%energy = this%energy+phonon_e(this%mat,3)-valley_offs(this%mat,X)
                this%v_subidx = new_valley_subidx(this,X,this%v_subidx)
                this%valley = X
            end if
            !print *, "G-valley: new", this%valley, this%energy
        else if (old_val == L) then
            !print *, "L-valley:", this%valley, this%energy
            if (type == 8) then
                this%energy = this%energy-phonon_e(this%mat,2)+valley_offs(this%mat,L)
                this%v_subidx = 0
                this%valley = G
            else if (type == 9) then
                this%energy = this%energy+phonon_e(this%mat,2)+valley_offs(this%mat,L)
                this%v_subidx = 0
                this%valley = G
            else if (type == 10) then
                this%energy = this%energy-phonon_e(this%mat,5)
                this%v_subidx = new_valley_subidx(this,L,this%v_subidx)
            else if (type == 11) then
                this%energy = this%energy+phonon_e(this%mat,5)
                this%v_subidx = new_valley_subidx(this,L,this%v_subidx)
            else if (type == 12) then
                this%energy = this%energy-phonon_e(this%mat,4)-diff_LX
                this%v_subidx = new_valley_subidx(this,X,this%v_subidx)
                this%valley = X
            else if (type == 13) then
                this%energy = this%energy+phonon_e(this%mat,4)-diff_LX
                this%v_subidx = new_valley_subidx(this,X,this%v_subidx)
                this%valley = X
            end if
            !print *, "L-valley: new", this%valley, this%energy
        else if (old_val == X) then
            !print *, "X-valley:", this%valley, this%energy
            if (type == 8) then
                this%energy = this%energy-phonon_e(this%mat,3)+valley_offs(this%mat,X)
                this%valley = G
                this%v_subidx = 0
            else if (type == 9) then
                this%energy = this%energy+phonon_e(this%mat,3)+valley_offs(this%mat,X)
                this%valley = G
                this%v_subidx = 0
            else if (type == 10) then
                this%energy = this%energy-phonon_e(this%mat,4)+diff_LX
                this%v_subidx = new_valley_subidx(this,L,this%v_subidx)
                this%valley = L
            else if (type == 11) then
                this%energy = this%energy+phonon_e(this%mat,4)+diff_LX
                this%v_subidx = new_valley_subidx(this,L,this%v_subidx)
                this%valley = L
            else if (type == 12) then
                this%energy = this%energy-phonon_e(this%mat,6)
                this%v_subidx = new_valley_subidx(this,X,this%v_subidx)
            else if (type == 13) then
                this%energy = this%energy+phonon_e(this%mat,6)
                this%v_subidx = new_valley_subidx(this,X,this%v_subidx)
            end if
        end if

        ! Randomize momentum in new valley
        call scatter_isotropic(this)
    end subroutine scatter_iv

    ! Get new valley index and subindex (different valleys of same type)
    function new_valley_subidx(this,new_valley,old_subidx) result (new_subidx)
        class(electron)             :: this
        integer                     :: new_valley, old_subidx, new_subidx, n
        !real(kind=s)                :: r

        new_subidx = 0

        if ((this%valley .ne. X) .and. (new_valley == X)) then
                new_subidx = int(rand(0)*6)+1
        end if
        if ((this%valley .ne. L) .and. (new_valley == L)) then
                new_subidx = int(rand(0)*8)+1
        end if

        n = 1
        if ((this%valley == X) .and. (new_valley == X)) then
            do while (n > 0)
                new_subidx = int(rand(0)*6)+1
                if (new_subidx .ne. old_subidx) then
                    n = 0
                end if
            end do
        end if
        n = 1
        if ((this%valley == L) .and. (new_valley == L)) then
            do while (n > 0)
                new_subidx = int(rand(0)*8)+1
                if (new_subidx .ne. old_subidx) then
                    n = 0
                end if
            end do
        end if
    end function new_valley_subidx

    ! Implementation of the isotropic scattering case, k-vector is randomized in current valley
    subroutine scatter_isotropic(this)
        class(electron)             :: this
        real(kind=s)                :: r1, r2, k_new_abs, ct

        if (this%valley .eq. G) then
        !if (this%valley .ne. X) then
            k_new_abs = sqrt(2.0_s*me(this%mat,this%valley)*this%energy*q*(1.0_s+nonp(this%mat,this%valley)*this%energy))/hb
        !end if
        !if ((this%valley .eq. X) .or. (this%valley .eq. L)) then
        else
            k_new_abs = sqrt(2.0_s*m0*this%energy*q*(1.0_s+nonp(this%mat,this%valley)*this%energy))/hb
        end if
        !call random_number(r1)
        r1 = rng_uniform(rng(this%t_id))
        !call random_number(r2)
        r2 = rng_uniform(rng(this%t_id))
        ct = (1.0_s-2.0_s*r1)
        if (ct >= 1.0_s) then
            print *, 'ct corrected iso', ct
            ct = 0.99999_s
        end if
        if (ct <= -1.0_s) then
            print *, 'ct corrected iso', ct
            ct = -0.99999_s
        end if
        if (ct == 0.0_s) then
            ct = 0.000001_s
            print *, 'ct corrected iso', ct
        end if
        this%k(1) = k_new_abs*sqrt(1.0_s-ct*ct)*cos(2.0_s*pi*r2)
        this%k(2) = k_new_abs*sqrt(1.0_s-ct*ct)*sin(2.0_s*pi*r2)
        this%k(3) = k_new_abs*ct
        !if (abs(ct) <= 1e-5_s) then
        !    print *, 'isotropic ct', ct, this%k
        !end if
        !if ((this%k(1) == 0.0_s) .or. (this%k(2) == 0.0_s) .or. (this%k(3) == 0.0_s)) then
        !    print *, 'isotropic k = 0', ct, this%k, r1, r2
        !end if
        !if (abs(this%k(1)) < 1e-10_s) then
        !    this%k(1) = 1e-10_s
        !end if
        !if (abs(this%k(2)) < 1e-10_s) then
        !    this%k(2) = 1e-10_s
        !end if
        !if (abs(this%k(3)) < 1e-10_s) then
        !    this%k(3) = 1e-10_s
        !end if
    end subroutine scatter_isotropic

    ! Small angle scattering needed for polar optical phonons, angle beta is calculated using the rejection method
    subroutine scatter_pop(this,scat_idx)
        class(electron)             :: this
        integer                     :: scat_idx, n
        real(kind=s)                :: r1, r2, k_new_abs, kxy, k, cth0, sth0, cfi0, sfi0, costheta, sb, e_new, e_old !, zeta
        real(kind=s)                :: Gamma_i, Gamma_n, beta, Ppop_max, Ppop
        real(kind=s),dimension(3)   :: x_new
        !real(kind=s),dimension(100) :: s1, s2, beta_s, Ppop

        n = 1
        e_old = this%energy
        if (scat_idx == 6) then
            e_new = e_old-phonon_e(this%mat,1)
        else
            e_new = e_old+phonon_e(this%mat,1)
        end if
        this%energy = e_new
        !if (this%energy <= 0.0_s) then
        !    this%energy = 1e-5_s
        !    print *, 'pop energy smaller 0', this%energy, this%r(1:2), E2idx(this%energy)
        !    print *, this%cidx, scat_idx, this%energy, this%k, this%vel, this%t_id
        !    print *, this%cidx, this%last_scat_idx(1), this%last_e(1), this%last_k(1,:)
        !    print *, this%cidx, this%last_scat_idx(2), this%last_e(2), this%last_k(2,:)
        !end if
        if (this%valley .eq. G) then
        !if (this%valley .ne. X) then
        !if (2.0_s*me(this%mat,this%valley)*this%energy*q*(1.0_s+nonp(this%mat,this%valley)*this%energy) <= 0) then
        !    print *, this%energy, this%valley, this%mat
        !end if
            k_new_abs = sqrt(2.0_s*me(this%mat,this%valley)*this%energy*q*(1.0_s+nonp(this%mat,this%valley)*this%energy))/hb
        !end if
        !if ((this%valley .eq. X) .or. (this%valley .eq. L)) then
        else
            k_new_abs = sqrt(2.0_s*m0*this%energy*q*(1.0_s+nonp(this%mat,this%valley)*this%energy))/hb
        end if
        ! Rotate coordinate system
        kxy = sqrt(this%k(1)*this%k(1)+this%k(2)*this%k(2))
        k = sqrt(kxy*kxy+this%k(3)*this%k(3))
        !if ((this%k(1) == 0.0_s) .or. (this%k(2) == 0.0_s) .or. (this%k(3) == 0.0_s)) then
        !    print *, this%cidx, scat_idx, this%energy, this%k, this%vel, this%t_id
        !    print *, this%cidx, this%last_scat_idx(1), this%last_e(1), this%last_k(1,:)
        !    print *, this%cidx, this%last_scat_idx(2), this%last_e(2), this%last_k(2,:)
        !end if
        !if (kxy == 0) then
        !    print *, 'kxy', kxy
        !end if
        !if (k == 0) then
        !    print *, 'k', k
        !end if
        cth0 = this%k(3)/k
        sth0 = kxy/k
        cfi0 = this%k(1)/kxy
        sfi0 = this%k(2)/kxy

        ! Determine scattering angle - using rejection method
        Gamma_i = e_old*(1.0_s+nonp(this%mat,this%valley)*e_old)
        Gamma_n = e_new*(1.0_s+nonp(this%mat,this%valley)*e_new)

        !if (e_new <= 0) then
        !    print *, this%energy, k_new_abs, this%k, Gamma_i, Gamma_n
        !end if

        ! Following the calculation of Fawcett
        Ppop_max = (sqrt(Gamma_i)*sqrt(Gamma_n)+nonp(this%MAT,this%valley)*e_old*e_new*1.0_s)**2.0_s/ &
                 & (Gamma_i+Gamma_n-2.0_s*sqrt(Gamma_i)*sqrt(Gamma_n)*1.0_s)
!        do while (n > 0)
!            call random_number(s1)
!            s1 = s1*2.0_s-1.0_s
!            call random_number(s2)
!            s2 = s2*Ppop_max
!            beta_s = acos(s1)
!            Ppop = (sqrt(Gamma_i)*sqrt(Gamma_n)+nonp(this%MAT,this%valley)*e_old*e_new*cos(beta_s))**2.0_s/ &
!                 & (Gamma_i+Gamma_n-2.0_s*sqrt(Gamma_i)*sqrt(Gamma_n)*cos(beta_s))
!            !print *, Ppop(1:2), beta_s(1:2)
!            do i=1,100
!                if (s2(i) .le. Ppop(i)) exit
!                !print *, i, s2(i), Ppop(i)
!            end do
!            if (i .le. 100) then
!                n = -1
!                costheta = s1(i)
!                !print *, i, s2(i), Ppop(i)
!            end if
!            !print *, n
!        end do
        !print *, costheta
        do while (n > 0)
            !call random_number(r1)
            r1 = rng_uniform(rng(this%t_id))
            r1 = r1*2.0_s-1.0_s
            !call random_number(r2)
            r2 = rng_uniform(rng(this%t_id))
            r2 = r2*Ppop_max
            beta = acos(r1)
            Ppop = (sqrt(Gamma_i)*sqrt(Gamma_n)+nonp(this%MAT,this%valley)*e_old*e_new*cos(beta))**2.0_s/ &
                 & (Gamma_i+Gamma_n-2.0_s*sqrt(Gamma_i)*sqrt(Gamma_n)*cos(beta))
            if (r2 <= Ppop) then
                n = -1
                costheta = r1
            end if
        end do
!       Evaluation of the scattering angle using a first guess of the PDF to find the maximum value
!       needed for the rejection method
!        tmp = calc_Ppop(this,e_old,e_new,Gamma_i,Gamma_n)
!        tmp(1) = maxval(tmp)
!        tmp(1) = tmp(1)+tmp(1)
!
!        do while (n > 0)
!            call random_number(r1)
!            call random_number(r2)
!            r2 = r2*tmp(1)
!            beta = r1*pi
!            Ppop = (sqrt(Gamma_i)*sqrt(Gamma_n)+nonp(this%MAT,this%valley)*e_old*e_new*cos(beta))**2.0_s/ &
!                 & (Gamma_i+Gamma_n-2.0_s*sqrt(Gamma_i)*sqrt(Gamma_n)*cos(beta))*sin(beta)
!            if (r2 < Ppop) then
!                n = -1
!                costheta = cos(beta)
!            end if
!        end do
!       Calculation without the rejection method, constant angle distribution function
!        call random_number(r1)
!        zeta = 2.0_s*sqrt(e_old*e_new)/((sqrt(e_old)-sqrt(e_new))*(sqrt(e_old)-sqrt(e_new)))
!        costheta = ((1.0_s+zeta)-(1.0_s+2.0_s*zeta)**r1)/zeta

        if (costheta >= 1.0_s) then
            print *, 'costheta corrected pop', costheta
            costheta = 0.99999_s
        end if
        if (costheta <= -1.0_s) then
            print *, 'costheta corrected pop', costheta
            costheta = -0.99999_s
        end if
        if (costheta == 0.0_s) then
            print *, 'costheta corrected pop', costheta
            costheta = 0.00001_s
        end if
        sb = sqrt(1.0_s-costheta*costheta)

        !call random_number(r1)
        r1 = rng_uniform(rng(this%t_id))
        x_new(1) = k_new_abs*sb*cos(2.0_s*pi*r1)
        x_new(2) = k_new_abs*sb*sin(2.0_s*pi*r1)
        x_new(3) = k_new_abs*costheta

        ! Transform back
        this%k(1) = x_new(1)*cfi0*cth0-x_new(2)*sfi0+x_new(3)*cfi0*sth0
        this%k(2) = x_new(1)*sfi0*cth0+x_new(2)*cfi0+x_new(3)*sfi0*sth0
        this%k(3) = -x_new(1)*sth0+x_new(3)*cth0
    end subroutine scatter_pop

    ! Small angle scattering needed for polar optical phonons, angle beta is calculated using the rejection method
    subroutine scatter_imp(this) !,costheta)
        class(electron)             :: this
        real(kind=s)                :: r1, k_new_abs, kxy, k, cth0, sth0, cfi0, sfi0, costheta, sb
        real(kind=s)                :: Gamma_i, beta2, e_beta
        real(kind=s),dimension(3)   :: x_new

        if (this%valley .eq. G) then
            k_new_abs = sqrt(2.0_s*me(this%mat,this%valley)*this%energy*q*(1.0_s+nonp(this%mat,this%valley)*this%energy))/hb
        !end if
        !if ((this%valley .eq. X) .or. (this%valley .eq. L)) then
        else
            k_new_abs = sqrt(2.0_s*m0*this%energy*q*(1.0_s+nonp(this%mat,this%valley)*this%energy))/hb
        end if

        ! Rotate coordinate system
        kxy = sqrt(this%k(1)*this%k(1)+this%k(2)*this%k(2))
        k = sqrt(kxy*kxy+this%k(3)*this%k(3))
        cth0 = this%k(3)/k
        sth0 = kxy/k
        cfi0 = this%k(1)/kxy
        sfi0 = this%k(2)/kxy

        ! Determine scattering angle - use only the Brooks-Herring angular distribution
        Gamma_i = this%energy*(1.0_s+nonp(this%mat,this%valley)*this%energy)

        ! Calculation without the rejection method, constant angle distribution function
        !call random_number(r1)
        r1 = rng_uniform(rng(this%t_id))
#if DIM == 0
        beta2 = donor_density(1)*q**2.0_s/(eps_static(this%MAT)*kB*T_lattice)
#elif DIM == 1
        beta2 = donor(this%r_idx)*q**2.0_s/(eps_static(this%MAT)*kB*T_lattice)
#elif DIM == 2
        beta2 = donor(this%r_idx(1),this%r_idx(2))*q**2.0_s/(eps_static(this%MAT)*kB*T_lattice)
#endif
        e_beta = hb**2.0_s*beta2/(2.0_s*me(this%MAT,this%valley))
        costheta = 1.0_s-(2.0_s*(1.0_s-r1)/(1.0_s+4.0_s*r1*Gamma_i/e_beta))

        sb = sqrt(1.0_s-costheta*costheta)
        !call random_number(r1)
        r1 = rng_uniform(rng(this%t_id))
        x_new(1) = k_new_abs*sb*cos(2.0_s*pi*r1)
        x_new(2) = k_new_abs*sb*sin(2.0_s*pi*r1)
        x_new(3) = k_new_abs*costheta

        ! Transform back
        this%k(1) = x_new(1)*cfi0*cth0-x_new(2)*sfi0+x_new(3)*cfi0*sth0
        this%k(2) = x_new(1)*sfi0*cth0+x_new(2)*cfi0+x_new(3)*sfi0*sth0
        this%k(3) = -x_new(1)*sth0+x_new(3)*cth0

        !if (abs(this%k(1)) < 1e-10_s) then
        !    print *, 'val equal to 0', this%k, r1, x_new, sb, costheta, cfi0, cth0, sfi0, sth0
        !    this%k(1) = 1e-10_s
        !end if
        !if (abs(this%k(2)) < 1e-10_s) then
        !    print *, 'val equal to 0', this%k, r1, x_new, sb, costheta, cfi0, cth0, sfi0, sth0
        !    this%k(2) = 1e-10_s
        !end if
        !if (abs(this%k(3)) < 1e-10_s) then
        !    print *, 'val equal to 0', this%k, r1, x_new, sb, costheta, cfi0, cth0, sfi0, sth0
        !    this%k(3) = 1e-10_s
        !end if
        !print *, 'New:', this%energy, sqrt(dot_product(this%k,this%k))
    end subroutine scatter_imp

    subroutine scatter_ac(this, proc)
        class(electron)             :: this
        integer                     :: proc, l_mat, n
        real(kind=s)                :: x1, x2, e_mu, C, l_nonp, gamma_i, ge_mu, pre, r1, r2
        real(kind=s)                :: P_ac, P_max, q_vec, q_e, f
        real(kind=s)                :: k_old_abs, k_new_abs, kxy, k, cth0, sth0, cfi0, sfi0, costheta, sb
        real(kind=s),dimension(3)   :: x_new

        l_nonp = nonp(this%MAT,this%valley)
        l_mat = this%mat
        e_mu = (me(this%MAT,this%valley)*sound_vel(this%MAT)**2.0_s)/2.0_s
        C = 4.0_s*sqrt(e_mu)/(kB*T_lattice*(1.0_s-4.0_s*l_nonp*e_mu))
        ge_mu = e_mu/(1.0_s-4.0_s*l_nonp*e_mu)
        gamma_i = this%energy*(1.0_s+l_nonp*this%energy)
        pre = sqrt(me(l_mat,this%valley))*(kB*T_lattice)**3.0_s*def_pot_ac(l_mat,1)**2.0_s/(2.0_s**(5.0_s/2.0_s)*pi*hb**4.0_s* &
            & sound_vel(l_mat)**4.0_s*rho(l_mat))*(this%energy*(1.0_s+l_nonp*this%energy))**(-1.0_s/2.0_s)
        if (gamma_i .lt. ge_mu) then
            x1 = C*(sqrt(e_mu)*(1.0_s+2.0_s*l_nonp*this%energy-sqrt(gamma_i)))
            x2 = C*(sqrt(e_mu)*(1.0_s+2.0_s*l_nonp*this%energy+sqrt(gamma_i)))
            P_max = pre*(1.0_s+2.0_s*l_nonp*this%energy+2.0_s*l_nonp*kB*T_lattice*x2)*x2**2.0_s*(1.0_s/(exp(x2)-1.0_s))
        else
            x1 = 0.0_s
            ! Absorption
            if (proc .eq. 1) then
                x2 = C*(sqrt(gamma_i)+sqrt(e_mu)*(1.0_s+2.0_s*l_nonp*this%energy))
                P_max = pre*(1.0_s+2.0_s*l_nonp*this%energy+2.0_s*l_nonp*kB*T_lattice*x2)*x2**2.0_s*(1.0_s/(exp(x2)-1.0_s))
            else
            ! Emission
                x2 = C*(sqrt(gamma_i)-sqrt(e_mu)*(1.0_s+2.0_s*l_nonp*this%energy))
                P_max = pre*(1.0_s+2.0_s*l_nonp*this%energy-2.0_s*l_nonp*kB*T_lattice*x2)*x2**2.0_s*((1.0_s/(exp(x2)-1.0_s))+1.0_s)
            end if
        end if

        ! Use rejection technique to find a q-vector (x first)
        n = 1
        do while (n > 0)
            !call random_number(r1)
            r1 = rng_uniform(rng(this%t_id))
            !call random_number(r2)
            r2 = rng_uniform(rng(this%t_id))
            if (gamma_i .lt. ge_mu) then
                f = (1.0_s/(x2-x1))
                r1 = (r1+x1*f)/f
                r2 = P_max*r2
            else
                r1 = r1*x2
                r2 = P_max*r2
            end if
            if (proc .eq. 1) then
                P_ac = pre*(1.0_s+2.0_s*l_nonp*this%energy+2.0_s*l_nonp*kB*T_lattice*r1)*r1**2.0_s*(1.0_s/(exp(r1)-1.0_s))
            else
                P_ac = pre*(1.0_s+2.0_s*l_nonp*this%energy-2.0_s*l_nonp*kB*T_lattice*r1)*r1**2.0_s*((1.0_s/(exp(r1)-1.0_s))+1.0_s)
            end if
            if (r2 <= P_ac) then
                q_vec = kB*T_lattice*r1/(hb*sound_vel(this%mat))
                q_e = hb*q_vec*sound_vel(this%mat)
                if (this%valley .gt. G) then
                    q_vec = kB*T_lattice*r1/(hb*sound_vel(this%mat))
                    q_e = hb*q_vec*sound_vel(this%mat)*sqrt(me(this%mat,this%valley)/m0)
                    !print *, this%energy, q_vec, q_e
                end if
                n = 0
!                print *, q_vec, q_e
            end if
        end do

        ! Length of old k-vector
        if (this%valley .eq. G) then
            k_old_abs = sqrt(2.0_s*me(this%mat,this%valley)*this%energy*q*(1.0_s+l_nonp*this%energy))/hb
        else
        !end if
        !if ((this%valley .eq. X) .or. (this%valley .eq. L)) then
            k_old_abs = sqrt(2.0_s*m0*this%energy*q*(1.0_s+l_nonp*this%energy))/hb
        end if

        ! Calculate new energy
        ! Absorption
        !print *, this%energy
        if (proc .eq. 1) then
            this%energy = this%energy+q_e
        else
        ! Emission
            this%energy = this%energy-q_e
        end if
        if (this%energy <= 0.0_s) then
            print *, 'energy smaller 0', this%energy
            this%energy = 1e-5_s
        end if
        !if (this%valley .eq. L) then
        !    print *, this%energy
        !end if
        ! Compute now the new k-vector
        ! Length of new k-vector
        if (this%valley .eq. G) then
            k_new_abs = sqrt(2.0_s*me(this%mat,this%valley)*this%energy*q*(1.0_s+l_nonp*this%energy))/hb
        else
        !end if
        !if (this%valley .gt. G) then! .or. (this%valley .eq. L)) then
            k_new_abs = sqrt(2.0_s*m0*this%energy*q*(1.0_s+l_nonp*this%energy))/hb
        !    print *, 'knew', k_old_abs, k_new_abs, q_vec
        end if


        ! Rotate coordinate system
        kxy = sqrt(this%k(1)*this%k(1)+this%k(2)*this%k(2))
        k = sqrt(kxy*kxy+this%k(3)*this%k(3))
        cth0 = this%k(3)/k
        sth0 = kxy/k
        cfi0 = this%k(1)/kxy
        sfi0 = this%k(2)/kxy

        ! Determine scattering angle by geometry - law of cosines
        costheta = (k_new_abs**2.0_s+k_old_abs**2.0_s-q_vec**2.0_s)/(2.0_s*k_old_abs*k_new_abs)
        if (costheta >= 1.0_s) then
            print *, 'costheta corrected ac', costheta
            costheta = 0.99999_s
        end if
        if (costheta <= -1.0_s) then
            print *, 'costheta corrected ac', costheta
            costheta = -0.99999_s
        end if
        !print *, 'AC Scattering', proc, this%cidx, this%valley
        !if (abs(costheta) .ge. 1.0_s) then
        !    print *, 'AC scat: Costheta greater than 1!', costheta, this%energy, this%valley, k_new_abs, k_old_abs
        !    if (costheta .lt. -1.0) then
        !        costheta = -0.9999_s
        !    else
        !        costheta = 0.9999_s
        !    end if
        !    print *, this%cidx, costheta
        !end if
        !if (this%valley .eq. L) then
        !    print *, costheta, q_vec, k_new_abs, k_old_abs, this%energy
        !end if

        sb = sqrt(1.0_s-costheta*costheta)
        !call random_number(r1)
        r1 = rng_uniform(rng(this%t_id))
        x_new(1) = k_new_abs*sb*cos(2.0_s*pi*r1)
        x_new(2) = k_new_abs*sb*sin(2.0_s*pi*r1)
        x_new(3) = k_new_abs*costheta

        ! Transform back
        this%k(1) = x_new(1)*cfi0*cth0-x_new(2)*sfi0+x_new(3)*cfi0*sth0
        this%k(2) = x_new(1)*sfi0*cth0+x_new(2)*cfi0+x_new(3)*sfi0*sth0
        this%k(3) = -x_new(1)*sth0+x_new(3)*cth0

        !if (abs(this%k(1)) < 1e-10_s) then
        !    print *, 'ac'
        !    print *, 'val equal to 0', this%k, r1, x_new, sb, costheta, cfi0, cth0, sfi0, sth0
        !    this%k(1) = 1e-10_s
        !end if
        !if (abs(this%k(2)) < 1e-10_s) then
        !    print *, 'ac'
        !    print *, 'val equal to 0', this%k, r1, x_new, sb, costheta, cfi0, cth0, sfi0, sth0
        !    this%k(2) = 1e-10_s
        !end if
        !if (abs(this%k(3)) < 1e-10_s) then
        !    print *, 'ac'
        !    print *, 'val equal to 0', this%k, r1, x_new, sb, costheta, cfi0, cth0, sfi0, sth0
        !    this%k(3) = 1e-10_s
        !end if
        !print *, 'New:', this%valley, this%energy, costheta, sqrt(dot_product(this%k,this%k)), this%k
    end subroutine scatter_ac

!    function calc_Ppop(this,e_old,e_new,Gamma_i,Gamma_n) result(Ppop)
!        class(electron)             :: this
!        real(kind=s)                :: e_old, e_new, Gamma_i, Gamma_n
!        real(kind=s),dimension(100) :: beta, Ppop
!        integer                     :: i
!
!        do i=1,100
!            beta(i) = pi/100.0_s*i
!        end do
!        Ppop = (sqrt(Gamma_i)*sqrt(Gamma_n)+nonp(this%MAT,this%valley)*e_old*e_new*cos(beta))**2.0_s/ &
!             & (Gamma_i+Gamma_n-2.0_s*sqrt(Gamma_i)*sqrt(Gamma_n)*cos(beta))*sin(beta)
!    end function calc_Ppop

    ! Initialize all carriers according to Maxwell-Boltzmann statistics - thermal equilibrium
    subroutine init(this, MAT)
         class(electron)            :: this
         integer                    :: MAT
         real(kind=s)               :: r,r2,k_abs,ct, l_nonp, l_me

         this%mat = MAT
         ! located in Gamma-valley
         this%valley = 1
         l_nonp = nonp(this%mat,this%valley)
         l_me = me(this%mat,this%valley)
         ! initial energy follows Maxwell-Boltzmann distribution
           !call random_number(r)
         !this%t_id = 1
         r = rng_uniform(rng(this%t_id))
         this%energy = -1.5_s*kB*t_lattice*log(r)!*4.0_s
         !this%energy = abs(log(1.0_s/r-1.0_s)*kB*T_lattice+0.15_s)
         call get_e_idx(this)

         ! select k
         k_abs = sqrt(2.0_s*me(MAT,1)*this%energy*(1.0_s+nonp(MAT,1)*this%energy))/hb
         !call random_number(r)
         r = rng_uniform(rng(this%t_id))
         !call random_number(r2)
         r2 = rng_uniform(rng(this%t_id))
         ct = (1.0_s-2.0_s*r)
         this%k(1) = k_abs*sqrt(1.0_s-ct*ct)*cos(2.0_s*pi*r2)
         this%k(2) = k_abs*sqrt(1.0_s-ct*ct)*sin(2.0_s*pi*r2)
         this%k(3) = k_abs*ct

         ! Give velocity
         this%vel = hb*this%k/(l_me*(1.0_s+2.0_s*l_nonp*this%energy))
    end subroutine init

#if DIM == 1
    subroutine init_1(this)
        class(electron)             :: this
        integer                     :: MAT
        integer,dimension(1)        :: ridx
        real(kind=s)                :: r, r1, rtmp

        call random_number(r)
        rtmp = r
        if (r == 1.0_s) then
            r = 0.999999999_s
        end if
        if (r == 0.0_s) then
            r = 0.000000001_s
        end if
        !r = rng_uniform(rng(this%t_id))
        ridx = minloc(abs(r-prob_init), MASK = ((r-prob_init .lt. 0) .and. (prob_init > 0.0_s)))
        if (r > prob_init(sum(num_nodes))) then
            ridx = sum(num_nodes)
        end if
                !print *, r, ridx(1)
        call random_number(r)
        !r = rng_uniform(rng(this%t_id))
        r1 = (r-0.5_s)*m_lx(ridx(1))+node_coord(ridx(1))
        !print *, r1
        if (r1 .lt. 0.0_s) then
            r1 = r1*(-1.0_s)
        end if
        if (r1 .gt. dev_x) then
            r1 = dev_x-(r1-dev_x)
        end if
        this%r = (/r1,0.0_s,0.0_s/)
        call get_r_idx(this%r,this%r_idx)
        !if (this%r_idx == 1) then
        !end if
        MAT = m_material(m_region(this%r_idx))
        !print *, MAT, this%r(1)
        call init(this,MAT)
        !print *, 'init', this%r(1), ridx, rtmp, this%energy, this%mat
    end subroutine init_1
#endif

#if DIM == 2
    subroutine init_2(this)
        class(electron)             :: this
        integer                     :: MAT
        integer,dimension(2)        :: ridx
        real(kind=s)                :: r, r1, r2

        call random_number(r)
        if (r == 1.0_s) then
            r = 0.99999999_s
        end if
        if (r == 0.0_s) then
            r = 0.000000001_s
        end if
        !ridx = maxloc(prob_init-r, MASK = (prob_init-r .lt. 0))
        ridx = minloc(abs(r-prob_init), MASK = ((r-prob_init .lt. 0) .and. (prob_init > 0.0_s)))
        !print *, r, ridx

        call random_number(r)
        r1 = (r-0.5_s)*m_lx(ridx(1))+node_coord(1,ridx(1))
        call random_number(r)
        r2 = (r-0.5_s)*m_ly(ridx(2))+node_coord(2,ridx(2))
        if (r1 .lt. 0.0_s) then
            r1 = r1*(-1.0_s)
        end if
        if (r1 .gt. dev_x) then
            r1 = dev_x-(r1-dev_x)
        end if
        if (r2 .lt. 0.0_s) then
            r2 = r2*(-1.0_s)
        end if
        if (r2 .gt. dev_y) then
            r2 = dev_y-(r2-dev_y)
        end if
        this%r = (/r1,r2,0.0_s/)
        call get_r_idx_2(this%r,this%r_idx)
        !print *, this%r_idx, this%r(1:2),m_region(this%r_idx(1),this%r_idx(2))
        MAT = m_material(m_region(this%r_idx(1),this%r_idx(2)))
        !print *, MAT
        call init(this,MAT)
        if (this%r(1) < 0) then
            print *, ridx, r
        end if
        if (this%r(2)  < 0) then
            print *, ridx, r
        end if
    end subroutine init_2
#endif

#if DIM == 3
    subroutine init_3(this)
        class(electron)             :: this
        integer                     :: MAT
        integer,dimension(3)        :: ridx
        real(kind=s)                :: r1, r2, r3, r

        r = rng_uniform(rng(this%t_id))
        ridx = maxloc(prob_init-r, MASK = (prob_init-r .lt. 0))
        !call random_number(r1)
        r = rng_uniform(rng(this%t_id))
        r1 = (r-0.5_s)*m_lx(ridx(1))+node_coord(1,ridx(1))
        !call random_number(r2)
        r = rng_uniform(rng(this%t_id))
        r2 = (r-0.5_s)*m_ly(ridx(2))+node_coord(2,ridx(2))
        !call random_number(r3)
        r = rng_uniform(rng(this%t_id))
        r3 = (r-0.5_s)*m_lz(ridx(3))+node_coord(3,ridx(3))
        !print *, r1, r2, r3
        if (r1 .lt. 0.0_s) then
            r1 = r1*(-1.0_s)
        end if
        if (r1 .gt. dev_x) then
            r1 = dev_x-(r1-dev_x)
        end if
        if (r2 .lt. 0.0_s) then
            r2 = r2*(-1.0_s)
        end if
        if (r2 .gt. dev_y) then
            r2 = dev_y-(r2-dev_y)
        end if
        if (r3 .lt. 0.0_s) then
            r3 = r3*(-1.0_s)
        end if
        if (r3 .gt. dev_z) then
            r3 = dev_z-(r3-dev_z)
        end if
        this%r = (/r1,r2,r3/)
        call get_r_idx_3(this%r,this%r_idx)
        !print *, this%r_idx, this%r(1:2),m_region(this%r_idx(1),this%r_idx(2))
        MAT = m_material(m_region(this%r_idx(1),this%r_idx(2),this%r_idx(3)))
        !print *, MAT
        call init(this,MAT)
    end subroutine init_3
#endif

    subroutine init_X(this, MAT)
         class(electron)            :: this
         integer                    :: MAT
         real(kind=s)               :: r,r2,k_abs,ct

         this%mat = MAT
         ! located in Gamma-valley
         this%valley = 2
         this%v_subidx = int(rand(0)*8)+1

         ! initial energy follows Maxwell-Boltzmann distribution
         !call random_number(r)
         r = rng_uniform(rng(this%t_id))
         this%energy = -1.5_s*kB*t_lattice*log(r)
         call get_e_idx(this)

         ! select k
         k_abs = sqrt(2.0_s*m0*this%energy*(1.0_s+nonp(MAT,1)*this%energy))/hb
         !call random_number(r)
         r = rng_uniform(rng(this%t_id))
         !call random_number(r2)
         r2 = rng_uniform(rng(this%t_id))
         ct = (1.0_s-2.0_s*r)
         this%k(1) = k_abs*sqrt(1.0_s-ct*ct)*cos(2.0_s*pi*r2)
         this%k(2) = k_abs*sqrt(1.0_s-ct*ct)*sin(2.0_s*pi*r2)
         this%k(3) = k_abs*ct
    end subroutine init_X

    subroutine time_2e(this,t_energy,n_idx)
        class(electron)                        :: this
        real(kind=s)                           :: l_nonp, l_me
        real(kind=s)                           :: E_lower, E_upper, r2_lower, r2_upper
        real(kind=s)                           :: dot_kE, dot_Efield, t_lower_1, t_lower_2, t_upper_1, t_upper_2, p
        real(kind=s)                           :: dot_k, sqrt_lower, sqrt_upper, sqrt_tmp
        real(kind=s),dimension(4)              :: t_arr
        integer,dimension(1),intent(out)       :: n_idx
        real(kind=s),intent(out)               :: t_energy
        real(kind=s),dimension(3)              :: e_field_tmp, e_field

        !e_field = el_field(this%r_idx)
        e_field = -this%loc_e_field
        l_nonp = nonp(this%mat,this%valley)
        l_me = me(this%mat,this%valley)
        ! Calculate the upper and lower energy boundary
        !e_idx = int((this%energy)/bin_dE)+1
        if (this%e_idx > 1) then
            E_lower = (this%e_idx-1)*bin_dE
            E_upper = this%e_idx*bin_dE
        else
            E_lower = 0.00001_s !1.0e-6_s
            E_upper = bin_dE
        end if
        !print *, this%energy, this%e_idx
        e_field_tmp = e_field
        ! Calculate the corresponding reciprocal vector length (already squared, see notes)
        r2_lower = (E_lower+l_nonp*E_lower**2.0_s)*2.0_s*l_me/hb**2.0_s
        r2_upper = (E_upper+l_nonp*E_upper**2.0_s)*2.0_s*l_me/hb**2.0_s
        if (this%valley == L) then
            r2_lower = (E_lower+l_nonp*E_lower**2.0_s)*2.0_s*m0/hb**2.0_s
            r2_upper = (E_upper+l_nonp*E_upper**2.0_s)*2.0_s*m0/hb**2.0_s
            if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
                e_field_tmp = matmul((T_hv_L(this%MAT,1,:,:)),e_field)
            else if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
                e_field_tmp = matmul((T_hv_L(this%MAT,2,:,:)),e_field)
            else if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
                e_field_tmp = matmul((T_hv_L(this%MAT,3,:,:)),e_field)
            else if ((this%v_subidx == 7) .or. (this%v_subidx == 8)) then
                e_field_tmp = matmul((T_hv_L(this%MAT,4,:,:)),e_field)
            end if
            !print *, e_field_tmp, e_field
        end if
        if (this%valley == X) then
            r2_lower = (E_lower+l_nonp*E_lower**2.0_s)*2.0_s*m0/hb**2.0_s
            r2_upper = (E_upper+l_nonp*E_upper**2.0_s)*2.0_s*m0/hb**2.0_s
            if ((this%v_subidx == 1) .or. (this%v_subidx == 2)) then
                e_field_tmp = matmul(T_hv(this%MAT,X,1,:,:),e_field)
            else if ((this%v_subidx == 3) .or. (this%v_subidx == 4)) then
                e_field_tmp = matmul(T_hv(this%MAT,X,2,:,:),e_field)
            else if ((this%v_subidx == 5) .or. (this%v_subidx == 6)) then
                e_field_tmp = matmul(T_hv(this%MAT,X,3,:,:),e_field)
            end if
            !print *, e_field_tmp, e_field
        end if

        ! Do some precalculations, only the length of the reciprocal space vector is changing
        dot_k = dot_product(this%k,this%k)
        dot_kE = dot_product(this%k,e_field_tmp)
        dot_Efield = dot_product(e_field_tmp,e_field_tmp)

        p = -dot_kE*hb/dot_Efield
        sqrt_tmp = p**2.0_s-(hb**2.0_s/dot_Efield*(dot_k-r2_lower))
        if (sqrt_tmp >= 0) then
            sqrt_lower = sqrt(sqrt_tmp) !sqrt(p**2.0_s-(hb**2.0_s/dot_Efield*(dot_k-r2_lower)))
        else
            sqrt_lower = -1e10_s
        end if
        sqrt_tmp = p**2.0_s-(hb**2.0_s/dot_Efield*(dot_k-r2_upper))
        if (sqrt_tmp >= 0) then
        sqrt_upper = sqrt(sqrt_tmp) !sqrt(p**2.0_s-(hb**2.0_s/dot_Efield*(dot_k-r2_upper)))
        else
            sqrt_upper = -1e10_s
        end if
        t_upper_1 = p+sqrt_upper
        t_upper_2 = p-sqrt_upper
        t_lower_1 = p+sqrt_lower
        t_lower_2 = p-sqrt_lower

        !t_energy = 0
        t_arr = (/t_upper_1, t_upper_2, t_lower_1, t_lower_2/)
        t_energy =  minval(t_arr, MASK = (t_arr .gt. 0) .and. (t_arr .gt. 1e-11_s))
        n_idx =  minloc(t_arr, MASK = (t_arr .gt. 0) .and. (t_arr .gt. 1e-11_s))
            !print *, this%energy, this%e_idx, this%k, t_energy
            !print *, n_idx, t_arr
        !if ((n_idx(1) == 1) .or. (n_idx(1) == 2)) then
        !    this%e_idx = this%e_idx+1
        !else if ((n_idx(1) == 3) .or. (n_idx(1) == 4)) then
        !    this%e_idx = this%e_idx-1
        !else if (n_idx(1) == 0) then
        !    this%e_idx = 1
        !    t_energy = 0
        !end if
        !if (this%valley == X) then
            !print *, E_lower, E_upper, this%energy, t_arr
        !end if
        !print *, t_arr, n_idx, t_energy, r2_lower, p**2, (hb**2.0_s/dot_Efield*(dot_k-r2_lower))
        !print *, E_lower, E_upper, t_energy, t_arr!, r2_lower, r2_upper
            !print *, this%energy, this%e_idx
    end subroutine time_2e

    subroutine get_e_idx(this)
        class(electron)         :: this
        this%e_idx = int(this%energy/bin_dE)+1
    end subroutine get_e_idx

!    function check_ppop() result(costheta)
!        integer                     :: scat_idx, n,j, mat, valley
!        real(kind=s)                :: r1, r2, k_new_abs, kxy, k, cth0, sth0, cfi0, sfi0, zeta, costheta, sb, e_new, e_old
!        real(kind=s)                :: Gamma_i, Gamma_n, beta, Ppop, Ppop_max
!        real(kind=s),dimension(3)   :: x_new
!        mat = 1
!        valley = 1
!        n = 1
!        e_old = 0.15
!        e_new = e_old-phonon_e(1,1)
!        ! Determine scattering angle - using rejection method
!        Gamma_i = e_old*(1.0_s+nonp(mat,valley)*e_old)
!        Gamma_n = e_new*(1.0_s+nonp(mat,valley)*e_new)
!
!        ! Following the calculation of Fawcett
!        Ppop_max = (sqrt(Gamma_i)*sqrt(Gamma_n)+nonp(MAT,valley)*e_old*e_new*1)**2.0_s/ &
!                 & (Gamma_i+Gamma_n-2.0_s*sqrt(Gamma_i)*sqrt(Gamma_n)*1)
!        do while (n > 0)
!            call random_number(r1)
!            r1 = r1*2.0_s-1.0_s
!            call random_number(r2)
!            r2 = r2*Ppop_max
!            beta = acos(r1)
!            Ppop = (sqrt(Gamma_i)*sqrt(Gamma_n)+nonp(MAT,valley)*e_old*e_new*cos(beta))**2.0_s/ &
!                 & (Gamma_i+Gamma_n-2.0_s*sqrt(Gamma_i)*sqrt(Gamma_n)*cos(beta))
!            if (r2 <= Ppop) then
!                n = -1
!                costheta = r1
!            end if
!        end do
!    end function check_ppop

    subroutine handle_heterostructure(this,border_idx,mat1,mat2)
        class(electron)                 :: this
        integer                         :: border_idx, mat1, mat2, visit
        real(kind=s)                    :: k_abs, kt, energy_new, delta_e

        this%visit = 0
        visit = -1
   !     if (this%cidx == 7508) then
        !if (this%valley > G) then
 !   print *, 'before heterjunction', this%cidx, this%energy, this%mat, mat1, mat2, this%r(1:2), this%r_old(1:2), this%vel(1:2)!, &
    !        & this%valley
   !     end if
   !print *, this%cidx, border_idx,mat1, mat2
        if ((border_idx == 1) .or. (border_idx == 2)) then
            kt = sqrt(this%k(2)**2.0_s+this%k(3)**2.0_s)
        end if
        if ((border_idx == 3) .or. (border_idx == 4)) then
            kt = sqrt(this%k(1)**2.0_s+this%k(3)**2.0_s)
        end if
        delta_e = emin(mat1)-emin(mat2)
        if (this%valley > G) then
            delta_e = (emin(this%mat)+valley_offs(this%mat,this%valley))-(emin(mat2)+valley_offs(mat2,this%valley))
        end if
        energy_new = this%energy+delta_e
        if (energy_new > 0.0_s) then
            if (this%valley .eq. G) then
                k_abs = sqrt(2.0_s*me(mat2,1)*energy_new*(1.0_s+nonp(mat2,1)*energy_new))/hb
            else
                k_abs = sqrt(2.0_s*m0*energy_new*(1.0_s+nonp(mat2,this%valley)*energy_new))/hb
            end if
        else
            k_abs = 0.0_s
        end if
            if (delta_e < 0) then
                !if ((this%energy+emin(mat1)+valley_offs(this%mat,this%valley) >= &
                !  & (emin(mat2)+valley_offs(mat2,this%valley))) .and. (k_abs > kt)) then
                if ((energy_new > 0) .and. (k_abs > kt)) then
                    this%visit = 1
                    ! subtract energy and conserve parallel momentum
                    this%energy = energy_new
                    this%mat = mat2

                    call cross_heterojunction(this, border_idx)

                    !if (this%cidx == 10000) then
                    !if (this%valley > G) then
   !                 print *, 'well -> barrier', this%cidx, this%energy, this%r(1:2)!, &
          !                  & this%k, this%valley, this%r_old(1:2), this%vel(1:2)!, k_abs, kt
                    !end if
                else
                    this%visit = 2
                    ! reflect
                  ! if (this%cidx == 7508) then
   !                 print *, 'particle reflection at barrier', this%cidx, this%energy, this%r(1:2)!, &
                  !         & this%k, this%valley
                   ! end if
                    !if (((this%vel(1) > 0.0_s) .and. (border_idx == 1)) .or. ((this%vel(1) < 0.0_s) .and. (border_idx == 2))) then
                    if ((border_idx == 1) .or. (border_idx == 2)) then
                        if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                            this%k(1) = this%k(1)*(-1.0_s)
                        end if
                        if (this%valley == L) then
                            call border_L_yz(this)
                        end if
                    end if
                    if (((this%vel(2) > 0.0_s) .and. (border_idx == 3)) .or. ((this%vel(2) < 0.0_s) .and. (border_idx == 4))) then
                    !if ((border_idx == 3) .or. (border_idx == 4)) then
                        if ((this%valley .eq. G) .or. (this%valley .eq. X)) then
                            this%k(2) = this%k(2)*(-1.0_s)
                            visit = -2
                        end if
                        if (this%valley == L) then
                            call border_L_xz(this)
                            visit = -3
                        end if
                    end if
                    if (border_idx == 1) then
                        this%r(1) = this%r(1)-1e-8_s
                    elseif (border_idx == 2) then
                        this%r(1) = this%r(1)+1e-8_s
                    elseif (border_idx == 3) then
                        this%r(2) = this%r(2)-1e-8_s
                    elseif (border_idx == 4) then
                        this%r(2) = this%r(2)+1e-8_s
                    end if
                    call compute_vel(this)
                    !if ((border_idx == 4) .and. (this%vel(2) > 0.0_s) .and. (this%r(2) > 55.0_s)) then
                    if ((border_idx == 3) .and. (this%vel(2) > 0.0_s) .and. (this%r(2) > 55.0_s)) then
                        print *, this%cidx, this%old_k(1:2), this%old_k2(1:2)
                        print *, this%cidx, this%k(1:2), this%vel(1:2)
                        print *, 'error in reflection part', this%cidx, this%k(1:2), this%mat, &
                     & this%valley, this%vel(1:2), this%r(1:2), visit, this%r_old(1:2),  this%r_old2(1:2)
                    end if
                end if
            else
                this%visit = 3
                ! add energy and conserve momentum
                !if (this%cidx == 10000) then
                !if (this%valley > G) then
                    !print *, 'barrier -> well', this%cidx,this%energy, this%r(1), this%r_old(1:2), this%mat, this%vel(1:2)
                !end if
                this%energy = energy_new !this%energy+delta_e !emin(mat1)-valley_offs(this%mat,this%valley)
                this%mat = mat2
                call cross_heterojunction(this, border_idx)

                !if (this%cidx == 10000) then
                !if (this%valley > G) then
    !                print *, 'barrier -> well'!, this%cidx, this%energy, this%r(1:2), k_abs, &
         !                  & this%k, this%valley!,this%r_old(1:2)
                !end if
            end if
        !end if
        call compute_vel(this)
        if (this%visit == 0) then
            print *, 'carrier not handled!'
        end if
        !if (this%cidx == 7508) then
        !if (this%valley > G) then
     !   print *, 'after heterjunction', this%cidx, this%energy, this%mat, mat1, mat2, this%r(1:2),this%r_old(1:2), this%vel(1:2)
        !end if
        !print *, this%vel(1)
    end subroutine handle_heterostructure

    subroutine cross_heterojunction(this,border_idx)
        class(electron)             :: this
        integer                     :: border_idx
        real(kind=s)                :: k_abs

        if (this%valley .eq. G) then
            k_abs = sqrt(2.0_s*me(this%mat,1)*this%energy*(1.0_s+nonp(this%mat,1)*this%energy))/hb
        else
            k_abs = sqrt(2.0_s*m0*this%energy*(1.0_s+nonp(this%mat,this%valley)*this%energy))/hb
        end if
        if (border_idx == 1) then
            this%k(1) = sqrt(k_abs**2.0_s-this%k(2)**2.0_s-this%k(3)**2.0_s)
        elseif (border_idx == 2) then
            this%k(1) = -sqrt(k_abs**2.0_s-this%k(2)**2.0_s-this%k(3)**2.0_s)
        end if
        if (border_idx == 3) then
            this%k(2) = sqrt(k_abs**2.0_s-this%k(1)**2.0_s-this%k(3)**2.0_s)
        elseif (border_idx == 4) then
            this%k(2) = -sqrt(k_abs**2.0_s-this%k(1)**2.0_s-this%k(3)**2.0_s)
        end if
        call compute_vel(this)
        !print *, this%cidx, this%k, this%vel, this%valley
        if (this%valley == L) then
            if ((border_idx == 1) .and. (this%vel(1) < 0.0_s)) then
                call border_L_yz(this)
            end if
            if ((border_idx == 2) .and. (this%vel(1) > 0.0_s)) then
                call border_L_yz(this)
            end if
        end if
        if (this%valley == L) then
            if ((border_idx == 3) .and. (this%vel(2) < 0.0_s)) then
                call border_L_xz(this)
            end if
            if ((border_idx == 4) .and. (this%vel(2) > 0.0_s)) then
                call border_L_xz(this)
            end if
        end if
        if (border_idx == 1) then
            this%r(1) = this%r(1)+1e-8_s
        elseif (border_idx == 2) then
            this%r(1) = this%r(1)-1e-8_s
        elseif (border_idx == 3) then
            this%r(2) = this%r(2)+1e-8_s
        elseif (border_idx == 4) then
            this%r(2) = this%r(2)-1e-8_s
        end if
    end subroutine cross_heterojunction

#if DIM == 1
    subroutine inject_carrier()
        integer                             :: i, contact_id, cell_id
        real(kind=s)                        :: r, t_left

        do contact_id=1,2
            if (contact_type(contact_id) .eq. OHMIC) then
            i = 0
            if (contact_id .eq. 1) then
                cell_id = 1
            else if (contact_id .eq. 2) then
                cell_id = sum(num_nodes)
            end if
            do while ((donor(cell_id)-contact_conc(contact_id)) .gt. (donor(cell_id)*0.001))
                    i = i+1
                    if (c(i)%valley .eq. 0) then
                        c(i)%t_id = 1
                        c(i)%cidx = i
                        call init(c(i),m_material(m_region(cell_id))) !c(i)%
                        call random_number(r)
                        c(i)%vel(1) = sqrt((-1.0_s*log(r)*2.0_s*kB*T_lattice)/(me(m_material(m_region(cell_id)),1)))
                        c(i)%k(1) = c(i)%vel(1)*(me(m_material(m_region(cell_id)),1)*(1.0_s+2.0_s* &
                                  & nonp(m_material(m_region(cell_id)),1)*c(i)%energy))/hb
                        if (contact_id .eq. 1) then
                            c(i)%r(1) = 0.0_s
                        else if (contact_id .eq. 2) then
                            c(i)%r(1) = dev_x
                            c(i)%k(1) = c(i)%k(1)*(-1.0_s)
                            c(i)%vel(1) = c(i)%vel(1)*(-1.0_s)
                        end if
                        call random_number(r)
                        t_left = r*t_step
                        if (t_left .gt. t_rest(contact_id)) then
                            t_left = t_rest(contact_id)
                        end if
                        c(i)%e_time = step*t_step-t_left
                        do while ((c(i)%synchronized < 0.1) .and. (c(i)%valley > 0) .and. (c(i)%op == 0))
                            call c(i)%drift(1)
                        end do
                        c(i)%synchronized = 0
                        call keep_ohmic_contacts()
                    end if
            end do
            end if
        end do
    end subroutine inject_carrier

    subroutine assign_force_1(this)
        class(electron)                 :: this
        integer                         :: i
        real(kind=s)                    :: x

            if ((this%valley > 0) .and. (this%op == 0)) then
                x = this%r(1)/m_dx(m_region(this%r_idx))
                i = int(x+1)
                this%loc_e_field(1) = el_field(i)*(1-(x-(i-1)))+ &
                        & el_field(i+1)*(x-(i-1))
                this%loc_e_field(1) = this%loc_e_field(1)
            end if
    end subroutine assign_force_1
#endif


#if DIM == 2
    subroutine assign_force_2(this)
        class(electron)                 :: this
        integer                         :: i, j
        real(kind=s)                    :: x, y

            if ((this%valley > 0) .and. (this%op == 0)) then
!                x = this%r(1)/m_dx
!                y = this%r(2)/m_dy
!                i = int(x+1)
!                j = int(y+1)
                i = int((this%r(1)+m_dx_half)/m_dx)+1
                j = int((this%r(2)+m_dy_half)/m_dy)+1
if (j < 0) then
    print *, this%cidx, this%r, this%r_old2, this%r_idx, this%energy, this%k, this%op
    this%valley = 0
else
!                this%loc_e_field(1) = el_field(1,i,j)*(1-(x-(i-1)))*(1-(y-(j-1)))+ &
!                        & el_field(1,i,j+1)*(1-(x-(i-1)))*(y-(j-1))+ &
!                        & el_field(1,i+1,j+1)*(x-(i-1))*(y-(j-1))+ &
!                        & el_field(1,i+1,j)*(x-(i-1))*(1-(y-(j-1)))
!                this%loc_e_field(2) = el_field(2,i,j)*(1-(x-(i-1)))*(1-(y-(j-1)))+ &
!                        & el_field(2,i,j+1)*(1-(x-(i-1)))*(y-(j-1))+ &
!                        & el_field(2,i+1,j+1)*(x-(i-1))*(y-(j-1))+ &
!                        & el_field(2,i+1,j)*(x-(i-1))*(1-(y-(j-1)))
                this%loc_e_field(1) = el_field(1,i,j)
                this%loc_e_field(2) = el_field(2,i,j)
    end if
            end if
    end subroutine assign_force_2
#endif

#if DIM == 3
    subroutine assign_force_3(this)
        class(electron)                 :: this
        integer                         :: k, i, j, m
        real(kind=s)                    :: x, y, z

                !if (this%cidx == 8081) then
                !print *, step, this%r
                !print *, x,y,z, i, j, k
                !print *, dev%g_el_field(step,1,i,j,k)
                !end if
            if ((this%valley > 0) .and. (this%op == 0)) then
                x = this%r(1)/m_dx(1)
                y = this%r(2)/m_dx(1)
                z = this%r(3)/m_dx(1)
                i = int(x+1)
                j = int(y+1)
                k = int(z+1)

                !if (this%cidx == 8081) then
                !print *, step, this%r
                !print *, x,y,z, i, j, k
                !print *, dev%g_el_field(step,1,i,j,k)
                !end if
                if ((j > 101) .or. (j<1)) then
                !    print *, 'assfdf'
                print *, this%cidx, this%loc_e_field, this%k, this%energy
                end if
                this%loc_e_field(1) = el_field(1,i,j,k)*(1-(x-(i-1)))*(1-(y-(j-1)))*(1-(z-(k-1)))+ &
                        & el_field(1,i,j+1,k)*(1-(x-(i-1)))*(y-(j-1))*(1-(z-(k-1)))+ &
                        & el_field(1,i,j+1,k+1)*(1-(x-(i-1)))*(y-(j-1))*(z-(k-1))+ &
                        & el_field(1,i,j,k+1)*(1-(x-(i-1)))*(1-(y-(j-1)))*(z-(k-1))+ &
                        & el_field(1,i+1,j+1,k+1)*(x-(i-1))*(y-(j-1))*(z-(k-1))+ &
                        & el_field(1,i+1,j,k)*(x-(i-1))*(1-(y-(j-1)))*(1-(z-(k-1)))+ &
                        & el_field(1,i+1,j+1,k)*(x-(i-1))*(y-(j-1))*(1-(z-(k-1)))+ &
                        & el_field(1,i+1,j,k+1)*(x-(i-1))*(1-(y-(j-1)))*(z-(k-1))
                this%loc_e_field(2) = el_field(2,i,j,k)*(1-(x-(i-1)))*(1-(y-(j-1)))*(1-(z-(k-1)))+ &
                        & el_field(2,i,j+1,k)*(1-(x-(i-1)))*(y-(j-1))*(1-(z-(k-1)))+ &
                        & el_field(2,i,j+1,k+1)*(1-(x-(i-1)))*(y-(j-1))*(z-(k-1))+ &
                        & el_field(2,i,j,k+1)*(1-(x-(i-1)))*(1-(y-(j-1)))*(z-(k-1))+ &
                        & el_field(2,i+1,j+1,k+1)*(x-(i-1))*(y-(j-1))*(z-(k-1))+ &
                        & el_field(2,i+1,j,k)*(x-(i-1))*(1-(y-(j-1)))*(1-(z-(k-1)))+ &
                        & el_field(2,i+1,j+1,k)*(x-(i-1))*(y-(j-1))*(1-(z-(k-1)))+ &
                        & el_field(2,i+1,j,k+1)*(x-(i-1))*(1-(y-(j-1)))*(z-(k-1))
                this%loc_e_field(3) = el_field(3,i,j,k)*(1-(x-(i-1)))*(1-(y-(j-1)))*(1-(z-(k-1)))+ &
                        & el_field(3,i,j+1,k)*(1-(x-(i-1)))*(y-(j-1))*(1-(z-(k-1)))+ &
                        & el_field(3,i,j+1,k+1)*(1-(x-(i-1)))*(y-(j-1))*(z-(k-1))+ &
                        & el_field(3,i,j,k+1)*(1-(x-(i-1)))*(1-(y-(j-1)))*(z-(k-1))+ &
                        & el_field(3,i+1,j+1,k+1)*(x-(i-1))*(y-(j-1))*(z-(k-1))+ &
                        & el_field(3,i+1,j,k)*(x-(i-1))*(1-(y-(j-1)))*(1-(z-(k-1)))+ &
                        & el_field(3,i+1,j+1,k)*(x-(i-1))*(y-(j-1))*(1-(z-(k-1)))+ &
                        & el_field(3,i+1,j,k+1)*(x-(i-1))*(1-(y-(j-1)))*(z-(k-1))

            end if
            !print *, this%cidx
    end subroutine assign_force_3
#endif
#if DIM == 2
    function par_tunnel(this, direction, mat2) result(did_tunnel)
        type(electron)  :: this
        type(tunnel)    :: tun_obj
        integer         :: did_tunnel, direction, i, mat2
        integer,dimension(1) :: e_idx
        real(kind=s)    :: r, delta_e, kt_old, kt_new, energy_new, energy_old

        did_tunnel = 0
        energy_old = this%energy
        delta_e = (emin(this%mat)+valley_offs(this%mat,this%valley))-(emin(mat2)+valley_offs(mat2,this%valley))
        !print *, delta_e
        !print *, this%r(1:2), this%r_idx, this%vel(2)
        do i=n_tunnel,1,-1
        if ((this%r_idx(1) >= tun(i)%t_x_idx(1)) .and. (this%r_idx(1) <= tun(i)%t_x_idx(2)) .and. &
         & ((this%r_idx(2) == tun(i)%t_y_idx(1)) .or. (this%r_idx(2) == tun(i)%t_y_idx(2))) .and. &
         & (tun(i)%direction == direction) .and. (this%energy < abs(delta_e)) .and. (step*t_step > 0.051_s)) then
            if ((direction == 2) .and. (this%vel(2) < 0.0_s)) then
                did_tunnel = i
                exit
                !print *, i, this%cidx, this%r, this%r_idx(1), this%r_idx(2), this%energy
            end if
            if ((direction == 1) .and. (this%vel(2) > 0.0_s)) then
                did_tunnel = i
                exit
                !print *, 'up', i, this%cidx, this%r, this%r_idx(1), this%r_idx(2), this%energy
            end if
        end if
        end do
            !if (this%valley > 1) then
            !    print *, 'carrier in L/X valley', this%valley, this%energy, this%r(1:2)
            !!    did_tunnel = 0
            !end if
        if (did_tunnel > 0) then
            tun_obj = tun(did_tunnel)
            !e_idx = minloc(abs(tun_obj%energy(:,this%valley)-this%energy))
            !print *, e_idx
            e_idx = int(this%energy/(tun_obj%energy(2,this%valley)-tun_obj%energy(1,this%valley))+0.5_s)+1
            !print *, e_idx
            if (this%energy > tun_obj%energy(1000,this%valley)) then
                e_idx(1) = 1000
            end if
            !energy_new = this%energy-tun_obj%energy(tun_obj%pot_min,this%valley) !(e_idx(1))
            energy_new = this%energy-tun_obj%epot_interp(tun_obj%pot_min,this%valley)
            !print *, e_idx, tun_obj%epot_interp(tun_obj%pot_min,this%valley), this%valley, did_tunnel
            if (energy_new < 0.0_s) then
                energy_new = 1e-5_s
            end if
            kt_old = sqrt(this%k(1)**2.0_s+this%k(3)**2.0_s)
            if (this%valley .eq. G) then
                kt_new = sqrt(2.0_s*me(mat2,1)*energy_new*(1.0_s+nonp(mat2,1)*energy_new))/hb
            else
                kt_new = sqrt(2.0_s*m0*energy_new*(1.0_s+nonp(mat2,this%valley)*energy_new))/hb
            end if

            !print *, 'tunneling', this%cidx, did_tunnel, this%energy, energy_new, kt_new, kt_old
            !print *, 'tunnel prob', tun_obj%tunnel_prob(e_idx(1),this%valley)
            r = rng_uniform(rng(this%t_id))
            if ((kt_new > kt_old) .and. (r < tun_obj%tunnel_prob(e_idx(1),this%valley))) then
!if (did_tunnel == 5) then
!            print *, this%r(1:2), this%vel(1:2), this%mat, this%valley
!end if
                    this%r(2) = tun_obj%x_interp(tun_obj%pot_min)!e_idx(1)
                    call get_r_idx_2(this%r,this%r_idx)
                    this%mat = m_material(m_region(tun_obj%t_x_idx(1)+1,tun_obj%t_y_idx(1)+1))
                    !if (tun(did_tunnel)%t_dest(tun(did_tunnel)%pot_min) <= tun(did_tunnel)%t_y(1)) then
                    !    this%mat = m_material(m_region(tun(did_tunnel)%t_x_idx(1)+1,tun(did_tunnel)%t_y_idx(1)-1))
                    !end if
                    this%energy = this%energy-tun_obj%epot_interp(tun_obj%pot_min,this%valley) !e_idx(1))
                    if (this%energy < 0.0_s) then
                        this%energy = 1e-5_s
                    end if
                    if (tun_obj%direction == 1) then
                        call cross_heterojunction(this,3)
                    end if
                    if (tun_obj%direction == 2) then
                        call cross_heterojunction(this,4)
                    end if
!if (did_tunnel == 5) then
!    print *, 'carrier tunneled', energy_old, this%energy, this%mat, this%r(1:2), this%k(1:3), this%vel, &
!    & tun_obj%epot_interp(tun_obj%pot_min,this%valley)
!end if
            else
                did_tunnel = 0
            end if
        end if
    end function par_tunnel
#endif
end module mc_core
