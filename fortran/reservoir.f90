module reservoir
    use config
    use mc_core
    use ohmic_contacts
    implicit none

    integer                 :: n_res
    real(kind=s)                     :: res_len
    integer             :: c_ent, c_left


#if DIM == 1
    type                    :: res
    real(kind=s)            :: doping, dop_charge
    real(kind=s),dimension(:),allocatable   :: conc, charge
    integer,dimension(2)    :: x
    integer                 :: loc, res_idx, max_c, del_carriers, add_carriers
    integer,dimension(:,:),allocatable    :: free_idx
    integer,dimension(:),allocatable      :: new_idx

    contains
    end type res

    type(res),dimension(:),allocatable   :: res_1

    contains
#endif

#if DIM == 2
    type                    :: res
    real(kind=s)            :: doping, dop_charge, charge=0.0_s, current, elconc_avg
    integer,dimension(2)    :: x, y
    integer                 :: loc, res_idx, max_c, del_carriers, add_carriers, test_vel_idx
    integer,dimension(:,:),allocatable    :: free_idx
    integer,dimension(:),allocatable      :: new_idx
    real(kind=s),dimension(:),allocatable      :: test_conc
    integer                 :: test_conc_idx, fill_sp
    real(kind=s),dimension(1000) :: vel = 0.0_s, vel_val=0.0_s
    real(kind=s),dimension(10000) :: test_vel

    contains
    procedure           :: fill_free_idx
    end type res

    type(res),dimension(:),allocatable   :: res_2

    contains
#endif

#if DIM == 1
    subroutine init_reservoir_1()
        integer                :: i, j, res_idx, max_c
        integer,dimension(2)             :: loc
        real(kind=s),dimension(2)        :: doping

        res_len = 32.45_s
        n_res = 0
        res_idx = 0
        loc = 0
        if (contact_type(1) == OHMIC) then
            res_idx = res_idx+1
            doping(res_idx) = donor(1)
            loc(res_idx) = LEFT
        end if
        if (contact_type(2) == OHMIC) then
            res_idx = res_idx+1
            doping(res_idx) = donor(sum(num_nodes))
            loc(res_idx) = RIGHT
        end if
        n_res = res_idx
        allocate(res_1(n_res))
        do i=1,res_idx
            res_1(i)%res_idx = i
            res_1(i)%doping = doping(i)
            res_1(i)%loc = loc(i)
            if (loc(i) == LEFT) then
                res_1(i)%x = (/1,8/)
                res_1(i)%dop_charge = (res_len)*res_1(i)%doping !22.5_s*res_1(i)%doping ! (node_coord(5)+m_dx(1)/2.0_s-node_coord(1))*res_1(i)%doping
            end if
            if (loc(i) == RIGHT) then
                res_1(i)%x = (/sum(num_nodes)-4,sum(num_nodes)/)
                res_1(i)%dop_charge = (node_coord(5))*res_1(i)%doping!(node_coord(sum(num_nodes))-node_coord(sum(num_nodes)-4))*res_1(i)%doping
            end if
            max_c = int(res_1(i)%dop_charge/p_charge)
            res_1(i)%max_c = max_c!*1.08_s*0.921_s
            allocate(res_1(i)%free_idx(NUM_THREADS,res_1(i)%max_c))
            res_1(i)%free_idx = 0
            res_1(i)%del_carriers = 0
            res_1(i)%add_carriers = 0
            allocate(res_1(i)%new_idx(NUM_THREADS))
            res_1(i)%new_idx = 1
            do j=1,res_1(i)%max_c
                call insert_particle_1(res_1(i),0)
            end do
        end do
        do i=1,res_idx
            allocate(res_1(i)%conc(int(res_len/m_dx(1)+2)))
            allocate(res_1(i)%charge(int(res_len/m_dx(1)+2)))
            call fill_free_idx(res_1(i))
        end do
        !print *, res_idx, res_1(1)%x, res_1(1)%doping, res_1(1)%dop_charge, res_1(1)%max_c
        !print *, res_idx, res_1(2)%x, res_1(2)%doping, res_1(2)%dop_charge, res_1(2)%max_c
    end subroutine init_reservoir_1

    subroutine fill_free_idx(this)
        type(res)           :: this
        integer             :: r, k, i, j

        j = 1
        do r=1,n_res
        do k=1,NUM_THREADS
            do i=1,res_1(r)%max_c
                do while (c(j)%valley .ne. 0)
                    j = j+1
                end do
                    res_1(r)%free_idx(k,i) = j
                    j = j+1
            end do
        end do
        end do
    end subroutine

    subroutine insert_particle_1(this,init)
        type(res)           :: this
        integer             :: init, k, omp_id, border

        if (init == 0) then
        k=1
        do while ((c(k)%valley .ne. 0))! .and. (c(k)%op == 1)) !(c(k)%valley > 0)
            k=k+1
        end do
        end if
        if (init == 1) then
            omp_id = omp_get_thread_num()+1
            k = this%free_idx(omp_id,this%new_idx(omp_id))
            c(k)%t_id = omp_id
            c(k)%cidx = k
            this%new_idx(omp_id) = this%new_idx(omp_id)+1
            !print *, this%new_idx(omp_id), omp_id
        end if
                !print *, 'k', k, omp_id
                c(k)%t_id = 1
                call init_res_1(this,c(k),init)
                c(k)%cidx = k
                c(k)%res_id = this%res_idx
                !print *, c(k)%res_id, c(k)%cidx, c(k)%r, c(k)%r_idx, c(k)%k, c(k)%energy, c(k)%mat, c(k)%op
                if (c(k)%r(1) < 20.0_s) then
                    border = 1
                else
                    border = 2
                end if
            !$omp atomic
            particle_in(border) = particle_in(border)+1
    end subroutine insert_particle_1

    subroutine init_res_1(this,cp,init)
        type(res)           :: this
        type(electron)      :: cp
        integer             :: init
        real(kind=s)        :: r, r2, r1, k_abs, ct, l_me, l_nonp

        cp%op = 1
        !if (this%loc == LEFT) then
        if (init == 0) then
            r = rng_uniform(rng(cp%t_id))
            if (this%loc == LEFT) then
                r2 = res_len*r !(node_coord(this%x(2))+m_dx(1))/2.0_s*r !(node_coord(this%x(2))+m_dx(1)/2.0_s-node_coord(this%x(1)))*r+node_coord(this%x(1))
                cp%r = (/r2,0.0_s,0.0_s/)
            end if
            if (this%loc == RIGHT) then
                r2 = dev_x-node_coord(6)*r !(node_coord(this%x(2))+m_dx(1)/2.0_s-node_coord(this%x(1)))*r+node_coord(this%x(1))
                cp%r = (/r2,0.0_s,0.0_s/)
            end if
            cp%ins_coord = r2
        else
            if (this%loc == LEFT) then
                cp%r = (/1e-5_s,0.0_s,0.0_s/)
            end if
            if (this%loc == RIGHT) then
                cp%r = (/dev_x-1e-5_s,0.0_s,0.0_s/)
            end if
        end if
        call get_r_idx(cp%r,cp%r_idx)
        cp%mat = m_material(m_region(cp%r_idx))
        ! located in Gamma-valley
        cp%valley = 1
        l_nonp = nonp(cp%mat,cp%valley)
        l_me = me(cp%mat,cp%valley)
        ! initial energy follows Maxwell-Boltzmann distribution
        r = rng_uniform(rng(cp%t_id))
        cp%energy = -1.5_s*kB*300.0_s*log(r)+0.008 !-1.5_s*kB*t_lattice*log(r)!+0.008_s
        call get_e_idx(cp)

        ! select k
        k_abs = sqrt(2.0_s*me(cp%mat,1)*cp%energy*(1.0_s+nonp(cp%mat,1)*cp%energy))/hb
        !call random_number(r)
        r = rng_uniform(rng(cp%t_id))
        !call random_number(r2)
        r2 = rng_uniform(rng(cp%t_id))
        ct = (1.0_s-2.0_s*r)
        cp%k(1) = k_abs*sqrt(1.0_s-ct*ct)*cos(2.0_s*pi*r2)
        cp%k(2) = k_abs*sqrt(1.0_s-ct*ct)*sin(2.0_s*pi*r2)
        cp%k(3) = k_abs*ct

        ! Give velocity
        cp%vel = hb*cp%k/(l_me*(1.0_s+2.0_s*l_nonp*cp%energy))

        if (this%loc == LEFT) then
            cp%k(1) = abs(cp%k(1))
            cp%vel(1) = abs(cp%vel(1))
        end if
        if (this%loc == RIGHT) then
            cp%k(1) = -abs(cp%k(1))
            cp%vel(1) = -abs(cp%vel(1))
        end if
        r = rng_uniform(rng(cp%t_id))
        cp%e_time = (step-0)*t_step-r*t_step
    end subroutine init_res_1

    subroutine res_drift_1()
        integer                 :: i, j
        real(kind=s)            :: tot_charge_conc, x, y
        real(kind=s)            :: test_conc

        ! Particle drift in reservoirs
        do j=1,n_res
        !$omp parallel private(i) shared(c)
        !$omp do
        !!schedule(static,100)
        do i=1,max_carriers
            if ((c(i)%op == 1) .and. (c(i)%res_id == j)) then
                c(i)%t_id = omp_get_thread_num()+1
                if (res_1(c(i)%res_id)%loc == LEFT) then
                    c(i)%loc_e_field(1) = (el_field(1))!+el_field(2)+el_field(3)+el_field(4))/4.0_s !+el_field(4))/3.0_s
                end if
                if (res_1(c(i)%res_id)%loc == RIGHT) then
                    c(i)%loc_e_field(1) = (el_field(sum(num_nodes)-1)+el_field(sum(num_nodes)-2))/2.0_s !-2
                end if
                do while ((c(i)%synchronized < 0.1) .and. (c(i)%valley > 0) .and. (c(i)%op == 1))
                    call c(i)%drift(1)
                    ! Left
                    if (res_1(c(i)%res_id)%loc == LEFT) then
                    if (c(i)%r(1) < 0) then
                        c(i)%r(1) = abs(c(i)%r(1))
                        if ((c(i)%valley .eq. G) .or. (c(i)%valley .eq. X)) then
                            c(i)%k(1) = c(i)%k(1)*(-1.0_s)
                        end if
                        if (c(i)%valley == L) then
                            call border_L_yz(c(i))
                        end if
                        !print *, 'reflected res', c(i)%cidx, c(i)%r(1:2)
                    end if
                    end if
                    if (res_1(c(i)%res_id)%loc == RIGHT) then
                    if (c(i)%r(1) > dev_x) then
                        c(i)%r(1) = dev_x+(dev_x-c(i)%r(1))
                        if ((c(i)%valley .eq. G) .or. (c(i)%valley .eq. X)) then
                            c(i)%k(1) = c(i)%k(1)*(-1.0_s)
                        end if
                        if (c(i)%valley == L) then
                            call border_L_yz(c(i))
                        end if
                      !  print *, 'reflected res', c(i)%cidx, c(i)%r(1:2)
                    end if
                    end if
                    call c(i)%compute_vel()
                    !print *, 'carrier res drift', c(i)%cidx, c(i)%r(1:2), c(i)%e_time

                    ! Carrier entered device region, inc counter and change particle type to device simulation
                    if (((c(i)%r(1) >= res_len) .and. (res_1(c(i)%res_id)%loc == LEFT)) .or. &
                     &  ((c(i)%r(1) < dev_x-res_len) .and. (res_1(c(i)%res_id)%loc == RIGHT))) then !m_dx(1)*4.0_s+m_dx(1)/2.0_s) then
                        !if (c(i)%r(1) < 0.0_s) then
                        !    c(i)%r(1) = abs(c(i)%r(1))
                        !end if
                        !if (c(i)%r(1) > dev_x) then
                        !    c(i)%r(1) = (dev_x)-(c(i)%r(1)-(dev_x))
                        !end if
                        if (res_1(c(i)%res_id)%loc == LEFT) then
                            c(i)%r(1) = c(i)%r(1)-res_len !1e-4_s
                        else
                            c(i)%r(1) = dev_x-1e-4_s
                        end if
                        call get_r_idx(c(i)%r,c(i)%r_idx)
                        c(i)%op = 0
                        !print *, 'carrier entered device', c(i)%cidx, c(i)%r(1:2), c(i)%e_time
                        !r = rng_uniform(rng(c(i)%t_id))
                        !c(i)%e_time = step*t_step-r*t_step
                        !    do while ((c(i)%synchronized < 0.1) .and. (c(i)%valley > 0))
                        !        call c(i)%drift(0)
                        !    end do
                        !c(i)%synchronized = 0
                        ! Inject a new carrier
                        if (res_1(c(i)%res_id)%del_carriers == 0) then
                            call insert_particle_1(res_1(j),1)
                        else
                            !$omp atomic
                            res_1(c(i)%res_id)%del_carriers = res_1(c(i)%res_id)%del_carriers-1
                        end if
                        if (res_1(c(i)%res_id)%add_carriers > 0) then
                            !print *, res_1(c(i)%res_id)%add_carriers
                            call insert_particle_1(res_1(j),1)
                            !$omp atomic
                            res_1(c(i)%res_id)%add_carriers = res_1(c(i)%res_id)%add_carriers-1
                        end if
                    end if

                end do
                c(i)%synchronized = 0
            end if
        end do
        !$omp end do
        !$omp end parallel
        end do
        do i=1,n_res
            call fill_free_idx(res_1(i))
            res_1(i)%new_idx = 1
            !call res_media_1
        end do
        if (n_res > 0) then
        if (mod(step,50) == 0) then
            !test_conc = (sum(dev%g_el_conc(step-49:step,1))/49_s)!+sum(dev%g_el_conc(step-49:step,2))/49_s)/2.0_s !+ &
                    !& sum(dev%g_el_conc(step-49:step,3))/49_s+sum(dev%g_el_conc(step-49:step,4))/49_s)/4.0_s
            test_conc = (sum(dev%g_el_conc(step-49:step,1))/49_s+sum(dev%g_el_conc(step-49:step,2))/49_s+ & !/2.0_s !+ &
                    & sum(dev%g_el_conc(step-49:step,3))/49_s+sum(dev%g_el_conc(step-49:step,4))/49_s)/4.0_s
            print *, 'avg. el conc', test_conc
            if (test_conc > 1e-4_s) then
                ! Find out how many carriers should be deleted
                res_1(1)%del_carriers = int(res_1(1)%max_c*(test_conc-1e-4_s)/test_conc)
                print *, 'del', res_1(1)%del_carriers
                !print *, int((test_conc-1e-4_s)*m_lx(1)/p_charge)
            end if
            if (test_conc < 1e-4_s) then
                ! Find out how many carriers should be inserted
                res_1(1)%add_carriers = int(res_1(1)%max_c*(1e-4_s-test_conc)/test_conc)
                print *, 'add', res_1(1)%add_carriers
                !print *, int((test_conc-1e-4_s)*m_lx(1)/p_charge)
            end if
        end if
        end if
    end subroutine res_drift_1

    subroutine res_media_1()
        integer         :: i,k
        real(kind=s)    :: x

        res_1(1)%charge = 0.0_s
        res_1(1)%conc = 0.0_s
        do k=1,max_carriers
            if ((c(k)%valley > 0) .and. (c(k)%op == 1)) then
                x = c(k)%r(1)/m_dx(1)
                i = int(x+1)
                !print *, c(k)%r(1), x, i
                res_1(1)%charge(i) = res_1(1)%charge(i)+(1-(x-(i-1)))
                res_1(1)%charge(i+1) = res_1(1)%charge(i+1)+(x-(i-1))
            end if
        end do
        res_1(1)%charge = res_1(1)%charge*p_charge
        res_1(1)%conc = res_1(1)%charge/m_lx(1:10)
        print *, res_1(1)%conc(1:9), sum(res_1(1)%conc(1:9))/9.0_s
    end subroutine res_media_1
#endif

#if DIM == 2
    subroutine init_reservoir_2()
        integer                         :: i, j, start_c, end_c, res_idx, max_c
        integer,dimension(10,3)         :: count_res
        real(kind=s),dimension(10)      :: res_dop

        res_len = 0.7_s-1e-8_s !0.0_s+m_dy/2.0_s-1e-8 !20.4999_s! 51.124_s
        ! count how many reservoirs we need
        res_idx = 0
        count_res = 0
        ! bottom contacts
        start_c = 0
        end_c = 0
        do i=1,num_nodes(1)
            if (i < num_nodes(1)) then
                if ((contact_type(i,1) == OHMIC) .and. (contact_type(i+1,1) == OHMIC) .and. (start_c == 0)) then
                    start_c = i
                end if
            end if
            if (start_c .ne. 0) then
                if (i > 1) then
                !print *, start_c, i
                    if ((contact_type(i,1) .ne. OHMIC) .and. (contact_type(i-1,1) == OHMIC))  then
                        end_c = i-1
                    elseif ((contact_type(i,1) == OHMIC) .and. (i == num_nodes(1))) then
                        end_c = i
                    end if
                end if
            end if
        end do
        if ((start_c > 0) .and. (end_c > 0)) then
            res_idx = res_idx+1
            count_res(res_idx,:) = (/BOTTOM,start_c,end_c/)
            res_dop(res_idx) = donor_density(m_region(int(start_c+(end_c-start_c)/2,8),1))
        end if
              !  print *, 'res start', start_c, end_c, res_dop(res_idx), count_res(res_idx,:)
               ! print *, contact_type(:,num_nodes(2))
        ! top contacts
        start_c = 0
        end_c = 0
        do i=1,int(num_nodes(1))!/3.0_s)
            if (i < num_nodes(1)) then
            if ((contact_type(i,num_nodes(2)) == OHMIC) .and. (contact_type(i+1,num_nodes(2)) == OHMIC) .and. (start_c == 0)) then
                    start_c = i
                end if
            end if
            if (start_c .ne. 0) then
                if (i > 1) then
                    if ((contact_type(i,num_nodes(2)) .ne. OHMIC) .and. (contact_type(i-1,num_nodes(2)) == OHMIC))  then
                        end_c = i-1
                    elseif ((contact_type(i,num_nodes(2)) == OHMIC) .and. (i == num_nodes(1))) then
                        end_c = i
                    end if
                end if
            end if
        end do
        if ((start_c > 0) .and. (end_c > 0)) then
            res_idx = res_idx+1
            count_res(res_idx,:) = (/TOP,start_c,end_c/)
            res_dop(res_idx) = donor_density(m_region(int(start_c+(end_c-start_c)/2,8),num_nodes(2)))
        end if
                !print *, 'res start', start_c, end_c
               ! print *, contact_type(:,num_nodes(2))

!if (n_out == 2) then
!        start_c = 0
!        end_c = 0
!        do i=int(num_nodes(1)/2.0_s),num_nodes(1)
!            if (i < num_nodes(1)) then
!            if ((contact_type(i,num_nodes(2)) == OHMIC) .and. (contact_type(i+1,num_nodes(2)) == OHMIC) .and. (start_c == 0)) then
!                 start_c = i
!            end if
!            end if
!            if (start_c .ne. 0) then
!                if ((contact_type(i,num_nodes(2)) .ne. OHMIC) .and. (contact_type(i-1,num_nodes(2)) == OHMIC))  then
!                    end_c = i-1
!                elseif ((contact_type(i,num_nodes(2)) == OHMIC) .and. (i == num_nodes(1))) then
!                    end_c = i
!                end if
!            end if
!        end do
!        if ((start_c > 0) .and. (end_c > 0)) then
!            res_idx = res_idx+1
!            count_res(res_idx,:) = (/TOP,start_c,end_c/)
!            res_dop(res_idx) = donor_density(m_region(int(start_c+(end_c-start_c)/2,8),num_nodes(2)))
!        end if
!end if

       ! left contacts
!        start_c = 0
!        end_c = 0
!        do i=1,num_nodes(2)
!            if ((contact_type(1,i) == OHMIC) .and. (contact_type(1,i+1) == OHMIC) .and. (start_c == 0)) then
!                start_c = i
!            end if
!            if (start_c .ne. 0) then
!                if ((contact_type(1,i) .ne. OHMIC) .and. (contact_type(1,i-1) == OHMIC))  then
!                    end_c = i-1
!                elseif ((contact_type(1,i) == OHMIC) .and. (i == num_nodes(2))) then
!                    end_c = i
!                end if
!            end if
!        end do
!        if ((start_c > 0) .and. (end_c > 0)) then
!            res_idx = res_idx+1
!            count_res(res_idx,:) = (/LEFT,start_c,end_c/)
!            res_dop(res_idx) = donor_density(m_region(1,int(start_c+(end_c-start_c)/2,8)))
!        end if
!        ! right contacts
!        start_c = 0
!        end_c = 0
!        do i=1,num_nodes(2)
!            if ((contact_type(num_nodes(1),i) == OHMIC) .and. (contact_type(num_nodes(1),i+1) == OHMIC) .and. (start_c == 0)) then
!                start_c = i
!            end if
!            if (start_c .ne. 0) then
!                if ((contact_type(num_nodes(1),i) .ne. OHMIC) .and. (contact_type(num_nodes(1),i-1) == OHMIC))  then
!                    end_c = i-1
!                elseif ((contact_type(num_nodes(1),i) == OHMIC) .and. (i == num_nodes(2))) then
!                    end_c = i
!                end if
!            end if
!        end do
!        if ((start_c > 0) .and. (end_c > 0)) then
!            res_idx = res_idx+1
!            count_res(res_idx,:) = (/RIGHT,start_c,end_c/)
!            res_dop(res_idx) = donor_density(m_region(num_nodes(1),int(start_c+(end_c-start_c)/2,8)))
!        end if

        allocate(res_2(res_idx))
        n_res = res_idx

        ! Get coordinates of ohmic contacts
        do i=1,res_idx
            res_2(i)%res_idx = i
            res_2(i)%doping = res_dop(i)
            res_2(i)%loc = count_res(i,1)
            !res_2(i)%new_idx = 1
            if ((count_res(i,1) == LEFT) .or. (count_res(i,1) == RIGHT)) then
                res_2(i)%x = (/1,8/)
                res_2(i)%y = (/count_res(i,2),count_res(i,3)/)
                res_2(i)%dop_charge = (node_coord(1,8)*(node_coord(2,res_2(i)%y(2))-node_coord(2,res_2(i)%y(1))))*res_2(i)%doping
            end if
            if ((count_res(i,1) == TOP) .or. (count_res(i,1) == BOTTOM)) then
                res_2(i)%y = (/1,8/)
                res_2(i)%x = (/count_res(i,2)+1,count_res(i,3)/)
                !res_2(i)%x(1) = res_2(i)%x(1)+1
                res_2(i)%dop_charge = (res_len*(node_coord(1,res_2(i)%x(2))-node_coord(1,res_2(i)%x(1))))*res_2(i)%doping
                !print *, (node_coord(2,5)),((node_coord(2,res_2(i)%x(2))-node_coord(2,res_2(i)%x(1)))), res_2(i)%doping
                print *, res_2(i)%loc, node_coord(1,res_2(i)%x(1)), node_coord(1,res_2(i)%x(2)), count_res(i,2:3),res_2(i)%doping
            end if
                ! Insert initial particle distribution in reservoir
                !!$omp parallel private(j) shared(c)
                !!$omp do schedule(static,100)
                max_c = int(res_2(i)%dop_charge/p_charge)
                res_2(i)%max_c = max_c*13.88_s
                print *, 'max c', max_c
                res_2(i)%fill_sp = int(res_2(i)%max_c/1.0_s)
                allocate(res_2(i)%free_idx(NUM_THREADS,res_2(i)%fill_sp))
                res_2(i)%free_idx = 0
                allocate(res_2(i)%new_idx(NUM_THREADS))
                res_2(i)%new_idx = 1
                allocate(res_2(i)%test_conc(int(0.006/t_step)))
                res_2(i)%test_conc = 0.0_s
                res_2(i)%test_conc_idx = 0
                res_2(i)%del_carriers = 0
                res_2(i)%add_carriers = 0
                res_2(i)%current = 0.0_s
                res_2(i)%test_vel_idx = 0
                res_2(i)%elconc_avg = 0.0_s
                do j=1,1000
                    res_2(i)%vel(j) = j*2.0_s
                    res_2(i)%vel_val(j) = res_2(i)%vel(j)*exp(-me(m_material(m_region(res_2(i)%x(1),res_2(i)%y(1))),1) &
                        & *(res_2(i)%vel(j))**2.0_s/(2.0_s*kB*T_lattice))
                end do
                if (load_carriers == 0) then
                do j=1,res_2(i)%max_c
                    call insert_particle_2(res_2(i),-1.0_s,0)
                end do
                end if

                !call plot(res_2(i)%vel,res_2(i)%vel_val)
                !!$omp end do
                !!$omp end parallel
        end do
        do i=1,n_res
            call res_2(i)%fill_free_idx
            res_2(i)%new_idx = 1
        end do

        allocate(particle_in(n_out))
        particle_in = 0
        allocate(particle_out(NUM_THREADS,n_out))
        particle_out = 0
        !allocate(dev%particle_flow_in(int(n_step/save_t),n_out))
        !dev%particle_flow_in = 0
        !allocate(dev%particle_flow_out(int(n_step/save_t),n_out))
        !dev%particle_flow_out = 0

   ! print *, res_2(1)%doping, node_coord(1,res_2(1)%x(1)), node_coord(1,res_2(1)%x(2)), &
   !         & res_2(1)%dop_charge, res_2(1)%dop_charge/p_charge, res_2(1)%x, res_2(1)%y
    !print *, res_2(2)%doping, node_coord(1,res_2(2)%x(1)), node_coord(1,res_2(2)%x(2)), &
    !        & res_2(2)%dop_charge, res_2(2)%dop_charge/p_charge, res_2(2)%x, res_2(2)%y
    !print *, res_2(1)%res_idx,res_2(1)%max_c!, res_2(2)%res_idx,res_2(2)%max_c
    end subroutine init_reservoir_2

    subroutine fill_free_idx(this)
        class(res)          :: this
        integer             :: r, k, i, j

        j = 1
        if (this%res_idx > 1) then
            j = res_2(this%res_idx-1)%free_idx(NUM_THREADS,this%fill_sp)
        end if
       ! !$omp parallel private(k) shared(this)
      !  !$omp do
        do k=1,NUM_THREADS
            do i=1,this%fill_sp
                do while (c(j)%valley .ne. 0)
                    j = j+1
                end do
                    this%free_idx(k,i) = j
                    j = j+1
            end do
        end do
     !   !$omp end do
      !  !$omp end parallel
    end subroutine

    subroutine insert_particle_2(this,r,init)
        type(res)           :: this
        integer             :: k, i
        real(kind=s)        :: r, r0, r1, r2
        integer             :: init, omp_id

        k = -1
        !r2 = rng_uniform(rng(omp_id))
        omp_id = 1
        if (init == 0) then
        k=1
        do while ((c(k)%valley .ne. 0))! .and. (c(k)%op == 1)) !(c(k)%valley > 0)
            k=k+1
        end do
        end if
        if (init == 1) then
            omp_id = omp_get_thread_num()+1
            if (this%new_idx(omp_id) == this%fill_sp) then
!                do i=1,n_res
!                    call res_2(i)%fill_free_idx
!                    res_2(i)%new_idx = 1
!                end do
            end if
            k = this%free_idx(omp_id,this%new_idx(omp_id))
            !c(k)%cidx = k
            !print *, k, omp_id
            this%new_idx(omp_id) = this%new_idx(omp_id)+1
        end if
                c(k)%t_id = omp_id
                c(k)%cidx = k
                c(k)%res_id = this%res_idx
                c(k)%res_scat_idx = (/this%x(1),num_nodes(2)/)
                c(k)%MAT = AIR
                do while (c(k)%MAT == AIR)
                    c(k)%op = 1
                    if (this%loc == LEFT) then
                        !c(k)%res_scat_idx = (/this%x(1),num_nodes(2)/)
                        r0 = rng_uniform(rng(c(k)%t_id))
                        r2 = (node_coord(2,this%y(2))-node_coord(2,this%y(1)))*r0+node_coord(2,this%y(1))
                        c(k)%r = (/m_dx+1e-5_s,r2,0.0_s/)
                        c(k)%ins_coord = r2
                        call get_r_idx_2(c(k)%r,c(k)%r_idx)
                        c(k)%mat = m_material(m_region(c(k)%r_idx(1),c(k)%r_idx(2)))
                    end if
                    if (this%loc == TOP) then
                        c(k)%res_scat_idx = (/this%x(1),num_nodes(2)/)
                        r0 = rng_uniform(rng(c(k)%t_id))
                    r2 = (node_coord(1,this%x(2))-node_coord(1,this%x(1))+m_dx-2e-4_s)*r0+node_coord(1,this%x(1))-m_dx_half+1e-4_s !node_coord(1,this%x(1))-m_dx/2.0_s+(node_coord(1,this%x(2))+m_dx-node_coord(1,this%x(1)))*r !
                        r1 = dev_y-1e-5_s
                        c(k)%r = (/r2,r1,0.0_s/)
                        call get_r_idx_2(c(k)%r,c(k)%r_idx)
                        if (c(k)%r_idx(1) < 30) then
                            !print *, c(k)%r_idx, c(k)%r
                        end if
                        c(k)%mat = m_material(m_region(c(k)%r_idx(1),c(k)%r_idx(2)))
                    end if
                    if (this%loc == BOTTOM) then
                        c(k)%res_scat_idx = (/this%x(1),1/)
                        r0 = rng_uniform(rng(c(k)%t_id))
                    r2 = (node_coord(1,this%x(2))-node_coord(1,this%x(1))+m_dx-2e-4_s)*r0+node_coord(1,this%x(1))-m_dx_half+1e-4_s !node_coord(1,this%x(1))-m_dx/2.0_s+(node_coord(1,this%x(2))+m_dx-node_coord(1,this%x(1)))*r !
                        r1 = 1e-5_s
                        c(k)%r = (/r2,r1,0.0_s/)
                        call get_r_idx_2(c(k)%r,c(k)%r_idx)
                        c(k)%mat = m_material(m_region(c(k)%r_idx(1),c(k)%r_idx(2)))
                    end if
                    call init_carrier(this,c(k))
                end do

                !call check_edge(k,this%res_idx)
    end subroutine insert_particle_2

!    subroutine init_res_2(this,cp,init)
!        type(res)           :: this
!        type(electron)      :: cp
!        integer             :: init, n
!        real(kind=s)        :: r, r2, r1, k_abs, ct, l_me, l_nonp
!        real(kind=s)        :: ra,rb, ra_f, kyz
!        real(kind=s),dimension(1) :: val_max, vel_max
!
!        ! located in Gamma-valley
!        cp%valley = 1
!        l_nonp = nonp(cp%mat,cp%valley)
!        l_me = me(cp%mat,cp%valley)
!        ! initial energy follows Maxwell-Boltzmann distribution
!        r = rng_uniform(rng(cp%t_id))
!        cp%energy = -1.5_s*kB*t_lattice*log(r)!*4.0_s!+0.15_s
!!        call get_e_idx(cp)
!!
!        ! select k
!        k_abs = sqrt(2.0_s*me(cp%mat,1)*cp%energy*(1.0_s+nonp(cp%mat,1)*cp%energy))/hb
!        !call random_number(r)
!        r = rng_uniform(rng(cp%t_id))
!        !call random_number(r2)
!        r2 = rng_uniform(rng(cp%t_id))
!        ct = (1.0_s-2.0_s*r)
!        cp%k(1) = k_abs*sqrt(1.0_s-ct*ct)*cos(2.0_s*pi*r2)
!        cp%k(2) = k_abs*sqrt(1.0_s-ct*ct)*sin(2.0_s*pi*r2)
!        cp%k(3) = k_abs*ct
!
!        ! Give velocity
!        cp%vel = hb*cp%k/(l_me*(1.0_s+2.0_s*l_nonp*cp%energy))
!
!        if (this%loc == LEFT) then
!            cp%k(2) = abs(cp%k(2))
!            cp%vel(2) = abs(cp%vel(2))
!        end if
!        if (this%loc == TOP) then
!            cp%k(2) = -abs(cp%k(2))
!            cp%vel(2) = -abs(cp%vel(2))
!        end if
!        if (this%loc == BOTTOM) then
!            cp%k(2) = abs(cp%k(2))
!            cp%vel(2) = abs(cp%vel(2))
!        end if
!
!        r = rng_uniform(rng(cp%t_id))
!        cp%e_time = (step-0)*t_step-r*t_step
!
!!print *, 'init', cp%vel, cp%k, cp%energy
!    end subroutine init_res_2

    subroutine init_carrier(this,cp)
        type(res)           :: this
        type(electron)      :: cp
        integer             :: init, n
        real(kind=s)        :: r, r2, r1, k_abs, l_me, l_nonp
        real(kind=s)        :: ra,rb, ra_f, kyz, A, B, C, test
        real(kind=s),dimension(1) :: val_max, vel_max

        ! located in Gamma-valley
        cp%valley = 1
        l_nonp = nonp(cp%mat,1)
        l_me = me(cp%mat,1)

        val_max = dble(maxval(this%vel_val))
        vel_max = dble(maxval(this%vel))
        n = 1
            do while (n > 0)
                ra = rng_uniform(rng(cp%t_id))
                rb = vel_max(1)*rng_uniform(rng(cp%t_id))
                ra_f = exp(-l_me*(rb)**2.0_s/(2.0_s*kB*T_lattice)) !
                if (ra <= ra_f) then
                    n = -1
                    ra = rng_uniform(rng(cp%t_id))
                    if (ra < 0.5_s) then
                        cp%vel(1) = rb
                    else
                        cp%vel(1) = -rb
                    end if
                end if
            end do
                ra = rng_uniform(rng(cp%t_id))
                rb = vel_max(1)*rng_uniform(rng(cp%t_id))
        n = 1
            do while (n > 0)
                ra = rng_uniform(rng(cp%t_id))
                rb = vel_max(1)*rng_uniform(rng(cp%t_id))
                ra_f = exp(-l_me*(rb)**2.0_s/(2.0_s*kB*T_lattice)) !
                if (ra <= ra_f) then
                    n = -1
                    ra = rng_uniform(rng(cp%t_id))
                    if (ra < 0.5_s) then
                        cp%vel(3) = rb
                    else
                        cp%vel(3) = -rb
                    end if
                end if
            end do
        if (this%loc == TOP) then
            n = 1
                do while (n > 0)
                    ra = val_max(1)*rng_uniform(rng(cp%t_id))
                    rb = vel_max(1)*rng_uniform(rng(cp%t_id))
                    ra_f = rb*exp(-me(cp%mat,1)*(rb-this%current/this%doping)**2.0_s/(2.0_s*kB*T_lattice)) !this%current/this%doping
                    if (ra <= ra_f) then
                        n = -1
                        cp%vel(2) = -rb
                    end if
                end do
        end if
        if (this%loc == BOTTOM) then
            n = 1
                do while (n > 0)
                    ra = val_max(1)*rng_uniform(rng(cp%t_id))
                    rb = vel_max(1)*rng_uniform(rng(cp%t_id))
                    ra_f = rb*exp(-me(cp%mat,1)*(rb-this%current/this%doping)**2.0_s/(2.0_s*kB*T_lattice)) !
                    if (ra <= ra_f) then
                        n = -1
                        cp%vel(2) = rb
                    end if
                end do
        end if

        C = (l_nonp-2.0_s*l_me*l_nonp**2.0_s*dot_product(cp%vel,cp%vel))
        A = (1.0_s-2.0_s*l_nonp*l_me*dot_product(cp%vel,cp%vel))/C
        B = -(dot_product(cp%vel,cp%vel)*l_me/2.0_s)/C
        if (A**2.0_s/4.0_s-B >= 0.0_s) then
            test = -A/2.0_s+sqrt(A**2.0_s/4.0_s-B)
        else
            cp%mat = AIR
        end if

        cp%energy = test
        cp%k = (l_me*(1.0_s+2.0_s*l_nonp*cp%energy)*cp%vel)/hb
        r = rng_uniform(rng(cp%t_id))
        cp%e_time = (step-0)*t_step-r*t_step
        if (E2idx(cp%energy) < 1) then
            cp%mat = AIR
        end if
        if ((abs(cp%vel(1)) < 1e-4_s) .or.  (abs(cp%vel(2)) < 1e-4_s) .or.  (abs(cp%vel(3)) < 1e-4_s)) then
            cp%mat = AIR
        end if
        if ((cp%energy > 0.4_s) .or. (cp%energy < 1e-5_s)) then
            cp%mat = AIR
        end if
       ! call compute_vel(cp)
    end subroutine init_carrier

    subroutine check_edge(i,j)
        integer             :: i, j
        real(kind=s)        :: edge

       if (c(i)%valley > 0) then
            if (node_coord(1,res_2(j)%x(1)) <= 1e-5_s) then
                edge = node_coord(1,res_2(j)%x(1))!-m_dx/2.0_s
            else
                edge = node_coord(1,res_2(j)%x(1))-m_dx/2.0_s+1e-5_s
            end if
            if (c(i)%r(1) < edge) then
                !print *, 'reflected res start', j, c(i)%cidx, c(i)%r(1:2), edge
                c(i)%r(1) =  edge+(edge-c(i)%r(1))
                if ((c(i)%valley .eq. G) .or. (c(i)%valley .eq. X)) then
                    c(i)%k(1) = c(i)%k(1)*(-1.0_s)
                end if
                if (c(i)%valley == L) then
                    call border_L_yz(c(i))
                end if
             !   call c(i)%compute_vel()
                !print *, 'reflected res', c(i)%cidx, c(i)%r(1:2), edge
            end if
            !else
                if ((node_coord(1,res_2(j)%x(2)) >= dev_x-1e-5_s)  .or. (node_coord(1,res_2(j)%x(2)) <= 1e-5_s)) then
                    edge = node_coord(1,res_2(j)%x(2))!+m_dx/2.0_s
                else
                    edge = node_coord(1,res_2(j)%x(2))+m_dx/2.0_s-1e-5_s
                end if
                if (c(i)%r(1) > edge) then
                    c(i)%r(1) = edge+(edge-c(i)%r(1)) !abs(c(i)%r(1))
                    if ((c(i)%valley .eq. G) .or. (c(i)%valley .eq. X)) then
                        c(i)%k(1) = c(i)%k(1)*(-1.0_s)
                    end if
                    if (c(i)%valley == L) then
                        call border_L_yz(c(i))
                    end if
                !    call c(i)%compute_vel()
                end if

            !end if
                if (c(i)%r(2) < 0.0_s) then
                    c(i)%r(2) = -c(i)%r(2)
                    if ((c(i)%valley .eq. G) .or. (c(i)%valley .eq. X)) then
                        c(i)%k(2) = c(i)%k(2)*(-1.0_s)
                    end if
                    if (c(i)%valley == L) then
                        call border_L_xz(c(i))
                    end if
                 !   call c(i)%compute_vel()
                end if

                if (c(i)%r(2) > dev_y) then
                    c(i)%r(2) = (dev_y)-(c(i)%r(2)-(dev_y))
                    if ((c(i)%valley .eq. G) .or. (c(i)%valley .eq. X)) then
                        c(i)%k(2) = c(i)%k(2)*(-1.0_s)
                    end if
                    if (c(i)%valley == L) then
                        call border_L_xz(c(i))
                    end if
                end if
            call c(i)%compute_vel()
            call get_r_idx_2(c(i)%r,c(i)%r_idx)
            end if
    end subroutine check_edge

    subroutine res_drift_2()
        integer                 :: j, i, c_id, action_idx

        do i=1,n_res
            do j=1,1000
                res_2(i)%vel(j) = j*2.0_s
                res_2(i)%vel_val(j) = res_2(i)%vel(j)*exp(-me(m_material(m_region(res_2(i)%x(1),res_2(i)%y(1))),1) &
                    & *(res_2(i)%vel(j)-res_2(i)%current/res_2(i)%doping)**2.0_s/(2.0_s*kB*T_lattice)) !-res_2(i)%current/res_2(i)%doping
            end do
        end do
        ! Particle drift in reservoirs
        !!$omp parallel private(i,j,c_id,action_idx) shared(c,res_2,particle_in,contact_id,res_len)
       ! !$omp parallel private(i,j) shared(c)
       ! !$omp do schedule(dynamic,100)
        !schedule(static,100)
        !do k=1,n_res
        do i=1,max_carriers
        !    if
            if ((c(i)%op == 1)) then
            if (res_2(c(i)%res_id)%add_carriers < 0) then
                res_2(c(i)%res_id)%add_carriers = 0
            endif
            if (res_2(c(i)%res_id)%del_carriers < 0) then
                res_2(c(i)%res_id)%del_carriers = 0
            endif
                c(i)%t_id = omp_get_thread_num()+1
                j = c(i)%res_id
                do while ((c(i)%synchronized < 0.1) .and. (c(i)%valley > 0))

                    !call check_edge(i,j)
                    if (c(i)%valley == 0) then
                        print *, 1, c(i)%energy, c(i)%r(1:2), c(i)%mat
                        exit
                    end if
                    if (c(i)%mat .ne. INGAAS) then
                        print *, 2, c(i)%energy, c(i)%r(1:2), c(i)%mat, c(i)%valley, me(c(i)%mat,c(i)%valley)
                        c(i)%valley = 0
                        exit
                    end if
                    if (E2idx(c(i)%energy) < 1) then
                        print *, 3, c(i)%energy, c(i)%r(1:2), c(i)%mat, c(i)%valley, me(c(i)%mat,c(i)%valley)
                        c(i)%valley = 0
                        exit
                    end if
                    if (res_2(j)%loc == TOP) then
                        c(i)%loc_e_field(1) = el_field(1,c(i)%r_idx(1),num_nodes(2))
                    c(i)%loc_e_field(2) = (el_field(2,c(i)%r_idx(1),num_nodes(2)-1)+el_field(2,c(i)%r_idx(1),num_nodes(2)-1))/2.0_s
                    end if
                    if (res_2(j)%loc == BOTTOM) then
                        c(i)%loc_e_field(1) = el_field(1,c(i)%r_idx(1),1)
                        c(i)%loc_e_field(2) = (el_field(2,c(i)%r_idx(1),2)+el_field(2,c(i)%r_idx(1),2))/2.0_s
                    end if
                    !c(i)%loc_e_field = 1e-9_s
                    call check_edge(i,j)
                    call c(i)%res_drift(action_idx)
                    call check_edge(i,j)
                        call get_r_idx_2(c(i)%r,c(i)%r_idx)
                        c(i)%mat = m_material(m_region(c(i)%r_idx(1),c(i)%r_idx(2)))
                        if (c(i)%mat == AIR) then
                            print *, 'injected in wrong mat', c(i)%cidx, c(i)%r(1:2), c(i)%mat,c(i)%energy
                            c(i)%valley = 0
                            exit
                        end if
                    ! Carrier entered device region, inc counter and change particle type to device simulation
                    if (((c(i)%r(2) < (dev_y-res_len)) .and. (res_2(j)%loc == TOP)) .or. &
                            & ((c(i)%r(2) > res_len) .and. (res_2(j)%loc == BOTTOM))) then
                        if (res_2(j)%loc == TOP) then
                            c(i)%r(2) = c(i)%r(2)+res_len
                        end if
                        if (res_2(j)%loc == BOTTOM) then
                            c(i)%r(2) = c(i)%r(2)-res_len
                        end if
                        call check_edge(i,j)
                        c_id = contact_id(c(i)%r_idx(1),c(i)%r_idx(2))
                        if (c_id == 0) then
                            c_id = 1
                            if (c(i)%r_idx(2) == 2) then
                                c_id = 1
                            end if
                            if (c(i)%r_idx(2) == num_nodes(2)-1) then
                                c_id = 2
                            end if
                        end if

                        if ((el_conc(c(i)%r_idx(1),c(i)%r_idx(2)) > res_2(j)%doping)) then
                            c(i)%valley = 0
                        else
                            !$omp atomic
                            particle_in(c_id) = particle_in(c_id)+1
                            c(i)%op = 0
                            do while ((c(i)%synchronized < 0.1) .and. (c(i)%valley > 0))
                                call check_edge(i,j)
                        call get_r_idx_2(c(i)%r,c(i)%r_idx)
                        c(i)%mat = m_material(m_region(c(i)%r_idx(1),c(i)%r_idx(2)))
                        if (c(i)%mat .ne. INGAAS) then
                            print *, 'injected in wrong mat'
                        end if
                                call c(i)%drift(0)
                                call check_edge(i,j)
                            end do
                        end if

                        ! Inject a new carrier
                        if ((res_2(j)%del_carriers == 0) .and. (res_2(j)%add_carriers == 0)) then
                            call insert_particle_2(res_2(j),c(i)%ins_coord,1)
                        elseif ((res_2(j)%add_carriers > 0) .and. (res_2(j)%del_carriers == 0)) then
                            call insert_particle_2(res_2(j),c(i)%ins_coord,1)
                            call insert_particle_2(res_2(j),c(i)%ins_coord,1)
                            !$omp atomic
                            res_2(j)%add_carriers = res_2(j)%add_carriers-1
                        elseif ((res_2(j)%add_carriers == 0) .and. (res_2(j)%del_carriers > 0)) then
                            !$omp atomic
                            res_2(j)%del_carriers = res_2(j)%del_carriers-1
                        else
                            print *, 'nothing matched', res_2(j)%del_carriers, res_2(j)%add_carriers
                        end if
                    end if
                end do
            end if
            c(i)%synchronized = 0
        end do
        !!$omp end do
        !!$omp end parallel
        !end do
        do i=1,n_res
            call fill_free_idx(res_2(i))
            res_2(i)%new_idx = 1
        end do
    end subroutine res_drift_2

    subroutine test_conc_2()
        integer             :: i, skip_c
        real(kind=s)        :: test_conc

        skip_c = 0
        do i=1,n_res
            res_2(i)%test_conc_idx = res_2(i)%test_conc_idx+1
            if (res_2(i)%loc == TOP) then
           res_2(i)%test_conc(res_2(i)%test_conc_idx) = &
                     & sum(el_conc(res_2(i)%x(1)+skip_c:res_2(i)%x(2)-skip_c,num_nodes(2)-1:num_nodes(2))) &
!                     & /((res_2(i)%x(2)-res_2(i)%x(1)+1-2*skip_c)*2.0_s)!*0.95
                     & /((res_2(i)%x(2)-res_2(i)%x(1)+1-2*skip_c)*2.0_s)!*0.95
            end if
            if (res_2(i)%loc == BOTTOM) then
           res_2(i)%test_conc(res_2(i)%test_conc_idx) = sum(el_conc(res_2(i)%x(1)+skip_c:res_2(i)%x(2)-skip_c,1:2)) &
!                     & /((res_2(i)%x(2)-res_2(i)%x(1)+1-2*skip_c)*2.0_s)!*0.95 !sum(el_conc(41:49,1:5))/9.0_s/5.0_s !
                     & /((res_2(i)%x(2)-res_2(i)%x(1)+1-2*skip_c)*2.0_s)!*0.95 !sum(el_conc(41:49,1:5))/9.0_s/5.0_s !
            end if
                if (res_2(i)%test_conc_idx == int(0.006/t_step)) then
                    test_conc = sum(res_2(i)%test_conc)/dble(int(0.006/t_step)) !50.0_s
                    print *, 'avg carrier conc', test_conc, int(0.006/t_step), ((res_2(i)%x(2)-res_2(i)%x(1)+1-2*skip_c)*2.0_s)
                    res_2(i)%test_conc_idx = 0
                    if (test_conc > res_2(i)%doping) then
                        ! Find out how many carriers should be deleted
                        !res_2(i)%max_c = res_2(i)%max_c+res_2(i)%del_carriers
                        res_2(i)%del_carriers = abs(int(res_2(i)%max_c*(test_conc-res_2(i)%doping)/test_conc ))
                        print *, 'del', res_2(i)%del_carriers, (test_conc-res_2(i)%doping)/test_conc, res_2(i)%max_c
                        !res_2(i)%max_c = res_2(i)%max_c-res_2(i)%del_carriers
                        !print *, int((test_conc-1e-4_s)*m_lx(1)/p_charge)
                    end if
                    if (test_conc  < res_2(i)%doping) then
                        ! Find out how many carriers should be inserted
                        !res_2(i)%max_c = res_2(i)%max_c-res_2(i)%add_carriers
                        res_2(i)%add_carriers = abs(int(res_2(i)%max_c*(res_2(i)%doping-test_conc )/test_conc))
                        print *, 'add', res_2(i)%add_carriers, (res_2(i)%doping-test_conc )/test_conc, res_2(i)%max_c
                        !print *, int((test_conc-1e-4_s)*m_lx(1)/p_charge)
                    end if
                    !res_2(i)%max_c = res_2(i)%max_c+res_2(i)%add_carriers
                    if (res_2(i)%add_carriers > res_2(i)%del_carriers) then
                        res_2(i)%add_carriers = res_2(i)%add_carriers-res_2(i)%del_carriers
                        res_2(i)%del_carriers = 0
                    end if
                    if (res_2(i)%add_carriers < res_2(i)%del_carriers) then
                        res_2(i)%del_carriers = res_2(i)%del_carriers-res_2(i)%add_carriers
                        res_2(i)%add_carriers = 0
                    end if
                else
                end if
        end do
    end subroutine test_conc_2

    subroutine floating_conc()
        integer             :: i, skip_c, start_step, l_step
        real(kind=s)        :: elconc_avg
        integer             :: del_carriers, add_carriers

        skip_c = 2
        i = 1
        start_step = 350
        l_step = step-start_step
        del_carriers = 0
        add_carriers = 0
        if (l_step > 0) then
        res_2(i)%elconc_avg = res_2(i)%elconc_avg*(l_step-1)/l_step+el_conc(45,1)/l_step! &
                        !& sum(el_conc(res_2(i)%x(1)+skip_c:res_2(i)%x(2)-skip_c,num_nodes(2)-0:num_nodes(2))) &
                        !& /((res_2(i)%x(2)-res_2(i)%x(1)+1-2*skip_c)*1.0_s)/step
        print *, l_step, res_2(i)%elconc_avg, res_2(i)%max_c
        else
            res_2(i)%elconc_avg = res_2(i)%doping
        end if
        res_2(i)%test_conc_idx = res_2(i)%test_conc_idx+1
        if (res_2(i)%test_conc_idx == int(0.05/t_step)+1) then
            res_2(i)%test_conc_idx = 1
            if (res_2(i)%elconc_avg > res_2(i)%doping) then
                ! Find out how many carriers should be deleted
                del_carriers = abs(int(res_2(i)%max_c*(res_2(i)%elconc_avg -res_2(i)%doping)/res_2(i)%elconc_avg ))
                print *, 'del', del_carriers
                !print *, int((test_conc-1e-4_s)*m_lx(1)/p_charge)
            end if
            if (res_2(i)%elconc_avg  < res_2(i)%doping) then
                ! Find out how many carriers should be inserted
                add_carriers = abs(int(res_2(i)%max_c*(res_2(i)%doping-res_2(i)%elconc_avg )/res_2(i)%elconc_avg ))
                print *, 'add', add_carriers
                !print *, int((test_conc-1e-4_s)*m_lx(1)/p_charge)
            end if
        res_2(i)%max_c = res_2(i)%max_c-del_carriers+add_carriers
        res_2(1)%add_carriers = add_carriers
        res_2(i)%del_carriers = del_carriers
        end if

        !print *, res_2(1)%add_carriers, res_2(1)%del_carriers
    end subroutine floating_conc

!    subroutine reset_contacts()
!        integer             :: i

!            contact_type(1:num_nodes(1),num_nodes(2)) = INSULATOR
!            contact_type(1:40,num_nodes(2)) = OHMIC
!            contact_pot(1:40,num_nodes(2)) = 0.0_s

            !contact_type((num_nodes(1)-40):num_nodes(1),num_nodes(2)) = OHMIC
!            contact_pot((num_nodes(1)-40):num_nodes(1),num_nodes(2)) = 0.8_s
!            deallocate(res_2)
!            do i=1,max_carriers
!                if (c(i)%op == 1) then
!                    c(i)%valley = 0
!                end if
!            end do
!            call init_reservoir_2()
!    end subroutine reset_contacts
#endif
end module reservoir
