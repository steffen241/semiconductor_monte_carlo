module tunneling
    use physconst
    use config
    use materialdef
    !use poisson
    use gnufor2
    use fgsl
    use rng_mod
    use device
    !use mc_core
    use, intrinsic :: iso_c_binding
    implicit none

    !real(kind=s)    :: barrier_l = 19.0_s

    type                            :: tunnel
    integer                         :: steps_compute_tunnel
    integer,dimension(2)            :: t_x_idx, t_y_idx
    real(kind=s),dimension(2)       :: t_x, t_y
    real(kind=s),dimension(1000)    :: x_interp, t_dest
    real(kind=s),dimension(:),allocatable :: epot, x
    real(kind=s),dimension(:,:),allocatable :: epot_org
!    real(kind=s),dimension(:,:),allocatable :: tunnel_prob, t_dest, energy
    real(kind=s),dimension(:,:,:),allocatable :: elpot_save
    real(kind=s),dimension(1000,3) :: tunnel_prob, energy, epot_interp
    integer                         :: direction, e_idx, x_idx, elpot_save_idx, idx, current_valley
    integer            :: pot_min
        integer,dimension(3)         :: save_tun = 1

    contains
!    procedure           :: fill_free_idx
    end type tunnel

    type(tunnel),dimension(:),allocatable   :: tun

    contains
#if DIM == 2
    subroutine init_tunneling()
        integer                 :: i, j
        integer                 :: s_idx, e_idx

        allocate(tun(n_tunnel))

        s_idx = 0
        e_idx = 0
        if (n_tunnel > 0) then
            do j=1,n_tunnel
            if (tunnel_x(1,j) > 0) then
                ! Find the start and end of barrier
                do i=1,num_nodes(2)-1
                if ((s_idx == 0) .and. (m_material(m_region(tunnel_x(1,j),i)) .ne. m_material(m_region(tunnel_x(1,j),i+1))) .and. &
                      & (emin(m_material(m_region(tunnel_x(1,j),i))) < emin(m_material(m_region(tunnel_x(1,j),i+1))))) then
                        s_idx = i
                    end if
                if ((e_idx == 0) .and. (m_material(m_region(tunnel_x(1,j),i)) .ne. m_material(m_region(tunnel_x(1,j),i+1))) .and. &
                      & (emin(m_material(m_region(tunnel_x(1,j),i))) > emin(m_material(m_region(tunnel_x(1,j),i+1))))) then
                        e_idx = i
                    end if
                end do
            end if
            print *, s_idx
            end do
        end if

        do i=1,n_tunnel
            tun(i)%steps_compute_tunnel = int(compute_tunnel/t_step)
            tun(i)%idx = i
            tun(i)%direction = tunnel_dir(i)
            tun(i)%t_x_idx = tunnel_x(:,i)
            tun(i)%t_y_idx = tunnel_y(:,i)
            tun(i)%t_y(1) = node_coord(2,tun(i)%t_y_idx(1)+1)-m_dy/2.0_s
            tun(i)%t_y(2) = node_coord(2,tun(i)%t_y_idx(2)-1)+m_dy/2.0_s
!            do j=1,1000
!                tun(i)%energy(j) = maxval(emin, MASK = (emin < 1.0_s))/1000.0_s*j
!            end do
            allocate(tun(i)%epot_org(tun(i)%t_y_idx(2)-tun(i)%t_y_idx(1)+1,3))
            allocate(tun(i)%epot(tun(i)%t_y_idx(2)-tun(i)%t_y_idx(1)+1))
            tun(i)%epot = 0.0_s
            allocate(tun(i)%x(tun(i)%t_y_idx(2)-tun(i)%t_y_idx(1)+1))
            tun(i)%x = node_coord(2,tun(i)%t_y_idx(1):tun(i)%t_y_idx(2))!-m_dy/2.0_s
            !allocate(tun(i)%tunnel_prob(tun(i)%t_y_idx(2)-tun(i)%t_y_idx(1)+2,1000))
            !allocate(tun(i)%t_dest(tun(i)%t_y_idx(2)-tun(i)%t_y_idx(1)+2,1000))
            !allocate(tun(i)%energy(tun(i)%t_y_idx(2)-tun(i)%t_y_idx(1)+2,1000))
            allocate(tun(i)%elpot_save(tun(i)%steps_compute_tunnel,tun(i)%t_x_idx(2)-tun(i)%t_x_idx(1)+1, &
                        & tun(i)%t_y_idx(2)-tun(i)%t_y_idx(1)+1))
            tun(i)%elpot_save_idx = 1
            tun(i)%tunnel_prob = 0.0_s
            print *, size(tun(i)%epot,1),size(tun(i)%elpot_save,3), tun(i)%t_x_idx
            !print *, tun(i)%t_x_idx, tun(i)%t_y
            !print *, m_material(m_region(tun(i)%t_x_idx(1),tun(i)%t_y_idx(2)))
            !print *, m_material(m_region(tun(i)%t_x_idx(1),tun(i)%t_y_idx(2)-1))
        end do
        !if (tunnel_hdf5 == 1) then
        !allocate(dev%tun_energy(n_tunnel,int(t_sim/0.02)+1,1000,3))
        !dev%tun_energy = 0.0_s
        !allocate(dev%tun_prob(n_tunnel,int(t_sim/0.02)+1,1000,3))
        !dev%tun_prob = 0.0_s
        !allocate(dev%tun_dest(n_tunnel,int(t_sim/0.02)+1,1000))
        !dev%tun_dest = 0.0_s
        !end if
    end subroutine init_tunneling

    subroutine calc_tunnelrate(this,valley)
        type(tunnel)                          :: this
        real(kind=s),dimension(20)            :: epot, x
        real(kind=s)                          :: t_capch, t_tfe, test, t_max, cb_offs
        !real(kind=s),dimension(1000)          :: e_eval, x_eval
        !real(kind=s),dimension(1000)          :: tunnel_prob, energy, t_dest
        integer                               :: valley, i, j, mat_dest, mat_org
        integer,dimension(1)                  :: tfe_idx, e_idx, xi_idx

        this%current_valley = valley
        ! Fill electric potential !el_pot2d
        !this%epot = el_pot2d(int(this%t_x_idx(2)-this%t_x_idx(1)/2.0_s),this%t_y_idx(1):this%t_y_idx(2))

        do i=1,this%steps_compute_tunnel
            this%epot = this%epot+this%elpot_save(i,int((this%t_x_idx(2)-this%t_x_idx(1))/2.0_s),:)
        end do
            this%epot = this%epot/dble(this%steps_compute_tunnel)
        mat_dest = m_material(m_region(this%t_x_idx(1)+1,this%t_y_idx(1)+2))
        if (this%direction == 1) then
            mat_org = m_material(m_region(this%t_x_idx(1)+1,this%t_y_idx(1)-1))
        elseif (this%direction == 2) then
            mat_org = m_material(m_region(this%t_x_idx(1)+1,this%t_y_idx(2)+1))
        end if
        !print *, mat_dest, mat_org, emin(mat_org), emin(mat_dest)
        !print *, emin(mat_dest)+valley_offs(mat_dest,L)-(emin(mat_org)+valley_offs(mat_org,L))
        !print *, this%epot
        !print *, size(sum(this%elpot_save(:,14,:), DIM = 1))
        !print *, this%epot(10), int(this%t_x_idx(2)-this%t_x_idx(1)/2.0_s)

        ! Normalize conduction band energy
        if (this%direction == 2) then
            this%epot = (-this%epot+(this%epot(size(this%epot))))
        elseif (this%direction == 1) then
            this%epot = (-this%epot+(this%epot(1)))
        end if
        cb_offs = emin(mat_dest)+valley_offs(mat_dest,valley)-(emin(mat_org)+valley_offs(mat_org,valley))
        this%epot_org(:,this%current_valley) = this%epot+cb_offs
        this%epot(2:size(this%epot)-1) = this%epot(2:size(this%epot)-1)+cb_offs
        print *, size(this%epot), this%epot
        !print *, this%x, this%epot_org(:,valley)
!        print *, size(this%x), this%x
if (this%idx == 5) then
        !call plot(this%x,this%epot_org(:,1))
end if
!        do i=1,1000
!            this%x_interp(i) = this%t_y(1)+(this%t_y(2)-(this%t_y(1)))/1000.0_s*dble(i-1)
!            this%epot_interp(i) = get_potential(this,this%x_interp(i))
!            print *, i, this%x_interp(i), this%epot_interp(i)
!        end do
!        do i=size(this%epot),2,-1
!            !print *, this%x(i), this%epot(i)
!            ! If slope is positive (no barrier in tunneling direction there is no tunneling
!            if (this%epot(i) > this%epot(i-1)) then
!                this%tunnel_prob(i,:) = 0.0_s
!            else
!                ! Find energy where direct tunneling turns to TFE
!                print *, 'tunneling possible'
!                xi_idx = minloc(abs(this%x_interp-this%x(i))) !, MASK=(this%x_interp-this%x(i) > 0.0_s))
!                print *, xi_idx, this%x(i), this%x_interp(xi_idx(1))
!                t_tfe = minval(this%epot_interp(1:xi_idx(1)))
!                tfe_idx = minloc(this%epot_interp(1:xi_idx(1)))
!                t_max = maxval(this%epot_interp(tfe_idx(1):xi_idx(1)))
!                print *, 'Max energy for direct tunneling', t_tfe, this%x_interp(xi_idx(1)), t_max
!                this%tunnel_prob(i,:) = 0.0_s
!                do j=1,1000
!                    this%energy(i,j) = t_max/1000.0_s*dble(j-1)
!                    this%e_idx = j
!                    this%x_idx = i
!                    ! Direct tunneling
!                    if (this%energy(i,j) < t_tfe) then
!                        !print *, this%energy(i,j)
!                        this%t_dest(i,j) = this%x_interp(xi_idx(1))!this%t_y(2)
!                        this%tunnel_prob(i,j) = exp(-2.0_s/hb*int_f_barrier(this,this%t_dest(i,j)))
!                   ! TFE, calculation until minimum energy in barrier
!        !            else
!        !                e_idx = minloc(this%epot_interp(tfe_idx(1):i)-this%energy(i), &
!        !                & MASK=((this%epot_interp(tfe_idx(1):i)-this%energy(i)) .gt. 0))+tfe_idx(1)
!        !                this%tunnel_prob(i) = exp(-2.0_s/hb*int_f_barrier(this,this%x_interp(e_idx(1))))
!                    end if
!                end do
!            end if
!        end do
!        call plot(this%energy(size(this%epot),:),this%tunnel_prob(size(this%epot),:))
!        call plot(this%energy(1,:),this%tunnel_prob(1,:))
        ! Minimum energy needed for direct tunneling: t_capch
        ! Maximum energy for direct tunneling (above TFE assumed: t_tfe
        do i=1,1000
            this%x_interp(i) = this%t_y(1)+(this%t_y(2)-this%t_y(1))/1000.0_s*dble(i)
            this%epot_interp(i,valley) = get_potential(this,this%x_interp(i))
        end do
        !call plot(this%x_interp,this%epot_interp)
        t_tfe = minval(this%epot_interp(:,valley))
        tfe_idx = minloc(this%epot_interp(:,valley))
        this%pot_min = tfe_idx(1)+1
        if (this%pot_min == 1001) then
            this%pot_min = 1000
        end if
        !print *, this%pot_min, t_tfe
        if (this%direction == 2) then
            t_max = maxval(this%epot_interp(tfe_idx(1):1000,valley))
        else
            t_max = maxval(this%epot_interp(1:tfe_idx(1),valley))
        end if
        !print *, t_max, t_tfe, tfe_idx
        !print *, get_potential(this,this%t_y(1)), get_potential(this,this%t_y(2))

        ! Calculate now tunnel probability for cap->channel / cap->barrier TFE
        !this%tunnel_prob = 0.0_s
        this%t_dest = 0.0_s
        do i=1,1000
            this%energy(i,valley) = t_max/1000.0_s*dble(i-1)
            this%e_idx = i
            ! Direct tunneling
            if (this%energy(i,valley) < t_tfe) then
                if (this%direction == 2) then
                    this%t_dest(i) = this%t_y(1)
                    this%tunnel_prob(i,valley) = exp(-2.0_s/hb*int_f_barrier(this,this%t_y(1)))
                elseif (this%direction == 1) then
                    this%t_dest(i) = this%t_y(2)
                    this%tunnel_prob(i,valley) = exp(-2.0_s/hb*int_f_barrier(this,this%t_y(2)))
                end if
           ! TFE, calculation until minimum energy in barrier
            else
                if (this%direction == 2) then
                    j = 1000
                    do while (this%epot_interp(j,valley) > this%energy(i,valley))
                        e_idx(1) = j
                        j = j-1
                    end do
                    this%t_dest(i) = this%x_interp(tfe_idx(1)) ! (e_idx(1))
                    this%tunnel_prob(i,valley) = exp(-2.0_s/hb*int_f_barrier(this,this%x_interp(e_idx(1))))
                elseif (this%direction == 1) then
                    j = 1
                    do while (this%epot_interp(j,valley) > this%energy(i,valley))
                            e_idx(1) = j
                            j = j+1
                    end do
                    this%t_dest(i) = this%x_interp(tfe_idx(1)) !(e_idx(1))
                    this%tunnel_prob(i,valley) = exp(-2.0_s/hb*int_f_barrier(this,this%x_interp(e_idx(1))))
                end if
            end if
        end do
        !if (tunnel_hdf5 == 1) then
        !dev%tun_prob(this%idx,this%save_tun(valley),:,valley) = this%tunnel_prob(:,valley)
        !dev%tun_energy(this%idx,this%save_tun(valley),:,valley) = this%energy(:,valley)
        !dev%tun_dest(this%idx,this%save_tun(valley),:) = this%t_dest
        !this%save_tun(valley) = this%save_tun(valley)+1
        !end if
        !call plot(this%energy(:,valley),this%tunnel_prob(:,valley))
        !call plot(this%energy,this%t_dest)
    end subroutine calc_tunnelrate

    function f_barrier( x, params ) bind(c)
        type(tunnel),pointer            :: this
        real( kind = c_double )        :: f_barrier, epot
        real( kind = c_double ), value :: x
        type( c_ptr ), value           :: params
        real(kind=fgsl_double),pointer    :: p

        params = params
        !epot = 0.52_s
        !energy = 0.3_s
        call c_f_pointer(params,this)
        if (get_potential(this,x)- this%energy(this%e_idx,this%current_valley) < 0.0_s) then
            print *, x, get_potential(this,x), this%energy(this%e_idx,this%current_valley)
            !print *, get_potential(x), x, get_potential(x)-p
        end if
        !print *, x, get_potential(this,x)
        f_barrier = sqrt(2.0_s*me(INALAS,this%current_valley)*(get_potential(this,x) &
            & -this%energy(this%e_idx,this%current_valley)))
    end function f_barrier

    function get_potential(this,x) result(el_pot)
        type(tunnel)                       :: this
        real(fgsl_double)                  :: el_pot, x
        !real(fgsl_double),dimension(20)    :: g_pot, xi
        integer                            :: i, direction

        ! Variables for GSL interpolation
        integer(fgsl_int)                  :: s
        integer(fgsl_size_t)               :: nmax
        type(fgsl_interp_accel)            :: acc
        type(fgsl_spline)                  :: spline

!        do i=1,20
!            xi(i) = dble(i-1)
!        end do
!        g_pot = (/0.35066277,  0.37853166,  0.40637618,  0.43426359,  0.44757858, &
!                       & 0.44631988,  0.43048924,  0.40008366,  0.36966932,  0.33924225, &
!                       & 0.30879721,  0.27832738,  0.24782376,  0.21727423,  0.18666238, &
!                       & 0.1559656 ,  0.12515233,  0.09521475,  0.06697019,  0.04261201 /)
        nmax = size(this%epot)

        spline = fgsl_spline_alloc(fgsl_interp_cspline,nmax)
        s = fgsl_spline_init(spline,this%x,this%epot_org(:,this%current_valley),nmax)
        el_pot = fgsl_spline_eval(spline,x,acc)
        call fgsl_spline_free(spline)
    end function get_potential
!
    function int_f_barrier(this, x) result(res)
        type(tunnel),target          :: this
        integer                      :: direction ! 1: Cap->channel; 2: Channel->Cap
        real( kind = fgsl_double )   :: res, error
        real(kind=s)                 :: x !1,x2
        integer( kind = fgsl_size_t) :: neval
        integer( kind = fgsl_int)    :: i
        type( fgsl_function )        :: func
        type(c_ptr)                  :: p
        real(kind=fgsl_double),target           :: energy

        !energy = 0.515_fgsl_double
        p = c_loc(this)
        func = fgsl_function_init( f_barrier, p )

        if (this%direction == 2) then
            i = fgsl_integration_qng ( func,                &
                                     x,                     &
                                     this%t_y(2),           &
                                     1e-2_fgsl_double,      &
                                     1e-2_fgsl_double,      &
                                     res, error, neval )
        elseif (this%direction == 1) then
            i = fgsl_integration_qng ( func,              &
                                     this%t_y(1),       &
                                     x,           &
                                     1e-2_fgsl_double,      &
                                     1e-2_fgsl_double,      &
                                     res, error, neval )
        end if
    end function int_f_barrier
#endif
end module tunneling
