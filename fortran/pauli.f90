module pauli
    use physconst
    use config
    use materialdef
    use gnufor2
    use fgsl
    use, intrinsic :: iso_c_binding
    implicit none

#if DIM == 0
    real(kind=s)        :: n_avg, k_drift, Ed != 1.0e-4_s
    type                            :: pep
    integer                         :: step, num_steps
    integer,dimension(:,:,:),allocatable :: num_particles
    ! We need variables for saving electron conc, average drift wave vector and average energy (each valley)
    real(kind=s),dimension(3)   :: n_conc, E_avg
    real(kind=s),dimension(3,3) :: k_drift
    real(kind=s),dimension(3)   :: fermi_lvl, el_temp
    integer                  :: mat, valley
    contains

    end type pep
#endif

#if DIM == 2
    real(kind=s)        :: n_avg, k_drift, Ed != 1.0e-4_s
    type                            :: pep
    integer                         :: steps_compute_pep, steps_save_pep, step, num_steps
    integer,dimension(:,:,:),allocatable :: num_particles
    ! We need variables for saving electron conc, average drift wave vector and average energy (each valley)
    real(kind=s),dimension(:,:,:),allocatable   :: n_conc, E_avg
    real(kind=s),dimension(:,:,:,:),allocatable :: k_drift
    real(kind=s),dimension(:,:,:,:),allocatable   :: step_n_conc, step_E_avg, test_energy
    real(kind=s),dimension(:,:,:,:,:),allocatable :: step_k_drift
    real(kind=s),dimension(:,:,:),allocatable     :: fermi_lvl, el_temp
    integer                                       :: mat, valley
    !real(kind=s)                    :: Ed
    contains

    end type pep
#endif

    type(pep)   :: p

    contains

#if DIM == 0
    subroutine init_pauli()
        p%step = 1
        p%num_particles = 0
        p%n_conc = 0.0_s
        p%k_drift = 0.0_s
        p%E_avg = 0.0_s
        p%fermi_lvl = -10.0_s
        p%el_temp = -10.0_s

        p%mat = 0
        p%valley = 0
    end subroutine init_pauli

    subroutine calc_fd()
        integer         :: i!, x_idx, y_idx, n_steps, j
        real(kind=s)    :: fermi, Te

        do i=1,3
            if (p%n_conc(i) > 0.0_s) then
                Ed = p%E_avg(i)
                call fd(i,fermi,Te)
                p%fermi_lvl(i) = fermi
                p%el_temp(i) = Te
            else
                p%fermi_lvl(i) = -10.0_s
                p%el_temp(i) = -2.0_s
            end if
        end do
    end subroutine calc_fd

    subroutine fd(valley,Ef,Te)
        integer(fgsl_int) :: status, i, valley
        integer,dimension(1)                   :: idx
        real(kind=s)                :: Ef, Te
        real(kind=s),dimension(100) :: x_guess, val_guess

        real(fgsl_double), parameter :: eps5 = 1.0d-5
        integer(fgsl_int), parameter :: itmax_root = 150
        type(c_ptr) :: ptr
        real(fgsl_double) :: ri, ra, xlo, xhi
        type(fgsl_function) :: func
        type(fgsl_root_fsolver) :: root_fslv

        ! Set variables
        n_avg = p%n_conc(valley)
        Ed = p%E_avg(valley)
        func = fgsl_function_init(fermi_level, ptr)
        root_fslv = fgsl_root_fsolver_alloc(fgsl_root_fsolver_bisection)
        val_guess = -1000.0_s
        do i=1,100
        x_guess = 6.0_s/100.0_s*i-3.0_s
        val_guess = fermi_level(x_guess(i),ptr)
        if (val_guess(i) > 1e-2) then
            exit
        end if
        end do
        func = fgsl_function_init(fermi_level, ptr)
        root_fslv = fgsl_root_fsolver_alloc(fgsl_root_fsolver_brent)
        status = fgsl_root_fsolver_set(root_fslv, func, -8.0_fgsl_double, &
               x_guess(i))
        i = 0
        do
            i = i + 1
            status = fgsl_root_fsolver_iterate(root_fslv)
            if (status /= fgsl_success .or. i > itmax_root) then
                exit
            end if
            ra = fgsl_root_fsolver_root(root_fslv)
            xlo = fgsl_root_fsolver_x_lower(root_fslv)
            xhi = fgsl_root_fsolver_x_upper(root_fslv)
            status = fgsl_root_test_interval (xlo, xhi, 0.0_fgsl_double, eps5)
            if (status == fgsl_success) exit
        end do
        call fgsl_root_fsolver_free(root_fslv)
        call fgsl_function_free(func)
        Ef = ra
        Te = electron_temp(ra)
        print *, Ef, Te
    end subroutine
#endif

#if DIM == 2
    subroutine init_pauli()
        allocate(p%num_particles(3,num_nodes(1),num_nodes(2)))
        allocate(p%n_conc(3,num_nodes(1),num_nodes(2)))
        allocate(p%E_avg(3,num_nodes(1),num_nodes(2)))
        allocate(p%k_drift(3,num_nodes(1),num_nodes(2),3))
        allocate(p%fermi_lvl(3,num_nodes(1),num_nodes(2)))
        allocate(p%el_temp(3,num_nodes(1),num_nodes(2)))

        p%step = 1
        p%num_particles = 0
        p%n_conc = 0.0_s
        p%k_drift = 0.0_s
        p%E_avg = 0.0_s
        p%fermi_lvl = -10.0_s
        p%el_temp = -10.0_s
        p%steps_save_pep = int(save_pep/t_step+0.5_s)
        p%steps_compute_pep = int(compute_pep/t_step+0.5_s)
        p%num_steps = p%steps_compute_pep/p%steps_save_pep
        print *, p%num_steps
        allocate(p%step_E_avg(p%num_steps,3,num_nodes(1),num_nodes(2)))
        allocate(p%step_n_conc(p%num_steps,3,num_nodes(1),num_nodes(2)))
        allocate(p%step_k_drift(p%num_steps,3,num_nodes(1),num_nodes(2),3))
        allocate(p%test_energy(p%num_steps,3,num_nodes(1),num_nodes(2)))
        p%test_energy = 0.0_s
        p%step_E_avg = 0.0_s
        p%step_n_conc = 0.0_s
        p%step_k_drift = 0.0_s
        print *, p%steps_save_pep, p%steps_compute_pep, p%num_steps
        p%mat = 0
        p%valley = 0
    end subroutine init_pauli

    subroutine calc_fd()
        integer         :: i, x_idx, y_idx, n_steps, j
        real(kind=s)    :: fermi, Te

        p%step = 1

        p%E_avg = 0.0_s
        p%n_conc = 0.0_s
        p%k_drift = 0.0_s
        n_steps = 0
        do i=1,1 !3

            do x_idx=1,num_nodes(1)
                do y_idx=1,num_nodes(2)
                    p%n_conc(i,x_idx,y_idx) = sum(p%step_n_conc(:,i,x_idx,y_idx),1)/p%num_steps*p_charge/(m_lx(x_idx)*m_ly(y_idx))
                end do
            end do
            !p%E_avg(i,:,:) = sum(p%step_E_avg(:,i,:,:),1)/p%num_steps

        do x_idx=1,num_nodes(1)
            do y_idx=1,num_nodes(2)
                do j=1,p%num_steps
                    if (p%step_E_avg(j,i,x_idx,y_idx) > 0.0_s) then
                        p%E_avg(i,x_idx,y_idx) = p%E_avg(i,x_idx,y_idx)+p%step_E_avg(j,i,x_idx,y_idx)
                        n_steps = n_steps+1
                    end if
                end do
                !p%k_drift(i,:,:,1) = sum(p%step_k_drift(:,i,:,:,1),1)/dble(n_steps)! p%num_steps
                !p%k_drift(i,:,:,2) = sum(p%step_k_drift(:,i,:,:,2),1)/dble(n_steps)!p%num_steps
                !p%k_drift(i,:,:,3) = sum(p%step_k_drift(:,i,:,:,3),1)/dble(n_steps)!p%num_steps
                p%E_avg(i,x_idx,y_idx) = p%E_avg(i,x_idx,y_idx)/dble(n_steps)
                n_steps = 0
            end do
        end do
        end do
        !print *, p%k_drift(1,30,77,:), p%n_conc(:,30,77), p%E_avg(:,30,77)
        !!$omp parallel private(x_idx)
        !!$omp do schedule(dynamic)
        do x_idx=1,num_nodes(1)
            do y_idx=1,num_nodes(2)
                p%mat = m_material(m_region(x_idx,y_idx))
                do i=1,3
                p%valley = i
        !        if (step*t_step > compute_pep) then
                if ((p%mat .ne. AIR) .and. (p%n_conc(i,x_idx,y_idx) > 0.0_s) .and. (p%E_avg(i,x_idx,y_idx) > 0.0_s) &
                   & .and. (p%E_avg(i,x_idx,y_idx) > 1e-5_s) .and. (i==1)) then
                    call fd(i,x_idx,y_idx,fermi,Te)
                    p%fermi_lvl(i,x_idx,y_idx) = fermi
                    p%el_temp(i,x_idx,y_idx) = Te
                    !print *, 'compl', x_idx, y_idx, fermi, Te
                else
                    !print *, i,x_idx,y_idx,p%step_E_avg(:,i,x_idx,y_idx)
                    p%fermi_lvl(i,x_idx,y_idx) = -2.0_s
                    p%el_temp(i,x_idx,y_idx) = -2.0_s
                    !print *, 'no calculation possible', x_idx, y_idx
                end if
         !       end if
                end do
            end do
        end do
        p%step_E_avg = 0.0_s
        p%step_n_conc = 0.0_s
        p%step_k_drift = 0.0_s
        p%test_energy = 0.0_s

        !!$omp end do
        !!$omp end parallel
    end subroutine calc_fd

    subroutine fd(valley,x_idx,y_idx,Ef,Te)
        integer(fgsl_int) :: status, i, valley, x_idx, y_idx
        integer,dimension(1)                   :: idx
        real(kind=s)                :: Ef, Te
        real(kind=s),dimension(100) :: x_guess, val_guess

        real(fgsl_double), parameter :: eps5 = 1.0d-4
        integer(fgsl_int), parameter :: itmax_root = 1500
        type(c_ptr) :: ptr
        real(fgsl_double) :: ri, ra, xlo, xhi
        type(fgsl_function) :: func
        type(fgsl_root_fsolver) :: root_fslv

        type(fgsl_error_handler_t)    :: error_h

        error_h = fgsl_set_error_handler_off()

        ! Set variables
        n_avg = p%n_conc(valley,x_idx,y_idx)
        Ed = p%E_avg(valley,x_idx,y_idx)

print *, x_idx, y_idx, n_avg, Ed
        func = fgsl_function_init(fermi_level, ptr)
        root_fslv = fgsl_root_fsolver_alloc(fgsl_root_fsolver_bisection)
        val_guess = -1000.0_s
        do i=1,100
        x_guess = 8.0_s/100.0_s*i-4.0_s
        val_guess = fermi_level(x_guess(i),ptr)
        !print *, i, x_guess(i), val_guess(i)
        if (val_guess(i) > 1e-2) then
            exit
        end if
        end do
        func = fgsl_function_init(fermi_level, ptr)
        root_fslv = fgsl_root_fsolver_alloc(fgsl_root_fsolver_brent)
        status = fgsl_root_fsolver_set(root_fslv, func, -4.0_fgsl_double, & !-3.0_fgsl_double, &
               x_guess(i))

        i = 0
        do
            i = i + 1
            status = fgsl_root_fsolver_iterate(root_fslv)
            if (status /= fgsl_success .or. i > itmax_root) then
                exit
            end if
            ra = fgsl_root_fsolver_root(root_fslv)
            xlo = fgsl_root_fsolver_x_lower(root_fslv)
            xhi = fgsl_root_fsolver_x_upper(root_fslv)
            status = fgsl_root_test_interval (xlo, xhi, 0.0_fgsl_double, eps5)
            if (status == fgsl_success) then
                exit
            else
                ra = -2.0_s
            end if
        end do
        !print *, ra, electron_temp(ra)
        call fgsl_root_fsolver_free(root_fslv)
        call fgsl_function_free(func)
        Ef = ra
        Te = electron_temp(ra)
        print *, Ef, Te
    end subroutine

    subroutine test_fd()
        integer(fgsl_int) :: status, i, valley, x_idx, y_idx
        integer,dimension(1)                   :: idx
        real(kind=s)                :: Ef, Te
        real(kind=s),dimension(1000) :: x_guess, val_guess
        real(fgsl_double), parameter :: eps5 = 1.0d-5
        integer(fgsl_int), parameter :: itmax_root = 1500
        type(c_ptr) :: ptr
        real(fgsl_double) :: ri, ra, xlo, xhi
        type(fgsl_function) :: func
        type(fgsl_root_fsolver) :: root_fslv
        type(fgsl_error_handler_t)    :: error_h

        error_h = fgsl_set_error_handler_off()
        ! Set variables
        n_avg = -1.3543333333333335e-3_s !p%n_conc(valley,x_idx,y_idx)
        Ed = 1.1528374236054982e-4_s !p%E_avg(valley,x_idx,y_idx)
!print *, x_idx, y_idx, n_avg, Ed
        ! initial guess
        p%mat = INGAAS
        p%valley = 1
        val_guess = -1000.0_s
        do i=1,100
        x_guess = 8.0_s/100.0_s*i-4.0_s
        val_guess = fermi_level(x_guess(i),ptr)
        print *, i, x_guess(i), val_guess(i)
        if (val_guess(i) > 1e-2) then
            exit
        end if
        end do
        func = fgsl_function_init(fermi_level, ptr)
        root_fslv = fgsl_root_fsolver_alloc(fgsl_root_fsolver_brent)
        status = fgsl_root_fsolver_set(root_fslv, func, -4.0_fgsl_double, & !-3.0_fgsl_double, &
               x_guess(i))

        i = 0
        do
            i = i + 1
            status = fgsl_root_fsolver_iterate(root_fslv)
            if (status /= fgsl_success .or. i > itmax_root) then
                exit
            end if
            ra = fgsl_root_fsolver_root(root_fslv)
            xlo = fgsl_root_fsolver_x_lower(root_fslv)
            xhi = fgsl_root_fsolver_x_upper(root_fslv)
            status = fgsl_root_test_interval (xlo, xhi, 0.0_fgsl_double, eps5)
            if (status == fgsl_success) then
                exit
            else
                print *, 'adf'
            end if
        end do
        !print *, ra, electron_temp(ra)
        call fgsl_root_fsolver_free(root_fslv)
        call fgsl_function_free(func)
        Ef = ra
        Te = electron_temp(ra)
        print *, Ef, Te
    end subroutine
#endif

    function fermi_level( x, params ) bind(c)
        implicit none
!        real(kind=fgsl_double),pointer :: arg
        real( kind = c_double )        :: fermi_level, Te
        real( kind = c_double ), value :: x
        type( c_ptr ), value           :: params
        integer(kind=fgsl_int)            :: i
        type(fgsl_error_handler_t)    :: error_h

        error_h = fgsl_set_error_handler_off()

     !   call c_f_pointer(params,arg)
        Te = electron_temp(x)
        fermi_level = int_n_conc(x,Te,i)-n_avg
        if (i .ne. 0) then
            x = -2.0_s
        end if
    end function fermi_level

    function energy_zero(Ef) result(E0)
        real(kind=s)        :: E0, Ef
        E0 = int_energy_zero(Ef)!/n_avg
    end function energy_zero

  function int_n_conc(Ef,temp,i) result(res)
        real( kind = fgsl_double )    :: res, Ef, temp
        real(kind=fgsl_double),dimension(2),target :: arg
        integer( kind = fgsl_int)     :: i
        type( fgsl_function )         :: func
        type(c_ptr)                   :: p_arg
        integer(fgsl_size_t), parameter  :: limit = 5000_fgsl_size_t
        type(fgsl_integration_workspace) :: integ_wk
        real(fgsl_double)                :: rda
        type(fgsl_error_handler_t)    :: error_h

        error_h = fgsl_set_error_handler_off()

        integ_wk = fgsl_integration_workspace_alloc(limit)
        arg(1) = Ef
        arg(2) = temp
        p_arg = c_loc(arg)
        func = fgsl_function_init(fermi_dos, p_arg)
        i = fgsl_integration_qagiu(func,0.0_fgsl_double,1.0d-5,1.0d-5,limit,integ_wk,res,rda)
        if (i .ne. 0) then
            print *, 'uups'
        end if
        call fgsl_function_free(func)
        call fgsl_integration_workspace_free(integ_wk)
    end function int_n_conc

  function int_energy_zero(Ef) result(res)
        real( kind = fgsl_double )    :: res, Ef
        real(kind=fgsl_double),target :: arg
        integer( kind = fgsl_int)     :: i
        type( fgsl_function )         :: func
        type(c_ptr)                   :: p_arg
        integer(fgsl_size_t), parameter  :: limit = 5000_fgsl_size_t
        type(fgsl_integration_workspace) :: integ_wk
        real(fgsl_double)                :: rda

        type(fgsl_error_handler_t)    :: error_h

        error_h = fgsl_set_error_handler_off()

        integ_wk = fgsl_integration_workspace_alloc(limit)
        arg = Ef
        p_arg = c_loc(arg)
        !print *, Ef
        func = fgsl_function_init(e_fermi_dos, p_arg)
        i = fgsl_integration_qagiu(func,0.0_fgsl_double,1.0d-5,1.0d-5,limit,integ_wk,res,rda)
        call fgsl_function_free(func)

        call fgsl_integration_workspace_free(integ_wk)
        !res = res
        !print *, res, res/n_avg
    end function int_energy_zero

    function e_fermi_dos( x, params ) bind(c)
        implicit none
        real(kind=fgsl_double),pointer :: arg
        real( kind = c_double )        :: e_fermi_dos
        real( kind = c_double ), value :: x
        type( c_ptr ), value           :: params

        call c_f_pointer(params,arg)
        !print *, arg, T_lattice, x, p%mat, p%valley
        !print *, x, dos(x,p%mat,p%valley)
        e_fermi_dos = x*fermi_dist(x,arg,T_lattice)*dos(x,p%mat,p%valley)
    end function e_fermi_dos

    function fermi_dos( x, params ) bind(c)
        implicit none
        real(kind=fgsl_double),dimension(:),pointer :: arg
        real( kind = c_double )        :: fermi_dos
        real( kind = c_double ), value :: x
        type( c_ptr ), value           :: params

        call c_f_pointer(params,arg,[2])
        !print *, x, arg
        fermi_dos = fermi_dist(x,arg(1),arg(2))*dos(x,p%mat,p%valley)
    end function fermi_dos

    function electron_temp(Ef) result(Te)
        real(kind=s)        :: Te, E0, Ef !Ed,

        !Ed = 0.29_s
        !Ed = p%Ed
        E0 = energy_zero(Ef)
        !print *, Ef, Ed, E0 !, E0*n_avg
        !print *, 2.0_s/(3.0_s*kB)*(Ed-E0), 2.0_s/(3.0_s*kB)*(Ed-E0)+T_lattice
        !if (Ed-E0 > 0) then
            Te = 2.0_s/(3.0_s*kB)*(Ed)!-E0)+T_lattice
        !else
        !    Te = T_lattice
        !end if
    end function electron_temp

    function fermi_dist(E,Ef,temp) result(fermi)
        real(kind=s)        :: fermi, E, Ef, temp

        fermi = 1.0_s/(exp((E-Ef)/(kB*temp))+1.0_s)
    end function fermi_dist

    function dos(energy,mat,valley) result(dd)
        real(kind=s)        :: dd, energy
        integer             :: mat, valley

        dd = (2.0_s*me(MAT,valley))**(3.0/2.0)/(2.0_s*pi**2.0_s*hb**3.0_s)* &
          & sqrt(energy*(1.0_s+nonp(MAT,valley)*energy))*(1.0_s+2.0_s*nonp(MAT,valley)*energy)
    end function dos
end module pauli
