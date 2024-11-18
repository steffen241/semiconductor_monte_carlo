module materialdef
    use physconst
    use config
    implicit none

    !integer,parameter                   :: GAAS_MODEL  = 0

    ! Herring-Vogt transformation matrix
    real(kind=s),dimension(NUM_MAT,3,4,3,3) :: T_hv, T_hv_inv
    real(kind=s),dimension(NUM_MAT,4,3,3)   :: T_hv_L, T_hv_L_inv, T_hv_L_t, T_hv_L_vel, D

    ! Material parameter
    real(kind=s),dimension(NUM_MAT,3)   :: me                         ! effective mass 1 Material, 2 Valley: G,L,X
    real(kind=s),dimension(NUM_MAT,3)   :: nonp                        ! nonparabolicity 1 Material, 2 Valley: G,L,X
    real(kind=s),dimension(NUM_MAT,3)   :: valley_offs                 ! energy offset of G, L and X-valley
    real(kind=s),dimension(NUM_MAT,6)   :: phonon_e                    ! Optical phonon energy: 1 Material Optical, 2 G-L, 3 G-X, 4 L-X, 5 L-L, 6 XX
    real(kind=s),dimension(NUM_MAT,6)   :: n_op                        ! Phonon occupation number


    real(kind=s),dimension(NUM_MAT,3)   :: def_pot_ac                  ! Acoustic deformation potential
    real(kind=s),dimension(NUM_MAT,3)   :: def_pot_nonp                ! Non-polar optical deformation potential
    real(kind=s),dimension(NUM_MAT)     :: def_pot_iv_GL, def_pot_iv_GX, def_pot_iv_LL, def_pot_iv_LX, def_pot_iv_XX  ! Intervalley deformation potential

    real(kind=s),dimension(NUM_MAT)     :: rho                         ! Material density
    real(kind=s),dimension(NUM_MAT)     :: sound_vel                   ! Longitudinal sound velocity

    real(kind=s),dimension(NUM_MAT)     :: eps_static, eps_hf          ! Static and high frequency relative permittivity
    real(kind=s),dimension(NUM_MAT)     :: F2                          ! Fröhlich coupling constant

    real(kind=s),dimension(NUM_MAT)     :: alloy_pot                   ! Alloy interaction potential
    real(kind=s),dimension(NUM_MAT)     :: lattice_const

    real(kind=s),dimension(NUM_MAT,3)   :: varshni_energy, varshni_alpha, varshni_beta, vall_offs_bow ! Parameter for used in the varshni model for temperature dependent valley offsets
    real(kind=s),dimension(NUM_MAT)     :: F, Ep, delta_so, me_G       ! Parameter for the computation of the new effective masses
    real(kind=s),dimension(NUM_MAT,3)   :: me_l, me_t
    real(kind=s),dimension(NUM_MAT)     :: F_bow, Ep_bow, delta_so_bow, me_bow
    real(kind=s),dimension(NUM_MAT)     :: mole_x, emin


    ! Array for scattering rates
    !real(kind=s),dimension(NUM_MAT,3,300000)     :: upper_bound = 0
    !real(kind=s),dimension(NUM_MAT,3,300000,15)  :: scat_rates = 0, n_scat_rates = 0
    real(kind=s),dimension(:,:,:),allocatable  :: upper_bound
    real(kind=s),dimension(:,:,:,:),allocatable  :: scat_rates, n_scat_rates
                                                                       ! scattering type: 1 - Self scattering
                                                                       !                  2 - Elastic ac phonon
                                                                       !                  3 - empty (later for inelastic ac phonon)
                                                                       !                  4 - Nonpolar phonon emission
                                                                       !                  5 - Nonpolar phonon absorption
                                                                       !                  6 - Polar phonon emission
                                                                       !                  7 - Polar phonon absorption
                                                                       ! Intervalley scattering now:
                                                                       !                  Gamma-valley
                                                                       !                  8 - Gamma->L emission
                                                                       !                  9 - Gamma->L absorption
                                                                       !                  10 - Gamma->X emission
                                                                       !                  11 - Gamma->X absorption
                                                                       !                  L-valley
                                                                       !                  8 - L->Gamma emission
                                                                       !                  9 - L->Gamma absorption
                                                                       !                  10 - L->L emission
                                                                       !                  11 - L->L absorption
                                                                       !                  12 - L->X emission
                                                                       !                  13 - L->X absorption
                                                                       !                  L-valley
                                                                       !                  8 - X->Gamma emission
                                                                       !                  9 - X->Gamma absorption
                                                                       !                  10 - X->L emission
                                                                       !                  11 - X->L absorption
                                                                       !                  12 - X->X emission
                                                                       !                  13 - X->X absorption
                                                                       !                  14 - Alloy scattering
                                                                       !                  15 - Impurity scattering

contains
    subroutine init_parameter
#if DIM > 0
        integer                                 :: diff_nodes
#endif
        real(kind=s),dimension(3,2)             :: me_tmp

        ! Allocate arrays for storing scattering rates
#if DIM == 0
        allocate(scat_rates(1,3,300000,15))
        allocate(n_scat_rates(1,3,300000,15))
        allocate(upper_bound(1,3,300000))
#elif DIM > 0
        !diff_nodes = size(m(1,:))
        diff_nodes = maxval(scat_m)
        !print *, diff_nodes
        allocate(scat_rates(diff_nodes,3,300000,15))
        allocate(n_scat_rates(diff_nodes,3,300000,15))
        allocate(upper_bound(diff_nodes,3,300000))
        !allocate(scat_rates(NUM_MAT,3,40000,15,diff_nodes))
        !allocate(n_scat_rates(NUM_MAT,3,40000,15,diff_nodes))
        !allocate(upper_bound(NUM_MAT,3,40000,diff_nodes))
#endif
        scat_rates = 0.0_s
        n_scat_rates = 0.0_s
        upper_bound = 0.0_s

        !
        ! Air
        !
        emin(AIR)             = 10.0_s
        valley_offs(AIR,:)    = (/0.0_s,0.0_s,0.0_s/)
        eps_static(AIR)      = 1.0_s*eps_vac !12.9_s*eps_vac !12.90_s*eps_vac
        eps_hf(AIR)          = 1.0_s*eps_vac
        !
        ! GaAs
        !
        ! - band structure parameter
        me(GAAS,:)           = [0.06223_s*m0,0.292_s*m0,0.471_s*m0]
        nonp(GAAS,:)         = [1.16_s,0.51_s,0.58_s] ![1.16_s,2.24_s,0.58_s] ![1.16_s,2.24_s,0.4_s] !0.85_s] ![1.16_s,1.5_s,1.1_s] ![1.16_s,2.81_s,0.58_s] ! !
        valley_offs(GAAS,:)  = [0.0_s,0.323_s,0.447_s]
        !print *, me(GAAS,L)

        ! - scattering rate parameters
        def_pot_ac(GAAS,1)    = 7.0_s !5.0_s
        def_pot_ac(GAAS,2)    = 7.0_s !7.0_s
        rho(GAAS)             = 5317.0_s*1e24_s*1e-18_s/q/1e27_s/q_si ! 5360
        sound_vel(GAAS)       = 5400.0_s*1e9_s/1e12_s
        phonon_e(GAAS,:)      = [35.36e-3_s,22e-3_s,20e-3_s,23e-3_s,19e-3_s,22e-3_s] ![35.36e-3_s, 20.0e-3_s,20.0e-3_s,20.0e-3_s,20.0e-3_s,20.0e-3_s] [35.36e-3_s, 22.69e-3_s,23.45e-3_s,21.85e-3_s,24.97e-3_s,24.31e-3_s] ![35.36e-3_s, 22.69e-3_s,23.45e-3_s,21.85e-3_s,24.97e-3_s,24.31e-3_s]
        n_op(GAAS,:)          = 1.0_s/(exp(phonon_e(GAAS,:)/(kB*T_lattice))-1.0_s)
        def_pot_iv_GL(GAAS)   = 55.0_s !55.0_s !52.5_s !53.0_s !
        def_pot_iv_GX(GAAS)   = 54.8_s !55.0_s !
        def_pot_iv_LL(GAAS)   = 59.4_s !60.0_s !
        def_pot_iv_LX(GAAS)   = 50.1_s ! 50.0_s !
        def_pot_iv_XX(GAAS)   = 29.9_s !30.0_s !
        def_pot_nonp(GAAS,1)  = 0.0_s;
        def_pot_nonp(GAAS,2)  = 21.0_s;
        eps_static(GAAS)      = 12.9_s*eps_vac !12.9_s*eps_vac !12.90_s*eps_vac
        eps_hf(GAAS)          = 10.9_s*eps_vac !10.8_s*eps_vac !10.8_s*eps_vac !
        F2(GAAS)              = phonon_e(GAAS,1)/(4.0_s)*(1.0_s/eps_hf(GAAS)-1.0_s/eps_static(GAAS))
        lattice_const(GAAS)   = 0.565325_s

        if (GAAS_MODEL .eq. 1) then
        ! - temperature dependent parameter
        ! - valley offset: Varshni model
        varshni_energy(GAAS,:)= [1.519_s,1.815_s,1.981_s]
        varshni_alpha(GAAS,:) = [0.5405e-3_s,0.46e-3_s,0.605e-3_s]
        varshni_beta(GAAS,:)  = [204.0_s,204.0_s,204.0_s]
        ! - parameters for effective masses
        F(GAAS)               = -1.94_s
        Ep(GAAS)              = 28.8_s
        delta_so(GAAS)        = 0.341_s
        !me_G(GAAS)            = 0.067_s
        me_l(GAAS,L)          = 1.9_s
        me_t(GAAS,L)          = 0.0754_s
        me_l(GAAS,X)          = 1.3_s
        me_t(GAAS,X)          = 0.23_s

        ! set the corresponding variables to the new parameter
        valley_offs(GAAS,:)   = varshni_model(GAAS, T_lattice)
        me(GAAS,G)            = effective_mass_model(GAAS, T_lattice)*m0
        me(GAAS,L)            = (me_l(GAAS,L)*me_t(GAAS,L)**2.0_s)**(1.0_s/3.0_s)*m0
        !me(GAAS,L)  = 1.25*m0
        !print *, me(GAAS,L), 1.25*m0
        me(GAAS,X)            = (me_l(GAAS,X)*me_t(GAAS,X)**2.0_s)**(1.0_s/3.0_s)*m0
        nonp(GAAS,G)          = 1.0_s/valley_offs(GAAS,G)*(1.0_s-me(GAAS,G)/m0)**2.0_s
        !Eg                    = valley_offs(GAAS,L)
        !nonp(GAAS,L)          = 1.0_s/Eg*(1.0_s-me(GAAS,G)/m0)**2.0_s*(1- &
        !                      & (Eg*delta_so(GAAS))/(3.0_s*(Eg+delta_so(GAAS))*(Eg+2.0_s*delta_so(GAAS)/3.0_s)))
        !print *, 'alpha L:', nonp(GAAS,L)
        !Eg                    = valley_offs(GAAS,X)
        !nonp(GAAS,X)          = 1.0_s/Eg*(1.0_s-me(GAAS,X)/m0)**2.0_s*(1- &
        !                      & (Eg*delta_so(GAAS))/(3.0_s*(Eg+delta_so(GAAS))*(Eg+2.0_s*delta_so(GAAS)/3.0_s)))
        !print *, 'alpha X:', nonp(GAAS,X)
        valley_offs(GAAS,:)   = valley_offs(GAAS,:)-valley_offs(GAAS,1)

        call init_herring_vogt(GAAS)
        end if
        !
        ! InAs
        !
        ! - temperature dependent parameter
        ! - valley offset: Varshni model
        varshni_energy(INAS,:)= [0.417_s,1.133_s,1.433_s]
        varshni_alpha(INAS,:) = [0.276e-3_s,0.276e-3_s,0.276e-3_s]
        varshni_beta(INAS,:)  = [93.0_s,93.0_s,93.0_s]
        F(INAS)               = -2.9_s
        Ep(INAS)              = 21.5_s
        delta_so(INAS)        = 0.39_s
        me_l(INAS,L)          = 0.64_s
        me_t(INAS,L)          = 0.05_s
        me_l(INAS,X)          = 1.13_s
        me_t(INAS,X)          = 0.16_s
        lattice_const(INAS)   = 0.60583_s
        sound_vel(INAS)       = 4410*1e9_s/1e12_s
        rho(INAS)             = 5680*1e24_s*1e-18_s/q/1e27_s/q_si
        eps_static(INAS)      = 15.15_s*eps_vac !12.9_s*eps_vac !12.90_s*eps_vac
        eps_hf(INAS)          = 12.25_s*eps_vac

        !
        ! AlAs
        !
        ! - temperature dependent parameter
        ! - valley offset: Varshni model
        varshni_energy(ALAS,:)= [3.099_s, 2.46_s, 2.24_s]
        varshni_alpha(ALAS,:) = [0.885e-3_s,0.605e-3_s,0.885e-3_s]
        varshni_beta(ALAS,:)  = [530.0_s, 204.0_s, 530.0_s]
        F(ALAS)               = -0.48_s
        Ep(ALAS)              = 21.1_s
        delta_so(ALAS)        = 0.28_s
        me_l(ALAS,L)          = 1.32_s
        me_t(ALAS,L)          = 0.15_s
        me_l(ALAS,X)          = 0.97_s
        me_t(ALAS,X)          = 0.22_s
        lattice_const(ALAS)   = 0.56611_s
        sound_vel(ALAS)       = 6480*1e9_s/1e12_s
        rho(ALAS)             = 3730*1e24_s*1e-18_s/q/1e27_s/q_si
        eps_static(ALAS)      = 9.4625*eps_vac
        eps_hf(ALAS)          = 8.1458*eps_vac

        !
        ! InGaAs
        !
        mole_x(INGAAS)        = dble(0.53) !dble(0.53) !mole_x_tmp
        !
        ! - band structure parameter
        emin(INGAAS)          = 0.0_s
        me(INGAAS,:)          = [0.0463_s*m0, 0.256_s*m0, 0.529_s*m0]
        nonp(INGAAS,:)        = [1.18, 0.588, 0.649] ![1.18_s, 0.588_s, 0.649_s]
        valley_offs(INGAAS,:) = [0.0_s, 0.55_s, 0.67_s]

        ! - scattering rate parameters, Indium mole fraction x=0.53
        def_pot_ac(INGAAS,1)  = 9.2_s
        rho(INGAAS)           = 5545_s*1e24_s*1e-18_s/q/1e27_s/q_si !5480
        sound_vel(INGAAS)     = 4756_s*1e9_s/1e12_s !4550
        phonon_e(INGAAS,:)    = [32.7e-3_s, 22.76e-3_s, 23.84e-3_s, 23.12e-3_s, 26.96e-3_s, 22.76e-3_s] ! ![32.8e-3_s, 25.4e-3_s, 25.7e-3_s, 30.2e-3_s, 24.8e-3_s, 28.4e-3_s]![32.7e-3_s, 22.76e-3_s, 23.84e-3_s, 23.12e-3_s, 26.96e-3_s, 22.76e-3_s] ![32.7e-3_s, 30.0e-3_s, 29.9e-3_s, 29.0e-3_s, 29.0e-3_s, 29.9e-3_s] !
        n_op(INGAAS,:)        = 1.0_s/(exp(phonon_e(INGAAS,:)/(kB*T_lattice))-1.0_s)
        def_pot_iv_GL(INGAAS) = 70.0_s !100.0_s !78.3_s  !70.0_s
        def_pot_iv_GX(INGAAS) = 70.0_s !100.0_s !113.2_s !70.0_s
        def_pot_iv_LL(INGAAS) = 70.0_s !100.0_s !64.0_s  !70.0_s
        def_pot_iv_LX(INGAAS) = 50.0_s !90.0_s !68.0_s  !50.0_s
        def_pot_iv_XX(INGAAS) = 58.0_s !90.0_s !85.4_s  !58.0_s
        def_pot_nonp(INGAAS,1)= 0.0_s;
        def_pot_nonp(INGAAS,2)= 30.0_s;
        eps_static(INGAAS)    = 14.092_s*eps_vac !13.85_s*eps_vac
        eps_hf(INGAAS)        = 11.6155_s*eps_vac !11.09_s*eps_vac
        F2(INGAAS)            = phonon_e(INGAAS,1)/(4.0_s)*(1.0_s/eps_hf(INGAAS)-1.0_s/eps_static(INGAAS))
        alloy_pot(INGAAS)     = 0.53 ! 0.6_s
        lattice_const(INGAAS) = 0.5867_s

        if (INGAAS_MODEL .eq. 1) then
        ! - temperature dependent parameter
        ! - valley offset: Varshni model; interpolation with bowing parameter
        vall_offs_bow(INGAAS,:) = [0.477_s, 0.33_s, 0.85_s]!0.8_s]!0.33_s, 0.08_s]!0.08_s]  0.5_s, 0.08_s] !
        ! - parameters for effective masses
        F_bow(INGAAS)         = 1.77_s
        Ep_bow(INGAAS)        = -1.48_s
        delta_so_bow(INGAAS)  = 0.15_s
        me_bow(INGAAS)        = 0.0091
        me_tmp                = effective_mass_model_alloys_LX(GAAS,INAS,mole_x(INGAAS))
        me_l(INGAAS,L)        = me_tmp(L,1)
        me_t(INGAAS,L)        = me_tmp(L,2)
        me_l(INGAAS,X)        = me_tmp(X,1)
        me_t(INGAAS,X)        = me_tmp(X,2)

        sound_vel(INGAAS)     = mole_x(INGAAS)*sound_vel(INAS)+(1.0_s-mole_x(INGAAS))*sound_vel(GAAS)
        lattice_const(INGAAS) = mole_x(INGAAS)*lattice_const(INAS)+(1.0_s-mole_x(INGAAS))*lattice_const(GAAS)
        rho(INGAAS)           = mole_x(INGAAS)*rho(INAS)+(1.0_s-mole_x(INGAAS))*rho(GAAS)
        eps_hf(INGAAS)        = mole_x(INGAAS)*eps_hf(INAS)+(1.0_s-mole_x(INGAAS))*eps_hf(GAAS)
        eps_static(INGAAS)    = mole_x(INGAAS)*eps_static(INAS)+(1.0_s-mole_x(INGAAS))*eps_static(GAAS)
        !print *, 'eps', eps_hf(INGAAS)/eps_vac, eps_static(INGAAS)/eps_vac

        ! set the corresponding variables to the new parameter
        valley_offs(INGAAS,:) = valley_offs_alloys(INGAAS,T_lattice,GAAS,INAS,mole_x(INGAAS))
        me(INGAAS,G)         = effective_mass_model_alloys(INGAAS,T_lattice,GAAS,INAS,mole_x(INGAAS))*m0
        me(INGAAS,L)         = (me_l(INGAAS,L)*me_t(INGAAS,L)**2.0_s)**(1.0_s/3.0_s)*m0
        me(INGAAS,X)         = (me_l(INGAAS,X)*me_t(INGAAS,X)**2.0_s)**(1.0_s/3.0_s)*m0
        !nonp(INGAAS,G)        = 1.0_s/valley_offs(INGAAS,G)*(1.0_s-me(INGAAS,G)/m0)**2.0_s
        print *, 'ingaas valley offs', valley_offs(INGAAS,:)
        valley_offs(INGAAS,:) = valley_offs(INGAAS,:)-valley_offs(INGAAS,1)

        ! set deformation potentials and phonon energy, according to mole fraction
        nonp(INGAAS,:) = abs(nonp(INGAAS,:)-nonp(GAAS,:))/0.53_s*mole_x(INGAAS)+nonp(GAAS,:)
        def_pot_ac(INGAAS,:) = abs(def_pot_ac(INGAAS,:)-def_pot_ac(GAAS,:))/0.53_s*mole_x(INGAAS)+def_pot_ac(GAAS,:)
        phonon_e(INGAAS,:) = abs(phonon_e(INGAAS,:)-phonon_e(GAAS,:))/0.53_s*mole_x(INGAAS)+phonon_e(GAAS,:)
        n_op(INGAAS,:) = 1.0_s/(exp(phonon_e(INGAAS,:)/(kB*T_lattice))-1.0_s)
        F2(INGAAS) = phonon_e(INGAAS,1)/(4.0_s)*(1.0_s/eps_hf(INGAAS)-1.0_s/eps_static(INGAAS))
        def_pot_iv_GL(INGAAS) = abs(def_pot_iv_GL(INGAAS)-def_pot_iv_GL(GAAS))/0.53_s*mole_x(INGAAS)+def_pot_iv_GL(GAAS)
        def_pot_iv_GX(INGAAS) = abs(def_pot_iv_GX(INGAAS)-def_pot_iv_GX(GAAS))/0.53_s*mole_x(INGAAS)+def_pot_iv_GX(GAAS)
        def_pot_iv_LL(INGAAS) = abs(def_pot_iv_LL(INGAAS)-def_pot_iv_LL(GAAS))/0.53_s*mole_x(INGAAS)+def_pot_iv_LL(GAAS)
        def_pot_iv_LX(INGAAS) = abs(def_pot_iv_LX(INGAAS)-def_pot_iv_LX(GAAS))/0.53_s*mole_x(INGAAS)+def_pot_iv_LX(GAAS)
        def_pot_iv_XX(INGAAS) = abs(def_pot_iv_XX(INGAAS)-def_pot_iv_XX(GAAS))/0.53_s*mole_x(INGAAS)+def_pot_iv_XX(GAAS)
        def_pot_nonp(INGAAS,:)= abs(def_pot_nonp(INGAAS,:)-def_pot_nonp(GAAS,:))/0.53_s*mole_x(INGAAS)+def_pot_nonp(GAAS,:)
        call init_herring_vogt(INGAAS)
        end if

        !
        ! InAlAs
        !
        emin(INALAS)          = 0.52_s
        mole_x(INALAS)        = dble(0.52_s)
        !
        ! - band structure parameter
        me(INALAS,:)         = [0.07_s*m0, 0.39_s*m0, 0.602_s*m0]
        nonp(INALAS,:)        = [0.571_s, 0.204_s, 0.204_s] ![0.843_s, 0.220_s, 0.17_s]!0.552_s, 0.204_s] !0.588_s]
        valley_offs(INALAS,:) = [0.0_s, 0.34_s, 0.6_s]

        ! - scattering rate parameters
        def_pot_ac(INALAS,1)  = 8.0_s
        !rho(INALAS)           = 4900_s*1e24_s*1e-18_s/q/1e27_s/q_si
        !sound_vel(INALAS)     = 4970_s*1e9_s/1e12_s
        phonon_e(INALAS,:)    = [41.0e-3_s, 29.0e-3_s, 29.0e-3_s, 29.0e-3_s, 29.0e-3_s, 29.0e-3_s] !, 30.8e-3_s, 29.10e-3_s, 36.0e-3_s, 31.1e-3_s, 29.1e-3_s]
        n_op(INALAS,:)        = 1.0_s/(exp(phonon_e(INALAS,:)/(kB*T_lattice))-1.0_s)
        def_pot_iv_GL(INALAS) = 100.0_s !50.0_s
        def_pot_iv_GX(INALAS) = 100.0_s!74.0_s
        def_pot_iv_LL(INALAS) = 100.0_s!50.0_s
        def_pot_iv_LX(INALAS) = 100.0_s!50.0_s
        def_pot_iv_XX(INALAS) = 100.0_s!74.0_s
        def_pot_nonp(INALAS,1)= 0.0_s;
        def_pot_nonp(INALAS,2)= 30.0_s;
        !eps_static(INALAS)    = 12.46_s*eps_vac
        !eps_hf(INALAS)        = 9.84_s*eps_vac

        sound_vel(INALAS)     = mole_x(INALAS)*sound_vel(INAS)+(1.0_s-mole_x(INALAS))*sound_vel(ALAS)
        lattice_const(INALAS) = mole_x(INALAS)*lattice_const(INAS)+(1.0_s-mole_x(INALAS))*lattice_const(ALAS)
        rho(INALAS)           = mole_x(INALAS)*rho(INAS)+(1.0_s-mole_x(INALAS))*rho(ALAS)
        eps_hf(INALAS)        = mole_x(INALAS)*eps_hf(INAS)+(1.0_s-mole_x(INALAS))*eps_hf(ALAS)
        eps_static(INALAS)    = mole_x(INALAS)*eps_static(INAS)+(1.0_s-mole_x(INALAS))*eps_static(ALAS)
        F2(INALAS)            = phonon_e(INALAS,1)/(4.0_s)*(1.0_s/eps_hf(INALAS)-1.0_s/eps_static(INALAS))
        alloy_pot(INALAS)     = 0.47_s
        !lattice_const(INALAS) = 0.5867_s

        if (INALAS_MODEL .eq. 1) then
        ! - temperature dependent parameter
        ! - valley offset: Varshni model; interpolation with bowing parameter
        vall_offs_bow(INALAS,:) = [0.7_s, 0.0_s, 0.0_s]
        ! - parameters for effective masses
        F_bow(INALAS)         = -4.44_s
        Ep_bow(INALAS)        = -4.81_s
        delta_so_bow(INALAS)  = 0.15_s
        me_bow(INALAS)        = 0.049
        me_tmp                = effective_mass_model_alloys_LX(ALAS,INAS,mole_x(INALAS))
        me_l(INALAS,L)        = me_tmp(L,1)
        me_t(INALAS,L)        = me_tmp(L,2)
        me_l(INALAS,X)        = me_tmp(X,1)
        me_t(INALAS,X)        = me_tmp(X,2)

        valley_offs(INALAS,:) = valley_offs_alloys(INALAS,T_lattice,ALAS,INAS,mole_x(INALAS))
        me(INALAS,G)         = effective_mass_model_alloys(INALAS,T_lattice,ALAS,INAS,mole_x(INALAS))*m0
        me(INALAS,L)         = (me_l(INALAS,L)*me_t(INALAS,L)**2.0_s)**(1.0_s/3.0_s)*m0
        me(INALAS,X)         = (me_l(INALAS,X)*me_t(INALAS,X)**2.0_s)**(1.0_s/3.0_s)*m0
        nonp(INALAS,G)        = 1.0_s/valley_offs(INALAS,G)*(1.0_s-me(INALAS,G)/m0)**2.0_s
        print *, 'inalas valley offs', valley_offs(INALAS,:)
        valley_offs(INALAS,:) = valley_offs(INALAS,:)-valley_offs(INALAS,1)
        call init_herring_vogt(INALAS)

        !print *, me(INALAS,:)/m0, nonp(INALAS,:)
        !print *, valley_offs(INALAS,:)
        !print *, valley_offs(INALAS,:),valley_offs(INALAS,:)-valley_offs(INALAS,1)
        end if
    end subroutine init_parameter

    subroutine init_herring_vogt(MAT)
        integer                        :: MAT
        real(kind=s)                   :: e_l, e_t, e_Ll, e_Lt, sq3, sq2, sq6
        real(kind=s),dimension(3,3)    :: T_hv_x, T_hv_x_inv, tmp, tmp_inv

        sq2 = 1.0_s/sqrt(2.0_s)
        sq3 = 1.0_s/sqrt(3.0_s)
        sq6 = 1.0_s/sqrt(6.0_s)
        e_l = sqrt(1.0_s/me_l(MAT,X))
        e_t = sqrt(1.0_s/me_t(MAT,X))
        e_Ll = sqrt(1.0_s/me_l(MAT,L))
        e_Lt = sqrt(1.0_s/me_t(MAT,L))
        T_hv_x = reshape((/e_Ll,0.0_s,0.0_s,0.0_s,e_Lt,0.0_s,0.0_s,0.0_s,e_Lt/), (/3,3/))
        T_hv_x_inv = reshape((/1.0_s/e_Ll,0.0_s,0.0_s,0.0_s,1.0_s/e_Lt,0.0_s,0.0_s,0.0_s,1.0_s/e_Lt/), (/3,3/))
        ! Matrix for L-valleys


        ! [111] and [_1_1_1] direction
        tmp = reshape((/sq3,-sq2,-sq6,sq3,sq2,-sq6,sq3,0.0_s,sqrt(2.0_s/3.0_s)/),(/3,3/))
        D(MAT,1,:,:) = tmp
        !tmp = reshape((/sq3,sq2,-sq6,sq3,sq2,-sq6,sq3,0.0_s,sqrt(2.0_s/3.0_s)/),(/3,3/))
        tmp_inv = reshape((/sq3,-sq2,-0.408248290463863_s,sq3,sq2,-0.408248290463863_s,sq3,0.0_s,sqrt(2.0_s/3.0_s)/), (/3,3/))
        !tmp = transpose(tmp)
        !print *, matmul(tmp,(/1,1,1/))
        T_hv_L(MAT,1,:,:) = (matmul(T_hv_x,tmp)) !reshape(T_hv_x, (/3,3/)) !matmul(T_hv_x,tmp) !reshape(T_hv_x, (/3,3/)) !matmul(T_hv_x,tmp)
        T_hv_L_inv(MAT,1,:,:) = transpose(matmul(T_hv_x_inv,tmp_inv))
        T_hv_L_vel(MAT,1,:,:) = matmul(T_hv_x_inv,tmp)
        !print *, T_hv_L(MAT,1,1,:)
        !print *, T_hv_L(MAT,1,2,:)
        !print *, T_hv_L(MAT,1,3,:)

        ! [_111] and [1_1_1] direction
        tmp = reshape((/-sq3,-sq2,sq6,sq3,-sq2,-sq6,sq3,0.0_s,sqrt(2.0_s/3.0_s)/),(/3,3/))
        D(MAT,2,:,:) = tmp
        T_hv_L(MAT,2,:,:) = (matmul(T_hv_x,tmp)) !reshape(T_hv_x, (/3,3/)) !matmul(T_hv_x,tmp)
        tmp_inv = reshape((/-sq3,-sq2,0.408248290463863_s,sq3,-sq2,-0.408248290463863_s,sq3,0.0_s,sqrt(2.0_s/3.0_s)/), (/3,3/))
        T_hv_L_inv(MAT,2,:,:) = transpose(matmul(T_hv_x_inv,tmp_inv))
        T_hv_L_vel(MAT,2,:,:) = matmul(T_hv_x_inv,tmp)

        ! [11_1] and [_1_11] direction
        tmp = reshape((/-sq3,sq2,sq6,-sq3,-sq2,sq6,sq3,0.0_s,sqrt(2.0_s/3.0_s)/),(/3,3/))
        D(MAT,3,:,:) = tmp
        T_hv_L(MAT,3,:,:) = (matmul(T_hv_x,tmp)) !reshape(T_hv_x, (/3,3/)) !matmul(T_hv_x,tmp)
        tmp_inv = reshape((/-sq3,sq2,0.408248290463863_s,-sq3,-sq2,0.408248290463863_s,sq3,0.0_s,sqrt(2.0_s/3.0_s)/), (/3,3/))
        T_hv_L_inv(MAT,3,:,:) = transpose(matmul(T_hv_x_inv,tmp_inv))
        T_hv_L_vel(MAT,3,:,:) = matmul(T_hv_x_inv,tmp)

        ! [1_1_1] and [_11_1] direction
        tmp = reshape((/sq3,sq2,-sq6,-sq3,sq2,sq6,sq3,0.0_s,sqrt(2.0_s/3.0_s)/),(/3,3/))
        D(MAT,4,:,:) = tmp
        T_hv_L(MAT,4,:,:) = (matmul(T_hv_x,tmp)) !reshape(T_hv_x, (/3,3/)) !matmul(T_hv_x,tmp)
        tmp_inv = reshape((/sq3,sq2,-0.408248290463863_s,-sq3,sq2,0.408248290463863_s,sq3,0.0_s,sqrt(2.0_s/3.0_s)/), (/3,3/))
        T_hv_L_inv(MAT,4,:,:) = (matmul(T_hv_x_inv,tmp_inv))
        T_hv_L_vel(MAT,4,:,:) = matmul(T_hv_x_inv,tmp)

        T_hv_L_t(MAT,1,:,:) = transpose(T_hv_L(MAT,1,:,:))
        T_hv_L_t(MAT,2,:,:) = transpose(T_hv_L(MAT,2,:,:))
        T_hv_L_t(MAT,3,:,:) = transpose(T_hv_L(MAT,3,:,:))
        T_hv_L_t(MAT,4,:,:) = transpose(T_hv_L(MAT,4,:,:))

        ! Matrix for x-direction in X-valleys
        T_hv(MAT,X,1,:,:) = reshape((/e_l,0.0_s,0.0_s,0.0_s,e_t,0.0_s,0.0_s,0.0_s,e_t/), (/3,3/))
        T_hv_inv(MAT,X,1,:,:) = reshape((/1.0_s/e_l,0.0_s,0.0_s,0.0_s,1.0_s/e_t,0.0_s,0.0_s,0.0_s,1.0_s/e_t/), (/3,3/))
        ! X-valley, y-direction
        T_hv(MAT,X,2,:,:) = reshape((/e_t,0.0_s,0.0_s,0.0_s,e_l,0.0_s,0.0_s,0.0_s,e_t/), (/3,3/))
        T_hv_inv(MAT,X,2,:,:) = reshape((/1.0_s/e_t,0.0_s,0.0_s,0.0_s,1.0_s/e_l,0.0_s,0.0_s,0.0_s,1.0_s/e_t/), (/3,3/))
        ! X-valley, z-direction
        T_hv(MAT,X,3,:,:) = reshape((/e_t,0.0_s,0.0_s,0.0_s,e_t,0.0_s,0.0_s,0.0_s,e_l/), (/3,3/))
        T_hv_inv(MAT,X,3,:,:) = reshape((/1.0_s/e_t,0.0_s,0.0_s,0.0_s,1.0_s/e_t,0.0_s,0.0_s,0.0_s,1.0_s/e_l/), (/3,3/))
    end subroutine init_herring_vogt

    function varshni_model(MAT,temp) result(varshni_offs)
        integer                        :: MAT
        real(kind=s)                   :: temp
        real(kind=s),dimension(3)      :: varshni_offs

        varshni_offs = varshni_energy(MAT,:)-varshni_alpha(MAT,:)*temp**2.0_s/(temp+varshni_beta(MAT,:))
    end function

    function valley_offs_alloys(MAT,temp,bin1,bin2,mole_x) result(vall_offs)
        integer                        :: MAT, bin1, bin2
        real(kind=s)                   :: temp, mole_x
        real(kind=s),dimension(3)      :: vall_offs, Eg_bin1, Eg_bin2 !, varshni_alpha_new, varshni_beta_new

        Eg_bin1 = varshni_model(bin1,temp)
        Eg_bin2 = varshni_model(bin2,temp)

        vall_offs = (1.0_s-mole_x)*Eg_bin1+mole_x*Eg_bin2-mole_x*(1.0_s-mole_x)*vall_offs_bow(MAT,:)

        !Eg_bin1 = varshni_energy(bin1,:)
        !Eg_bin2 = varshni_energy(bin2,:)

        !varshni_alpha(MAT,:) = (1.0_s-mole_x)*varshni_alpha(bin1,:)+mole_x*varshni_alpha(bin2,:)
        !varshni_beta(MAT,:) = (1.0_s-mole_x)*varshni_beta(bin1,:)+mole_x*varshni_beta(bin2,:)

        !varshni_energy(MAT,:) = (1.0_s-mole_x)*Eg_bin1+mole_x*Eg_bin2-mole_x*(1.0_s-mole_x)*vall_offs_bow(MAT,:)
        !vall_offs = varshni_model(MAT,temp)
    end function

    function effective_mass_model(MAT,temp) result(me_G)
        integer                        :: MAT
        real(kind=s)                   :: temp
        real(kind=s)                   :: me_G
        real(kind=s),dimension(3)      :: Eg

        Eg = varshni_model(MAT,temp)
        me_G = 1/((1.0_s+2.0_s*F(MAT))+Ep(MAT)*(Eg(1)+2.0_s/3.0_s*delta_so(MAT))/(Eg(1)*(Eg(1)+delta_so(MAT))))
        print *, MAT, Eg(1), me_G
    end function effective_mass_model

    function effective_mass_model_alloys(MAT,temp,bin1,bin2,mole_x) result(me_G)
        integer                        :: MAT, bin1, bin2
        real(kind=s)                   :: temp, me_G, mole_x, me_bin1, me_bin2
        !real(kind=s),dimension(3)      :: Eg

        !F(MAT) = (1.0_s-mole_x)*F(bin1)+mole_x*F(bin2)-mole_x*(1.0_s-mole_x)*F_bow(MAT)
        !Ep(MAT) = (1.0_s-mole_x)*Ep(bin1)+mole_x*Ep(bin2)-mole_x*(1.0_s-mole_x)*Ep_bow(MAT)
        !delta_so(MAT) = (1.0_s-mole_x)*delta_so(bin1)+mole_x*delta_so(bin2)-mole_x*(1.0_s-mole_x)*delta_so_bow(MAT)
        !Eg = valley_offs_alloys(MAT,temp,bin1,bin2,dble(0.53))
        !me_G = 1.0_s/((1.0_s+2.0_s*F(MAT))+Ep(MAT)*(Eg(1)+2.0_s/3.0_s*delta_so(MAT))/(Eg(1)*(Eg(1)+delta_so(MAT))))
        !print *, MAT, me_G
        me_bin1 = effective_mass_model(bin1,temp)
        me_bin2 = effective_mass_model(bin2,temp)
        me_g = (1.0_s-mole_x)*me_bin1+mole_x*me_bin2-mole_x*(1.0_s-mole_x)*me_bow(MAT)
        !print *, MAT, me_G
        !print *, me_G, F(INGAAS), Ep(INGAAS), delta_so(INGAAS)
    end function effective_mass_model_alloys

    function effective_mass_model_alloys_LX(bin1,bin2,mole_x) result(me)
        integer                         :: bin1, bin2
        real(kind=s)                    :: mole_x
        real(kind=s),dimension(3,2)     :: me               ! 1: longitudinal effective mass; 2: transversal effective mass

        me(L,1) = (1.0_s-mole_x)*me_l(bin1,L)+mole_x*me_l(bin2,L)
        me(X,1) = (1.0_s-mole_x)*me_l(bin1,X)+mole_x*me_l(bin2,X)
        me(L,2) = (1.0_s-mole_x)*me_t(bin1,L)+mole_x*me_t(bin2,L)
        me(X,2) = (1.0_s-mole_x)*me_t(bin1,X)+mole_x*me_t(bin2,X)
    end function effective_mass_model_alloys_LX
end module materialdef
