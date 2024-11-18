module scattering_rates
    use physconst
    use materialdef
    use config
    use gnufor2
    use integrals
    implicit none

    contains

    subroutine calc_scatrates
        !call calc_scatrates_GAAS()
        call calc_scatrates_INGAAS()
        call calc_scatrates_INALAS()
    end subroutine calc_scatrates

!    subroutine calc_scatrates_GAAS()
!        integer                              :: i
!        real(kind=s),dimension(3,size_E)     :: ph_ac, ph_ac_tmp
!        real(kind=s),dimension(2,size_E)     :: ph_tmp, ph_iv_GL, ph_iv_GX, &
!                                              & ph_iv_LG, ph_iv_LL, ph_iv_LX, &
!                                              & ph_iv_XG, ph_iv_XL, ph_iv_XX
!        real(kind=s),dimension(3,2,size_E)   :: ph_nonp, ph_pol, ph_in_ac
!        real(kind=s),dimension(:,:,:),allocatable :: self_scat
!
!        allocate(self_scat(3,size(m(1,:)),size_E))
!        ph_ac = 0;
!        ph_ac_tmp = 0;
!
!        ! Gamma valley
!        ph_ac_tmp(G,:) = ph_ac_scatrate(GAAS,1,1)
!        ph_ac(G,1:4000) = ph_ac_tmp(GAAS,1:4000)
!        ph_ac_tmp(G,:) = ph_ac_scatrate(GAAS,1,2)
!        ph_ac(G,4001:size_E) = ph_ac_tmp(GAAS,4001:size_E)
!        ph_tmp = ph_nonp_scatrate(GAAS,G,1)
!        ph_nonp(G,1:2,1:4000) = ph_tmp(1:2,1:4000)
!        ph_tmp = ph_nonp_scatrate(GAAS,G,2)
!        ph_nonp(G,1:2,4001:size_E) = ph_tmp(1:2,4001:size_E)
!        ph_iv_GL = ph_iv_scatrate(GAAS,G,L)
!        ph_iv_GX = ph_iv_scatrate(GAAS,G,X)
!        ph_pol(G,1:2,:) = ph_pol_scatrate(GAAS,G)
!        ph_in_ac(G,1:2,:) = ph_inelastic_ac(GAAS,G)
!
!        ! L valley
!        ph_ac(L,:) = ph_ac_scatrate(GAAS,L,2)
!        ph_nonp(L,1:2,:) = ph_nonp_scatrate(GAAS,L,2)
!        ph_iv_LG(1:2,:) = ph_iv_scatrate(GAAS,L,G)
!        ph_iv_LL(1:2,:) = ph_iv_scatrate(GAAS,L,L)
!        ph_iv_LX(1:2,:) = ph_iv_scatrate(GAAS,L,X)
!        ph_pol(L,1:2,:) = ph_pol_scatrate(GAAS,L)
!        ph_in_ac(L,1:2,:) = ph_inelastic_ac(GAAS,L)
!
!        ! X valley
!        ph_ac(X,:) = ph_ac_scatrate(GAAS,X,2)
!        ph_nonp(X,1:2,:) = ph_nonp_scatrate(GAAS,X,2)
!        ph_iv_XG(1:2,:) = ph_iv_scatrate(GAAS,X,G)
!        ph_iv_XL(1:2,:) = ph_iv_scatrate(GAAS,X,L)
!        ph_iv_XX(1:2,:) = ph_iv_scatrate(GAAS,X,X)
!        ph_pol(X,1:2,:) = ph_pol_scatrate(GAAS,X)
!        ph_in_ac(X,1:2,:) = ph_inelastic_ac(GAAS,X)
!
!        ! Sort all the scattering rates in the final array
!        ! scat_rates(MAT,valley,energy,scattering type): see materialdef.f90 for description
!        do i=1,size(m(1,:))
!        !scat_rates(GAAS,G,:,2,i)   = ph_ac(G,:)
!        scat_rates(GAAS,G,:,2,i)   = ph_in_ac(G,1,:)
!        scat_rates(GAAS,G,:,3,i)   = ph_in_ac(G,2,:)
!        scat_rates(GAAS,G,:,4,i)   = ph_nonp(G,1,:)
!        scat_rates(GAAS,G,:,5,i)   = ph_nonp(G,2,:)
!        scat_rates(GAAS,G,:,6,i)   = ph_pol(G,1,:)
!        scat_rates(GAAS,G,:,7,i)   = ph_pol(G,2,:)
!        scat_rates(GAAS,G,:,8,i)   = ph_iv_GL(1,:)
!        scat_rates(GAAS,G,:,9,i)   = ph_iv_GL(2,:)
!        scat_rates(GAAS,G,:,10,i)  = ph_iv_GX(1,:)
!        scat_rates(GAAS,G,:,11,i)  = ph_iv_GX(2,:)
!
!        ! Compute the impurity scattering rate only for the Gamma-valley, neglected at high electric fields
!        if (IMPURITY_SCAT .eq. 1) then
!            if (m_material(i) .eq. GAAS) then
!                scat_rates(GAAS,G,:,15,i)  = impurity_ridley(GAAS,G,donor_density(i))
!            end if
!        end if
!
!        !scat_rates(GAAS,L,:,2,i)   = ph_ac(L,:)
!        scat_rates(GAAS,L,:,2,i)   = ph_in_ac(L,1,:)
!        scat_rates(GAAS,L,:,3,i)   = ph_in_ac(L,2,:)
!        scat_rates(GAAS,L,:,4,i)   = ph_nonp(L,1,:)
!        scat_rates(GAAS,L,:,5,i)   = ph_nonp(L,2,:)
!        scat_rates(GAAS,L,:,6,i)   = ph_pol(L,1,:)
!        scat_rates(GAAS,L,:,7,i)   = ph_pol(L,2,:)
!        scat_rates(GAAS,L,:,8,i)   = ph_iv_LG(1,:)
!        scat_rates(GAAS,L,:,9,i)   = ph_iv_LG(2,:)
!        scat_rates(GAAS,L,:,10,i)  = ph_iv_LL(1,:)
!        scat_rates(GAAS,L,:,11,i)  = ph_iv_LL(2,:)
!        scat_rates(GAAS,L,:,12,i)  = ph_iv_LX(1,:)
!        scat_rates(GAAS,L,:,13,i)  = ph_iv_LX(2,:)
!
!        !scat_rates(GAAS,X,:,2,i)   = ph_ac(X,:)
!        scat_rates(GAAS,X,:,2,i)   = ph_in_ac(X,1,:)
!        scat_rates(GAAS,X,:,3,i)   = ph_in_ac(X,2,:)
!        scat_rates(GAAS,X,:,4,i)   = ph_nonp(X,1,:)
!        scat_rates(GAAS,X,:,5,i)   = ph_nonp(X,2,:)
!        scat_rates(GAAS,X,:,6,i)   = ph_pol(X,1,:)
!        scat_rates(GAAS,X,:,7,i)   = ph_pol(X,2,:)
!        scat_rates(GAAS,X,:,8,i)   = ph_iv_XG(1,:)
!        scat_rates(GAAS,X,:,9,i)   = ph_iv_XG(2,:)
!        scat_rates(GAAS,X,:,10,i)  = ph_iv_XL(1,:)
!        scat_rates(GAAS,X,:,11,i)  = ph_iv_XL(2,:)
!        scat_rates(GAAS,X,:,12,i)  = ph_iv_XX(1,:)
!        scat_rates(GAAS,X,:,13,i)  = ph_iv_XX(2,:)
!        end do
!
!        ! Find boundary for the upper bound calculations
!        call calc_upper_bound(GAAS,G)
!        call calc_upper_bound(GAAS,L)
!        call calc_upper_bound(GAAS,X)
!        ! With this we can proceed to compute the self scattering rate, basically the difference
!        ! between the real scattering rates and the upper bound
!        self_scat(G,:,:) = self_scatrate(GAAS,G)
!        self_scat(L,:,:) = self_scatrate(GAAS,L)
!        self_scat(X,:,:) = self_scatrate(GAAS,X)
!        do i=1,size(m(1,:))
!        scat_rates(GAAS,G,:,1,i) = self_scat(G,i,:)
!        scat_rates(GAAS,L,:,1,i) = self_scat(L,i,:)
!        scat_rates(GAAS,X,:,1,i) = self_scat(X,i,:)
!        end do
!        call normalize_scatrates(GAAS,G)
!        call normalize_scatrates(GAAS,L)
!        call normalize_scatrates(GAAS,X)
!    end subroutine calc_scatrates_GAAS

    subroutine calc_scatrates_INGAAS()
        integer                              :: i, idx
        real(kind=s),dimension(2,size_E)     :: ph_iv_GL, ph_iv_GX, &
                                              & ph_iv_LG, ph_iv_LL, ph_iv_LX, &
                                              & ph_iv_XG, ph_iv_XL, ph_iv_XX
        real(kind=s),dimension(3,2,size_E)   :: ph_nonp, ph_pol, ph_in_ac
        real(kind=s),dimension(:,:,:),allocatable :: self_scat

        !allocate(self_scat(3,size(m(1,:)),size_E))
        allocate(self_scat(3,maxval(scat_m),size_E))

        ph_nonp(G,1:2,:) = ph_nonp_scatrate(INGAAS,G,1)
        ph_nonp(L,1:2,:) = ph_nonp_scatrate(INGAAS,L,2)
        ph_nonp(X,1:2,:) = ph_nonp_scatrate(INGAAS,X,1)

        ph_pol(G,1:2,:) = ph_pol_scatrate(INGAAS,G)
        ph_pol(L,1:2,:) = ph_pol_scatrate(INGAAS,L)
        ph_pol(X,1:2,:) = ph_pol_scatrate(INGAAS,X)

        ph_iv_GL = ph_iv_scatrate(INGAAS,G,L)
        ph_iv_GX = ph_iv_scatrate(INGAAS,G,X)

        ph_iv_LG(1:2,:) = ph_iv_scatrate(INGAAS,L,G)
        ph_iv_LL(1:2,:) = ph_iv_scatrate(INGAAS,L,L)
        ph_iv_LX(1:2,:) = ph_iv_scatrate(INGAAS,L,X)

        ph_iv_XG(1:2,:) = ph_iv_scatrate(INGAAS,X,G)
        ph_iv_XL(1:2,:) = ph_iv_scatrate(INGAAS,X,L)
        ph_iv_XX(1:2,:) = ph_iv_scatrate(INGAAS,X,X)

        ph_in_ac(G,1:2,:) = ph_inelastic_ac(INGAAS,G)
        ph_in_ac(L,1:2,:) = ph_inelastic_ac(INGAAS,L)
        ph_in_ac(X,1:2,:) = ph_inelastic_ac(INGAAS,X)

        ! Sort all the scattering rates in the final array
        ! scat_rates(MAT,valley,energy,scattering type): see materialdef.f90 for description
        !do i=1,size(m(1,:))
        do idx=1,size(scat_m)
            if ((scat_rate_comp(idx) == 1) .and. (m_material(idx) == INGAAS)) then
            !print *, material, i
        i = scat_m(idx)
        !scat_rates(INGAAS,G,:,2,i)   = ph_ac_scatrate(INGAAS,G,1)
        scat_rates(i,G,:,2)   = ph_in_ac(G,1,:)
        scat_rates(i,G,:,3)   = ph_in_ac(G,2,:)
        scat_rates(i,G,:,4)   = ph_nonp(G,1,:)
        scat_rates(i,G,:,5)   = ph_nonp(G,2,:)
        scat_rates(i,G,:,6)   = ph_pol(G,1,:)
        scat_rates(i,G,:,7)   = ph_pol(G,2,:)
        scat_rates(i,G,:,8)   = ph_iv_GL(1,:)
        scat_rates(i,G,:,9)   = ph_iv_GL(2,:)
        scat_rates(i,G,:,10)  = ph_iv_GX(1,:)
        scat_rates(i,G,:,11)  = ph_iv_GX(2,:)
        scat_rates(i,G,:,14)  = alloy_scatrate(INGAAS,G,mole_x(INGAAS))

        if ((IMPURITY_SCAT .eq. 1) .and. (donor_density(i) > 0.0_s)) then
            !if (m_material(i) .eq. INGAAS) then
                scat_rates(i,G,:,15)  = impurity_ridley(INGAAS,G,donor_density(i))
            !end if
        end if
        !scat_rates(INGAAS,L,:,2,i)   = ph_ac_scatrate(INGAAS,L,1)
        scat_rates(i,L,:,2)   = ph_in_ac(L,1,:)
        scat_rates(i,L,:,3)   = ph_in_ac(L,2,:)
        scat_rates(i,L,:,4)   = ph_nonp(L,1,:)
        scat_rates(i,L,:,5)   = ph_nonp(L,2,:)
        scat_rates(i,L,:,6)   = ph_pol(L,1,:)
        scat_rates(i,L,:,7)   = ph_pol(L,2,:)
        scat_rates(i,L,:,8)   = ph_iv_LG(1,:)
        scat_rates(i,L,:,9)   = ph_iv_LG(2,:)
        scat_rates(i,L,:,10)   = ph_iv_LL(1,:)
        scat_rates(i,L,:,11)  = ph_iv_LL(2,:)
        scat_rates(i,L,:,12)  = ph_iv_LX(1,:)
        scat_rates(i,L,:,13)  = ph_iv_LX(2,:)
        scat_rates(i,L,:,14)  = alloy_scatrate(INGAAS,L,mole_x(INGAAS))

        !scat_rates(INGAAS,X,:,2,i)   = ph_ac_scatrate(INGAAS,X,1)
        scat_rates(i,X,:,2)   = ph_in_ac(X,1,:)
        scat_rates(i,X,:,3)   = ph_in_ac(X,2,:)
        scat_rates(i,X,:,4)   = ph_nonp(X,1,:)
        scat_rates(i,X,:,5)   = ph_nonp(X,2,:)
        scat_rates(i,X,:,6)   = ph_pol(X,1,:)
        scat_rates(i,X,:,7)   = ph_pol(X,2,:)
        scat_rates(i,X,:,8)   = ph_iv_XG(1,:)
        scat_rates(i,X,:,9)   = ph_iv_XG(2,:)
        scat_rates(i,X,:,10)  = ph_iv_XL(1,:)
        scat_rates(i,X,:,11)  = ph_iv_XL(2,:)
        scat_rates(i,X,:,12)  = ph_iv_XX(1,:)
        scat_rates(i,X,:,13)  = ph_iv_XX(2,:)
        scat_rates(i,X,:,14)  = alloy_scatrate(INGAAS,X,mole_x(INGAAS))
!        end do

        ! find borders for the upper bound calculations
        call calc_upper_bound(i,G)
        call calc_upper_bound(i,L)
        call calc_upper_bound(i,X)

        ! With this we can proceed to compute the self scattering rate, basically the difference
        ! between the real scattering rates and the upper bound
        self_scat(G,:,:) = self_scatrate(i,G)
        self_scat(L,:,:) = self_scatrate(i,L)
        self_scat(X,:,:) = self_scatrate(i,X)
        !do i=1,size(m(1,:))
        scat_rates(i,G,:,1) = self_scat(G,i,:)
        scat_rates(i,L,:,1) = self_scat(L,i,:)
        scat_rates(i,X,:,1) = self_scat(X,i,:)

        call normalize_scatrates(i,G)
        call normalize_scatrates(i,L)
        call normalize_scatrates(i,X)
        end if
        end do
    end subroutine calc_scatrates_INGAAS

    subroutine calc_scatrates_INALAS()
        integer                              :: i, idx
        real(kind=s),dimension(2,size_E)     :: ph_iv_GL, ph_iv_GX, &
                                              & ph_iv_LG, ph_iv_LL, ph_iv_LX, &
                                              & ph_iv_XG, ph_iv_XL, ph_iv_XX
        real(kind=s),dimension(3,2,size_E)   :: ph_nonp, ph_pol, ph_in_ac
        real(kind=s),dimension(:,:,:),allocatable :: self_scat

        allocate(self_scat(3,maxval(scat_m),size_E))

        ph_nonp(G,1:2,:) = ph_nonp_scatrate(INALAS,G,1)
        ph_nonp(L,1:2,:) = ph_nonp_scatrate(INALAS,L,2)
        ph_nonp(X,1:2,:) = ph_nonp_scatrate(INALAS,X,1)

        ph_pol(G,1:2,:) = ph_pol_scatrate(INALAS,G)
        ph_pol(L,1:2,:) = ph_pol_scatrate(INALAS,L)
        ph_pol(X,1:2,:) = ph_pol_scatrate(INALAS,X)

        ph_iv_GL = ph_iv_scatrate(INALAS,G,L)
        ph_iv_GX = ph_iv_scatrate(INALAS,G,X)

        ph_iv_LG(1:2,:) = ph_iv_scatrate(INALAS,L,G)
        ph_iv_LL(1:2,:) = ph_iv_scatrate(INALAS,L,L)
        ph_iv_LX(1:2,:) = ph_iv_scatrate(INALAS,L,X)

        ph_iv_XG(1:2,:) = ph_iv_scatrate(INALAS,X,G)
        ph_iv_XL(1:2,:) = ph_iv_scatrate(INALAS,X,L)
        ph_iv_XX(1:2,:) = ph_iv_scatrate(INALAS,X,X)

        ph_in_ac(G,1:2,:) = ph_inelastic_ac(INALAS,G)
        ph_in_ac(L,1:2,:) = ph_inelastic_ac(INALAS,L)
        ph_in_ac(X,1:2,:) = ph_inelastic_ac(INALAS,X)

        ! Sort all the scattering rates in the final array
        ! scat_rates(MAT,valley,energy,scattering type): see materialdef.f90 for description
        do idx=1,size(scat_m)
            if ((scat_rate_comp(idx) == 1) .and. (m_material(idx) == INALAS)) then
        i = scat_m(idx)
        !scat_rates(INGAAS,G,:,2,i)   = ph_ac_scatrate(INGAAS,G,1)
        scat_rates(i,G,:,2)   = ph_in_ac(G,1,:)
        scat_rates(i,G,:,3)   = ph_in_ac(G,2,:)
        scat_rates(i,G,:,4)   = ph_nonp(G,1,:)
        scat_rates(i,G,:,5)   = ph_nonp(G,2,:)
        scat_rates(i,G,:,6)   = ph_pol(G,1,:)
        scat_rates(i,G,:,7)   = ph_pol(G,2,:)
        scat_rates(i,G,:,8)   = ph_iv_GL(1,:)
        scat_rates(i,G,:,9)   = ph_iv_GL(2,:)
        scat_rates(i,G,:,10)  = ph_iv_GX(1,:)
        scat_rates(i,G,:,11)  = ph_iv_GX(2,:)
        scat_rates(i,G,:,14)  = alloy_scatrate(INALAS,G,mole_x(INALAS))

        if ((IMPURITY_SCAT .eq. 1) .and. (donor_density(i) > 0.0_s)) then
            !if (m_material(i) .eq. INALAS) then
                scat_rates(i,G,:,15)  = impurity_ridley(INALAS,G,donor_density(i))
            !end if
        end if
        !scat_rates(INGAAS,L,:,2,i)   = ph_ac_scatrate(INGAAS,L,1)
        scat_rates(i,L,:,2)   = ph_in_ac(L,1,:)
        scat_rates(i,L,:,3)   = ph_in_ac(L,2,:)
        scat_rates(i,L,:,4)   = ph_nonp(L,1,:)
        scat_rates(i,L,:,5)   = ph_nonp(L,2,:)
        scat_rates(i,L,:,6)   = ph_pol(L,1,:)
        scat_rates(i,L,:,7)   = ph_pol(L,2,:)
        scat_rates(i,L,:,8)   = ph_iv_LG(1,:)
        scat_rates(i,L,:,9)   = ph_iv_LG(2,:)
        scat_rates(i,L,:,10)  = ph_iv_LL(1,:)
        scat_rates(i,L,:,11)  = ph_iv_LL(2,:)
        scat_rates(i,L,:,12)  = ph_iv_LX(1,:)
        scat_rates(i,L,:,13)  = ph_iv_LX(2,:)
        scat_rates(i,L,:,14)  = alloy_scatrate(INALAS,L,mole_x(INALAS))

        !scat_rates(INGAAS,X,:,2,i)   = ph_ac_scatrate(INGAAS,X,1)
        scat_rates(i,X,:,2)   = ph_in_ac(X,1,:)
        scat_rates(i,X,:,3)   = ph_in_ac(X,2,:)
        scat_rates(i,X,:,4)   = ph_nonp(X,1,:)
        scat_rates(i,X,:,5)   = ph_nonp(X,2,:)
        scat_rates(i,X,:,6)   = ph_pol(X,1,:)
        scat_rates(i,X,:,7)   = ph_pol(X,2,:)
        scat_rates(i,X,:,8)   = ph_iv_XG(1,:)
        scat_rates(i,X,:,9)   = ph_iv_XG(2,:)
        scat_rates(i,X,:,10)  = ph_iv_XL(1,:)
        scat_rates(i,X,:,11)  = ph_iv_XL(2,:)
        scat_rates(i,X,:,12)  = ph_iv_XX(1,:)
        scat_rates(i,X,:,13)  = ph_iv_XX(2,:)
        scat_rates(i,X,:,14)  = alloy_scatrate(INALAS,X,mole_x(INALAS))
        !end do
        ! find borders for the upper bound calculations
        call calc_upper_bound(i,G)
        call calc_upper_bound(i,L)
        call calc_upper_bound(i,X)

        ! With this we can proceed to compute the self scattering rate, basically the difference
        ! between the real scattering rates and the upper bound
        self_scat(G,:,:) = self_scatrate(i,G)
        self_scat(L,:,:) = self_scatrate(i,L)
        self_scat(X,:,:) = self_scatrate(i,X)
        !do i=1,size(m(1,:))
        scat_rates(i,G,:,1) = self_scat(G,i,:)
        scat_rates(i,L,:,1) = self_scat(L,i,:)
        scat_rates(i,X,:,1) = self_scat(X,i,:)
        !end do

        call normalize_scatrates(i,G)
        call normalize_scatrates(i,L)
        call normalize_scatrates(i,X)
        end if
        end do
    end subroutine calc_scatrates_INALAS
!
!    subroutine calc_scatrates_INALAS()
!        integer                              :: i
!        real(kind=s),dimension(2,size_E)     :: ph_iv_GL, ph_iv_GX, &
!                                              & ph_iv_LG, ph_iv_LL, ph_iv_LX, &
!                                              & ph_iv_XG, ph_iv_XL, ph_iv_XX
!        real(kind=s),dimension(3,2,size_E)   :: ph_nonp, ph_pol
!        real(kind=s),dimension(:,:,:),allocatable :: self_scat
!
!        ph_nonp(G,1:2,:) = ph_nonp_scatrate(INALAS,G,1)
!        ph_nonp(L,1:2,:) = ph_nonp_scatrate(INALAS,L,2)
!        ph_nonp(X,1:2,:) = ph_nonp_scatrate(INALAS,X,1)
!
!        ph_pol(G,1:2,:) = ph_pol_scatrate(INALAS,G)
!        ph_pol(L,1:2,:) = ph_pol_scatrate(INALAS,L)
!        ph_pol(X,1:2,:) = ph_pol_scatrate(INALAS,X)
!
!        ph_iv_GL = ph_iv_scatrate(INALAS,G,L)
!        ph_iv_GX = ph_iv_scatrate(INALAS,G,X)
!
!        ph_iv_LG(1:2,:) = ph_iv_scatrate(INALAS,L,G)
!        ph_iv_LL(1:2,:) = ph_iv_scatrate(INALAS,L,L)
!        ph_iv_LX(1:2,:) = ph_iv_scatrate(INALAS,L,X)
!
!        ph_iv_XG(1:2,:) = ph_iv_scatrate(INALAS,X,G)
!        ph_iv_XL(1:2,:) = ph_iv_scatrate(INALAS,X,L)
!        ph_iv_XX(1:2,:) = ph_iv_scatrate(INALAS,X,X)
!
!        ! Sort all the scattering rates in the final array
!        ! scat_rates(MAT,valley,energy,scattering type): see materialdef.f90 for description
!        do i=1,size(m(:,1))
!        scat_rates(INALAS,G,:,2,i)   = ph_ac_scatrate(INALAS,G,1)
!        scat_rates(INALAS,G,:,4,i)   = ph_nonp(G,1,:)
!        scat_rates(INALAS,G,:,5,i)   = ph_nonp(G,2,:)
!        scat_rates(INALAS,G,:,6,i)   = ph_pol(G,1,:)
!        scat_rates(INALAS,G,:,7,i)   = ph_pol(G,2,:)
!        scat_rates(INALAS,G,:,8,i)   = ph_iv_GL(1,:)
!        scat_rates(INALAS,G,:,9,i)   = ph_iv_GL(2,:)
!        scat_rates(INALAS,G,:,10,i)  = ph_iv_GX(1,:)
!        scat_rates(INALAS,G,:,11,i)  = ph_iv_GX(2,:)
!        scat_rates(INALAS,G,:,14,i)  = alloy_scatrate(INALAS,G,0.53_s)
!
!        if (m_material(i) .eq. INALAS) then
!            scat_rates(INALAS,G,:,15,i)  = impurity_ridley(INALAS,G,donor_density(i))
!        end if
!
!        scat_rates(INALAS,L,:,2,i)   = ph_ac_scatrate(INALAS,L,1)
!        scat_rates(INALAS,L,:,4,i)   = ph_nonp(L,1,:)
!        scat_rates(INALAS,L,:,5,i)   = ph_nonp(L,2,:)
!        scat_rates(INALAS,L,:,6,i)   = ph_pol(L,1,:)
!        scat_rates(INALAS,L,:,7,i)   = ph_pol(L,2,:)
!        scat_rates(INALAS,L,:,8,i)   = ph_iv_LG(1,:)
!        scat_rates(INALAS,L,:,9,i)   = ph_iv_LG(2,:)
!        scat_rates(INALAS,L,:,10,i)  = ph_iv_LL(1,:)
!        scat_rates(INALAS,L,:,11,i)  = ph_iv_LL(2,:)
!        scat_rates(INALAS,L,:,12,i)  = ph_iv_LX(1,:)
!        scat_rates(INALAS,L,:,13,i)  = ph_iv_LX(2,:)
!        scat_rates(INALAS,L,:,14,i)  = alloy_scatrate(INALAS,L,0.53_s)
!
!        scat_rates(INALAS,X,:,2,i)   = ph_ac_scatrate(INALAS,X,1)
!        scat_rates(INALAS,X,:,4,i)   = ph_nonp(X,1,:)
!        scat_rates(INALAS,X,:,5,i)   = ph_nonp(X,2,:)
!        scat_rates(INALAS,X,:,6,i)   = ph_pol(X,1,:)
!        scat_rates(INALAS,X,:,7,i)   = ph_pol(X,2,:)
!        scat_rates(INALAS,X,:,8,i)   = ph_iv_XG(1,:)
!        scat_rates(INALAS,X,:,9,i)   = ph_iv_XG(2,:)
!        scat_rates(INALAS,X,:,10,i)   = ph_iv_XL(1,:)
!        scat_rates(INALAS,X,:,11,i)  = ph_iv_XL(2,:)
!        scat_rates(INALAS,X,:,12,i)  = ph_iv_XX(1,:)
!        scat_rates(INALAS,X,:,13,i)  = ph_iv_XX(2,:)
!        scat_rates(INALAS,X,:,14,i)  = alloy_scatrate(INALAS,X,0.53_s)
!        end do
!
!        ! Find borders for the upper bound calculations
!        call calc_upper_bound(INALAS,G)
!        call calc_upper_bound(INALAS,L)
!        call calc_upper_bound(INALAS,X)
!
!        ! With this we can proceed to compute the self scattering rate, basically the difference
!        ! between the real scattering rates and the upper bound
!        self_scat(G,:,:) = self_scatrate(INALAS,G)
!        print *, 'test3'
!        self_scat(L,:,:) = self_scatrate(INALAS,L)
!        self_scat(X,:,:) = self_scatrate(INALAS,X)
!        print *, 'test5'
!        do i=1,size(m(:,1))
!        scat_rates(INALAS,G,:,1,i) = self_scat(G,i,:)
!        scat_rates(INALAS,L,:,1,i) = self_scat(L,i,:)
!        scat_rates(INALAS,X,:,1,i) = self_scat(X,i,:)
!        end do
!
!        call normalize_scatrates(INALAS,G)
!        call normalize_scatrates(INALAS,L)
!        call normalize_scatrates(INALAS,X)
!    end subroutine calc_scatrates_INALAS

    ! Normalize tables for the scattering rates
    subroutine normalize_scatrates(scat_i,valley)
        integer                           :: valley, i, j, k, scat_i !MAT,

        !do k=1,size(m(1,:))
        k = scat_i
            do i=1,size_E
                do j=1,15
                    n_scat_rates(k,valley,i,j) = sum(scat_rates(k,valley,i,1:j))/upper_bound(k,valley,i)
                end do
            end do
        !end do
    end subroutine normalize_scatrates

    ! Compute the constant (for dE) upper bound needed for the free flight determination
    subroutine calc_upper_bound(scat_i,valley) ! MAT,valley
        real(kind=s)                :: max_scatrate
        integer                     :: dE, valley, i, j, k, start_idx, end_idx, last_idx, scat_i
        dE = E2idx(bin_dE)
        last_idx = int(E2idx(3.0_s)/dE)
        end_idx = 0

        !do k=1,size(m(1,:))
        k = scat_i
            do i=1,last_idx
                start_idx = dE*(i-1)+1
                end_idx = dE*i
                max_scatrate = maxval(sum(scat_rates(k,valley,start_idx:end_idx,:),2))
                do j=start_idx,end_idx
                    upper_bound(k,valley,j) = max_scatrate
                end do
            end do
            if (end_idx < size_E) then
                start_idx = (dE*last_idx)+1
                end_idx = size_E
                max_scatrate = maxval(sum(scat_rates(k,valley,start_idx:end_idx,:),2))
                do j=start_idx,end_idx
                    upper_bound(k,valley,j) = max_scatrate
                end do
            end if
        !end do
    end subroutine calc_upper_bound

    ! Calculate the self scattering rate, needed for free flight calculation
    function self_scatrate(scat_i,valley) result(self_scat_new)
        integer                            :: valley, i, k, scat_i ! MAT,
        real(kind=s),dimension(:,:),allocatable     :: self_scat_new
        real(kind=s)                       :: tot_scatrate

        !allocate(self_scat_new(size(m(1,:)),size_E))
        !do k=1,size(m(1,:))
        allocate(self_scat_new(maxval(scat_m),size_E))
        k = scat_i
            do i=1,size_E
                !tot_scatrate = sum(scat_rates(MAT,valley,i,:,k))
                !self_scat_new(k,i) = upper_bound(MAT,valley,i,k)-tot_scatrate
                tot_scatrate = sum(scat_rates(k,valley,i,:))
                self_scat_new(k,i) = upper_bound(k,valley,i)-tot_scatrate
            end do
        !end do
    end function self_scatrate

    ! This function calculates the elastic acoustic phonon
    ! scattering rate. Arguments are the material index, valley and the index for chosing the deformation potential
    function ph_ac_scatrate(MAT,valley,def_idx) result(ph_ac)
        integer                            :: MAT, valley, def_idx
        real(kind=s),dimension(size_E)     :: ph_ac

        ph_ac = pi*kB*T_lattice*def_pot_ac(MAT,def_idx)**2.0_s/(hb*sound_vel(MAT)**2.0_s*rho(MAT))*dos(MAT,valley,E)
    end function ph_ac_scatrate

    ! Function for non-polar optical phonon scattering, arguments
    ! see acoustic phonon scattering rate
    function ph_nonp_scatrate(MAT,valley,def_idx) result(ph_nonp)
        integer                            :: MAT, valley, def_idx
        real(kind=s),dimension(2,size_E)   :: ph_nonp
        real(kind=s),dimension(size_E)     :: E_em, E_abs

        ! Calculate energy array for emission/absorption process
        E_em = build_energy_arr(-phonon_e(MAT,1))
        E_abs = build_energy_arr(phonon_e(MAT,1))

        ! Emission
        ph_nonp(1,:) = def_pot_nonp(MAT,def_idx)**2.0_s*pi/(2.0_s*rho(MAT)*phonon_e(MAT,1)/hb)*(n_op(MAT,1)+1.0_s)* &
                     & dos(MAT,valley,E_em)
        ! Absorption
        ph_nonp(2,:) = def_pot_nonp(MAT,def_idx)**2.0_s*pi/(2.0_s*rho(MAT)*phonon_e(MAT,1)/hb)*(n_op(MAT,1))* &
                     & dos(MAT,valley,E_abs)
    end function ph_nonp_scatrate

    ! Calculation of the polar optical phonon scattering rate. Parameters are Material and valley
    function ph_pol_scatrate(MAT,valley) result(ph_pol)
        integer                            :: MAT, valley
        real(kind=s),dimension(2,size_E)   :: ph_pol
        real(kind=s),dimension(size_E)     :: E_em, E_abs, D_q2_em, D_q2_abs

        ! Calculate energy array for emission/absorption process
        E_em = build_energy_arr(-phonon_e(MAT,1))
        E_abs = build_energy_arr(phonon_e(MAT,1))

        ! Get direction-weighted DOS
        D_q2_em = dw_dos(MAT,valley,E_em)
        D_q2_abs = dw_dos(MAT,valley,E_abs)

        ! Calculate now the final scattering rate for polar optical phonon scattering
        ! Emission
        ph_pol(1,:) = 1.0_s/(2.0_s*hb)*F2(MAT)*(n_op(MAT,1)+1.0_s)*D_q2_em
        ! Absorption
        ph_pol(2,:) = 1.0_s/(2.0_s*hb)*F2(MAT)*n_op(MAT,1)*D_q2_abs

!       !Conventional calculation, but this works follow direction-weighted DOS
!        real(kind=s),dimension(size_E)     :: E_em, E_abs, Gamma_i, Gamma_em, Gamma_abs, &
!                                            & A_em, A_abs, B_em, B_abs, C_em, C_abs, F0_em, F0_abs
!        real(kind=s)                       :: alpha
!
!        ! Do precalculations
!        alpha = nonp(MAT,valley)
!        E_em = build_energy_arr(-phonon_e(MAT,1))
!        E_abs = build_energy_arr(phonon_e(MAT,1))
!        Gamma_i = E*(1+alpha*E)
!        Gamma_em = E_em*(1+alpha*E_em)
!        Gamma_abs = E_abs*(1+alpha*E_abs)
!        A_em = (2.0_s*(1.0_s+alpha*E)*(1.0_s+alpha*E_em)+alpha*(Gamma_i+Gamma_em))**2.0_s
!        A_abs = (2.0_s*(1.0_s+alpha*E)*(1.0_s+alpha*E_abs)+alpha*(Gamma_i+Gamma_abs))**2.0_s
!        B_em = -2.0_s*alpha*sqrt(Gamma_i*Gamma_em)*(4*(1+alpha*E)*(1+alpha*E_em)+ &
!             & alpha*(Gamma_i+Gamma_em))
!        B_abs = -2.0_s*alpha*sqrt(Gamma_i*Gamma_abs)*(4*(1+alpha*E)*(1+alpha*E_abs)+ &
!              & alpha*(Gamma_i+Gamma_abs))
!        C_em = 4.0_s*(1.0_s+alpha*E)*(1.0_s+alpha*E_em)*(1.0_s+2.0_s*alpha*E)*(1.0_s+2.0_s*alpha*E_em)
!        C_abs = 4.0_s*(1.0_s+alpha*E)*(1.0_s+alpha*E_abs)*(1.0_s+2.0_s*alpha*E)*(1.0_s+2.0_s*alpha*E_abs)
!        F0_em = 1.0_s/C_em*(A_em*log(abs((sqrt(Gamma_i)+sqrt(Gamma_em))/(sqrt(Gamma_i)-sqrt(Gamma_em))))+B_em)
!        F0_abs = 1.0_s/C_abs*(A_abs*log(abs((sqrt(Gamma_i)+sqrt(Gamma_abs))/(sqrt(Gamma_i)-sqrt(Gamma_abs))))+B_abs)
!
!        ! Emission
!        ph_pol(1,:) = sqrt(me(MAT,valley))*phonon_e(MAT,1)/hb/(4.0_s*pi*sqrt(2.0_s)*hb)* &
!                    & (1.0_s/eps_hf(MAT)-1.0_s/eps_static(MAT))*(n_op(MAT,1)+1.0_s)*F0_em*(1.0_s+2.0_s*alpha*E_em)/sqrt(Gamma_i)
!        ! Absorption
!        ph_pol(2,:) = sqrt(me(MAT,valley))*phonon_e(MAT,1)/hb/(4.0_s*pi*sqrt(2.0_s)*hb)* &
!                    & (1.0_s/eps_hf(MAT)-1.0_s/eps_static(MAT))*n_op(MAT,1)*F0_abs*(1.0_s+2.0_s*alpha*E_abs)/sqrt(Gamma_i)
    end function ph_pol_scatrate

    ! This function computes the invervalley transition rates, input parameters
    ! are only the material, initial and final valley
    function ph_iv_scatrate(MAT,valley,final_valley) result(ph_iv)
        integer                            :: MAT, valley, final_valley, Zf
        real(kind=s),dimension(2,size_E)   :: ph_iv
        real(kind=s),dimension(size_E)     :: E_em, E_abs
        real(kind=s)                       :: def_pot_iv, ph_e, nop_tmp

        ph_e = 0
        Zf = 0
        nop_tmp = 0
        def_pot_iv = 0
        if ((valley == G) .and. (final_valley == L)) then
            def_pot_iv = def_pot_iv_GL(MAT)
            Zf = 4
            ph_e = phonon_e(MAT,2)
            nop_tmp = n_op(MAT,2)
            E_em = build_energy_arr(-ph_e-valley_offs(MAT,L))
            E_abs = build_energy_arr(ph_e-valley_offs(MAT,L))
        else if ((valley == G) .and. (final_valley == X)) then
            def_pot_iv = def_pot_iv_GX(MAT)
            Zf = 3
            ph_e = phonon_e(MAT,3)
            nop_tmp = n_op(MAT,3)
            E_em = build_energy_arr(-ph_e-valley_offs(MAT,X))
            E_abs = build_energy_arr(ph_e-valley_offs(MAT,X))
        else if ((valley == L) .and. (final_valley == G)) then
            def_pot_iv = def_pot_iv_GL(MAT)
            Zf = 1
            ph_e = phonon_e(MAT,2)
            nop_tmp = n_op(MAT,2)
            E_em = build_energy_arr(-ph_e+valley_offs(MAT,L))
            E_abs = build_energy_arr(ph_e+valley_offs(MAT,L))
        else if ((valley == L) .and. (final_valley == L)) then
            def_pot_iv = def_pot_iv_LL(MAT)
            Zf = 3
            ph_e = phonon_e(MAT,5)
            nop_tmp = n_op(MAT,5)
            E_em = build_energy_arr(-ph_e)
            E_abs = build_energy_arr(ph_e)
        else if ((valley == L) .and. (final_valley == X)) then
            def_pot_iv = def_pot_iv_LX(MAT)
            Zf = 3
            ph_e = phonon_e(MAT,4)
            nop_tmp = n_op(MAT,4)
            E_em = build_energy_arr(-ph_e-(valley_offs(MAT,X)-valley_offs(MAT,L)))
            E_abs = build_energy_arr(ph_e-(valley_offs(MAT,X)-valley_offs(MAT,L)))
        else if ((valley == X) .and. (final_valley == G)) then
            def_pot_iv = def_pot_iv_GX(MAT)
            Zf = 1
            ph_e = phonon_e(MAT,3)
            nop_tmp = n_op(MAT,3)
            E_em = build_energy_arr(-ph_e+valley_offs(MAT,X))
            E_abs = build_energy_arr(ph_e+valley_offs(MAT,X))
        else if ((valley == X) .and. (final_valley == L)) then
            def_pot_iv = def_pot_iv_LX(MAT)
            Zf = 4
            ph_e = phonon_e(MAT,4)
            nop_tmp = n_op(MAT,4)
            E_em = build_energy_arr(-ph_e+(valley_offs(MAT,X)-valley_offs(MAT,L)))
            E_abs = build_energy_arr(ph_e+(valley_offs(MAT,X)-valley_offs(MAT,L)))
        else if ((valley == X) .and. (final_valley == X)) then
            def_pot_iv = def_pot_iv_XX(MAT)
            Zf = 2
            ph_e = phonon_e(MAT,6)
            nop_tmp = n_op(MAT,6)
            E_em = build_energy_arr(-ph_e)
            E_abs = build_energy_arr(ph_e)
        end if

        ! Emission
        ph_iv(1,:) = Zf*def_pot_iv**2.0_s*pi/(2.0_s*rho(MAT)*ph_e/hb)*(nop_tmp+1.0_s)* &
                   & dos(MAT,final_valley,E_em)
        ! Absorption
        ph_iv(2,:) = Zf*def_pot_iv**2.0_s*pi/(2.0_s*rho(MAT)*ph_e/hb)*(nop_tmp+0.0_s)* &
                   & dos(MAT,final_valley,E_abs)
    end function ph_iv_scatrate

    function alloy_scatrate(MAT,valley,x_mole) result(alloy)
        integer                              :: MAT, valley
        real(kind=s),dimension(size_E)       :: alloy
        real(kind=s)                         :: x_mole

        alloy = 3*pi**3.0_s/(16.0_s*hb)*x_mole*(1.0_s-x_mole)*lattice_const(MAT)**3.0_s/4.0_s* &
              & alloy_pot(MAT)**2.0_s*dos(MAT,valley,E)
    end function alloy_scatrate

    function impurity_ridley(MAT,valley,Nd) result(imp_r)
        integer                              :: MAT, valley
        real(kind=s)                         :: Nd, d
        real(kind=s),dimension(size_E)       :: imp_r, imp_bh, k, v

        ! Average distance between impurities
        d = (2.0_s*pi*Nd)**(-1.0_s/3.0_s)
        ! Be careful, we take the abs of k- and v-vector as the value in x-direction (other components = 0)
        k = sqrt(2.0_s*me(MAT,valley)*E*(1.0_s+nonp(MAT,valley)*E))/hb
        imp_bh = impurity_bh(MAT,valley,Nd)
        ! Calculate velocity to corresponding k-vector
        v = hb*k/(me(MAT,valley)*(1.0_s+2.0_s*nonp(MAT,valley)*E))
        imp_r = v/d*(1.0_s-exp(-d*imp_bh/v))
        !call plot(E,imp_bh,E,imp_r)
    end function impurity_ridley

    function impurity_bh(MAT,valley,Nd) result(imp)
        integer                              :: MAT, valley
        real(kind=s),dimension(size_E)       :: imp
        real(kind=s)                         :: Nd, beta2, ebeta

        beta2 = Nd*q**2.0_s/(eps_static(MAT)*kB*T_lattice)
        !print *, beta2
        !print *, eps_static(MAT)
        ebeta = hb**2.0_s*beta2/(2.0_s*me(MAT,valley))
        !imp = (2.0_s*sqrt(2.0_s)*pi*Nd*sqrt(me(MAT,valley))/((4.0_s*pi*eps_static(MAT))**2.0_s*hb**2.0_s*beta2)) &
        !    & *(E*(1.0_s+nonp(MAT,valley)*E))**(-1.0_s/2.0_s)*(1.0_s+2.0_s*nonp(MAT,valley)*E)
        !print *, 2.0_s*sqrt(2.0_s)*pi*Nd*sqrt(me(MAT,valley))/((4.0_s*pi*eps_static(MAT))**2.0_s*hb**2.0_s*beta2)
        imp = (2.0_s**(5.0_s/2.0_s)*pi*Nd)/((4.0_s*pi*eps_static(MAT))**2.0_s*ebeta**(4.0_s/2.0_s)*sqrt(me(MAT,valley))) &
            & *sqrt(E*(1.0_s+nonp(MAT,valley)*E))*(1.0_s+2.0_s*nonp(MAT,valley)*E)/(1.0_s+4.0_s*(E*(1.0_s+nonp(MAT,valley)*E) &
            & /ebeta))
        !print *, E(10), imp(10)
    end function impurity_bh

    function ph_inelastic_ac(MAT, valley) result(ph_ine_ac)
        integer                              :: MAT, valley, i
        real(kind=s),dimension(2,size_E)     :: ph_ine_ac
        real(kind=s)                         :: C, gamma_i, ge_mu
        real(kind=s)                         :: e_mu, x1a, x2a, x1e, x2e, pre

        ph_ine_ac = 0.0_s
        e_mu = (me(MAT,valley)*sound_vel(MAT)**2.0_s)/2.0_s
        C = 4.0_s*sqrt(e_mu)/(kB*T_lattice*(1.0_s-4.0_s*nonp(MAT,valley)*e_mu))
        ge_mu = e_mu/(1.0_s-4.0_s*nonp(MAT,valley)*e_mu)
        do i=1,size_E
            pre = sqrt(me(MAT,valley))*(kB*T_lattice)**3.0_s*def_pot_ac(MAT,1)**2.0_s/(2.0_s**(5.0_s/2.0_s)*pi*hb**4.0_s* &
                & sound_vel(MAT)**4.0_s*rho(MAT))*(E(i)*(1.0_s+nonp(MAT,valley)*E(i)))**(-1.0_s/2.0_s)
            gamma_i = E(i)*(1.0_s+nonp(MAT,valley)*E(i))
            if (gamma_i .lt. ge_mu) then
                x1a = C*(sqrt(e_mu)*(1.0_s+2.0_s*nonp(MAT,valley)*E(i))-sqrt(gamma_i))
                x2a = C*(sqrt(e_mu)*(1.0_s+2.0_s*nonp(MAT,valley)*E(i))+sqrt(gamma_i))
            else
                x1a = 0.0_s
                x2a = C*(sqrt(gamma_i)+sqrt(e_mu)*(1.0_s+2.0_s*nonp(MAT,valley)*E(i)))
                x1e = 0.0_s
                x2e = C*(sqrt(gamma_i)-sqrt(e_mu)*(1.0_s+2.0_s*nonp(MAT,valley)*E(i)))
            end if

            ! Emission
            if (gamma_i .ge. ge_mu) then
            ph_ine_ac(1,i) = pre*((1.0_s+2.0_s*nonp(MAT,valley)*E(i))*(int_G1(x2e)-int_G1(x1e)) &
                           & -2.0_s*nonp(MAT,valley)*kB*T_lattice*(int_G2(x2e)-int_G2(x1e)))
            !ph_ine_ac(1,i) = pre*(int_G1(x2e)-int_G1(x1e))
            end if
            ! Absorption
            ph_ine_ac(2,i) = pre*((1.0_s+2.0_s*nonp(MAT,valley)*E(i))*(int_F1(x2a)-int_F1(x1a)) &
                           & +2.0_s*nonp(MAT,valley)*kB*T_lattice*(int_F2(x2a)-int_F2(x1a)))
            !ph_ine_ac(2,i) = pre*(int_F1(x2a)-int_F1(x1a))

        end do
    end function ph_inelastic_ac

    function i_ac_F1(x) result(res)
        real(kind=s)                        :: x, res, x_a

        x_a = 2.6_s
        if (x .le. x_a) then
            res = x**2.0_s/2.0_s-x**3.0_s/6.0_s+x**4.0_s/48.0_s-x**6.0_s/4320.0_s+ &
                & x**8.0_s/241920.0_s-x**10.0_s/1209600.0_s+x**12.0_s/622702080.0_s
        else
            !res = 0.0_s
            res = x_a**2.0_s/2.0_s-x_a**3.0_s/6.0_s+x_a**4.0_s/48.0_s-x_a**6.0_s/4320.0_s+ &
                & x_a**8.0_s/241920.0_s-x_a**10.0_s/1209600.0_s+x_a**12.0_s/622702080.0_s+exp(-x_a)*(x_a**2.0_s+2.0_s*x_a+2.0_s)- &
                & exp(-x)*(x**2.0_s+2.0_s*x+2.0_s)
        end if
        !tmp = int_F1(x)
        !print *, x, res, tmp
    end function i_ac_F1

    function i_ac_F2(x) result(res)
        real(kind=s)                        :: x, res, x_a

        x_a = 5.5_s
        if (x .le. x_a) then
            res = x**3.0_s/3.0_s-x**4.0_s/8.0_s+x**5.0_s/60.0_s-x**7.0_s/5040.0_s+ &
                & x**9.0_s/272160.0_s-x**11.0_s/143305600.0_s+x**12.0_s/622702080.0_s
        else
            !res = 0.0_s
            res = x**3.0_s/3.0_s-x**4.0_s/8.0_s+x**5.0_s/60.0_s-x**7.0_s/5040.0_s+ &
        & x**9.0_s/272160.0_s-x**11.0_s/143305600.0_s+x**12.0_s/622702080.0_s+exp(-x_a)*(x_a**3.0_s+3*x_a**2.0_s+6.0_s*x_a+6.0_s) &
        & -exp(-x)*(x**3.0_s+3.0_s*x**2.0_s+6.0_s*x+6.0_s)
            !    & x_a**8.0_s/241920_s-x_a**10.0_s/1209600_s+x_a**12.0_s/622702080_s+exp(-x_a)*(x_a**2.0_s+2.0_s*x_a+2.0_s)- &
            !    & exp(-x)*(x**2.0_s+2.0_s*x+2.0_s)
        end if
    end function i_ac_F2

    function i_ac_G1(x) result(res)
        real(kind=s)                        :: x, res, x_e

        x_e = 4.0_s
        if (x .le. x_e) then
            res = x**2.0_s/2.0_s+x**3.0_s/6.0_s+x**6.0_s/4320.0_s+x**8.0_s/241920.0_s- &
                & x**10.0_s/12096000.0_s+x**12.0_s/622702080.0_s
        else
            !res = 0.0_s
            res = x_e**2.0_s/2.0_s+x_e**3.0_s/6.0_s+x_e**6.0_s/4320.0_s+x_e**8.0_s/241920.0_s- &
                & x_e**10.0_s/12096000.0_s+x_e**12.0_s/622702080.0_s+exp(-x_e)*(x_e**2.0_s+2*x_e+2.0_s) &
              & -x_e**3.0_s/3.0_s-exp(-x)*(x**2.0_s+2.0_s*x+2.0_s)+x**3.0_s/3.0_s
        end if

        !print *, x
        !tmp = int_G1(x)
        !print *, x, res, tmp
    end function i_ac_G1

    function i_ac_G2(x) result(res)
        real(kind=s)                        :: x, res, x_e

        x_e = 4.0_s
        if (x .le. x_e) then
            res = x**3.0_s/3.0_s+x**4.0_s/8.0_s+x**5.0_s/60.0_s-x**7.0_s/5040.0_s+ &
                & x**9.0_s/272160.0_s-x**11.0_s/143305600.0_s+x**13.0_s/622702080.0_s
        else
            !res = 0.0_s
            res = x_e**3.0_s/3.0_s+x_e**4.0_s/8.0_s+x_e**5.0_s/60.0_s-x_e**7.0_s/5040.0_s+ &
   & x_e**9.0_s/272160.0_s-x_e**11.0_s/143305600.0_s+x_e**13.0_s/622702080.0_s+exp(-x_e)*(x_e**2.0_s+3*x_e**2.0_s+6.0_s*x_e+6.0_s) &
   & -x_e**4.0_s/4.0_s-exp(-x)*(x**3.0_s+3.0_s*x**2.0_s+6.0_s*x+6.0_s)+x**4.0_s/4.0_s
        end if
    end function i_ac_G2

    function dos(MAT,valley,E_tmp) result(D)
        !implicit none
        integer                              :: MAT, valley
        real(kind=s),dimension(size_E)       :: E_tmp, D

        ! calculate DOS for given MATERIAL and valley
        D = (2.0_s*me(MAT,valley))**(3.0/2.0)/(2.0_s*pi**2.0_s*hb**3.0_s)* &
          & sqrt(E_tmp*(1.0_s+nonp(MAT,valley)*E_tmp))*(1.0_s+2.0_s*nonp(MAT,valley)*E_tmp)
    end function dos

    function dw_dos(MAT,valley,final_energy) result(D_q2)
        integer                              :: MAT, valley
        real(kind=s)                         :: alpha
        real(kind=s),dimension(size_E)       :: final_energy, Gamma_i, Gamma_n, A, B, C, F0, D_q2

        alpha = nonp(MAT,valley)
        Gamma_i = E*(1+alpha*E)
        Gamma_n = final_energy*(1+alpha*final_energy)
        A = (2.0_s*(1.0_s+alpha*E)*(1.0_s+alpha*final_energy)+alpha*(Gamma_i+Gamma_n))**2.0_s
        B = -2.0_s*alpha*sqrt(Gamma_i*Gamma_n)*(4*(1+alpha*E)*(1+alpha*final_energy)+ &
             & alpha*(Gamma_i+Gamma_n))
        C = 4.0_s*(1.0_s+alpha*E)*(1.0_s+alpha*final_energy)*(1.0_s+2.0_s*alpha*E)*(1.0_s+2.0_s*alpha*final_energy)
        F0 = 1.0_s/C*(A*log(abs((sqrt(Gamma_i)+sqrt(Gamma_n))/(sqrt(Gamma_i)-sqrt(Gamma_n))))+B)
        D_q2 = sqrt(me(MAT,valley))*phonon_e(MAT,1)/hb/(sqrt(2.0_s)*hb)*(1.0_s/eps_hf(MAT)-1.0_s/eps_static(MAT))* &
             & (1+2*alpha*final_energy)/sqrt(Gamma_i)*F0/(2.0_s*pi/hb*F2(MAT))
    end function dw_dos

    function build_energy_arr(energy_diff) result(E_result)
        integer                              :: i
        real(kind=s),dimension(size_E)       :: E_result
        real(kind=s)                         :: energy_diff

        do i=1,size_E
            E_result(i) = E(i)+energy_diff
            if (E_result(i) < 0) then
                E_result(i) = 0
            end if
        end do
    end function build_energy_arr

    function E2idx(E) result(idx)
        integer         :: idx
        real(kind=s)    :: E
        idx = int(E*100000.0_s) ! 0
        if (idx == 0) then
            idx = 1
        end if
    end function E2idx
end module scattering_rates
