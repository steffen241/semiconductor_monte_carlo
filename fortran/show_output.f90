module show_output
    use gnufor2
    use physconst
    use materialdef
    use scattering_rates
    use config
    use mc_core
    use poisson
    use hdf5
    use device
    implicit none
    real(kind=s),dimension(200,200)           :: mean_e = 0, mean_vel = 0, bulk_time=0, mean_vel_G=0, mean_vel_L=0,mean_vel_X=0
    real(kind=s),dimension(200,200)           :: mean_el_temp_G=0, mean_el_temp_L=0, mean_el_temp_X=0
    real(kind=s),dimension(200,200)           :: mean_el_fermi_G=0, mean_el_fermi_L=0, mean_el_fermi_X=0
    real(kind=s),dimension(200)               :: mean_vel_steps, mean_e_steps, mean_vel_tot, mean_e_tot
    integer,dimension(200,200,3)              :: valley_occ
    real(kind=s),dimension(3,200,200)         :: mean_energy = 0

    contains

#if DIM == 2
    subroutine count_carriers()
        integer                :: count_c, i
        integer,dimension(2)   :: count_r

        count_c = 0
        count_r = 0

        !$omp parallel do reduction(+:count_c,count_r)
        do i=1,max_carriers
            if ((c(i)%valley > 0) .and. (c(i)%op == 0)) then
                count_c = count_c+1
            elseif ((c(i)%valley > 0) .and. (c(i)%op == 1)) then
                if (c(i)%res_id == 1) then
                    count_r(1) = count_r(1)+1
                end if
                if (c(i)%res_id == 2) then
                    count_r(2) = count_r(2)+1
                end if
            end if
        end do
        !$omp end parallel do
        print *, char(9), char(9), 'Number of carriers in device/reservoir:', count_c, count_r
        res_2(1)%max_c = count_r(1)
        !res_2(2)%max_c = count_r(2)
        if ((n_res >= 1)) then
        if ((res_2(1)%max_c > 0) .and. (count_r(1) < res_2(1)%max_c))then
!            print *, res_2(1)%max_c-count_r(1)
            res_2(1)%add_carriers = res_2(1)%max_c-count_r(1)
        end if
        end if
!
        if ((n_res >= 2)) then
        if ((res_2(2)%max_c > 0) .and. (count_r(2) < res_2(2)%max_c))then
            res_2(2)%add_carriers = res_2(2)%max_c-count_r(2)
        end if
        end if
        !print *, res_2(1)%add_carriers,res_2(2)%add_carriers
!        print *, char(9), char(9), el_pot2d(50,86)
    end subroutine count_carriers
#endif

    subroutine save_materialdef()
        integer(HID_T)          :: file_id       ! File identifier
        integer(HID_T)          :: dset_id       ! Dataset identifier
        integer(HID_T)          :: dspace_id     ! Dataspace identifier
        integer(HID_T)          :: group_id      ! Group identifier
        CHARACTER(LEN=40)       :: filename = "/tmp/matdef.h5" ! File name
        integer                 :: error

        integer                               :: i,j
        real(kind=s),dimension(301)           :: temp, em_G
        real(kind=s),dimension(101)           :: conc
        real(kind=s),dimension(301,3)         :: v_offs
        real(kind=s),dimension(301,101,3)     :: v_offs_a
        real(kind=s),dimension(301,101)       :: me_G_a
        real(kind=s),dimension(3,2,101)       :: me_LX
        real(kind=s),dimension(29)            :: x_bow
        !real(kind=s),dimension(29,3)          :: X_offs
        ! GaAs valley offsets
        do i=1,301
            temp(i) = i-1
            v_offs(i,:) = varshni_model(GAAS, dble(i))
            em_G(i) = effective_mass_model(GAAS, dble(i))
        end do
        !call plot(temp,v_offs(:,1),temp,v_offs(:,2),temp,v_offs(:,3))


        call h5open_f(error)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)


        call h5screate_simple_f(1, (/int(301,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "temperature", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, temp, (/int(301,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5gcreate_f(file_id, "gaas", group_id, error)
        call h5screate_simple_f(2, (/int(301,8),int(3,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "gaas/valley_offs", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, v_offs, (/int(301,8),int(3,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(301,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "gaas/em_G", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, em_G, (/int(301,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(2,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "gaas/em_L", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, (/me_l(GAAS,L),me_t(GAAS,L)/), (/int(2,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(2,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "gaas/em_X", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, (/me_l(GAAS,X),me_t(GAAS,X)/), (/int(2,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5gclose_f(group_id, error)

        ! InAs valley offsets
        do i=1,301
            temp(i) = i-1
            v_offs(i,:) = varshni_model(INAS, dble(i))
            em_G(i) = effective_mass_model(INAS, dble(i))
        end do

        call h5gcreate_f(file_id, "inas", group_id, error)
        call h5screate_simple_f(2, (/int(301,8),int(3,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "inas/valley_offs", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, v_offs, (/int(301,8),int(3,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(301,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "inas/em_G", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, em_G, (/int(301,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(2,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "inas/em_L", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, (/me_l(INAS,L),me_t(INAS,L)/), (/int(2,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(2,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "inas/em_X", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, (/me_l(INAS,X),me_t(INAS,X)/), (/int(2,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5gclose_f(group_id, error)

        ! InGaAs valley offsets
        do i=1,301
        do j=1,101
            temp(i) = i-1
            conc(j) = (j-1)/100.0_s
            v_offs_a(i,j,:) = valley_offs_alloys(INGAAS,temp(i),GAAS,INAS,dble(conc(j)))
            me_G_a(i,j) = effective_mass_model_alloys(INGAAS,temp(i),GAAS,INAS,dble(conc(j)))
            me_LX(:,:,j) = effective_mass_model_alloys_LX(GAAS,INAS,dble(conc(j)))
        end do
        end do

        ! X-valley offset for different X-bowing parameter
        x_bow = (/0.0_s,0.05_s,0.1_s,0.15_s,0.2_s,0.25_s,0.3_s,0.35_s,0.4_s,0.45_s, &
              &  0.5_s,0.55_s,0.6_s,0.65_s,0.7_s,0.75_s,0.8_s,0.85_s,0.9_s,0.95_s, &
              & 1.0_s,1.05_s,1.1_s,1.15_s,1.2_s,1.25_s,1.3_s,1.35_s,1.4_s/)
        !do i=1,size(x_bow)
        !    print *, 'xbow', x_bow(i)
        !    vall_offs_bow(INGAAS,3) = x_bow(i)
        !    X_offs(i,:) = valley_offs_alloys(INGAAS,300.0_s,GAAS,INAS,0.53_s)
        !end do
        !call plot(conc,v_offs_a(301,:,1),conc,v_offs_a(301,:,2),conc,v_offs_a(301,:,3))

        call h5screate_simple_f(1, (/int(101,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "conc", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, conc, (/int(101,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5gcreate_f(file_id, "ingaas", group_id, error)
        call h5screate_simple_f(3, (/int(301,8),int(101,8),int(3,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "ingaas/valley_offs", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, v_offs_a, (/int(301,8),int(101,8),int(3,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(301,8),int(101,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "ingaas/em_G", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, me_G_a, (/int(301,8),int(101,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(2,8),int(101,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "ingaas/em_L", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, me_LX(2,:,:), (/int(2,8),int(101,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(2,8),int(101,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "ingaas/em_X", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, me_LX(3,:,:), (/int(2,8),int(101,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        !call h5screate_simple_f(2, (/int(29,8),int(3,8)/), dspace_id, error)
        !call h5dcreate_f(file_id, "ingaas/X_offs", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        !call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, X_offs, (/int(29,8),int(3,8)/), error)
        !call h5dclose_f(dset_id, error)
        !call h5sclose_f(dspace_id, error)

        call h5gclose_f(group_id, error)

        ! AlAs valley offsets
        do i=1,301
            temp(i) = i-1
            v_offs(i,:) = varshni_model(ALAS, dble(i))
            em_G(i) = effective_mass_model(ALAS, dble(i))
        end do

        call h5gcreate_f(file_id, "alas", group_id, error)
        call h5screate_simple_f(2, (/int(301,8),int(3,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "alas/valley_offs", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, v_offs, (/int(301,8),int(3,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(301,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "alas/em_G", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, em_G, (/int(301,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(2,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "alas/em_L", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, (/me_l(ALAS,L),me_t(ALAS,L)/), (/int(2,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(2,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "alas/em_X", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, (/me_l(ALAS,X),me_t(ALAS,X)/), (/int(2,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5gclose_f(group_id, error)

        ! InAlAs valley offsets
        do i=1,301
        do j=1,101
            temp(i) = i-1
            conc(j) = (j-1)/100.0_s
            v_offs_a(i,j,:) = valley_offs_alloys(INALAS,temp(i),ALAS,INAS,dble(conc(j)))
            me_G_a(i,j) = effective_mass_model_alloys(INALAS,temp(i),ALAS,INAS,dble(conc(j)))
            me_LX(:,:,j) = effective_mass_model_alloys_LX(ALAS,INAS,dble(conc(j)))
        end do
        end do
        !call plot(conc,v_offs_a(301,:,1),conc,v_offs_a(301,:,2),conc,v_offs_a(301,:,3))

        call h5gcreate_f(file_id, "inalas", group_id, error)
        call h5screate_simple_f(3, (/int(301,8),int(101,8),int(3,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "inalas/valley_offs", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, v_offs_a, (/int(301,8),int(301,8),int(3,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(301,8),int(101,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "inalas/em_G", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, me_G_a, (/int(301,8),int(101,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(2,8),int(101,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "inalas/em_L", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, me_LX(2,:,:), (/int(2,8),int(101,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(2,8),int(101,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "inalas/em_X", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, me_LX(3,:,:), (/int(2,8),int(101,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5gclose_f(group_id, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)
    end subroutine save_materialdef

    subroutine calc_mean_val(field_step,step,MAT)
        integer                               :: step, field_step, MAT, i, num
        real(kind=s),dimension(num_carriers)  :: mean_e_tmp
        real(kind=s)                          :: tmp
        real(kind=s)                          :: tmp_e

        !bulk_time(field_step,step) = time_tmp
        !print *, bulk_time(field_step,step)
        ! Get the average carrier velocity; save for each step in global array
        mean_vel(field_step,step) = sum(c(1:num_carriers)%vel(3))/num_carriers
        ! Anisotropic, [111]-direction
        !mean_vel(field_step,step) = sqrt((sum(c(1:num_carriers)%vel(1))/num_carriers)**2.0_s+ &
        !                          &  (sum(c(1:num_carriers)%vel(2))/num_carriers)**2.0_s+ &
        !                          &  (sum(c(1:num_carriers)%vel(1))/num_carriers)**2.0_s)
        tmp = 0
        num = 0
        tmp_e = 0
        do i=1,num_carriers
           if (c(i)%valley == X) then
               tmp = tmp+c(i)%vel(3)
               num = num+1
               tmp_e = tmp_e+c(i)%energy
           end if
        end do
        if (num > 0) then
            mean_energy(3,field_step,step) = tmp_e/num !+valley_offs(MAT,X)
            mean_vel_X(field_step,step) = tmp/num
            print *, char(9), char(9), 'Mean velocity in X-valley:', mean_vel_X(field_step,step)
        end if
        num = 0
        tmp = 0
        tmp_e = 0
        do i=1,num_carriers
           if ((c(i)%valley == L)) then! .and. (c(i)%v_subidx == 2)) then
               tmp = tmp+c(i)%vel(3)
               num = num+1
               tmp_e = tmp_e+c(i)%energy
           end if
        end do
        if (num > 0) then
            mean_energy(2,field_step,step) = tmp_e/num !+valley_offs(MAT,L)
            mean_vel_L(field_step,step) = tmp/num
            print *, char(9), char(9), 'Mean velocity in L-valley:', mean_vel_L(field_step,step)
        end if
        num = 0
        tmp = 0
        tmp_e = 0
        do i=1,num_carriers
           if ((c(i)%valley == G)) then! .and. (c(i)%v_subidx == 2)) then
               tmp = tmp+c(i)%vel(3)
               tmp_e = tmp_e+c(i)%energy
               num = num+1
           end if
        end do
        if (num > 0) then
            mean_energy(1,field_step,step) = tmp_e/num
            mean_vel_G(field_step,step) = tmp/num
            print *, char(9),char(9),'Mean velocity in G-valley:', mean_vel_G(field_step,step)
        end if

        ! Compute the mean energy, adding valley offsets for higher valleys
        where (c(1:num_carriers)%valley == G) mean_e_tmp = c(1:num_carriers)%energy
        where (c(1:num_carriers)%valley == L) mean_e_tmp = c(1:num_carriers)%energy+valley_offs(MAT,L)
        where (c(1:num_carriers)%valley == X) mean_e_tmp = c(1:num_carriers)%energy+valley_offs(MAT,X)
        mean_e(field_step,step) = sum(mean_e_tmp)/num_carriers
        !print *, sum(c(:)%vel(1))/num_carriers, sum(c(:)%vel(2))/num_carriers, sum(c(:)%vel(3))/num_carriers
        print *, char(9), 'Mean velocity:', mean_vel(field_step,step), 'Mean energy:', mean_e(field_step,step)
        !print *, char(9), 'Mean energy:', mean_energy(:,field_step,step)

        !print *, 'Mean energy:', sum(mean_e_tmp)/num_carriers
    end subroutine

    subroutine calc_valley_occ(field_step,step)
        integer                               :: step, field_step
        integer                               :: i, cg, cl, cx

        cg = 0
        cl = 0
        cx = 0
        do i=1,num_carriers
            if (c(i)%valley == G) then
                cg = cg+1
            else if (c(i)%valley == L) then
                cl = cl+1
            else if (c(i)%valley == X) then
                cx = cx+1
            end if
        end do
        valley_occ(field_step,step,:) = (/cg, cl, cx/)
        print *, char(9), '#carriers in valley: ', 'G:', cg, 'L:', cl, 'X:', cx
    end subroutine calc_valley_occ

    subroutine save_scatrates_hdf5()
        integer(HID_T)          :: file_id       ! File identifier
        integer(HID_T)          :: dset_id       ! Dataset identifier
        integer(HID_T)          :: dspace_id     ! Dataspace identifier
        integer(HID_T)          :: group_id     ! Dataspace identifier
        CHARACTER(LEN=40)       :: filename = "/tmp/gaas_scatrates.h5" ! File name
        character(len=7)        :: groupname
        integer                 :: error, i, j, k
        real(kind=s),dimension(int(maxval(scat_m)),3,size_E/100,15)    :: scat_rates_tmp
        real(kind=s),dimension(size_E)             :: e_val

        do i=1,size_E/100
            e_val(i) = E(i*100)
        end do

        do k=1,int(maxval(scat_m))
        do j=1,3
        do i=1,size_E/100
            scat_rates_tmp(k,j,i,:) = scat_rates(k,j,i*100,:)
        end do
        end do
        end do

        !filename = scatrate_filename(1:scatrate_filename_len)
        !print *, filename

        call h5open_f(error)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

        call h5screate_simple_f(1, (/int(size_E/100,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "energy", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, e_val, (/int(size_E/100,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)


        do k=1,int(maxval(scat_m))
        write (groupname, "(A6,I1)") "region", k
        call h5gcreate_f(file_id, groupname, group_id, error)
        !call h5gopen_f(file_id, groupname, group_id, error)
        !print *, groupname//"/G"
!        call h5gcreate_f(file_id, "region1", group_id, error)
!        call h5gopen_f(file_id, "region1", group_id, error)
!        print *, scat_rates_tmp(G,5,:,k)
        call h5screate_simple_f(2, (/int(size_E/100,8),int(15,8)/), dspace_id, error)
        call h5dcreate_f(file_id, groupname//"/G", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, scat_rates_tmp(k,G,:,:), (/int(size_E/100,8),int(15,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
        call h5screate_simple_f(2, (/int(size_E/100,8),int(15,8)/), dspace_id, error)
        call h5dcreate_f(file_id, groupname//"/L", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, scat_rates_tmp(k,L,:,:), (/int(size_E/100,8),int(15,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
        call h5screate_simple_f(2, (/int(size_E/100,8),int(15,8)/), dspace_id, error)
        call h5dcreate_f(file_id, groupname//"/X", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, scat_rates_tmp(k,X,:,:), (/int(size_E/100,8),int(15,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
        call h5gclose_f(group_id, error)
        end do

        call h5fclose_f(file_id, error)
        call h5close_f(error)
    end subroutine save_scatrates_hdf5

    subroutine save_bulk_hdf5()
        integer(HID_T)          :: file_id       ! File identifier
        integer(HID_T)          :: dset_id       ! Dataset identifier
        integer(HID_T)          :: dspace_id     ! Dataspace identifier
        integer(HID_T)          :: group_id      ! Group identifier
        CHARACTER(LEN=40)       :: filename = "/tmp/gaas.h5" ! File name
        integer                 :: error

        filename = bulk_filename(1:bulk_filename_len)

        call h5open_f(error)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

        call h5gcreate_f(file_id, "valley_occ", group_id, error)
        call h5gopen_f(file_id, "valley_occ", group_id, error)
        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "valley_occ/G", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, valley_occ(1:k_max,1:n_step,G), &
           & (/int(k_max,8),int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "valley_occ/L", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, valley_occ(1:k_max,1:n_step,L), &
           & (/int(k_max,8),int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "valley_occ/X", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, valley_occ(1:k_max,1:n_step,X), &
           & (/int(k_max,8),int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)
        call h5gclose_f(group_id, error)

        call h5gcreate_f(file_id, "energy", group_id, error)
        call h5gopen_f(file_id, "energy", group_id, error)
        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "energy/G", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_energy(G,1:k_max,1:n_step), &
           & (/int(k_max,8),int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "energy/L", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_energy(L,1:k_max,1:n_step), &
           & (/int(k_max,8),int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "energy/X", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_energy(X,1:k_max,1:n_step), &
           & (/int(k_max,8),int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)
        call h5gclose_f(group_id, error)

        call h5gcreate_f(file_id, "velocity", group_id, error)
        call h5gopen_f(file_id, "velocity", group_id, error)
        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "velocity/G", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_vel_G(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "velocity/L", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_vel_L(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "velocity/X", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_vel_X(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "vel", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_vel(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5dclose_f(dset_id, error)
        !call h5dcreate_f(file_id, "energy", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        !call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_e(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        !call h5dclose_f(dset_id, error)
        print *, bulk_time(1,:)
        call h5dcreate_f(file_id, "time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, bulk_time(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(k_max,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_field", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, k_array(1:k_max), (/int(k_max,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5gcreate_f(file_id, "el_temp", group_id, error)
        call h5gopen_f(file_id, "el_temp", group_id, error)
        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_temp/G", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_el_temp_G(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_temp/L", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_el_temp_L(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_temp/X", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_el_temp_X(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5gcreate_f(file_id, "el_fermi", group_id, error)
        call h5gopen_f(file_id, "el_fermi", group_id, error)
        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_fermi/G", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_el_fermi_G(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_fermi/L", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_el_fermi_L(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5screate_simple_f(2, (/int(k_max,8),int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_fermi/X", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mean_el_fermi_X(1:k_max,1:n_step), (/int(k_max,8), int(n_step,8)/), error)
        call h5sclose_f(dspace_id, error)
        call h5dclose_f(dset_id, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)
    end subroutine save_bulk_hdf5

#if DIM == 0
    subroutine calc_pauli()
        integer                     :: i, mat
        integer,dimension(3)        :: valley
        real(kind=s)                :: l_nonp, l_me
        real(kind=s),dimension(3)   :: k_tmp, E_tmp, E_total

        p%k_drift = 0.0_s
        p%n_conc = 0.0_s
        valley = 0
        E_total = 0.0_s

        ! Determine k_drift
        do i=1,num_carriers
            valley(c(i)%valley) = valley(c(i)%valley)+1
            p%k_drift(c(i)%valley,:) = p%k_drift(c(i)%valley,:)+c(i)%k
        end do
        p%k_drift(1,:) = p%k_drift(1,:)/valley(1)
        p%k_drift(2,:) = p%k_drift(2,:)/valley(2)
        p%k_drift(3,:) = p%k_drift(3,:)/valley(3)

        ! Compute <E(k-k_d)>
        do i=1,num_carriers
            mat = m_material(1)
            l_nonp = nonp(mat,c(i)%valley)
            l_me = me(mat,c(i)%valley)
            k_tmp = c(i)%k-p%k_drift(c(i)%valley,:)
            if ((c(i)%valley == X) .or. (c(i)%valley == L)) then
            E_tmp(c(i)%valley) = -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(k_tmp,k_tmp)/q)/ &
                            & (2.0_s*m0*l_nonp))
            else
            E_tmp(c(i)%valley) = -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(k_tmp,k_tmp)/q)/ &
                            & (2.0_s*l_me*l_nonp))
            end if
            E_total(c(i)%valley) = E_total(c(i)%valley)+E_tmp(c(i)%valley)
        end do
        E_total = E_total/valley
        p%E_avg = E_total
        p%n_conc = (donor_density(1)/num_carriers)*valley
        !print *, p%k_drift(1,:), E_total, valley, p%n_conc
!    print *, -1.0_s/(2.0_s*l_nonp)+sqrt(1.0_s/(4.0_s*(l_nonp**2.0_s))+(hb**2.0_s*dot_product(p%k_drift(1,:),p%k_drift(1,:))/q)/ &
!                            & (2.0_s*l_me*l_nonp))
    end subroutine calc_pauli

    subroutine save_pauli(field_step,step)
        integer             :: field_step, step

        mean_el_fermi_G(field_step,step) = p%fermi_lvl(1)
        mean_el_fermi_L(field_step,step) = p%fermi_lvl(2)
        mean_el_fermi_X(field_step,step) = p%fermi_lvl(3)

        mean_el_temp_G(field_step,step) = p%el_temp(1)
        mean_el_temp_L(field_step,step) = p%el_temp(2)
        mean_el_temp_X(field_step,step) = p%el_temp(3)
    end subroutine save_pauli
#endif

#if DIM == 1
    subroutine save_device()
        integer(HID_T)          :: file_id       ! File identifier
        integer(HID_T)          :: dset_id       ! Dataset identifier
        integer(HID_T)          :: dspace_id     ! Dataspace identifier
        integer(HID_T)          :: group_id      ! Group identifier
        CHARACTER(LEN=25)       :: filename = "/home/ss/Python/device.h5" ! File name
        integer                 :: error, i
        real(kind=s),dimension(:),allocatable   :: time

        filename = dev2_filename(1:dev2_filename_len)

        call h5open_f(error)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

        call h5gcreate_f(file_id, "setup", group_id, error)
        call h5gopen_f(file_id, "setup", group_id, error)

        call h5screate_simple_f(1, (/int(sum(num_nodes),8)/), dspace_id, error)
        call h5dcreate_f(file_id, "setup/nodes", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%nodes), (/int(sum(num_nodes),8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        allocate(time(n_step))
        do i=1,int(n_step,4)
            time(i) = i*t_step
        end do

        call h5screate_simple_f(1, (/int(n_step,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "setup/time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(time), (/int(n_step,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(n_step,8),int(sum(num_nodes),8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_conc", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_el_conc), (/int(n_step,8),int(sum(num_nodes),8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(n_step,8),int(sum(num_nodes),8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_pot", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_el_pot), (/int(n_step,8),int(sum(num_nodes),8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(n_step,8),int(sum(num_nodes),8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_field", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_el_field), (/int(n_step,8),int(sum(num_nodes),8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(n_step,8),int(2,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "particles_out", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dev%particle_flow_out, (/int(n_step,8),int(2,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(n_step,8),int(2,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "particles_in", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dev%particle_flow_in, (/int(n_step,8),int(2,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(n_step,8),int(sum(num_nodes),8)/), dspace_id, error)
        call h5dcreate_f(file_id, "number_carriers", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dev%g_num, (/int(n_step,8),int(sum(num_nodes),8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(n_step,8),int(sum(num_nodes),8)/), dspace_id, error)
        call h5dcreate_f(file_id, "velocity", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_vel), (/int(n_step,8),int(sum(num_nodes),8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(2, (/int(n_step,8),int(sum(num_nodes),8)/), dspace_id, error)
        call h5dcreate_f(file_id, "energy", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_energy), (/int(n_step,8),int(sum(num_nodes),8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)
    end subroutine save_device!        character(len=1)        :: dname

#endif

#if DIM == 2
!    subroutine save_device()
!        integer(HID_T)          :: file_id       ! File identifier
!        integer(HID_T)          :: dset_id       ! Dataset identifier
!        integer(HID_T)          :: dspace_id     ! Dataspace identifier
!        integer(HID_T)          :: group_id      ! Group identifier
!        CHARACTER(LEN=25)       :: filename = "/home/ss/Python/device.h5" ! File name
!        integer                 :: error, i, j
!        character(len=1)        :: dname
!        real(kind=s),dimension(:),allocatable   :: time
!        real(kind=s),dimension(:,:),allocatable :: cb_offs
!
!        allocate(cb_offs(num_nodes(1),num_nodes(2)))
!        do i=1,num_nodes(1)
!            do j=1,num_nodes(2)
!                cb_offs(i,j) = emin(m_material(m_region(i,j)))
!            end do
!        end do
!
!        filename = dev2_filename(1:dev2_filename_len)
!
!        call h5open_f(error)
!        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
!
!        call h5gcreate_f(file_id, "setup", group_id, error)
!        call h5gopen_f(file_id, "setup", group_id, error)
!
!        !call h5screate_simple_f(2, (/int(num_cells,8),int(3,8)/), dspace_id, error)
!        !call h5dcreate_f(file_id, "setup/idx2sub", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
!        !call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, A_idx2sub, (/int(num_cells,8),int(3,8)/), error)
!        !call h5dclose_f(dset_id, error)
!        !call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(1, (/int(l_save_x,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "setup/node_x", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(node_coord(1,save_x(1):save_x(2))), (/int(l_save_x,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(1, (/int(l_save_y,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "setup/node_y", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(node_coord(2,save_y(1):save_y(2))), (/int(l_save_y,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(1, (/int(n_step/save_t,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "setup/time", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%time), (/int(n_step/save_t,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(num_nodes(1),8),int(num_nodes(2),8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "setup/cb_offs", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(cb_offs), (/int(num_nodes(1),8),int(num_nodes(2),8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        if (tunnel_hdf5 == 1) then
!        call h5gcreate_f(file_id, "tunnel", group_id, error)
!        call h5gopen_f(file_id, "tunnel", group_id, error)
!
!        do i=1,n_tunnel
!        !write (groupname, "(A6,I1)") "region", k
!        write(dname, '(i1)' ) i
!        call h5gcreate_f(file_id, "tunnel/"//dname, group_id, error)
!        call h5gopen_f(file_id, "tunnel/"//dname, group_id, error)
!        print *, "tunnel/"//dname//"/prob"
!        call h5screate_simple_f(2, (/int(t_sim/0.02+1,8),int(1000,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "tunnel/"//dname//"/G_prob", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%tun_prob(i,:,:,G)), (/int(t_sim/0.02+1,8),int(1000,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(t_sim/0.02+1,8),int(1000,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "tunnel/"//dname//"/G_energy", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%tun_energy(i,:,:,G)), (/int(t_sim/0.02+1,8),int(1000,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(t_sim/0.02+1,8),int(1000,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "tunnel/"//dname//"/L_prob", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%tun_prob(i,:,:,L)), (/int(t_sim/0.02+1,8),int(1000,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(t_sim/0.02+1,8),int(1000,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "tunnel/"//dname//"/L_energy", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%tun_energy(i,:,:,L)), (/int(t_sim/0.02+1,8),int(1000,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(t_sim/0.02+1,8),int(1000,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "tunnel/"//dname//"/X_prob", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%tun_prob(i,:,:,X)), (/int(t_sim/0.02+1,8),int(1000,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(t_sim/0.02+1,8),int(1000,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "tunnel/"//dname//"/X_energy", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%tun_energy(i,:,:,X)), (/int(t_sim/0.02+1,8),int(1000,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(t_sim/0.02+1,8),int(1000,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "tunnel/"//dname//"/dest", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%tun_dest(i,:,:)), (/int(t_sim/0.02+1,8),int(1000,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        end do
!        end if
!        !allocate(time(n_step))
!        !do i=1,int(n_step,4)
!        !    time(i) = i*t_step
!        !end do
!
!        !call h5screate_simple_f(1, (/int(n_step,8)/), dspace_id, error)
!        !call h5dcreate_f(file_id, "setup/time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
!        !call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time, (/int(n_step,8)/), error)
!        !call h5dclose_f(dset_id, error)
!        !call h5sclose_f(dspace_id, error)int(n_step/save_t)
!
!        call h5screate_simple_f(3, (/int(n_step/save_t,8),int(l_save_x,8),int(l_save_y,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "el_conc", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_el_conc), (/int(n_step/save_t,8), &
!                                   & int(l_save_x,8),int(l_save_y,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(3, (/int(n_step/save_t,8),int(l_save_x,8),int(l_save_y,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "el_pot", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_el_pot), (/int(n_step/save_t,8), &
!                                    & int(l_save_x,8),int(l_save_y,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(3, (/int(n_step/save_t,8),int(l_save_x,8),int(l_save_y,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "energy", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_energy), (/int(n_step/save_t,8), &
!                                    & int(l_save_x,8),int(l_save_y,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(3, (/int(n_step/save_t,8),int(l_save_x,8),int(l_save_y,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "number_carriers", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dev%g_num, (/int(n_step/save_t,8), &
!                                    & int(l_save_x,8),int(l_save_y,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(4, (/int(n_step/save_t,8),int(2,8),int(l_save_x,8),int(l_save_y,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "el_field", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_el_field), (/int(n_step/save_t,8),int(2,8), &
!                                    & int(l_save_x,8),int(l_save_y,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(4, (/int(n_step/save_t,8),int(2,8),int(l_save_x,8),int(l_save_y,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "velocity", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_vel), (/int(n_step/save_t,8),int(2,8), &
!                                   & int(l_save_x,8),int(l_save_y,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        !call h5screate_simple_f(1, (/int(1,8)/), dspace_id, error)
!        !call h5dcreate_f(file_id, "force", H5T_NATIVE_REAL, dspace_id, dset_id, error)
!        !call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(sqrt(dot_product(c(2)%loc_e_field,c(2)%loc_e_field))), (/int(1,8)/), error)
!        !call h5dclose_f(dset_id, error)
!        !call h5sclose_f(dspace_id, error)
!!
!        call h5screate_simple_f(2, (/int(n_step/save_t,8),int(n_out,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "particles_out", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dev%particle_flow_out, (/int(n_step/save_t,8),int(n_out,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(n_step/save_t,8),int(n_out,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "particles_in", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dev%particle_flow_in, (/int(n_step/save_t,8),int(n_out,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!!
!!        call h5screate_simple_f(2, (/int(n_step,8),int(sum(num_nodes),8)/), dspace_id, error)
!!        call h5dcreate_f(file_id, "vel", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
!!        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dev%g_vel, (/int(n_step,8),int(sum(num_nodes),8)/), error)
!!        call h5dclose_f(dset_id, error)
!!        call h5sclose_f(dspace_id, error)
!!
!!        call h5screate_simple_f(2, (/int(n_step,8),int(sum(num_nodes),8)/), dspace_id, error)
!!        call h5dcreate_f(file_id, "energy", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
!!        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dev%g_energy, (/int(n_step,8),int(sum(num_nodes),8)/), error)
!!        call h5dclose_f(dset_id, error)
!!        call h5sclose_f(dspace_id, error)
!
!        call h5fclose_f(file_id, error)
!        call h5close_f(error)
!    end subroutine save_device
#endif

#if DIM == 3
    subroutine save_device()
        integer(HID_T)          :: file_id       ! File identifier
        integer(HID_T)          :: dset_id       ! Dataset identifier
        integer(HID_T)          :: dspace_id     ! Dataspace identifier
        integer(HID_T)          :: group_id      ! Group identifier
        CHARACTER(LEN=25)       :: filename = "/home/ss/Python/device.h5" ! File name
        integer                 :: error!, step !, i
        real(kind=s),dimension(:),allocatable   :: time

        !step = 1

        filename = dev2_filename(1:dev2_filename_len)

        call h5open_f(error)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

        call h5gcreate_f(file_id, "setup", group_id, error)
        call h5gopen_f(file_id, "setup", group_id, error)

        !call h5screate_simple_f(2, (/int(sum(num_cells),8),int(3,8)/), dspace_id, error)
        !call h5dcreate_f(file_id, "setup/idx2sub", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
        !call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, A_idx2sub, (/int(sum(num_cells),8),int(3,8)/), error)
        !call h5dclose_f(dset_id, error)
        !call h5sclose_f(dspace_id, error)

        !call h5screate_simple_f(1, (/int(sum(num_nodes),8)/), dspace_id, error)
        !call h5dcreate_f(file_id, "setup/nodes", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        !call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dev%nodes, (/int(sum(num_nodes),8)/), error)
        !call h5dclose_f(dset_id, error)
        !call h5sclose_f(dspace_id, error)

        !allocate(time(n_step))
        !do i=1,int(n_step,4)
        !    time(i) = i*t_step
        !end do

        !call h5screate_simple_f(1, (/int(n_step,8)/), dspace_id, error)
        !call h5dcreate_f(file_id, "setup/time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        !call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time, (/int(n_step,8)/), error)
        !call h5dclose_f(dset_id, error)
        !call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(4, (/int(n_step,8),int(num_nodes(1),8),int(num_nodes(2),8), &
                            & int(num_nodes(3),8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_conc", H5T_NATIVE_REAL, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_el_conc), (/int(n_step,8),int(num_nodes(1),8),int(num_nodes(2),8), &
                            & int(num_nodes(3),8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
!
        call h5screate_simple_f(4, (/int(n_step,8),int(num_nodes(1),8),int(num_nodes(2),8), &
                            & int(num_nodes(3),8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_pot", H5T_NATIVE_REAL, dspace_id, dset_id, error)
        !call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dev%g_el_pot, (/int(n_step,8),int(sum(num_cells),8)/), error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_el_pot), (/int(n_step,8),int(num_nodes(1),8),int(num_nodes(2),8), &
                            & int(num_nodes(3),8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(5, (/int(n_step,8),int(3,8),int(num_nodes(1),8),int(num_nodes(2),8), &
                            & int(num_nodes(3),8)/), dspace_id, error)
        call h5dcreate_f(file_id, "el_field", H5T_NATIVE_REAL, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(dev%g_el_field), (/int(n_step,8),int(3,8),int(num_nodes(1),8), &
                            & int(num_nodes(2),8),int(num_nodes(3),8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        call h5screate_simple_f(1, (/int(1,8)/), dspace_id, error)
        call h5dcreate_f(file_id, "force", H5T_NATIVE_REAL, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sngl(sqrt(dot_product(c(2)%loc_e_field,c(2)%loc_e_field))), (/int(1,8)/), error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(n_step,8),int(2,8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "particles_out", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dev%particle_flow_out, (/int(n_step,8),int(2,8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(n_step,8),int(sum(num_nodes),8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "vel", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dev%g_vel, (/int(n_step,8),int(sum(num_nodes),8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)
!
!        call h5screate_simple_f(2, (/int(n_step,8),int(sum(num_nodes),8)/), dspace_id, error)
!        call h5dcreate_f(file_id, "energy", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
!        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dev%g_energy, (/int(n_step,8),int(sum(num_nodes),8)/), error)
!        call h5dclose_f(dset_id, error)
!        call h5sclose_f(dspace_id, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)
    end subroutine save_device
#endif
end module show_output
