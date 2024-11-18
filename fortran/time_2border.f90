module time_2border
    use config

    implicit none
    contains

#if DIM == 1
    subroutine time2r_1(t_real, r_idx, r, v, idx, cidx)
        real(kind=s),dimension(3)   :: r, v
        integer                     :: idx, idx_l, idx_r, cidx
        real(kind=s)                :: t_real, t1, t2, t3, t4
        real(kind=s),dimension(4)   :: t_arr
        integer,dimension(1)        :: r_idx

        t1 = 0.0_s
        t2 = 0.0_s
        t3 = 0.0_s
        t4 = 0.0_s
        idx_l = idx-1
        if (idx <= 1) then
            idx_l = 1
        end if
        idx_r = idx+1
        if (idx >= sum(num_nodes)) then
            idx_r = sum(num_nodes)
        end if
        !print *, r_idx, r, v, idx, cidx
        ! Time to positive coordinate border
        !print *, cidx, v, r
        t1 = ((node_coord(idx)+node_dist(2,idx))-r(1))/v(1)
        !t2 = ((node_coord(idx)+3.0_s/2.0_s*m_dx(1))-r(1))/v(1)
        !print *, cidx, t2, ((node_coord(idx)+3.0_s/2.0_s*m_dx(1))-r(1))/v(1)

        t2 = ((node_coord(idx)+node_dist(2,idx)+2.0_s*node_dist(2,idx_r))-r(1))/v(1)
        !if ((r(1) .gt. (dev_x-node_dist(1,sum(num_nodes))-1e-5_s)) .and. (contact_type(2) .ne. INSULATOR)) then
        !    t1 = (dev_x-r(1))/v(1)
        !    t2 = t1
        !end if
        ! Time to negative coordinate border
        t3 = ((node_coord(idx)-node_dist(1,idx))-r(1))/v(1)
        t4 = ((node_coord(idx)-node_dist(1,idx)-2.0_s*node_dist(1,idx_l))-r(1))/v(1)
        !if ((r(1) .lt. (node_dist(2,1)+1e-5_s)) .and. (contact_type(1) .ne. INSULATOR)) then
        !    t3 = -r(1)/v(1)
        !    t4 = t3
        !end if

        !print *, cidx, r(1), r(1)+((node_coord(idx)+3.0_s/2.0_s*m_dx(1))-r(1)), (node_coord(idx)+3.0_s/2.0_s*m_dx(1))
        !print *, cidx, r(1), t1, t2, t3, t4, r(1)+t1*v(1), r(1)+t2*v(1), r(1)+t3*v(1), r(1)+t4*v(1), v(1)
        !print *, cidx,
        t_arr = (/t1,t2,t3,t4/)
        t_real = minval(t_arr, MASK = (t_arr .gt. 0) .and. (t_arr .gt. 1e-15_s))
        r_idx = minloc(t_arr, MASK = (t_arr .gt. 0) .and. (t_arr .gt. 1e-15_s))
    end subroutine time2r_1

    subroutine get_r_idx(r, idx)
        real(kind=s),dimension(3),intent(in)   :: r
        integer,dimension(1)        :: r_idx
        integer                     :: i
        integer,intent(out)         :: idx

        r_idx = minloc(abs(r(1)-node_coord))
        !if (((node_coord(r_idx(1)-1)+node_dist(1,r_idx(1)-1)) .gt. r(1)) .and. (r_idx(1) .gt. 1)) then
        !    r_idx = r_idx-1
        !end if
        !if (((node_coord(r_idx(1)+1)-node_dist(2,r_idx(1)+1)) .lt. r(1)) .and. (r_idx(1) .lt. sum(num_nodes))) then
        !    r_idx = r_idx+1
        !end if
        !print *, r_idx, node_coord(r_idx-1),node_coord(r_idx),node_coord(r_idx+1)
        idx = r_idx(1)
        !print *, r(1), idx, r(1)-node_coord(idx-1), r(1)-node_coord(idx), r(1)-node_coord(idx+1)
    end subroutine get_r_idx

    subroutine get_efield_idx(r, v, idx)
        real(kind=s),dimension(3),intent(in)   :: r, v
        integer,dimension(1)        :: r_idx
        integer,intent(out)         :: idx

        r_idx = minloc(abs(r(1)-node_coord))
        if (((node_coord(r_idx(1)-1)+node_dist(1,r_idx(1)-1)) .gt. r(1)) .and. (r_idx(1) .gt. 1)) then
            r_idx = r_idx-1
        end if
        if (((node_coord(r_idx(1)+1)-node_dist(2,r_idx(1)+1)) .lt. r(1)) .and. (r_idx(1) .lt. sum(num_nodes))) then
            r_idx = r_idx+1
        end if
        !print *, r_idx, node_coord(r_idx-1),node_coord(r_idx),node_coord(r_idx+1)
        idx = r_idx(1)
        !print *, r(1), idx, r(1)-node_coord(idx-1), r(1)-node_coord(idx), r(1)-node_coord(idx+1)
    end subroutine get_efield_idx
#endif

#if DIM == 2
    subroutine time2r_2(t_real,b_idx,r,r_idx,vel)
        real(kind=s),dimension(3)     :: r, vel
        integer,dimension(2)          :: r_idx
        real(kind=s),dimension(8)     :: t
        real(kind=s)                  :: t_real
        integer,dimension(1)          :: b_idx
!        real(kind=s)                  :: i

        t = 0.0_s
        !i = 1.0_s

!        ! See notes, x lines
!        t(1) = (node_coord(1,r_idx(1))+m_dx/2.0_s-r(1))/vel(1)
!        t(2) = (node_coord(1,r_idx(1))+3.0_s/2.0_s*m_dx-r(1))/vel(1)
!        t(5) = (node_coord(1,r_idx(1))-m_dx/2.0_s-r(1))/vel(1)
!        t(6) = (node_coord(1,r_idx(1))-3.0_s/2.0_s*m_dx-r(1))/vel(1)
!        ! y lines
!        t(3) = (node_coord(2,r_idx(2))+m_dy/2.0_s-r(2))/vel(2)
!        t(4) = (node_coord(2,r_idx(2))+3.0_s/2.0_s*m_dy-r(2))/vel(2)
!        t(7) = (node_coord(2,r_idx(2))-m_dy/2.0_s-r(2))/vel(2)
!        t(8) = (node_coord(2,r_idx(2))-3.0_s/2.0_s*m_dy-r(2))/vel(2)
!
!        !if ((r(1) <= 40.0_s) .and. (r(2) > 239.0_s)) then!(contact_type(1,20) == OHMIC)) then
!        !    t(5) = -r(1)/vel(1)
!        !    t(6) = t(5)
!        !end if
!        !if ((r(2) >= (dev_y-m_ly(num_nodes(2))))) then! .and. (r(1) < 40.0_s)) then!(contact_type(num_nodes(1),20) == OHMIC)) then
!        !    t(3) = (dev_y-r(2))/vel(2)
!        !    t(4) = t(3)
!        !end if
!!        if ((r(1) <= m_lx(1)) .and. (r(2) > 239.0_s)) then!(contact_type(1,20) == OHMIC)) then
!!            t(5) = -r(1)/vel(1)
!!            t(6) = t(5)
!!        end if
!!        if ((r(1) >= (dev_x-m_lx(num_nodes(1)))) .and. (r(2) > 199.0_s)) then!(contact_type(num_nodes(1),20) == OHMIC)) then
!!            t(1) = (dev_x-r(1))/vel(1)
!!            t(2) = t(1)
!!        end if
!        !print *, t, r_idx, m_dx, vel(1:2)
!
!        t_real = minval(t, MASK = (t .gt. 0) .and. (t .gt. 1e-10_s))
!        b_idx = minloc(t, MASK = (t .gt. 0) .and. (t .gt. 1e-10_s))


!        do while ((t(1) <= 1e-10_s) .and. (t(3) <= 1e-10_s) .and. (t(5) <= 1e-10_s) .and. (t(7) <= 1e-10_s))
            !print *, i, t(1), t(5), t(3), t(7)
           ! if (abs(vel(1)) < 1e-4_s) then
          !      print *, 'asdf', vel(1:2), r
          !  end if
            t(1) = (node_coord(1,r_idx(1))+m_dx_half-r(1))/vel(1)
            t(5) = (node_coord(1,r_idx(1))-m_dx_half-r(1))/vel(1)
            t(3) = (node_coord(2,r_idx(2))+m_dy_half-r(2))/vel(2)
            t(7) = (node_coord(2,r_idx(2))-m_dy_half-r(2))/vel(2)
!if (i > 1_s) then
!    print *, r(1:2), vel(1:2), t(1), t(5), t(3), t(7)
!end if
!            i = i+2.0_s
!        end do

!        print *, node_coord(2,r_idx(2))
!        print *, r(1), t(1), t(2), t(5), t(6), r(1)+t(1)*vel(1), r(1)+t(2)*vel(1), r(1)+t(5)*vel(1), r(1)+t(6)*vel(1), vel(1)
     !   print *, r(1), t(1), t(5), r(1)+t(1)*vel(1), r(1)+t(5)*vel(1), vel(1)
     !   print *, r(2), t(3), t(7), r(2)+t(3)*vel(2), r(2)+t(7)*vel(2), vel(2)

!        print *, r(2), t(3), t(4), t(7), t(8), r(2)+t(3)*vel(2), r(2)+t(4)*vel(2), r(2)+t(7)*vel(2), r(2)+t(8)*vel(2), vel(2)
!        print *, t_real, b_idx
        t_real = minval(t, MASK = (t .gt. 0))! .and. (t .gt. 1e-10_s))
        b_idx = minloc(t, MASK = (t .gt. 0))! .and. (t .gt. 1e-10_s))
!        if ((b_idx(1) == 2) .or. (b_idx(1) == 4) .or. (b_idx(1) == 6) .or. (b_idx(1) == 8)) then
!            print *, 'time 2border', t, b_idx, r
!        end if
    end subroutine time2r_2

    subroutine get_r_idx_2(r, idx)
        real(kind=s),dimension(3),intent(in)     :: r
!        integer,dimension(1)                     :: r_idx
        integer,dimension(2),intent(out)         :: idx
!        integer,dimension(2)            :: idx_new

!        integer                         :: i

        !r_idx = minloc(abs(r(1)-node_coord(1,1:num_nodes(1))))
        !idx(1) = r_idx(1)
        !r_idx = minloc(abs(r(2)-node_coord(2,1:num_nodes(2))))
        !idx(2) = r_idx(1)
        idx(1) = int((r(1)+m_dx_half)/m_dx)+1
        idx(2) = int((r(2)+m_dy_half)/m_dy)+1

        if (idx(1)>num_nodes(1)) then
            idx(1) = num_nodes(1)
        end if
        if (idx(2)>num_nodes(2)) then
            idx(2) = num_nodes(2)
        end if
        if (idx(1)<1) then
            idx(1) = 1
        end if
        if (idx(2)<1) then
            idx(2) = 1
        end if

        !if ((idx(1) .ne. idx_new(1)) .or. (idx(2) .ne. idx_new(2))) then
        !    print *, r, idx, idx_new
        !end if
    end subroutine get_r_idx_2
#endif

#if DIM == 3
    subroutine time2r_3(t_real,b_idx,r,r_idx,vel)
        real(kind=s),dimension(3)     :: r, vel
        integer,dimension(3)          :: r_idx
        real(kind=s),dimension(12)    :: t
        real(kind=s)                  :: t_real
        integer,dimension(1)          :: b_idx

        ! See notes, xy-planes, z-direction
        t(1) = (node_coord(3,r_idx(3))+m_dx(1)/2.0_s-r(3))/vel(3)
        t(2) = (node_coord(3,r_idx(3))+3.0_s/2.0_s*m_dx(1)-r(3))/vel(3)
        t(3) = (node_coord(3,r_idx(3))-m_dx(1)/2.0_s-r(3))/vel(3)
        t(4) = (node_coord(3,r_idx(3))-3.0_s/2.0_s*m_dx(1)-r(3))/vel(3)
        ! xz-planes, y-direction
        t(5) = (node_coord(2,r_idx(2))+m_dx(1)/2.0_s-r(2))/vel(2)
        t(6) = (node_coord(2,r_idx(2))+3.0_s/2.0_s*m_dx(1)-r(2))/vel(2)
        t(7) = (node_coord(2,r_idx(2))-m_dx(1)/2.0_s-r(2))/vel(2)
        t(8) = (node_coord(2,r_idx(2))-3.0_s/2.0_s*m_dx(1)-r(2))/vel(2)
        ! yz-planes, x-direction
        t(9) = (node_coord(1,r_idx(1))+m_dx(1)/2.0_s-r(1))/vel(1)
        t(10) = (node_coord(1,r_idx(1))+3.0_s/2.0_s*m_dx(1)-r(1))/vel(1)
        t(11) = (node_coord(1,r_idx(1))-m_dx(1)/2.0_s-r(1))/vel(1)
        t(12) = (node_coord(1,r_idx(1))-3.0_s/2.0_s*m_dx(1)-r(1))/vel(1)

        if ((r(1) <= m_lx(1)) .and. (contact_type(LEFT,10,10) == OHMIC)) then
            t(11) = -r(1)/vel(1)
            t(12) = t(11)
        end if
!        if ((r(1) >= (dev_x-m_lx(num_nodes(1)))) .and. (contact_type(num_nodes(1),20) == OHMIC)) then
!            t(1) = (dev_x-r(1))/vel(1)
!            t(2) = t(1)
!        end if
        t_real = minval(t, MASK = (t .gt. 0) .and. (t .gt. 1e-15_s))
        b_idx = minloc(t, MASK = (t .gt. 0) .and. (t .gt. 1e-15_s))
        !print *, r, r_idx, r(1)+t(9)*vel(1),r(1)+t(10)*vel(1),r(1)+t(11)*vel(1),r(1)+t(12)*vel(1)
        !print *, t, t_real, b_idx
    end subroutine time2r_3

    subroutine get_r_idx_3(r, idx)
        real(kind=s),dimension(3),intent(in)     :: r
        integer,dimension(1)                     :: r_idx
        integer,dimension(3),intent(out)         :: idx

        r_idx = minloc(abs(r(1)-node_coord(1,1:num_nodes(1))))
        idx(1) = r_idx(1)
        r_idx = minloc(abs(r(2)-node_coord(2,1:num_nodes(2))))
        idx(2) = r_idx(1)
        r_idx = minloc(abs(r(3)-node_coord(3,1:num_nodes(3))))
        idx(3) = r_idx(1)
    end subroutine get_r_idx_3
#endif
end module time_2border
