!!!!! Diffusion near arbitrary Koch
!!!!! Reflection after first arrval
!!!!! Yilin YE @ 05/2023


Program main
    use MATHS
    implicit none
    include "para.h"
    integer :: i,j,k,l,jj ! Loop variables
    integer :: outpoint(2), intpoint(3), out_gap, int_gap, rec_old1, rec_old2
    integer :: qmax, g_act, g_test, case, infvl, supvl, search_length, k1, k2
    integer :: cross, inside, touch !! Variables for Reflection
    integer :: traj_max
    real*8 :: angle_ini, tta, new_x, new_y, radius_min, radius, u, l1, l2, ell, ell_hat, time_diff
    real*8 :: avg_seuil, avg_local, avg_temps
    real*8,allocatable :: position(:,:), coord_x(:), coord_y(:), account(:), zeta(:)
    integer,allocatable :: rec_index(:,:), search_index(:) !! Arrays to record indices
    integer,external:: noloop !! Function to return reasonnable segment index
    real*8,external :: span !! Function to calculate distance between two points.
	
    integer :: date(8)
    real*8 :: temps_debut, temps_fin, temps_provisoire!, time_begin, time_end !! Define variables to record time consumed
    character*128 :: rue_total !! Define path name for outputs



    call cpu_time(temps_debut)
    continue
        62 format(10X,'X coord',10X,'Y coord',10X,'Generation',10X,'Radius',10X,'Time',10X,'Local Time')
        63 format(12X,'X coord',12X,'Y coord',12X,'Time',12X,'Local Time',12X,'Threshold')
        64 format('#',3X,'Segment Index',3X,'Hitting Probability')
        71 format(2(3X,ES20.10))
        72 format(2(3X,f16.8),2X,I6,3(2X,f16.8))
        73 format(5(4X,f16.8))
        74 format(2X,I6,10X,ES10.4)
        !76 format(8(2X,I6),f14.5) !! format for file(36)=debug.txt
        91 format('# Polygon',I4,', Gmax',I4,', Length',f8.2,', Thickness ',ES8.2)
        92 format('# No. Particle',I9,', Rho ratio',f6.2,', Reactivity ',ES8.2,', Diff coef ',ES8.2)
        94 format('PSD @ YYE. Version 2.0.2.. Simulation on ',I4,2('-',I2),'   ',I2,'h',I2)
        95 format("Done for particle",I9,", Time used (s)",ES11.3)
        99 format('It spent',f10.4,' seconds for the whole program.')
    continue



    call obtaindata
    call date_and_time(values=date)
    rue_total = "rfl-"
    open(unit=30,file=TRIM(rue_total)//"record.txt")
        write(30,94) date(1), date(2), date(3), date(5), date(6)
    open(unit=31,file=TRIM(rue_total)//"coord.txt")
        write(31,94) date(1), date(2), date(3), date(5), date(6)
    open(unit=32,file=TRIM(rue_total)//"particle_traj.txt")
        write(32,94) date(1), date(2), date(3), date(5), date(6)
    write(32,62); write(32,*)
    open(unit=33,file=TRIM(rue_total)//"particle_final.txt")
        write(33,94) date(1), date(2), date(3), date(5), date(6)
        write(33,63); write(33,*)
    open(unit=34,file=TRIM(rue_total)//"pbb.txt")
        write(34,94) date(1), date(2), date(3), date(5), date(6)
        write(34,91) polygon_shape, Generation_max, Length, thickness
        write(34,92) num_particle, rho_ratio, q_rect, diff_coef; write(34,*)
        write(34,64); write(34,*)
    !open(unit=35,file=TRIM(rue_total)//"zeta.txt"); write(35,94) date(1), date(2), date(3), date(5), date(6)
    open(unit=36,file=TRIM(rue_total)//"debug.txt")
        write(36,94) date(1), date(2), date(3), date(5), date(6); !write(36,*)
    open(unit=37,file=TRIM(rue_total)//"temps.txt")
        write(37,94) date(1), date(2), date(3), date(5), date(6); 
        write(37,91) polygon_shape, Generation_max, Length, thickness
        write(37,92) num_particle, rho_ratio, q_rect, diff_coef; write(37,*)
    num_points = polygon_shape * 4**Generation_max + 1 !! + 1 to form a circle
    allocate(position(num_particle,2)); allocate(coord_x(num_points)); allocate(coord_y(num_points))
    allocate(account(num_points-1)) !! account = count @ Python
    coord_x = 0.0d0; coord_y = 0.0d0; account = 0.0d0
    traj_max = 5
    avg_seuil = 0.0d0; avg_local = 0.0d0; avg_temps = 0.0d0


	
	!!! Generate initial polygon coordinates
    angle_ini = twopi / polygon_shape; tta = pi * (1.0d0 - 2.0d0/polygon_shape) 
	!!!!! Note, all integers in Python should be converted as real values for +-*/, else mod taken...
    new_x = - Length/2.0d0; new_y = new_x * tan(tta/2.0d0)
    do i = 0, polygon_shape - 1 
        j = i * 4**Generation_max + 1
        coord_x(j) = new_x; coord_y(j) = new_y !coordinate.append((new_x,new_y))
        new_x = new_x + Length * cos(i * angle_ini)
        new_y = new_y + Length * sin(i * angle_ini)
    end do
    coord_x(num_points) = coord_x(1); coord_y(num_points) = coord_y(1) ! To form a circle.
	!!! Insert fractal points inside
    do i = 0, Generation_max - 1
        k = 4**i * polygon_shape
        out_gap = 4**(Generation_max - i); int_gap = out_gap / 4
        do j = 0, k-1
            l = out_gap * j
            outpoint(1) = 1 + l
            outpoint(2) = outpoint(1) + out_gap
            intpoint(1) = outpoint(1) + int_gap
            intpoint(2) = intpoint(1) + int_gap
            intpoint(3) = intpoint(2) + int_gap
            call koch_line(coord_x(outpoint(1)), coord_y(outpoint(1)), coord_x(outpoint(2)), coord_y(outpoint(2)), &
			&pi/alpha_angle_ratio, direct, coord_x(intpoint(1)), coord_y(intpoint(1)), &
			&coord_x(intpoint(2)), coord_y(intpoint(2)), coord_x(intpoint(3)), coord_y(intpoint(3)))
        end do
    end do
    !u = sqrt( 2 * (1 - cos(pi/alpha_angle_ratio)) ) + 2
    !seg_l = Length / u**Generation_max !! 22/05/23 added. Give the minimal segment length.

    !goto 1102
	!!! Output fractal points
    do i=1,num_points
        write(31,71) coord_x(i), coord_y(i)
    end do
    !1102 continue
	
	
    allocate(rec_index(Generation_max + 2, 2))
    out_gap = 4**Generation_max * polygon_shape !! 28/04/23 Added. (g=0) -> 3, (g=1) -> 12, (g=2) -> 48
    !!! Diffusion for each particle
    do i = 1, num_particle
        1100 continue
        
		!! Initialization
        ell = 0.0d0; call random_number(u)
        ell_hat = -log(u) / q_rect !! PDF needed. psi(l) = q*exp(-q*l); 
        time_diff = 0.0d0

        position(i,:) = 0.0d0; g_act = 1
        if (i <= traj_max) then
            write(32,*) "# ",i !! Furnish particle index in traj.
            write(32,72) position(i,1), position(i,2), g_act, 0.0d0, time_diff, ell
        end if
        !write(36,*) "Position Initialized!"; !write(36,*) "***** i =", i, "*****"
        !write(36,*) "Local Time Threshold", ell_hat; !write(36,*)
        
        rec_index(1,1) = 1; rec_index(1,2) = num_points
        do j = 2, Generation_max + 2
            rec_index(j,1) = rec_index(1,1) ! = 1
            rec_index(j,2) = rec_index(1,2) ! = num_points
        end do


		!### [1.1] **Initialize `g_act` ** (for start-up)
        g_act = 2 !! g_act = {1, 2, 3, ..., Gmax+2} ~ g = {0, 0, 1, ..., Gmax}
        radius_min = Length / 4.0d0 !! Arbitrary diffusion at first step but still inside the boundary.
        call diffusion(position(i,1), position(i,2), radius_min)
        time_diff = time_diff + radius_min**2 / (4.0d0 * diff_coef)
        if (i <= traj_max) write(32,72) position(i,1), position(i,2), g_act, radius_min, time_diff, ell

		
		!! Iteration for GAFRW
        do while (g_act < Generation_max + 3) !! 等价于 (1.gt.0)
            
            int_gap = 4**(Generation_max - g_act + 2)
            !write(36,*) "[1.2] Determine g & interval."
            !write(36,*) "g_act / int_gap =", g_act, int_gap
            
            
            continue
				!### [2] **Find 2 nearest points** (not for start-up)
                !! 取上一级端点往左右额外延伸各一点
                if (g_act == 2) then
                    k1 = rec_index(g_act-1,1) ! = 1 
                    k2 = rec_index(g_act-1,2) - int_gap ! = num_points - int_gap
                    infvl = k1 + int_gap; supvl = k2 - int_gap
                    search_length = polygon_shape - 2
                else
                    infvl = rec_index(g_act-1,1)
                    supvl = rec_index(g_act-1,2)
                    call cherche(infvl, supvl, int_gap, k1, k2, search_length) !! subroutine cherche(niveau, connu1, connu2, gap, nouveau1, nouveau2, combien)
                end if
                l1 = span(coord_x(k1), coord_y(k1), position(i,1), position(i,2))
                l2 = span(coord_x(k2), coord_y(k2), position(i,1), position(i,2))
                

                allocate(search_index(search_length))
                do j = 1, search_length
                    search_index(j) = noloop(infvl + (j-1) * int_gap)
                end do
                
                !write(36,*) "[2] Determine search space"
                do j = 1, search_length !! 为何放弃 do j = infvl, supvl, int_gap? 答曰考虑区间跨过首尾分界点不严格递增
                    k = search_index(j)
                    u = span(coord_x(k), coord_y(k), position(i,1), position(i,2))
                    !write(36,*) "search index / distance =", k, u
                    call ordre(k1, k2, l1, l2, k, u) !! subroutine ordre(k1,k2,l1,l2, k,l)
                end do
                rec_index(g_act,1)= k1; rec_index(g_act,2)= k2
                deallocate(search_index)
            continue

			
			
			!### [3] **Determine diffusion radius** (able for start-up)
            search_length = noloop(rec_index(g_act,2) - rec_index(g_act,1)) / int_gap + 1
            !write(36,*) "[3] Determine diffusion radius"
            !write(36,*) "endpoint index", rec_index(g_act,1), rec_index(g_act,2), " search length =", search_length
            allocate(search_index(search_length))
            do j = 1, search_length
                search_index(j) = noloop(rec_index(g_act,1) + (j-1) * int_gap)
            end do
            !write(36,*) "Search index", search_index(:)
            radius_min = Length!/4**(g_act-1)
            k = search_index(1) !! 04/05/2023 添加，防自检报错
            do jj = 1, search_length - 1
                j = search_index(jj)
                l = search_index(jj+1) !! 04/05/2023 修改search_length大小，减少循环中加减运算。
                call distance(coord_x(j),coord_y(j), coord_x(l),coord_y(l), position(i,1), position(i,2), case, radius)
                if (case.ne.-1) then
                    if (radius_min > radius) then
                        radius_min = radius
                        k = j
                    end if
                else
                    deallocate(search_index) !! 25/05/23 added. Avoid < Attempting to allocate already allocated variable 'search_index' >
                    goto 1100
                end if
            end do
            deallocate(search_index)



			!### [4] **Diffusion** with given `radius_min`.
            !write(36,*) "[4] Diffusion uniform"
            !write(36,*) "R_min =", radius_min!, "(last) radius =", radius
            !write(36,*) "g before", g_act
            if (radius_min == Length) goto 1100
            if (case.ne.+2) then
            !if (radius_min > 0.0d0) then
                call diffusion(position(i,1), position(i,2), radius_min)
                time_diff = time_diff + radius_min**2 / (4.0d0 * diff_coef)
                if (i <= traj_max) write(32,72) position(i,1), position(i,2), g_act, radius_min, time_diff, ell !! Record trajs of each particle
                !write(36,*) "Position afterwards.", position(i,1), position(i,2)
            else !! radius_min == 0.0d0
                if (g_act == Generation_max+2) then
                    cross = 0; inside = 1; touch = 0

                    do while (inside == 1) !! Repete the process once the particle enter the layer.
                        !call reflection(xy1234,xy0,cross,inside,touche)
                        call reflection(coord_x(k),coord_y(k), coord_x(noloop(k+1)),coord_y(noloop(k+1)), &
                        &position(i,1),position(i,2), coord_x(noloop(k-1)),coord_y(noloop(k-1)), & 
                        &coord_x(noloop(k+2)),coord_y(noloop(k+2)), cross, inside, touch, radius_min, time_diff, ell, ell_hat)
                        if (abs(cross) == 1) k = noloop(k + cross) !! Perhaps hitting segments change.
                        !write(36,*) "Call Reflection + 1", i
                        !write(36,*) "Position afterwards.", position(i,1), position(i,2)
                        !write(36,*) "cross, inside, touch = ", cross, inside, touch
                        !write(36,*) "Time =", time_diff, "Local Time =", ell
                        if (i <= traj_max) write(32,72) position(i,1), position(i,2), g_act, radius_min, time_diff, ell
                    end do
                
                    if ((ell >= ell_hat).or.(inside == 2)) then
                        if (i <= traj_max) then
                            write(32,*) "Final Time Used =", time_diff
                            write(32,*) "Final Local Time =", ell
                            write(32,*) "Local Time Threshold =", ell_hat
                            write(32,*)
                        end if
                        if (i <= 20) write(33,73) position(i,1), position(i,2), time_diff, ell, ell_hat !! Print final positions, diffusion time, local time
                        account(k) = account(k) + 1
                        !write(36,*) "***** Particle Attached! *****"; !write(36,*) "i = ",i; !write(36,*); !write(36,*); !write(36,*)
                        !write(37,'(ES14.6)') time_diff
                        avg_temps = avg_temps + time_diff
                        avg_local = avg_local + ell
                        avg_seuil = avg_seuil + ell_hat
                        goto 1101 !! Jump the out while loop.
                    end if
                else
                    g_act = g_act + 1
                end if
            end if
            !write(36,*) "g after ", g_act



            !### [5] **Determine if go out** (not for start up)
            !write(36,*) "[5] Verify correct region."
            !do g_test = 2, g_act !! Search all passed level g. Make sure the particle not out of region at all level.
            g_test = 2
            do while (g_test <= g_act)
                int_gap = 4**(Generation_max - g_test + 2) ! 定义该级别搜索步长

                if (g_test.eq.2) then
                    k1 = rec_index(g_test-1,1) ! = 1 
                    k2 = rec_index(g_test-1,2) - int_gap ! = num_points - int_gap
                    infvl = k1 + int_gap; supvl = k2 - int_gap
                    search_length = polygon_shape - 2
                else
                    infvl = rec_index(g_test-1,1) !! Define search inf
                    supvl = rec_index(g_test-1,2) !! Define search sup
                    call cherche(infvl, supvl, int_gap, k1, k2, search_length) !! subroutine cherche(niveau, connu1, connu2, gap, nouveau1, nouveau2, combien)
                end if
                l1 = span(coord_x(k1), coord_y(k1), position(i,1), position(i,2))
                l2 = span(coord_x(k2), coord_y(k2), position(i,1), position(i,2))

                
                
                allocate(search_index(search_length))
                do j = 1, search_length
                    search_index(j) = noloop(infvl + (j-1) * int_gap)
                end do
                rec_old1 = rec_index(g_test,1); rec_old2 = rec_index(g_test,2) ! 该级别下原来最近两端点位置

                !write(36,*) "g_test / int_gap =", g_test, int_gap
                !write(36,*) "Endpoint index", rec_index(g_test-1,:)
                !write(36,*) "k1/k2 & l1/l2", k1, k2, l1, l2
                !write(36,*) "search length = ", search_length
                do j = 1, search_length
                    k = search_index(j)
                    u = span(coord_x(k), coord_y(k), position(i,1), position(i,2))
                    !write(36,*) "search_index / distance", k, u
                    call ordre(k1,k2, l1,l2, k,u) !! subroutine ordre(k1,k2,l1,l2, k,l)
                end do
                rec_index(g_test,1) = k1; rec_index(g_test,2) = k2
                deallocate(search_index)
                !write(36,*) "k1/k2 & l1/l2", k1, k2, l1, l2
                
                
                
                if ((noloop(rec_old1).ne.noloop(rec_index(g_test,1))).or.(noloop(rec_old2).ne.noloop(rec_index(g_test,2)))) then
                    !write(36,*) "Break @ g_test =", g_test
                    g_act = g_test !! 之后g要重新算，固定为首次出现异常的g_test
                    !if (g_act > Generation_max+2) g_act = Generation_max
                end if

                g_test = g_test + 1
            end do
            

            if (noloop(rec_index(g_act,2) - rec_index(g_act,1)) == int_gap) then
                if (g_act < Generation_max + 2) g_act = g_act + 1
            end if

            
            !write(36,*) "Done once for while loop."
            !write(36,*) "i =",i,"g_act =",g_act
            !write(36,*) "#rec(1)", rec_index(:,1)
            !write(36,*) "#rec(2)", rec_index(:,2)
            !write(36,*); !write(36,*)

	
        end do 

        1101 continue
        if (MOD(i,ceiling(num_particle/5.0d0)).eq.0) then
            call cpu_time(temps_provisoire)
            write(*,95) i,temps_provisoire-temps_debut
            write(30,95) i,temps_provisoire-temps_debut
        end if
    end do
	


	!!! Normalize account
    do i = 1, num_points - 1
        account(i) = account(i) / num_particle
        write(34,74) i, account(i)
    end do
    
    write(30,*); write(30,*)
    avg_seuil = avg_seuil / num_particle; write(30,*) "Average Threshold  ", avg_seuil
    avg_local = avg_local / num_particle; write(30,*) "Average Local Time ", avg_local
    avg_temps = avg_temps / num_particle; write(30,*) "Average Diff Time  ", avg_temps
    write(30,*); write(30,*); write(30,*)
    


    goto 1104
        qmax = 10
        allocate(zeta(qmax))
        do i = 1, qmax
            u = 0.0d0
            do k = 1, num_points - 1
                u = u + account(k)**i
            end do
            zeta(i) = u
            write(35,*) i,zeta(i)
        end do
        deallocate(zeta)
    1104 continue



    deallocate(position); deallocate(coord_x); deallocate(coord_y); deallocate(account)
    deallocate(rec_index)
    close(31); close(32); close(33); close(34); !close(35); 
    close(36); close(37)


    call cpu_time(temps_fin)
    write(*,99) temps_fin-temps_debut !! Write the total time consumed to the screen.
    write(30,99) temps_fin-temps_debut
    close(30)
	! close files?
end Program main





MODULE MATHS
    implicit none
    real(kind=8),parameter :: pi = 4.0d0*atan(1.0d0), twopi = 2.0d0*pi
    real(kind=8),parameter :: sqrt2 = sqrt(2.0d0), sqrt3 = sqrt(3.0d0), septds3 = 7.0d0/sqrt3, cinqds3 = 5.0d0/sqrt3
    CONTAINS
    function normaldist(mean,std,n) result(r)
        implicit none
        real(kind=8),intent(in) :: mean,std
        integer,intent(in) :: n
        real(kind=8) :: r(n,2)
        real(kind=8),dimension(n,2) :: zeta
        !call random_seed()
        call random_number(zeta)
        r(:,1) = dsqrt(-2.0d0*log(zeta(:,1)))*cos(twopi*zeta(:,2))
        r(:,2) = dsqrt(-2.0d0*log(zeta(:,1)))*sin(twopi*zeta(:,2))
        r = mean + std * r
    end function normaldist
END MODULE MATHS





subroutine obtaindata
    implicit none
	include "para.h"
    integer :: status1
    character*128 :: msg,ruein
	
    ruein="input_rfl.txt"
    open(unit=201,file=TRIM(ruein),form='formatted',status='old',action='read',iostat=status1,iomsg=msg)
	
    read(201,*) polygon_shape !! 3 = triangular; 4 = square; 5 = star; etc.
    read(201,*) Generation_max !! The maximum value of generation / recursion
    read(201,*) num_particle !! Number of particles
    read(201,*) Length !! Edge length for initial polygon
    read(201,*) alpha_angle_ratio !! Fractal angle 
    read(201,*) direct
    read(201,*) thickness
    read(201,*) rho_ratio
    read(201,*) q_rect
    read(201,*) diff_coef
	
    close(201)
end subroutine
subroutine koch_line(x1,y1,x2,y2,alpha,direction, b1,b2,c1,c2,d1,d2)
    use MATHS
    implicit none
    real*8, intent(in) :: x1,y1,x2,y2,alpha ! coordinates for start/end points & fractal angle
    integer, intent(in) :: direction ! 凹 +1 = concave, 凸 -1 = convex 
    real*8, intent(out) :: b1,b2,c1,c2,d1,d2 ! coordinates for inserted points
    real*8 :: deltax,deltay,l,coef,segm,beta,theta,changex,changey,degree
    
    deltax = x2 - x1; deltay = y2 - y1 
    l = sqrt((deltax)**2 + (deltay)**2) ! the length of the line
    coef = sqrt( 2 * (1 - cos(alpha)) ) + 2; segm = l / (coef); beta = (pi - alpha)/2
    if (x1.eq.x2) then
        if (y1.lt.y2) theta = + pi / 2
        if (y1.gt.y2) theta = - pi / 2
    else
        theta = atan(deltay/deltax)
    end if
    if (x1.gt.x2) theta = theta + pi
    
    changex = deltax / coef; changey = deltay / coef
    b1 = x1 + changex; b2 = y1 + changey ! second point: one third in each direction from the first point
    degree = theta + beta * direction
    c1 = b1 + segm * cos(degree); c2 = b2 + segm * sin(degree) ! third point: rotation for multiple of 60 degrees
    d1 = x2 - changex; d2 = y2 - changey ! fourth point: two thirds in each direction from the first point
end subroutine
subroutine distance(x1, y1, x2, y2, x0, y0, dehors, dmin)
    use MATHS
    implicit none
	include "para.h"
    real*8, intent(in) :: x1, y1, x2, y2, x0, y0 !! Given 2 points: (x1,y1) & (x2,y2); Particle: (x0,y0)
    integer,intent(out):: dehors
    real*8, intent(out):: dmin
    real*8 :: deltax, deltay, dist_xy, eps, theta, st, ct, xdist, ydist, dxnew, dynew, d1, d2
    real*8,external :: span

    deltax = x2 - x1; deltay = y2 - y1
    !dist_xy = sqrt(deltax**2 + deltay**2)
    dist_xy = span(x1,y1, x2,y2)
    !eps = dist_xy * thickness !eps = dist_xy / 10**3, 19/05/23 modified
    eps = thickness !! 20/05/23 modified. Independent variable, not a ratio.
    
    if (x1.eq.x2) then
        if (y1.le.y2) then
            theta = + pi / 2.0d0
        else
            theta = - pi / 2.0d0
        end if
    else
        theta = atan(deltay/deltax)
    end if
    if (x1.gt.x2) theta = theta + pi
    st = sin(theta); ct = cos(theta)
    
    xdist = x0 - x1; ydist = y0 - y1
    dxnew = + ct * xdist + st * ydist; dynew = - st * xdist + ct * ydist
    
    if (dynew.lt.0) then !# To avoid errors while g>2. 或出现某点距离某边负距离，未抵达边界便停止扩散。
        dmin = Length; dehors = -1
    else
        if ((dxnew.ge.0).and.(dxnew.le.dist_xy)) then
            dmin = dynew; dehors = +1
            if (dmin.lt.eps) then
                !dmin = 0.0d0
                dehors = +2
            end if
        else
            dehors = 0
            d1 = sqrt(xdist**2 + ydist**2)
            d2 = sqrt((x0-x2)**2 + (y0-y2)**2)
            if (d1.le.d2) then
                dmin = d1
            else
                dmin = d2
            end if
        end if
    end if
end subroutine
real*8 function span(x1, y1, x2, y2)
    implicit none
    real*8 :: x1, y1, x2, y2
    span = sqrt((x1-x2)**2 + (y1-y2)**2)
    return
end function span




subroutine cherche(connu1, connu2, gap, nouveau1, nouveau2, combien)
    !! Here is a procedure to extend the search space
    implicit none !! "niveau" = g_act/g_test; "gap" = distance of two point indices
    integer,intent(in) :: connu1, connu2, gap !! "connu1/2" refer to obvious search range, namely #rec(g_act-1,:)
    integer,intent(out):: nouveau1, nouveau2, combien !! "nouveau1/2" = two new point indices after extension.
    integer,external :: noloop

    combien = noloop(connu2 - connu1) / gap + 1
    nouveau1 = noloop(connu1 - gap)
    nouveau2 = noloop(connu2 + gap)
end subroutine
integer function noloop(index)
    implicit none
    include "para.h"
    integer :: index

    if (index >= num_points) then
        noloop = mod(index, num_points-1)
    else
        if (index <= 0) index = index + num_points - 1
        noloop = index
    end if
    return
end function noloop



subroutine ordre(k1,k2,l1,l2, k,l)
    !! Here is a procedure to find 2 nearest points. Added 03/05/2023
    implicit none
    include "para.h"
    integer,intent(in) :: k
    integer :: k1, k2, kint
    real*8,intent(in) :: l
    real*8 :: l1, l2, u, lint

    !! 首先判断已知两点是否等距，若是则随机放弃一个
    if (l1 == l2) then
        if (l < l1) then
            call random_number(u)
            if (u < 0.5d0) then
                k1 = k; l1 = l
            else
                k2 = k; l2 = l
            end if
        end if
    else
        if (l1 > l2) then !! 调整顺序确保靠前的更近
            kint = k1; k1 = k2; k2 = kint
            lint = l1; l1 = l2; l2 = lint
        end if

        if (l < l1) then
            k2 = k1; k1 = k
            l2 = l1; l1 = l
        else
            if (l < l2) then
                k2 = k; l2 = l
                goto 3101
            end if

            if (l == l2) then !! 若二者相等，则随机放弃其一
                call random_number(u)
                if (u < 0.5d0) then
                    k2 = k; l2 = l
                end if
                goto 3101
            end if

            3101 continue
        end if
    end if

    !! 确保小序号在前。同时考虑绕过首尾分界点情况
    if (((k1 > k2).and.((k1 - k2) < num_points/2)).or.&
    & ((k1 < k2).and.((k2 - k1) > num_points/2))) then
        lint = l1; l1 = l2; l2 = lint
        kint = k1; k1 = k2; k2 = kint
    end if
end subroutine



subroutine diffusion(x0, y0, rayon)
    !! Here is a procedure to realize uniform diffusion on the circle of radius = "rayon". Added 03/05/23.
    use MATHS
    implicit none
    real*8,intent(in) :: rayon
    real*8 :: x0, y0
    real*8 :: u, random_angle, move_x, move_y

    call random_number(u); random_angle = u * twopi
    move_x = rayon * cos(random_angle); move_y = rayon * sin(random_angle)
    x0 = x0 + move_x; y0 = y0 + move_y
end subroutine



subroutine reflection(x1,y1, x2,y2, x0,y0, x3,y3, x4,y4, croise, interieure, touche, rayon, time, local_time, threshold)
    !! Here is a procedure to realize reflection when particle inside boundary layer. Added 19/05/23.
    use MATHS
    implicit none
    include "para.h"
    real*8, intent(in) :: x1,y1, x2,y2, x3,y3, x4,y4 !! Boundary points with order 3-1-2-4
    real*8,intent(out) :: rayon
    real*8 :: x0, y0, time, local_time !! Particle coordinates, inside 1-2. time t; local_time ell_t
    real*8, intent(in) :: threshold
    integer :: croise, interieure, touche !! encounter: how many times reflection on the boundary; interieure: inside (=1) layer or not (=0); touche: once touch boundary (=1) or not (=0);
    integer :: case1, case2, case_new !! Left/Right boundary point environments. Condition after diffusion.
    real*8 :: d1, d2, radius_min, d14, d23 !! Left/Right boundary point distance. 
    real*8 :: deltax, deltay, dist_xy, theta, st, ct, xdist, ydist, dxnew, dynew, u !! Variables for Rotations.
    real*8 :: eps, rho, delta !! epaisseur, rayon, temps.
    real*8,external :: span
    !!! 1. Rotation; 2. Determine case & radius; 3. Diffusion; 4. Reflection; 5. Determine inside or not; 6. Anti-rotation


    !!! [1] Rotation
    deltax = x2 - x1; deltay = y2 - y1
    dist_xy = sqrt(deltax**2 + deltay**2) !! Do not understand. But this expression run faster.
    !dist_xy = span(x1,y1,x2,y2)
    !dist_xy = seg_l
    
    if (x1.eq.x2) then
        if (y1.le.y2) then
            theta = + pi / 2.0d0
        else
            theta = - pi / 2.0d0
        end if
    else
        theta = atan(deltay/deltax)
    end if
    if (x1.gt.x2) theta = theta + pi
    st = sin(theta); ct = cos(theta)
    xdist = x0 - x1; ydist = y0 - y1
    dxnew = + ct * xdist + st * ydist; dynew = - st * xdist + ct * ydist
    
    
    !!! [2] Determine case & radius
    eps = thickness
    rho = rho_ratio * eps
    case1 = 6; case2 = 6 !! Help to determine whether we call "distance", and then which point is closer.

    if (abs(croise) == 2) then
        if (croise == 2) then
            d2 = span(x2,y2, x0,y0)
            case2 = 0
            radius_min = min(rho, d2)
        else
            d1 = span(x1,y1, x0,y0)
            case1 = 0
            radius_min = min(rho, d1)
        end if
        goto 3202
    end if

    if (dxnew < septds3 * eps) then
        d23 = span(x2,y2,x3,y3)
        !d23 = sqrt((x2-x3)**2 + (y2-y3)**2)
        if (d23 > sqrt2 * dist_xy) then
            d1 = span(x1,y1,x0,y0)
        else
            !call distance(x3,y3, x1,y1, x0,y0, case1, d1)
            d1 = abs(sqrt3 * dxnew - dynew) / 2.0d0
            case1 = 0 
        end if
        radius_min = min(rho, d1); croise = -1
    else if (dxnew > dist_xy - septds3 * eps) then
        d14 = span(x1,y1,x4,y4)
        !d14 = sqrt((x1-x4)**2 + (y1-y4)**2)
        if (d14 > sqrt2 * dist_xy) then
            d2 = span(x2,y2,x0,y0)
        else
            !call distance(x2,y2, x4,y4, x0,y0, case2, d2)
            d2 = abs(sqrt3 * dxnew + dynew - dist_xy) / 2.0d0
            case2 = 0
        end if
        radius_min = min(rho, d2); croise = +1
    else
        radius_min = rho
    end if
    3202 continue
    rayon = radius_min


    if (radius_min < 1.0d-4 * eps) then !! Once the particle is trapped to a point, stop the simulation directly.
        interieure = 2
        call random_number(u)
        if ((case1==0).and.(u<0.5d0)) croise = -1
        if ((case2==0).and.(u<0.5d0)) croise = +1
        goto 3201
    end if
    
    
    !!! [3] Diffusion
    call diffusion(dxnew, dynew, radius_min)
    !!! [4] Reflection ? 
    if (dynew < 0.0d0) then 
        touche = touche + 1
        dynew = - dynew
    end if


    !!! [6] Anti-Rotation
    xdist = ct * dxnew - st * dynew; ydist = st * dxnew + ct * dynew
    x0 = x1 + xdist; y0 = y1 + ydist
    

    !!! [7] Add Local Time
    if (touche > 0) then !! Only the reflection case can increase local time.
        delta = radius_min**2 / (4.0d0 * diff_coef)
        time = time + delta
        local_time = local_time + sqrt(pi * diff_coef * delta / 2.0d0)
        if (local_time >= threshold) then
            interieure = 2 !! 25/05/23 added. Stop "Reflection" even inside the diffusion layer.
            goto 3201
        end if    
    end if    
    

    !!! [5] Inside ?
    if (croise == 0) then !! Assure that even after the diffusion, the particle doesn't enter other boundary. septds3!
        if (dynew > eps) interieure = 0
    else
        !if (case1 > case2) then !! Particle closer to right boundary point
        if (croise > 0) then
            call distance(x2,y2, x4,y4, x0,y0, case_new, d2)
            select case (case_new)
                case(0) !! angle = 120
                    if (dynew > eps) then
                        interieure = 0
                        goto 3201
                    end if
                    if (dxnew > dist_xy) then
                        if (d2 > eps) then
                            interieure = 0
                        else
                            croise = +2
                        end if
                    end if
                case(1) !! angle = 60, but not inside layer
                    if (dynew > eps) interieure = 0
                case(2) !! angle = 60, but inside layer
                    if (d2 < dynew) then
                        croise = +1 !! Attach to right segment
                    else
                        croise = 0
                    end if
            end select
        else !! Particle closer to left boundary point
            call distance(x3,y3, x1,y1, x0,y0, case_new, d1)
            select case (case_new)
                case(0) !! angle = 120
                    if (dynew > eps) then
                        interieure = 0
                        goto 3201
                    end if
                    if (dxnew < 0) then
                        if (d1 > eps) then
                            interieure = 0
                        else
                            croise = -2
                        end if
                    end if
                case(1) !! angle = 60, but not inside layer
                    if (dynew > eps) interieure = 0
                case(2) !! angle = 60, but inside layer
                    if (d1 < dynew) then
                        croise = -1 !! Attach to left segment
                    else
                        croise = 0
                    end if
            end select
        end if
    end if
    3201 continue
    
end subroutine
    
