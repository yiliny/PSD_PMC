!!!!! Diffusion near arbitrary Koch
!!!!! Yilin YE @ 04/2023


Program main
    use MATHS
    implicit none
    include "para.h"
    integer :: i,j,k,l,jj ! Loop variables
    integer :: outpoint(2), intpoint(3), out_gap, int_gap, rec_old1, rec_old2
    integer :: qmax, g_act, g_test, case, infvl, supvl, search_length, k1, k2
    real*8 :: angle_ini, tta, new_x, new_y, radius_min, radius, u, l1, l2
    real*8,allocatable :: position(:,:), coord_x(:), coord_y(:), account(:), zeta(:)
    integer,allocatable :: rec_index(:,:), search_index(:)
    integer,external:: noloop
    real*8,external :: span
	
    integer :: date(8)
    real*8 :: temps_debut, temps_fin!, time_begin, time_end !! Define variables to record time consumed
    character*128 :: rue_total !! Define path name for outputs



    call cpu_time(temps_debut)
    continue
        !76 format(8(2X,I6),f14.5) !! format for file(36)=debug.txt
        94 format('PSD @ YYE. Simulation done on ',I4,2('-',I2),'   ',I2,'h',I2)
        99 format('It spent',f10.4,' seconds for the whole program.')
    continue



    call obtaindata
    call date_and_time(values=date)
    rue_total = "diff-"
    open(unit=31,file=TRIM(rue_total)//"coord.txt"); write(31,94) date(1), date(2), date(3), date(5), date(6)
    open(unit=32,file=TRIM(rue_total)//"particle_traj.txt"); write(32,94) date(1), date(2), date(3), date(5), date(6)
    open(unit=33,file=TRIM(rue_total)//"particle_final.txt"); write(33,94) date(1), date(2), date(3), date(5), date(6)
    open(unit=34,file=TRIM(rue_total)//"sgm_pbb.txt"); write(34,94) date(1), date(2), date(3), date(5), date(6)
    open(unit=35,file=TRIM(rue_total)//"zeta.txt"); write(35,94) date(1), date(2), date(3), date(5), date(6)
    open(unit=36,file=TRIM(rue_total)//"debug_v3.txt"); write(36,94) date(1), date(2), date(3), date(5), date(6)
    num_points = polygon_shape * 4**Generation_max + 1 !! + 1 to form a circle
    allocate(position(num_particle,2)); allocate(coord_x(num_points)); allocate(coord_y(num_points))
    allocate(account(num_points-1)) !! account = count @ Python
    coord_x = 0.0d0; coord_y = 0.0d0; account = 0.0d0


	
	!!! Generate initial polygon coordinates
    angle_ini = twopi / polygon_shape; tta = pi * (1.0d0 - 2.0d0/polygon_shape) 
	!!!!! Note, all integers in Python should be converted as real values for +-*/, else mod taken...
    new_x = - Length/2.0d0; new_y = new_x * tan(tta/2)
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

    goto 1102
	!!! Output fractal points
    do i=1,num_points
        write(31,*) coord_x(i), coord_y(i)
    end do
    1102 continue
	
	
    allocate(rec_index(Generation_max + 2, 2))
    out_gap = 4**Generation_max * polygon_shape !! 28/04/23 Added. (g=0) -> 3, (g=1) -> 12, (g=2) -> 48
    !!! Diffusion for each particle
    do i = 1, num_particle
        1100 continue
        
		!! Initialization
        if (i.le.50) write(32,*) "# ",i !! Furnish particle index in traj.
        position(i,:) = 0.0d0; g_act = 1
        if (i.le.50) write(32,*) position(i,1), position(i,2), g_act
        !write(36,*) "Position Initialized!"; !write(36,*) "***** i = *****", i; !write(36,*)
        
        rec_index(1,1) = 1; rec_index(1,2) = num_points
        do j = 2, Generation_max + 2
            rec_index(j,1) = rec_index(1,1) ! = 1
            rec_index(j,2) = rec_index(1,2) ! = num_points
        end do

		!### [1.1] **Initialize `g_act` ** (for start-up)
        g_act = 2 !! g_act = {1, 2, 3, ..., Gmax+2} ~ g = {0, 0, 1, ..., Gmax}
        radius_min = Length / 4.0d0 !! Arbitrary diffusion at first step but still inside the boundary.
        call diffusion(position(i,1), position(i,2), radius_min)
        if (i.le.50) write(32,*) position(i,1), position(i,2), g_act

		
		!! Iteration for GAFRW
        do while (g_act < Generation_max + 3) !! 等价于 (1.gt.0)
            
            int_gap = 4**(Generation_max - g_act + 2)
            !write(36,*) "[1.2] g_act =", g_act, "int_gap =", int_gap
            
            
            continue
				!### [2] **Find 2 nearest points** (not for start-up)
                !! 取上一级端点往左右额外延伸各一点
                if (g_act.eq.2) then
                    k1 = rec_index(g_act-1,1) ! = 1 
                    k2 = rec_index(g_act-1,2) - int_gap ! = num_points - int_gap
                    infvl = k1 + int_gap; supvl = k2 - int_gap
                    search_length = polygon_shape - 2
                else
                    infvl = rec_index(g_act-1,1)
                    supvl = rec_index(g_act-1,2)
                    call cherche(g_act, infvl, supvl, int_gap, k1, k2, search_length) !! subroutine cherche(niveau, connu1, connu2, gap, nouveau1, nouveau2, combien)
                end if
                l1 = span(coord_x(k1), coord_y(k1), position(i,1), position(i,2))
                l2 = span(coord_x(k2), coord_y(k2), position(i,1), position(i,2))
                

                allocate(search_index(search_length))
                do j = 1, search_length
                    search_index(j) = noloop(infvl + (j-1) * int_gap)
                end do
                
                do j = 1, search_length !! 为何放弃 do j = infvl, supvl, int_gap? 答曰考虑区间跨过首尾分界点不严格递增
                    k = search_index(j)
                    u = span(coord_x(k), coord_y(k), position(i,1), position(i,2))
                    !write(36,*) "[2] g_act =", g_act, "search_index & span =",k,u
                    call ordre(k1, k2, l1, l2, k, u) !! subroutine ordre(k1,k2,l1,l2, k,l)
                end do
                rec_index(g_act,1)= k1; rec_index(g_act,2)= k2
                deallocate(search_index)
            continue

			
			
			!### [3] **Determine diffusion radius** (able for start-up)
            search_length = noloop(rec_index(g_act,2) - rec_index(g_act,1) + num_points - 1) / int_gap + 1
            !write(36,*) "[3]: #rec(g_act,:) =", rec_index(g_act,1), rec_index(g_act,2), "search length =", search_length
            allocate(search_index(search_length))
            do j = 1, search_length
                search_index(j) = noloop(rec_index(g_act,1) + (j-1) * int_gap)
            end do
            !write(36,*) search_index(:)
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
                end if
            end do
            deallocate(search_index)



			!### [4] **Diffusion** with given `radius_min`.
            !write(36,*) "[4] starts w/ rmin =", radius_min, "(last) radius =", radius
            if (radius_min == Length) goto 1100
            if (radius_min > 0.0d0) then
                call diffusion(position(i,1), position(i,2), radius_min)
                if (i.le.50) write(32,*) position(i,1), position(i,2), g_act !! Record trajs of each particle
                !write(36,*) "Diffusion->", position(i,1), position(i,2)
            else !! radius_min == 0.0d0
                if (g_act == Generation_max+2) then
                    if (i.le.50) write(32,*) position(i,1), position(i,2), g_act !! Print final positions
                    if (i < 500) write(33,*) position(i,1), position(i,2) !! Print final positions
                    account(k) = account(k) + 1
                    !write(36,*) "***** Particle Attached! *****"; !write(36,*) "i = ",i; !write(36,*); !write(36,*); !write(36,*)
                    goto 1101 !! Jump the out while loop.
                else
                    g_act = g_act + 1
                end if
            end if


            !### [5] **Determine if go out** (not for start up)
            !write(36,*) "section[5]"
            do g_test = 2, g_act !! Search all passed level g. Make sure the particle not out of region at all level.
                int_gap = 4**(Generation_max - g_test + 2) ! 定义该级别搜索步长

                if (g_test.eq.2) then
                    k1 = rec_index(g_test-1,1) ! = 1 
                    k2 = rec_index(g_test-1,2) - int_gap ! = num_points - int_gap
                    infvl = k1 + int_gap; supvl = k2 - int_gap
                    search_length = polygon_shape - 2
                else
                    infvl = rec_index(g_test-1,1) !! Define search inf
                    supvl = rec_index(g_test-1,2) !! Define search sup
                    call cherche(g_test, infvl, supvl, int_gap, k1, k2, search_length) !! subroutine cherche(niveau, connu1, connu2, gap, nouveau1, nouveau2, combien)
                end if
                l1 = span(coord_x(k1), coord_y(k1), position(i,1), position(i,2))
                l2 = span(coord_x(k2), coord_y(k2), position(i,1), position(i,2))

                
                
                allocate(search_index(search_length))
                do j = 1, search_length
                    search_index(j) = noloop(infvl + (j-1) * int_gap)
                end do
                rec_old1 = rec_index(g_test,1); rec_old2 = rec_index(g_test,2) ! 该级别下原来最近两端点位置

                !write(36,*) "[5]: #rec(g_test-1,:) & int_gap =", rec_index(g_test-1,:), int_gap
                !write(36,*) "k1/k2 & l1/l2", k1, k2, l1, l2
                !write(36,*) "search length = ", search_length
                do j = 1, search_length
                    k = search_index(j)
                    u = span(coord_x(k), coord_y(k), position(i,1), position(i,2))
                    !write(36,*) "g_test =", g_test," search_index & span", k, u
                    call ordre(k1,k2, l1,l2, k,u) !! subroutine ordre(k1,k2,l1,l2, k,l)
                end do
                rec_index(g_test,1) = k1; rec_index(g_test,2) = k2
                deallocate(search_index)
                
                
                
                !if ((rec_old1 > rec_index(g_test,1)).or.(rec_old2 < rec_index(g_test,2))) then
                if ((noloop(rec_old1 - rec_index(g_test,1) + num_points - 1) > 0).or.&
                & (noloop(rec_index(g_test,2) - rec_old2 + num_points - 1) > 0)) then
                    !write(36,*) "Break @ g_act =", g_act
                    g_act = g_test !! 之后g要重新算，固定为首次出现异常的g_test
                    if (g_act > Generation_max+2) g_act = Generation_max
                else 
                    !! 范围可能缩小，于是重新赋值
                    rec_index(g_test,1) = rec_old1; rec_index(g_test,2) = rec_old2
                end if

            end do
            

            if (noloop(rec_index(g_act,2) - rec_index(g_act,1) + num_points - 1) == int_gap) then
                if (g_act < Generation_max + 2) g_act = g_act + 1
            end if
            
            
            !write(36,*) "Done once for while loop."
            !write(36,*) "i =",i,"g_act =",g_act
            !write(36,*) "#rec(1)", rec_index(:,1)
            !write(36,*) "#rec(2)", rec_index(:,2)
            !write(36,*); !write(36,*)

	
        end do 
        1101 continue
    end do
	


	!!! Normalize account
    do i = 1, num_points - 1
        account(i) = account(i) / num_particle
        write(34,*) i, ",", account(i), ","
    end do



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





    deallocate(position); deallocate(coord_x); deallocate(coord_y); deallocate(account); deallocate(zeta)
    deallocate(rec_index)
    close(31); close(32); close(33); close(34); close(35); close(36)


    call cpu_time(temps_fin)
    write(*,99) temps_fin-temps_debut !! Write the total time consumed to the screen.
	! close files?
end Program main





MODULE MATHS
    implicit none
    real(kind=8),parameter :: pi=4.0d0*atan(1.0d0),twopi=2.0d0*pi
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
	
    ruein="diff_input_v3.txt"
    open(unit=201,file=TRIM(ruein),form='formatted',status='old',action='read',iostat=status1,iomsg=msg)
	
    read(201,*) polygon_shape !! 3 = triangular; 4 = square; 5 = star; etc.
    read(201,*) Generation_max !! The maximum value of generation / recursion
    read(201,*) num_particle !! Number of particles
    read(201,*) Length !! Edge length for initial polygon
    read(201,*) alpha_angle_ratio !! Fractal angle 
    read(201,*) direct
	
    close(201)
end subroutine
subroutine koch_line(x1,y1,x2,y2,alpha,direction, b1,b2,c1,c2,d1,d2)
    use MATHS
    implicit none
	include "para.h"
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
    real*8 :: deltax, deltay, ell, eps, theta, st, ct, xdist, ydist, dxnew, dynew, d1, d2

    deltax = x2 - x1; deltay = y2 - y1
    ell = sqrt(deltax**2 + deltay**2); eps = ell / 10**3
    
    if (x1.eq.x2) then
        if (y1.le.y2) then
            theta = + pi / 2
        else
            theta = - pi / 2
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
        if ((dxnew.ge.0).and.(dxnew.le.ell)) then
            dmin = dynew; dehors = +1
            if (dmin.lt.eps) dmin = 0.0d0
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




subroutine cherche(niveau, connu1, connu2, gap, nouveau1, nouveau2, combien)
    !! Here is a procedure to extend the search space
    implicit none
    include "para.h"
    integer,intent(in) :: niveau, connu1, connu2, gap
    !! "niveau" = g_act/g_test;
    !! "connu1/2" refer to obvious search range, namely #rec(g_act-1,:)
    !! "gap" = distance of two point indices
    integer,intent(out):: nouveau1, nouveau2, combien
    !! "nouveau1/2" = two new point indices after extension.
    integer,external :: noloop

    if (niveau.eq.2) then
        nouveau1 = connu1 !+ gap
        nouveau2 = connu2 !- gap
        combien = (nouveau2 - nouveau1) / gap - 2 !! 不考虑首尾两点
    else
        !combien = (connu2 - connu1) / gap + 1
        combien = noloop(connu2 - connu1 + num_points - 1) / gap + 1
        nouveau1 = noloop(connu1 - gap + num_points - 1) !mod(connu1 - gap + num_points - 1, num_points - 1)
        nouveau2 = noloop(connu2 + gap + num_points - 1) !mod(connu2 + gap + num_points - 1, num_points - 1)
    end if
end subroutine
integer function noloop(index)
    implicit none
    include "para.h"
    integer :: index
    !noloop = mod(index,num_points-1)
    !if (noloop.eq.0) noloop = num_points - 1
    if (index >= num_points) then
        noloop = mod(index, num_points-1)
    else
        noloop = index
    end if
    return
end function noloop



subroutine ordre(k1,k2,l1,l2, k,l)
    !! Here is a procedure to find 2 nearest points. Added 03/05/2023
    implicit none
    include "para.h"
    integer,intent(in) :: k
    integer :: k1, k2
    integer :: kint
    real*8,intent(in) :: l
    real*8 :: l1, l2
    real*8 :: u, lint
    real*8,external :: noloop

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

    call random_number(u)
    random_angle = u * twopi
    move_x = rayon * cos(random_angle)
    move_y = rayon * sin(random_angle)

    x0 = x0 + move_x
    y0 = y0 + move_y
end subroutine