!!!!! Diffusion near arbitrary Koch
!!!!! Yilin YE @ 04/2023


Program main
	use MATHS
	implicit none
	include "para.h"
	integer :: i,j,k,l ! Loop variables
	integer :: outpoint(2), intpoint(3), out_gap, int_gap, rec_index(6)
	integer :: num_points, qmax, g_last, g_act
	real*8 :: angle_ini, tta, new_x, new_y, radius_min, radius, u, random_angle, move_x, move_y
	real*8 :: rec_dist(6), dist_old, dist_new
	integer,allocatable :: search_coord(:)
	real*8,allocatable :: position(:,:), coord_x(:), coord_y(:), account(:), zeta(:)
	real*8,external :: distance, span
	
	integer :: date(8)
	real*8 :: time_begin, time_end, temps_debut, temps_fin !! Define variables to record time consumed
	character*128 :: rue_total !! Define path name for outputs
	




	call cpu_time(temps_debut)
	continue
		94 format('PSD @ YYE. Simulation done on ',I4,2('-',I2),'   ',I2,'h',I2)
		99 format('It spent',f10.4,' seconds for the whole program.')
	continue



	call obtaindata
	!write(*,*) polygon_shape, Generation_max, num_particle, Length, alpha_angle_ratio
	call date_and_time(values=date)
	rue_total = "diff-"
	open(unit=31,file=TRIM(rue_total)//"coord.txt"); write(31,94) date(1), date(2), date(3), date(5), date(6)
	open(unit=32,file=TRIM(rue_total)//"particle_traj.txt"); write(32,94) date(1), date(2), date(3), date(5), date(6)
	open(unit=33,file=TRIM(rue_total)//"particle_final.txt"); write(33,94) date(1), date(2), date(3), date(5), date(6)
	open(unit=34,file=TRIM(rue_total)//"sgm_pbb.txt"); write(34,94) date(1), date(2), date(3), date(5), date(6)
	open(unit=35,file=TRIM(rue_total)//"zeta.txt"); write(35,94) date(1), date(2), date(3), date(5), date(6)
	num_points = polygon_shape * 4**Generation_max + 1
	allocate(position(num_particle,2)); allocate(coord_x(num_points)); allocate(coord_y(num_points))
	allocate(account(num_points-1)) !! account = count @Python
	position = 0.0d0; coord_x = 0.0d0; coord_y = 0.0d0; account = 0.0d0


	
	!!! Generate initial polygon coordinates
	angle_ini = twopi / polygon_shape; tta = pi * (1.0d0 - 2.0d0/polygon_shape) 
	!!!!! Note, all integers in Python should be converted as real values for +-*/, else mod taken...
	new_x = - Length/2.0d0; new_y = new_x * tan(tta/2)
	do i = 0, polygon_shape - 1
		!coordinate.append((new_x,new_y))
		j = i * 4**Generation_max + 1
		coord_x(j) = new_x; coord_y(j) = new_y
        new_x = new_x + Length * cos(i * angle_ini)
        new_y = new_y + Length * sin(i * angle_ini)
	end do
	coord_x(num_points) = coord_x(1); coord_y(num_points) = coord_y(1) !coordinate.append(coordinate[0]) # To form a circle.
	
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
	
	!!! Output fractal points
	do i=1,num_points
		write(31,*) coord_x(i), coord_y(i)
	end do
	


	!!! Diffusion for each particle
	do i = 1, num_particle
		!! Initialization
		if (i.le.50) write(32,*) "#    ",i !! Furnish particle index in traj.
		g_act = 0; g_last = -1; rec_dist = Length
		rec_index(1) = 1; rec_index(2) = num_points
		rec_index(3) = 1; rec_index(4) = num_points
		rec_index(5) = 1; rec_index(6) = num_points
		!! Iteration for GAFRW
		do while (1.gt.0)
			!! Determine the search space by start/end points' indices
			if (g_last.eq.-1) then
				!! Find 2 nearest points
			else
				if (g_act.lt.g_last) then
					rec_dist(3) = rec_dist(1); rec_index(3) = rec_index(1)
					rec_dist(4) = rec_dist(2); rec_index(4) = rec_index(2)
					goto 1102
				else
					!! Find 2 nearest points
				end if
			end if

			!! Find 2 nearest points in given g
			out_gap = 4**(Generation_max - g_act); int_gap = out_gap / 4
			k = (rec_index(4) - rec_index(3)) / out_gap
			allocate(search_coord(k))
			do j = rec_index(3), rec_index(4), out_gap !0, k-1
				!l = out_gap * j
				!search_coord(j) = rec_index(3) + l
				u = span(coord_x(search_coord(j)), coord_y(search_coord(j)), position(i,1), position(i,2))
				!! Define rec_index(4): (g_last argmin1, g_last argmin2, g_act argmin1, g_act argmin2)
				!! Define rec_dist(4): (g_last min1, g_last min2, g_act min1, g_act min2)
				call vertex_update(u, search_coord(j), rec_dist(3), rec_dist(4), rec_index(3), rec_index(4))
				call petitavant(rec_index(3), rec_index(4))
			end do
			1102 continue


			radius_min = Length
			do j = 1, num_points-1
				radius = distance(coord_x(j), coord_y(j), coord_x(j+1), coord_y(j+1), position(i,1), position(i,2))
				if (radius_min.ge.radius) then
					radius_min = radius
					k = j
				end if
			end do
			
			
			
			
			if ((g_act.eq.Generation_max).and.(radius_min.eq.0)) then
				account(k) = account(k) + 1.0d0
				goto 1101
			else
				if (i.le.50) write(32,*) position(i,1), position(i,2) !! Record trajs of each particle
				!! Diffusion: GAFRW
				call random_number(u)
				random_angle = u * twopi
				move_x = radius_min * cos(random_angle)
				move_y = radius_min * sin(random_angle)

				dist_old = span(0.0d0,0.0d0,position(i,1),position(i,2))
				position(i,1) = position(i,1) + move_x
				position(i,2) = position(i,2) + move_y
				dist_new = span(0.0d0,0.0d0,position(i,1),position(i,2))
			end if
			
			!! Find 2 nearest points. Make sure particles stay in the same region.
			!! -> rec(5,6)
			g_last = g_act

			!! Update g_act 
			if ((rec_index(3).le.rec_index(5)).and.(rec_index(6).le.rec_index(4))) then
				if (dist_new.gt.dist_old) then
					g_act = g_act + 1
					rec_dist(1) = rec_dist(3); rec_index(1) = rec_index(3)
					rec_dist(2) = rec_dist(4); rec_index(2) = rec_index(4)
					if (g_act.gt.Generation_max) g_act = Generation_max
				end if
				if (dist_new.lt.dist_old) then
					g_act = g_act - 1
					if (g_act.lt.0) g_act = 0
				end if
				! if (dist_new.eq.dist_old) g_act = g_act
			! else g_act = g_act
			end if	


			
		end do 
		1101 continue
		!write(*,*) "Done for particle",i
	end do
	
	!!! Print final positions
	do i = 1, num_particle
		write(33,*) position(i,1), position(i,2)
	end do

	!!! Normalize account
	do i = 1, num_points - 1
		account(i) = account(i) / num_particle
		write(34,*) i, account(i)
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
	close(31); close(32); close(33); close(34); close(35)


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
	
	ruein="diff_input.txt"
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





real*8 function distance(x1, y1, x2, y2, x0, y0)
	use MATHS
	implicit none
	include "para.h"
	real*8 :: x1, y1, x2, y2, x0, y0 !! Given 2 points: (x1,y1) & (x2,y2); Particle: (x0,y0)
	real*8 :: deltax, deltay, ell, eps, theta, st, ct, xdist, ydist, dxnew, dynew, dmin, d1, d2

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
        dmin = Length
    else
        if ((dxnew.ge.0).and.(dxnew.le.ell)) then
            dmin = dynew
			if (dmin.lt.eps) dmin = 0.0d0
        else
            !d1 = sqrt((x0-x1)**2 + (y0-y1)**2)
			d1 = sqrt(xdist**2 + ydist**2)
            d2 = sqrt((x0-x2)**2 + (y0-y2)**2)
            if (d1.le.d2) then
                dmin = d1
            else
                dmin = d2
			end if
        end if
	end if
	distance = dmin
	return
end function distance





real*8 function span(x1, y1, x2, y2)
	implicit none
	real*8 :: x1, y1, x2, y2
	span = sqrt((x1-x2)**2 + (y1-y2)**2)
	return
end function span





subroutine vertex_update(valeur, arg, min1, min2, argmin1, argmin2)
	!! Here is a procedure to update possible change of searching range.
	implicit none
	real*8, intent(in) :: valeur
	integer, intent(in) :: arg
	real*8 :: min1, min2
	integer:: argmin1, argmin2, argint
	if (valeur.le.min1) then
		min2 = min1; argmin2 = argmin1
		min1 = valeur; argmin1 = arg
	else
		if (valeur.le.min2) then
			min2 = valeur; argmin2 = arg
		end if
	end if
end subroutine
subroutine petitavant(v1,v2)
	!! Here is a procedure to make sure the first value not large than the second one.
	implicit none
	integer :: v1,v2,vint
	if (v1.gt.v2) then
		vint = v1
		v1 = v2
		v2 = vint
	end if
end subroutine