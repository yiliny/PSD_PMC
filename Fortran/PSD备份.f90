!!!!! Diffusion near arbitrary Koch
!!!!! Yilin YE @ 04/2023


Program main
	use MATHS
	implicit none
	include "para.h"
	integer :: i,j,k,l ! Loop variables
	integer :: outpoint(2), intpoint(3), out_gap, int_gap
	integer :: num_points, qmax, g_last, g_act, g_test, case, k1,k2,kint,supvl
	real*8 :: angle_ini, tta, new_x, new_y, radius_min, radius, u, random_angle, move_x, move_y, l1,l2,lint
	real*8,allocatable :: position(:,:), coord_x(:), coord_y(:), account(:), zeta(:)
	integer,allocatable :: rec_index(:,:)
	real*8,allocatable :: rec_dist(:,:)
	real*8 :: rec_old1, rec_old2
	real*8,external :: span
	
	integer :: date(8)
	real*8 :: time_begin, time_end, temps_debut, temps_fin !! Define variables to record time consumed
	character*128 :: rue_total !! Define path name for outputs
	




	call cpu_time(temps_debut)
	continue
		76 format(8(2X,I6),f14.5) !! format for file(36)=debug.txt
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
	open(unit=36,file=TRIM(rue_total)//"debug.txt"); write(36,94) date(1), date(2), date(3), date(5), date(6)
	num_points = polygon_shape * 4**Generation_max + 1 !! + 1 to form a circle
	allocate(position(num_particle,2)); allocate(coord_x(num_points)); allocate(coord_y(num_points))
	allocate(account(num_points-1)) !! account = count @Python
	position = 0.0d0; coord_x = 0.0d0; coord_y = 0.0d0; account = 0.0d0


	
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
	
	!!! Output fractal points
	do i=1,num_points
		write(31,*) coord_x(i), coord_y(i)
	end do
	
	
	
	allocate(rec_index(Generation_max + 2, 2)); allocate(rec_dist(Generation_max + 2, 2))
	!!! Diffusion for each particle
	do i = 1, num_particle
		
		!! Initialization
		if (i.le.50) write(32,*) "# ",i !! Furnish particle index in traj.
		do j = 1, Generation_max + 2
			rec_index(j,1) = 1
			rec_index(j,2) = num_points
			rec_dist(j,1) = span(coord_x(rec_index(j,1)), coord_y(rec_index(j,1)), position(i,1), position(i,2))
			rec_dist(j,2) = span(coord_x(rec_index(j,2)), coord_y(rec_index(j,2)), position(i,1), position(i,2))
		end do
		!### [1.1] **Initialize `g_last` & `g_act`** (for start-up)
		g_act = 1; g_last = 0 !! @ level g, g_act = g+1!
		int_gap = 4**(Generation_max - g_act + 1)
		
		
		!! Iteration for GAFRW
		do while (1.gt.0)
			if (g_last.ne.0) then
				!### [1.2] **Initialize `g_last` & `g_act`** (not for start-up)
				if (g_last.lt.Generation_max + 2) then
					g_last = g_last + 1
				else
					g_last = Generation_max + 2
				end if
				int_gap = 4**(Generation_max - g_last + 1)
				write(36,*) g_last, int_gap


				!### [2] **Find 2 nearest points** (not for start-up)
				rec_dist(g_last,:) = rec_dist(g_last-1,:)
				rec_index(g_last,:)=rec_index(g_last-1,:)
				l1 = rec_dist(g_last-1,1); l2 = rec_dist(g_last-1,2)
				k1 =rec_index(g_last-1,1); k2 =rec_index(g_last-1,2)
				if (l1.gt.l2) then
					lint = l1; kint = k1;
					l1 = l2; k1 = k2;
					l2 = lint; k2 = kint
				end if
				supvl = rec_index(g_last-1,2)

				do j = rec_index(g_last-1,1) + int_gap, supvl-1, int_gap
					u = span(coord_x(j), coord_y(j), position(i,1), position(i,2))
					write(36,*) "g_last=",g_last,"j&u=",j,u
					if (u.lt.l1) then
						l2 = l1; k2 = k1
						l1 = u; k1 = j
					else
						if (u.lt.l2) then
							l2 = u; k2 = j
						end if
					end if
				end do
				call petitavant(k1,k2)
				rec_index(g_last,1) = k1; rec_index(g_last,2) = k2
				
				if (g_last.eq.1) then
					if (abs(rec_index(g_last-1,2) - rec_index(g_last-1,1)).gt.int_gap) then
						rec_index(g_last,1) = rec_index(g_last-1,2) - int_gap
						rec_index(g_last,2) = rec_index(g_last-1,2)
					end if
				end if
			end if

			
			
			!### [3] **Determine diffusion radius** (able for start-up)
			radius_min = Length
			do j = rec_index(g_last,1), rec_index(g_last,2)-1, int_gap
				l = j + int_gap
				call distance(coord_x(j),coord_y(j), coord_x(l),coord_y(l), position(i,1), position(i,2), case, radius)
				if (case.ne.-1) then
					if (radius_min.gt.radius) then
						radius_min = radius
						k = j
					end if
				end if
			end do



			!### [4] **Diffusion** with given `radius_min`.
			if (radius_min.gt.0.0d0) then
				call random_number(u); random_angle = u * twopi
				move_x = radius_min * cos(random_angle); move_y = radius_min * sin(random_angle)
				position(i,1) = position(i,1) + move_x; position(i,2) = position(i,2) + move_y
			else
				if (g_last.eq.Generation_max+2) then
					account(k) = account(k) + 1
				else
					g_last = g_last + 1
				end if
				goto 1101
			end if

			do j = 1, Generation_max+2
				rec_dist(j,1) = span(coord_x(rec_index(j,1)),coord_y(rec_index(j,1)),position(i,1),position(i,2))
				rec_dist(j,2) = span(coord_x(rec_index(j,2)),coord_y(rec_index(j,2)),position(i,1),position(i,2))
			end do



			if (g_last.eq.0) goto 1101
			!### [5] **Determine if go out** (not for start up)
			write(36,*) "section[5]"
			do g_test = 1, g_last
				supvl = rec_index(g_test-1,2)
				rec_old1 = rec_index(g_test,1)
				rec_old2 = rec_index(g_test,2)
				int_gap = 4**(Generation_max - g_test + 1)
				l1 = rec_dist(g_test-1,1); l2 = rec_dist(g_test-1,2)
				k1 =rec_index(g_test-1,1); k2 =rec_index(g_test-1,2)
				if (l1.gt.l2) then
					lint = l1; kint = k1
					l1 = l2; k1 = k2
					l2 = lint; k2 = kint
				end if

				do j = rec_index(g_test-1,1)+int_gap, supvl-1, int_gap
					u = span(coord_x(j), coord_y(j), position(i,1), position(i,2))
					write(36,*) "g_test=", g_test,"j&u", j,u
					if (u.lt.l1) then
						l2 = l1; k2 = k1
						l1 = u; k1 = j
					else
						if (u.lt.l2) then
							l2 = u; k2 = j
						end if
					end if
				end do
				call petitavant(k1,k2)
				rec_index(g_test,1) = rec_old1; rec_index(g_test,2) = rec_old2

				if (g_test.eq.1) then
					if (abs(rec_index(g_test,2) - rec_index(g_test,1)) > int_gap) then
						rec_index(g_test,1) = rec_index(g_test-1,2) - int_gap
						rec_index(g_test,2) = rec_index(g_test-1,2)
					end if
				end if

				if ((rec_old1.gt.rec_index(g_test,1)).or.(rec_old2.lt.rec_index(g_test,2))) then
					g_last = g_test
					write(36,*) "Break @ g_last =", g_last
					goto 1101
				end if

			end do


			
			!! Determine the search space by start/end points' indices
			if ((g_last.eq.0).and.(g_act.eq.1)) then
				int_gap = 4**(Generation_max - g_last)
				!g_last = g_last + 1
				goto 1102 !! At first, search all g=0 points
			else
				int_gap = 4**(Generation_max - g_last + 1)
			end if
			!! Find 2 nearest points @ g, inside the interval vertex @ g-1, ∆ = 4^{Gmax - g_last}
			if (g_last.ge.Generation_max+1) goto 1102
			do j = rec_index(g_last,1), rec_index(g_last,2), int_gap
				!!!!! Note we start from rec(1,1/2) @ g=0
				u = span(coord_x(j), coord_y(j), position(i,1), position(i,2))
				!! Define rec(g,2) for each generation g = 0, 1, ..., gmax. 1/2 for interval left/right point.
				call vertex_update(u, j, rec_dist(g_act,1), rec_dist(g_act,2), rec_index(g_act,1), rec_index(g_act,2))
				call petitavant(rec_index(g_act,1), rec_index(g_act,2))
				write(*,*) i,j
			end do
			1102 continue

			!! Find the minimal radius for diffusion
			radius_min = Length/4**g_last
			do j = rec_index(g_act,1), rec_index(g_act,2)-1, int_gap ! distance for segments, thus w/ -1
				l = j + int_gap
				!radius = distance(coord_x(j), coord_y(j), coord_x(l), coord_y(l), position(i,1), position(i,2))
				call distance(coord_x(j), coord_y(j), coord_x(l), coord_y(l), position(i,1), position(i,2), case, radius)
				if (radius_min.ge.radius) then
					radius_min = radius
					k = j
				end if
			end do

			!! If not at gmax, then diffuse.
			if (radius_min.gt.0) then
				if (i.le.50) write(32,*) position(i,1), position(i,2) !! Record trajs of each particle
				call random_number(u); random_angle = u * twopi !! Diffusion: GAFRW
				move_x = radius_min * cos(random_angle); move_y = radius_min * sin(random_angle)
				position(i,1) = position(i,1) + move_x; position(i,2) = position(i,2) + move_y
			else
				if (g_last.eq.Generation_max+1) then !! End of the whole diffusion.
					if (i.le.50) write(32,*) position(i,1), position(i,2) !! Print final positions
					write(33,*) position(i,1), position(i,2) !! Print final positions
					account(k) = account(k) + 1.0d0
					goto 1101 !! Get out of Loop for Particle i
				else	
					g_last = g_act; g_act = g_act + 1
					goto 1103 !! End this while loop. Already attached to boundary, thus we consider the next g.
				end if
			end if

			!! At first, jump the return block
			if (g_last.eq.0) then
				g_last = g_act
				g_act = g_act + 1
				goto 1103
			end if

			!! Find 2 nearest points. Make sure particles stay in the same region @ all levels
			l = 1 !! Determine if #rec @ g_last are still the nearest points @ g_last - 1.
			do while (l.lt.g_act) !! Perhaps more than 1 step needed for g--. Loop until no change
				g_test = l; g_last = g_test + 1; out_gap = 4**(Generation_max - g_test + 1)
				rec_old1 = rec_index(g_last,1); rec_old2 = rec_index(g_last,2)
				do j = rec_index(g_test,1), rec_index(g_test,2)-1, out_gap
					u = span(coord_x(j), coord_y(j), position(i,1), position(i,2))
					call vertex_update(u, j, rec_dist(g_last,1), rec_dist(g_last,2), rec_index(g_last,1), rec_index(g_last,2))
					call petitavant(rec_index(g_last,1), rec_index(g_last,2))
				end do
				if ((rec_old1.le.rec_index(g_last,1)).and.(rec_old2.ge.rec_index(g_last,2))) then
					rec_index(g_last,1) = rec_old1; rec_index(g_last,2) = rec_old2
					l = l + 1
				else
					g_last = g_test
					g_act = g_last + 1
					goto 1103
				end if
			end do
			
			!! Update g_act if pass checks @ all levels.
			if (g_last.lt.Generation_max+1) then
				g_last = g_act; g_act = g_act + 1
			else
				g_last = Generation_max+1
			end if
			1103 continue

			write(36,76) i, g_last, rec_index(1,1), rec_index(1,2), rec_index(2,1), rec_index(2,2), &
			&rec_index(3,1), rec_index(3,2), radius_min
			write(36,*) rec_dist(1,1), rec_dist(1,2), rec_dist(2,1), rec_dist(2,2), rec_dist(3,1), rec_dist(3,2)

	
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
	deallocate(rec_index); deallocate(rec_dist)
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
        dmin = Length
		dehors = -1
    else
        if ((dxnew.ge.0).and.(dxnew.le.ell)) then
			dehors = +1
            dmin = dynew
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





subroutine vertex_update(valeur, arg, min1, min2, argmin1, argmin2)
	!! Here is a procedure to update possible change of searching range.
	implicit none
	real*8, intent(in) :: valeur
	integer, intent(in) :: arg
	real*8 :: min1, min2, minint
	integer:: argmin1, argmin2, argint

	if (min1.gt.min2) then
		minint = min1; min1 = min2; min2 = minint
		argint = argmin1; argmin1 = argmin2; argmin2 = argint
	end if

	if (valeur.lt.min1) then
		min2 = min1; argmin2 = argmin1
		min1 = valeur; argmin1 = arg
	else
		if (valeur.lt.min2) then
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
