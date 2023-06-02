!!!!! Diffuion inside circle for validation.
!!!!! Yilin YE @ 06/2023

Program main
    use MATHS
    implicit none
    integer :: i,j,k
    integer :: num_points, num_particle, traj_max, final_max, inside, touch, n_iter
    real*8 :: angle_ini, tta, taille, radius_min, u, eps, ell, ell_hat, time_diff
    real*8 :: rho_ratio, q_rect, diff_coef, d, delta, rint, xini, yini, theta0, r0
    real*8,allocatable :: position(:,:), coord_x(:), coord_y(:), account(:), pbb_ana(:)
    integer,external:: noloop !! Function to return reasonnable segment index
    real*8,external :: span,pente !! Function to calculate distance between two points.

    !! Define system varaibles
    integer :: date(8)
    real*8 :: temps_debut, temps_fin, temps_provisoire
    character*128 :: rue_total


    call cpu_time(temps_debut)
    continue !! formats
        62 format(10X,'X coord',10X,'Y coord',10X,'Radius',10X,'Time',10X,'Local Time')
        63 format(12X,'X coord',12X,'Y coord',12X,'Time',12X,'Local Time',12X,'Threshold')
        64 format('#',3X,'Segment Index',3X,'Hitting Probability')
        71 format(2(3X,ES16.6))
        72 format(2(3X,f16.8),3(2X,f16.8))
        73 format(5(4X,f16.8))
        74 format(2X,I6,10X,ES16.8)
        91 format('# Polygon',I10,', Circle Radius',f8.2,', Thickness ',ES8.2)
        92 format('# No. Particle',I9,', Rho ratio',f6.2,', Reactivity ',ES8.2,', Diff coef ',ES8.2)
        94 format('Circle Test @ YYE. Version 1.0.0.. Simulation on ',I4,2('-',I2),'   ',I2,'h',I2)
        95 format('Done for particle',I9,', Time used (s)',ES11.3)
        99 format('It spent',f10.4,' seconds for the whole program.')
    continue


    
    
    !! Define boundary parameters
    num_points = 61
    num_particle = 5000000
    traj_max = 5
    final_max = 100
    allocate(position(num_particle,2))
    allocate(coord_x(num_points))
    allocate(coord_y(num_points))
    allocate(account(num_points-1))
    coord_x = 0.0d0
    coord_y = 0.0d0
    account = 0.0d0
    
    
    !! Determine diffusion/reflection parameters
    taille = 1.0d0
    eps = 1.0d-3
    rho_ratio = 2.0d0
    q_rect = 1.0d1
    diff_coef = 1.0d0
    
    
    !! Generate data files.
    call date_and_time(values=date)
    rue_total = "cc-"
    open(unit=31,file=TRIM(rue_total)//"coord.txt")
        write(31,94) date(1), date(2), date(3), date(5), date(6)
    open(unit=32,file=TRIM(rue_total)//"particle_traj.txt")
        write(32,94) date(1), date(2), date(3), date(5), date(6)
        write(32,62); write(32,*)
    open(unit=33,file=TRIM(rue_total)//"particle_final.txt")
        write(33,94) date(1), date(2), date(3), date(5), date(6)
        write(33,63); write(33,*)
    open(unit=34,file=TRIM(rue_total)//"pbnum.txt")
        write(34,94) date(1), date(2), date(3), date(5), date(6)
        write(34,91) num_points-1, taille, eps
        write(34,92) num_particle, rho_ratio, q_rect, diff_coef; write(34,*)
        write(34,64); write(34,*)
    continue


    !! Generate boundary points
    angle_ini = twopi/(num_points - 1)
    do i = 1, num_points
        j = i - 1
        tta = angle_ini * j
        coord_x(i) = taille * cos(tta)
        coord_y(i) = taille * sin(tta)
        write(31,71) coord_x(i), coord_y(i)
    end do


    !! Compute analytical probability.
    xini = 0.2d0; yini = 0.0d0
    theta0 = pente(xini,yini); r0 = span(0.0d0,0.0d0,xini,yini)
    open(unit=35,file=TRIM(rue_total)//"pbana.txt")
        write(35,94) date(1), date(2), date(3), date(5), date(6)
        write(35,92) num_particle, rho_ratio, q_rect, diff_coef
        write(35,*) "r0, θ0 = ", r0, theta0; write(35,*)

        write(34,*) "r0, θ0 = ", r0, theta0; write(34,*)
    allocate(pbb_ana(num_points-1))
    n_iter = 50
    do k = 1, num_points-1
        pbb_ana(k) = 1.0d0 / (taille * (num_points-1))
        do j = 1, n_iter
            pbb_ana(k) = pbb_ana(k) + 1.0d0/(pi * taille) * (r0/taille)**j * (sin(j * (twopi*k/(num_points-1)-theta0)) &
            & - sin(j * (twopi*(k-1)/(num_points-1)-theta0))) / (1.0d0/j+1.0d0/(q_rect*taille))
        end do
        write(35,74) k, pbb_ana(k)
    end do
    deallocate(pbb_ana)
    close(35)




    !! Diffusion for all particles
    do i = 1, num_particle
        
		!! Initialization
        call random_number(u)
        ell_hat = -log(u) / q_rect !! PDF needed. psi(l) = q*exp(-q*l); 
        ell = 0.0d0; time_diff = 0.0d0
        inside = 0; touch = 0
        position(i,1) = xini; position(i,2) = yini
        if (i <= traj_max) then
            write(32,*); write(32,*) "#    ",i
            write(32,72) position(i,1), position(i,2), radius_min, time_diff, ell
        end if


        !! Non-stop diffusion
        do while (1 > 0)
            radius_min = taille - span(0.0d0,0.0d0,position(i,1),position(i,2))
            if (radius_min > eps) then !! Not enter diffusion layer
                call diffusion(position(i,1), position(i,2), radius_min)
                if (radius_min > taille) then
                    write(32,*) "Leave the circle !!"
                    goto 1101
                end if
            else
                inside = 1
                radius_min = eps * rho_ratio
                call diffusion(position(i,1), position(i,2), radius_min)

                d = span(0.0d0,0.0d0,position(i,1),position(i,2))
                if (d < taille - eps) then !! Out of diffusion layer
                    inside = 2
                    !touch = 0
                else if (d > taille) then !! Need reflection!
                    touch = touch + 1

                    if (position(i,1) == 0.0d0) then                      
                        if (position(i,2) > 0.0d0) then
                            position(i,2) = +2.0d0 * taille - position(i,2)
                        else
                            position(i,2) = -2.0d0 * taille - position(i,2)
                        end if
                    else !! 不解方程，直接旋转到x轴正向，计算距离后旋转回去
                        tta = pente(position(i,1), position(i,2))
                        !call rotation(position(i,1), position(i,2), -tta, xint, yint)
                        rint = 2.0d0 * taille - d !span(0.0d0,0.0d0,position(i,1), position(i,2))
                        position(i,1) = rint * cos(tta)
                        position(i,2) = rint * sin(tta)
                    end if
                end if


                if (touch > 0) then
                    delta = radius_min**2 / (4.0d0 * diff_coef)
                    time_diff = time_diff + delta
                    ell = ell + sqrt(pi * diff_coef * delta / 2.0d0)
                    if (ell >= ell_hat) then
                        if (i <= traj_max) write(32,72) position(i,1), position(i,2), radius_min, time_diff, ell
                        goto 1101
                    end if
                end if    
                
                !! If leave diffusion layer, then clear "touch" and reset "inside"
                if (inside == 2) then
                    touch = 0
                    inside = 0
                end if
            end if


            if (i <= traj_max) write(32,72) position(i,1), position(i,2), radius_min, time_diff, ell
        end do
        

        1101 continue

        tta = pente(position(i,1),position(i,2))
        u = (num_points-1) * tta / twopi
        k = ceiling(u) !! 向上取整，角度趋近0则第一边概率++
        account(k) = account(k) + 1


        if (i <= final_max) then
            write(33,73) position(i,1), position(i,2), time_diff, ell, ell_hat
            write(32,*) "Local Time Threshold =", ell_hat
            write(32,*)
        end if
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





    deallocate(position); deallocate(coord_x); deallocate(coord_y); deallocate(account)
    close(31)


    call cpu_time(temps_fin)
    write(*,99) temps_fin-temps_debut !! Write the total time consumed to the screen.
    !write(30,99) temps_fin-temps_debut
    !close(30)


end Program main






MODULE MATHS
    implicit none
    real(kind=8),parameter :: pi = 4.0d0*atan(1.0d0), twopi = 2.0d0*pi
END MODULE MATHS



real*8 function span(x1,y1,x2,y2)
    implicit none
    real*8 :: x1,y1,x2,y2
    span = sqrt((x1-x2)**2 + (y1-y2)**2)
    return
end function span



integer function noloop(index,num_points)
    implicit none
    integer :: index, num_points
    if (index >= num_points) then
    noloop = mod(index, num_points - 1)
    else
    if (index <= 0) index = index + num_points - 1
    noloop = index
    end if
    return
end function noloop



subroutine diffusion (x0,y0,rayon)
    use MATHS
    implicit none
    real*8,intent(in) :: rayon
    real*8 :: x0,y0
    real*8 :: u, random_angle, move_x, move_y

    call random_number(u)
    random_angle = u * twopi
    move_x = rayon * cos(random_angle)
    move_y = rayon * sin(random_angle)
    x0 = x0 + move_x
    y0 = y0 + move_y
end subroutine



real*8 function pente(x0,y0)
    use MATHS
    implicit none
    real*8 :: x0, y0, theta

    if (x0 == 0.0d0) then
        if (y0 > 0.0d0) then
            theta = pi / 2.0d0
        else
            theta = pi / 2.0d0 * 3.0d0
        end if
    else
        theta = atan(y0/x0)
        if (theta < 0.0d0) then
            if (x0 < 0.0d0) then
                theta = theta + pi
            else
                theta = theta + twopi
            end if
        else if (theta == 0.0d0) then
            if (x0 < 0.0d0) theta = pi
        else
            if (x0 < 0.0d0) theta = theta + pi
        end if
    end if

    pente = theta
    return
end function



subroutine rotation(x0,y0,angle,xnew,ynew)
    use MATHS
    implicit none
    real*8, intent(in) :: x0,y0,angle
    real*8,intent(out) :: xnew,ynew
    real*8,external :: span,pente
    real*8 :: theta, radius

    radius = span(0.0d0,0.0d0,x0,y0)
    theta = pente(x0,y0)

    xnew = radius * cos(theta + angle)
    ynew = radius * sin(theta + angle)
end subroutine
