Program test
    implicit none
    !integer :: i,j(5)
    integer :: k1,k2,k
    real*8 :: l1,l2,l
    real*8 :: temps_debut, temps_fin!,u
    !integer :: date(8)

    call cpu_time(temps_debut)

    k1 = 1; k2 = 33; l1 = 496.39983773319386; l2 = 590.13500642954523
    k = 17; l = 656.47831490427484
    write(*,*) k1,k2,l1,l2, k,l
    call ordre(k1,k2,l1,l2, k,l)
    write(*,*) k1,k2,l1,l2, k,l
    
    
    
    
    call cpu_time(temps_fin)
    write(*,*) temps_fin - temps_debut

end Program test



recursive subroutine Gamma(x,v)
    implicit none
    integer :: x,v,w

    if (x == 0) v = 1
    if (x == 1) then
        v = 1
    else
        call Gamma(x-1,w)
        v = x * w
    end if
end subroutine





subroutine ordre(k1,k2,l1,l2, k,l)
    !! Here is a procedure to find 2 nearest points. Added 03/05/2023
    implicit none
    !include "para.h"
    integer :: k1, k2, k, kint
    real*8 :: l1, l2, l, u, lint

    !! 首先判断已知两点是否等距，若是则随机放弃一个
    if (l1 == l2) then
        if (l < l1) then
            call random_number(u)
            if (u < 0.5d0) then
                k1 = k; l1 = l
            else
                k2 = k1; l2 = l
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

            if (l == l2) then
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
    !if (((k1 > k2).and.((k1-k2) < num_points/2)).or.((k1<k2).and.((k2-k1) > num_points/2))) then
    !    kint = k1; k1 = k2; k2 = kint
    !    lint = l1; l1 = l2; l2 = lint
    !end if
end subroutine