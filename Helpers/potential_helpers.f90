program potential_helpers
    !
    ! Test Program for S and T Integrals.
    !
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    real(dp) :: a_p, time, s, t

    a_p = 5.202545d0
    time = 3.0d0

    call s_and_t(a_p, time, s, t)

    write(*,*)a_p,time,s,t

end program


subroutine s_and_t(a_p_in, time, s, t)
    !
    ! Integrates the S and T terms from the gas disk potential. 
    ! Given the Gas Model in Morisihima+ (2010).
    !
    ! Cf. Nagasawa+ (2000), Appendix, Eq. (A3)
    !     http://adsabs.harvard.edu/abs/2000AJ....119.1480N
    ! Cf. Morishima+ (2010), Section 2.2
    !     http://adsabs.harvard.edu/abs/2010Icar..207..517M
    ! 
    ! Generate Python Binding:
    ! $ f2py -c -m potential_helpers potential_helpers.f90
    ! $ ipython
    ! > import potential_helpers as ph
    ! > ph.s_and_t(1.0, 10.0)
    !
    ! Test w/ Jupiter:
    ! $ gfortran -o potential_helpers potential_helpers.f90
    ! $ ./potential_helpers
    !
    ! Or use the Makefile:
    ! $ make ph
    !
    ! Original Code by Ryuiji Morishima.
    ! 
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    real(dp), intent(in)  :: a_p_in ! Planet Location (AU)
    real(dp), intent(in)  :: time   ! Decay Time (Units of Tau)
    real(dp), intent(out) :: s, t

    real(dp) :: sigma_0, h_0, alpha, r_in, r_out ! Gas Stuff
    real(dp) :: a_p ! Planet Parameters
    real(dp) :: a_1, a_2

    ! Parameters
    real(dp), parameter :: PI = 3.14159265359d0
    real(dp), parameter :: G  = 6.674d-8          ! cm3/g/s2
    real(dp), parameter :: AU = 1.49598d13        ! cm
    
    ! Grid Stuff
    integer :: nr1max = 100, nr2max = 100, nr3max = 200
    integer :: nz1max = 10, nz2max = 100
    integer :: np1max = 10, np2max = 100
    integer :: nr, nz, np 
    integer :: nrmax, npmax
    integer :: nzmax
    integer :: nr1m, nr2m, nr3m 
    real(dp) :: nnr, nnz, nnp, nnr1m, nnr2m, nnr3m, nnp1m, nnp2m
    real(dp) :: nnz1m, nnz2m
    real(dp) :: dnr, dnr1, dnr2, dnr3, dnp, dnp1, dnp2
    real(dp) :: dnz, dnz1, dnz2, z_1
    real(dp) :: frho0, frho, rho
    real(dp) :: r, p, z, h, p_1

    ! Summing Stuff
    real(dp) :: s1, s2, t1, t2
    real(dp) :: rpa, ar2, arz, cp, c1, c32, c52

    ! Gas Stuff
    real(dp) :: sigma

    ! Gas Stuff
    sigma_0 = 2000.0d0 ! g/cm2
    r_in = 0.1d0       ! AU
    r_out = 35.0d0     ! AU
    alpha = 1.0d0
    h_0 = 0.03358d0    ! Scale Height @ 1 AU [AU]

    ! Conversions
    r_in = r_in * AU
    r_out = r_out * AU
    a_p = a_p_in * AU
    h_0 = h_0 * AU

    ! Decay Gas
    sigma_0 = sigma_0 * exp(-time)

    ! Some Gas Constant
    ! Factor 2*2 comes from integrating a quarter of the space only.
    ! (0 to 180 degree and 0 to zmax; x4 for total integral)
    frho0 = 4.0d0/sqrt(2.0*PI)

    ! Grid Stuff
    nr1m = nr1max
    nr2m = nr2max
    nr3m = nr3max

    nnr1m = nr1max
    nnr2m = nr2max
    nnr3m = nr3max

    nnz1m = nz1max
    nnz2m = nz2max

    nnp1m = np1max
    nnp2m = np2max

    ! High-Res Grid on Planet
    a_1 = 0.98d0 * a_p
    a_2 = 1.02d0 * a_p

    ! Is the planet actually further out than the inner disk edge?
    if(r_in .lt. a_1)then
        dnr1 = (a_1 - r_in) / nnr1m 
    else
        nr1m  = 0
        nnr1m = 0
        dnr1  = 0.0d0
    endif

    ! Is the planet away from the edges?
    if(r_out .gt. a_1 .and. r_in .lt. a_2) then
        dnr2 = (a_2 - a_1) / nnr2m
    else
        nr2m  = 0
        nnr2m = 0
        dnr2  = 0.0d0
    endif

    ! Is the outer edge further out than the planet?
    if(r_out .gt. a_2)then
        dnr3 = (r_out - a_2) / nnr3m
    else
        nr3m  = 0
        nnr3m = 0
        dnr3  = 0.0d0
    endif
 
    nrmax = nr1m + nr2m + nr3m
    nzmax = nz1max + nz2max

    ! Two Grids in Angle
    p_1 = 0.03d0 * PI
    dnp1 = p_1 / nnp1m         ! 0 to 0.03 PI 
    dnp2 = (PI - p_1) / nnp2m  ! 0.03 PI to PI
    npmax = np1max + np2max

    t = 0.d0
    s = 0.d0

    ! Loop Radius
    do nr = 1, nrmax

        ! Compute r (cell-centered), dr
        nnr = nr
        if(nr.le.nr1m)then
            ! 0.1 AU to 0.98 * a_planet
            r = r_in + (nnr - 0.5) * dnr1
            dnr = dnr1
        elseif(nr .gt. nr1m .and. nr .le. (nr1m + nr2m))then
            ! 0.98 * a_planet to 1.02 * a_planet
            r = a_1 + ((nnr - nnr1m) - 0.5) * dnr2
            dnr = dnr2
        elseif(nr.gt.(nr1m+nr2m))then
            ! 1.02 * a_planet to 35.0 AU
            r = a_2 + ((nnr - (nnr1m + nnr2m)) - 0.5) * dnr3
            dnr = dnr3
        endif

        ! Disk Density, Vertically Isothermal (h_0 * r * r**0.25)
        sigma = sigma_0*(r/AU)**(-alpha)    ! g/cm2
        h = h_0*(r/AU)**(1.25D0)            ! cm

        ! write(*,*)nr, dnr, r/au, sigma, h/AU, h_0

        ! Two Grids in Vertical
        z_1 = 0.03D0 * 0.5D0 * r
        dnz1 = z_1 / nnz1m              ! z = 0 to 0.03 * 0.5 r 
        dnz2 = (0.5D0*r - z_1) / nnz2m  ! z = 0.03 * r to 0.5 r 

        ! Gas Density
        frho = frho0 * dnr * r * sigma/h
        ! write(*,*)nr,dnr,r,sigma/h,frho

        ar2 = 2.0d0 * a_p * r
        rpa = r/a_p

        ! Loop Vertical
        do nz = 1, nzmax

            ! Compute z (cell centered), dz
            nnz = nz
            if(nz .le. nz1max)then
                ! z = 0 to 0.03 * 0.5 * r
                z = (nnz - 0.5D0) * dnz1
                dnz = dnz1
            else
                ! z = 0.03 * 0.5 * r to 0.5 * r
                z = ((nnz - nnz1m) - 0.5D0) * dnz2 + z_1
                dnz = dnz2
            endif
            ! write(*,*)nr,nz,z

            ! Compute Density
            rho = frho * dnz * exp(-0.5d0 * (z/h)**2.0d0)
            ! write(*,*)nr,nz,rho

            arz = a_p**2.0d0 + r**2.0d0 + z**2.0d0
            ! write(*,*)nr,nz,arz

            ! Loop Angle
            do np = 1, npmax

                ! Compute phi (cell centered), nphi
                nnp = np
                if(np .le. np1max)then
                    ! 0 to 0.03 PI 
                    p = (nnp - 0.5d0) * dnp1
                    dnp = dnp1
                else
                    ! 0.03 PI to PI
                    p = p_1 + ((nnp -nnp1m) - 0.5d0) * dnp2 
                    dnp = dnp2
                endif

                cp = cos(p)
                c1 = (arz - ar2*cp)/a_p**2.0d0
                ! c1 = (arz - ar2*cp)
                c32 = 1.0 / (c1 * sqrt(c1))
                c52 = c32 / c1
                cp = cp * rpa

                s = s + (-cp*c32 + 3.D0*(z/a_p)**2*c52)*(rho*dnp) 
                t = t + ((-3.D0 + 2.D0*cp)*c32 + &
                         3.D0*(1.D0-cp)**2*c52)*(rho*dnp)

            end do

        end do

        !write(*,*)nr,r/AU,-s*0.25*G
        
    end do

    ! write(*,*)s,t
    s = -0.25d0 * G * s / a_p
    t =  0.25d0 * G * t / a_p
    ! write(*,*)s,t

end subroutine
