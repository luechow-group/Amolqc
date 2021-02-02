! Copyright (C) 2020 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module coulombDensity_m
    use kinds_m, only: r8
    use electronDensity_m, only: rho
    use wfData_m, only: separated_electrons
    use normal_m, only: normal,  normal_ran
    use random_m, only: myran
    use wfData_m, only: atoms
    use globalUtils_m, only: iul, logmode
    use parsing_m, only: getinta
    use error_m, only: assert

    implicit none
    private
    public :: get_density_repulsion, is_coulomb_density_initialized, coulomb_density

    logical :: is_initialized = .False.
    integer :: mc_steps = 100
contains
    function get_density_repulsion(position) result(repulsion)
        real(r8), intent(in) :: position(3)
        real(r8) :: repulsion
        integer :: se_distribution(separated_electrons) ! separated_electrons distribution to cores

        real(r8) :: f, r(3), rr(3)
        integer :: step, d, index, i, j
        real(r8) :: g, gg

        repulsion = 0._r8

        if (mc_steps > 0) then
            j=1
            do i=1,SIZE(atoms)
                se_distribution(j:j+atoms(i)%pa-1) = i
                j = j + atoms(i)%pa
            end do
            do step = 1, mc_steps
                ! chosing atom randomly, weighted by separated electrons
                index = se_distribution(INT(separated_electrons * myran()) + 1)

                ! random 3x normal distributed vector in proximity to the atom
                r = atoms(index)%Get_position()
                do d=1,3
                    r(d) = r(d) + normal_ran()
                end do

                ! calculating the value of the probability distribution
                g = 0.0
                do i=1,SIZE(atoms)
                    gg = (1._r8 * atoms(i)%pa) / separated_electrons  ! be careful with integer division!
                    rr = atoms(i)%Get_position()
                    do d=1,3
                        gg = gg * normal(r(d)-rr(d))
                    end do
                    g = g + gg
                end do

                ! doing the actual sampling
                call rho(r, f, orbitals=separated_electrons/2)
                repulsion = repulsion + f/(g * NORM2(position-r))
            end do

            repulsion = repulsion / mc_steps

            ! to reset wfData entries
            call rho(position, f, orbitals=separated_electrons/2)
        end if
    end function get_density_repulsion

    subroutine coulomb_density(lines, nl)
        character(len=120), intent(in) :: lines(nl)
        integer, intent(in) :: nl
        integer :: in_mc_steps = -1, iflag

        call getinta(lines,nl,'mc_steps=',in_mc_steps,iflag)
        call assert(in_mc_steps >= 0, 'coulomb_density: invalid mc_steps number')

        if (logmode >=2) then
            write(iul,'(A)') "Coulomb density calculation is initialized"
            write(iul,'(A,I15)') "Setting Monte Carlo integration mc_steps to = ", in_mc_steps
        end if

        is_initialized = .True.

        mc_steps = in_mc_steps
    end subroutine coulomb_density

    function is_coulomb_density_initialized() result(bool)
        logical :: bool
        bool = is_initialized
    end function is_coulomb_density_initialized

end module coulombDensity_m