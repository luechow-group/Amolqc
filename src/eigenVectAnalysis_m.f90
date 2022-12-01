! Copyright (C) 2022 Michel V. Heinz
!
! SPDX-License-Identifier: GPL-3.0-or-later

module eigenVectAnalysis_m
    use kinds_m, only: r8, i4
    use parsing_m, only: getinta, getloga
    use global_m
    use fPsi2_m
    use atom_m
    use mpiInterface_m, only: myMPIReduceSumInteger
    use rwSample_m, only: rWSample, getSampleSize, getFirst, isNext, getNext, getWalker
    use psiMax_m, only: psimax
    use randomWalker_m, only: RandomWalker, pos, resetTo
    use singularityParticles_m, only: singularity_particles

    implicit none

    private
    public :: eigenVectAnalysis

contains
    subroutine eigenVectAnalysis(smpl, psimax_obj)
        type(RWSample), intent(inout) :: smpl
        type(psimax), intent(inout)   :: psimax_obj
        type(RandomWalker), pointer   :: rw
        integer                       :: i, lwork, info, yml, k
        real(r8)                      :: xx(3*getNElec())
        real(r8)                      :: H(SIZE(xx),SIZE(xx))
        real(r8)                      :: lambda(SIZE(xx)), work2(3*SIZE(xx)-1)

        !initialize variables
        rw => getFirst(smpl)

        !create output document
        open(newunit=yml, file = 'eigenvect_analysis.yml')
        write(yml, '(a)') '---'
        write(yml, '(a)') 'Structures:'

        !eigenvector analysis
        do
            !get sample
            call pos(rw, xx)

            !calculate Hessian
            call psimax_obj%calculateHessian(xx, H, 2E+5_r8)

            !get eigenvalues and -vectors
            lwork = 3*SIZE(xx)-1
            call DSYEV('V', 'U', SIZE(xx), H, SIZE(xx), lambda, work2, lwork, info)
            !call assert(info == 0, 'eigenVectAnalysis: Inversion failed!')

            !write eigenvalues to document
            write(yml, '(a)', ADVANCE = 'no') '  - Eigenvalues: ['
            do i = 1, SIZE(lambda)
                write(yml, '(es16.6e3)', ADVANCE = 'no') lambda(i)
                if (i /= SIZE(lambda)) write(yml, '(a)', ADVANCE = 'no') ','
            end do
            write(yml, '(a)')']'

            !write eigenvectors to document
            write(yml, '(a)') '    Eigenvectors: ['
            do i = 1, SIZE(lambda)
                do k = 1, SIZE(lambda)/3
                    write(yml,'(a, es16.6e3, a)', ADVANCE = 'no') '      [', H(3*k-2, i), ','
                    write(yml,'(es16.6e3, a)', ADVANCE = 'no') H(3*k-1, i), ','
                    write(yml,'(es16.6e3, a)') H(3*k, i), '],'
                    end do
                end do
            write(yml, '(a)')'      ]'

            !get next sample
            if (.not. isNext(smpl)) exit
            rw => getNext(smpl)
        end do

        !close document
        write(yml, '(a)') '...'
        close(yml)

    end subroutine eigenVectAnalysis

end module eigenVectAnalysis_m