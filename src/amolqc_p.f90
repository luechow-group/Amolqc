! Copyright (C) 1996-2021 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

!
!               ----------------------
!                  A  M  O  L  Q  C
!               ----------------------
!
! Atomic and Molecular Quantum Monte Carlo Calculations
!
! initial version:
!  Arne Luechow, Penn State University, Feb. 1996
!
! main author:
!  Arne Luechow, RWTH Aachen University
!
! with contributions from:
! Sebastian Manten, Christian Diedrich, Annika Bande, Tony C. Scott,
! Annett Schwarz, Rene Petz, Raphael Berner, Alexander Sturm,
! Marko Hermsen, Kaveh Haghighi Mood, Christoph Schulte,
! Leonard Reuter, Michael A. Heuer, Jil Ludovicy
!
!
! main program



program amolqc_p

   use global_m
   use mainLoop_m
   use init_m, only: initAmolqc, finalizeAmolqc

   implicit none


   call initAmolqc()

   call mainloop()

   call finalizeAmolqc()

end program amolqc_p


