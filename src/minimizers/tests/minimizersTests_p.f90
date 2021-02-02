! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

program minimizersTests_p
use minimizer_tm
use minimizer_ws_tm
use myfctn_tm

implicit none

call minimizer_test()
call minimizer_ws_test()
call myfctn_test()

end program minimizersTests_p
