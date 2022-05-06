$begin_subroutine(name=jvs1)
$sample(change_size,new_size=1,last)
x$qmc(accept_ratio=.5,
vmc,discard_all,steps=1500,block_len=500)
x$qmc(step_stride=10,acc_size=1k,block_len=200,accept_ratio=.5,discard=1,move=umr,
vmc,accumulate)
$sample(remove_outliers,no_replace)
$end_subroutine()


$begin_subroutine(name=jvs2)
$sample(change_size,new_size=1,last)
x$qmc(accept_ratio=.5,
vmc,discard_all,steps=1500,block_len=500)
x$qmc(step_stride=10,acc_size=1k,block_len=200,accept_ratio=.5,discard=1,move=umr,
vmc,accumulate)
$save_result()
$sample(remove_outliers,no_replace)
$end_subroutine()


$sample(change_size,new_size=1,last)
x$qmc(step_stride=10,acc_size=1k,block_len=200,accept_ratio=.5,discard=1,move=umr,
vmc,accumulate)
$sample(remove_outliers,no_replace)
x$optimize_parameters(E_ref=REQUIRED,eq_iter=5,optmode=1,
params=jastrow,variance_min,method=varmin,E_ref_adp,NL2SOL_D_mode=1,max_iter=1,eq_call=jvs1)
$sample(change_size,new_size=1,last)
x$qmc(step_stride=10,acc_size=1k,block_len=200,accept_ratio=.5,discard=1,move=umr,
vmc,accumulate)
x$optimize_parameters(E_ref=REQUIRED,eq_iter=5,optmode=1,
params=jastrow,variance_min,method=varmin,E_ref_adp,NL2SOL_D_mode=0,max_iter=1,write_wf,eq_call=jvs2)
$sample(change_size,new_size=1,last)
x$qmc(accept_ratio=.5,
vmc,discard_all,steps=1500,block_len=500)
x$qmc(step_stride=10,acc_size=1k,block_len=200,accept_ratio=.5,discard=1,move=umr,
vmc,accumulate)
$save_result()
$print_results()

