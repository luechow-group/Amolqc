$begin_subroutine(name=jvf)
$sample(change_size,new_size=1,last)
$qmc(vmc,steps=1500,block_len=500,accept_ratio=0.5,discard_all)
x$qmc(vmc,accumulate,step_stride=10,acc_size=1k,block_len=200,
accept_ratio=0.5,discard=1,move=umr)
$save_result()
$sample(remove_outliers,no_replace)
$end_subroutine()


$sample(change_size,new_size=1,last)
x$qmc(vmc,accumulate,step_stride=10,acc_size=1k,block_len=200,
accept_ratio=0.5,discard=1,move=umr)
$sample(remove_outliers,no_replace)
x$optimize_parameters(params=jastrow,method=varmin,E_ref=REQUIRED,eq_iter=5
,variance_min,write_wf,eq_call=jvf,E_ref_adp)
$sample(change_size,new_size=1,last)
$qmc(vmc,steps=1500,block_len=500,accept_ratio=0.5,discard_all)
x$qmc(vmc,accumulate,step_stride=10,acc_size=1k,block_len=200,
accept_ratio=0.5,discard=1,move=umr)
$save_result()
$print_results()
