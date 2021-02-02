$begin_subroutine(name=jenw)
$sample(change_size,new_size=1,last)
$qmc(vmc,steps=1500,block_len=500,accept_ratio=0.5,discard_all)
x$qmc(vmc,accumulate,step_stride=10,acc_size=1k,block_len=200
accept_ratio=0.5,discard=1)
$save_result()
$sample(remove_outliers,no_replace)
$end_subroutine()


$sample(change_size,new_size=1,last)
x$qmc(vmc,accumulate,step_stride=10,move=umr,acc_size=1k,
accept_ratio=0.5,discard=1)
$sample(remove_outliers,no_replace)
x$optimize_parameters(params=jastrow,energy_min,method=newton,eq_iter=5
,max_var=50.0,eq_call=jenw,write_wf)
$sample(change_size,new_size=1,last)
$qmc(vmc,steps=1500,block_len=500,accept_ratio=0.5,discard_all)
x$qmc(vmc,accumulate,step_stride=10,move=umr,acc_size=1k,accept_ratio=0.5,
discard=1)
$save_result()
$print_results()
