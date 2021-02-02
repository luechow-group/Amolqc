$begin_subroutine(name=jelm)
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
x$optimize_parameters(params=jastrow,energy_min,method=lm_newton,eq_iter=5,
nu=0.001,delta_f_min=0.0001,max_var=50.0,eq_call=jelm,write_wf)
$sample(change_size,new_size=1,last)
$qmc(vmc,steps=1500,block_len=500,accept_ratio=0.5,discard_all)
x$qmc(vmc,accumulate,step_stride=10,acc_size=1k,block_len=200,
accept_ratio=0.5,discard=1,move=umr)
$save_result()
$print_results()
