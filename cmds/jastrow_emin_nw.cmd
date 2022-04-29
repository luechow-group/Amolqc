$begin_subroutine(name=jenr)
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
x$optimize_parameters(eq_iter=5,max_var=50.,optmode=1,
params=jastrow,energy_min,method=nr,eq_call=jenr,write_wf)
$sample(change_size,new_size=1,last)
x$qmc(accept_ratio=.5,
vmc,discard_all,steps=1500,block_len=500)
x$qmc(step_stride=10,acc_size=1k,block_len=200,accept_ratio=.5,discard=1,move=umr,
vmc,accumulate)
$save_result()
$print_results()

