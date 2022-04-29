$begin_subroutine(name=jelin)
$sample(change_size,new_size=1,last)
x$qmc(accept_ratio=.5,
vmc,discard_all,steps=1500,block_len=500)
x$qmc(step_stride=10,acc_size=1k,block_len=200,accept_ratio=.5,discard=1,move=umr,
vmc,accumulate)
$save_result()
x$stop_if(variance>50.0)
$sample(remove_outliers,no_replace)
$end_subroutine()


$sample(change_size,new_size=1,last)
x$qmc(step_stride=10,acc_size=1k,block_len=200,accept_ratio=.5,discard=1,move=umr,
vmc,accumulate)
$sample(remove_outliers,no_replace)
x$optimize_parameters(target_E=REQUIRED,target_var=0.,eq_iter=5,max_ev=5,ev_sample_size=100,cffac=0.,optmode=1,
params=jastrow,energy_min,method=lin,eq_call=jelin,write_wf)
$sample(change_size,new_size=1,last)
x$qmc(accept_ratio=.5,
vmc,discard_all,steps=1500,block_len=500)
x$qmc(step_stride=10,acc_size=1k,block_len=200,accept_ratio=.5,discard=1,move=umr,
vmc,accumulate)
$save_result()
$print_results()

