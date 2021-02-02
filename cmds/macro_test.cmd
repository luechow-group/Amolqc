! testing
$qmc(vmc,steps=20,blocks=3,tau=0.01,discard=3)
$change_jastrow(new_jastrow=sm3)
! testing
$wf(write,file='ne-t.wf')
qmc(vmc,steps=20,blocks=3,tau=0.01,discard=3)

