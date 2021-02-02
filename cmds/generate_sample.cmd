x$sample(create,start=density,generate=random,size=100)
$sample(remove_outliers)
$qmc(vmc,accept_ratio=0.5,steps=30,block_len=10,persist=9,discard_all)
x$qmc(vmc,accept_ratio=0.5,steps=400,block_len=10,discard_all)

