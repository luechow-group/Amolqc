x$sample(start=density,generate=random,size=100,
create)
$sample(remove_outliers)
$qmc(steps=30,block_len=10,accept_ratio=.5,persist=9,
vmc,discard_all)
x$qmc(steps=400,block_len=10,accept_ratio=.5,
vmc,discard_all)

