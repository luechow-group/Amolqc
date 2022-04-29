x$sample(start=density,size=10,
create)
$qmc(steps=300,block_len=100,accept_ratio=.5,persist=9,
vmc,discard_all)
$sample(remove_outliers)
$sample(change_size,new_size=1)
x$qmc(steps=700,block_len=100,accept_ratio=.5,
vmc,discard_all)

