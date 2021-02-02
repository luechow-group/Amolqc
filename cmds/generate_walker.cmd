X$sample(create, start=density, size=10)
$qmc(vmc, steps=300, block_len=100, persist=9, accept_ratio=0.5, discard_all)
$sample(remove_outliers)
$sample(change_size, new_size=1)
x$qmc(vmc, steps=700, block_len=100, accept_ratio=0.5, discard_all)

