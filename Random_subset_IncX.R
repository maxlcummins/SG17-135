IncX <- df3[df3$`IncX1_1 76/82 93%` == 1,] 

IncX_names <- rownames(IncX)

set.seed(1)

base::sample(IncX_names, size = 10) %>% sort()



