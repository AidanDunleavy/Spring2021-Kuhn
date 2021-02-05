D <- expand.grid( flavor = c(1, 2, 3),
                  sugarlevel = c("reg", "sug"),
                  block = c(1,2,3)) 
D[] <- lapply(D, factor)