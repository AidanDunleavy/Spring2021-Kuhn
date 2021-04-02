x1 <- seq(-1, 1, by=0.01); x1
x2 <- seq(-1, 1, by=0.01); x2
X1 <- c(rep(0,40401)); X1
X2 <- c(rep(0,40401)); X2
c <- 0
for (i in 1:201) {
  for (j in 1:201) {
    if ((x2[j] >= -2*x1[i] - 2) & (x2[j] <= -2*x1[i] + 1)) {
      c <- c + 1
      X1[c] <- x1[i]  
      X2[c] <- x2[j]
    }
  }
}
grid <- data.frame(cbind(X1,X2))
grid.constr <- grid[-(c+1):-nrow(grid),]; grid.constr # candidate points

# creating a D-optimal design
# for irregular experimental region
library(AlgDesign)
d <- data.frame(var=paste("x",1:2,sep=""),  
                low=-1,
                high=1,
                center=0,
                nLevels=21,
                round=1,
                factor=FALSE); d
constFcn <- function(x){x[2] >= -2*x[1] - 2 && x[2] <= -2*x[1] + 1}
desCon <- optMonteCarlo(~quad(.), 
                        d, 
                        constraint=constFcn,
                        criterion = "D"); desCon

# finer grid for irregular experimental region
library(AlgDesign)
d <- data.frame(var=paste("x",1:2,sep=""),  
                low=-1,
                high=1,
                center=0,
                nLevels=201,
                round=2,
                factor=FALSE); d
constFcn <- function(x){x[2] >= -2*x[1] - 2 && x[2] <= -2*x[1] + 1}
desCon <- optMonteCarlo(~quad(.), 
                        d, 
                        constraint=constFcn,
                        criterion = "D"); desCon






A <- factor(rep(c(1,1,2,2,3,3), each = 4)); A # 36 runs for whole factor A
B <- factor(rep(c(1,2,3,1,2,3), each = 4)); B # whole factor B
wp <- factor(rep(c(1,2,3,4,5,6,7,8,9), each = 4)); wp
C <- factor(rep(rep(c(1,2), each = 2)), times = 9); C # subplot factor C
D <- factor(rep(c(1,2), times = 16)); D # subplot factor D
sp <- factor(rep(c(1,2), times = 16)); sp
splitP <- data.frame(cbind(wp,sp,A,B,C,D)); splitP # split plot design
