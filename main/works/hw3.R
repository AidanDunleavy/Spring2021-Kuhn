
# random assignment for the block
library(agricolae)
treat <- c("flav1-reg", "flav2-reg", "flav3-reg", "flav1-sug",
           "flav2-sug", "flav3-sug") 
# seed 11 gives a single randomization for ROC data
outdesign <- design.rcbd(treat, 3, seed = 11); rcb 
rcb <- outdesign$book; rcb
# ROC of mice for different temps, block over mice
levels(rcb$block) <- c("subject1", "subject2", "subject3"); rcb




D <- expand.grid( flavor = c(1, 2, 3),
                  sugarlevel = c("reg", "sug"),
                  block = c(1,2,3)) 
D[] <- lapply(D, factor)

D
y = c(10.3, 18.8, 11.2, 12.4, 16.1, 15.1, 9.1, 12.1, 6.5, 16.1, 7.2, 9.8, 4.1, 8.1, 3.5, 6.1, 2.2, 5.8)
gum <- data.frame( D, y ); gum
set.seed(2591) # ensures everyone has same randomization
D <- D[order(sample(1:18)), ]; D

modelgum = lm(y~ block + flavor*sugarlevel, data = gum)

anova(modelgum)





D <- expand.grid(barley=c("manchuria","svansota","velvet","trebi","peatland"),
                 year=c(1931,1932),
                 place=c(1,2,3,4,5,6))

D[] = lapply(D, factor)

growthlevels=c(81,105.4,119.7,109.7,98.3,80.7,82.3,80.4,87.2,84.2,146.6,142,150.7,
               191.5,145.7,100.4,115.5,112.2,147.7,108.1,82.3,77.3,78.4,131.3,89.6,
               103.1,105.1,116.5,139.9,129.6,119.8,121.4,124,140.8,124.8,98.9,61.9,
               96.2,125.5,75.7,98.9,89,69.1,89.3,104.1,66.4,49.9,96.7,61.9,80.3,86.9,
               77.1,78.9,101.8,96,67.7,66.7,67.4,91.8,94.1)

GLdata <- data.frame(cbind(D,growthlevels)); GLdata
m.gl <- aov(growthlevels ~ place + year + barley, data=GLdata) 
summary(m.gl)




m.gl <- aov(growthlevels ~ place*year + barley, data=GLdata)
summary(m.gl)




# D <- expand.grid(barley=c("manchuria","svansota","velvet","trebi","peatland"),year_place=c("1931-1","1932-1","1931-2","1932-2","1931-3","1932-3","1931-4","1932-4","1931-5","1932-5","1931-6","1932-6"))
# 
# D[] = lapply(D, factor)
# 
# growthlevels=c(81,105.4,119.7,109.7,98.3,80.7,82.3,80.4,87.2,84.2,146.6,142,150.7,191.5,145.7,100.4,115.5,112.2,147.7,108.1,82.3,77.3,78.4,131.3,89.6,103.1,105.1,116.5,139.9,129.6,119.8,121.4,124,140.8,124.8,98.9,61.9,96.2,125.5,75.7,98.9,89,69.1,89.3,104.1,66.4,49.9,96.7,61.9,80.3,86.9,77.1,78.9,101.8,96,67.7,66.7,67.4,91.8,94.1)
# 
# GLdata <- data.frame(cbind(D,growthlevels)); GLdata
# barley year_place growthlevels
# m.gl <- aov(growthlevels ~ barley + year.place, data=GLdata) 
# summary(m.gl) # RCB design





D <- expand.grid(barley=c("manchuria","svansota","velvet","trebi","peatland"),year.place=c("1931-1","1932-1","1931-2","1932-2","1931-3","1932-3","1931-4","1932-4","1931-5","1932-5","1931-6","1932-6"))

D[] = lapply(D, factor)

growthlevels=c(81,105.4,119.7,109.7,98.3,80.7,82.3,80.4,87.2,84.2,146.6,142,150.7,191.5,145.7,100.4,115.5,112.2,147.7,108.1,82.3,77.3,78.4,131.3,89.6,103.1,105.1,116.5,139.9,129.6,119.8,121.4,124,140.8,124.8,98.9,61.9,96.2,125.5,75.7,98.9,89,69.1,89.3,104.1,66.4,49.9,96.7,61.9,80.3,86.9,77.1,78.9,101.8,96,67.7,66.7,67.4,91.8,94.1)

GLdata <- data.frame(cbind(D,growthlevels)); GLdata
m.gl <- aov(growthlevels ~ barley + year.place, data=GLdata) 
summary(m.gl) # RCB design

model.tables( m.gl , type = "means" )
# tukey pairwise comparison of ROC means
TukeyHSD( m.gl, "barley")
