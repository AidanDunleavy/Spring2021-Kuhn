#' ---
#' title: Test 1
#' author: Aidan Dunleavy
#' date: 01/21/2021
#' ---
head(dough)
data <- dough
names(data) <- c("rt", "Group")
data$Group <- ordered(data$Group,
                      levels = c("0.25", "0.5", "0.75", "1"))
plot(data)
data
res.aov <- aov(rt ~ Group, data = data)
summary(res.aov)
plot(res.aov)

