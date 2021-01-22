#' ---
#' title: HW 1
#' author: Aidan Dunleavy
#' date: 01/20/2021
#' ---
#' PROBLEM 1
data = read.csv("https://raw.githubusercontent.com/AidanDunleavy/Fall-2020-STAT-40001/master/paper.csv")
data$Group <- ordered(data$Group,
                         levels = c("4", "4.75", "5.5", "6"))
names(data) <- c("time", "Group")

levels(data$Group)
library("ggpubr")
plot.new()
ggboxplot(data, x = "Group", y = "time",
          color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#00bb57"),
          order = c("4", "4.75", "5.5", "6"),
          ylab = "time", xlab = "Treatment")

plot.new()
boxplot(time ~ Group, data = data,
        xlab = "Treatment", ylab = "time",
        frame = FALSE, col = c("#00AFBB", "#E7B800", "#FC4E07", "#00bb57"))

res.aov <- aov(time ~ Group, data = data)
summary(res.aov)

plot.new()
par(mfrow=c(2,2))
plot(res.aov)
par(mfrow = c(1,1))
data$Group
plot.new()
plot(c("4", "4.75", "5.5", "6"), by(data$time, data$Group, mean),
     col = c("#00AFBB", "#E7B800", "#FC4E07", "#00bb57"))

plot(data$Group, data$time)
model = lm(time ~ Group, data = data)
summary(model)

