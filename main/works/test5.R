#9.2
square <- factor(c(
  "I", "I", "I", "I", "I", "I", "I", "I", 
  "I", "I", "I", "I", "I", "I", "I", "I",
  "I", "I", "I", "I", "I", "I", "I", "I",
  "II", "II", "II", "II", "II", "II", "II", "II",
  "II", "II", "II", "II", "II", "II", "II", "II",
  "II", "II", "II", "II", "II", "II", "II", "II"
)); square
group <- factor(c(
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3,
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3
)); group
cow <- factor(c(
  1, 1, 1, 1, 2, 2, 2, 2,
  3, 3, 3, 3, 4, 4, 4, 4,
  5, 5, 5, 5, 6, 6, 6, 6,
  7, 7, 7, 7, 8, 8, 8, 8,
  9, 9, 9, 9, 10, 10, 10, 10,
  11, 11, 11, 11, 12, 12, 12, 12
)); cow
period <- factor(c(
  1, 2, 3, 4, 1, 2, 3, 4,
  1, 2, 3, 4, 1, 2, 3, 4,
  1, 2, 3, 4, 1, 2, 3, 4,
  1, 2, 3, 4, 1, 2, 3, 4,
  1, 2, 3, 4, 1, 2, 3, 4,
  1, 2, 3, 4, 1, 2, 3, 4
)); period
diet <- factor(c(
  "A", "B", "C", "C", "A", "B", "C", "C",
  "B", "C", "A", "A", "B", "C", "A", "A",
  "C", "A", "B", "B", "C", "A", "B", "B",
  "A", "C", "B", "B", "A", "C", "B", "B",
  "B", "A", "C", "C", "B", "A", "C", "C",
  "C", "B", "A", "A", "C", "B", "A", "A"
)); diet
carry <- factor(c(
  "0", "A", "B", "C", "0", "A", "B", "C", 
  "0", "B", "C", "A", "0", "B", "C", "A",
  "0", "C", "A", "B", "0", "C", "A", "B",
  "0", "A", "C", "B", "0", "A", "C", "B",
  "0", "B", "A", "C", "0", "B", "A", "C",
  "0", "C", "B", "A", "A", "0", "C", "B"
)); carry
fcm <- c(
  38.66, 37.43, 34.39, 31.30, 25.72, 26.13, 23.35, 18.69,
  48.85, 46.88, 41.99, 39.61, 30.80, 29.29, 26.41, 23.16,
  34.64, 32.27, 28.50, 27.13, 25.35, 26.00, 23.86, 19.92,
  35.19, 33.50, 28.41, 25.12, 21.80, 23.91, 21.69, 17.55,
  32.90, 33.12, 27.52, 25.10, 21.37, 21.97, 19.38, 16.57,
  30.40, 29.50, 26.70, 23.09, 22.84, 20.97, 18.59, 16.10
); fcm
data <- data.frame(square, group, cow, period, diet, carry, fcm); data
m <- lm(fcm ~ cow + period + diet + carry, 
        data = data, 
        contrasts = list(cow = contr.sum, period = contr.sum, 
                         diet = contr.sum, carry = contr.sum))
library(car)
summary.aov(m, type = "I", singular.ok = TRUE)
summary.aov(m, type = "III", singular.ok = TRUE)

library(lsmeans)
lsm(m)


with(data, tapply(fcm, diet, mean))
sqrt(with(data, tapply(fcm, diet, var)))
with(data, tapply(fcm, carry, mean))
sqrt(with(data, tapply(fcm, carry, var)))
m <- mean(data$fcm); m # grand mean roc


aggregate(data$fcm, by=list(Group=data$diet), FUN = sum)
