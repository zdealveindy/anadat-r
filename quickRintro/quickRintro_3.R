# The plotting figures in R

cars

cars$speed
cars$dist

plot (x = cars$speed, y = cars$dist)

plot (dist ~ speed, data = cars)  # formula interface: y ~ x, data

plot (dist ~ speed, data = cars, main = 'Dependence of distance on speed',
      xlab = 'Speed (mph)', ylab = 'Stopping distance (ft)',
      xlim = c(0, 30))

hist (cars$dist, col = 'red')

plot (dist ~ speed, data = cars, type = 'n', axes = FALSE, ann = FALSE)
points (dist ~ speed, data = cars, col = 'red')
# text (dist ~ speed, data = cars)
axis (1)
axis (2)
box ()
title (main = 'Dependence of distance on speed', xlab = 'Speed (mph)', ylab = 'Distance (ft)')

#Orange
jpeg (filename = 'Orange_scatterplot.jpeg', 
      width = 8, height = 6, units = 'in', res = 300)
plot (circumference ~ age, data = Orange, col = Orange$Tree)
dev.off ()
