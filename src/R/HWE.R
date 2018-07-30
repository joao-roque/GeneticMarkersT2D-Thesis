library(deFinetti)

data <- read.csv(file = "../../data/plots/hwe/frequencies.csv")
data <- data[-1]

deFinetti.plot(data, with_F_color = FALSE, markerpos = 2)

