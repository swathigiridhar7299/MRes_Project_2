#######################################################################
## Necessary Libraries
library(npreg)
library(tidyverse)
#######################################################################
df <- read.csv("file_name.csv")
##Choosing the treatment
data <- df[df$Treatment == "Control", ]
glimpse(data)
data[apply(data, 1, function(x) all(is.finite(x))), ]
# Remove the leading "D" from each value in the Day column
data$Day <- as.numeric(sub("D", "", data$Day))
# Verify the updated Day column
print(data$Day)
#############################################################################
# Construct the smoothing splines and compare with a linear fit 
## Fit the ss
mod.ss <- ss(data$Day,data$Early.Granulocytes.Count, method = "GCV")
mod.ss
summary(mod.ss)
## Plot the ss fit
plot(mod.ss,xlab='Day',ylab='Normalized Early Granulocytes Cell Count',col = "red")
# Fit linear regression model
lm_model <- lm(Early.Granulocytes.Count ~ Day, data = data)
abline(lm_model, lty = 2,col="blue")
points(data$Day, data$Early.Granulocytes.Count)
# Add a legend
legend("topleft", legend = c("ss", "lm"), lty = c(1, 2), col = c("red", "blue"), lwd = 2, bty = "n")
#############################################################################
## Note - looped through all treatment types and doses