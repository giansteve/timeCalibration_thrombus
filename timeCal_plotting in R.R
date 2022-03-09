# setting job ---------------------
# Setting the environment
library(factoextra)
library(readxl)
library(devtools)
library(ggbiplot)
library(graphics)
library(corrplot)
library(ggplot2)

# set directory
setwd("M:/IFM/User/melito/Server/Projects/TimeCalibration")

# load data ---------------------------------------------------------------------------
mat = read_excel("MRIdata.xlsx", range = "B20:J26", na="")

# plot dtata -----------------
ggplot(mat, aes(x = mat[1], y = mat[6])) +
  geom_line() +
  geom_ribbon(aes(ymin = mat[6] - mat[7],
                  ymax = mat[6] + mat[7]), alpha = 0.2)