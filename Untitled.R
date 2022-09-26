setwd("/Users/mollycui/Desktop/R script/Research 10-Epid PMCMC")

v1 <- read.csv("andre_estimates_21_02.txt", sep  = "\t") %>%
  rowSums()
y1 <- data.frame(value = v1) %>%
 mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)




v2 <- read.csv("covid.csv", header=FALSE, stringsAsFactors=FALSE) %>%
  rowSums()
  v2<-v2[1:37]
y2 <- data.frame(value = v2) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)
#y2<- y2[1:37,]



plot(y1,type='l',col='red')
plot(y2,type='o')
lines(y1,type='o',col='red')

```{r load_data, echo = T}
library(rbi)
library(rbi.helpers)
# Load the data
v <- read.csv("covid.csv", header=FALSE, stringsAsFactors=FALSE) %>%
  rowSums()
y <- data.frame(value = v) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)
y<- y[1:30,]
ncores <- 8
minParticles <- max(ncores, 16)
```