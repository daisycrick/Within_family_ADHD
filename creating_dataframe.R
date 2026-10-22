rm(list = ls())
install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")
library(dplyr)
library(purrr)
library(tidyr)
install.packages("readr")
library(readr)

setwd ("~/OneDrive - University College London/Documents/ADHD_brisbane_paper/shrinkage/update/dat_s0.4_se0.5/MR/") # Set the file path to your local directory


a<-read_csv("ieu-b-4763_presteiger.csv")
b<-read_csv("ieu-b-4815_presteiger.csv")
c<-read_csv("ieu-b-4819_presteiger.csv")
d<-read_csv("ieu-b-4821_presteiger.csv")
f<-read_csv("ieu-b-4823_presteiger.csv")
g<-read_csv("ieu-b-4825_presteiger.csv")
h<-read_csv("ieu-b-4829_presteiger.csv")
i<-read_csv("ieu-b-4833_presteiger.csv")
j<-read_csv("ieu-b-4835_presteiger.csv")
k<-read_csv("ieu-b-4837_presteiger.csv")
l<-read_csv("ieu-b-4839_presteiger.csv")
m<-read_csv("ieu-b-4841_presteiger.csv")
n<-read_csv("ieu-b-4843_presteiger.csv")
o<-read_csv("ieu-b-4845_presteiger.csv")
p<-read_csv("ieu-b-4847_presteiger.csv")
q<-read_csv("ieu-b-4849_presteiger.csv")

a<-subset(a, select=c(3,6,7,8,9,10,11))
b<-subset(b, select=c(3,6,7,8,9,10,11))
c<-subset(c, select=c(3,6,7,8,9,10,11))
d<-subset(d, select=c(3,6,7,8,9,10,11))
f<-subset(f, select=c(3,6,7,8,9,10,11))
g<-subset(g, select=c(3,6,7,8,9,10,11))
h<-subset(h, select=c(3,6,7,8,9,10,11))
i<-subset(i, select=c(3,6,7,8,9,10,11))
j<-subset(j, select=c(3,6,7,8,9,10,11))
k<-subset(k, select=c(3,6,7,8,9,10,11))
l<-subset(l, select=c(3,6,7,8,9,10,11))
m<-subset(m, select=c(3,6,7,8,9,10,11))
n<-subset(n, select=c(3,6,7,8,9,10,11))
o<-subset(o, select=c(3,6,7,8,9,10,11))
p<-subset(p, select=c(3,6,7,8,9,10,11))
q<-subset(q, select=c(3,6,7,8,9,10,11))


dat<-rbind(a,b,c,d,f,g,h,i,j,k,l,m,n,o,p,q)
head(dat)

dat$s<-0.4

#dat2 <- dat[, c(1, 2, 8, 9, 10)]

write.csv(dat, "dat_s0.4_se0.5_combined.csv", row.names = FALSE)


a<-read_csv("~/OneDrive - University College London/Documents/ADHD_brisbane_paper/shrinkage/update/dat_s0.05_se0.5/MR/dat_s0.05_se0.5_combined.csv")
b<-read_csv("~/OneDrive - University College London/Documents/ADHD_brisbane_paper/shrinkage/update/dat_s0.1_se0.5/MR/dat_s0.1_se0.5_combined.csv")
c<-read_csv("~/OneDrive - University College London/Documents/ADHD_brisbane_paper/shrinkage/update/dat_s0.2_se0.5/MR/dat_s0.2_se0.5_combined.csv")
d<-read_csv("~/OneDrive - University College London/Documents/ADHD_brisbane_paper/shrinkage/update/dat_s0.3_se0.5/MR/dat_s0.3_se0.5_combined.csv")
e<-read_csv("~/OneDrive - University College London/Documents/ADHD_brisbane_paper/shrinkage/update/dat_s0.4_se0.5/MR/dat_s0.4_se0.5_combined.csv")


combined<-rbind(a,b,c,d,e)
write.csv(combined, "MR_shirnkage_no_delta.csv", row.names = FALSE)
