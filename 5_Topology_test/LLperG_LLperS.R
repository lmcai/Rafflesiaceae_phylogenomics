y=read.table('LLperS.csv',header = T)
x=read.table('LLperG.csv',header = T)
s<-x[order(x$Rate), ]
plot(x$tree_median_root2tip,x$Tree1-x$Tree2,pch=20,cex=0.3,xlim=c(0,2),ylim=c(-50,25))
abline(h=0)
plot(x$tree_median_root2tip,x$Tree1-x$Tree3,pch=20,cex=0.3,xlim=c(0,2),ylim=c(-50,25))
abline(h=0)
plot(x$tree_median_root2tip,x$Tree1-x$Tree4,pch=20,cex=0.3,xlim=c(0,2),ylim=c(-50,25))
abline(h=0)
plot(x$tree_median_root2tip,x$Tree1-x$Tree5,pch=20,cex=0.3,xlim=c(0,2),ylim=c(-50,25))
abline(h=0)
plot(x$tree_median_root2tip,x$Tree1-x$Tree6,pch=20,cex=0.3,xlim=c(0,2),ylim=c(-50,25))
abline(h=0)
plot(x$tree_median_root2tip,x$Tree1-x$Tree7,pch=20,cex=0.3,xlim=c(0,2),ylim=c(-50,25))
abline(h=0)
plot(x$tree_median_root2tip,x$Tree1-x$Tree8,pch=20,cex=0.3,xlim=c(0,2),ylim=c(-50,25))
abline(h=0)

y$row_num=1:3099966
y$class <- as.factor(y$row_num %% 3 + 1)
a=x[y$class==1,]
b=x[y$class==2,]
c=x[y$class==3,]

##########################################
# Define the number of reorders
num_reorders <- 1000

# Initialize a matrix to store cumulative values
cumulative_values <- matrix(0, nrow = 2134, ncol = num_reorders)

# Perform reordering and compute cumulative sums
for (i in 1:num_reorders) {
  shuffled_dataset <- x[sample(1:nrow(x)), ]
  cumulative_values[, i] <- cumsum(shuffled_dataset$Tree3 - shuffled_dataset$Tree5)
}

# Calculate the averaged curve
average_curve <- rowMeans(cumulative_values)

# Plot the curves
plot(1:2134, cumulative_values[, 1], ylim=c(-200,500),type = "l", col = rgb(0.8, 0.8, 0.8, 0.1),
     xlab = "Number of genes", ylab = "Cumulative delta LL")
for (i in 2:num_reorders) {
  lines(1:2134, cumulative_values[, i], col = rgb(0.8, 0.8, 0.8, 0.1))
}
lines(1:2134, average_curve, col = "blue", lwd = 2)  # Add averaged curve in blue

h=0
for (i in 1:2134){
	h=h+s$Tree3[i]-s$Tree5[i]
	s$temp[i]=h
}
lines(s$temp,type='l',col='red')


sum(y$Tree1[y$Site %% 3 ==0])
[1] -19256079
> sum(y$Tree2[y$Site %% 3 ==0])
[1] -19256084
> sum(y$Tree3[y$Site %% 3 ==0])
[1] -19256139
> sum(y$Tree4[y$Site %% 3 ==0])
[1] -19256252
> sum(y$Tree5[y$Site %% 3 ==0])
[1] -19255965
> sum(y$Tree6[y$Site %% 3 ==0])
[1] -19256128
> sum(y$Tree7[y$Site %% 3 ==0])
[1] -19256323
> sum(y$Tree8[y$Site %% 3 ==0])
[1] -19256358

sum(y$Tree1[y$Site %% 3 ==1])
[1] -9762863
> sum(y$Tree2[y$Site %% 3 ==1])
[1] -9762903
> sum(y$Tree3[y$Site %% 3 ==1])
[1] -9762884
> sum(y$Tree4[y$Site %% 3 ==1])
[1] -9763647
> sum(y$Tree5[y$Site %% 3 ==1])
[1] -9763129
> sum(y$Tree6[y$Site %% 3 ==1])
[1] -9763031
> sum(y$Tree7[y$Site %% 3 ==1])
[1] -9762969
> sum(y$Tree8[y$Site %% 3 ==1])
[1] -9762945

sum(y$Tree1[y$Site %% 3 ==2])
[1] -7779921
> sum(y$Tree2[y$Site %% 3 ==2])
[1] -7779918
> sum(y$Tree3[y$Site %% 3 ==2])
[1] -7779897
> sum(y$Tree4[y$Site %% 3 ==2])
[1] -7780473
> sum(y$Tree5[y$Site %% 3 ==2])
[1] -7780020
> sum(y$Tree6[y$Site %% 3 ==2])
[1] -7780029
> sum(y$Tree7[y$Site %% 3 ==2])
[1] -7779916
> sum(y$Tree8[y$Site %% 3 ==2])
[1] -7779900