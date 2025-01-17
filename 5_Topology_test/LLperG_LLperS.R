y=read.csv('LLperS.csv',header = T)
x=read.csv('LLperG.csv',header = T)
s<-x[order(x$Rate), ]

#################
pdf("LLperG.pdf", width = 8.5, height = 11)

# Calculate the plotting region to fit within 1-inch margins
# The inner area is 6.5 inches wide and 9 inches high after leaving 1-inch margins
layout(matrix(1:8, nrow = 4, byrow = TRUE)) # Arrange plots in a 4x2 grid
par(omi = c(1, 1, 1, 1), # Outer margins: 1 inch on all sides
    mai = c(0.5, 0.5, 0.4, 0.2)) # Inner margins for individual panels

# Generate the plots
plot(x$Tree1 - x$Tree2, type = 'l', main = "Tree1 - Tree2")
plot(x$Tree1 - x$Tree3, type = 'l', main = "Tree1 - Tree3")
plot(x$Tree1 - x$Tree4, type = 'l', main = "Tree1 - Tree4")
plot(x$Tree1 - x$Tree5, type = 'l', main = "Tree1 - Tree5")
plot(x$Tree1 - x$Tree6, type = 'l', main = "Tree1 - Tree6")
plot(x$Tree1 - x$Tree7, type = 'l', main = "Tree1 - Tree7")
plot(x$Tree1 - x$Tree8, type = 'l', main = "Tree1 - Tree8")
# If there are fewer than 8 plots, add blank space
plot.new()  # Adds an empty plot to fill the layout if necessary
# Close the PDF device
dev.off() # Ensures all plots are written to the file


y$row_num=1:3099966
y$class <- as.factor(y$row_num %% 3 + 1)
a=x[y$class==1,]
b=x[y$class==2,]
c=x[y$class==3,]

##########################################
# Support per gene examination
# Define the number of reorders
num_reorders <- 1000

# Initialize a matrix to store cumulative values
cumulative_values <- matrix(0, nrow = 2135, ncol = num_reorders)

# Perform reordering and compute cumulative sums
for (i in 1:num_reorders) {
  shuffled_dataset <- x[sample(1:nrow(x)), ]
  cumulative_values[, i] <- cumsum(shuffled_dataset$Tree1 - shuffled_dataset$Tree4)
}

# Calculate the averaged curve
average_curve <- rowMeans(cumulative_values)

# Calculate 95% confidence interval
lower_bound <- apply(cumulative_values, 1, quantile, probs = 0.05)
upper_bound <- apply(cumulative_values, 1, quantile, probs = 0.95)

# Plot the curves
plot(1:2135, cumulative_values[, 1], ylim=c(-100,1500),type = "l", col = rgb(0.1, 0.1, 0.1, 0.05),
     xlab = "Number of genes", ylab = "Cumulative delta LL")
for (i in 2:num_reorders) {
  lines(1:2135, cumulative_values[, i], col = rgb(0.1, 0.1, 0.1, 0.05))
}

lines(1:2135, average_curve, col = "yellow", lwd = 2)  # Add averaged curve in blue
# Add 95% confidence interval band
polygon(c(1:2135, 2135:1), c(lower_bound, rev(upper_bound)), 
        col = rgb(0.2, 0, 0.7, 0.3), border = NA)

h=0
for (i in 1:2135){
	h=h+s$Tree1[i]-s$Tree4[i]
	s$temp[i]=h
}
lines(s$temp,type='l',col='red')

###########################################
#Test if genes grouping Raff and Apo is evolving faster

##############
#DNA data
z=read.csv('monophyly_raff_apo.csv')
merged_df <- merge(x, z, by = "Gene_ID")
t_test_result <- t.test(merged_df$Rate[merged_df$Monophyly_raff_apo=='N'], merged_df$Rate[merged_df$Monophyly_raff_apo=='Y'], alternative = "greater")
print(t_test_result)

	Welch Two Sample t-test

data:  merged_df$Rate[merged_df$Monophyly_raff_apo == "N"] and merged_df$Rate[merged_df$Monophyly_raff_apo == "Y"]
t = -0.81793, df = 839.49, p-value = 0.7932
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 -0.03098754         Inf
sample estimates:
mean of x mean of y 
0.9996147 1.0098986 

#the other direction:
t_test_result <- t.test(merged_df$Rate[merged_df$Monophyly_raff_apo=='Y'], merged_df$Rate[merged_df$Monophyly_raff_apo=='N'], alternative = "greater")
print(t_test_result)

	Welch Two Sample t-test

data:  merged_df$Rate[merged_df$Monophyly_raff_apo == "Y"] and merged_df$Rate[merged_df$Monophyly_raff_apo == "N"]
t = 0.81793, df = 839.49, p-value = 0.2068
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 -0.01041975         Inf
sample estimates:
mean of x mean of y 
1.0098986 0.9996147 

########################
###Protein data
a=read.csv('aa.LLperG.csv')
merged_df <- merge(a, z, by = "Gene_ID")
print(t_test_result)

	Welch Two Sample t-test

data:  merged_df$Rate[merged_df$Monophyly_raff_apo == "Y"] and merged_df$Rate[merged_df$Monophyly_raff_apo == "N"]
t = 0.45399, df = 895.49, p-value = 0.325
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 -0.03037397         Inf
sample estimates:
mean of x mean of y 
0.9948065 0.9832435 

#the other direction
t_test_result <- t.test(merged_df$Rate[merged_df$Monophyly_raff_apo=='N'], merged_df$Rate[merged_df$Monophyly_raff_apo=='Y'], alternative = "greater")
print(t_test_result)

	Welch Two Sample t-test

data:  merged_df$Rate[merged_df$Monophyly_raff_apo == "N"] and merged_df$Rate[merged_df$Monophyly_raff_apo == "Y"]
t = -0.45399, df = 895.49, p-value = 0.675
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 -0.05349986         Inf
sample estimates:
mean of x mean of y 
0.9832435 0.9948065 
#############################################

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