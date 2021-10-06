# Pseudocodes
# Limit to first 500 iterations 
  ## Set reps in ex1_3 to 500
  ## Set burns in ex1_3 to 0
## plot raw Gibbs sample output
par(mfrow=c(3,2))
df = as.data.frame(out2)
p1 <- ggplot(data= df,
            aes(x=as.numeric(row.names(df)), y=V1)) +
  geom_line() +
  labs(x = "iterations",
       title = "intercept")

p2 <- ggplot(data= df,
            aes(x=as.numeric(row.names(df)), y=V2)) +
  geom_line() +
  labs(x = "iterations",
       title = "B_1")

p3 <- ggplot(data= df,
            aes(x=as.numeric(row.names(df)), y=V3)) +
  geom_line() +
  labs(x = "iterations",
       title = "B_2")

df = as.data.frame(out3)
p4 <- ggplot(data= df,
            aes(x=as.numeric(row.names(df)), y=out3)) +
  geom_line() +
  labs(x = "iterations",
       title = "rho")

df = as.data.frame(out4)
p5 <- ggplot(data= df,
            aes(x=as.numeric(row.names(df)), y=out4)) +
  geom_line() +
  labs(x = "iterations",
       title = "sigma^2")


library(patchwork)
(p1|p3)/(p2|p4)/(p5|plot_spacer())


## plot recursive means of Gibbs sample output
df = as.data.frame(out2)
df1=as.data.frame(colMeans(matrix(df[,1], 10)))
df2=as.data.frame(colMeans(matrix(df[,2], 10)))
df3=as.data.frame(colMeans(matrix(df[,3], 10)))
df=cbind(df1,df2,df3)
names(df)=c("intercept","B_1","B_2")

p1 <- ggplot(data= df,
             aes(x=as.numeric(row.names(df1)), y=intercept)) +
  geom_line() +
  labs(x = "iterations",
       title = "intercept")
p1


df1=as.data.frame(colMeans(matrix(df[,2], 10)))
p2 <- ggplot(data= df,
             aes(x=as.numeric(row.names(df1)), y=B_1)) +
  geom_line() +
  labs(x = "iterations",
       title = "B_1")
p2


df1=as.data.frame(colMeans(matrix(df[,3], 10)))
p3 <- ggplot(data= df,
             aes(x=as.numeric(row.names(df1)), y=B_2)) +
  geom_line() +
  labs(x = "iterations",
       title = "B_2")

df = as.data.frame(out3)
df1=as.data.frame(colMeans(matrix(df[,1], 10)))
names(df1)[1]="out3"
p4 <- ggplot(data= df1,
             aes(x=as.numeric(row.names(df1)), y=out3)) +
  geom_line() +
  labs(x = "iterations",
       title = "rho")
p4

df = as.data.frame(out4)
df1=as.data.frame(colMeans(matrix(df[,1], 10)))
names(df1)[1]="out4"
p5 <- ggplot(data= df1,
             aes(x=as.numeric(row.names(df1)), y=out4)) +
  geom_line() +
  labs(x = "iterations",
       title = "sigma^2")
p5

(p1|p3)/(p2|p4)/(p5|plot_spacer())

## plot acf
par(mfrow=c(2,3))

df = as.data.frame(out2)
acf(df$V1)
acf(df$V2)
acf(df$V3)
df = as.data.frame(out3)
acf(df$out3)
df = as.data.frame(out4)
acf(df$out4)
## plot raw Gibbs sample output

# Do again but for reps = 25000 iterations
