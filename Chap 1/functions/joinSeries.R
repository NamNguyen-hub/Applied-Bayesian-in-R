if(!"devtools" %in% rownames(installed.packages())) install.packages("devtools")

# Install the stable development verions from GitHub
devtools::install_github("crsh/papaja")
library("papaja")
library(ggplot2)


#Pseudo Codes
##Load series
setwd("D:/GitHub/Applied-Bayesian-in-R/Chap 1")

Credit_filepath = sprintf("Credit_HPfilter_%s.txt",country)
df2 <- read.table(Credit_filepath, header=TRUE, sep=",")
df2 <- na.omit(df2[-c(2)]) #Remove country name column because redundancy


HP_filepath = sprintf("HPindex_HPfilter_%s.txt",country)

df1 <- read.table(HP_filepath, header=TRUE, sep=",")

### observed series
Y
### Forecast Series
EB1
### Rearrange and combine
T_combined = nrow(Y)+12
Y_combined=matrix(NA,T_combined,10)
Y_combined[1:nrow(Y),1]=Y
Y_combined[(nrow(Y)+1):(nrow(Y)+12),2:10]=t(EB1)
##Convert to time series
##Set date and variable

# Example codes: https://www.statmethods.net/advstats/timeseries.html
# # from Jan 2009 to Dec 2014 as a time series object
# myts <- ts(myvector, start=c(2009, 1), end=c(2014, 12), frequency=12)
# 
# # subset the time series (June 2014 to December 2014)
# myts2 <- window(myts, start=c(2014, 6), end=c(2014, 12))
# 
# # plot series
# plot(myts)

dates = seq(as.Date("1948-01-01"), by = "quarter", length.out = 262)
dates = as.character(dates)

Y_combined = cbind(Y_combined,dates)
Y_combined = as.data.frame(Y_combined)
myts <- ts(Y_combined, start=c(1948, 1), frequency=12)

#Plot
##Load variables

# Ex1
# #Graph raw data for each country
# ggplot(rates_plot, aes(date, obs_value, color = borrowers_country)) +
#   geom_hline(yintercept = 0, linetype = "dashed",
#              color = "grey70", size = 0.02) +
#   geom_line(show.legend = FALSE) +
#   facet_wrap(~borrowers_country) +
#   theme_light() +
#   theme(panel.grid = element_blank()) +
#   labs(x = NULL, y = NULL,
#        title = "Credit to household",
#        subtitle = "as percentage of GDP")

# Ex2
# library(crsh/papaja)
# 
# ## Not run: 
# # Copied from ?ggtheme
library("reshape2")
library("ggplot2")

for (i in 1:10){
  Y_combined[,i]=as.numeric(Y_combined[,i])
}

Y_combined[,11]=as.Date(Y_combined[,11])

Y_combined=subset(Y_combined, dates>as.Date("2000-01-01"))
test_data_long <- melt(Y_combined, id="dates")  # convert to long format

p <- ggplot(data=test_data_long,
       aes(x=dates, y=value, colour=variable)) +
  geom_line()
p
p + theme_apa()


## End(Not run)


