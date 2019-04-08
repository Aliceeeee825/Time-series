library("quantmod")
sp500<-new.env()
dj30 <- c("MMM","AXP","AAPL","BA","CAT","CVX","CSCO","KO","DIS","DWDP",
          "XOM","GS","HD","IBM","INTC","JNJ","JPM","MCD","MRK","MSFT",
          "NKE","PFE","PG","TRV","UTX","UNH","VZ","V","WMT","WBA")


#download data for all 30 stocks from the webstie. Change the data to monthly data using apply monthly and extract the data out by using coredata. 
#all the NA from the data are removed
# monthly return is also calculated by getting the difference between each month after taking log of each data

ne <- new.env()
getSymbols('MMM', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_MMM <- apply.monthly(ne$MMM$MMM.Adjusted, last)
data1=coredata(diff(log(data_MMM)))
data1 = na.omit(data1)

getSymbols('AXP', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_AXP <- apply.monthly(ne$AXP$AXP.Adjusted, last)
data2=coredata(diff(log(data_AXP)))
data2 = na.omit(data2)

getSymbols('AAPL', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_AAPL <- apply.monthly(ne$AAPL$AAPL.Adjusted, last)
data3 = coredata(diff(log(data_AAPL)))
data3 = na.omit(data3)

getSymbols('BA', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_BA <- apply.monthly(ne$BA$BA.Adjusted, last)
data4 = coredata(diff(log(data_BA)))
data4 = na.omit(data4)

getSymbols('CAT', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_CAT <- apply.monthly(ne$CAT$CAT.Adjusted, last)
data5 = coredata(diff(log(data_CAT)))
data5 = na.omit(data5)

getSymbols('CVX', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_CVX <- apply.monthly(ne$CVX$CVX.Adjusted, last)
data6 = coredata(diff(log(data_CVX)))
data6 = na.omit(data6)

getSymbols('CSCO', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_CSCO <- apply.monthly(ne$CSCO$CSCO.Adjusted, last)
data7 = coredata(diff(log(data_CSCO)))
data7= na.omit(data7)

getSymbols('KO', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_KO <- apply.monthly(ne$KO$KO.Adjusted, last)
data8 = coredata(diff(log(data_KO)))
data8 = na.omit(data8)

getSymbols('DIS', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_DIS <- apply.monthly(ne$DIS$DIS.Adjusted, last)
data9 = coredata(diff(log(data_DIS)))
data9 = na.omit(data9)

getSymbols('DWDP', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_DWDP <- apply.monthly(ne$DWDP$DWDP.Adjusted, last)
data10 = coredata(diff(log(data_DWDP)))
data10 = na.omit(data10)

getSymbols('XOM', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_XOM <- apply.monthly(ne$XOM$XOM.Adjusted, last)
data11= coredata(diff(log(data_XOM)))
data11 = na.omit(data11)

getSymbols('GS', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_GS <- apply.monthly(ne$GS$GS.Adjusted, last)
data12 = coredata(diff(log(data_GS)))
data12 = na.omit(data12)

getSymbols('HD', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_HD <- apply.monthly(ne$HD$HD.Adjusted, last)
data13 = coredata(diff(log(data_HD)))
data13 = na.omit(data13)

getSymbols('IBM', env = ne, src = 'yahoo', from = '1999-12-01', to ='2018-12-31', adjust = TRUE)
data_IBM <- apply.monthly(ne$IBM$IBM.Adjusted, last)
data14 = coredata(diff(log(data_IBM)))
data14 = na.omit(data14)

s15 = new.env()
getSymbols("INTC", env=s15, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
data15 = apply.monthly(s15$INTC$INTC.Adjusted, last)
data15 = coredata(diff(log(data15)))
data15 = na.omit(data15)

s16 = new.env()
getSymbols("JNJ", env=s16, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
data16 = apply.monthly(s16$JNJ$JNJ.Adjusted, last)
data16 = coredata(diff(log(data16)))
data16 = na.omit(data16)

s17 = new.env()
getSymbols("JPM", env=s17, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
data17 = apply.monthly(s17$JPM$JPM.Adjusted, last)
data17 = coredata(diff(log(data17)))
data17 = na.omit(data17)

s18 = new.env()
getSymbols("MCD", env=s18, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
data18 = apply.monthly(s18$MCD$MCD.Adjusted, last)
data18 = coredata(diff(log(data18)))
data18 = na.omit(data18)

s19 = new.env()
getSymbols("MRK", env=s19, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
data19 = apply.monthly(s19$MRK$MRK.Adjusted, last)
data19 = coredata(diff(log(data19)))
data19 = na.omit(data19)

s20 = new.env()
getSymbols("MSFT", env=s20, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
data20 = apply.monthly(s20$MSFT$MSFT.Adjusted, last)
data20 = coredata(diff(log(data20)))
data20 = na.omit(data20)

s21 = new.env()
getSymbols("NKE", env=s21, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
data21 = apply.monthly(s21$NKE$NKE.Adjusted, last)
data21 = coredata(diff(log(data21)))
data21 = na.omit(data21)

s22 = new.env()
getSymbols("PFE", env=s22, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
data22 = apply.monthly(s22$PFE$PFE.Adjusted, last)
data22 = coredata(diff(log(data22)))
data22 = na.omit(data22)

s22 = new.env()
getSymbols("PG", env=s22, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
data23 = apply.monthly(s22$PG$PG.Adjusted, last)
data23 = coredata(diff(log(data23)))
data23 = na.omit(data23)

getSymbols("TRV", env=s22, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
data24 = apply.monthly(s22$TRV$TRV.Adjusted, last)
data24 = coredata(diff(log(data24)))
data24 = na.omit(data24)

getSymbols("UTX", env=s22, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
dataUTX = apply.monthly(s22$UTX$UTX.Adjusted, last)
data25 = coredata(diff(log(dataUTX)))
data25 = na.omit(data25)


getSymbols("UNH", env=s22, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
dataUNH = apply.monthly(s22$UNH$UNH.Adjusted, last)
data26 = coredata(diff(log(dataUNH)))
data26 = na.omit(data26)

getSymbols("VZ", env=s22, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
dataVZ = apply.monthly(s22$VZ$VZ.Adjusted, last)
data27 = coredata(diff(log(dataVZ)))
data27 = na.omit(data27)

s22 = new.env()
getSymbols("V", env=s22, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
dataV = apply.monthly(s22$V$V.Adjusted, last)
data28 = coredata(diff(log(dataV)))
data28 = na.omit(data28)

getSymbols("WMT", env=s22, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
dataWMT = apply.monthly(s22$WMT$WMT.Adjusted, last)
data29 = coredata(diff(log(dataWMT)))
data29 = na.omit(data29)

getSymbols("WBA", env=s22, src="yahoo", from="1999-12-01", to ="2018-12-31", adjust = TRUE)
dataWBA = apply.monthly(s22$WBA$WBA.Adjusted, last)
data30 = coredata(diff(log(dataWBA)))
data30 = na.omit(data30)



# Calculate the variance of forecaster using quadratic form
# d: vector of dj coefficients (j=0,..., m-2)
# X: log returns

muF<-function(d,X){mean(X)*sum(d)}

varF<-function(d,X){
  M<-length(d)-1
  acfs<- acf(X, plot=F, type="covariance", lag.max=M)$acf
  Gamma<-toeplitz(as.vector(acfs))
  d%*%Gamma%*%as.vector(d)
}

rhoF<-function(d,X){
  M<-length(d)-1
  acfs<- acf(X, plot=F, type ="covariance", lag.max=M+2)$acf
  temp<-d%*%matrix(acfs[abs(outer(0:M,1:(M+1), "-")) +1,,1],
                   M+1, M+1) %*% as.vector(d)
  temp/varF(d,X)
}

corXF<-function(d,X){
  Mp<-length(d)
  acfs<- acf(X, plot=F, type= "covariance", lag.max=Mp)$acf
  sum(d*acfs[-1])/sqrt(acfs[1]*varF(d,X))
}

Hold<-function(rho){pi/acos(rho)}
# m > r >=1
d<-function(m,r){ c((m-r)*((0:(r-1))+1), r*(m-(r: (m-1))-1))}
# retX: log asset return
# m: long-term MA
# r: short-term MA

#return the input dataset's expected return, holding period, correlation, covariance, mean
ruleReturn<-function(retX, m, r){
  vX<-sd(retX)
  mX<-mean(retX)
  mF<-muF(d(m,r),retX)
  vF<-sqrt(varF(d(m,r),retX))
  rXF<-corXF(d(m,r),retX)
  rF<-rhoF(d(m,r),retX)
  ER<-sqrt(2/pi)*vX*rXF*exp(-mF*mF/(2*vF*vF))+mX*(1-2*pnorm(-mF/vF))
  H<-Hold(rF)
  list("ER"=ER, "H"=H, "rhoF"=rF, "VF"=vF, "muF"=mF,
       "corXF"=rXF)
}

# referece: last year's STA457 assignment 1 answer code 

#combine all the data in to one matrix 
library('rowr')
returns <- cbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12, data13, data14, data15,data16, data17, data18, data19, data20, data21, data22, data23,
                      data24, data25, data26, data27, data28, data29, data30, fill= 0)

#define a function 'cal' to calculate optimal m and r for a single stock 
cal <- function(stock){
  holder <- numeric(0)
  m <- numeric(0)
  r <- numeric(0)
  for (i in 2:11){
    for (j in (i+1):12){
      if (i<j){
        holder <- c(holder, ruleReturn(stock, m=j, r=i) [[1]])
        m <- c(m, j)
        r <- c(r, i)
      }
    }
  }
  m_optimal <- m[which.max(holder)]
  r_optimal <- r[which.max(holder)]
  list(optimal_m = m_optimal, optimal_r = r_optimal)
}

#create an empty list to hold all the optimal m and r result for all 30 stocks
optimal <- c()
#using loop to apply cal function to each of the stock and add the result to the optimal list; add the name of the
#stock to the end of its optimal m and r combination.
for (i in 1:30){
  stocks <- return[[i]]
  optimal <- c(optimal, cal(stocks), name = dj30[i])
}
optimal

# Input data are the stock price, we got X_t from stock price by using x_t = ln(p_t/p_t-1). X_t follows time series. 
# F_t are calculated from prices. When F_t is greater than 0, B_t = 1. When F_t is less than 0, B_t = -1.
#The sign of B_t decides either buy (+1) or sell (-1) the stock based on the change in price.
# From B_t-1 and X_t, we got R_t, which is the ruled return of the stock at time t. 
# Holding period can be used as anthoer standard of stock return as well. Longer holding period, aka larger H, indicates higher return.
# According to the formula given in the appendix last year, double MA = (m-r)(j+1) when 0 <=j <=r-1
# double MA = r(m-j-1) when r<=j<=m-1. If we want to calculate d_j, we need to know m and r 
# Expectation of ruled return is calculated for each time point t.
# The optimal combinations of m and r for all the stocks are calculated by choosing the m and r that generating the highest expectated r.



#QA 2) 
library('rowr')
returns <- cbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12, data13, data14, data15,data16, data17, data18, data19, data20, data21, data22, data23,
                      data24, data25, data26, data27, data28, data29, data30, fill= 0)

#Equal Weighted
#To calculate the equal weighted portfolio, we first imported all the data from yahoo.
#We constructed a matrix "returns" that contains all the monthly-returns per stock from the year 1999 to
#the year 2018 (228 months in total) for all 30 stocks (so "returns" contains 30 * 228 months total for
#all 30 stocks).
#since all the stocks are equal weighted,
#we constructed Br to be the vector that contains the mean of each stocks monthly returns over 
#228 months.
#We can then use the build in function mean() and var() to calculate the mean and variance 
#of Br to see the mean and variance of all 30 stocks.
#note that the number is very small, it is because data27 (Visa) stock only went public after 2018.
Br <- numeric(nrow(returns))
for(i in 1:nrow(returns)){
  Br[i] <- 1/nrow(returns) * sum(returns[,i], na.rm = TRUE)
}

#Risk Parity 
#We used a similar approach to calculate the RP portfolio as we used for calculating EW.
#we calculated the standard deviation (sds) of the monthly return for all 30 stocks.
#given that we have the standard deviation for all 30, we can use it to calculate the
#weighted sum by using the Hint b) given to us in the assignment 
#sheet.
#we iterate the weighted sum over all monthly return for all 228 months for all 30 stocks.
#we assign Rp to be a list of risk-parity allocation for all 30 stocks calculated from
#their respective weighted sums.
#We can then use the build in function mean() and var() to calculate the mean and variance 
#of Rp to see the mean and variance of all 30 stocks.
#note that the number is normal because data27 (Visa) stock's allocation has been changed
#due to the weighted sum calculation. Therefore, this portfolio is better than using
#equal weighted and gives a better number.
Rp <- numeric(nrow(returns))
sds <- apply(returns, 2, sd, na.rm = TRUE)
ws <- (1/sds)/(sum(1/sds))
for (i in 1:nrow(returns)) {
  Rp[i] <- sum(returns[i,] * ws, na.rm = TRUE)
}


#annualize
#We use sharpe ratio to calculate the annualized performance, which is the ratio of 
#the average return earned in excess of the risk-free rate per unit of volatility/total risk.
#the function will use annualized expected returns (mean), volatilities(sqr(var))
#and the riskfree rate of 0.02 to calculate the performance of the portfolio.
#You can call either annualize(Br) or annualize(Rp) to see their respective
#performances.
annualize <- function (strat){
  riskFreeRate <- 0.02
  return((12 * mean(strat) - riskFreeRate)/(sqrt(12) * sqrt(var(strat))))
}


#QB 1)

#sigma_t is the ex-ante volatility estimate
#Volatility is the measurement of the variation in trading price over time. Ex-ante means the prediction before the future events take place.
#Therefore, sigma_t can be used to predict future variability in prices.

# delta = 0.2 is given by the professor 
delta = 0.2

# Sigma squared function is given in the assignment and we define a function called calSigma based on that.
# CalSigma can calcualte sigma_t for each time slot and return a vector of sigma_t for that dataset.
# Two for loops are involved to go through each line of the data and for each month
# t starts from 13 since t-1-i cannt less than 0.
calSigma <- function(data){
  sig <- c()
  for (t in 13:nrow(returns)){
    acc <- c()
    for (i in 0:11){
      acc <- c(acc, ((1-delta)*(delta^i)*(data[t-1-i]-mean(data))^2))
    }
    sig <- c(sig, sqrt(12*sum(acc)))
  }
  return(sig)
}
#calculate sigma for all 30 companies
calAll30Sig <- function(){
  stocks <- c('data1', 'data2', 'data3', 'data4', 'data5', 'data6', 'data7', 'data8', 'data9', 'data10', 'data11', 'data12', 'data13', 'data14', 'data15','data16', 'data17', 'data18', 'data19', 'data20', 'data21', 'data22', 'data23',
              'data24', 'data25', 'data26', 'data27', 'data28', 'data29', 'data30')
  list1 <- list()
  for(i in 1:30)
    list1[[i]]<- calSigma(get(stocks[i]))
  return(list1)
}

#calculate mean sigma for all 30 stocks
calAll30SigMean <- function(){
  all30 <- calAll30Sig()
  data <- c()
  datamean <- 0
  result <- c()
  for (i in 1:30) {
    data <- all30[[i]]
    datamean <- mean(data)
    result <- c(result, datamean)
  }
  return(result)
}

#B(2)
#define x and y for the linear regression model
#y = r_(s,t)/ sigma_(s,t-1)
#x = r_(s,t-h)/ sigma_(s,t-h)
#y is the excess return in month t.
#x is ys return lagged h months.
#we learned in previous courses (Sta302, Sta303) that the R-Squared
#calculated in a linear regression model measures how close the data
#are to the fitted line. In this question, when we try to find the
#optimal h for the predicative regressions for all 30 DJ constituents,
#we are basically trying to maximize R-Squared. In which, will give us
#the h that maximize the excess return.
# h+1 to make x and y in the same length

#calculate max r-squared
calmaxRsquared <- function(data,i){
  rst_sigmast <- NULL
  list2 = calAll30()
  rst_sigmast <- c(rst_sigmast, data[13:nrow(data)]/list2[[i]])
  rsqr <- c()
  m <- c(1:12)
  for(h in 1:12){
    y <- rst_sigmast[(h+1):length(rst_sigmast)]
    x <- sign(rst_sigmast[(h+1): length(rst_sigmast)])
    model <- lm(y ~ x)
    
    rsqr <- c(rsqr,summary(model)$r.squared)
  }
  final <- m[which.max(rsqr)]
  return(final)
}

#calculate for all 30
calmaxRsquared30 <- function(){
  stocks <- c('data1', 'data2', 'data3', 'data4', 'data5', 'data6', 'data7', 'data8', 'data9', 'data10', 'data11', 'data12', 'data13', 'data14', 'data15','data16', 'data17', 'data18', 'data19', 'data20', 'data21', 'data22', 'data23',
              'data24', 'data25', 'data26', 'data27', 'data28', 'data29', 'data30')
  rSquared30 <- c()
  for (i in 1:30){
    rSquared30 <- c(rSquared30, calmaxRsquared(get(stocks[i]), i))
  }
  return(rSquared30)
}


#b 3)
#list3 is created to hold all the sigma_t for all 30 stocks.
list3 <- calAll30()
#a is a holder for the result of sign(r_s,t-h:t)*40%/sigma_s,t*r_s,t:t+1
a <- numeric(nrow(returns)-12)
#tsmom is a empty holder of 30 elements waiting for the tsmom result of each stock
tsmom <- numeric(30)
# two loops of s and t are used to go over each time t in each stock, s indicates which stock we are dealing with and
# t indicates which time slot we are working on.
# h_s represent the holding period and r symbolize the return.
# for each time t, we begin with calculate position using sign(r_s, t-h:t)*40%/sigma_s
# After calculating the position, we multiply the position vector with the rule return of stock s at time t to t+1
# At the end of 30 loops, we sum all the result together without NA in the vector and divide it by 30.
for (s in 1:30){
  for (t in 13:nrow(returns)){
    position_s <- sign(returns[t-12:t, s])*(40/100)/(list3[[s]][t])
    a[t] <- sum(position_s * returns[t:t+1, s])
  }
  tsmom[s] <- sum(na.omit(a))
}
total_tsmom <- sum(tsmom)/30
total_tsmom

#rTSMOM_t,t+1 represents the time series momentum. Previous data is used to produce time series momentum so that 
#we can make prediction in that stock's performace and earn more profit as well as avoid loss. 


#C 1)
# The function ER is defined based on the formula provided in the announcement by professor Lin
# m represents long term moving average
# h represents lagged time
# d is a vector consists of dj coefficents 
# retx is the input dataset we want to work on

ER <- function(m, h, d, retx){
  sum1 <- c()
  sum2 <- c()
  sum3 <- c()
  result <- 0
  for (t in 13:nrow(returns)){
    for (i in 1:h-1){
      for (j in 1: m-2){
        rule <- ruleReturn(retx, j, i)
        sum1 <- c(sum1, sum(d * as.vector(corXF(t, t+i-j) - rule[[1]]^2)))
        print(sum1)
      }
      sum2 <- c(sum2, sum(sum1))
    }
    sum3 <- c(sum3, sum(sum2))
  }
  result <- sum(sum3)
  return(result)
}


#C 2)
# define a function called cal_ER to extract ER for each stock
# similar steps are taken as the part b in question A
# the function return the optimal m and h for each stock 
cal_ER <- function(stock){
  result_c <- numeric(0)
  m <- numeric(0)
  r <- numeric(0)
  for (i in 2:11){
    for (j in (i+1):12){
      if (i<j){
        holder <- c(result_c,  ER(m =j, h=i, d=d, get(stocks[i])))
        m <- c(m, i)
        h <- c(h, j)
      }
    }
  }
  m_optimal <- m[which.max(holder)]
  h_optimal <- h[which.max(holder)]
  list(optimal_m = m_optimal, optimal_h = h_optimal)
}

#create an empty list to hold all the optimal m and h result for all 30 stocks
optimal_c <- c()
#using loop to apply cal_ER function to each of the stock and add the result to the optimal_c list; add the name of the
#stock to the end of its optimal m and h combination.
for (i in 1:30){
  stocks <- return[[i]]
  optimal_c <- c(optimal_c, cal_ER(stocks), name = dj30[i])
}
optimal_c
