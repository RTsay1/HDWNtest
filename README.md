# HDWNtest
# The R code: HDWNtest.R contains the R code for running the high-dimensional white noise test of Tsay (2020)
# The input includes:
# (a) xt: a T-by-k data matrix of time series, each column is a time series. Thus, T is the sample size and k is the number of series
# (b) lag: The number of lags of autocorrelation matrices used in the test (default is 10)
# (c) alpha: type-I error of the white-noise test (default are 10%, 5% and 1%)
# (d) subfr: the proportion of series to be used in test when k >= T. That is, when the number of series is greater than or equal to the sample size,
#      the test uses the subfr times k to perform the test. Default is 0.75.
# (e) output: a logical input to control the output. Default is TRUE.
#
# Output includes:
# (a) Tst: the test statistics for each lag
# (b) Test: The overall test statistics
# (c) pv: p-value of the joint tests
# (d) CVm: critical values of the joint tests
# (e) pvi: p-value of individual lag test
# (f) CVi: critical values of individual lag test.
