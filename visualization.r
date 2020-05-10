library(ggplot2)

# fake data

plot_I_intervals = function(expected_I, bootstrap_samples, confidence) {
  "
  Plots the daily expected number of new cases with confidence intervals.
  
  Parameters
  ----------
  expected_I : numeric vector
    Daily expected number of new cases.
  bootstrap_samples : DataFrame
    Bootstrapped samples of the time series. Each row must be a time series
    of the same length as expected_I.
  confidence : numeric
    Confidence of the intervals to be drawn.
  "
  # extract confidence intervals
  lower = function(c) quantile(c, (1 - confidence) / 2)
  upper = function(c) quantile(c, confidence + ((1 - confidence) / 2))
  lower_bound = apply(bootstrap_samples, MARGIN=2, FUN=lower)
  upper_bound = apply(bootstrap_samples, MARGIN=2, FUN=upper)
  # plot
  df = data.frame(expected_I)
  colnames(df) = "I(t)"
  df["t"] = 1:nrow(df)
  ggplot(df, aes(x=t, y=`I(t)`)) +
    geom_line() +
    geom_ribbon(aes(ymin=lower_bound, ymax=upper_bound), alpha=0.2, fill="blue")
}