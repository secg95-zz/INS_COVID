bootstrap_samples = function( expected_I, observed_I, number_samples, window_size=3){
  "
  Once the expected_I vector is computed is posible to generate the bootsrapped
  samples for  a confidence interval cosntruction.
  
  Parameters
  ----------
  Expected_I : numeric vector
    Number of sintomatic people on a given day (estimated).
  Observed_I : numeric vector
    Number of sintomatic people on a given day.
  number_samples : numeric
    number of samples.
  window_size : numeric
    smoothing factor.
  
  Returns
  -------
  standarized_residual : data_frame
    bootstrapped sample
  "
  standarized_residual <- standarized_residuals(expected_I, observed_I, window_size)
  samples = list()
  for( i in 1:number_samples)
  {
    sample <- sample(standarized_residual$SR, length(expected_I), replace = TRUE)
    sample = sample*standarized_residual$DE
    sample_I = expected_I + sample
    samples <- list.append(samples, i=sample_I)
  }
  
  df <- data.frame(matrix(unlist(samples), nrow=length(samples), byrow=T))
  return(df)
}

standarized_residuals = function(expected_I, observed_I, window_size=3){
  "
  Residuals normalization for an umbaised bootstrap.
  
  Parameters
  ----------
  Expected_I : numeric vector
    Number of sintomatic people on a given day (estimated).
  Observed_I : numeric vector
    Number of sintomatic people on a given day.
  window_size : numeric
    smoothing factor.
  
  Returns
  -------
  standarized_residual : list
    SR : residuals
    DE : standard deviation
  "
  standarized_residual <- observed_I - expected_I
  DES = 1:length(expected_I)
  for( i in seq(1, length(expected_I) - (length(expected_I) %% window_size), window_size))
  {
    DE = sd(observed_I[i:(i + window_size-1)] - expected_I[i:(i + window_size-1)])
    standarized_residual[i:(i + window_size-1)] = standarized_residual[i:(i + window_size-1)]/DE 
    DES[i:(i + window_size-1)] = DE
  }
  
  if (length(expected_I) %% window_size != 0)
  {
    aux = length(expected_I) %% window_size
    DE = sd(observed_I[(length(expected_I) - aux):length(expected_I)] - expected_I[(length(expected_I) - aux):length(expected_I)])
    DES[(length(expected_I) - aux):length(expected_I)] = DE
    standarized_residual[(length(expected_I) - aux):length(expected_I)] = standarized_residual[(length(expected_I) - aux):length(expected_I)]/DE
  }
  return(list("SR"=standarized_residual, "DE"=DES))
}