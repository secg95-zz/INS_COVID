bootstrap_samples = function( expected_I, observed_I, number_samples, window_size=3 , use_wa = FALSE, beta = 0.9){
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

      if(use_wa){
        DES = 1:length(expected_I)
        temp_de = standarized_residual$DE
        for(j in 2:length(observed_I))
        { 
          temp_de[j] = beta*temp_de[j - 1] + (1-beta)*temp_de[j]
        }
        sample = sample*temp_de
      }else{
        sample = sample*standarized_residual$DE
      }
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
  half_window_size = floor(window_size / 2)
  standarized_residual <- observed_I - expected_I
  DES = 1:length(expected_I)
  for( i in seq(half_window_size + 1, length(expected_I) - half_window_size + 1))
  {
    indices = (i-half_window_size):(i + half_window_size-1)
    DE = sd(observed_I[indices] - expected_I[indices])
    standarized_residual[i] = standarized_residual[i]/DE 
    DES[i] = DE
  }

  first_window_idexes = 1:half_window_size
  last_window_indexes = (length(expected_I) - half_window_size + 1):length(expected_I)
  DE_in = sd(observed_I[first_window_idexes] - expected_I[first_window_idexes])
  DE_out = sd(observed_I[last_window_indexes] - expected_I[last_window_indexes])
  DES[first_window_idexes] = DE_in
  DES[last_window_indexes] = DE_out
  standarized_residual[first_window_idexes] = standarized_residual[first_window_idexes]/DE_in
  standarized_residual[last_window_indexes] = standarized_residual[last_window_indexes]/DE_out 
  return(list("SR"=standarized_residual, "DE"=DES))
}


compute_variances = function(observed_I, window_size=3){
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
  variance_observed = 1:length(observed_I)
  half_window_size = floor(window_size / 2)
  for( i in seq(half_window_size + 1, length(observed_I) - half_window_size + 1))
  {
    indices = (i-half_window_size):(i + half_window_size-1)
    variance = var(observed_I[indices])
    variance_observed[i] = variance
  }

  first_window_idexes = 1:half_window_size
  last_window_indexes = (length(observed_I) - half_window_size + 1):length(observed_I)
  DE_in = var(observed_I[first_window_idexes])
  DE_out = var(observed_I[last_window_indexes])
  variance_observed[first_window_idexes] = DE_in
  variance_observed[last_window_indexes] = DE_out
  return(variance_observed)
}