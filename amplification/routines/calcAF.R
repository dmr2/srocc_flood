get_AF <- function(rp,z,Ne_hist,histCurveSamps,slr) {
  
  q <- 1/rp
  
  # The height of the historical flood with frequency q
  hist_height <- round(z[which.min(abs(Ne_hist-q))],2)
  
  mi <- NULL
  ElogNprojshift <- NULL
  
  randi <- sample(1:1000, 1000, replace=T) # random GPD parameter samples
  
  for(iii in 1:1000){
    mi <- which.min( abs(z-(hist_height-slr)) )
    ElogNprojshift[iii] = log(histCurveSamps[randi[iii],mi])
  }
  
  af = exp(ElogNprojshift-log(q))
  
  return ( list(AF=af,HistHeight=hist_height)  )
}