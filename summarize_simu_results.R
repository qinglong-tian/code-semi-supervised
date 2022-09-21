compare_mse <- function(summary_simulation, index)
{
  mean(sapply(summary_simulation, function(x) {
    x$mse3[index]
  })) -> mse3
  mean(sapply(summary_simulation, function(x) {
    x$mse4[index]
  })) -> mse4
  mean(sapply(summary_simulation, function(x) {
    x$mse34[index]
  })) -> mse34
  
  valVec <- c(mse3, mse4, mse34)
  names(valVec) <- c("col3", "col4", "comb")
  
  return(valVec)
}

compare_cp <- function(summary_simulation, index)
{
  mean(sapply(summary_simulation, function(x) {
    x$cp3[index]
  })) -> cp3
  mean(sapply(summary_simulation, function(x) {
    x$cp4[index]
  })) -> cp4
  mean(sapply(summary_simulation, function(x) {
    x$cp34[index]
  })) -> cp34
  
  cpVec <- c(cp3, cp4, cp34)
  names(cpVec) <- c("col3", "col4", "comb")
  return(cpVec)
}