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
  mean(sapply(summary_simulation, function(x) {
    x$mse_add[index]
  })) -> mseadd
  mean(sapply(summary_simulation, function(x) {
    x$mse_debias[index]
  })) -> msedebias
  
  valVec <- c(mse3, mse4, mse34, mseadd, msedebias)
  names(valVec) <- c("m1", "m2", "comb", "additive", "debias")
  
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
  mean(sapply(summary_simulation, function(x) {
    x$cp_add[index]
  })) -> cpadd
  mean(sapply(summary_simulation, function(x) {
    x$cp_debias[index]
  })) -> cpdebias
  
  cpVec <- c(cp3, cp4, cp34, cpadd, cpdebias)
  names(cpVec) <- c("col3", "col4", "comb", "additive", "debias")
  return(cpVec)
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
  mean(sapply(summary_simulation, function(x) {
    x$cp_add[index]
  })) -> cpadd
  mean(sapply(summary_simulation, function(x) {
    x$cp_debias[index]
  })) -> cpdebias
  
  cpVec <- c(cp3, cp4, cp34, cpadd, cpdebias)
  names(cpVec) <- c("col3", "col4", "comb", "additive", "debias")
  return(cpVec)
}

compare_se <- function(summary_simulation, index)
{
  mean(sapply(summary_simulation, function(x) {
    x$se3[index]
  })) -> cp3
  mean(sapply(summary_simulation, function(x) {
    x$se4[index]
  })) -> cp4
  mean(sapply(summary_simulation, function(x) {
    x$se34[index]
  })) -> cp34
  mean(sapply(summary_simulation, function(x) {
    x$se_add[index]
  })) -> cpadd
  mean(sapply(summary_simulation, function(x) {
    x$se_debias[index]
  })) -> cpdebias
  
  cpVec <- c(cp3, cp4, cp34, cpadd, cpdebias)
  names(cpVec) <- c("col3", "col4", "comb", "additive", "debias")
  return(cpVec)
}
