
MECHA_R <- function(path = "granar_examples/MECHAv4_light.py", Kx = F){
  
  tic <- Sys.time()
  message("Launch MECHA")
  MECHA <- try(source_python(path), silent = T)
  
  kr_tot <- py_to_r_2(kr_tot)
  if(is.null(MECHA)){MECHA <- c(0)}
  if(str_detect(MECHA[1], "Error ")){
    message("error NaN: kr_tot <- NA")
    kr_tot <- c(NA, NA, NA)
  }
  if(MECHA[1] == "Error in py_run_file_impl(file, local, convert) : \n  IndexError: index 25 is out of bounds for axis 0 with size 25\n"){
    message("Undersized matrix of connected cell.")
    kr_tot <- c(NA, NA, NA)
  }
  print(Sys.time()-tic)
  if(Kx){
    K_xyl_spec <- py_to_r_2(K_xyl_spec)
    K <- data.frame(kr_tot,K_xyl_spec)
    return(K)
  }else{
    K <- kr_tot
    return(K)
    }
}

py_to_r_2 <- function(x){
  
  ele <- str_split(paste0(x), " ")[[1]]
  kr <- as.numeric(unlist(regmatches(ele,gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?",ele))))
  return(kr)
}