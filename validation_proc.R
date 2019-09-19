
library(tidyverse)
library(granar)
library(readxl)

dir = getwd()

RS <- read_excel("./RootScan_data/Burton_2013.xlsx")
RS <- RS%>%
  mutate(r_CT = sqrt(RXSA/pi),
         r_stele = sqrt(TSA/pi),
         d_stele = r_stele*2,
         r_cortex = r_CT - r_stele,
         one_cell_area = (CCA/CC)*0.75,
         cell_dia_cortex = (r_CT-r_stele)/round(CF),
         n_layers_cortex = round(CF)-2, # endodermis and exodermis should not be associate to cortex cell
         cell_dia_stele = 0.011, # between 9.4 and 13.5 micrometer (Steudle et al. 1993; Gao et al. 2015)
         n_X = round(r_stele*33), # olometric link between the stele radius and the metaxylem size (Moreno-ortega 2016)
         one_x = XVA/n_X,
         max_size = 2*sqrt(one_x/pi))

RS_tot <- rbind(RS %>% mutate(id = 1), RS %>%  mutate(id = 2), RS %>%  mutate(id = 3))
RS_tot <- RS_tot%>%
  mutate(repet = id,
         id_sim = c(1:nrow(RS_tot)),
         name = "Burton et al. 2013",
         assession = paste0(asses,sion))

params <- read_param_xml("./input/Zea_mays_Heymans_2019.xml")
sim <- create_anatomy(parameters = params)
plot_anatomy(sim)
k <- 1
sim <- list()
for(id_sim in k:nrow(RS_tot)){

  temp <- RS_tot[id_sim,]
  message(paste0(">>>>>> ", id_sim , " out of  ", nrow(RS_tot), " <<<<<<<"))
  # Change the parameters --------------
  #Stele
  params$value[params$name == "stele" & params$type == "cell_diameter"] <- temp$cell_dia_stele
  params$value[params$name == "stele" & params$type == "layer_diameter"] <- temp$d_stele-2*temp$cell_dia_stele
  params$value[params$name == "xylem" & params$type == "max_size"] <- temp$max_size - 0.011/2
  params$value[params$name == "xylem" & params$type == "n_files"] <- temp$n_X
  params$value[params$name == "phloem" & params$type == "max_size"] <- temp$cell_dia_stele
  params$value[params$name == "pericycle" & params$type == "cell_diameter"] <- temp$cell_dia_stele
  params$value[params$name == "endodermis" & params$type == "cell_diameter"] <- temp$cell_dia_cortex
  
  #Cortex
  params$value[params$name == "cortex" & params$type == "cell_diameter"] <- temp$cell_dia_cortex
  params$value[params$name == "cortex" & params$type == "n_layers"] <- temp$n_layers_cortex
  params$value[params$name == "exodermis" & params$type == "cell_diameter"] <- temp$cell_dia_cortex
  params$value[params$name == "epidermis" & params$type == "cell_diameter"] <- temp$cell_dia_cortex*0.33
  params$value[params$name == "aerenchyma" & params$type == "proportion"] <- temp$p_A/100
  
  t_before <- Sys.time()
  # Run GRANAR ------------
  sim <- list()
  n_try <- 1
  while (TRUE) {
    did_it <- 1
    message("GRANAR in process")
    sim <- try(create_anatomy(parameters = params, verbatim=F), silent = TRUE)
    n_try <- n_try + 1
    if(is.null(sim)){
      sim <- try(error("nerror = ..."), silent = T)
      did_it <- 0
    }
    if(n_try > 100){
      sim <- list()
      sim$output$value <- NA
      did_it <- 0
    }
    if(!is(sim,'try-error')) break
  }
  if(did_it == 1){
    write_anatomy_xml(sim = sim,
                      path = "./MECHA/cellsetdata/current_root.xml")
    fc <- file.copy(from = "./MECHA/cellsetdata/current_root.xml",
                    to = paste0("./MECHA/cellsetdata/root_",temp$assession,"_",temp$repet,".xml"),
                    overwrite = T)
    output <- output_as_RootScan(sim)
  }else{
    output[,1:12] <- NA
  }
  
  RS_tot[id_sim, "model_area"] <- output$m_RXSA
  RS_tot[id_sim, "model_TSA"] <- output$m_TSA
  RS_tot[id_sim, "model_TSA_2"] <- output$m_TSA_2
  RS_tot[id_sim, "model_TCA"] <- output$m_TCA
  RS_tot[id_sim, "model_CCA"] <- output$m_CCA
  RS_tot[id_sim, "model_XVA"] <- output$m_XVA
  RS_tot[id_sim, "model_CC"] <- output$m_CC
  RS_tot[id_sim, "model_AA"] <- output$m_AA
  RS_tot[id_sim, "model_pA"] <- output$m_pA
  RS_tot[id_sim, "sim_time_G"] <- output$granar_time
  
  k <- k + 1
}

RS_tot <- RS_tot%>%
  mutate(model_pA = (1-model_CCA/model_TCA)*100,
         model_cortex_w = sqrt(model_area/pi)-sqrt(model_TSA/pi))

# Working with Python in R ------------------
# Download Anaconda for python 
# Set a environment for MECHA with numpy 1.15.2, scipy 1.1.0, networkx 2.3, lxml 4.2.5, maplotlib 3.0.0
library(reticulate)
use_python("../AppData/Local/Continuum/envs/MECHA")
os <- import("os")
numpy <- import("numpy")
scipy <- import("scipy")
networkx <- import("networkx")
lxml <- import("lxml")
matplotlib <- import("matplotlib")
source("MECHA_R.R")

# MECHA loop on the generate cross-section -------------
k <- 1
for(id_sim in k:nrow(RS_tot)){
  
  temp <- RS_tot[id_sim,]
  message(paste0(">>>>>> ", k , " out of  ", nrow(RS_tot),". In process: ",temp$assession," ",temp$repet, " <<<<<<<"))
  
  fc <- try(file.copy(from = paste0("MECHA/cellsetdata/root_",temp$assession,"_",temp$repet,".xml"),
                      to = "MECHA/cellsetdata/current_root.xml",
                      overwrite = T), silent = T)
  if(!str_detect(fc[1], "Error :")){
    K <- MECHA_R(path = "MECHA/MECHAv4_light.py", Kx = F)
    
    fc <- file.copy(from = "MECHA/Projects/granar/out/GRANAR/Project_Test/baseline/Macro_prop_1,0.txt",
                    to = paste0("MECHA/cellsetdata/Burton/macro/Macro_prop_1,0_",id_sim,".txt"),
                    overwrite = T)
    fc <- file.copy(from = "MECHA/Projects/granar/out/GRANAR/Project_Test/baseline/Macro_prop_2,1.txt",
                    to = paste0("MECHA/cellsetdata/Burton/macro/Macro_prop_2,1_",id_sim,".txt"),
                    overwrite = T)
    fc <- file.copy(from = "MECHA/Projects/granar/out/GRANAR/Project_Test/baseline/Macro_prop_4,2.txt",
                    to = paste0("MECHA/cellsetdata/Burton/macro/Macro_prop_4,2_",id_sim,".txt"),
                    overwrite = T)
    
  }else{
    K <- c(NA,NA,NA)
  }
  
  RS_tot[id_sim, "model_kr0"] <- K[1]
  RS_tot[id_sim, "model_kr1"] <- K[2]
  RS_tot[id_sim, "model_kr2"] <- K[3]
  
  k <- k + 1
}

write.csv(RS_tot, "RS_Burton.csv")
# Burton et al. 2013 processed ----


RS_Gao <- read_excel("./RootScan_data/Gao_2015.xlsx")
RS_Gao <- rbind(RS_Gao %>% mutate(id = 1), RS_Gao %>%  mutate(id = 2), RS_Gao %>%  mutate(id = 3))
RS_Gao <- RS_Gao%>%
  mutate(repet = id,
         id_sim = c(1:nrow(RS_Gao)),
         name = "Gao et al. 2015")

k <- 1
for(id_sim in 1:nrow(RS_Gao)){

  temp <- RS_Gao[id_sim,]
  message(paste0(">>>>>> ", id_sim , " out of  ", nrow(RS_Gao), " <<<<<<<"))
  # Change the parameteRS_Gao --------------
  #Stele
  params$value[params$name == "stele" & params$type == "cell_diameter"] <- temp$diam_stele/1000
  params$value[params$name == "stele" & params$type == "layer_diameter"] <- temp$stele/1000-3*temp$diam_stele/1000
  
  params$value[params$name == "xylem" & params$type == "max_size"] <- temp$dia_xyl/100
  params$value[params$name == "xylem" & params$type == "n_files"] <- temp$n_xylem
  
  params$value[params$name == "phloem" & params$type == "max_size"] <- temp$diam_stele/1000
  
  params$value[params$name == "pericycle" & params$type == "cell_diameter"] <- temp$diam_stele/1000
  
  params$value[params$name == "endodermis" & params$type == "cell_diameter"] <- temp$diam_cortex/1000
  
  #Cortex
  params$value[params$name == "cortex" & params$type == "cell_diameter"] <- temp$diam_cortex/1000
  params$value[params$name == "cortex" & params$type == "n_layers"] <- temp$n_layer_cortex-2
  
  params$value[params$name == "exodermis" & params$type == "cell_diameter"] <- temp$diam_cortex/1000
  
  params$value[params$name == "epidermis" & params$type == "cell_diameter"] <- (temp$diam_cortex/1000)*0.33
  
  params$value[params$name == "aerenchyma" & params$type == "proportion"] <- temp$aerenchyma
  params$value[params$name == "aerenchyma" & params$type == "n_files"] <- 10
  
  sim <- list()
  granar_time <- NULL
  
  # Run the model -----------------
  
  # Run GRANAR ------------
  sim <- list()
  n_try <- 1
  while (TRUE) {
    did_it <- 1
    message("GRANAR in process")
    sim <- try(create_anatomy(parameters = params, verbatim=F), silent = TRUE)
    n_try <- n_try + 1
    if(is.null(sim)){
      sim <- try(error("nerror = ..."), silent = T)
      # sim$output$value <- NA
      did_it <- 0
    }
    # if(sim[1] == "Error in 1:nrow(cell_size) : argument of length 0\n" ){
    #   sim <- list()
    #   sim$output$value <- NA
    #   did_it <- 0
    # }
    if(n_try > 5){
      sim <- list()
      sim$output$value <- NA
      did_it <- 0
    }
    if(!is(sim,'try-error')) break
  }
  if(did_it == 1){
    write_anatomy_xml(sim = sim, 
                      path = "MECHA/cellsetdata/current_root.xml")
    fc <- file.copy(from = "MECHA/cellsetdata/current_root.xml",
                    to = paste0("MECHA/cellsetdata/Zhengdan958_",temp$ID[1],"_",temp$repet[1], ".xml"),
                    overwrite = T)
    
    output <- output_as_RootScan(sim)
  }else{
    output[,1:12] <- NA
  }
  
  RS_Gao[id_sim, "model_area"] <- output$m_RXSA
  RS_Gao[id_sim, "model_TSA"] <- output$m_TSA
  RS_Gao[id_sim, "model_TSA_2"] <- output$m_TSA_2
  RS_Gao[id_sim, "model_TCA"] <- output$m_TCA
  RS_Gao[id_sim, "model_CCA"] <- output$m_CCA
  RS_Gao[id_sim, "model_XVA"] <- output$m_XVA
  RS_Gao[id_sim, "model_CC"] <- output$m_CC
  RS_Gao[id_sim, "model_AA"] <- output$m_AA
  RS_Gao[id_sim, "model_pA"] <- output$m_pA
  RS_Gao[id_sim, "sim_time_G"] <- output$granar_time
}


# MECHA loop on the generate cross-section -------------
k <- 1
for(id_sim in k:nrow(RS_Gao)){
  
  temp <- RS_Gao[id_sim,]
  message(paste0(">>>>>> ", k , " out of  ", nrow(RS_Gao),". In process: ",temp$ID[1]," ",temp$repet, " <<<<<<<"))
  
  fc <- try(file.copy(from = paste0("MECHA/cellsetdata/root_",temp$ID,"_",temp$repet,".xml"),
                      to = "MECHA/cellsetdata/current_root.xml",
                      overwrite = T), silent = T)
  if(!str_detect(fc[1], "Error :")){
    K <- MECHA_R(path = "MECHA/MECHAv4_light.py", Kx = F)
    
    fc <- file.copy(from = "MECHA/Projects/granar/out/GRANAR/Project_Test/baseline/Macro_prop_1,0.txt",
                    to = paste0("MECHA/cellsetdata/Gao/macro/Macro_prop_1,0_",id_sim,".txt"),
                    overwrite = T)
    fc <- file.copy(from = "MECHA/Projects/granar/out/GRANAR/Project_Test/baseline/Macro_prop_2,1.txt",
                    to = paste0("MECHA/cellsetdata/Gao/macro/Macro_prop_2,1_",id_sim,".txt"),
                    overwrite = T)
    fc <- file.copy(from = "MECHA/Projects/granar/out/GRANAR/Project_Test/baseline/Macro_prop_4,2.txt",
                    to = paste0("MECHA/cellsetdata/Gao/macro/Macro_prop_4,2_",id_sim,".txt"),
                    overwrite = T)
    
  }else{
    K <- c(NA,NA,NA)
  }
  
  RS_Gao[id_sim, "model_kr0"] <- K[1]
  RS_Gao[id_sim, "model_kr1"] <- K[2]
  RS_Gao[id_sim, "model_kr2"] <- K[3]
  
  k <- k + 1
}

write.csv(RS_Gao, "RS_Gao.csv")

# Gao et al. 2015 processed ----

RS_Chimungu <- read_excel("RootScan_data/Chimungu_2014.xlsx")

RS_t <- rbind(RS_Chimungu %>% mutate(id = 1), RS_Chimungu %>%  mutate(id = 2), RS_Chimungu %>%  mutate(id = 3))
RS_Chimungu <- RS_t%>%
  mutate(repet = id,
         id_sim = c(1:nrow(RS_t)),
         name = "Chimungu et al.2014")
k <- 1
for(id_sim in k:nrow(RS_Chimungu)){
  
  temp <- RS_Chimungu[id_sim,]
  message(paste0(">>>>>> ", id_sim , " out of  ", nrow(RS_Chimungu), " <<<<<<<"))
  # Change the parameteRS_Chimungu --------------
  #Stele
  params$value[params$name == "stele" & params$type == "cell_diameter"] <- 0.015
  params$value[params$name == "stele" & params$type == "layer_diameter"] <- temp$SD-2*0.015
  params$value[params$name == "xylem" & params$type == "max_size"] <- temp$SD*0.0762 + 0.0325
  params$value[params$name == "xylem" & params$type == "n_files"] <- round((temp$SD/2)*33)
  params$value[params$name == "phloem" & params$type == "max_size"] <- 0.015
  params$value[params$name == "pericycle" & params$type == "cell_diameter"] <- 0.015
  params$value[params$name == "endodermis" & params$type == "cell_diameter"] <- temp$cortex_dia
  #Cortex
  params$value[params$name == "cortex" & params$type == "cell_diameter"] <- temp$cortex_dia
  params$value[params$name == "cortex" & params$type == "n_layers"] <- round(temp$CF)-2
  params$value[params$name == "exodermis" & params$type == "cell_diameter"] <- temp$cortex_dia
  params$value[params$name == "epidermis" & params$type == "cell_diameter"] <- (temp$cortex_dia)*0.33
  params$value[params$name == "aerenchyma" & params$type == "proportion"] <- temp$p_A
  params$value[params$name == "aerenchyma" & params$type == "n_files"] <- 10
  
  t_before <- Sys.time()
  # Run GRANAR ------------
  sim <- list()
  n_try <- 1
  while (TRUE) {
    did_it <- 1
    message("GRANAR in process")
    sim <- try(create_anatomy(parameters = params, verbatim=F), silent = TRUE)
    n_try <- n_try + 1
    if(is.null(sim)){
      sim <- try(error("nerror = ..."), silent = T)
      did_it <- 0
    }
    if(n_try > 100){
      sim <- list()
      sim$output$value <- NA
      did_it <- 0
    }
    if(!is(sim,'try-error')) break
  }
  if(did_it == 1){
    write_anatomy_xml(sim = sim,
                      path = "MECHA/cellsetdata/current_root.xml")
    fc <- file.copy(from = "MECHA/cellsetdata/current_root.xml",
                    to = paste0("MECHA/cellsetdata/",temp$RIL,"_",temp$Treatment,"_",temp$repet,".xml"),
                    overwrite = T)
    output <- output_as_RootScan(sim)
  }else{
    output[,1:12] <- NA
  }
  
  RS_Chimungu[id_sim, "model_area"] <- output$m_RXSA
  RS_Chimungu[id_sim, "model_TSA"] <- output$m_TSA
  RS_Chimungu[id_sim, "model_TSA_2"] <- output$m_TSA_2
  RS_Chimungu[id_sim, "model_TCA"] <- output$m_TCA
  RS_Chimungu[id_sim, "model_CCA"] <- output$m_CCA
  RS_Chimungu[id_sim, "model_XVA"] <- output$m_XVA
  RS_Chimungu[id_sim, "model_CC"] <- output$m_CC
  RS_Chimungu[id_sim, "model_AA"] <- output$m_AA
  RS_Chimungu[id_sim, "model_pA"] <- output$m_pA
  RS_Chimungu[id_sim, "sim_time_G"] <- output$granar_time
}


# MECHA loop on the generate cross-section -------------
k <- 1
for(id_sim in k:nrow(RS_Chimungu)){
  
  temp <- RS_Chimungu[id_sim,]
  message(paste0(">>>>>> ", k , " out of  ", nrow(RS_Chimungu),". In process: ",temp$RIL," ",temp$Treatment," ",temp$repet, " <<<<<<<"))
  
  fc <- try(file.copy(from = paste0("MECHA/cellsetdata/",temp$RIL,"_",temp$Treatment,"_",temp$repet,".xml"),
                      to = "MECHA/cellsetdata/current_root.xml",
                      overwrite = T), silent = T)
  if(!str_detect(fc[1], "Error :")){
    K <- MECHA_R(path = "MECHA/MECHAv4_light.py", Kx = F)
    
    fc <- file.copy(from = "MECHA/Projects/granar/out/GRANAR/Project_Test/baseline/Macro_prop_1,0.txt",
                    to = paste0("MECHA/cellsetdata/Chimungu/macro/Macro_prop_1,0_",id_sim,".txt"),
                    overwrite = T)
    fc <- file.copy(from = "MECHA/Projects/granar/out/GRANAR/Project_Test/baseline/Macro_prop_2,1.txt",
                    to = paste0("MECHA/cellsetdata/Chimungu/macro/Macro_prop_2,1_",id_sim,".txt"),
                    overwrite = T)
    fc <- file.copy(from = "MECHA/Projects/granar/out/GRANAR/Project_Test/baseline/Macro_prop_4,2.txt",
                    to = paste0("MECHA/cellsetdata/Chimungu/macro/Macro_prop_4,2_",id_sim,".txt"),
                    overwrite = T)
    
  }else{
    K <- c(NA,NA,NA)
  }
  
  RS_Chimungu[id_sim, "model_kr0"] <- K[1]
  RS_Chimungu[id_sim, "model_kr1"] <- K[2]
  RS_Chimungu[id_sim, "model_kr2"] <- K[3]
  
  k <- k + 1
}

write.csv(RS_Chimungu, "RS_Chimungu.csv")

# Once the script have run
# All data can be combine in a excel file

# Thne stored data can be ploted
Sum_up <- read_excel("Thesis/2019-01 GRANAR/Validation/Sum_up_03.xlsx",
                     col_types = c("text", rep("numeric", 22)))
Sum_up <- Sum_up%>%
  mutate(cortex_w = sqrt(RXSA/pi)-sqrt(TSA/pi),
         CC_size = model_cortex_w/round(CF),
         model_stelaire = model_TSA/model_area,
         stelaire = TSA/RXSA,
         cort_stele = model_TCA/model_TSA,
         model_kr0 = model_kr0*0.001157407,
         model_kr1 = model_kr1*0.001157407,
         model_kr2 = model_kr2*0.001157407)%>%
  filter(!is.na(model_area))

err <- Sum_up%>%
  transmute(n = na.omit(n()),
            RXSA_sd = sqrt(sum(na.omit(RXSA - model_area)^2)/(n-1)),
            RXSA_err = 100*sum(na.omit(abs(RXSA - model_area)/RXSA))/n,
            RXSA_r = 1- (sum(na.omit(model_area-RXSA)^2)/(sum(na.omit(model_area - mean(model_area))^2))),
            TSA_r = 1- (sum((model_TSA-TSA)^2)/(sum((model_TSA - mean(model_TSA))^2))),
            TSA_sd = sqrt(sum(na.omit(TSA - model_TSA)^2)/(n-1)),
            TSA_err = 100*sum(na.omit(abs(TSA - model_TSA)/TSA))/n,
            TCA_sd = sqrt(sum(na.omit(TCA - model_TCA)^2)/(n-1)),
            TCA_r = 1- (sum((model_TCA-TCA)^2)/(sum((model_TCA - mean(model_TCA))^2))),
            TCA_err = 100*sum(na.omit(abs(TCA - model_TCA)/TCA))/n,
            XVA_sd = sqrt(sum(na.omit(XVA - model_XVA)^2)/(n-1)),
            XVA_err = 100*sum(na.omit(abs(XVA - model_XVA)/XVA))/n,
            XVA_r = 1- (sum((model_XVA-XVA)^2)/(sum((model_XVA - mean(model_XVA))^2))),
            cortex_w_r = 1- (sum((model_cortex_w-cortex_w)^2)/(sum((model_cortex_w - mean(model_cortex_w))^2))),
            cortex_err = 100*sum(na.omit(abs(cortex_w - model_cortex_w)/cortex_w))/n,
            AA_sd = sqrt(sum(na.omit(AA - model_AA)^2)/(n-1)),
            AA_err = 100*sum(na.omit(abs(AA - model_AA)/AA))/n,
            AA_r = 1- (sum((model_AA-AA)^2)/(sum((model_AA - mean(model_AA))^2)))
            #CC_sd = sqrt(sum(na.omit(CC - model_CC)^2)/(n-1)),
            #CC_err = sum(na.omit(abs(CC - model_CC)/CC))/n
  )


Sum_up%>%
  ggplot()+
  geom_point(aes(stelaire, model_stelaire, colour = name, shape = name), alpha = 0.5, size = 2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2)+
  #geom_text(x= 0.6, y = 2.5, label = paste0("R²: ",round(err$RXSA_r,3)))+
  theme_cowplot()+
  labs(colour = "Data from:", 
       shape = "Data from:")+
  xlab("ratio stele area over the cross section area")+
  ylab("GRANAR ratio stele area over the cross section area")

Sum_up%>%
  ggplot()+
  geom_point(aes(RXSA, model_area, colour = name, shape = name), alpha = 0.5, size = 2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2)+
  #geom_text(x= 0.6, y = 2.5, label = paste0("R²: ",round(err$RXSA_r,3)))+
  theme_cowplot()+
  labs(colour = "Data from:", 
       shape = "Data from:")+
  xlab("Observed area of the cross-section (mm²)")+
  ylab("GRANAR total area (mm²)")

Sum_up%>%
  ggplot()+
  geom_point(aes(TSA, model_TSA, colour = name, shape = name), alpha = 0.5, size = 2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2)+
  geom_text(x= 0.2, y = 0.6, label = paste0("R²: ",round(err$TSA_r,3)))+
  theme_cowplot()+
  labs(colour = "Data from:", 
       shape = "Data from:")+
  xlab("Observed area of the stele (mm²)")+
  ylab("GRANAR stele area (mm²)")

Sum_up%>%
  ggplot()+
  geom_point(aes(TCA, model_TCA, colour = name, shape = name), alpha = 0.5, size = 2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2)+
  #geom_text(x= 0.5, y = 2, label = paste0("R²: ",round(err$TCA_r,3)))+
  theme_cowplot()+
  labs(colour = "Data from:", 
       shape = "Data from:")+
  xlab("Observed area of the total cortical area (mm²)")+
  ylab("GRANAR total cortical area (mm²)")
  


Sum_up%>%
  ggplot()+
  geom_point(aes(XVA, model_XVA, colour = name, shape = name), alpha = 0.5, size = 2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2)+
  #geom_text(x= 0.03, y = 0.13, label = paste0("R²: ",round(err$XVA_r,3)))+
  theme_cowplot()+
  labs(colour = "Data from:", 
       shape = "Data from:")+
  xlab("Observed area of the xylem area (mm²)")+
  ylab("GRANAR xylem area (mm²)")
  

Sum_up%>%
  ggplot()+
  geom_point(aes(AA, model_AA, colour = name, shape = name), alpha = 0.5, size = 2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2)+
  #geom_text(x= 0.1, y = 0.35, label = paste0("R²: ",round(err$AA_r,3)))+
  theme_cowplot()+
  labs(colour = "Data from:", 
       shape = "Data from:")+
  xlab("aerenchyma area (mm²)")+
  ylab("GRANAR aerenchyma area (mm²)")


Sum_up%>%
  mutate(cortex_w = sqrt(RXSA/pi)-sqrt(TSA/pi))%>%
  ggplot()+
  geom_point(aes(cortex_w, model_cortex_w, colour = name, shape = name), alpha = 0.5, size = 2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2)+
  theme_cowplot()+
  labs(colour = "Data from:", 
       shape = "Data from:")+
  xlab("Observed cortex width (mm)")+
  ylab("GRANAR cortex width (mm)")


Sum_up%>%
  ggplot()+
  geom_histogram(aes(model_kr0), fill = "seagreen", alpha = 0.8, bins = 60)+
  geom_histogram(aes(model_kr1), bins = 60, fill = "blue", alpha = 0.3)+
  geom_histogram(aes(model_kr2), fill = "red", alpha = 0.3, bins = 60)+
  xlab("radial conductivity (Kr)")+
  theme_cowplot()

cortex <- lm(model_kr0 ~ model_cortex_w, data = Sum_up)
summary(cortex)

Sum_up%>%
  ggplot(aes(model_cortex_w, model_kr0*10^8))+
  geom_point(aes(colour = model_XVA, shape = name), size = 4, alpha = 0.4)+
  scale_colour_viridis()+
  xlab("cortex width (mm)")+
  ylab(bquote('Simulated kr '*10^8~ (m~ MPa^-1~ s^-1)*''))+
  labs(colour = "Xylem area (mm²)", shape = "Data from:")+
  theme_cowplot()+
  guides(shape = F)


Sum_up%>%
  ggplot(aes(model_area, model_kr0*10^8))+
  #geom_smooth(aes(group = factor(round(model_XVA*20)/20)), method = "lm", linetype = 2)+
  geom_point(aes(colour = model_XVA, shape = name), size = 4, alpha = 0.3)+
  # geom_smooth(method = "lm", se = F, colour = "red", linetype = 1)+
  scale_colour_viridis()+
  xlab('Total cortical area (mm²)')+
  ylab(bquote('Simulated Kr '*10^8~ (m~ MPa^-1~ s^-1)*''))+
  labs(colour = "Xylem area (mm²)", shape = "Data from:")+
  theme_cowplot()+
  guides(shape = F)

Sum_up%>%
  ggplot(aes(model_area, model_kr0*10^5))+
  #geom_smooth(aes(group = factor(round(model_XVA*20)/20)), method = "lm", linetype = 2)+
  geom_point(aes(colour = model_XVA, shape = name), size = 4, alpha = 0.3)+
  # geom_smooth(method = "lm", se = F, colour = "red", linetype = 1)+
  scale_colour_viridis()+
  xlab('Total cortical area (mm²)')+
  ylab(bquote('Simulated Kr '*10^5~ (cm~ hPa^-1~ d^-1)*''))+
  labs(colour = "Xylem area (mm²)", shape = "Data from:")+
  theme_cowplot()+
  guides(shape = F)

Sum_up%>%
  ggplot(aes(model_area, model_kr0*10^5))+
  #geom_smooth(aes(group = factor(round(model_XVA*20)/20)), method = "lm", linetype = 2)+
  geom_point(aes(colour = model_XVA, shape = name), size = 4, alpha = 0.3)+
  # geom_smooth(method = "lm", se = F, colour = "red", linetype = 1)+
  scale_colour_viridis()+
  xlab('Total area (mm²)')+
  ylab(bquote('Simulated Kr '*10^5~ (cm~ hPa^-1~ d^-1)*''))+
  labs(colour = "Xylem area (mm²)", shape = "Data from:")+
  theme_cowplot()+
  guides(shape = F)

Sum_up%>%
  ggplot(aes(cort_stele, model_kr0*10^5))+
  #geom_smooth(aes(group = factor(round(model_XVA*20)/20)), method = "lm", linetype = 2)+
  geom_point(aes(colour = model_XVA, shape = name), size = 4, alpha = 0.3)+
  # geom_smooth(method = "lm", se = F, colour = "red", linetype = 1)+
  scale_colour_viridis()+
  xlab('Total cortical area (mm²)')+
  ylab(bquote('Simulated Kr '*10^5~ (cm~ hPa^-1~ d^-1)*''))+
  labs(colour = "Xylem area (mm²)", shape = "Data from:")+
  theme_cowplot()+
  guides(shape = F)

Sum_up%>%
  ggplot(aes(round(CF), model_kr0*10^5))+
  #geom_smooth(aes(group = factor(round(model_XVA*20)/20)), method = "lm", linetype = 2)+
  geom_point(aes(colour = model_cortex_w, shape = name), size = 4, alpha = 0.3)+
  # geom_smooth(method = "lm", se = F, colour = "red", linetype = 1)+
  scale_colour_viridis()+
  xlab('n cortical layer')+
  ylab(bquote('Simulated Kr '*10^5~ (cm~ hPa^-1~ d^-1)*''))+
  labs(colour = "cortex width (mm)", shape = "Data from:")+
  theme_cowplot()+
  guides(shape = F)

Sum_up%>%
  ggplot(aes(CC_size, model_kr0*10^5))+
  #geom_smooth(aes(group = factor(round(model_XVA*20)/20)), method = "lm", linetype = 2)+
  geom_point(aes(colour = model_cortex_w, shape = name), size = 4, alpha = 0.3)+
  # geom_smooth(method = "lm", se = F, colour = "red", linetype = 1)+
  scale_colour_viridis()+
  xlab('mean cortex cell diameter')+
  ylab(bquote('Simulated Kr '*10^5~ (cm~ hPa^-1~ d^-1)*''))+
  labs(colour = "cortex width (mm)", shape = "Data from:")+
  theme_cowplot()+
  guides(shape = F)

Sum_up%>%
  ggplot(aes(model_stelaire, model_kr0*10^8))+
  #geom_smooth(aes(group = factor(round(model_XVA*20)/20)), method = "lm", linetype = 2)+
  geom_point(aes(colour = model_XVA, shape = name), size = 4, alpha = 0.3)+
  # geom_smooth(aes(group = factor(n_X)),method = "lm", se = F, linetype = 1)+
  scale_colour_viridis()+
  xlab('ratio stele:CT area')+
  ylab(bquote('Simulated kr '*10^8~ (m~ MPa^-1~ s^-1)*''))+
  labs(colour = "Xylem area (mm²)", shape = "Data from:")+
  theme_cowplot()+
  guides(shape = F)

Sum_up%>%
  filter(!is.na(model_kr0))%>%
  ggplot(aes(model_pA, model_kr0*10^8))+
 # ggplot(aes(model_pA, norm_kr))+
  geom_point(aes(colour = model_cortex_w, shape = name), size = 5,alpha = 0.3)+
  scale_colour_viridis()+
  xlab("Aerenchyma proportion (%)")+
  ylab(bquote('Simulated Kr '*10^8~ (m~ MPa^-1~ s^-1)*''))+
  labs(colour = "cortex width (mm)", shape = "Data from:")+
  theme_cowplot()

Sum_up%>%
  filter(!is.na(model_kr0))%>%
  ggplot(aes(model_pA, model_kr0*10^8))+
  # ggplot(aes(model_pA, norm_kr))+
  geom_point(aes(shape = name, colour = model_cortex_w), size = 5,alpha = 0.3)+
  scale_colour_viridis()+
  xlab("Aerenchyma proportion (%)")+
  ylab(bquote('Simulated kr '*10^8~ (m~ MPa^-1~ s^-1)*''))+
  labs(colour = "cortex width (mm)", shape = "Data from:")+
  theme_cowplot()+
  theme(legend.position = "bottom")

Sum_up%>%
  filter(model_cortex_w > 0.269 & model_cortex_w < 0.274)%>%
  ggplot(aes(model_cortex_w, model_kr0*10^5))+
  geom_point(aes(shape = name, colour = model_CCA/model_CC), size = 5,alpha = 1)+
  scale_colour_viridis()+
  theme_cowplot()+
  theme(legend.position = "bottom")


for(i in 2:ncol(Sum_up)){
  print(colnames(Sum_up)[i])
  a <- Sum_up%>%
    transmute(Kr = model_kr0,
              v = c(Sum_up[[i]]))
  fit <- aov(Kr ~ v, data = a)
  print(summary(fit))
  message(paste0("RMSE value of the ",colnames(Sum_up)[i]," is : ", rmse(fit, a), " and cor :",  round(cor(a$Kr, a$v,  method = "pearson", use = "complete.obs"),3)))
  message(paste0("R sq value of the ",colnames(Sum_up)[i]," is : ",round(rsquare(fit, a),3) ))
}


cort <- Sum_up%>%
  filter(CC_size > mean(CC_size),
         CF < mean(round(CF)))
cort%>%
  ggplot(aes(CC_size, model_kr0*10^5))+
  #geom_smooth(aes(group = factor(round(model_XVA*20)/20)), method = "lm", linetype = 2)+
  geom_point(aes(colour = factor(round(CF)), shape = name), size = 4, alpha = 0.3)+
  # geom_smooth(method = "lm", se = F, colour = "red", linetype = 1)+
  xlab('mean cortex cell diameter')+
  ylab(bquote('Simulated Kr '*10^5~ (cm~ hPa^-1~ d^-1)*''))+
  labs(colour = "cortex width (mm)", shape = "Data from:")+
  theme_cowplot()+
  guides(shape = F)


Sum_up%>%
  mutate(rs = model_XVA/model_TSA)%>%
  ggplot(aes(model_stelaire, model_kr0*10^5))+
  #geom_smooth(aes(group = factor(round(model_XVA*20)/20)), method = "lm", linetype = 2)+
  geom_point(aes(colour = rs, shape = name), size = 4, alpha = 0.3)+
  # geom_smooth(method = "lm", se = F, colour = "red", linetype = 1)+
  scale_colour_viridis()+
  xlab('ratio')+
  ylab(bquote('Simulated Kr '*10^5~ (cm~ hPa^-1~ d^-1)*''))+
  labs(colour = "cortex width (mm)", shape = "Data from:")+
  theme_cowplot()+
  guides(shape = F)


ggplot()+
  geom_hline(yintercept = c(27, 2.2), linetype = 2, colour = "red")+ # Steudle 1993
  geom_hline(yintercept = c(22), linetype = 1, colour = "red")+ # Frensch
  geom_hline(yintercept = c(10.5,21.8), linetype = 2, colour = "blue")+ # Schambil
  geom_hline(yintercept = c(10,20), linetype = 2, colour = "green")+ # Miller 1985
  geom_hline(yintercept = c(10.3), linetype = 1, colour = "green")+ # Shone
  geom_hline(yintercept = c(7), linetype = 1, colour = "blue")+ # Ye
  geom_hline(yintercept = c(12,9), linetype = 2, colour = "grey")+ # Steaudle 1987
  geom_hline(yintercept = c(22, 1.85), linetype = 4, colour = "red")+ # Doussan
  geom_hline(yintercept = c(8.1, 21.1), linetype = 3, size = 2, colour = "magenta")




  

