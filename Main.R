source('alra.R')
source('Functions.R')
library(tidyverse)
library(rsvd)
library(ggplot2)
library(scImpute)
library(splatter)
library(scater)
library(reshape2)



# Set numbers of Genes and Cells

num_genes <- 1000
num_cells <- 300

# Set Parameters (without dropout)

  params <- newSplatParams()
  params <- setParam(params, "nGenes", num_genes)
  params <- setParam(params, "batchCells", num_cells)
  params <- setParam(params, "mean.shape", 0.11)
  params <- setParam(params, "seed", 200000)
  params
  
# Simulation (without dropout)
  
sim1 <- splatSimulate(params)
  
# Set Dropout Parameters
  
params <- setParam(params, "dropout.shape", -1)
params <- setParam(params, "dropout.mid", 5.5)
params <- setParam(params, "dropout.type", "experiment")
  
# Simulation (with dropout)
  
sim1_D <- splatSimulate(params)
  
# Caculate Pecentages of Biological Zeros and Technical Zeros
  
sim1_C <- counts(sim1)
sim1_D_C <- counts(sim1_D)
Count_B_T(sim1_C, sim1_D_C)
  
# ALRA Imputation
  
sim1_D_C_norm <- N_Data(sim1_D_C)
k_choice <- choose_k((sim1_D_C_norm))
sim1_D_C_norm_completed <- alra(sim1_D_C_norm,k=k_choice$k)[[3]]
sim1_C_norm <- N_Data(sim1_C)
print("ALRA")
Count_I(sim1_C_norm, sim1_D_C_norm, sim1_D_C_norm_completed)
SCCz(sim1_C_norm, sim1_D_C_norm, sim1_D_C_norm_completed)
SCCn(sim1_C_norm, sim1_D_C_norm, sim1_D_C_norm_completed)
SCCal(sim1_C_norm, sim1_D_C_norm, sim1_D_C_norm_completed)
  
# ScImpute Imputation
# Please customize file.path if used on another machine
  
write.csv(sim1_D_C, 'sim1_D_C.csv')
scimpute(# full path to raw count matrix
  count_path = file.path("C:", "Users", "fate_", "Desktop", "中国科学院","数学 毕业设计","test","sim1_D_C.csv"), 
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = "./",           # full path to output directory
  labeled = FALSE,          # cell type labels not available
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 2,             # 2 cell subpopulations
  ncores = 1)              # number of cores used in parallel computation
  
sim1_scImpute <- read.csv("scimpute_count.csv")
sim1_scImpute <- data.matrix(sim1_scImpute)
sim1_scImpute <- sim1_scImpute[,-1]
sim1_scImpute_norm <- N_Data(sim1_scImpute)
print("scImpute")
Count_I(sim1_C_norm, sim1_D_C_norm, sim1_scImpute_norm)
SCCz(sim1_C_norm, sim1_D_C_norm, sim1_scImpute_norm)
SCCn(sim1_C_norm, sim1_D_C_norm, sim1_scImpute_norm)
SCCal(sim1_C_norm, sim1_D_C_norm, sim1_scImpute_norm)

# Plot

#Plot_EM(sim1_C_norm,0,5.5)
#Plot_EM(sim1_D_C_norm,0,5.5)
#Plot_EM(sim1_D_C_norm_completed,0,5.5)
#Plot_EM(sim1_scImpute_norm,0,5.5)

