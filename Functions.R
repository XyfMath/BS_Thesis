### Functions used

# Normalize Data
# Same as the normalization function in ALRA, without removing lines where all elements
# are zeros.

N_Data <- function(A){
  M <- rowSums(A)
  B <- A
  for (i in 1:(nrow(A))) {
    if(M[i] != 0){
      for(j in 1:(ncol(A))){
        B[i,j] <- (A[i,j] / M[i])
      }
    }
  }
  B <- B * 10E3
  B <- log(B +1)
  return(B)
}

# Calculate percentages of Biological Zeros and Technical Zeros
# A: True expression matrix, B: Expression matrix with dropout

Count_B_T <- function(A, B){
  n <- 0
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      if (A[i,j] == 0){
        n <- (n+1)
      }
    }
  }
  m <- 0
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      if ((A[i,j] != 0) && (B[i,j] == 0)){
        m <- (m+1)
      }
    }
  }
  bzr <- n / (nrow(A) * ncol(A))
  tzr <- m / (nrow(A) * ncol(A))
  return (matrix(c("biological zero rate", bzr, "techical zero rate",tzr),nrow=2))
}

# Calculate Imputation Parameters i.e. biological zeros presevation rate
# and technical zeros completion rate

Count_I <- function(A, B, C){
  a <- 0
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      if ((near(A[i,j], 0) == TRUE)){
        a <- (a+1)
      }
    }
  }
  b <- 0
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      if ((near(A[i,j], 0) == FALSE) && (near(B[i,j], 0) == TRUE)){
        b <- (b+1)
      }
    }
  }
  n <- 0
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      if ((near(A[i,j], 0) == TRUE) && (near(C[i,j], 0) == TRUE)){
        n <- (n+1)
      }
    }
  }
  p <- 0
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      if ((near(A[i,j], 0) == FALSE) && (near(B[i,j], 0) == TRUE) && (near(C[i,j], 0) == FALSE)){
        p <- (p+1)
      }
    }
  }
  bpr <- n / a
  bmr <- (1 - bpr)
  tcr <- p / b
  tmr <- (1 - tcr)
  return (matrix(c("biological preservation rate", bpr, "biological mistake rate", bmr, "technical complete rate",tcr, "technical mistake rate",tmr),nrow=2))
}

# Plot a expression matrix

Plot_EM <- function(A,a,b){
  colnames(A) <- NULL
  rownames(A) <- NULL
  melted_A <- melt(A)
  head(melted_A)
  ggplot(data = melted_A, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + 
    scale_fill_gradient(low = "white", high = "red", limits=c(a,b))
}


# Calculate sample correlation coefficient of imputed values of zeros in the expression matrix with dropout
# and their true expression values
# A: True expression matrix, B: Expression matrix with dropout, C: Imputed expression matrix
SCCz <- function(A,B,C){
  n <- 0
  aveI <- 0
  aveT <- 0
  cov <- 0
  varA <- 0
  varC <- 0
  s <- 0
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      if (B[i,j] == 0){
        n <- (n+1)
        aveI <- aveI + C[i,j]
        aveT <- aveT + A[i,j]
      }
    }
  }
  aveI <- aveI/n
  aveT <- aveT/n
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      if (B[i,j] == 0){
        cov <- cov + (A[i,j]- aveT)*(C[i,j]- aveI)
        varA <- varA + (A[i,j]- aveT)^2
        varC <- varC + (C[i,j]- aveI)^2
      }
    }
  }
  s <- cov/sqrt(varA*varC)
  s <- as.numeric(s)
  print("Calculate sample correlation coefficient -- zeros: ")
  return(s)
}

# Calculate sample correlation coefficient of imputed non-dropout values in the expression matrix with dropout
# and their true expression values
# A: True expression matrix, B: Expression matrix with dropout, C: Imputed expression matrix
SCCn <- function(A,B,C){
  n <- 0
  aveI <- 0
  aveT <- 0
  cov <- 0
  varA <- 0
  varC <- 0
  s <- 0
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      if ((A[i,j] == 0)||(B[i,j] != 0)){
        n <- (n+1)
        aveI <- aveI + C[i,j]
        aveT <- aveT + A[i,j]
      }
    }
  }
  aveI <- aveI/n
  aveT <- aveT/n
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      if ((A[i,j] == 0)||(B[i,j] != 0)){
        cov <- cov + (A[i,j]- aveT)*(C[i,j]- aveI)
        varA <- varA + (A[i,j]- aveT)^2
        varC <- varC + (C[i,j]- aveI)^2
      }
    }
  }
  s <- cov/sqrt(varA*varC)
  s <- as.numeric(s)
  print("Calculate sample correlation coefficient -- non-dropout: ")
  return(s)
}

# Calculate sample correlation coefficient of imputed values in the expression matrix with dropout
# and their true expression values
# A: True expression matrix, B: Expression matrix with dropout, C: Imputed expression matrix
SCCal <- function(A,B,C){
  n <- 0
  aveI <- 0
  aveT <- 0
  cov <- 0
  varA <- 0
  varC <- 0
  s <- 0
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
        n <- (n+1)
        aveI <- aveI + C[i,j]
        aveT <- aveT + A[i,j]
      
    }
  }
  aveI <- aveI/n
  aveT <- aveT/n
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
        cov <- cov + (A[i,j]- aveT)*(C[i,j]- aveI)
        varA <- varA + (A[i,j]- aveT)^2
        varC <- varC + (C[i,j]- aveI)^2
      
    }
  }
  s <- cov/sqrt(varA*varC)
  s <- as.numeric(s)
  print("Calculate sample correlation coefficient -- all: ")
  return(s)
}

###