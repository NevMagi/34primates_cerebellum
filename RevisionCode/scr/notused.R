# Function to apply random transformations to a dataframe
transform_df <- function(df) {
  # Randomly transform MedianCerebellum and MedianCerebrum columns
  transformed_df <- df
  transformed_df$MedianCerebellum <- transformed_df$MedianCerebellum * runif(1, 0.5, 2.0)
  transformed_df$MedianCerebrum <- transformed_df$MedianCerebrum * runif(1, 0.5, 2.0)
  return(transformed_df)
}

# Function to apply random transformations to a dataframe
# Constrained version that assumes relatively minute differential shrinkage across tissues.
transform_df_constraint <- function(df) {
  # Generate random transformation factors for MedianCerebellum and MedianCerebrum columns
  transformation_cerebellum <- runif(nrow(df), 0.5, 2.0)
  transformation_cerebrum <- runif(nrow(df), 0.5, 2.0)
  
  # Generate random ratios between MedianCerebellum and MedianCerebrum columns
  ratio <- runif(nrow(df), 0.91, 1.1) ## Constrained to apply transformations on the MedianCerebellum and MedianCerebrum colums so that the ratio of transformations per row is 0.91 < x > 1.1.
  
  # Apply transformations to MedianCerebellum and MedianCerebrum columns using the random ratios and transformation factors
  transformed_df <- df
  transformed_df$MedianCerebellum <- transformed_df$MedianCerebellum * (transformation_cerebellum / ratio)
  transformed_df$MedianCerebrum <- transformed_df$MedianCerebrum * (transformation_cerebrum * ratio)
  
  return(transformed_df)
}


# Generate 1000 transformed dataframes and compute PGLS and PIC for each
bm <- corBrownian(1, tree)
for (i in 1:10) {
  # Apply random transformations to the original dataframe
  transformed_df <- transform_df(phenoMedian)
  
  # Compute PGLS under Brownian Motion
  pgls_model <- gls(MedianCerebellum ~ MedianCerebrum, data=transformed_df, correlation=bm)
  pgls_summary <- summary(pgls_model)
  pgls_coefs <- coef(pgls_model)
  
  # Compute PIC regression
  pic_model <- lm(pic.C.med ~ pic.V.med + 0, data=transformed_df)
  pic_summary <- summary(pic_model)
  pic_coefs <- coef(pic_model)
  
  # Save results to file
  sink("3.Robustness_variability/3.1.PGLS-PIC-C-V-1000.txt", append=TRUE)
  cat(paste0("Randomly transformed dataframe ", i, ":\n"))
  cat(paste0("PGLS summary:\n", capture.output(pgls_summary), "\n"))
  cat(paste0("PGLS coefficients:\n", pgls_coefs, "\n"))
  cat("---------------------------\n")
  cat(paste0("PIC summary:\n", capture.output(pic_summary), "\n"))
  cat(paste0("PIC coefficients:\n", pic_coefs, "\n\n"))
  sink()
}



### STORING ALL RESULTS IN A SINGLE CSV ###

library(nlme)

# set number of transformations
n_transformations <- 10

# create a data frame to store results
results_df <- data.frame(matrix(ncol = 4, nrow = n_transformations))
colnames(results_df) <- c("pgls_summary", "pgls_coefs", "pic_summary", "pic_coefs")

# loop through transformations
for(i in 1:n_transformations) {
  # make a copy of original dataframe
  phenoMedian_transformed <- phenoMedian
  
  # create a random multiplier between 0.5 and 2.0
  multiplier <- runif(1, 0.5, 2)
  
  # create a random ratio between 0.91 and 1.1 for MedianCerebellum and MedianCerebrum
  ratio <- runif(1, 0.91, 1.1)
  
  # apply random transformation to MedianCerebellum column
  phenoMedian_transformed$MedianCerebellum <- phenoMedian_transformed$MedianCerebellum * ratio * multiplier
  
  # apply random transformation to MedianCerebrum column
  phenoMedian_transformed$MedianCerebrum <- phenoMedian_transformed$MedianCerebrum * (1/ratio) * multiplier
  
  # perform PGLS under Brownian Motion
  bm <- corBrownian(1, tree)
  pglsModel <- gls(MedianCerebellum ~ MedianCerebrum, data = phenoMedian_transformed, correlation = bm)
  
  # extract PGLS summary and coefficients
  pgls_summary <- capture.output(summary(pglsModel))
  pgls_coefs <- as.numeric(coef(pglsModel))
  
  # perform PIC regression forced through 0
  PICmodel <- lm(pic.C.med ~ pic.V.med + 0, data = phenoMedian_transformed)
  
  # extract PIC summary and coefficients
  pic_summary <- capture.output(summary(PICmodel))
  pic_coefs <- as.numeric(coef(PICmodel))
  
  # store results in data frame
  results_df[i, 1] <- paste(pgls_summary, collapse = "\n")
  results_df[i, 2] <- paste(pgls_coefs, collapse = ",")
  results_df[i, 3] <- paste(pic_summary, collapse = "\n")
  results_df[i, 4] <- paste(pic_coefs, collapse = ",")
  
}

# write results to CSV
write.csv(results_df, "3.Robustness_variability/3.2.PGLS-PIC-C-V-1000.csv", row.names = FALSE)