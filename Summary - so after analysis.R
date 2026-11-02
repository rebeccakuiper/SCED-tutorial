################################################################################
# Summary of the person-specific GORICA results (run after the ENET_Pk scripts)
#
# Works with restriktor 0.6-x (CRAN) and 0.7+: columns are selected by name,
# and the ratio matrix is matched to the hypotheses by position (row/column
# names differ between versions, e.g. 0.7 adds ' (best)' to a row name).
################################################################################

# Participants (objects output_P1, output_P3, ... must exist in the workspace)
P_names <- c("P1", "P3", "P4", "P5", "P6", "P7", "P9", "P10", "P11")

outputs <- mget(paste0("output_", P_names), envir = globalenv(), ifnotfound = list(NULL))
missing_P <- P_names[vapply(outputs, is.null, logical(1))]
if (length(missing_P) > 0) {
  stop("Run the analysis scripts first; not found: ",
       paste0("output_", missing_P, collapse = ", "))
}
names(outputs) <- P_names


## Helper: one tidy table per person ----------------------------------------------
goric_table <- function(out, participant) {
  res  <- as.data.frame(out$result)
  cols <- c("model", "loglik", "penalty", "gorica", "gorica.weights")
  if (!all(cols %in% names(res))) {
    stop("Unexpected goric output for ", participant, "; columns found: ",
         paste(names(res), collapse = ", "))
  }
  # Ratios of GORICA weights (row hypothesis versus column hypothesis);
  # rows and columns are in the same order as the hypotheses in 'result'.
  ratios <- as.matrix(out$ratio.gw)
  if (!all(dim(ratios) == nrow(res))) stop("Unexpected ratio.gw for ", participant)
  ratios <- as.data.frame(ratios)
  names(ratios) <- paste("vs.", res$model)

  data.frame(Participant = participant,
             Hypothesis  = res$model,
             res[, cols[-1]],
             ratios,
             row.names = NULL, check.names = FALSE)
}

mydata <- do.call(rbind, Map(goric_table, outputs, P_names))
rownames(mydata) <- NULL
num_cols   <- vapply(mydata, is.numeric, logical(1))
ratio_cols <- grepl("^vs\\. ", names(mydata))
mydata[num_cols] <- round(mydata[num_cols], 3)


## Person-level results -------------------------------------------------------------
for (p in P_names) {
  cat("\n==", p, "==\n")
  print(mydata[mydata$Participant == p, -1], row.names = FALSE)
}


## Preferred hypothesis per person ---------------------------------------------------
# Preferred: has lowest GORICA value (i.e., highest GORICA weight).
# An informative hypothesis is 'weak' if it is not supported more than the
# unconstrained hypothesis.
overview <- do.call(rbind, lapply(P_names, function(p) {
  d     <- mydata[mydata$Participant == p, ]
  unc_w <- d$gorica.weights[d$Hypothesis == "unconstrained"]
  inf   <- d[d$Hypothesis != "unconstrained", ]
  data.frame(Participant = p,
             Preferred   = d$Hypothesis[which.min(d$gorica)],
             Weight      = max(d$gorica.weights),
             Weak        = paste(inf$Hypothesis[inf$gorica.weights < unc_w], collapse = ", "))
}))
overview


## Save -------------------------------------------------------------------------------
# Tab-separated text with a decimal comma 
# (opens in Excel; Excel may warn that the format does not match the .xls extension).
write.table(mydata, "Results.xls", row.names = FALSE, sep = "\t", dec = ",")

# Same, with ratios above 1000 shown as '>1000' (as in Table 4 of the paper)
mydata_1000 <- mydata
mydata_1000[num_cols] <- lapply(mydata[num_cols], formatC, format = "f", digits = 3,
                                decimal.mark = ",")
for (j in which(ratio_cols)) {
  mydata_1000[mydata[[j]] > 1000, j] <- ">1000"
}
mydata_1000
write.table(mydata_1000, "Results_1000.xls", row.names = FALSE, sep = "\t", quote = FALSE)

######################
