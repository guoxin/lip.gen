load("reports.Rdata")
load("fixAlleleData")
load("BLOSUM62.2.Rdata")
source("internals.R")
for(nm in names(reports)){
    predict <- reports[[nm]][["predict"]]
    predict.matrix <- divisor <- matrix(0, ncol = length(predict), nrow = length(predict))

    predict.matrix <- sweep(predict.matrix, MARGIN = 1, STATS = predict, FUN = "+")
    predict.matrix <- sweep(predict.matrix, MARGIN = 2, STATS = predict, FUN = "-")
    predict.matrix <- abs(predict.matrix)

    divisor <- stringToK3(string.set = fixAlleleData[[nm]][["peptides"]], K1 = BLOSUM62.2, beta = 0.11387)
    divisor <- sqrt(2 - 2 * divisor)
    diag(divisor <- 1)
    quotient <- predict.matrix / divisor
    diag(quotient) <- 0
    gpgx.time.stamp("for Allele", nm, "lip value is ", max(quotient))
}
