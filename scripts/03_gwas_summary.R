#!/usr/bin/env Rscript
# Genomic inflation (lambda), QQ and Manhattan plots for fastGWA results.
#
# Usage: Rscript 03_gwas_summary.R [gwas_dir]

args     <- commandArgs(trailingOnly = TRUE)
gwas_dir <- if (length(args) >= 1) args[1] else "gwas"
files    <- list.files(gwas_dir, pattern = "\\.fastGWA$", full.names = TRUE)
if (length(files) == 0) stop("No .fastGWA files found in ", gwas_dir)

for (f in files) {
  label <- sub("\\.fastGWA$", "", basename(f))
  res   <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
  res   <- res[!is.na(res$P), ]

  lambda <- median(qchisq(res$P, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
  n_sig  <- sum(res$P < 5e-8)
  message(sprintf("%-8s SNPs = %d, lambda = %.3f, P < 5e-8: %d",
                  label, nrow(res), lambda, n_sig))

  png(file.path(gwas_dir, paste0(label, "_qq.png")), width = 1600, height = 1600, res = 300)
  obs <- -log10(sort(res$P))
  exp <- -log10(ppoints(length(obs)))
  plot(exp, obs, pch = 20, cex = 0.4,
       xlab = expression(Expected ~ -log[10](P)), ylab = expression(Observed ~ -log[10](P)),
       main = sprintf("%s (lambda = %.3f)", label, lambda))
  abline(0, 1, col = "red")
  dev.off()

  res  <- res[order(res$CHR, res$POS), ]
  ends <- tapply(res$POS, res$CHR, max)
  offs <- setNames(c(0, cumsum(as.numeric(ends))[-length(ends)]), names(ends))
  x    <- res$POS + offs[as.character(res$CHR)]
  mids <- offs + ends / 2

  png(file.path(gwas_dir, paste0(label, "_manhattan.png")), width = 3000, height = 1400, res = 300)
  plot(x, -log10(res$P), pch = 20, cex = 0.4, xaxt = "n",
       col = c("grey30", "steelblue")[res$CHR %% 2 + 1],
       xlab = "Chromosome", ylab = expression(-log[10](P)), main = label)
  axis(1, at = mids, labels = names(ends), cex.axis = 0.6)
  abline(h = -log10(5e-8), col = "red", lty = 2)
  dev.off()
}
