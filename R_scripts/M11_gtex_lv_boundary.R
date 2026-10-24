# M11 · GTEx v8 left-ventricle expression boundary of the 332 LRGs
# Added in v9.2 (response to pre-submission revision R1). Output = Supplementary Table S8.
# Input: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
#   Download: https://gtexportal.org/home/downloads/adult-gtex/bulk-gex ("GTEx v8 median TPM by tissue")
# Output: S8_GTEx_LV_expression.csv (332 genes x median TPM / detectable / in-model)
# Expected values (from the reference run; see manuscript Suppl. Table S8):
#   mapped 307/332 (92.5%) | TPM>=1 285 (85.8%) | median TPM 18.56 (IQR 6.50-44.43)
#   13 model genes all detectable, TPM 4.76 (JMJD1C) - 103.44 (DDX17)
suppressMessages(library(data.table))

gct_path <- 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz'
stopifnot(file.exists(gct_path))

gct <- fread(gct_path, skip = 2)
lv_col <- grep('Heart - Left Ventricle', names(gct), value = TRUE)
stopifnot(length(lv_col) == 1)
message('LV column: ', lv_col)

lrg <- readLines('lrg332.txt')            # in data/; keep it next to the script or adjust the path
lrg <- trimws(lrg[lrg != ''])

# H2AZ1 is queried under its current symbol H2AFZ (renamed in GTEx v8)
lrg_query <- lrg
lrg_query[lrg_query == 'H2AZ1'] <- 'H2AFZ'

res <- data.frame(Gene = lrg, TPM = NA_real_, stringsAsFactors = FALSE)
sym <- gct[['Description']]; names(sym) <- gct[['Name']]
tpm <- gct[[lv_col]]; names(tpm) <- sym
for (i in seq_along(lrg)) {
  q <- lrg_query[i]
  res$TPM[i] <- if (q %in% names(tpm)) tpm[[q]] else NA_real_
}

res$Detectable <- ifelse(is.na(res$TPM), 'not mapped', ifelse(res$TPM >= 1, 'Yes', 'No'))
write.csv(res, 'S8_GTEx_LV_expression.csv', row.names = FALSE)

# ---- Self-verification ----
mapped   <- sum(!is.na(res$TPM))
detect   <- sum(res$TPM >= 1, na.rm = TRUE)
cat(sprintf('mapped %d/%d (%.1f%%)\n', mapped, length(lrg), 100*mapped/length(lrg)))
cat(sprintf('TPM>=1 %d (%.1f%%)\n', detect, 100*detect/length(lrg)))
cat(sprintf('median TPM %.2f (IQR %.2f-%.2f)\n',
    median(res$TPM, na.rm=TRUE), quantile(res$TPM, .25, na.rm=TRUE), quantile(res$TPM, .75, na.rm=TRUE)))
cat('JMJD1C:', res$TPM[res$Gene=='JMJD1C'], ' DDX17:', res$TPM[res$Gene=='DDX17'], '\n')
cat('== Expected: 307 (92.5%) / 285 (85.8%) / 18.56 / JMJD1C 4.76 / DDX17 103.44 ==\n')
