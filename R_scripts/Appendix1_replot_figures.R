# ============================================================
# Appendix 1: replot three figures (Figure 4B / 5C / 7, base versions)
# ============================================================
# [What it does]
#   Figure 4B: DLRG cross-module enrichment OR bar chart (DLRG_module_enrichment.csv)
#   Figure 5C: GSEA-KEGG bubble plot, Top 15 (Table5_GSEA_KEGG_results.csv)
#   Figure 7 : MR forest plot, base version, 26 genes (MR_results/Table3_MR_main_results.csv)
#   Each figure is written as PNG (screen preview) + TIFF (300 dpi, LZW,
#   submission-ready).
#
# [Expected values (automatically compared after the run; on mismatch check
#  the console output)]
#   Fig4B: module rows = 18 (grey is excluded from enrichment); green module
#          OR=7.47, P=6.5e-07
#   Fig5C: significant KEGG pathways (p.adjust<0.05) = 80; the plot takes the
#          Top 15 by |NES|
#   Fig7 : testable genes = 26; nominally significant (red) = 5
#          (SOD1/H2AZ1/ACAT2/LDHB/EIF4G1); all q>=0.132 (min 0.132)
#
# [Dependencies] ggplot2 only (installed by M0). If missing:
#   install.packages("ggplot2")   # or use a CRAN mirror near you
#
# [Paths] edit only the DIR line below. Use forward slashes on Windows
#   (e.g. "./GEO_data"; backslashes cause escape errors).
# ============================================================

DIR <- "./GEO_data"   # set to the folder holding the result CSVs

F_FIG4B <- file.path(DIR, "DLRG_module_enrichment.csv")
F_FIG5C <- file.path(DIR, "Table5_GSEA_KEGG_results.csv")
F_FIG7  <- file.path(DIR, "MR_results", "Table3_MR_main_results.csv")
OUTDIR  <- file.path(DIR, "replot_output")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

suppressMessages(library(ggplot2))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUTDIR, paste0(name, ".png")), p, width=w, height=h, dpi=300)
  ggsave(file.path(OUTDIR, paste0(name, ".tiff")), p, width=w, height=h,
         dpi=300, compression="lzw")
  cat("[saved]", name, ".png/.tiff ->", OUTDIR, "\n")
}

# ============================================================
# Part 1: Figure 4B — DLRG cross-module enrichment OR bar chart
# ============================================================
enr <- read.csv(F_FIG4B, check.names=FALSE, fileEncoding="UTF-8-BOM")
cat("[Fig4B] module rows =", nrow(enr), "(expected 18)\n")
g <- enr[enr$module=="green", ]
cat("[Fig4B] green: OR =", g$odds_ratio, ", P =", g$p_value, "(expected 7.47 / 6.5e-07)\n")

enr$log2OR <- log2(as.numeric(enr$odds_ratio))
enr$module <- factor(enr$module, levels=enr$module[order(enr$log2OR)])
# Fill palette: green module = green; other p<0.05 = dark blue; ns = light grey
enr$fill <- ifelse(enr$module=="green", "#2e7d32",
            ifelse(as.numeric(enr$p_value)<0.05, "#1e3a8a", "#c3cbd6"))
# significance stars
enr$star <- ifelse(as.numeric(enr$p_value)<0.05, "*", "")
# OR value labels (only OR>=2, to avoid crowding)
enr$orlab <- ifelse(as.numeric(enr$odds_ratio)>=2, sprintf("%.1f", as.numeric(enr$odds_ratio)), "")

p4b <- ggplot(enr, aes(x=module, y=log2OR, fill=fill)) +
  geom_col(width=0.72, color="grey25", linewidth=0.25) +
  geom_text(aes(label=paste0(orlab, star)),
            hjust=ifelse(enr$log2OR>=0, -0.25, 1.25), size=3.2, color="grey15") +
  geom_hline(yintercept=0, linewidth=0.4) +
  scale_fill_identity() +
  labs(x=NULL, y="log2(Odds Ratio)",
       title="DLRG enrichment across WGCNA modules") +
  coord_flip() +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold", size=12),
        axis.text=element_text(color="black"))
save_fig(p4b, "Figure4B_DLRG_module_enrichment_OR", 6.5, 5.5)

# ============================================================
# Part 2: Figure 5C — GSEA-KEGG bubble plot
# ============================================================
kk <- read.csv(F_FIG5C, check.names=FALSE, fileEncoding="UTF-8-BOM")
cat("[Fig5C] total rows =", nrow(kk), "\n")
sig <- kk[!is.na(as.numeric(kk$p.adjust)) & as.numeric(kk$p.adjust)<0.05, ]
cat("[Fig5C] significant pathways =", nrow(sig), "(expected 80)\n")

top <- sig[order(-abs(as.numeric(sig$NES))), ][1:15, ]
top$Description <- factor(top$Description, levels=top$Description[order(as.numeric(top$NES))])
top$size <- as.numeric(top$setSize)

p5c <- ggplot(top, aes(x=NES, y=Description)) +
  geom_point(aes(size=size, fill=NES), shape=21, color="grey30", stroke=0.3) +
  geom_vline(xintercept=0, linetype=2, linewidth=0.35, color="grey55") +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0,
                       name="NES") +
  scale_size_continuous(range=c(2.2, 6), name="Gene set size") +
  labs(x="Normalized Enrichment Score (NES)", y=NULL,
       title="GSEA-KEGG: HF vs NF (Top 15 by |NES|)") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold", size=12),
        axis.text=element_text(color="black"),
        legend.box="vertical", legend.key.size=grid::unit(0.8,"lines"))
save_fig(p5c, "Figure5C_GSEA_KEGG_bubble", 8, 6)

# ============================================================
# Part 3: Figure 7 — MR forest plot, base version (26 genes)
# ============================================================
# Note: this is a clean base version; fonts/palette/spacing are left free to
# customize. Data come straight from Table 3 and need no edits.
# ============================================================
mr <- read.csv(F_FIG7, check.names=FALSE, fileEncoding="UTF-8-BOM")
ok <- mr[mr$Status=="OK", ]
cat("[Fig7] testable genes =", nrow(ok), "(expected 26)\n")
red <- ok[as.numeric(ok$P)<0.05, ]$Gene
cat("[Fig7] nominally significant (red) =", length(red), ":", paste(red, collapse="/"), "\n")
cat("[Fig7] min q =", min(as.numeric(ok$q_BH)), "(expected 0.132)\n")

ci <- do.call(rbind, strsplit(gsub("\\[|\\]| ", "", ok$CI95), ","))
ok$beta <- as.numeric(ok$beta); ok$lo <- as.numeric(ci[,1]); ok$hi <- as.numeric(ci[,2])
ok <- ok[order(as.numeric(ok$P)), ]
ok$Gene <- factor(ok$Gene, levels=rev(ok$Gene))
ok$col <- ifelse(as.numeric(ok$P)<0.05, "#c0392b", "grey35")
ok$qlab <- sprintf("%.2f", as.numeric(ok$q_BH))

xr <- max(abs(c(ok$lo, ok$hi))) * 1.18
p7 <- ggplot(ok, aes(x=beta, y=Gene)) +
  geom_vline(xintercept=0, linetype=2, linewidth=0.4, color="grey55") +
  geom_pointrange(aes(xmin=lo, xmax=hi, color=col),
                  size=0.55, linewidth=0.5, fatten=2.2, show.legend=FALSE) +
  geom_text(aes(x=xr*0.99, label=qlab), hjust=1, size=2.8, color="grey30") +
  scale_color_identity() +
  scale_x_continuous(limits=c(-xr, xr), breaks=scales::pretty_breaks(n=7)) +
  labs(x="log(OR) per SD increase in genetically predicted expression  (95% CI)",
       y=NULL, title="cis-eQTL MR: 26 testable genes vs heart failure") +
  annotate("text", x=xr*0.99, y=Inf, label="q (BH)", vjust=1.8, hjust=1,
           size=2.8, color="grey30", fontface="bold") +
  theme_classic(base_size=10.5) +
  theme(plot.title=element_text(face="bold", size=11.5),
        axis.text=element_text(color="black"),
        axis.text.y=element_text(size=9))
save_fig(p7, "Figure7_MR_forest_base", 6.8, 8.2)

cat("\n========== Appendix 1: all three figures done ==========\n")
cat("output directory:", OUTDIR, "\n")
cat("verify the console anchor lines (18 module rows / green 7.47 / 80 KEGG terms / 26 genes / 5 red / q 0.132)\n")
