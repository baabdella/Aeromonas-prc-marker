if (!require("phytools")) install.packages("phytools")
library(ape)
library(phytools)

# 1. Load and Clean
t1 <- read.tree("Reference_tree.nwk")
t2 <- read.tree("MLST_tree.nwk")

t1$tip.label <- gsub("[._-]", " ", t1$tip.label)
t2$tip.label <- gsub("[._-]", " ", t2$tip.label)

common <- intersect(t1$tip.label, t2$tip.label)
p1 <- keep.tip(t1, common)
p2 <- keep.tip(t2, common)

# 2. Global Congruence Stats
rf_val <- as.numeric(dist.topo(unroot(multi2di(p1)), unroot(multi2di(p2)), method="PH85"))
max_rf <- 2 * (length(common) - 3)
score <- 1 - (rf_val / max_rf)

# 3. Stability & Color Assignment
d1 <- cophenetic(p1)
d2 <- cophenetic(p2)[rownames(d1), colnames(d1)]
stability <- sapply(rownames(d1), function(x) cor(d1[x,], d2[x,]))

# Define threshold and colors
threshold <- 0.85
unstable_strains <- names(stability[stability < threshold])

# Create color map for the plot
edge_cols <- rep(make.transparent("#2E86C1", 0.4), length(common)) # Professional Blue
names(edge_cols) <- common
edge_cols[unstable_strains] <- make.transparent("#E74C3C", 0.8)    # Sharp Red

# 4. Save Stability Scores to CSV
# Adding a "Status" column to make the CSV easy to filter
stability_df <- data.frame(
  Strain = names(stability),
  Stability_Score = round(as.numeric(stability), 4),
  Status = ifelse(stability < threshold, "Unstable (Red)", "Stable (Blue)")
)
write.csv(stability_df, "Strain_Stability_Results_MLST.csv", row.names = FALSE)

# 5. Export Figure
pdf("Phylogenetic_Stability_Analysis_MLST.pdf", width = 12, height = 13)

# oma = c(bottom, left, top, right)
par(oma = c(12, 1, 6, 1), mar = c(2, 1, 2, 1)) 

obj <- cophylo(p1, p2, rotate = TRUE, iterations = 1000)

plot(obj, 
     link.type = "curved", 
     link.col = edge_cols, 
     link.lwd = 2, 
     fsize = 0.8, 
     ftype = "i", 
     pts = FALSE)

# --- VISIBLE LABELS ---
mtext("Phylogenetic Stability & Congruence Analysis", side = 3, line = 3, cex = 1.5, font = 2, outer = TRUE)
mtext("Reference Phylogeny", side = 3, line = 0, at = 0.2, cex = 1.2, font = 3, outer = TRUE)
mtext("MLST Phylogeny", side = 3, line = 0, at = 0.8, cex = 1.2, font = 3, outer = TRUE)

# Stats Line
stats_text <- paste0("Common Strains: ", length(common), 
                     "  |  Congruence Score: ", round(score, 4),
                     "  |  Unstable Strains (Red): ", length(unstable_strains))
mtext(stats_text, side = 1, line = 2, cex = 1.1, font = 1, outer = TRUE)

# --- THE LEGEND (FORCED VISIBILITY) ---
par(xpd = NA) 
legend(x = mean(par("usr")[1:2]), 
       y = par("usr")[3] - (diff(par("usr")[3:4]) * 0.1), # Lowered to 20% below plot
       legend = c("Stable Topology (Corr > 0.85)", "Unstable/Misplaced (Corr < 0.85)"), 
       lwd = 8, 
       col = c("#2E86C1", "#E74C3C"), 
       bty = "n", 
       horiz = TRUE, 
       xjust = 0.5, 
       cex = 1.2)

dev.off()

cat("✅ CSV Created: Strain_Stability_Results_MLST.csv\n")
cat("✅ PDF Created: Phylogenetic_Stability_Analysis_MLST.pdf\n")