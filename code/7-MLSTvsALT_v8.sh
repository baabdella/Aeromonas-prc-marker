#!/bin/bash

# ==============================================================================
# 🧬 MASTER PIPELINE: MLST vs ALTERNATIVE MARKERS vs cgSNP (v7 - Robust)
# Handling Single/Multiple Genes & Missing Data
# ==============================================================================

# --- 1. Configuration ---
MLST_GENES="gyrB groL gltA metG recA ppsA" 
ALT_GENES="bamC" # Can be one or multiple    corC fliF hpt



ROARY_CSV="gene_presence_absence.csv"
SNP_TREE="cgSNP_reference.nwk"
BAKTA_DIR_NAME="1_bakta_result" 
THREADS=16

# Workspace Setup
MASTER_DIR="PanPhylogeny_$(date +%Y%m%d_%H%M)"
mkdir -p "$MASTER_DIR/MLST/1_raw" "$MASTER_DIR/MLST/2_aligned"
mkdir -p "$MASTER_DIR/ALT/1_raw" "$MASTER_DIR/ALT/2_aligned"
mkdir -p "$MASTER_DIR/Comparison"

# Path Detection
if [ -d "../$BAKTA_DIR_NAME" ]; then 
    BAKTA_PATH="$(dirname "$(pwd)")/$BAKTA_DIR_NAME"
else 
    BAKTA_PATH="$(pwd)/$BAKTA_DIR_NAME"
fi

# Check for cgSNP tree before starting
if [ ! -f "$SNP_TREE" ]; then
    echo "❌ Error: $SNP_TREE not found. Please provide the reference tree."
    exit 1
fi

# --- 2. System Check (Optional - handled gracefully) ---
FORTRAN_LIB=$(find /usr/lib/x86_64-linux-gnu/ -name "libgfortran.so.[45]" | head -n 1)
if [ ! -z "$FORTRAN_LIB" ] && [ ! -L /usr/lib/x86_64-linux-gnu/libgfortran.so ]; then 
    echo "Found gfortran lib at $FORTRAN_LIB. Creating symlink..."
    sudo ln -sf "$FORTRAN_LIB" /usr/lib/x86_64-linux-gnu/libgfortran.so 2>/dev/null || echo "Note: Could not create symlink, continuing anyway."
fi

# --- 3. Processing Function ---
run_pipeline() {
    local SET_NAME=$1
    local GENE_LIST=$2
    local TARGET_DIR="$MASTER_DIR/$SET_NAME"

    echo "🚀 Starting $SET_NAME extraction (Query: $GENE_LIST)..."
    python3 <<EOF
import pandas as pd
import os, pathlib, re
from Bio import SeqIO

def find_best_gene_match(target, available_genes):
    if target in available_genes: return target
    # Match exact or with Roary suffixes like _1, _2
    pattern = re.compile(rf"^{target}(_[0-9]+)?$", re.IGNORECASE)
    matches = [g for g in available_genes if pattern.match(str(g))]
    return matches[0] if matches else None

active_query = "$GENE_LIST".split()
df = pd.read_csv("$ROARY_CSV", low_memory=False)
all_roary_genes = df['Gene'].astype(str).tolist()
strains = [c for c in df.columns if c not in ['Gene','Non-unique Gene name','Annotation','No. isolates','No. sequences','Avg sequences per isolate','Genome Fragment','Order within Fragment','Accessory Fragment','Accessory Order with Fragment','QC','Min group size nuc','Max group size nuc','Avg group size nuc']]

tag_map = {}
found_any = False
for q in active_query:
    best_match = find_best_gene_match(q, all_roary_genes)
    if best_match:
        found_any = True
        subset = df[df['Gene'] == best_match][['Gene'] + strains].dropna(axis=1, how='all')
        for s in strains:
            if s in subset.columns:
                val = subset[s].values[0]
                if pd.notna(val):
                    # Handle multiple tags in one cell
                    for t in str(val).replace('\t', ' ').split():
                        tag_map[t.strip()] = (s, best_match)

if not found_any:
    print(f"⚠️ No genes from list {active_query} found in $ROARY_CSV")
    exit(0)

ffn_files = list(pathlib.Path("$BAKTA_PATH").rglob("*.ffn"))
for ffn in ffn_files:
    try:
        for rec in SeqIO.parse(str(ffn), "fasta"):
            if rec.id in tag_map:
                s_name, g_name = tag_map[rec.id]
                rec.id = s_name.replace(" ", "_").replace(".", "_").replace("-", "_")
                rec.description = ""
                with open(os.path.join("$TARGET_DIR/1_raw", f"{g_name}.fasta"), "a") as f:
                    SeqIO.write(rec, f, "fasta")
    except Exception as e: pass
EOF

    VALID_ALIGNMENTS=()
    for f in "$TARGET_DIR/1_raw"/*.fasta; do
        if [ -f "$f" ] && [ "$(grep -c ">" "$f")" -gt 2 ]; then
            gene_name=$(basename "$f" .fasta)
            echo "  - Aligning $gene_name..."
            mafft --auto --thread $THREADS "$f" > "$TARGET_DIR/2_aligned/${gene_name}_aligned.fasta" 2>/dev/null
            VALID_ALIGNMENTS+=("$TARGET_DIR/2_aligned/${gene_name}_aligned.fasta")
        fi
    done

    # --- Single vs Multiple Gene Logic ---
    if [ ${#VALID_ALIGNMENTS[@]} -eq 1 ]; then
        echo "✅ Single gene detected for $SET_NAME. Proceeding to tree..."
        cp "${VALID_ALIGNMENTS[0]}" "$TARGET_DIR/concatenated.fasta"
        fasttree -nt -gtr < "$TARGET_DIR/concatenated.fasta" > "$TARGET_DIR/${SET_NAME}_tree.nwk" 2>/dev/null
    elif [ ${#VALID_ALIGNMENTS[@]} -gt 1 ]; then
        echo "✅ ${#VALID_ALIGNMENTS[@]} genes detected for $SET_NAME. Concatenating..."
        seqkit concat "${VALID_ALIGNMENTS[@]}" -o "$TARGET_DIR/concatenated.fasta"
        fasttree -nt -gtr < "$TARGET_DIR/concatenated.fasta" > "$TARGET_DIR/${SET_NAME}_tree.nwk" 2>/dev/null
    else
        echo "❌ Error: No valid alignments for $SET_NAME. Tree will not be built."
    fi
}

run_pipeline "MLST" "$MLST_GENES"
run_pipeline "ALT" "$ALT_GENES"

# --- 4. R Analysis with Improved Checks ---
echo "📊 Running Triple Congruence Analysis..."

Rscript -e "
if (!require('phytools')) install.packages('phytools', repos='http://cran.us.project.org')
library(ape); library(phytools)

t_mlst <- if(file.exists('$MASTER_DIR/MLST/MLST_tree.nwk')) read.tree('$MASTER_DIR/MLST/MLST_tree.nwk') else NULL
t_alt  <- if(file.exists('$MASTER_DIR/ALT/ALT_tree.nwk')) read.tree('$MASTER_DIR/ALT/ALT_tree.nwk') else NULL
t_snp  <- if(file.exists('$SNP_TREE')) read.tree('$SNP_TREE') else NULL

analyze_pair <- function(tree1, tree2, n1, n2, genes1, genes2) {
    if(is.null(tree1) || is.null(tree2)) {
        message(paste('Skipping', n1, 'vs', n2, ': Tree not found'))
        return(NULL)
    }
    
    # Standardize labels
    tree1\$tip.label <- gsub('_', ' ', tree1\$tip.label)
    tree2\$tip.label <- gsub('_', ' ', tree2\$tip.label)
    common <- intersect(tree1\$tip.label, tree2\$tip.label)
    
    if(length(common) < 4) {
        message('Too few common taxa for comparison'); return(NULL)
    }

    p1 <- multi2di(unroot(keep.tip(tree1, common)))
    p2 <- multi2di(unroot(keep.tip(tree2, common)))
    
    d1 <- cophenetic(p1); d2 <- cophenetic(p2)[rownames(d1), colnames(d1)]
    stability <- sapply(rownames(d1), function(x) cor(d1[x,], d2[x,]))
    misplaced <- names(stability[stability < 0.90])
    
try({
        obj <- cophylo(p1, p2, rotate=TRUE)
        pdf(file.path('$MASTER_DIR/Comparison', paste0(n1, '_vs_', n2, '_final.pdf')), width=12, height=10)
        
        # --- ADD THIS LINE TO CREATE SPACE FOR TEXT ---
        # mar = c(bottom, left, top, right)
        par(mar = c(5, 1, 5, 1)) 
        
        cols <- rep(make.transparent('#5D6D7E', 0.2), length(common)); names(cols) <- common
        cols[misplaced] <- 'red'
        
        plot(obj, link.col=cols, link.lwd=2, fsize=0.8, ftype='i', pts=FALSE)
        
        # Header labels
        mtext(paste(n1, 'Tree'), side=3, line=2, at=0.1, cex=1.3, font=2, col='blue')
        mtext(paste('Genes:', genes1), side=3, line=1, at=0.1, cex=0.7, col='blue')
        mtext(paste(n2, 'Tree'), side=3, line=2, at=0.9, cex=1.3, font=2, col='darkgreen')
        
        # Footer Score
        rf_val <- as.numeric(dist.topo(p1, p2, method='PH85'))
        score <- 1 - (rf_val / (2 * (length(common) - 3)))
        mtext(paste0('Congruence (1-RF): ', round(score, 4)), side=1, line=2, font=3)
        
        dev.off()
        
        write.csv(data.frame(Strain=names(stability), Stability_Score=stability, Genes_Analyzed=genes1), 
                  file.path('$MASTER_DIR/Comparison', paste0(n1, '_vs_', n2, '_stability.csv')), row.names=F)
        return(data.frame(Pair=paste(n1, 'vs', n2), Congruence=score, Genes_Used=genes1))
    })
}

summary_stats <- rbind(
    analyze_pair(t_mlst, t_snp, 'MLST', 'cgSNP', '$MLST_GENES', 'Genome'),
    analyze_pair(t_alt, t_snp, 'ALT', 'cgSNP', '$ALT_GENES', 'Genome'),
    analyze_pair(t_alt, t_mlst, 'ALT', 'MLST', '$ALT_GENES', '$MLST_GENES')
)
if(!is.null(summary_stats)) write.csv(summary_stats, '$MASTER_DIR/Comparison/CONGRUENCE_SUMMARY.csv', row.names=F)
"

# --- 5. Final Report with File Handling ---
python3 <<EOF
import pandas as pd
import os

comp_dir = "$MASTER_DIR/Comparison"
mlst_file = f"{comp_dir}/MLST_vs_cgSNP_stability.csv"
alt_file = f"{comp_dir}/ALT_vs_cgSNP_stability.csv"

if os.path.exists(mlst_file) and os.path.exists(alt_file):
    mlst_snp = pd.read_csv(mlst_file)
    alt_snp  = pd.read_csv(alt_file)

    df = mlst_snp.rename(columns={'Stability_Score': 'MLST_vs_SNP', 'Genes_Analyzed': 'MLST_Genes'})
    df = df.merge(alt_snp.rename(columns={'Stability_Score': 'ALT_vs_SNP', 'Genes_Analyzed': 'ALT_Genes'}), on='Strain')

    def judge(row):
        m, a = row['MLST_vs_SNP'], row['ALT_vs_SNP']
        status = "ALT better" if a > m + 0.05 else ("MLST better" if m > a + 0.05 else "Comparable")
        if a < 0.8 and m < 0.8: status += " (Poor)"
        return status

    df['Judgment'] = df.apply(judge, axis=1)
    df.to_csv(f"{comp_dir}/FINAL_STRAIN_JUDGMENT.csv", index=False)
    print("✅ Final strain judgment created.")
else:
    print("⚠️ Skipping final judgment: One or more comparison files are missing (likely because a tree failed).")
EOF

echo "🏁 Process Complete. Results: $MASTER_DIR/Comparison"
