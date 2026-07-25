This repository contains the dataset, custom scripts, and phylogenomic results supporting our manuscript: "A Highly Congruent Single-Locus Core Marker for Resolution of Aeromonas Species Dynamics"

.
├── code/
│   ├── 0-Genomics V1.ipynb
│   ├── 1-ext_gene_locus_tag.py
│   ├── 2-ext_genes_seq.sh
│   ├── 3-rename_fasta_header.py
│   ├── 4-align_MAFFT.sh
│   ├── 5-FastTree.sh
│   ├── 6-best_cgSNPs_fit.R
│   ├── 7-MLSTvsALT_v8.sh
│   ├── 8-Phylogenetic_Stability_Analysis.R
│   └── 9-extract_genes.py
├── data/
│   ├── Discovery_set_checkM_filtration.xlsx
│   ├── Discovery_set_rawdata_accessions.xlsx
│   ├── Validation_set_checkM_filtration.xlsx
│   ├── Validation_set_pairwise_distance_prc.xlsx
│   ├── Validation_set_rawdata_accessions.xlsx
│   └── Validation_set_skani_ANI_matrix.xlsx
└── results/
    ├── alignments/
    │   ├── Discovery_mlst_concatenated_aligned.fasta
    │   ├── Discovery_prc_aligned.fasta
    │   ├── Validation_prc_alignment.fasta
    │   ├── Validation_query_prc_seq.fasta
    │   └── MLST genes/
    │       ├── gltA_aligned.fasta
    │       ├── groL_aligned.fasta
    │       ├── gyrB_aligned.fasta
    │       ├── metG_aligned.fasta
    │       ├── ppsA_aligned.fasta
    │       └── recA_aligned.fasta
    ├── tang_MLST/
    │   ├── gltA_genes_renamed.nwk
    │   ├── groL_genes_renamed.nwk
    │   ├── gyrB_genes_renamed.nwk
    │   ├── metG_genes_renamed.nwk
    │   ├── ppsA_genes_renamed.nwk
    │   ├── prc_genes_renamed.nwk
    │   ├── recA_genes_renamed.nwk
    │   ├── Phylogenetic_Stability_Analysis_gltA.pdf
    │   ├── Phylogenetic_Stability_Analysis_groL.pdf
    │   ├── Phylogenetic_Stability_Analysis_gyrB.pdf
    │   ├── Phylogenetic_Stability_Analysis_metG.pdf
    │   ├── Phylogenetic_Stability_Analysis_ppsA.pdf
    │   ├── Phylogenetic_Stability_Analysis_recA.pdf
    │   ├── Strain_Stability_Results_gltA.csv
    │   ├── Strain_Stability_Results_groL.csv
    │   ├── Strain_Stability_Results_gyrB.csv
    │   ├── Strain_Stability_Results_metG.csv
    │   ├── Strain_Stability_Results_ppsA.csv
    │   └── Strain_Stability_Results_recA.csv
    └── trees/
        ├── MLST_tree.nwk
        ├── Phylogenetic_Stability_Analysis_MLST.pdf
        ├── Phylogenetic_Stability_Analysis_MLSTvsPRC.pdf
        ├── Phylogenetic_Stability_Analysis_prc.pdf
        ├── Reference_tree.nwk
        ├── Strain_Stability_Results_MLST.csv
        ├── Strain_Stability_Results_MLSTvsPRC.csv
        ├── Strain_Stability_Results_prc.csv
        └── prc_tree.nwk
