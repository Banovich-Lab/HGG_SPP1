This repository contains the custom code to reproduce the results presented in:

## Overcoming myeloid-driven resistance to CAR T therapy by targeting SPP1 

Sharareh Gholamin<sup>1,2,3</sup>*<sup>†</sup>, Heini M Natri<sup>4</sup>*, Yuqi Zhao<sup>5</sup>, Shengchao Xu<sup>1</sup>, Maryam Aftabizadeh<sup>1</sup></sup>, Begonya Comin-Anduix<sup>6</sup>, Supraja Saravanakumar<sup>1</sup>, Christian Masia<sup>1</sup>, Robyn A. Wong<sup>1</sup>, Lance Peter<sup>4</sup>, Mei-i Chung<sup>4</sup>, Evan D Mee<sup>4</sup>, Brenda Aguilar<sup>1</sup>, Renate Starr<sup>1</sup>, Davis Y. Torrejon<sup>7</sup>, Darya Alizadeh<sup>1</sup>, Xiwei Wu<sup>5</sup>, Anusha Kalbasi<sup>8</sup>, Antoni Ribas<sup>7</sup>, Stephen Forman<sup>1</sup>, Behnam Badie<sup>9</sup>, Nicholas E Banovich<sup>4</sup>**, Christine E Brown<sup>1</sup>**<sup>†</sup>

1. City of Hope Beckman Research Institute and Medical Center, Duarte, CA 91010, USA 
2. California Institute of Technology, Division of Biology and Biological Engineering, Pasadena, CA 91125, USA 
3. Department of Radiation Oncology, City of Hope, Duarte, CA 91010, USA 
4. Division of Bioinnovation and Genome Sciences, The Translational Genomics Research Institute (TGen), Phoenix, AZ 85004, USA 
5. Department of Computational and Quantitative Medicine within Beckman Research Institute, City of Hope, Duarte, CA 91010, USA 
6. Division of Surgical Oncology, Department of Surgery, UCLA, Los Angeles, CA 90095, USA. 
7. Division of Hematology-Oncology, Department of Medicine, UCLA, Los Angeles, CA 90095, USA 
8. Division of Radiation Oncology, Department of Medicine, Stanford, Palo Alto, CA 94305, USA 
9. Division of Neurosurgery, City of Hope, Duarte, CA 91010, USA 

 
&ast; indicates shared first authorship. 
&ast;&ast; indicates shared senior authorship. 
<sup>†</sup> indicates shared correspondence. 

Manuscript correspondence to: cbrown@coh.org, sgholamin@coh.org
Code/Github repository contact: hnatri@tgen.org

Preprint URL pending.

## Single cell RNA-sequencing data processing and analysis

Raw and processed GBM/HGG scRNA-seq data are deposited at the Gene Expression Omnibus (GEO) under accession number GSE290291 and mouse data under accession GSE289327. Code for the processing and initial analysis of the mouse data can be found in a separate repository at [https://github.com/zhaoyuqi616/SPP1-Targeted_CAR_T_Therapy](https://github.com/zhaoyuqi616/SPP1-Targeted_CAR_T_Therapy)).

### Pre-processing and integration

1. Processing batches 1-19 and samples run without multiplexing: 13384_first_batches_preprocessing.R
2. Processing batches 42-1, 42-2, 43-1, 43-2, 44-2, 45, 46, and 47: 13384_tumor_preprocessing.R
3. Preparing all sample data for integration: 13384_tumor_prepare_integration.R
4. Detecting ambient RNA using SoupX and integrating: 13384_tumor_soupX_integrate.R
5. Doublet detection using DoubletFinder: 13384_tumor_DoubletFinder.R
6. Filtering doublets: 13384_filter_doublets.R

### Analysis and figures

* Figure 1 and Extended Figures 1-2
⋅⋅* Fig1_13384_GBM_scRNAseq.R
⋅⋅* Fig1_GBMHGG_CellChat.R
* Figure 2
⋅⋅* Fig2_GBMHGG_SPP1_survival_wholetumor.R
⋅⋅* Fig2_GBMHGG_SPP1_HR.R
⋅⋅* Fig2_GBMHGG_SPP1_survival.R
⋅⋅* Fig2_GBMHGG_SPP1_survival_plot.R
⋅⋅* Fig2_GBMHGG_GSEA.R
⋅⋅* Fig2_GBMHGG_DEGs.R
* Figure 3
⋅⋅* Fig3_GBMHGG_JAK1KO_integrated.R
⋅⋅* Fig3_JAK1KO_WT_CellChat.R
* Figure 4
⋅⋅* Fig4_13384_GBM_JAK1KO_WT_compare.R
* Figure 5
⋅⋅* Fig5_bulkRNAseq.R
⋅⋅* Fig5_bulkRNAseq_ctdeconv.R
* Extended Figure 1, panel c
⋅⋅* ExtendedFig_GBM_inferCNV.R
⋅⋅* ExtendedFig_GBM_inferCNV_result.R
* Quality metrics
⋅⋅* ExtendedFig_GBMHGG_QC.R

