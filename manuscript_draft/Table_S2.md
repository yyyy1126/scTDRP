# Table S2 | GSE130116 B-ALL and Healthy Donor Sample Information

## B-ALL Diagnostic Samples

| GEO Accession | Sample Name | Genetic Subtype | Age (years) | Sex | Source Institute | Cell Input | QC-Pass Cells | B-Cell Selected |
|:---|:---|:---|:---:|:---:|:---|:---|:---:|:---:|
| GSM3732336 | ETV001_NYU_DIAGNOSIS | ETV-RUNX1 | 11 | Male | NYU Langone Health | 1:5 CD19+:CD19−CD45+ BM | 3,833 | 2,641 |
| GSM3732338 | ETV002_COG_DIAGNOSIS | ETV-RUNX1 | 2 | Male | Children's Oncology Group (COG) | 1:5 CD19+:CD19−CD45+ BM | 1,648 | 454 |
| GSM3732341 | ETV003_SJ_DIAGNOSIS | ETV-RUNX1 | 9 | Female | St. Jude Children's Research Hospital | 1:5 CD19+:CD19−CD45+ BM | 668 | 236 |
| GSM3732344 | ETV004_COG_DIAGNOSIS | ETV-RUNX1 | 2 | Male | Children's Oncology Group (COG) | 1:5 CD19+:CD19−CD45+ BM | 650 | 151 |
| GSM3732347 | ETV005_COG_DIAGNOSIS | ETV-RUNX1 | 2 | Male | Children's Oncology Group (COG) | 1:5 CD19+:CD19−CD45+ BM | 248 | 52 |
| GSM3732354 | PH001_COG_DIAGNOSIS | Ph+ (BCR-ABL1) | 12 | Male | Children's Oncology Group (COG) | 1:5 CD19+:CD19−CD45+ BM | 36 | 17 |
| GSM3732357 | PH002_SJ_DIAGNOSIS | Ph+ (BCR-ABL1) | 8 | Female | St. Jude Children's Research Hospital | 1:5 CD19+:CD19−CD45+ BM | 56 | 18 |
| **Total** | — | — | — | — | — | — | **7,139** | **3,569** |

## Healthy Donor Samples

| GEO Accession | Sample Name | Age (years) | Sex | Source | QC-Pass Cells |
|:---|:---|:---:|:---:|:---|:---:|
| GSM3732350 | HEALTHY001_LONZA_CD45 | 21 | Male | Lonza commercial BM | 22 |
| GSM3732351 | HEALTHY002_WC4_CD45 | 42 | Male | WC4 donor BM | 16 |
| GSM3732352 | HEALTHY003_SC_CD45 | 25 | Female | SC donor BM | 16 |
| GSM3732353 | HEALTHY004_SC_CD45 | 20 | Male | SC donor BM | 3 |
| **Total** | — | — | — | — | **57** |

## Notes

1. **Data source:** Witkowski et al., *Cancer Cell* 37(6):867–882.e12 (2020); GEO GSE130116.
2. **Platform:** 10x Genomics Chromium v2; shared barcode whitelist (737,280 barcodes).
3. **QC criteria:** n_counts ≥ 200, n_genes ≥ 100, pct_mt < 20%.
4. **B-cell selection:** B-cell score > median (score_genes with CD19, MS4A1, CD79A, CD79B, PAX5, VPREB1, RAG1, RAG2).
5. Healthy donor cells were not used as the normal reference due to low cell numbers after QC; the normal B-cell trajectory was constructed from 16,878 B cells curated from the blood map dataset (Pro-B VDJ, Large Pre-B, Immature B, Mature B).
