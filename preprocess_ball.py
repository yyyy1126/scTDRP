#!/usr/bin/env python3
"""
GSE130116 (Witkowski 2019 B-ALL) 预处理脚本

数据特点:
- 24个样本共享同一个 10x Genomics v2 barcode whitelist (737,280 barcodes)
- 每个样本的 matrix.mtx.gz 都是完整 33538 genes x 737280 barcodes 的raw unfiltered矩阵
- 需要根据non-zero columns识别每个样本的真实细胞

输出:
- GSE130116_combined_preprocessed.h5ad: 原始合并数据 (all assigned barcodes)
- GSE130116_qc_filtered_loose.h5ad: QC过滤后的数据
"""

import os
import scipy.io
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy.sparse import vstack

DATA_DIR = '../../1.data/raw/scRNAseq/Witkowski2019_BALL/'

# ========================== 读取注释文件 ==========================
print('Reading features...')
features = pd.read_csv(os.path.join(DATA_DIR, 'GSE130116_features.tsv.gz'),
                       sep='\t', header=None, compression='gzip')
features.columns = ['ensembl_id', 'gene_symbol', 'feature_type']
print(f'Features: {features.shape}')
print(f'Duplicate gene symbols: {features["gene_symbol"].duplicated().sum()}')

# 去重基因symbol
features['gene_symbol_unique'] = features['gene_symbol'].astype(str)
duplicates = features['gene_symbol'].duplicated(keep=False)
if duplicates.any():
    dup_counts = {}
    new_symbols = []
    for sym in features['gene_symbol']:
        if sym not in dup_counts:
            dup_counts[sym] = 0
            new_symbols.append(sym)
        else:
            dup_counts[sym] += 1
            new_symbols.append(f"{sym}.{dup_counts[sym]}")
    features['gene_symbol_unique'] = new_symbols
    print(f'Made {duplicates.sum()} duplicate symbols unique')

print('Reading mapping...')
mapping = pd.read_csv(os.path.join(DATA_DIR, 'cell_to_sample_mapping.csv'))
mapping = mapping[mapping['sample_id'] != 'unassigned'].copy()
mapping = mapping.reset_index(drop=True)
print(f'Mapping assigned: {len(mapping)}')

# ========================== 合并样本矩阵 ==========================
diagnosis_samples = sorted([s for s in mapping['sample_id'].unique() if 'DIAGNOSIS' in s])
healthy_samples = sorted([s for s in mapping['sample_id'].unique() if 'HEALTHY' in s])
print(f'Diagnosis samples: {len(diagnosis_samples)}')
print(f'Healthy samples: {len(healthy_samples)}')

all_barcodes = []
all_sample_ids = []
all_conditions = []
matrices = []

for sample in diagnosis_samples + healthy_samples:
    print(f'Processing {sample}...')
    mtx_path = os.path.join(DATA_DIR, f'{sample}.matrix.mtx.gz')
    with gzip.open(mtx_path, 'rt') as f:
        mat = scipy.io.mmread(f).tocsr()

    sample_indices = mapping[mapping['sample_id'] == sample].index.values
    print(f'  Sample cells: {len(sample_indices)}')

    sub_mat = mat[:, sample_indices].tocsr().T  # cells x genes
    sub_barcodes = mapping.iloc[sample_indices]['barcode'].values
    print(f'  Extracted shape: {sub_mat.shape}')

    matrices.append(sub_mat)
    all_barcodes.extend([f'{sample}_{b}' for b in sub_barcodes])
    all_sample_ids.extend([sample] * len(sub_barcodes))
    all_conditions.extend(['B-ALL_Diagnosis' if 'DIAGNOSIS' in sample else 'Healthy'] * len(sub_barcodes))
    del mat

print('\nConcatenating matrices...')
X_all = vstack(matrices, format='csr')
print(f'Combined matrix: {X_all.shape}')

obs = pd.DataFrame({
    'sample_id': all_sample_ids,
    'condition': all_conditions,
}, index=all_barcodes)

var = pd.DataFrame({
    'ensembl_id': features['ensembl_id'].values,
}, index=features['gene_symbol_unique'].values)

adata = sc.AnnData(X=X_all, obs=obs, var=var)
print(f'AnnData: {adata.shape}')
print(adata.obs['condition'].value_counts().to_string())
print('\nBy sample:')
print(adata.obs['sample_id'].value_counts().to_string())

output_path = os.path.join(DATA_DIR, 'GSE130116_combined_preprocessed.h5ad')
adata.write_h5ad(output_path)
print(f'Saved: {output_path}')

# ========================== QC过滤 ==========================
print('\n' + '='*60)
print('Quality Control Filtering')
print('='*60)

ncounts = np.array(adata.X.sum(axis=1)).flatten()
ngenes = np.array((adata.X > 0).sum(axis=1)).flatten()
mt_genes = adata.var_names.str.startswith('MT-')
pct_mt = np.nan_to_num(np.array(adata[:, mt_genes].X.sum(axis=1)).flatten() / ncounts * 100)

adata.obs['n_counts'] = ncounts
adata.obs['n_genes'] = ngenes
adata.obs['pct_mt'] = pct_mt

# Loose QC: n_counts >= 200, n_genes >= 100, pct_mt < 20%
mask = (ncounts >= 200) & (ngenes >= 100) & (pct_mt < 20)
print(f'QC pass (n_counts>=200, n_genes>=100, pct_mt<20): {mask.sum()} / {len(mask)}')

adata_qc = adata[mask].copy()
print(f'QC adata: {adata_qc.shape}')
print(adata_qc.obs['condition'].value_counts().to_string())
print('\nBy sample:')
print(adata_qc.obs['sample_id'].value_counts().to_string())

qc_path = os.path.join(DATA_DIR, 'GSE130116_qc_filtered_loose.h5ad')
adata_qc.write_h5ad(qc_path)
print(f'Saved: {qc_path}')
