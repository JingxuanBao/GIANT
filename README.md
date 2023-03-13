# Genetically informed brain atlas (GIANT) for enhancing neuroimaging genetics studies

This repository includes the **Genetically Informed brAiN Ttlas (GIANT)** with its implementation, and all imaging-genetics subsequent analysis summary statistics for the GIANT. 

## GIANT - A genetically informed brain atlas

We introduce GIANT, a genetically informed brain atlas. Our 3D clustering algorithm was applied separately to gray matter and white matter, and the best-tuned brain parcellations were integrated to define GIANT.  Our annotation results divided GIANT into 7 anatomical sub-structures: cerebellum (**A**), deep gray matter and white matter structure (**B**), frontal structure (**C** and **D**), occipital structure (**F**), parietal structure (**G**), temporal structure (**H**), and others (**E**).  

The NIFTI for GIANT and its specifications: https://github.com/JingxuanBao/GIANT/tree/main/GIANT

<p>
	<img src="figures/Annotation_1_8.png" width="864">
</p>

## 3D heritability-aware brain parcellation model

We propose a 3D heritability-aware brain parcellation model that integrates voxel-wise heritability and spatial proximity to cluster brain voxels into genetically informed regions.

You may find the implementation and example: https://github.com/JingxuanBao/GIANT/blob/main/scripts/Her_Atlas.py

## Voxel-level heritability estimates

The voxel-level heritability estiamtes for UK Biobank (UKBB): https://github.com/JingxuanBao/GIANT/tree/main/example_data/voxel_heritability_UKBB

The voxel-level heritability estimates for Alzheimer's disease neuroimaging initiative (ADNI): https://github.com/JingxuanBao/GIANT/tree/main/example_data/voxel_heritability_ADNI

## Regional-level heritability estimates for GIANT

The regional-level heritability estimates for UKBB and ADNI using GIANT brian parcellation: https://github.com/JingxuanBao/GIANT/tree/main/GIANT/GIANT_regional_heritability

## Genome-wide association study summary statistics

Genome builder: GRCh37



### Polygenetic risk scores

box folder link

### Regional-level genetic correlations

box folder link
