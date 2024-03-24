## Usage: python BrainRegionNaming.py --target_atlas "D:/Box/02-Research/220311-BayesSpace-BrainAtlas/221019-Updated/3-3-out/CH_atlas.nii.gz" --ref_atlas "D:/Box/02-Research/220311-BayesSpace-BrainAtlas/1-Data/MUSE_atlas/BLSA_SPGR+MPRAGE_averagetemplate_muse_222.nii.gz" --ref_atlas_specification "D:/Box/02-Research/220311-BayesSpace-BrainAtlas/Manuscript/BrainRegionNaming/MUSE_ROI_specification.csv" --out "D:/Box/02-Research/220311-BayesSpace-BrainAtlas/Manuscript/BrainRegionNaming/CH_HABP_Atlas.csv" --threshold 0.05
import skimage
import pandas as pd
import numpy as np
import nibabel as nib
import nilearn
from datetime import datetime
import scipy.io
import sys
import os
from itertools import compress

import math
import sys
import re

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--target_atlas', type=str, required=True,help='Path to atlas file you want to annotate.')
parser.add_argument('--ref_atlas', type=str, required=True,help='Path to the reference atlas file.')
parser.add_argument('--ref_atlas_specification', type=str, required=True,help='Path to the specification of the reference file.')
parser.add_argument('--out', type=str, required=True,help='Path to the specification of the reference file.')
parser.add_argument('--threshold', type=float, required=True,help='Ignore ROI with percentage lower than the threshold.')
args = parser.parse_args()

target_ori = nib.load(args.target_atlas)
target_atlas = target_ori.get_fdata()
ref_ori = nib.load(args.ref_atlas)
ref_atlas = ref_ori.get_fdata()
specification = pd.read_csv(args.ref_atlas_specification).values

target_ROI = np.unique(target_atlas)[1:] # exclude the background

threshold = args.threshold

rois = []
annotation_rois = []
annotation_high_rois = []
annotation_high_rois_one = []

xsorted=np.argsort(specification[:,0])
for roi in target_ROI:
    print(roi)
    rois.append(roi.astype(int))
    mask = target_atlas==roi
    ref_roi = ref_atlas[mask]
    label,count = np.unique(ref_roi,return_counts=True)
    percentage = count/np.sum(mask)
    label = label[percentage>=threshold]
    percentage = percentage[percentage>=threshold]
    percentage = (percentage.round(2)*100).astype(int)
    
    ypos = np.searchsorted(specification[xsorted,0], label)
    index = xsorted[ypos]
    
    roi_specification = specification[index,:]
    annotation1 = [str(percentage[i])+" "+str(roi_specification[i,1]) for i in range(percentage.shape[0])]
    annotation1 = ' + '.join(annotation1)
    annotation_rois.append(annotation1)
    # print(percentage)
    # print(roi_specification)
    
    label, idx, _ = np.unique(roi_specification[:,2], return_counts=True, return_inverse=True)
    percentage = np.bincount(idx, percentage).astype(int) 
    annotation2 = [str(percentage[i])+" "+str(label[i]) for i in range(percentage.shape[0])]
    annotation2 = ' + '.join(annotation2)
    annotation_high_rois.append(annotation2)
    # print(annotation2)

    annotation_high_rois_one.append(label[np.argmax(percentage)])
    
df = pd.DataFrame(list(zip(rois, annotation_rois,annotation_high_rois,annotation_high_rois_one)),columns =['INDEX', 'NAME','GROUP','GROUP_UNIQ'])
df.to_csv(args.out,index=False)
