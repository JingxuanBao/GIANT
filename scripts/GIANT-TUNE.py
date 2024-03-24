import argparse

## Take commandline input
parser = argparse.ArgumentParser(description='GIANT-tune: Genetically Informed brAiN aTlas (Version 2)\n Please refer to the GitHub for more information: https://github.com/JingxuanBao/GIANT\n Contact:jingxuan.bao@pennmedicine.upenn.edu.',
    usage='use "python %(prog)s --help" for more information',
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--heritability', type=str, required=True, help='Path to heritability brainmap NIFTI file.')
parser.add_argument('--ref', type=str, required=True, help='Path to reference atlas NIFTI file.')
parser.add_argument('--c', type=float, default=1, help='Biological validaty hyperparameter.')
parser.add_argument('--atlas', type=str, required=True, nargs='+', help='Path to GIANT atlas NIFTI file.')
parser.add_argument('--mask', type=str, required=True, help='Path to mask NIFTI file.')
parser.add_argument('--steps', type=float, default=6, help='The radius (level=sphere) or number of steps (level=step).')
parser.add_argument('--neighbor', type=str, choices=['sphere','step'], default='sphere', help='How to define the neighbor? sphere or step.')
parser.add_argument('--out_prefix', type=str, required=True, help='Prefix path to your output (parcellations, MCMC chain, log-likelihood chain).')
args = parser.parse_args()

import skimage
import pandas as pd
import numpy as np
import nibabel as nib
import nilearn
from nilearn import datasets, plotting
from datetime import datetime
import scipy.io
import sys
import matplotlib.pyplot as plt
import os
from itertools import compress
from sklearn.metrics import calinski_harabasz_score,davies_bouldin_score
import math
import networkx as nx
import sys
import seaborn as sns
import matplotlib.pyplot as plt
 
def Generate_Neighbor_IdxMat(mask,level=6,neighbor="sphere"):

    # This function returns a neighbor matrix for each non-background voxel
    # The matrix has number of rows equals to number of non-background voxels
    # The matrix has number of columns equals to number of neighborhood voxels for given voxel
    # All negative values in the matrix means the corresponding neighborhood voxel is masked out
    # The order of voxel (rows of matrix) is equal to img.ravel()[mask.ravel()]

    ## neighbor defines how the neighborhood voxels are obtained
    ### "sphere" means all the voxels within the "level"(default is 3) Euclidean distance of center voxel
    ### "step" means all the voxels that can be reached by walk "level" number of steps from the center voxel

    within_subset_idx_mat = -np.ones(mask.shape)
    within_subset_idx_mat[mask] = np.arange(np.sum(mask))

    within_linear_idx_mat = np.arange(np.prod(mask.shape)).reshape(mask.shape)
    within_linear_idx_mat[mask==False] = -1

    R,C,P = np.meshgrid(np.arange(2*level+1)-level, np.arange(2*level+1)-level,np.arange(2*level+1)-level,indexing='ij')
    if neighbor=="step":
        connectivity = np.add.reduce(np.fabs(np.indices([level*2+1,level*2+1,level*2+1])-level), 0)
        connectivity = (connectivity<=level)
    elif neighbor=="sphere":
        connectivity = (R**2+C**2+P**2<=level**2)
    else:
        raise SystemExit("neighbor must be step or sphere!")
    #connectivity=connectivity+connectivity-1
    connectivity[(level,level,level)]=False
    ncol = np.sum(connectivity)
    nrow = np.sum(mask)
    neighbor_mat = np.zeros((nrow,ncol))
    mask_flat = mask.ravel()
    y_idx = np.arange(np.prod(mask.shape))[mask_flat]
    for curr_idx in range(y_idx.shape[0]):
        if curr_idx==0 or curr_idx%10000==0:
            print(curr_idx,end=" ")
            sys.stdout.flush()
        r,c,p = np.unravel_index(y_idx[curr_idx], mask.shape)
        r = (r+R)[connectivity].ravel()
        c = (c+C)[connectivity].ravel()
        p = (p+P)[connectivity].ravel()
        idx_greater0 = (r>0) & (c>0) & (p>0)
        if np.sum(idx_greater0)!=r.shape[0]:
            # If some voxels have negative coordinate
            neighbor_before_limitcheck = -np.ones(r.shape[0])
            r = r[idx_greater0]
            c = c[idx_greater0]
            p = p[idx_greater0]
            neighbor_before_limitcheck[idx_greater0] = np.array([within_subset_idx_mat[(r[i],c[i],p[i])] for i in range(r.shape[0])])
            neighbor_mat[curr_idx,:] = neighbor_before_limitcheck
        else:
            # If no negative coordinate
            neighbor_mat[curr_idx,:] = np.array([within_subset_idx_mat[(r[i],c[i],p[i])] for i in range(r.shape[0])])
    return (within_subset_idx_mat.ravel()[mask_flat]).astype(int),neighbor_mat.astype(int),(within_linear_idx_mat.ravel()[mask_flat]).astype(int)

def Neighbor_NumClust(subset_idx_connectivity,atlas,mask):
    atlas_inear = atlas[mask]
    num_neighbor_clust = np.zeros(subset_idx_connectivity.shape[0])
    for idx in range(subset_idx_connectivity.shape[0]):
        good_point_idx = subset_idx_connectivity[idx,:]!=-1
        num_neighbor_clust[idx] = np.unique(atlas_inear[subset_idx_connectivity[idx,good_point_idx]]).shape[0]
    return np.mean(num_neighbor_clust),np.std(num_neighbor_clust)

def Create_Network_from_Connectivity(connectivity):
    G = nx.Graph()
    G.add_nodes_from(list(range(connectivity.shape[0])))
    for node in range(connectivity.shape[0]):
        if node==0 or node%10000==0:
            print(node,end=" ")
            sys.stdout.flush()
        G.add_edges_from([(node,x.astype(float)) for x in connectivity[node,:][connectivity[node,:]>=0]])
    return G

def Check_Connected(network,atlas):
    within_subset_idx_mat = -np.ones(atlas.shape)
    within_subset_idx_mat[atlas>0] = np.arange(np.sum(atlas>0))
    rois = np.unique(atlas[atlas>0])
    count = np.zeros(rois.shape[0])
    for roi in range(rois.shape[0]):
        #print(roi,end=" ")
        sys.stdout.flush()
        sub_idx = within_subset_idx_mat[atlas==rois[roi]]
        sub_network = network.subgraph(sub_idx)
        count[roi] = len(list(nx.connected_components(sub_network)))
    print(" ")
    sys.stdout.flush()
    return rois,count

img = nib.load(args.mask) 
mask = img.get_fdata()
mask = mask ==1

_,subset_idx_connectivity,_ = Generate_Neighbor_IdxMat(mask,level=args.steps,neighbor=args.neighbor)

her_3d = nib.load(args.heritability)
her_3d = her_3d.get_fdata()

network_masked = Create_Network_from_Connectivity(subset_idx_connectivity)

c=args.c

old_stdout = sys.stdout
log_file = open(args.out_prefix+".gtune","w")
sys.stdout = log_file

### reference atlas
img = nib.load(args.ref)
data = img.get_fdata()

CHscore = calinski_harabasz_score(her_3d[mask].reshape((-1,1)),data[mask])
DBscore = davies_bouldin_score(her_3d[mask].reshape((-1,1)),data[mask])
aveNumClust_bl,sdNumClust_bl = Neighbor_NumClust(subset_idx_connectivity,data,mask)
rois,count = Check_Connected(network_masked,data)

print("All regions with same index are connected (ref)? "+str(all(count==1)))
if (all(count==1)==False):
    print("Not connected regions (ref): "+str(rois[count>1]))
print("Number of unique neighborhood regions per voxel (ref): "+str(aveNumClust_bl)+"+-"+str(sdNumClust_bl))
print("CH score (ref): "+str(CHscore))
print("DB score (ref): "+str(DBscore))
print(" ")

### target atlases
CH_best = 0.0
DB_best = math.inf
for element in args.atlas:
    print(element) 
    img_target = nib.load(element)
    atlas = img_target.get_fdata()

    CHscore = calinski_harabasz_score(her_3d[mask].reshape((-1,1)),atlas[mask])
    DBscore = davies_bouldin_score(her_3d[mask].reshape((-1,1)),atlas[mask])
    aveNumClust,sdNumClust = Neighbor_NumClust(subset_idx_connectivity,atlas,mask)
    rois,count = Check_Connected(network_masked,atlas)

    print("All regions with same index are connected (target)?"+str(all(count==1)))
    if (all(count==1)==False):
        print("Not connected regions (target): "+str(rois[count>1]))
    print("Number of unique neighborhood regions per voxel (target): "+str(aveNumClust)+"+-"+str(sdNumClust))   
    print("CH score (target): "+str(CHscore))
    print("DB score (target): "+str(DBscore))
    
    if aveNumClust>(aveNumClust_bl-c*sdNumClust_bl) and aveNumClust<(aveNumClust_bl+c*sdNumClust_bl):
        if CHscore>CH_best:
            CH_best = CHscore.copy()
            CH_atlas = element
        if DBscore<DB_best:
            DB_best = DBscore.copy()
            DB_atlas = element
    else:
        print("This atlas didn't pass the cluster number check!")
        
    print(" ")

print("The best tunning atlas:")
print("CH score (best): "+str(CH_best))
print("CH atlas (best): "+CH_atlas)
print("DB score (best): "+str(DB_best))
print("DB atlas (best): "+DB_atlas)

sys.stdout = old_stdout
log_file.close()

# python ./GIANT-TUNE.py --heritability /gpfs/fs001/cbica/home/baoji/230626-GIANT-Atlas/8-out/gm_in.nii.gz --ref /gpfs/fs001/cbica/home/baoji/230626-GIANT-Atlas/Data/Atlas/BLSA_SPGR+MPRAGE_averagetemplate_muse.nii.gz --atlas /gpfs/fs001/cbica/home/baoji/230626-GIANT-Atlas/9-out/GIANT_gm_in_0.5-atlas-24344.0.nii.gz /gpfs/fs001/cbica/home/baoji/230626-GIANT-Atlas/9-out/GIANT_gm_in_1-atlas-20412.0.nii.gz /gpfs/fs001/cbica/home/baoji/230626-GIANT-Atlas/9-out/GIANT_gm_in_20-atlas-53.0.nii.gz /gpfs/fs001/cbica/home/baoji/230626-GIANT-Atlas/9-out/GIANT_gm_in_30-atlas-53.0.nii.gz --mask /gpfs/fs001/cbica/home/baoji/230626-GIANT-Atlas/9-out/gm_mask.nii.gz --out_prefix /gpfs/fs001/cbica/home/baoji/230626-GIANT-Atlas/11-out/test