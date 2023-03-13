## Usage: python Her_Atlas.py gm 3 10 150 100 python Her_Atlas.py gm 3 10 150 100 ../example_data/conventional_atlas/convention_2x2x2mm3_atlas.nii.gz

######################
## Script for GIANT ##
######################
import numpy as np
#import numpy_indexed as npi
import nibabel as nib
import nilearn
from nilearn import plotting
import scipy
import skimage
import pandas as pd
import scipy.io
import os
from datetime import datetime
from sklearn.cluster import KMeans
import networkx as nx
import sys

def nums(first_number, last_number, step=1):
    return range(first_number, last_number+1, step)

def Create_Network_from_Connectivity(connectivity):
    G = nx.Graph()
    G.add_nodes_from(list(range(connectivity.shape[0])))
    for node in range(connectivity.shape[0]):
        if node==0 or node%10000==0:
            print(node,end=" ")
            sys.stdout.flush()
        G.add_edges_from([(node,x.astype(float)) for x in connectivity[node,:][connectivity[node,:]>=0]])
    return G

def ReAssignRegion(network,atlas):
    within_subset_idx_mat = -np.ones(atlas.shape)
    within_subset_idx_mat[atlas>0] = np.arange(np.sum(atlas>0))
    rois = np.unique(atlas[atlas>0])
    count = np.zeros(rois.shape[0])
    result_atlas = np.zeros(atlas.shape)
    for roi in range(rois.shape[0]):
        print(roi,end=" ")
        sys.stdout.flush()
        sub_idx = within_subset_idx_mat[atlas==rois[roi]]
        sub_network = network.subgraph(sub_idx)
        NumOfIsolated = len(list(nx.connected_components(sub_network)))
        cons_modification = 1.0/NumOfIsolated
        for sub_roi in range(NumOfIsolated):
            result_atlas[np.isin(within_subset_idx_mat,np.array(list(list(nx.connected_components(sub_network))[sub_roi])))] = rois[roi]+sub_roi*cons_modification
    result_atlas[result_atlas>0] = scipy.stats.rankdata(result_atlas[result_atlas>0],"dense")
    print(" ")
    sys.stdout.flush()
    return result_atlas

def Generate_Neighbor_IdxMat(mask,level=3,neighbor="sphere"):

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

def Her_Atlas(img,mask,affine,save_nifti_path_prefix,init_path,
              level=3,neighbor="sphere",
              chain_length=500,burn_in=200,
              alpha=1,beta=0.01,gamma=None,
              save_chain=None):

    # Define neighborhood voxels
    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Start to generate connectivity matrix ...")
    sys.stdout.flush()

    # Define neighborhood voxels
    subset_idx,subset_idx_connectivity,subsetL_idx = Generate_Neighbor_IdxMat(mask,level=level,neighbor=neighbor)
    voxel_with_neighbor = np.sum(subset_idx_connectivity>=0,axis=1)>0
    subset_idx_connectivity = subset_idx_connectivity[voxel_with_neighbor,:]
    neighbor_voxel_index = []
    for neigh_idx in range(subset_idx_connectivity.shape[0]):
        neighbor_idx_nonnegative = subset_idx_connectivity[neigh_idx,:]>=0
        neighbor_voxel_index.append(subset_idx_connectivity[neigh_idx,neighbor_idx_nonnegative])
    neighbor_num = np.array([x.shape[0] for x in neighbor_voxel_index])
    # np.sum(img[mask][subset_idx]==np.take(img,subsetL_idx)) # Uncomment to check if indeces are matching
    # Check if any voxel does not have neighbor
    voxel_with_neighbor = np.sum(subset_idx_connectivity>=0,axis=1)>0
    network = Create_Network_from_Connectivity(subset_idx_connectivity)

    print("")
    sys.stdout.flush()
    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Done!")
    sys.stdout.flush()
    print("There are "+str(subset_idx_connectivity.shape[0])+" voxels where "+str(np.sum(voxel_with_neighbor==False))+" of them are isolated.")
    sys.stdout.flush()


    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Read and preprocess data ...")
    sys.stdout.flush()
    y = img[subset_idx].reshape(-1)

    nrep = chain_length
    n=y.shape[0] # number of voxels
    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Done!")
    sys.stdout.flush()

    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Start initialization ...")
    sys.stdout.flush()
    #gamma = subset_idx_connectivity.shape[1]/4
    if gamma == None:
        gamma = np.sqrt(subset_idx_connectivity.shape[1])

    init_z = nib.load(init_path)
    init_z = init_z.get_fdata()
    q = np.unique(init_z[mask]).shape[0] # number of clusters

    df_sim_z = np.zeros((n,nrep))
    df_simu_mu = np.zeros((q,nrep))
    df_simu_sigma2 = np.zeros(nrep)

    ## Initialize parameters
    #kmeans = KMeans(n_clusters=q, random_state=23).fit(y.reshape(-1, 1))
    #df_sim_z[:,0] = kmeans.labels_+1
    ##df_sim_z[:,0] = init_z[mask]
    #_,_,count = np.unique(df_sim_z[:,0],return_inverse=True,return_counts=True)
    #df_simu_mu[:,0] = np.bincount(df_sim_z[:,0].astype(int), weights=y,minlength=q+1)[1:]/count
    df_sim_z[:,0] = scipy.stats.rankdata(init_z[mask],"dense")
    _,_,count = np.unique(df_sim_z[:,0],return_inverse=True,return_counts=True)
    df_simu_mu[:,0] = np.bincount(df_sim_z[:,0].astype(int), weights=y,minlength=q+1)[1:]/count

    mu0 = np.mean(y)
    #lambda0 = 0.01
    lambda0 = np.std(df_simu_mu[:,0])

    #_, sigma2 = npi.group_by(df_sim_z[:,0].astype(int)).std(y)
    df_simu_sigma2[0] = np.nanmean(np.sqrt(np.bincount(df_sim_z[:,0].astype(int), weights=(y-df_simu_mu[df_sim_z[:,0].astype(int)-1,0])**2,minlength=q+1)[1:]/count))
    #df_simu_sigma2[0] = 0.1
    alpha_iv_gamma=alpha+np.sum(voxel_with_neighbor)/2

    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Done!")
    sys.stdout.flush()

    ## MCMC
    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Start MCMC process ...")
    sys.stdout.flush()
    iterN = 0
    while iterN <(nrep-1):

        iterN=iterN+1
        print(iterN, end =" ")
        sys.stdout.flush()

        # simulate muk
        groupSumY = np.bincount(df_sim_z[:,iterN-1].astype(int), weights=y,minlength=q+1)[1:]
        unique_cluster,_,groupCount = np.unique(df_sim_z[:,iterN-1],return_inverse=True,return_counts=True)
        var_mu_k = np.zeros(q)
        var_mu_k[(unique_cluster-1).astype(int)] = 1/(groupCount/df_simu_sigma2[iterN-1]+1/lambda0)
        mean_mu_k = (groupSumY/df_simu_sigma2[iterN-1]+mu0/lambda0)*var_mu_k
        df_simu_mu[:,iterN] = np.random.normal(mean_mu_k,np.sqrt(var_mu_k))

        # simulate sigma2
        beta_IG = beta+0.5*np.sum((y-df_simu_mu[df_sim_z[:,iterN-1].astype(int)-1,iterN])**2)
        df_simu_sigma2[iterN] = 1/np.random.gamma(alpha_iv_gamma, 1/beta_IG) #alpha_IG = alpha+0.5*n = 42047.0

        # simulate zi
        df_sim_z[:,iterN] = df_sim_z[:,iterN-1]
        temp_z = df_sim_z[voxel_with_neighbor,iterN]
        for neigh_idx in range(len(neighbor_voxel_index)):
            ## pdf of y_i,mu_k,sigma**2
            z_i_j = temp_z[neighbor_voxel_index[neigh_idx]]
            unique_z, counts_z = np.unique(z_i_j, return_counts=True)
            classLabel = np.zeros(q)
            classLabel[(unique_z-1).astype(int)] = (2*counts_z*gamma/np.sum(counts_z))
            f_i_k = scipy.stats.norm.pdf(y[neigh_idx], df_simu_mu[:,iterN], np.sqrt(df_simu_sigma2[iterN])) * np.exp(classLabel)
            prob_label = f_i_k/np.sum(f_i_k)
            temp_z[neigh_idx] = np.random.choice(np.arange(q)+1, size=1,p=prob_label)
        df_sim_z[voxel_with_neighbor,iterN] = temp_z
    print("")
    sys.stdout.flush()
    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Done! ")
    sys.stdout.flush()

    if save_chain !=None:
        with open(save_chain, 'wb') as f:
            np.save(f, df_simu_mu)
            np.save(f, df_simu_sigma2)
            #np.save(f, df_sim_w)
            np.save(f, df_sim_z)

    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Start to generate heritability aware brain atlas using posterior samples ...")
    sys.stdout.flush()

    print("Number of burn in is "+str(burn_in))
    sys.stdout.flush()

    clusterGroup = scipy.stats.mode(df_sim_z[:,burn_in:],axis=1,nan_policy='omit')
    clusterGroup = clusterGroup[0].reshape(-1)

    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Done!")
    sys.stdout.flush()

    print("[ "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+" ] Start to save atlas to NIFTI file ...")
    sys.stdout.flush()
    atlas = (np.zeros(mask.shape).astype(int)).ravel()
    atlas[subsetL_idx] = clusterGroup
    atlas = atlas.reshape(mask.shape)
    atlas[atlas>0] = scipy.stats.rankdata(atlas[atlas>0],"dense")
    atlas = ReAssignRegion(network,atlas)
    atlas[atlas>0] = scipy.stats.rankdata(atlas[atlas>0],"dense").astype(int)
    # np.sum(np.take(atlas,subsetL_idx)==clusterGroup) # Uncomment to check if voxels are correctly assigned
    new_image = nib.Nifti1Image(atlas, affine)
    nib.save(new_image,save_nifti_path_prefix+"-clust-"+str(np.max(atlas[atlas>0]))+".nii.gz")

    return atlas

tissue_name = sys.argv[1] ## gm or wm
num_level = sys.argv[2]
num_gamma = sys.argv[3]
num_chain_length =  sys.argv[4]
num_burn_in = sys.argv[5]
init_path = sys.argv[6]

data_her = pd.read_csv("../example_data/voxel_heritability_UKBB/"+str(tissue_name)+"_heritability.csv", index_col=0)
y = data_her[['heritability']].values
img = y
img1 = nib.load("../example_data/conventional_atlas/convention_2x2x2mm3_gm_wm_mask.nii.gz")
data = img1.get_fdata()
if str(tissue_name)=="wm":
    mask = data==250
elif str(tissue_name)=="gm":
    mask = data==150

atlas = Her_Atlas(y,mask,img1.affine,init_path=init_path,
                  save_nifti_path_prefix="../example_data/example_output/"+str(tissue_name)+"_"+str(num_gamma),
                  level=int(num_level),
                  neighbor="sphere",
                  chain_length=int(num_chain_length),
                  burn_in=int(num_burn_in),
                  gamma=float(num_gamma),
                  save_chain="../example_data/example_output/MCMC_Chain_"+str(tissue_name)+"_"+str(num_gamma)+".npy"
                  )
