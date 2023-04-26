def plot_thesis_dendrogram(fmap_x,fmap_y,fmap_flux):
    """
    Plots figures for dendrogram analysis
    Inputs:
        fmap_x = fmap_x (Footprint-weighted Flux Map data x-coordinates)
        fmap_y = fmap_y (Footprint-weighted Flux Map data y-coordinates)
        fmap_flux = matched_ffp['ch4'] (Footprint-weighted Flux Map flux values)
        * all input data have been coarsened to 30 m resolution
    Outputs:
        Figure a) Dendrogram
        Figure b) Footprint-weighted Flux Map plot
        Figure c) Clustered flux map
    Return:
        fig (subplot with Figures a,b,c)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram, linkage
    from sklearn.cluster import AgglomerativeClustering

    

    # DENDROGRAM
    fig = plt.figure(figsize=(8, 15))
    plt.subplot(3,1,1)
    plt.scatter(fmap_x,fmap_y,c=fmap_flux,marker = 's',s=350)
    plt.axvline(0,c='k',linestyle='--',alpha=0.7)
    plt.axhline(0,c='k',linestyle='--',alpha=0.7)
    plt.title('Footprint-weighted Flux Map',fontsize=20,fontweight='bold')
    plt.xlabel('Distance from tower [m]',fontsize=16)
    plt.ylabel('Distance from tower [m]',fontsize=16)
    plt.text(-162,137,'a)',fontsize = 22)

    plt.subplot(3,1,2)
    dend_data = np.transpose([np.unique(fmap_flux),np.unique(fmap_flux)])
    linked = linkage(dend_data,'ward')

    dendrogram(linked, 
            orientation='top', 
            distance_sort='descending',
            truncate_mode='lastp',
            p=30)
    # plt.axhline(0.2,c='k',linestyle='--',alpha=0.5)
    plt.title('Dendrogram',fontsize=20,fontweight='bold')
    plt.ylabel('distance',fontsize=16)
    plt.xlabel('Data clusters',fontsize=16)
    plt.text(7,805,'b)',fontsize = 22)

    plt.subplot(3,1,3)
    dend_data = np.transpose([fmap_flux,fmap_flux])
    linked = linkage(dend_data,'ward')
    #now cluster
    n_clusters = 5
    cluster = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage='ward')

    cn_PC = cluster.fit_predict(dend_data)
    #find mean pattern of each cluster
    cluster_pattern_PC = np.empty((n_clusters,np.shape(dend_data)[1]))

    # For 2 clusters
    for cluster_num in range(n_clusters):
        inds_PC = np.argwhere(cn_PC==cluster_num)
        cluster_pattern_PC[cluster_num,:] = np.mean(dend_data[inds_PC,:],axis=0)
    # Plotting    
    for i in range(len(np.unique(cn_PC))):
        clus_num = np.where(cn_PC==i)[0]
        plt.scatter(np.array(fmap_x)[clus_num],np.array(fmap_y)[clus_num],label=f'cluster {i+1}',marker = 's',s=350)
    plt.axvline(0,c='k',linestyle='--',alpha=0.7)
    plt.axhline(0,c='k',linestyle='--',alpha=0.7)
    plt.ylabel('metres')
    plt.xlabel('metres')
    plt.title('Footprint pixel by cluster',fontsize = 20, fontweight='bold')
    plt.legend(loc='upper right',fontsize='small',markerscale=0.5)
    plt.text(-170,148,'c)',fontsize = 22)

    plt.tight_layout()

    # plt.show()

    return fig

# ---------------------------------------------------------------------------------------------------------

def compiled_regression(FARF_list, correlation="Spearman"):
    import numpy as np
    import pickle 
    from scipy import stats
    import matplotlib.pyplot as plt
    import os

    """
    Compiles FARF analyses from /MSc_Geography/MSc Thesis/Data/Compilation

    Takes in: 
            1) List of cluster-completed FARF data. 
            E.g:
            ['Young_FFP=202105027-202106029_L8=20210527',
             'Young_FFP=202106007-202106021_L8=20210612',
             'Young_FFP=202107020-202108020_L8=20210808']

            2) Correlation approach. Either "Spearman" or "Pearson". Default is set to "Spearman"

    Outputs 2 figures:
        Fig 1) Subplot of each observation 
        Fig 2) Regression plot of all observations collapsed into a single dataset.
    Returns:
        Fig 1
        Fig 2
        compiled data
    
    """
    # Folder that all .p datafiles are stored in
    dataPath = '/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data/Compilation'

    indices = ['NDWI','NDMI','NDVI','temp']

    # Storing all runs into single dict, while keeping runs separate
    collected_regression = {'period':[],'r_val':[],'p_val':[],'landsat':[],'ffp':[]}

    # Collapsing all runs into single lists
    collapsed_landsat = {'NDWI':[],'NDMI':[],'NDVI':[],'temp':[]}
    collapsed_ffp = []

    # Storing all linear equation data:
    yfit = {'slope': [],'offset':[]}

    # Loading data
    for this_FARF in FARF_list:
        os.chdir(dataPath)
        with open(f'{this_FARF}.p', 'rb') as fp:
            data = pickle.load(fp)
        collected_regression['period'].append(this_FARF)
        collected_regression['landsat'].append(data['mean_clustered_landsat'])
        collected_regression['ffp'].append(data['mean_clustered_ffp'])
        collapsed_ffp.extend(data['mean_clustered_ffp'])
        for key in collapsed_landsat.keys():
            collapsed_landsat[key].extend(data['mean_clustered_landsat'][key])
    
    # Defining Spearman permutation test
    def statistic(x):  # permute only `x`
        return stats.spearmanr(x, y).statistic

    # Plotting figure 1: Subplot of each observation handled individually
    fig1 = plt.figure(figsize=(25,3*len(FARF_list)))
    plot_count = 1
    for ii in range(len(collected_regression['period'])):
        for kk,idx in enumerate(indices):
            plt.subplot(len(FARF_list),4,plot_count)
            
            # Plot data
            x = collected_regression['landsat'][ii][idx]
            y = collected_regression['ffp'][ii]

            # Calculating correlation coefficients
            if correlation == 'Pearson':
                r_val, p_val = stats.pearsonr(x,y)
            if correlation == 'Spearman':
                r_val, p_val = stats.spearmanr(x,y)
                perm = stats.permutation_test((x,), statistic, permutation_type='pairings')
                p_val = perm.pvalue
            corr = f'r: {np.round(r_val,6)} \np-value: {np.round(perm.pvalue,4)}'

            plt.scatter(x,y,label=corr)

            # Plotting regression line
            m, b = np.polyfit(x,y,1)
            yfit_line = m*np.array(x)+b
            plt.plot(x,yfit_line,'orange',label=f'Slope: {np.round(m,4)} \nOffset: {np.round(b,4)}')

            plt.ylabel('CH4 flux\n[nmol/m$^2$/s]',fontsize = 14)
            plt.xlabel(idx,fontsize = 16,fontweight='bold')
            plt.legend(loc='lower left', bbox_to_anchor=(1, 0.6))

            if kk%4==0: # Title only at the start of each row
                title_name = f"{collected_regression['period'][ii][-8:-4]}\n{collected_regression['period'][ii][-4:-2]}-{collected_regression['period'][ii][-2:]}"
                plt.annotate(title_name, xy=(0, 0.5), xycoords='axes fraction', xytext=(-0.6, 0.5),
                             textcoords='axes fraction', va='center', ha='center',
                             fontsize=20, fontweight='bold')
                
            # Storing regression data:
            collected_regression['p_val'].append(p_val)
            collected_regression['r_val'].append(r_val)
            # Storing linear equation data:
            yfit['slope'].append(m)
            yfit['offset'].append(b)
            
            plot_count += 1

    plt.tight_layout()

    # Plotting Figure 2: Compiled Regression
    fig2 = plt.figure(figsize=(15,3))
    for idx,key in enumerate(collapsed_landsat.keys()):
        plt.subplot(1,4,idx+1)
        plt.scatter(collapsed_landsat[key],collapsed_ffp,s=10)
        plt.title(key,fontsize=16,fontweight='bold')
        plt.ylabel('CH4 Flux  [µmol/m$^2$/s]')
        plt.xlabel(key)
    plt.tight_layout()

    compiled_data = {'collected_regression':collected_regression,'collapsed_landsat':collapsed_landsat,'collapsed_ffp':collapsed_ffp,'yfit':yfit}

    return fig1, fig2, compiled_data

def ffp_clustering(data,GHG_var='ch4',n_clusters = 5):
    """
    Performs Agglomerative Clustering with Ward's linkage.
    Input data: 
        data: data dict created from sector_plot.py (Relies on having to call sector_plot.py first before running this)
        GHG_var: flux variable to be clustered. Default set to 'ch4'
        n_clusters: number of cluster groups. Default set to 5
    return:
        clustered_landsat (dict holding 5-datapoint-clusters for each Landsat 8 product)
        clustered_ffp (dict holding 5-datapoint-cluster for ch4)
    """
    #now cluster
    from sklearn.cluster import AgglomerativeClustering
    import numpy as np

    dend_data = np.transpose([data['matched_ffp'][GHG_var],data['matched_ffp'][GHG_var]])

    # n_clusters = 5 # parameter passed to function input
    cluster = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage='ward')

    cn_PC = cluster.fit_predict(dend_data)
    #find mean pattern of each cluster
    cluster_pattern_PC = np.empty((n_clusters,np.shape(dend_data)[1]))

    # For 2 clusters
    for cluster_num in range(n_clusters):
        inds_PC = np.argwhere(cn_PC==cluster_num)
        cluster_pattern_PC[cluster_num,:] = np.mean(dend_data[inds_PC,:],axis=0)

    print(f'Confirming current Gas: {GHG_var}')

    # Compiling clustered data into dicts.
    clustered_ffp = {}
    clustered_landsat = {}
    for clus_idx in range(n_clusters):
        clus_num = np.where(cn_PC == clus_idx)[0]
        clustered_ffp[f'xr_{clus_idx+1}'] = np.array(matched_ffp['xr'])[clus_num]
        clustered_ffp[f'yr_{clus_idx+1}'] = np.array(matched_ffp['yr'])[clus_num]
        clustered_ffp[f'ch4_{clus_idx+1}'] = np.array(matched_ffp[GHG_var])[clus_num]
        clustered_landsat[f'NDVI_{clus_idx+1}'] = np.array(landsat['NDVI'])[clus_num]
        clustered_landsat[f'NDWI_{clus_idx+1}'] = np.array(landsat['NDWI'])[clus_num]
        clustered_landsat[f'NDMI_{clus_idx+1}'] = np.array(landsat['NDMI'])[clus_num]
        clustered_landsat[f'MNDWI2_{clus_idx+1}'] = np.array(landsat['MNDWI2'])[clus_num]
        clustered_landsat[f'temp_{clus_idx+1}'] = np.array(landsat['temp'])[clus_num]
    mean_clustered_ffp = []
    mean_clustered_landsat = {'NDVI':[],'NDWI':[],'NDMI':[],'MNDWI2':[],'temp':[]}
    for i in range(n_clusters):
        mean_clustered_ffp.append(np.mean(clustered_ffp[f'ch4_{i+1}']))
        for idx,key in enumerate(spatialData):
            mean_clustered_landsat[key].append(np.mean(clustered_landsat[f'{key}_{i+1}']))
    
    return mean_clustered_landsat, mean_clustered_ffp


def save_cluster(data,L8_filename,save_folder):
    """
    Saves clustered figures and data
    Input:
        data: data dict created from sector_plot.py (Relies on having to call sector_plot.py first before running this)
        L8_filename: same filename used for sector_plot.py; imported from Google Earth Engine
            e.g. 'Young_spatial_indices_2021_May_2023_April_Large.csv'
        save_folder: folder name.
    Return:
        fig1: Regression
        fig2: Clustering map
        data: Stored data in a .pickle file
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy as stats
    import pandas as pd

    # Getting clustered data
    from thesis_plots import ffp_clustering
    mean_clustered_landsat, mean_clustered_ffp = ffp_clustering(data=data,GHG_var='ch4',n_clusters=5)

    saved_indices = ['NDVI','NDWI','NDMI','temp']

    # Figure 1: Regression plot
    fig1 = plt.figure(figsize = (len(saved_indices)*4,4.5))
    for idx,key in enumerate(saved_indices):
        x_data = mean_clustered_landsat[key]
        y_data = mean_clustered_ffp

        plt.subplot(1,len(saved_indices),idx+1)

        r_val, p_val = stats.spearmanr(x_data,y_data)
        corr = f'r = {np.round(r_val,5)} \np = {np.round(p_val,5)}'

        plt.scatter(x_data,y_data,label=corr)
        plt.title(f'CH4 vs. {key}',fontweight='bold',fontsize=14)
        m, b = np.polyfit(x_data,y_data,1)
        yfit = m*np.array(x_data)+b
        plt.plot(x_data,yfit,'orange')
        
        plt.legend()
        plt.ylabel('CH4 flux [µmol/m$^2$/s]',fontsize=14)
        plt.xlabel(f'{key}',fontsize=14)
        
    fig1.tight_layout()

    # Figure 2: Cluster plot
    from thesis_plots import plot_thesis_dendrogram
    fig2 = plot_thesis_dendrogram(fmap_x = data['matched_ffp']['xr'],fmap_y = data['matched_ffp']['yr'],fmap_flux = data['matched_ffp']['ch4'])

    # Adding clustered data to final data
    data['mean_clustered_landsat'] = mean_clustered_landsat
    data['mean_clustered_ffp'] = mean_clustered_ffp
    data['corr'] =[r_val,p_val]
    return fig1, fig2, data
