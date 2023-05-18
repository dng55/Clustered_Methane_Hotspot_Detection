def preprocess_maps(siteName,fluxmap_filename,L8_date,ffp_datapath):
    """
    Map data preprocessing to create comparable datasets
    Step 1: Remove Landsat 8 data outside the flux footprint
    Step 2: Coarsen Flux Maps to 30 m resolution
    
    Inputs:
        Sitename = Name of flux site {string}
            e.g. 'Young', 'Hogg', 'US-Myb', 'US-WPT'
        fluxmap = Flux map output file name {string}
            e.g. '202207017-202208017.csv'
        L8_date = Landsat 8 image date
            e.g. 20220523
        ffp_datapath = path name to folder holding ffp outputs.
            
    Returns dict holding all step 1 & 2 processed data
    """
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
    import os
    # Loading paths and setup data
    from thesis_setup import thesis_setup
    coordinates,L8_filename,workPath,thesisPath,figurePath,dataPath,L8Path,FluxMapPath,savePath = thesis_setup(siteName)
    FluxMapPath = ffp_datapath

    # Step 1: Remove Landsat 8 data outside the flux footprint
    os.chdir(L8Path)
    satdata = pd.read_csv(L8_filename,delimiter = ',',header = 1) # Loading L8 data
    
    # List of remote sensing indices to be used in this analysis
    indices = ['MNDWI2','MNDWI','NDWI','NDVI','temp']
    
    # Landsat ID suffix corresponding to analysis date. Refer to script description for list of dates.
    ANALYSIS_DATE = L8_date

    os.chdir(workPath)
    # Importing sub-function that grabs Landsat 8 data
    from get_spatial import get_spatial
    spatialData = {'MNDWI2':[],'NDMI':[],'NDVI':[],'NDWI':[],'temp':[]}
    lonData, latData, spatialData['MNDWI2'] = get_spatial(ANALYSIS_DATE, 'MNDWI_SW2',satdata,'daily',coordinates)
    lonData, latData, spatialData['NDVI'] = get_spatial(ANALYSIS_DATE, 'NDVI',satdata,'daily',coordinates)
    lonData, latData, spatialData['NDWI'] = get_spatial(ANALYSIS_DATE, 'NDWI',satdata,'daily',coordinates)
    lonData, latData, spatialData['NDMI'] = get_spatial(ANALYSIS_DATE, 'MNDWI_SW1',satdata,'daily',coordinates)
    lonData, latData, spatialData['temp'] = get_spatial(ANALYSIS_DATE, 'CELSIUS',satdata,'daily',coordinates)
    
    # Loading fluxmap data from FluxMapPath
    os.chdir(FluxMapPath)
    ffp = {}
    ffp['xr'] = pd.read_csv(siteName+'_fluxMap_x_'+fluxmap_filename+'.csv',header = None) # x-coordinates
    ffp['yr'] = pd.read_csv(siteName+'_fluxMap_y_'+fluxmap_filename+'.csv',header = None) # y-coordinates
    ffp['co2'] = pd.read_csv(siteName+'_fluxMap_co2_'+fluxmap_filename+'.csv',header = None) # CO2 spatial data
    ffp['ch4'] = pd.read_csv(siteName+'_fluxMap_ch4_'+fluxmap_filename+'.csv',header = None) # CH4 spatial data
    ffp['h'] = pd.read_csv(siteName+'_fluxMap_h_'+fluxmap_filename+'.csv',header = None) # Sensible heat spatial data
    
    os.chdir(workPath)
    # Importing sub-function that cuts out the landsat pixels found in the flux footprint area
    from landsat_footprint import landsat_footprint
    
    # Storing landsat data in dict called landsat
    landsat = {'MNDWI2':[],'NDMI':[],'NDVI':[],'NDWI':[],'temp':[],'lonData':[],'latData':[]}
    
    # Each iteration of landsat_footprint sub-function deals with only 1 spatial index. Do this for each index of interest
    sat = landsat_footprint(lonData,latData,spatialData['MNDWI2'], ffp)
    landsat['MNDWI2'] = sat['spatialData']
    sat = landsat_footprint(lonData,latData,spatialData['NDMI'], ffp)
    landsat['NDMI'] = sat['spatialData']
    sat = landsat_footprint(lonData,latData,spatialData['NDWI'], ffp)
    landsat['NDWI'] = sat['spatialData']
    sat = landsat_footprint(lonData,latData,spatialData['NDVI'], ffp)
    landsat['NDVI'] = sat['spatialData']
    sat = landsat_footprint(lonData,latData,spatialData['temp'], ffp)
    landsat['temp'] = sat['spatialData']
    landsat['lonData'] = sat['lonData']
    landsat['latData'] = sat['latData']
    
    #Matching FFP resolution to landsat resolution
    from landsat_footprint import ffp_matched_to_landsat
    matched_ffp = ffp_matched_to_landsat(landsat,ffp) # dict keys are the same for "matched_ffp" as for "landsat"
    
    allData = {'siteName':siteName,'landsat':landsat, 'matched_ffp':matched_ffp, 'ffp':ffp,'lonData':lonData,'latData':latData,
                         'spatialData':spatialData,'run':[L8_date,fluxmap_filename]}
    
    return allData

# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------
def ffp_clustering(data,GHG_var='ch4',n_clusters = 5):
    """
    Performs Agglomerative Clustering with Ward's linkage.
    Input data: 
        data: data dict created from preprocessing_maps function in thesis.py (Perform preprocessing first!)
        GHG_var: flux variable to be clustered. Default set to 'ch4'
        n_clusters: number of cluster groups. Default set to 5
    return:
        data: data dict from input but with additions:
            'clustered_landsat' = list holding 5-datapoint-clusters for each Landsat 8 product)
            'clustered_ffp' = list holding 5-datapoint-cluster for ch4)
            'GHG_var' = Gas flux being measured (use to make sure axis labels are correct)
        fig = Linear regression plot with all available L8 products.
    """
    #now cluster
    from sklearn.cluster import AgglomerativeClustering
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats


    dend_data = np.transpose([data['matched_ffp'][GHG_var],data['matched_ffp'][GHG_var]])

    # n_clusters = 5 # parameter passed to function input
    # cluster = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage='ward')
    cluster = AgglomerativeClustering(n_clusters=n_clusters, metric='euclidean', linkage='ward')

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
        clustered_ffp[f'xr_{clus_idx+1}'] = np.array(data['matched_ffp']['xr'])[clus_num]
        clustered_ffp[f'yr_{clus_idx+1}'] = np.array(data['matched_ffp']['yr'])[clus_num]
        clustered_ffp[f'ch4_{clus_idx+1}'] = np.array(data['matched_ffp'][GHG_var])[clus_num]
        clustered_landsat[f'NDVI_{clus_idx+1}'] = np.array(data['landsat']['NDVI'])[clus_num]
        clustered_landsat[f'NDWI_{clus_idx+1}'] = np.array(data['landsat']['NDWI'])[clus_num]
        clustered_landsat[f'NDMI_{clus_idx+1}'] = np.array(data['landsat']['NDMI'])[clus_num]
        clustered_landsat[f'MNDWI2_{clus_idx+1}'] = np.array(data['landsat']['MNDWI2'])[clus_num]
        clustered_landsat[f'temp_{clus_idx+1}'] = np.array(data['landsat']['temp'])[clus_num]
    mean_clustered_ffp = []
    mean_clustered_landsat = {'NDVI':[],'NDWI':[],'NDMI':[],'MNDWI2':[],'temp':[]}
    for i in range(n_clusters):
        mean_clustered_ffp.append(np.mean(clustered_ffp[f'ch4_{i+1}']))
        for idx,key in enumerate(data['spatialData']):
            mean_clustered_landsat[key].append(np.mean(clustered_landsat[f'{key}_{i+1}']))
    
    # Fig: Calculating and plotting linear regression

    # Defining Spearman permutation test
    def statistic(x_data):  # permute only `x`
        return stats.spearmanr(x_data, y_data).statistic
    
    L8_products = ['NDVI','NDWI','NDMI','temp','MNDWI2']

    fig = plt.figure(figsize = (24,4))
    for idx,key in enumerate(L8_products):
        x_data = mean_clustered_landsat[key]
        y_data = mean_clustered_ffp

        plt.subplot(1,5,idx+1)
        # Plotting linear regression line
        m, b = np.polyfit(x_data,y_data,1)
        yfit = m*np.array(x_data)+b
        plt.plot(x_data,yfit,'orange',label=f'Slope: {np.round(m,4)} \nOffset: {np.round(b,4)}')

        # Calculating correlation coefficients with permuted Spearman
        r_val, p_val = stats.spearmanr(x_data,y_data)
        perm = stats.permutation_test((x_data,), statistic, permutation_type='pairings')
        p_val = perm.pvalue
        corr = f'r: {np.round(r_val,6)} \np-value: {np.round(perm.pvalue,4)}'

        plt.scatter(x_data,y_data,label = corr)
        plt.title(f'{key}',fontweight='bold',fontsize=18)
        plt.xlabel(f'{key}',fontsize=14)
        if GHG_var == 'ch4':
            plt.ylabel(f'{GHG_var} flux \n[nmol/m$^2$/s]',fontsize=14)
        elif GHG_var == 'co2':
            plt.ylabel(f'{GHG_var} flux \n[µmol/m$^2$/s]',fontsize=14)
        plt.legend()
        
    plt.tight_layout()
    # plt.suptitle(f'Landsat image:\n{data["run"][0]}\n\nFlux Map period:\n{data["run"][1]}',y=1.45,fontsize=20,fontweight = 'bold')
    plt.suptitle(f'Landsat image:              Flux Map period:\n          {data["run"][0]}              {data["run"][1]}',y=1.2,fontsize=20,fontweight = 'bold')

    data['clustered_landsat'] = mean_clustered_landsat
    data['clustered_ffp'] = mean_clustered_ffp
    data['GHG_var'] = GHG_var
    return data, fig

# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------
def save_cluster(data,save_folder):
    """
    Saves clustered figures and data to save_folder path.
    Regression statistics for each L8 index are added to data, and then saved as a .pickle.
    Input:
        data: data dict created from {preprocessing_maps} and {ffp_clustering} functions in thesis.py (Preprocessing & cluster first!)
        L8_filename: same filename used for sector_plot.py; imported from Google Earth Engine
            e.g. 'Young_spatial_indices_2021_May_2023_April_Large.csv'
        save_folder: folder name.
    Return:
        fig1: Regression (fig2, cluster plot, is saved, but not returned.)
        data: Stored data in a .pickle file
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
    import pandas as pd
    import pickle
    import os
    from thesis import plot_thesis_dendrogram

    mean_clustered_landsat = data['clustered_landsat']
    mean_clustered_ffp = data['clustered_ffp']

    saved_indices = ['NDVI','NDWI','NDMI','temp']

    # Creating dict key to store regression data
    data['regression'] = {'NDVI':{'rval':[],'pval':[],'slope':[],'offset':[]},
                          'NDWI':{'rval':[],'pval':[],'slope':[],'offset':[]},
                          'NDMI':{'rval':[],'pval':[],'slope':[],'offset':[]},
                          'temp':{'rval':[],'pval':[],'slope':[],'offset':[]}}

    # Figure 1: Regression plot
    fig1 = plt.figure(figsize = (len(saved_indices)*6,4))

    def statistic(x_data):  # Defining Spearman permutation test
        return stats.spearmanr(x_data, y_data).statistic

    for idx,key in enumerate(saved_indices):
        x_data = data['clustered_landsat'][key]
        y_data = data['clustered_ffp']
        GHG_var = data['GHG_var']

        plt.subplot(1,len(saved_indices),idx+1)
        
        # Plotting linear regression line
        m, b = np.polyfit(x_data,y_data,1)
        yfit = m*np.array(x_data)+b
        plt.plot(x_data,yfit,'orange',label=f'Slope: {np.round(m,4)} \nOffset: {np.round(b,4)}')
        
        # Calculating correlation coefficients with permuted Spearman
        r_val, p_val = stats.spearmanr(x_data,y_data)
        perm = stats.permutation_test((x_data,), statistic, permutation_type='pairings')
        p_val = perm.pvalue
        corr = f'r: {np.round(r_val,6)} \np-value: {np.round(perm.pvalue,4)}'

        plt.scatter(x_data,y_data,label = corr)
        plt.title(f'CH$_4$ vs. {key}',fontweight='bold',fontsize=14)

        plt.legend(loc='lower left', bbox_to_anchor=(1, 0.7))
        plt.xlabel(f'{key}',fontsize=14)
        
        if GHG_var == 'ch4':
            plt.ylabel(f'CH$_4$ flux [nmol/m$^2$/s]',fontsize=14)
        elif GHG_var == 'co2':
            plt.ylabel(f'CO$_2$ flux [µmol/m$^2$/s]',fontsize=14)
        
        # Adding regression statistics to data
        data['regression'][key]['rval'].append(r_val)
        data['regression'][key]['pval'].append(p_val)
        data['regression'][key]['slope'].append(m)
        data['regression'][key]['offset'].append(b)
        
    fig1.tight_layout()

    # Figure 2: Cluster plot
    fig2 = plot_thesis_dendrogram(fmap_x = data['matched_ffp']['xr'],fmap_y = data['matched_ffp']['yr'],fmap_flux = data['matched_ffp']['ch4'])
    # fig2.clf()

    # Saving figures and data to save_folder path
    save_filename = f'{data["siteName"]}_FFP={data["run"][1]}_L8={str(data["run"][0])}'
    save_filepath = f'{save_folder}/{save_filename}'
    # Making sub-folder for this run if it does not already exist
    if not os.path.exists(save_filepath):
        os.makedirs(save_filepath)
    # Saving figures
    fig1.savefig(f'{save_filepath}/Fig1_Regression.png')
    fig2.savefig(f'{save_filepath}/Fig2_Cluster.png')
    fig1.savefig(f'{save_filepath}.png') # Saving this figure outside the folder to easily see run results
    # Saving datafile
    with open(save_filepath+'/'+save_filename+'.p', 'wb') as fp:
        pickle.dump(data, fp, protocol=pickle.HIGHEST_PROTOCOL)

    print(f'Figures and data saved to: \n{save_filepath}')

    return fig1, data

# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------
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
    # plt.text(-162,137,'a)',fontsize = 22)

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
    # plt.text(7,805,'b)',fontsize = 22)

    plt.subplot(3,1,3)
    dend_data = np.transpose([fmap_flux,fmap_flux])
    linked = linkage(dend_data,'ward')
    #now cluster
    n_clusters = 5
    cluster = AgglomerativeClustering(n_clusters=n_clusters, metric='euclidean', linkage='ward')

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
    # plt.text(-170,148,'c)',fontsize = 22)

    plt.tight_layout()

    # plt.show()

    return fig

# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------
def compiled_regression(FARF_list, dataPath, correlation="Spearman"):
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

            2) dataPath = folder path of saved data for compilation
            3) Correlation approach. Either "Spearman" or "Pearson". Default is set to "Spearman"

    Outputs 2 figures:
        Fig 1) Subplot of each observation 
        Fig 2) Regression plot of all observations collapsed into a single dataset.
    Returns:
        Fig 1
        Fig 2
        compiled data
    
    """
    # Folder that all .p datafiles are stored in

    indices = ['NDWI','NDMI','NDVI','temp']

    # Storing all runs into single dict, while keeping runs separate
    collected_regression = {'period':[],'r_val':[],'p_val':[],'landsat':[],'ffp':[]}

    collapsed_regression = {'slope':[],'offset':[],'r_val':[],'p_val':[]}

    # Collapsing all runs into single analysis
    collapsed_landsat = {'NDWI':[],'NDMI':[],'NDVI':[],'temp':[],'period':[]}
    collapsed_ffp = []

    # Storing all linear equation data:
    yfit = {'slope': [],'offset':[]}

    print('Collecting data...', end='')

    # Loading data
    for this_FARF in FARF_list:
        os.chdir(dataPath)
        with open(f'{this_FARF}.p', 'rb') as fp:
            data = pickle.load(fp)
        # collected_regression for aggregating observation-separated data
        collected_regression['period'].append(this_FARF)
        collected_regression['landsat'].append(data['clustered_landsat'])
        collected_regression['ffp'].append(data['clustered_ffp'])
        collapsed_ffp.extend(data['clustered_ffp'])

        # collapsed_landsat for all across-observation data joined together
        for key in list(collapsed_landsat.keys())[:-1]: # Skipping the last key which is 'period'
            collapsed_landsat[key].extend(data['clustered_landsat'][key])
        # Giving each datapoint a date label.
        this_period = [int(this_FARF.split('L8=')[-1])]*len(data['clustered_landsat']['NDVI']) # Arbitratily calling NDVI just to get the length
        collapsed_landsat['period'].extend(this_period) 
    
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

            # Plotting regression line (doing this first to prioritize order in legend)
            m, b = np.polyfit(x,y,1)
            yfit_line = m*np.array(x)+b
            plt.plot(x,yfit_line,'orange',label=f'Slope: {np.round(m,4)} \nOffset: {np.round(b,4)}')

            # Calculating correlation coefficients
            if correlation == 'Pearson':
                r_val, p_val = stats.pearsonr(x,y)
            elif correlation == 'Spearman':
                r_val, p_val = stats.spearmanr(x,y)
                perm = stats.permutation_test((x,), statistic, permutation_type='pairings')
                p_val = perm.pvalue
            corr = f'r: {np.round(r_val,6)} \np-value: {np.round(p_val,4)}'

            plt.scatter(x,y,label=corr)

            plt.ylabel('CH4 flux\n[nmol/m$^2$/s]',fontsize = 14)
            plt.xlabel(idx,fontsize = 16)
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
    print('done!')
    print('Collapsing data to single analysis...', end = '')

    # Plotting Figure 2: Compiled Regression
    fig2 = plt.figure(figsize=(24,4))
    for idx,key in enumerate(list(collapsed_landsat.keys())[:-1]):
        plt.subplot(1,4,idx+1)

        # Plot data
        x = collapsed_landsat[key]
        y = collapsed_ffp

        # Plotting regression line (doing this first to prioritize order in legend)
        m, b = np.polyfit(x,y,1)
        yfit_line = m*np.array(x)+b
        plt.plot(x,yfit_line,'orange',label=f'Slope: {np.round(m,4)} \nOffset: {np.round(b,4)}')

        # Regression stats
        r_val, p_val = stats.spearmanr(x,y)
        perm = stats.permutation_test((x,), statistic, permutation_type='pairings')
        p_val = perm.pvalue
        corr = f'r: {np.round(r_val,4)} \np-value: {np.round(perm.pvalue,4)}'
        
        plt.scatter(x,y,label=corr)
        
        plt.title(key,fontsize=18,fontweight='bold')
        plt.legend(loc='lower left', bbox_to_anchor=(1, 0.7))
        plt.ylabel('CH4 Flux  [nmol/m$^2$/s]',fontsize = 14)
        plt.xlabel(key,fontsize = 14)

        # Storing regression data:
        collapsed_regression['slope'].append(m)
        collapsed_regression['offset'].append(b)
        collapsed_regression['r_val'].append(r_val)
        collapsed_regression['p_val'].append(p_val)
        if idx == 1:
            print('half way there...',end='')
    plt.tight_layout()
    print('done!')

    compiled_data = {'collected_regression':collected_regression,'collapsed_regression':collapsed_regression,'collapsed_landsat':collapsed_landsat,'collapsed_ffp':collapsed_ffp,'yfit':yfit}

    return fig1, fig2, compiled_data

# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------

def check_L8(site,L8_date,product='NDVI',pixelSize=100):
    """
    Plots raw L8 map exported from Google Earth Engine
    Inputs:
        site: 'Young', 'Hogg', 'US-Myb', or 'US-WPT'
        L8_date: date of L8 image (e.g. 20210527)
        product: remote sensing product to be plotted. Default set to NDVI
        pixelSize: size of plotted pixel. Default set to 100.
    Outputs:
        Plot of L8 map.
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from thesis_setup import thesis_setup
    import os

    # Importing paths
    tower_coordinates,L8_filename,workPath,thesisPath,figurePath,dataPath,L8Path,FluxMapPath,savePath = thesis_setup(site)

    # Loading L8 data
    os.chdir(L8Path)
    satdata = pd.read_csv(L8_filename,delimiter = ',',header = 1)

    N = len(satdata['id'])

    imageCollection = {'latitude':[],'longitude':[],'NDVI':[],'NDWI':[],'MNDWI_SW1':[],'MNDWI_SW2':[],'CELSIUS':[]}

    for i in range(N):
        id_date = satdata['id'][i][12:]
        if id_date == str(L8_date):
            for key in imageCollection.keys():
                imageCollection[key].append(satdata[key][i])
    plt.scatter(imageCollection['longitude'],imageCollection['latitude'],c=imageCollection[product],marker='s',s=pixelSize)
    xzoom = -0.009
    yzoom = -0.005
    # plt.xlim([min(imageCollection['longitude'])-xzoom,max(imageCollection['longitude'])+xzoom])
    # plt.ylim([min(imageCollection['latitude'])-yzoom,max(imageCollection['latitude'])+yzoom])
    plt.axvline(tower_coordinates[0],c='k')
    plt.axhline(tower_coordinates[1],c='k')
    os.chdir(thesisPath)

