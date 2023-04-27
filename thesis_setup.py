def thesis_setup(Sitename):
    """
    Loads setup data for Thesis scripts
    
    Input:
        Sitename = Name of flux site {string}
            e.g. 'Young', 'Hogg', 'US-Myb', 'US-WPT'
    Returns:
        tower_coordinates {float} [longitude,latitude]
        L8_filename {string} Google Earth Engine output file name
        
        workPath {string} path for working Flux map python
        figurePath {string} folder path to store figures in
        dataPath {string} folder path to save data ouputs
        L8Path {string} path where Landsat 8 data from Google Earth Engine are saved
        FluxMapPath {string} path where Flux Map data are saved
        savePath {string} path to save runs to.
    """
    # Paths - Change as needed
    workPath = '/Users/darianng/Documents/Msc_Geography/Methane_Hotspot/FARF_Code/python_code'
    figurePath = '/Users/darianng/Documents/MSc_Geography/MSc Thesis/Figures'
    dataPath = '/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data'
    L8Path = '/Volumes/GoogleDrive/My Drive/Micromet_GEE'
    FluxMapPath = '/Users/darianng/Documents/MSc_Geography/Methane_Hotspot/FARF_Code/python_code/data'
    savePath = '/Users/darianng/Documents/Msc_Geography/Methane_Hotspot/FARF_Code/python_code/Saved_Data_Cluster'
    
    # Getting site data
    if Sitename == 'Hogg':
        tower_coordinates = [-100.534,50.371] # Hogg
        L8_filename = 'Hogg_spatial_indices_2021_May_2022_Nov.csv'
        
    elif Sitename == 'Young':
        tower_coordinates = [-100.20242,50.3623] # Young updated coordinates
        L8_filename = 'Young_spatial_indices_2021_May_2023_April_Large.csv'
        
    elif Sitename == 'US-Myb':
        tower_coordinates = [-121.7651,38.0498] # Mayberry
        L8_filename ='US-Myb_spatial_indices_2020_Jan_2022_June.csv'
        
    elif Sitename == 'US-WPT':
        tower_coordinates = [-82.9962,41.4646] # WPT
        L8_filename = 'US-WPT_spatial_indices_2012_Jan_2013_Dec.csv'
    else:
        print('Error: No data setup for tower yet.')
    
    return tower_coordinates,L8_filename,workPath,figurePath,dataPath,L8Path,FluxMapPath,savePath