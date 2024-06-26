import numpy as np
import time
import astropy.constants as const
import astropy.units as u

def galdtype_darksage(Nannuli=30):
    floattype = np.float32
    Galdesc_full = [
                    ('Galaxy_Classification'        , np.int32),
                    ('GalaxyIndex'                  , np.int64),
                    ('HaloIndex'                    , np.int32),
                    ('SimulationHaloIndex'          , np.int32),
                    ('TreeIndex'                    , np.int32),
                    ('SnapNum'                      , np.int32),
                    ('CentralGalaxyIndex'           , np.int64),
                    ('Central_Galaxy_Mvir'          , floattype),
                    ('mergeType'                    , np.int32),
                    ('mergeIntoID'                  , np.int32),
                    ('mergeIntoSnapNum'             , np.int32),
                    ('dT'                           , floattype),
                    ('Pos'                          , (floattype, 3)),
                    ('Vel'                          , (floattype, 3)),
                    ('halospin'                         , (floattype, 3)),
                    ('Len'                          , np.int32),
                    ('LenMax'                       , np.int32),
                    ('Mvir'                         , floattype),
                    ('Rvir'                         , floattype),
                    ('Vvir'                         , floattype),
                    ('Vmax'                         , floattype),
                    ('VelDisp'                      , floattype),                # velocity dispersion of dark matter only
                    ('DiscRadii'                    , (floattype, Nannuli+1)),   # radius of each annulus for every galaxy (every galaxy has 31 radii that pertain to each annulus)
                    ('Cold_Gas_Mass'                      , floattype),
                    ('Total_Stellar_Mass'                  , floattype),
                    ('Mergerdriven_Bulge_Mass'              , floattype),
                    ('Instabilitydriven_Bulge_Mass'          , floattype),
                    ('Hot_Gas_Mass'                       , floattype),
                    ('Ejected_Gas_Mass'                  , floattype),
                    ('Black_Hole_Mass'                , floattype),
                    ('IntraCluster_Stars'            , floattype),
                    ('DiscGas'                      , (floattype, Nannuli)),
                    ('DiscStars'                    , (floattype, Nannuli)),
                    ('j_Stellar_Disk'                    , (floattype, 3)),
                    ('j_Cold_Gas'                      , (floattype, 3)),
                    ('SpinClassicalBulge'           , (floattype, 3)),
                    ('StarsInSitu'                  , floattype),
                    ('StarsInstability'             , floattype),
                    ('StarsMergeBurst'              , floattype),
                    ('DiscHI'                       , (floattype, Nannuli)),
                    ('DiscH2'                       , (floattype, Nannuli)),
                    ('DiscSFR'                      , (floattype, Nannuli)),
                    ('Metals_Cold_Gas'                , floattype),
                    ('Metals_Stellar_Mass'            , floattype),
                    ('Classical_Metals_Bulge_Mass'     , floattype),
                    ('Secular_Metals_Bulge_Mass'       , floattype),
                    ('Metals_Hot_Gas'                 , floattype),
                    ('Metals_Ejected_Mass'            , floattype),
                    ('Metals_IntraCluster_Stars'      , floattype),
                    ('DiscGasMetals'                , (floattype, Nannuli)),
                    ('DiscStarsMetals'              , (floattype, Nannuli)),
                    ('SfrFromH2'                    , floattype),
                    ('SfrInstab'                    , floattype),
                    ('SfrMergeBurst'                , floattype),
                    ('SfrDiskZ'                     , floattype),
                    ('SfrBulgeZ'                    , floattype),
                    ('Disk_Scale_Radius'              , floattype),
                    ('Cool_Scale_Radius'              , floattype),
                    ('Stellar_Disc_Scale_Radius'       , floattype),
                    ('Cooling'                      , floattype),
                    ('Heating'                      , floattype),
                    ('QuasarEnergy'                      , floattype),
                    ('BHaccreted'                       , floattype),
                    ('BondiBHaccreted'                       , floattype),
                    ('RadioBHaccreted'                      , floattype),
                    ('QuasarBHaccreted'                      , floattype),
                    ('InstaBHaccreted'                      , floattype),
                    ('MergerBHaccreted'                      , floattype),
                    ('RadioBlackHoleMass'                      , floattype),
                    ('QuasarBlackHoleMass'                      , floattype),
                    ('InstaBlackHoleMass'                      , floattype),
                    ('MergerBlackHoleMass'                      , floattype),
                    ('TimeofFirstAccretionRadio'    ,  floattype),
                    ('TimeofLastAccretionRadio'    ,  floattype),
                    ('TimeofFirstAccretionQuasar'    ,  floattype),
                    ('TimeofLastAccretionQuasar'    ,  floattype),
                    ('TimeofFirstAccretionInsta'    ,  floattype),
                    ('TimeofLastAccretionInsta'    ,  floattype),
                    ('TimeofFirstAccretionMerger'    ,  floattype),
                    ('TimeofLastAccretionMerger'    ,  floattype),
                    ('Last_Major_Merger'              , floattype),
                    ('Last_Minor_Merger'              , floattype),
                    ('OutflowRate'                  , floattype),
                    ('Subhalo_Mvir_at_Infall'                   , floattype),
                    ('Subhalo_Vvir_at_Infall'                   , floattype),
                    ('Subhalo_Vmax_at_Infall'                   , floattype)
                    ]
    names = [Galdesc_full[i][0] for i in range(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in range(len(Galdesc_full))]
    Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)
    return Galdesc

"""
field = [ 'Galaxy_Classification', 'GalaxyIndex', 'HaloIndex', 'Pos', 'Vel', 'Mvir', 'Rvir',
    'Vvir', 'Vmax', 'Central_Galaxy_Mvir', 'VelDisp', 'DiscRadii', 'Cold_Gas_Mass', 'Total_Stellar_Mass',
    'Mergerdriven_Bulge_Mass', 'Instabilitydriven_Bulge_Mass', 'Hot_Gas_Mass', 'Ejected_Gas_Mass',
    'Black_Hole_Mass', 'j_Stellar_Disk', 'j_Cold_Gas', 'DiscStars', 'DiscHI',
    'DiscH2', 'DiscSFR', 'SfrFromH2', 'SfrInstab', 'SfrMergeBurst', 'Disk_Scale_Radius', 'Cool_Scale_Radius',
    'Cooling', 'Heating', 'QuasarEnergy', 'BHaccreted', 'BondiBHaccreted', 'RadioBHaccreted', 'QuasarBHaccreted', 'InstaBHaccreted',
    'MergerBHaccreted','RadioBlackHoleMass', 'QuasarBlackHoleMass', 'InstaBlackHoleMass', 'MergerBlackHoleMass',
    'TimeofFirstAccretionRadio', 'TimeofLastAccretionRadio', 'TimeofFirstAccretionQuasar', 'TimeofLastAccretionQuasar',
    'TimeofFirstAccretionInsta', 'TimeofLastAccretionInsta', 'TimeofFirstAccretionMerger', 'TimeofLastAccretionMerger',
    'Last_Major_Merger', 'Last_Minor_Merger', 'OutflowRate', 'Subhalo_Mvir_at_Infall', 'Subhalo_Vvir_at_Infall',
    'Subhalo_Vmax_at_Infall', 'SnapNum' ]
"""

def galdtype_darksage_old(Nannuli=30):
    floattype = np.float32
    Galdesc_full = [
                    ('Galaxy_Classification'                         , np.int32),
                    ('GalaxyIndex'                  , np.int64),
                    ('HaloIndex'                    , np.int32),
                    ('SimulationHaloIndex'          , np.int32),
                    ('TreeIndex'                    , np.int32),
                    ('SnapNum'                      , np.int32),
                    ('CentralGalaxyIndex'           , np.int64),
                    ('Central_Galaxy_Mvir'                  , floattype),
                    ('mergeType'                    , np.int32),
                    ('mergeIntoID'                  , np.int32),
                    ('mergeIntoSnapNum'             , np.int32),
                    ('dT'                           , floattype), 
                    ('Pos'                          , (floattype, 3)),
                    ('Vel'                          , (floattype, 3)),
                    ('halospin'                         , (floattype, 3)),
                    ('Len'                          , np.int32),
                    ('LenMax'                       , np.int32),
                    ('Mvir'                         , floattype),
                    ('Rvir'                         , floattype),
                    ('Vvir'                         , floattype),
                    ('Vmax'                         , floattype),
                    ('VelDisp'                      , floattype),                # velocity dispersion of dark matter only
                    ('DiscRadii'                    , (floattype, Nannuli+1)),   # radius of each annulus for every galaxy (every galaxy has 31 radii that pertain to each annulus)
                    ('Cold_Gas_Mass'                      , floattype),
                    ('Total_Stellar_Mass'                  , floattype),
                    ('Mergerdriven_Bulge_Mass'              , floattype),
                    ('Instabilitydriven_Bulge_Mass'          , floattype),
                    ('Hot_Gas_Mass'                       , floattype),
                    ('Ejected_Gas_Mass'                  , floattype),
                    ('Black_Hole_Mass'                , floattype),
                    ('IntraCluster_Stars'            , floattype),
                    ('DiscGas'                      , (floattype, Nannuli)),
                    ('DiscStars'                    , (floattype, Nannuli)),
                    ('j_Stellar_Disk'                    , (floattype, 3)),
                    ('j_Cold_Gas'                      , (floattype, 3)),
                    ('SpinClassicalBulge'           , (floattype, 3)),
                    ('StarsInSitu'                  , floattype),
                    ('StarsInstability'             , floattype),
                    ('StarsMergeBurst'              , floattype),
                    ('DiscHI'                       , (floattype, Nannuli)),
                    ('DiscH2'                       , (floattype, Nannuli)),
                    ('DiscSFR'                      , (floattype, Nannuli)), 
                    ('Metals_Cold_Gas'                , floattype),
                    ('Metals_Stellar_Mass'            , floattype),
                    ('Classical_Metals_Bulge_Mass'     , floattype),
                    ('Secular_Metals_Bulge_Mass'       , floattype),
                    ('Metals_Hot_Gas'                 , floattype),
                    ('Metals_Ejected_Mass'            , floattype),
                    ('Metals_IntraCluster_Stars'      , floattype),
                    ('DiscGasMetals'                , (floattype, Nannuli)),
                    ('DiscStarsMetals'              , (floattype, Nannuli)),
                    ('SfrFromH2'                    , floattype),
                    ('SfrInstab'                    , floattype),
                    ('SfrMergeBurst'                , floattype),
                    ('SfrDiskZ'                     , floattype),
                    ('SfrBulgeZ'                    , floattype),
                    ('Disk_Scale_Radius'              , floattype),
                    ('Cool_Scale_Radius'              , floattype), 
                    ('Stellar_Disc_Scale_Radius'       , floattype),
                    ('Cooling'                      , floattype),
                    ('Heating'                      , floattype),
                    ('Last_Major_Merger'              , floattype),
                    ('Last_Minor_Merger'              , floattype),
                    ('OutflowRate'                  , floattype),
                    ('Subhalo_Mvir_at_Infall'                   , floattype),
                    ('Subhalo_Vvir_at_Infall'                   , floattype),
                    ('Subhalo_Vmax_at_Infall'                   , floattype)
                    ]
    names = [Galdesc_full[i][0] for i in range(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in range(len(Galdesc_full))]
    Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)
    return Galdesc

def darksage_snap(fpre, filelist, fields=[], Nannuli=30):
    # Read full Dark Sage snapshot, going through each file and compiling into 1 array
    # fpre is the name of the file up until the _ before the file number
    # filelist contains all the file numbers you want to read in
    
    Galdesc = galdtype_darksage()
    if len(fields)==0: fields=list(Galdesc.names)
    NtotGalsSum = 0
    
    # First calculate the total number of galaxies that will fill the array
    for i in filelist:
        fin = open(fpre+'_'+str(i), 'rb')
        Ntrees = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file
        NtotGals = np.fromfile(fin,np.dtype(np.int32),1)[0]
        NtotGalsSum += NtotGals
        fin.close()

    G = np.empty(NtotGalsSum, dtype=Galdesc)[fields] # Intialise the galaxy array
    NtotGalsSum = 0 # reset for next loop

    # Loop through files to fill in galaxy array
    for i in filelist:
        print('reading file', i)
        fin = open(fpre+'_'+str(i), 'rb')
        Ntrees = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file
        NtotGals = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.
        GalsPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1) # Read the number of gals in each tree
        G1 = np.fromfile(fin, Galdesc, NtotGals) # Read all the galaxy data
        fin.close()
        G[NtotGalsSum:NtotGalsSum+NtotGals] = G1[fields]
        NtotGalsSum += NtotGals

    return G


def galaxy_sample( sample, data ):

    import pandas as pd

    gal = {}
    gal_dic = {}

    # Here we choose whether our sample should contain all galaxies, central galaxies only, or satellite galaxies only
    start = time.time()
    print('Loading data...')
    
    key=0
    for key in data.dtype.names:
        gal_dic[key] = list( data[key] )

    gal = pd.DataFrame( gal_dic )
    print('Time taking turning dataset to pandas dataframe: %f s.' % (time.time() - start) )

    if (sample == "centrals"):
        gal = gal.loc[ ( gal.Galaxy_Classification == 0 ) ]
        print( "You have chosen a sample of " + sample + "galaxies." )

    elif (sample == "satellites"):
        gal = gal.loc[ ( gal.Galaxy_Classification == 1 ) ]
        print( "You have chosen a sample of " + sample + "galaxies.")  

    else:
        print( "You have chosen a sample of " + sample + "galaxies.")

    return gal

def param_dicts(sample, data ):

    # Defining parameters
    """
    See text file for units.
    For the following, 1e10 has been applied and h has already been divided:
    Total Stellar Mass
    Bulge Stellar Mass
    Disk Stellar Mass
    HI Mass
    H2 Mass
    Mvir
    Cold Gas Mass
    Hot Gas Mass
    Black_Hole_Mass
    """

    h = 0.7
    gal = galaxy_sample( sample, data )
    Vvir = gal.Vvir   
    Rvir = gal.Rvir

    gal['Bulge_Stellar_Mass'] = (gal.Mergerdriven_Bulge_Mass + 
        gal.Instabilitydriven_Bulge_Mass) / h
    gal['Disk_Stellar_Mass'] = (gal.Total_Stellar_Mass - 
        gal.Bulge_Stellar_Mass) / h

    gal = gal.reset_index(drop=True)

    ## Summing up annuli properties
    start = time.time()
    print('Summing up annuli properties...')
    Total_Star_Formation_Rate = []
    HI_Mass = []
    H2_Mass = []
    for i in range(len(gal.DiscSFR)):
        Total_DiscSFR = np.sum(gal.DiscSFR[i])
        Total_HI_Mass = np.sum(gal.DiscHI[i])
        Total_H2_Mass = np.sum(gal.DiscH2[i])

        Total_Star_Formation_Rate.append(Total_DiscSFR)
        HI_Mass.append(Total_HI_Mass)
        H2_Mass.append(Total_H2_Mass)

        
    gal['Total_Star_Formation_Rate'] = Total_Star_Formation_Rate
    gal['HI_Mass'] = HI_Mass
    gal['H2_Mass'] = H2_Mass
    print('Time: %f s.' % (time.time() - start) )

    # delete columns we don't need anymore to save space
    gal = gal.drop(['DiscSFR', 'DiscHI', 'DiscH2', 'DiscRadii', 'DiscStars'], axis=1)

    ##### Adjustments. Multiplying by 1e10 for some parameters
    start = time.time()
    print('Lastly, converting to log form and multiplying by 1e10s...')

    gal.loc[( gal.Total_Stellar_Mass == 0.0 ), 'Total_Stellar_Mass'] = 1e-10

    gal.Total_Stellar_Mass = 1e10 * gal.Total_Stellar_Mass / h
    gal.Disk_Stellar_Mass = 1e10 * gal.Disk_Stellar_Mass
    gal.HI_Mass = 1e10 * gal.HI_Mass / h
    gal.H2_Mass = 1e10 * gal.H2_Mass / h

    gal.Mvir = 1e10 * gal.Mvir / h
    gal.Cold_Gas_Mass = 1e10 * gal.Cold_Gas_Mass / h
    gal.Hot_Gas_Mass = 1e10 * gal.Hot_Gas_Mass / h


    # Taking the log10 of these values
    gal['logSFR'] = np.log10( gal.Total_Star_Formation_Rate )
    gal['logsSFR'] = np.log10( gal.Total_Star_Formation_Rate / gal.Total_Stellar_Mass) 
    # this is to have all logsSFR values less than -12 to just be -12.0
    gal.loc[( gal.logsSFR < -12.0 ), 'logsSFR'] = -12.0    

    # Calculating morphology
    gal['morph'] = gal.Disk_Stellar_Mass / gal.Total_Stellar_Mass

    # Cold gas fraction
    gal['CGM_SM'] = gal.Cold_Gas_Mass / gal.Total_Stellar_Mass

    # Cold baryonic fraction
    gal['coldbar_frac'] = ( gal.HI_Mass + 
        gal.H2_Mass ) / ( gal.Total_Stellar_Mass + 
        gal.HI_Mass + gal.H2_Mass + 
        gal.Hot_Gas_Mass )
    
    gal['logcoldbar_frac'] = np.log10( gal.coldbar_frac )
    gal['logHeating'] = np.log10( 1e10 * gal.Heating )

    logTSMvir = np.log10( gal.Total_Stellar_Mass / 
        gal.Mvir )

    # Accretion rate
    gal['logBHL'] = np.log10( (0.1*gal.BHaccreted.to_numpy()*u.Msun/u.yr *const.c**2).to(u.erg/u.s).value )
    #gal.loc[( gal.BHaccreted == 0.0 ), 'logBHL'] = 1e-7

    # Black hole mass
    #gal.loc[( gal.Black_Hole_Mass == 0.0 ), 'Black_Hole_Mass'] = 1e-7

    gal.Black_Hole_Mass = 1e10 * gal.Black_Hole_Mass / h
    gal['logBHM'] = np.log10( gal.Black_Hole_Mass )
    gal['t_dyn'] = 2. * np.pi * 1e6 * (gal.Rvir / gal.Vvir )

    logTSM = np.log10( gal.Total_Stellar_Mass )
    logMvir = np.log10( gal.Mvir )
    logVmax = np.log10( gal.Vmax )
    logCGM = np.log10( gal.Cold_Gas_Mass )
    logHGM = np.log10( gal.Hot_Gas_Mass )
    gal['logTSM'] = logTSM
    gal['logMvir'] = logMvir                                          
    gal['logVmax'] = logVmax
    gal['logCGM'] = logCGM
    gal['logHGM'] = logHGM        
    gal['logTSMvir'] = logTSMvir

    if sample == "satellites":
        #subhalo information (only do this for the satellite population)
        gal = gal.loc[( gal.Subhalo_Mvir_at_Infall > 0 )]

        logsubTSMvir_infall = np.log10( gal.Total_Stellar_Mass /
            (1e10*gal.Subhalo_Mvir_at_Infall ) )
        logSubMvir_infall = np.log10( 1e10 * gal.Subhalo_Mvir_at_Infall )
        gal['logSubMvir_infall'] = logSubMvir_infall
        gal['logsubTSMvir_infall'] = logsubTSMvir_infall

        gal['q'] = (1e10 * gal.Central_Galaxy_Mvir) / gal.Mvir
        gal['qinfall'] = gal.Central_Galaxy_Mvir / gal.Subhalo_Mvir_at_Infall

    print('Time: %f s.' % (time.time() - start) ) 

    return gal