def getnav(segyfile,navfile):
    '''
    VJS 1/2016
    Extract nav from a segy
    Input:
       segyfile:   String with segy file to extract nav
       navfile:    String with output nav file
    Output:
       To navfile
       Shot     Lon     Lat
    
    '''
    
    import numpy as np
    import obspy as obs
    from struct import unpack
    
    #Import segy data:
    sdat=obs.read(segyfile)
    
    #Extract header:
    #Set to empty arrays:
    lon=np.zeros((len(sdat),1))
    lat=np.zeros((len(sdat),1))
    shot=np.zeros((len(sdat),1))
    
    #run through each trace:
    for i in range(len(sdat)):
        header=sdat[i].stats.segy.trace_header.unpacked_header
        #Unpack the source lat/lon, format is arcseconds:
        lonh=unpack(">i",header[72:76])[0]
        lath=unpack(">i",header[76:80])[0]
        shot[i]=unpack(">i",header[8:12])[0]
        #Convert to decimal degrees:
        lon[i]=lonh/(60.0*60.0)
        lat[i]=lath/(60.0*60.0)
        
    #Print to file:
    out=np.c_[shot,lon,lat]
    np.savetxt(navfile,out,fmt='%6i\t%12.8f\t%10.8f')
    
    
def nav_ll2utm(navfile_ll,navfile_utm):
    '''
    Convert a navfile made above to a utm file for kingdom
    VJS 1/2016
    Input:  
        navfile_ll:     String of path to navfile, format: shot lon lat
        navfile_utm:    String of path to utm navfile, format: shot X Y
    Output:
        To navfile_utm:
            Shot/RP  Lon  Lat  
    '''
    
    import numpy as np
    from pyproj import Proj
    
    #Read in ll navfile:
    lldat=np.genfromtxt(navfile_ll)
    
    #Sort:
    shot=lldat[:,0]
    lon=lldat[:,1]
    lat=lldat[:,2]
    
    #Make projection:
    p=Proj(proj='utm',zone='11S',ellps='WGS84')
    
    #Project:
    UTMx,UTMy=p(lon,lat)

    #Print out and save:
    out=np.c_[shot,UTMx,UTMy]
    np.savetxt(navfile_utm,out,fmt='%6i\t%12.5f\t%12.5f')
    