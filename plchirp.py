def plchirp(segyfile,fshot,lshot,stime,etime,cscale):
    #PLOT chirp data
    
    import sys
    import numpy as np
    import obspy as obs
    import matplotlib.pyplot as plt
    from matplotlib import rc
    from struct import unpack
    sys.path.append('/Users/sahakian/Software/PY')
    
    #Read segy file:
    seg=obs.read(segyfile)
    
    ##Set Parameters##
    t0=stime
    tend=etime
    