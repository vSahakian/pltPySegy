#Script to plot all OBS
#VJS 9/2014
#

import numpy as np
import plotShots

#Do for all OBS in here (line 7)
obs=np.array([70, 71, 72, 73, 74, 28, 75, 76, 77])
#obs=np.array([70, 72])


#Shot XYZ for line 7
linefile='/Users/dmelgar/Downloads/canopystuff/salton_line7_landshot.ixyz'

#Set other params:
throw=1          #throw out land shots
shot_beg=70000   #First shot to use
shot_end=70090   #Last shot to use
savepdf=1        #Save pdf versions of plots?
trscale=0.09     #Number by which to scale/normalize amplitudes
tlims=np.array([0,5]) #Time limits

#Now make plots:
for i in range(len(obs)):
    obsnum=obs[i]
    print obsnum
    obsfile='/Users/dmelgar/Downloads/canopystuff/obs'+np.str(obsnum)+'_Salton2011_land_L28ZFagc.segy'
    ktpfile='/Users/dmelgar/Downloads/canopystuff/Site'+np.str(obsnum)+'_picks.ktp'
    plotShots.plot_shot_gather_byRange(obsnum,obsfile,shot_beg,shot_end,tlims,throw,ktpfile,linefile,trscale,savepdf)
    
    

