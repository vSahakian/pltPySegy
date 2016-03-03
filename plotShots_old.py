def plot_shot_gather_byRange(obsnum,obsfile,shot_beg,shot_end,tlims,throw,ktpfile,linefile,trscale,savepdf):
    #VJS 9/2014
    #plot_shot_gather.py
    #Plot shot gather for one OBS, including picks, by along-track range (not by shot)
    #Input:
    #       obsnum:         	    OBS number for shot gather
    #       obsfile:                String with full path for obs file
    #       shot_beg, shot_end:     First and last shot to plot
    #       tlims:                  Array: [t0 tend]
    #       throw:                  if ==1, throw out landshots. if == 0, keep land shots.
    #       ktpfile:                String with full path for ktp file
    #       linefile:               String with full path for shot file (with shot, x, y, z: along-track.)
    #       savepdf:                If ==1, save pdf form.  If ==0, do not save pdf.
    
    import sys
    import numpy as np
    import obspy as obs
    import matplotlib.pyplot as plt
    from matplotlib import rc
    from struct import unpack
    sys.path.append('/Users/vjsahakian/Software/PY')
    import tred
    import read_ktp
    
    #Define fonts for plot (default for axes, etc):
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':11})
    
    #Read OBS file
    obsn=obs.read(obsfile)
    #obsfile='/Users/sahakian/SaltonSea/OLD_segy/SSIP_OBS_segy/obs74/segy/obs74_Salton2011_7_L28ZFagc.segy'
    
    #Enter ktp (pick) file
    #ktpfile='/Users/sahakian/SaltonSea/Modeling/PkFiles/ttNewpicks/Site72_picks.ktp'
    
    #Throw out land shots option
    #throw=1
    
    #Enter beginning and last shot for plotting/correlation with picks
    #shot_beg=614
    #shot_end=747
    
    #Enter time limits of plot:
    t0=tlims[0]
    tend=tlims[1]
    
    #Read in shot info
    shotinfo=np.loadtxt(linefile)
    shotnums=shotinfo[:,0]
    
    #Find shots within range of shot_beg to shot_end
    shotind=np.where((shotnums>=shot_beg) & (shotnums<= shot_end))
    #Pull out x position of shots (shotx)
    shotx=shotinfo[shotind,1][0]
    
    shots=np.arange(shot_beg,shot_end+1,1)
    
    #Get source to receiver range
    rng=np.zeros(len(obsn))
    for k in range(len(obsn)):
        header=obsn[k].stats.segy.trace_header.unpacked_header
        rng[k]=unpack(">i",header[36:40])[0]
    shot_rang_match=[shots,rng]
    
    #Assign reduction velocity  
    vred = 3000  #m/s
    
    #Convert ktp file to lists/arrays with the various pick segments for each branch
    [pg, pn, pb]=read_ktp.read_ktp(ktpfile,throw,shot_end)
    pg,pn,pb=np.array(pg),np.array(pn),np.array(pb)
                    
    #Find ranges for each in order to get time with reduction velocity
    pgin,pnin,pbin=[],[],[]
    pgout,pnout,pbout=[],[],[]
    pgr,pnr,pbr=[],[],[]
    
    #For Pg
    for i in range(pg.shape[1]):
        pgin=[]
        for j in range(len(pg[0][i])):
            snum=pg[0][i][j]
            ind=np.where(shot_rang_match[0]==snum)[0]
            pgin.append(shot_rang_match[1][ind][0])
        pgout.append(pgin)
    pgr=np.array(pgout)
            
    #For Pn
    for i in range(pn.shape[1]):
        pnin=[]
        for j in range(len(pn[0][i])):
            snum=pn[0][i][j]
            ind=np.where(shot_rang_match[0]==snum)[0]
            pnin.append(shot_rang_match[1][ind][0])
        pnout.append(pnin)
    pnr=np.array(pnout)      
            
    #For Pb
    for i in range(pb.shape[1]):
        pbin=[]
        for j in range(len(pb[0][i])):
            snum=pb[0][i][j]
            ind=np.where(shot_rang_match[0]==snum)[0]
            pbin.append(shot_rang_match[1][ind][0])
        pbout.append(pbin)
    pbr=np.array(pbout)        
    
    #Plot!
    plt.close('all')
    f=plt.figure(num=1, figsize=(6.5,3.5))
    
    #Set figure size
    #f.set_size_inches(6,3)
    
    for i in range(len(obsn)):
        t=obsn[i].times()
        tnew=tred.tred(t,rng[i],vred)
        amp=obsn[i].data
        #Scale by 0.18 - arbitrary, based on visual.  This is to keep it normalized between the average shot range difference, ~0.08 for obs72
        amp=trscale*(amp/max(abs(amp)))
        #Get positive part for shading
        ineg=amp<=0
        amp_pos=amp.copy()
        amp_pos[ineg]=0
        #Get negative part for shading
        ipos=amp>=0
        amp_neg=amp.copy()
        amp_neg[ipos]=0
        tind=np.where((tnew>=t0) & (tnew<=tend))
        plt.fill_betweenx(tnew[tind],shotx[i],amp_pos[tind]+shotx[i],facecolor='gray') #Shade positive
        #plt.fill_betweenx(tnew,amp_neg+i,i,facecolor='red') #Shade negative
        plt.plot(amp[tind]+shotx[i],tnew[tind],color='0.3') #Plot unshaded wiggle
        
        #Plot travel-time picks - number of segments is in the 1 index of *.shape
        #pg: 
        ttmp=None
        for i in range(pg.shape[1]):
            ttmp=tred.tred(pg[1][i],pgr[i],vred)
            #Turn the shots into an array of integers to use as indices
            pg_shotind=np.array(pg[0][i]).astype(int)
            #subtract the first shot number from each called shot here so the index 
            #used for shotx starts at 0 if the first shot called is shot_beg, and after
            #that will be "offset" by shot_beg 
            plt.plot(shotx[pg_shotind-shot_beg],ttmp,color='k',alpha=0.7,linewidth=2)
        ##pn: - this has been finagled, need to fix read_ktp to get rid of empty array
        ttmp=None
        for i in range(pn.shape[1]):
            ttmp=tred.tred(pn[1][i],pnr[i],vred)
            pn_shotind=np.array(pn[0][i]).astype(int)
            plt.plot(shotx[pn_shotind-shot_beg],ttmp,color='m',alpha=0.7,linewidth=2)
        ##pb:
        ttmp=None
        for i in range(pb.shape[1]):
            ttmp=tred.tred(pb[1][i],pbr[i],vred)
            pb_shotind=np.array(pb[0][i]).astype(int)            
            plt.plot(shotx[pb_shotind-shot_beg],ttmp,color='b',alpha=0.7,linewidth=2)
            
        
        #plt.ylim([0,1.7])
        #plt.show()
        #plt.plot(amp+i,tnew,color='k')
    plt.ylim([0,1.5])
    #plt.xlim([min(obsn[0].data),len(obsn)+max(obsn[len(obsn)-1].data)])
    plt.xlim([shotx[0], shotx[-1]])
    plt.xlabel('Along-track Distance (km)',fontsize=12, fontname='Helvetica')
    plt.ylabel('Time (s) - X/'+str(vred/1000)+'km/s',fontsize=12, fontname='Helvetica')
    plt.title('OBS '+str(obsnum)+' Shot Gather',fontsize=12, fontname='Helvetica')
    
    #Save figure
    plt.savefig('obs'+str(obsnum)+'.png')
    
    if savepdf==1:
        plt.savefig('obs'+str(obsnum)+'.svg',format='svg')
        plt.show()
    else:
        plt.show()
    
    #f=plt.figure(1)
    #obsn[40].plot(show=True)
    ##plt.plot(trace70)
    #plt.savefig('obstest.png',format='png')
    
    