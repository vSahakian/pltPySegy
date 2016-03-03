def plot_shot_gather_byRange(obsnum,obsfile,shot_beg,shot_end,throw,ktpfile,linefile):
    #VJS 9/2014
    #plot_shot_gather.py
    #Plot shot gather for one OBS, including picks, by along-track range (not by shot)
    #Input:
    #       obsnum:         	    OBS number for shot gather
    #       obsfile:                String with full path for obs file
    #       shot_beg, shot_end:     First and last shot to plot
    #       throw:                  if ==1, throw out landshots. if == 0, keep land shots.
    #       ktpfile:                String with full path for ktp file
    #       linefile:               String with full path for shot file (with shot, x, y, z: along-track.)
    
    import sys
    import numpy as np
    import obspy as obs
    import matplotlib.pyplot as plt
    from struct import unpack
    sys.path.append('/Users/sahakian/Software/PY')
    import tred
    import read_ktp
    
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
    
    #pgtnew=[]
    ##ttmp=[]
    #for i in range(pgr.shape[0]):
    #    ttmp=tred.tred(pg[1][i],pgr[i],vred)
    #    pgtnew.append(ttmp)
    
    #Plot!
    plt.close('all')
    f=plt.figure(1)
    for i in range(len(obsn)):
        t=obsn[i].times()
        tnew=tred.tred(t,rng[i],vred)
        amp=obsn[i].data
        amp=4*(amp/max(abs(amp)))
        #Get positive part for shading
        ineg=amp<=0
        amp_pos=amp.copy()
        amp_pos[ineg]=0
        #Get negative part for shading
        ipos=amp>=0
        amp_neg=amp.copy()
        amp_neg[ipos]=0
        plt.fill_betweenx(tnew,i+shot_beg,amp_pos+i+shot_beg,facecolor='gray') #Shade positive
        #plt.fill_betweenx(tnew,amp_neg+i,i,facecolor='red') #Shade negative
        plt.plot(amp+i+shot_beg,tnew,color='0.3') #Plot unshaded wiggle
        
        #Plot travel-time picks - number of segments is in the 1 index of *.shape
        #pg: 
        ttmp=None
        for i in range(pg.shape[1]):
            ttmp=tred.tred(pg[1][i],pgr[i],vred)
            plt.plot(pg[0][i],ttmp,color='k',alpha=0.7,linewidth=2)
        #pn: - this has been finagled, need to fix read_ktp to get rid of empty array
        ttmp=None
        for i in range(pn.shape[1]):
            ttmp=tred.tred(pn[1][i],pnr[i],vred)
            plt.plot(pn[0][i],ttmp,color='m',alpha=0.7,linewidth=2)
        #pb:
        ttmp=None
        for i in range(pb.shape[1]):
            ttmp=tred.tred(pb[1][i],pbr[i],vred)
            plt.plot(pb[0][i],ttmp,color='b',alpha=0.7,linewidth=2)
            
        
        #plt.ylim([0,1.7])
        #plt.show()
        #plt.plot(amp+i,tnew,color='k')
    plt.ylim([0,1.5])
    #plt.xlim([min(obsn[0].data),len(obsn)+max(obsn[len(obsn)-1].data)])
    plt.xlim([shot_beg, shot_end])
    plt.xlabel('Shot Number')
    plt.ylabel('Time (s) - X/'+str(vred/1000)+'km/s')
    plt.savefig('obsn.png')
    plt.show()
    
    #f=plt.figure(1)
    #obsn[40].plot(show=True)
    ##plt.plot(trace70)
    #plt.savefig('obstest.png',format='png')
    
    