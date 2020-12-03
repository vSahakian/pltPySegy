def read_ktp(ktpfile,throw,last_shot):
    '''
    Convert ktp file used for pltsegy into digestible format for plotting.
    Usage:
        [pg,pn,pb]=read_ktp(ktpfile)
    Input:
        ktpfile:       String with full path of ktpfile
        throw:         if throw==1, throw out land shots.  if throw==0, keep land shots.
        last_shot:     last shot to compute
    Output:
        pg:            list of first-arrival pick segments
        pn:            list of second-arrival pick segments
        pb:            list of basement or reflector pick segments 
    Note on format:
        p*[0] contains lists of shots for each segment
        p*[1] contains lists of pick times for each segment   
    '''
    
    import numpy as np
    
    pg_t, pn_t, pb_t=np.array([]),np.array([]),np.array([])
    pg_s, pn_s, pb_s=np.array([]),np.array([]),np.array([])
    
    f=open(ktpfile,'r')
    ktp=f.readlines()
    k=0
    print k
    while k<len(ktp)-1:
        if ktp[k]=='Pg\n':
            k=k+6
            while ktp[k]!='comment\n':
                str=ktp[k].split()
                shot=float(str[0])
                time=float(str[1])
                pg_t=np.append(pg_t,time)
                pg_s=np.append(pg_s,shot)
                k=k+1
        if ktp[k]=='Pn\n':
            k=k+6
            while ktp[k]!='comment\n':
                str=ktp[k].split()
                shot=float(str[0])
                time=float(str[1])
                pn_t=np.append(pn_t,time)
                pn_s=np.append(pn_s,shot)
                k=k+1
        if ktp[k]=='Pb\n':
            #Need to do k+5 instead of k+6 because the next loop needs the counter at the beginning not the end, so the while loop will end on the last k
            k=k+5
            while ktp[k]!='comment\n' and k<len(ktp)-1:  #<len(ktp)-1 because the counter updates right after this - so going into it, it must be less than one less than end of file
                k=k+1
                str=ktp[k].split()
                shot=float(str[0])
                time=float(str[1])
                pb_t=np.append(pb_t,time)
                pb_s=np.append(pb_s,shot)
                
        else:
            k=k+1
            
    
    #Find where there are NaN's        
    iNan_pgs=np.where(np.isnan(pg_s))[0]
    iNan_pns=np.where(np.isnan(pn_s))[0]
    iNan_pbs=np.where(np.isnan(pb_s))[0]
    
    #Initialize new lists, we'll append each "line" (before nan) to this
    pg_times, pn_times, pb_times=[],[],[]
    pg_shots, pn_shots, pb_shots=[],[],[]
    pgs_tmp,pns_tmp,pbs_tmp=[],[],[]
    pgt_tmp,pnt_tmp,pbt_tmp=[],[],[]
    
    #If there are non nan's, what do you do (just take teh value of the original shots/times array):
    if len(iNan_pgs)==0:
        pg_shots=pg_s
        pg_times=pg_t
    #If there are nan's, we want to make separate lines, so append.  
    if len(iNan_pgs)>0:
         for i in range(len(iNan_pgs)):
             #First, if it's the first nan, take everything before that nan until that nan
             if i==0:
                 pgs_tmp=pg_s[0:iNan_pgs[0]]
                 pgt_tmp=pg_t[0:iNan_pgs[0]]
                 pg_shots.append(pgs_tmp)
                 pg_times.append(pgt_tmp)
             #Then go from the last nan to the current one
             else:
                 pgs_tmp=pg_s[iNan_pgs[i-1]+1:iNan_pgs[i]]
                 pgt_tmp=pg_t[iNan_pgs[i-1]+1:iNan_pgs[i]]
                 pg_shots.append(pgs_tmp)
                 pg_times.append(pgt_tmp)
    #Then turn into arrays:  
    pg_shots=np.array(pg_shots)
    pg_times=np.array(pg_times)
     
#Now for the Pn picks 
    if len(iNan_pns)==0:
        pn_shots=pn_s
        pn_times=pn_t
    if len(iNan_pns)>0:
         for i in range(len(iNan_pns)):
             if i==0:
                 pns_tmp=pn_s[0:iNan_pns[0]]
                 pnt_tmp=pn_t[0:iNan_pns[0]]
                 pn_shots.append(pns_tmp)
                 pn_times.append(pnt_tmp)
             else:
                 pns_tmp=pn_s[iNan_pns[i-1]+1:iNan_pns[i]]
                 pnt_tmp=pn_t[iNan_pns[i-1]+1:iNan_pns[i]]
                 pn_shots.append(pns_tmp)
                 pn_times.append(pnt_tmp)
    pn_shots=np.array(pn_shots)
    pn_times=np.array(pn_times)
                     
#And the Pb picks (reflections)
    if len(iNan_pbs)==0:
        pb_shots=pb_s
        pb_times=pb_t
    if len(iNan_pbs)>0:
         for i in range(len(iNan_pbs)):
             if i==0:
                 pbs_tmp=pb_s[0:iNan_pbs[0]]
                 pbt_tmp=pb_t[0:iNan_pbs[0]]
                 pb_shots.append(pbs_tmp)
                 pb_times.append(pbt_tmp)
             else:
                 pbs_tmp=pb_s[iNan_pbs[i-1]+1:iNan_pbs[i]]
                 pbt_tmp=pb_t[iNan_pbs[i-1]+1:iNan_pbs[i]]
                 pb_shots.append(pbs_tmp)
                 pb_times.append(pbt_tmp)
    pb_shots=np.array(pb_shots)
    pb_times=np.array(pb_times)
 
  #Give option to throw out land shots 
    pstemp_in,pstemp_out=[],[]  
    pttemp_in,pttemp_out=[],[] 
    
    #Throw out for Pg
    if throw==1:
        for i in range(pg_shots.shape[0]):
            pstemp_in,pttemp_in=[],[]
            for j in range(len(pg_shots[i])):
                if pg_shots[i][j]<=last_shot:                  
                    pstemp_in.append(pg_shots[i][j])
                    pttemp_in.append(pg_times[i][j])
            pstemp_out.append(pstemp_in)
            pttemp_out.append(pttemp_in)
    pgs_noland=np.array(pstemp_out)
    pgt_noland=np.array(pttemp_out)
        
    #Do for Pn...        
    pstemp_in,pstemp_out=[],[]  
    pttemp_in,pttemp_out=[],[] 
    
    if throw==1:
        for i in range(pn_shots.shape[0]):
            pstemp_in,pttemp_in=[],[]
            for j in range(len(pn_shots[i])):
                if pn_shots[i][j]<=last_shot:                  
                    pstemp_in.append(pn_shots[i][j])
                    pttemp_in.append(pn_times[i][j])
            pstemp_out.append(pstemp_in)
            pttemp_out.append(pttemp_in)
    pns_noland=np.array(pstemp_out)
    pnt_noland=np.array(pttemp_out)            
    
    #And for Pb
    pstemp_in,pstemp_out=[],[]  
    pttemp_in,pttemp_out=[],[]  
    
    if throw==1:
        for i in range(pb_shots.shape[0]):
            pstemp_in,pttemp_in=[],[]
            for j in range(len(pb_shots[i])):
                if pb_shots[i][j]<=last_shot:                  
                    pstemp_in.append(pb_shots[i][j])
                    pttemp_in.append(pb_times[i][j])
            pstemp_out.append(pstemp_in)
            pttemp_out.append(pttemp_in)
    pbs_noland=np.array(pstemp_out)
    pbt_noland=np.array(pttemp_out)     
        
        
    if throw==0:
        pg,pn,pb=[pg_shots,pg_times],[pn_shots, pn_times],[pb_shots, pb_times]
        pg,pn,pb=np.array(pg),np.array(pn),np.array(pb)
    if throw==1:
        pg,pn,pb=[pgs_noland,pgt_noland],[pns_noland,pnt_noland],[pbs_noland,pbt_noland]
        pg,pn,pb=np.array(pg),np.array(pn),np.array(pb)
    return pg, pn, pb
         
    
