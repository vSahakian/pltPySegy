###############
## FUNCTIONS ##
###############


############
def uniform_gain(stream_obj,gain):
    '''
    Apply a uniform gain to the data, amp(t)*gain
    Input:
        stream_obj:             Stream object with traces of segy data
        gain:                   Float with gain to apply to data
    Output:
        gain_stream:            Stream object with gain applied
    '''
    from numpy import require,float32
    
    ## Make copy of stream object:
    gain_stream = stream_obj.copy()
    
    ## For each trace, apply the exponential function:
    for i_trace in range(len(stream_obj)):
        gain_stream[i_trace].data = require(stream_obj[i_trace].data * gain,dtype=float32)
    
    ## Return the new stream:
    return gain_stream


############
def sqrt_gain(stream_obj,gain):
    '''
    Apply a sqrt gain function to the data, amp(t)*sqrt(gain*t)
    Input:
        stream_obj:             Stream object with traces of segy data
        gain:                   Float with gain to apply to data
    Output:
        gain_stream:            Stream object with gain applied
    '''
    from numpy import require,float32,sqrt
    
    ## Make copy of stream object:
    gain_stream = stream_obj.copy()
    
    ## For each trace, apply the exponential function:
    for i_trace in range(len(stream_obj)):
        i_times = stream_obj[i_trace].times()
        gain_stream[i_trace].data = require(stream_obj[i_trace].data * sqrt(gain*i_times),dtype=float32)
    
    ## Return the new stream:
    return gain_stream



############
def exp_gain(stream_obj,exp_power):
    '''
    Apply an exponential gain function to the data, amp(t)*exp(exp_power*t)
    Input:
        stream_obj:             Stream object with traces of segy data
        exp_power:              Scalar with exponential power to apply to data
    Output:
        gain_stream:            Stream object with gain applied
    '''
    from numpy import exp,require,float32
    
    ## Make copy of stream object:
    gain_stream = stream_obj.copy()
    
    ## For each trace, apply the exponential function:
    for i_trace in range(len(stream_obj)):
        i_times = stream_obj[i_trace].times()
        gain_stream[i_trace].data = require(stream_obj[i_trace].data * exp(exp_power*i_times),dtype=float32)
    
    ## Return the new stream:
    return gain_stream


###########
def gain(stream,option1,parameters,option2):
    '''
    GAIN: Gain a group of traces.
    
      gain(d,dt,option1,parameters,option2);
    
      Input   
          stream:    Obspy stream object to gain
          option1:   = 'time' parameters = [a,b],  gain = t.^a . * exp(-bt)
                     = 'agc' parameters = [agc_gate], length of the agc gate in secs
          option2:   = 0  No normalization
                     = 1  Normalize each trace by amplitude
                     = 2  Normalize each trace by rms value
    
      OUT  dout(nt,nx): traces after application of gain function
    '''
    
    import numpy as np
    import math
    from scipy.signal.windows import triang
    from scipy.signal import convolve2d as conv2
    
    ## Get number of samples and traces:
    nx = len(stream)
    nt = len(stream[0].times())
    
    dt = stream.stats.binary_file_header.sample_interval_in_microseconds * 1e-6

    gained_stream = stream.copy()
    
    if option1 == 'time':
        a = parameters[0]
        b = parameters[1]
        t = [x*dt for x in range(nt)]
        tgain = [(x**a)*math.exp(x*b) for x in t]

        for k in range(nx):
            gained_stream[k].data = stream[k].data*tgain

    elif option1 == 'agc':
        L = parameters/dt+1
        L = np.floor(L/2)
        h = triang(2*L+1)
        shaped_h  = h.reshape(len(h),1)

        for k in range(nx):
            aux = stream[k].data
            e = aux**2
            shaped_e = e.reshape(len(e),1)
            
            rms = np.sqrt(conv2(shaped_e,shaped_h,"same"))
            epsi = 1e-10*max(rms)
            op = rms/(rms**2+epsi)
            op = op.reshape(len(op),)

            gained_stream[k].data = stream[k].data*op

    #Normalize by amplitude 
    if option2==1:
        for k in range(nx):
            aux = stream[k].data
            amax = max(abs(aux))
            gained_stream[k].data = gained_stream[k].data/amax 

    #Normalize by rms 
    if option2==2:
        for k in range(nx):
            aux = stream[k].data
            amax = np.sqrt(sum(aux**2)/nt)
            gained_stream[k].data = gained_stream[k].data/amax


    return gained_stream
    

############
def get_alongtrackdistance(x,y,x0=0,y0=0,distance_format='fromfirstpoint'):
    '''
    compute alongtrack distance for UTM measurments
    Input:
        x:              Array with UTM x values
        y:              Array with UTM y values
        x0:             Float
    '''
    from numpy import sqrt,diff,append,cumsum
    
    ## If it's distance from the first point:
    if distance_format == 'fromfirstpoint':
        ## Get along track distance - subtract the first value of x and y from the array to get distance along line:
        alongtrack_distance_utm = sqrt((x - x0)**2 + (y - y0)**2)
        ## Divide by 1000 for km
        alongtrack_distance = alongtrack_distance_utm/1000
    elif distance_format == 'string':
        ## Get the differences betwen them, in km
        dx = diff(x/1000)
        dy = diff(y/1000)
        dl = sqrt(dx**2 + dy**2)
        alongtrack_distance = append(0,cumsum(dl))
        
    return alongtrack_distance




############
def extract_data(segy_st,plotheight,kmperinch,distance_format='fromfirstpoint'):
    '''
    Plot all shots, don't filter by any shot
    Input:
        segy_st:                Stream object with segy data to plot
        plotheight:             Plot height in inches
        kmperinch:              Number of kilometers per inch to plot
        distance_format:        Format to use for computing alongtrack distance.
                                    Default: 'fromfirstpoint' - distance from the first point.
                                    Also could use 'string' - cumulative distance
    Output:
        alongtrack_distance:    1D array with the alongtrack distance of the shots on this line, in km
        times:                  1D array with the times on this line in seconds
        DISTANCE:               2d array from meshgrid with the distances, shape len(distance),len(times)
        TIMES:                  2d array from meshgrid with the times, shape len(distance),len(times)
        AMP:                    2d array from meshgrid with the seismic amplitudes, shape len(distance),len(times)
        plotwidth:              Float with the width in inches of the line
        trace_x:                Array with the UTM x values of the traces
        trace_y:                Array with the UTM y values of the traces
    '''
    
    from numpy import array,append,meshgrid,zeros_like
    
    ## Get times for plotting from first trace/shot:
    times = segy_st[0].times()
    
    print('sorting and getting along track distances')
    ## Get source x and y:
    trace_x = array([])
    trace_y = array([])
    for i_trace in range(len(segy_st)):
        i_trace_x = segy_st[i_trace].stats.segy.trace_header['source_coordinate_x']
        i_trace_y = segy_st[i_trace].stats.segy.trace_header['source_coordinate_y']
        i_trace_coordscalar = segy_st[i_trace].stats.segy.trace_header['scalar_to_be_applied_to_all_coordinates']
        
        ## Append to arrays, but first divide by the coordinate scalar that puts it in meters:
        trace_x = append(trace_x,i_trace_x/i_trace_coordscalar)
        trace_y = append(trace_y,i_trace_y/i_trace_coordscalar)
        
    ## Get distance:
    alongtrack_distance = get_alongtrackdistance(trace_x,trace_y,x0=trace_x[0],y0=trace_y[0],distance_format=distance_format)
    
#    ## If it's distance from the first point:
#    if distance_format == 'fromfirstpoint':
#        ## Get along track distance - subtract the first value of x and y from the array to get distance along line:
#        alongtrack_distance_utm = sqrt((trace_x - trace_x[0])**2 + (trace_y - trace_y[0])**2)
#        ## Divide by 1000 for km
#        alongtrack_distance = alongtrack_distance_utm/1000
#    elif distance_format == 'string':
#        ## Get the differences betwen them, in km
#        dx = diff(trace_x/1000)
#        dy = diff(trace_y/1000)
#        dl = sqrt(dx**2 + dy**2)
#        alongtrack_distance = append(0,cumsum(dl))
        
    ## GEt the plot dimensions, given the km per inch, and the total number of km on the line:
    plotwidth = alongtrack_distance[-1] / kmperinch
    
    print('getting meshgrid')
    ## Get a meshgrid with distance as the x to plot and time as the y:
    DISTANCE,TIMES = meshgrid(alongtrack_distance,times)
    ## Initiate the amplitude array:
    AMP = zeros_like(TIMES)
    
    print('getting amplitudes')
    ## Grab amplitudes out:
    for i_trace in range(len(segy_st)):
        AMP[:,i_trace] = segy_st[i_trace].data
    
    return alongtrack_distance,times,DISTANCE,TIMES,AMP,plotwidth,trace_x,trace_y


#############
#def 


############
def interpolate_data(AMP,distance,time,distance_newsamplenum,time_newsamplenum):
    '''
    Given amplitude, distance, and time arryas, interpolate with a new sampling "decimation" interval
    Input:
        AMP:                    Array with shape (len(distance),len(times)) for seismic amplitudes
        distance:               1D array with distances that come from segy stream
        time:                   1D array with times that come from segy stream
        distance_newsamplenum:  Float with how frequently to interpolate new distance on this array
        time_newsamplenum:      Float with how frequently to interpolate new time on this array
    Output:
        DISTANCEinterp:         Array with new distance, meshgrid shapes from distance_newsamplenum and time_newsamplenum     
        TIMEinterp:             Array with new times, meshgrid shapes from distance_newsamplenum and time_newsamplenum  
        amp_interp:             Array with new amplitudes, meshgrid shapes from distance_newsamplenum and time_newsamplenum  
        distance_new:           1D array with new distance values
        time_new:               1D array with new time values
    
    '''
    from scipy.interpolate import RectBivariateSpline
    from numpy import unique,mean,diff,meshgrid,arange

    ## Get unique values of distance (since three shots a second) - this is needed
    ##   for rect bivariate spline, strictly increasing in x and y:
    unique_distance, udistind= unique(distance,return_index=True)
        

    ## Get interpolant on old data - but the dimensions of the Z array
    ##   must be (len(x),len(y)). The dimensions of AMP to match meshgrid
    ##   are (len(y),len(x)). So find the *columns* where AMP matches the first
    ##   of the unique value in distance (udistind), then transpose it so its shape
    ##   matches (len(unique_distance),len(times)):
    splinef = RectBivariateSpline(unique_distance,time,AMP[:,udistind].T)
    
    ## Get new distance and time values, with spacing multiplied by new sample number.
    
    ## GEt average spacing of the unique distance values and of all the time values:
    dist_spacing = diff(unique_distance)
    mean_dist_spacing = mean(dist_spacing)
    mean_time_spacing = mean(diff(time))
    
    ## Get new sampled array with the number of samples 
    distance_new = arange(distance[0],distance[-1]+(mean_dist_spacing/distance_newsamplenum),mean_dist_spacing/distance_newsamplenum)
    time_new = arange(time[0],time[-1]+(mean_time_spacing/time_newsamplenum),mean_time_spacing/time_newsamplenum)

    print('interpolating')
    amp_interp = splinef(distance_new,time_new)
    DISTANCEinterp, TIMEinterp = meshgrid(distance_new, time_new)    

    ## return - return the transpose of the interpolated amplitudes to match
    ##   meshgrid of distance and time
    return DISTANCEinterp,TIMEinterp,amp_interp.T,distance_new,time_new

############
def subsample_grids(alongtrack_distance,times,AMP,DISTANCE,TIMES,starttime=None,endtime=None,startdistance=None,enddistance=None):
    '''
    Subsample 
    Input: 
        alongtrack_distance:    1D array with the alongtrack distance, presampling
        times:                  1D array with the times, presampling
        AMP:                    2D array with the seismic amplitudes, presampling, shape len(distance),len(times)
        DISTANCE:               2D array with the meshgrid distances, presampling, shape len(distance),len(times)
        TIMES:                  2D array with the meshgrid times, presampling, shape len(distance),len(times)
        starttime:              Float with the start time for sampling. Default: None
        endtime:                Float with the end time for sampling. Default: None
        startdistance:          Float with the start distance for sampling. Default: None
        enddistance:            Float with the end distance for sampling. Default: None
    Output:
        DISTANCE_sampled:       2D array with the meshgrid distances, sampled, shape len(sampdistance),len(samptimes)
        TIMES_sampled:          2D array with the meshgrid times, sampled, shape len(sampdistance),len(samptimes)
        AMP_sampled:            2D array with the meshgrid amplitudes, sampled, shape len(sampdistance),len(samptimes)
    '''
    from numpy import where
    
    print('subsampling')
    
    ## If it's only time being subsampled:
    if (startdistance == None) & (enddistance == None) & (starttime != None) & (endtime != None):
        ## get plot index for times and amps:
        time_plot_ind = where((times >= starttime) & (times <=endtime))[0]
        TIMES_sampled = TIMES[time_plot_ind,:]
        DISTANCE_sampled = DISTANCE[time_plot_ind,:]
        AMP_sampled = AMP[time_plot_ind,:]
    
    elif (startdistance != None) & (enddistance != None) & (starttime == None) & (endtime == None):
        ## get plot index for distance and amps:
        dist_plot_ind = where((alongtrack_distance >= startdistance) & (alongtrack_distance <= enddistance))[0]
        TIMES_sampled = TIMES[:,dist_plot_ind]
        DISTANCE_sampled = DISTANCE[:,dist_plot_ind]
        AMP_sampled = AMP[:,dist_plot_ind]
        
    elif (startdistance != None) & (enddistance != None) & (starttime != None) & (endtime != None):
        ## get plot index for distance, times and amps:
        time_plot_ind = where((times >= starttime) & (times <=endtime))[0]
        dist_plot_ind = where((alongtrack_distance >= startdistance) & (alongtrack_distance <= enddistance))[0]
        TIMES_sampled = TIMES[time_plot_ind,dist_plot_ind]
        DISTANCE_sampled = DISTANCE[time_plot_ind,dist_plot_ind]
        AMP_sampled = AMP[time_plot_ind,dist_plot_ind]
    
#    ## INterpolate the range you want to plot
    ## get new distance and time to use.
    return DISTANCE_sampled,TIMES_sampled,AMP_sampled
    
    

    
############  
def plot_chirp(DISTANCEp,TIMESp,AMPp,plotcmap,plotheight,plotwidth,pcmin,pcmax,line4title,watervel,trace_x,trace_y,horizons=None,horizoncolors=None,faults=None,faultstyles=None,depthshift=0):
    '''
    Plot given distance/time/amplitude grids from chirp data
    Input:
        DISTANCEp:          Array with distances in grid form from meshgrid (shape: len(time) x len(alongtrack_distance))
        TIMESp:             Array with times in grid form from meshgrid (shape: len(time) x len(alongtrack_distance))
        AMPp:               Array with amplitudes from chirp data, interpolated/gridded/subsampled, any combination thereof or none. (shape: len(time) x len(alongtrack_distance))
        plotcmap:           String with colormap to use
        plotheight:         Float with plot height in inches
        plotwidth:          Float with plot width in inches
        pcmin:              Minimum for plot color
        pcmax:              Maximum for plot color
        line4title:         String with line number for title
        watervel:           Float with the water velocity to use for righ taxis in m/s
        trace_x:            Array with the UTM x values of the traces
        trace_y:            Array with the UTM y values of the traces
        horizons:           Dictionary with horizon info for the line. If None (default), no horizon plotted.
        horizoncolors:      Pandas dataframe with columns: 'horizon', 'color'. Color in hex.
        faults:             Dictionary with fault info for the line. If None (default), no faults plotted.
        faultstyles:        Pandas dataframe with columns: 'fault', 'color', 'linestyle'. Color in hex.
        depthshift:         Float with the value in meters to shift all converted depths (i.e., assume 0 TWTT is @ depthshift meters)
    Output:
        allshot_figure
    '''
    print('plotting')
    
    import matplotlib.pyplot as plt
    from numpy import str,array
    import pandas as pd
    
    ## Function to conevrt TWTT to meters:
    depth_meters = lambda depth_twtt: ((depth_twtt/2)*watervel) + depthshift
    depth_twtt = lambda depth_meters: ((depth_meters - depthshift)/watervel)*2
    
    ## Make figure
    ## the plot height and width will actulaly depend on if there is a legend,
    ##   or not. 
    if ((horizons == None) and (faults == None)):
        ## Make figure
        print('no legend')
        allshots_figure, twtt_ax = plt.subplots(figsize=(plotwidth/0.8,plotheight))
    else:
        ## Then, want the height to be the same, but need the variable "plotwidth"
        ##   to be 30% larger. i.e., plotwidth = 0.7 of the figure width,
        ##   so figure width = plotwidth (specified) / 0.7.
        allshots_figure, twtt_ax = plt.subplots(figsize=(plotwidth/0.65,plotheight))
        print('legend')
    
    
    #allshots_figure, twtt_ax = plt.subplots(figsize=(plotwidth,plotheight))
    ## Set labels on a second y axis on right
    #rightax = allshots_axis.twinx()
    m_ax = twtt_ax.secondary_yaxis("right", functions=(depth_meters,depth_twtt))

    twtt_ax.pcolormesh(DISTANCEp,TIMESp,AMPp,cmap=plotcmap,vmin=pcmin,vmax=pcmax,shading='gouraud')
    
    
    ## If there were horizons specified:
    if horizons != None:
        ## then for each horizon in the dictionary, plot all segments, separately:
        for i_horizon_ind in range(len(horizons['horizons'])):
            ## loop through segments:
            i_horizon = horizons['horizons'][i_horizon_ind]
            i_horizon_color = horizoncolors.loc[horizoncolors.horizon == i_horizon].color.values[0]
            for j_segment_ind in range(len(horizons[i_horizon]['X'])):
                ij_x = array(horizons[i_horizon]['X'][j_segment_ind])
                ij_y = array(horizons[i_horizon]['Y'][j_segment_ind])
                ij_z = array(horizons[i_horizon]['Z'][j_segment_ind])/1000 ## convert from ms to s
                
                ## Get alongtrack distance:
                ij_dist = get_alongtrackdistance(ij_x,ij_y,x0=trace_x[0],y0=trace_y[0])
                
                ## Plot - if it's the first, add a label:
                if j_segment_ind == 0:
                    twtt_ax.plot(ij_dist,ij_z,linewidth=0.8,color=i_horizon_color,alpha=0.75,label=i_horizon)
                else:
                    twtt_ax.plot(ij_dist,ij_z,linewidth=0.8,color=i_horizon_color,alpha=0.75)
                    
    ## If there were faults given:
    if faults != None:
        ## make a list for hte labels
        faultplot_counter = []
        ## then for each fault segment in the dictionary, plot:
        for i_fault_ind in range(len(faults['X'])):
            ## plot each segment.
            i_fault = faults['faultnames'][i_fault_ind]
            i_fault_color = faultstyles.loc[faultstyles.fault == i_fault].color.values[0]
            i_fault_linestyle = faultstyles.loc[faultstyles.fault == i_fault].linestyle.values[0]
            
            ## Get X, Y, Z
            ij_x = array(faults['X'][i_fault_ind])
            ij_y = array(faults['Y'][i_fault_ind])
            ij_z = array(faults['Z'][i_fault_ind])/1000 # convert from ms to s
            
            ## Get alongtrack distance:
            ij_dist = get_alongtrackdistance(ij_x,ij_y,x0=trace_x[0],y0=trace_y[0])
            
            ## Plot.
            ## If this fault is not in the list yet (hasn't been plotted on this line),
            ##   then add it and label it:
            if i_fault not in faultplot_counter:
                faultplot_counter.append(i_fault)
                twtt_ax.plot(ij_dist,ij_z,linewidth=1,color=i_fault_color,linestyle=i_fault_linestyle,alpha=0.75,label=i_fault)
            else:
                twtt_ax.plot(ij_dist,ij_z,linewidth=1,color=i_fault_color,linestyle=i_fault_linestyle,alpha=0.75)
                
    if ((faults != None) or (horizons !=None)):
        ## First, adjust subplots to give room for legend:
        plt.subplots_adjust(left=0.1,right=0.75)
        twtt_ax.legend(loc='center left',bbox_to_anchor=(1.10,0.5),ncol=1)
        print('added legend')
    elif ((faults == None) and (horizons == None)):
        ##don'ta djust right hand side, just left...
        plt.subplots_adjust(left=0.1,right=0.9)
        
    ## Set axis limits for time, in case faults are deeper:
    twtt_ax.set_ylim([TIMESp[0,0],TIMESp[-1,0]])     
    
    ## invert axis
    twtt_ax.invert_yaxis()

    ## Add label to this axis:
    m_ax.set_ylabel('Depth (m) \n v = ' + str(watervel) + ' m/s')
    
    plt.xlabel('Alongtrack Distance (km)')
    twtt_ax.set_ylabel('TWTT (s)')
    plt.title('Line ' + line4title)
            
    
    return allshots_figure
