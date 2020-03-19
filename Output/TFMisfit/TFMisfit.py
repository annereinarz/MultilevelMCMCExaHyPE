#!/home/anne/anaconda3/bin/python3
import sys,getopt
import pandas as pd
import obspy as obs
import numpy as np
import timeFrequencyMisfit_2 as tf_misfit

def processExaData(num,var,var_names):
    processed = pd.DataFrame()
    processed['time'] = pd.to_numeric(num[num.columns.values[1]])
    for i in range(0,len(var)):
        processed[var_names[i]]    = pd.to_numeric(num[num.columns.values[2+var[i]]])
    processed.sort_values(by=["time"])
    return processed

# Read seismograms into obspy format
# station = station name + id
# synth = new var name for the signal in obspy format
def CreateObspyTraceFromSeismogram(station, synth,variables):
    #load the synthetics into an obspy dataformat
    st_syn= obs.read()
    st_syn.clear()
    xyz='xyz'
    for var in variables:
        tr = obs.Trace()
        tr.stats.station = 'receiver.'+var
        tr.data = synth[var].as_matrix()
        
        tr.stats.delta = synth["time"][1]-synth["time"][0]
        tr.stats.starttime = obs.UTCDateTime(2015, 5, 12, 7, 5, 19.00)  #put whatever date you want, eg. your Date of Birth
        st_syn.append(tr)
    return st_syn

def interpolateSignal(st,t1,t2):
    x=t2
    xvals=t1
    _sample=np.zeros((1,len(t1)))
    for i in range(0,1):
        y=st[i]
        yinterp=np.interp(xvals,x,y)
        _sample[i]=yinterp
    return _sample

def main(argv):
    arguments_short="hr:s:o:l:n:e:p"
    arguments_long =["help","reference=","simulation=","output=","low_frequency","high_frequency","end_time","show"]
    arguments_span =[ ("--"+arguments_long[i].replace("=",""), "-"+list(arguments_short.replace(":",""))[i] ) for i in range(0,len(arguments_long)) ]
    output="output"
    arg_show=False

    end_time=False

    f1 = 0.001
    f2 = 300
    help_string="./TFMisfit.py -r <reference> -s <simulation> -o <output> -l <low_frequency> -n <high_frequency> -e <end_time> -p"

    try:
        opts, args = getopt.getopt(argv,arguments_short,arguments_long)
    except getopt.GetoptError as error:
        print(error)
        print(help_string)
        sys.exit(2);
    for opt, arg in opts:
        if opt in arguments_span[0] : #help
            print("Help:")
            print(help_string)
            sys.exit()
        elif opt in arguments_span[1] :#reference
            ref_seissmogram=arg
        elif opt in arguments_span[2] :#simulation
            sim_seissmogram=arg
        elif opt in arguments_span[3] :#output
            output=arg
        elif opt in arguments_span[4] :#low_frequency
            f1=float(arg)
        elif opt in arguments_span[5] :#high_frequency
            f2=float(arg)
        elif opt in arguments_span[6] :#end_time
            end_time=float(arg)
        elif opt in arguments_span[7] :#show
            arg_show=True


    

    #TODO make this an argument
    var_indeces=[1]
    var_names=["h"]


    df_sim = pd.read_csv(sim_seissmogram)
    df_sim = processExaData(df_sim,var_indeces,var_names)
    
    df_ref = pd.read_csv(ref_seissmogram)
    df_ref = processExaData(df_ref,var_indeces,var_names)

    #End time defined by reference data
    if(end_time):
        max_time=end_time
    else:
        max_time=df_ref['time'][df_ref['time'].count()-1]
        
    min_time=df_ref['time'][ 0 ]
    df_sim = df_sim[(df_sim['time'] <= max_time) & (df_sim['time'] >= min_time)]
    df_ref = df_ref[(df_ref['time'] <= max_time) & (df_ref['time'] >= min_time)]

    ref_data= CreateObspyTraceFromSeismogram('reference',df_ref,var_names)

    #apply bandpass filter
    ref_flt = ref_data.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)

    sim_data= CreateObspyTraceFromSeismogram('simerence',df_sim,var_names)
    #apply bandpass filter
    sim_flt = sim_data.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)

    #interpolate ref data
    ref_flt_int    = interpolateSignal(ref_flt,df_sim["time"],df_ref["time"])

    dt_sim=df_sim["time"][1]-df_sim["time"][0]
    tf_misfit.plot_tf_misfits(sim_flt,ref_flt_int,output,".",dt=dt_sim,t0=0, fmin=f1, fmax=f2, st2_isref=True,show=arg_show)

if __name__ == "__main__":
    main(sys.argv[1:])

