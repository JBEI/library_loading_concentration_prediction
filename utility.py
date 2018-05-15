import pandas as pd
from IPython.display import display
import os,glob
import numpy as np

#Peak Detection
from scipy.interpolate import interp1d
from scipy.integrate import quad
import peakutils

#Machine Learning
from tpot import TPOTRegressor
from sklearn.model_selection import cross_val_predict
from sklearn.ensemble import RandomForestRegressor

#Plotting Tools
#if __name__ == "__main__":
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns

test=False

def aproximate_sequence(df,n=100,start=0,stop=15000):
    """Aproximate a discrete sequence using a truncated expansion."""
    
    #Create even spacing
    int_lims = np.linspace(start,stop,n+1)
    
    #Create Interpolation function
    abundance_curve = interp1d(df['base pairs'],df['Value'],fill_value='extrapolate')
    
    #Aproximate Function using an integral
    state = [quad(abundance_curve,a,b,limit=2000)[0] for a,b in zip(int_lims[:-1],int_lims[1:])]
    bps = int_lims[:-1]
    
    return state,bps

def extract_features(df):
    pass
    

def find_peaks(df,test=False):
    """Find the peaks of a ladder gel. Return a list of peak times and the 
    base pairs that correspond to those times.
    """
    
    #display(df)
    peaks = peakutils.indexes(df['Value'].values, thres=0.15, min_dist=10)
    #print(peaks)
    #plt.figure(figsize=(10,10))
    #df['Value'].plot()
    #plt.scatter(peaks,np.zeros(len(peaks)))
    #plt.show()
    peak_bps = [35,50,150,300,400,500,600,700,1000,2000,3000,7000,10380]
    peak_times = df['Time'].values[peaks[1:-1]]
    #print(peak_times)
    
    if test:
        plt.figure()
        df.plot('Time','Value')
        plt.scatter(peak_times,np.zeros(len(peak_times)))
        plt.show()
        print(len(peaks),len(peak_bps))
    
    return peak_times,peak_bps

def normalize_peak(df,start_bp=500,end_bp=11000):
    
    #Subtract out Noise Floor
    noise_floor = np.mean(df['Value'].tail(20))
    df['Value'] -= noise_floor
    df.loc[df['Value']<0,'Value']=0
    
    #Select only important segment
    LIBRARY_SEGMENT = (df['base pairs'] > start_bp) & (df['base pairs'] < end_bp)
    df = df.loc[LIBRARY_SEGMENT]
    
    #Normalize Values to have an area under the curve of 1
    area = np.trapz(df['Value'],x=df['base pairs'])
    df.loc[:,'Value'] /= area
    
    return df

def load_data(file_id,plot=False):
    """load bioanalyzer data from a set of two csvs. Return a pandas series."""
    for file in glob.glob(file_id+'*'):
        #Load Ladder Data
        if 'Ladder' in file:
            ladder_df = pd.read_csv(file,skiprows=17)[:-1]
            ladder_df = ladder_df.apply(np.vectorize(float))
            peak_times,peak_bps = find_peaks(ladder_df,test=test)
            #print(peak_times,peak_bps)
            #print(len(peak_times),len(peak_bps))
            if len(peak_times) != len(peak_bps):
                return None
            
            time_to_bps = interp1d(peak_times,peak_bps,fill_value='extrapolate',kind='slinear')
            
        #Load Sample Data
        elif 'Sample' in file:
            sample_df = pd.read_csv(file,skiprows=17)[:-2]
            sample_df = sample_df.apply(np.vectorize(float))
            sample_df['base pairs'] = sample_df['Time'].apply(time_to_bps)
            for column in ['Library Loading Concentration (pM)','Cluster Density (K/mm^2)']:
                if column not in sample_df.columns:
                    sample_df[column] = float('NaN')
        
        #plt.scatter(peak_times,peak_bps)
        #plt.plot(np.arange(0,120),np.vectorize(time_to_bps)(np.arange(0,120)))
        #plt.show()
    
    
    llc = sample_df['Library Loading Concentration (pM)'].values[0]
    cluster_density = sample_df['Cluster Density (K/mm^2)'].values[0]
    sample_df = sample_df.loc[:,sample_df.columns.isin(['Time','Value','base pairs'])]
    
    #Normalize Peak 
    sample_df = normalize_peak(sample_df)
    
    #Plot Base Pairs Vs Values
    #sample_df['logbp'] = np.log10(sample_df['base pairs'])
    if plot:
        sample_df.plot('base pairs','Value')
        plt.show()
    
    #sample_df.plot('logbp','Value')
    #plt.show()
    
    #sample_df.plot('Time','Value')
    #plt.show()
    
    #sample_df.plot('Time','base pairs')
    #plt.show()
    
    #display(sample_df)
    state,bps = aproximate_sequence(sample_df,n=10,stop=12000)    
    
    #Create df
    columns = ['Library Loading Concentration','Cluster Density'] + [str(int(bp)) for bp in bps]
    data = [[llc,cluster_density] + state]
    
    state_df = pd.DataFrame(data,columns=columns)
    
    return state_df


def predict_loading_concentration(spectra,model,low=5,high=20,n=1000,output_file=None):
    cluster_densities = []
    loading_concentrations = np.linspace(low,high,n)
    for lc in loading_concentrations:
        X = np.append(lc,spectra).reshape(1,-1)
        y = model.predict(X)
        cluster_densities.append(y[0])
    
    #print(cluster_densities)
    plt.plot(loading_concentrations,cluster_densities)
    plt.title('Loading Concentration vs Predicted Cluster Density')
    plt.xlabel('Loading Concentration [pM]')
    plt.ylabel('Predicted Cluster Density [k/mm^2]')
    plt.xlim([low,high])
    
    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file,dpi=600)