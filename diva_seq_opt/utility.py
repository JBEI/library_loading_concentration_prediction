#Supress Warnings from Integrator (This could be made more specific)
import warnings
warnings.filterwarnings("ignore")


import pandas as pd
import numpy as np
import os,glob
import re

#Machine Learning
from tpot import TPOTRegressor
from sklearn.model_selection import cross_val_predict
from sklearn.ensemble import RandomForestRegressor

#Plotting Tools
import matplotlib.pyplot as plt

#Peak Detection
from scipy.interpolate import interp1d
from scipy.integrate import quad
import peakutils


TEST=False


def aproximate_sequence(df,n=100,start=0,stop=15000):
    """Reduce the Feature Space by Estimating a Curve with Its Integral"""
    
    #Create even spacing
    int_lims = np.linspace(start,stop,n+1)
    
    #Create Interpolation function
    abundance_curve = interp1d(df['base pairs'],df['Value'],fill_value='extrapolate')
    
    #Aproximate Function using an integral
    state = [quad(abundance_curve,a,b,limit=2000)[0] for a,b in zip(int_lims[:-1],int_lims[1:])]
    bps = int_lims[:-1]
    
    return state,bps
    

def find_peaks(df,test=False,peak_bps=[35,50,150,300,400,500,600,700,1000,2000,3000,7000,10380]):
    """Find the peaks of a ladder gel. Return a list of peak times and the 
    base pairs that correspond to those times.
    """
    
    for threshold in [0.30,0.20,0.15,0.10,0.05]:
        peaks = peakutils.indexes(df['Value'].values, thres=threshold, min_dist=10)
        peak_times = df['Time'].values[peaks[1:-1]]
    
        if test:
            print('testing!')
            plt.figure()
            df.plot('Time','Value')
            plt.scatter(peak_times,np.zeros(len(peak_times)))
            plt.show()
            print(len(peaks),len(peak_bps))
            print(len(peak_times),len(peak_bps))
            
        if len(peak_times) == len(peak_bps):
            break

    return peak_times,peak_bps


def normalize_peak(df,start_bp=500,end_bp=11000):
    '''Normalize BioAnalyzer Curve so its total area is 1.'''
    
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


def prepare_data(ladder_file,bioanalyzer_file):
    '''Convert a Ladder File and BioAnalyzer File into a machine learning data point'''
    
    #Load & Parse Ladder  
    ladder_df = pd.read_csv(ladder_file,skiprows=17)[:-1]
    ladder_df = ladder_df.apply(np.vectorize(float))
    peak_times,peak_bps = find_peaks(ladder_df)
    
    #Check to make sure the right number of peaks were found in the ladder
    if len(peak_times) != len(peak_bps):
        Warning('Could Not Find Peaks in Ladder File! Bug Zak to Increase the Robustness of the Peak Finding Algorithm')
        return None
    
    time_to_bps = interp1d(peak_times,peak_bps,fill_value='extrapolate',kind='slinear')
    
    #Read Bioanalyzer File
    sample_df = pd.read_csv(bioanalyzer_file,skiprows=17)[:-2]
    
    #Only Grab Relevant Columns
    sample_df = sample_df.loc[:,sample_df.columns.isin(['Time','Value'])]
    
    #Convert Values to Floats
    sample_df = sample_df.apply(np.vectorize(float))
    
    #Rescale From Retention Time to Base Pairs
    sample_df['base pairs'] = sample_df['Time'].apply(time_to_bps)
    
    #Normalize Curve
    sample_df = normalize_peak(sample_df)
    
    return sample_df


def get_sample_metadata(bioanalyzer_file):
    '''Extracts libary loading concentration and cluster density from the bioanalyzer file'''
    df = pd.read_csv(bioanalyzer_file,skiprows=17)[:-2]
    metadata = []
    for column in ['Library Loading Concentration (pM)','Cluster Density (K/mm^2)']:
        if column in df.columns:
            meta_val = df[column].values[0]
            if isinstance(meta_val,str):
                meta_val = float(meta_val.replace(',',''))
            metadata.append(meta_val)
        else:
            return None,None
    return metadata


def get_file_pairs(data_dir,
                   file_id_re='.*[?=\_]',
                   ladder_file_re='.*ladder',
                   bioanalyzer_fila_re='sample'
                  ):
    
    #Get list of all files in the training directory
    file_paths = glob.glob(data_dir+'/*')
    
    #Find File Pairs
    file_pairs = {}
    for file_path in file_paths:
        match = re.match(file_id_re,file_path).group(0)
        file_pairs.setdefault(match,[]).append(file_path)
    
    #create training file pairs
    training_files = []
    for _,file_pair in file_pairs.items():
        #Check if the First one is a ladder file
        match = re.match(ladder_file_re,file_pair[0],flags=re.I)
        if match is None:
            training_files.append([file_pair[1],file_pair[0]])
        else:
            training_files.append(file_pair)
    return training_files


def load_training_data(training_dir,plot=False):
    """load bioanalyzer training data from training directory."""
    training_files = get_training_files(training_dir)
    
    df = pd.DataFrame()
    for ladder_file,bioanalyzer_file in training_files:
        sample_df = prepare_data(ladder_file,bioanalyzer_file)
        llc,cluster_density = get_sample_metadata(bioanalyzer_file)
    
        if (sample_df is None) or (llc is None):
            continue
            
        #Plot Base Pairs Vs Values
        if plot:
            sample_df.plot('base pairs','Value')
            plt.show()
        
        state,bps = aproximate_sequence(sample_df,n=10,stop=12000)    
    
        #Create df
        columns = ['Library Loading Concentration','Cluster Density'] + [str(int(bp)) for bp in bps]
        data = [[llc,cluster_density] + state]
        df = pd.concat([df,pd.DataFrame(data,columns=columns)]) 
    return df


def fit_model(training_df,output_file):
    X = training_df.loc[:,df.columns != 'Cluster Density']
    y = training_df['Cluster Density']
    model = RandomForestRegressor().fit(X,y)
    
    with open(output_file,'wb') as fp:
        pickle.dump(model,fp)
    
    return model


def plot_model_performance(model,training_df):
    '''Create a Plot of the Models Cross Validated Performance'''
    
    X = training_df.loc[:,df.columns != 'Cluster Density']
    y = training_df['Cluster Density']
    y_p = cross_val_predict(model,X,y)

    plt.figure(figsize=(15,6))
    plt.subplot(1,2,1)
    plt.scatter(y,y_p)
    axes = plt.gca()
    ymin, ymax = axes.get_ylim()
    xmin, xmax = axes.get_xlim()
    min_val = min(xmin,ymin)
    max_val = max(ymax,xmax)
    plt.plot([min_val,max_val],[min_val,max_val],'k--')
    plt.xlim([min_val,max_val])
    plt.ylim([min_val,max_val])
    plt.title('Actual Vs Predicted Cluster Density')
    plt.xlabel('Actual Cluster Density')
    plt.ylabel('Predicted Cluster Density')

    plt.subplot(1,2,2)
    y_err = [ye -ye_p for ye,ye_p in zip(y,y_p)]
    print(min(y),max(y),max(y)-min(y))
    sns.distplot(y_err)
    plt.title('Prediction Error Histogram (Mean Error: {:0.0f})'.format(np.mean(np.abs(y_err))))
    plt.xlabel('Error in Predicting Cluster Density')

    plt.tight_layout()

    plt.savefig('../figures/ModelAccuracy.pdf',dpi=600)
    plt.show()