#!/usr/bin/python

import pickle
import pandas as pd
import numpy as np
from utility import predict_loading_concentration,find_peaks,normalize_peak,aproximate_sequence
from scipy.interpolate import interp1d
import argparse

#Allow Matplotlib to be used on commandline
import matplotlib as mpl
mpl.use('Agg')

#Supress Warnings from Integrator
import warnings
warnings.filterwarnings("ignore")


#Handle inputs
if __name__ == "__main__":
        
    #Parse Inputs
    parser = argparse.ArgumentParser(description='Use BioAnalyzer Data to Predict Library Loading Concentration for Library Barcoding')
    
    parser.add_argument('BAF',type=str, help='BioAnalyzer CSV')
    parser.add_argument('LAD',type=str, help='Ladder CSV')
    parser.add_argument('OF' ,type=str, help='Prediction Plot output file, *.png')
    args = parser.parse_args()
    
    #Load & Parse Ladder  
    ladder_df = pd.read_csv(args.LAD,skiprows=17)[:-1]
    ladder_df = ladder_df.apply(np.vectorize(float))
    peak_times,peak_bps = find_peaks(ladder_df)
    time_to_bps = interp1d(peak_times,peak_bps,fill_value='extrapolate',kind='slinear')

    #Check to make sure the right number of peaks were found in the ladder
    if len(peak_times) != len(peak_bps):
        raise('Could Not Find Peaks in Ladder File! Bug Zak to Increase the Robustness of the Peak Finding Algorithm')
    
    #Load & Parse Bioanalyzer Run
    sample_df = pd.read_csv(args.BAF,skiprows=17)[:-2]
    sample_df = sample_df.apply(np.vectorize(float))
    sample_df['base pairs'] = sample_df['Time'].apply(time_to_bps)
    sample_df = sample_df.loc[:,sample_df.columns.isin(['Time','Value','base pairs'])]
    sample_df = normalize_peak(sample_df)
    
    #========I Can Probably Remove This ============
    #for column in ['Library Loading Concentration (pM)','Cluster Density (K/mm^2)']:
    #    if column not in sample_df.columns:
    #        sample_df[column] = float('NaN')
            
    #llc = sample_df['Library Loading Concentration (pM)'].values[0]
    #cluster_density = sample_df['Cluster Density (K/mm^2)'].values[0]
    #===============================================
    
    #Create Features
    state,bps = aproximate_sequence(sample_df,n=10,stop=12000)    
    
    #Create df
    #state_df = pd.DataFrame(state,columns=[str(int(bp)) for bp in bps])
    
    #Load Model
    with open('model/model.pkl','rb') as fp:
        model = pickle.load(fp)
        
    predict_loading_concentration(state,model,output_file=args.OF)