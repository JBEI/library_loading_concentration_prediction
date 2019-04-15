#!/usr/bin/python
from diva_seq_opt.utility import aproximate_sequence,prepare_data
from diva_seq_opt import predict_loading_concentration
import diva_seq_opt
import argparse
import pickle
import numpy as np
import os

#Allow Matplotlib to be used on commandline
import matplotlib as mpl
mpl.use('Agg')



#Handle inputs
def main():    
    #Parse Inputs
    parser = argparse.ArgumentParser(description='Use BioAnalyzer Data to Predict Library Loading Concentration for Library Barcoding')
    
    parser.add_argument('BAF',type=str, help='BioAnalyzer CSV')
    parser.add_argument('LAD',type=str, help='Ladder CSV')
    parser.add_argument('OF' ,type=str, help='Prediction Plot output file, *.png')
    args = parser.parse_args()
    
    #Load Model
    model_path = os.path.dirname(os.path.realpath(diva_seq_opt.__file__))
    model_path = os.path.join(model_path,'model/model30.pkl')
    with open(model_path,'rb') as fp:
        model = pickle.load(fp)
    
    #Prepare Data
    sample_df = prepare_data(args.LAD,args.BAF)
    
    #Create Features
    state,bps = aproximate_sequence(sample_df,n=10,stop=12000)    
    
    #Predict
    predict_loading_concentration(state,model,output_file=args.OF)