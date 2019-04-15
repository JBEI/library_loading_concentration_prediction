#!/usr/bin/python
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

from diva_dna_seq.utility import aproximate_sequence,prepare_data
from diva_dna_seq import predict_loading_concentration
import argparse
import pickle



#Allow Matplotlib to be used on commandline
import matplotlib as mpl
mpl.use('Agg')

#Handle inputs
#Parse Inputs
parser = argparse.ArgumentParser(description='Use BioAnalyzer Data to Predict Library Loading Concentration for Library Barcoding')

parser.add_argument('BAF',type=str, help='BioAnalyzer CSV')
parser.add_argument('LAD',type=str, help='Ladder CSV')
parser.add_argument('OF' ,type=str, help='Prediction Plot output file, *.png')
args = parser.parse_args()

#Load Model
with open('../model/model30.pkl','rb') as fp:
    model = pickle.load(fp)

#Prepare Data
sample_df = prepare_data(args.LAD,args.BAF)

#Create Features
state,bps = aproximate_sequence(sample_df,n=10,stop=12000)    

#Predict
predict_loading_concentration(state,model,output_file=args.OF)