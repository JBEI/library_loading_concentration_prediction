import matplotlib.pyplot as plt
import numpy as np

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