#!/usr/bin/env python

import sys, os, random
import pickle
import matplotlib.pyplot as plt

def main():
    
    target_dir = '/home/rsamant2/scratch16-jgray21/rsamant2/pHLIP_structures'
    list_of_targets = target_dir + '/targets1.list'
    
    with open( list_of_targets, 'rt' ) as f: 
        test_cases = f.readlines()
        test_cases = [ x.strip() for x in test_cases ]
    print(test_cases)
    
    for name in test_cases:
        print(name)
        plt.figure()
        for i in range(5):
            filename = target_dir + "/" + name + "/result_model_" + str(i+1) + ".pkl"
            #print(filename)
            with open(filename,'rb') as f:
                data = pickle.load(f)
                
            X = list(range(0,len(data['plddt'])))
            Y = data['plddt']
            
            plt.plot(X, Y, label = "model_"+str(i))
            plt.legend()
            plt.xlabel('Residue #')
            plt.ylabel('predicted lddt')
            plt.savefig(name+'plddt_plot.png')
    
if __name__ == "__main__": main()
