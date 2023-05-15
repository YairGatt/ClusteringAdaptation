#!/usr/bin/env python

"""
This script corrects for multiple hypotheses using the Benjamini-Hochberg method
"""

import os,sys
from statsmodels.stats.multitest import multipletests

def correct_multiple_hypotheses(p_value_dict, alpha=0.05):
	"""
	Correct a dictionary of p_values in the format of {object:p_value} using the Benjamini Hochberg method
	"""
	#initialize adjusted p values dict
	padjust_dict = {}
	#run multiple hypothesis correction in python using the Benjamini-Hochberg method
	truths, padjust, correctedp1, correctedp2 = multipletests(p_value_dict.values(), alpha*2, method = 'fdr_bh')
	#define p_adjust dict
	for index, pathway in enumerate(p_value_dict.keys()):
		padjust_dict[pathway] = padjust[index]
	#return results
	return padjust_dict

def correct_multiple_hypotheses_dict_in_dict(p_value_dict, alpha=0.05):
        """
        Correct a dictionary of p_values in the format of {object1:object2:p_value} using the Benjamini Hochberg method
        """
        #initialize adjusted p values dict
        padjust_dict = {}
        #run multiple hypothesis correction in python using the Benjamini-Hochberg method
        vector = []
        for i in p_value_dict:
            for j in p_value_dict[i]: vector.append(p_value_dict[i][j])
        #correct
        truths, padjust, correctedp1, correctedp2 = multipletests(vector, alpha*2, method = 'fdr_bh')
        #define p_adjust dict
        index = 0
        for i in p_value_dict:
            padjust_dict[i] = {}
            for j in p_value_dict[i]:
                padjust_dict[i][j] = padjust[index]
                index += 1
        #return results
        return padjust_dict

def main(argv=None):
	#This should not be called from teh command line
	raise Exception("This script is not intended to be called from the command line")

if __name__ == "__main__":
		exit(main())
