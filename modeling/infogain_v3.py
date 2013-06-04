#
#   Module: infogain_v3.py
#	Author: Jie Shen  @ NCTR/FDA
#	Date: Jul. 12, 2012.
#
#       Due to the low efficiency of the old version, this one is optimized. it combines the t-test and fc in one
#       function and using the Xposi and Xnega lists as input.
#
#
##############################################
#   Module: infogain_v2.py
#	Author: Jie Shen  @ NCTR/FDA
#	Date: Jun. 27, 2012.
#   Modified: Jul 3. 2012
#
#   Added a method to return the T test P value of the feature
#
#	Description:
#		Module to define the class of Feature (gene, transcripts,...)
#       The purpose is to calculate the INFORMATION GAIN of a giving feature
#
#	Dependence: math.py, numpy.py,
#
#	Usage: f=Feature(No., Name, Y, X, Bins, NOTES)
#           where No. is the order number of the feature, Name is the name of the feature
#           Y is the ordered list of the classification labels, X is the ordered list of the feature values
#           the order of Y and X are exactly mapped. Bins is the range to category the distribution of
#           feature values.
#
#          Information Gain: f.ig()
#
#
###################Old Description for infogain.py###################
#	Module: selected.py
#	Author: Jie Shen
#	Date: Constructed on Dec. 16, 2009; Modified on Jan. 4, 2010.
#	
#	Description: 
#		Module to define classes used in INFORMATION GAIN 
#		of features of molecular fingerprint
#
#	Dependence: math.py, numpy.py,
#	
#	Usage: f=Feature(No., N_of_Class_1, N_of_Class_1_containing_feature_t, \
#				N_of_Class_0, N_of_Class_0_containing_feature_t, NOTES)
#
#	# Pay Attention!!!!
# 	This module is wrong !!! until it was found on Jan. 20, 2010
#	Using the wrong formulations!!!!!!!!!
#
#	Corrected by Jie Shen
#		Jan. 20, 2010
#
#######################################################################3
#
#!/bin/python	




#declare class
from math import *
from numpy import histogram

#define a feature

class Feature:
    def __init__(self, no, name, Y, X, BinArray, notes):

        self.no=no						#Public. Feature Order Number
        self.name=name					#Feature Name
        self.Y=Y                        #sample class labels
        self.X=X                        #feature values of samples
        self.BinArray=BinArray          #bins
        self.notes=notes                #notes
        Xp=[]
        Xn=[]
        c1=0;c0=0
        for i in range(0,len(Y)):
            if Y[i]==1:
                c1=c1+1
                Xp.append(X[i])
            else:
                c0=c0+1
                Xn.append(X[i])
        self.n_c1=c1
        self.n_c0=c0
        self.Xposi=Xp
        self.Xnega=Xn


#    def hc(self):
#        p_c1=self.n_c1*1.0/(self.n_c1+self.n_c0)	#Probability of the first category
#        p_c0=1-p_c1									#Probability of the second category
#        hc=-p_c1*log(p_c1)/log(2)-p_c0*log(p_c0)/log(2)
#        return hc								#Entropy of the classification system (data set)
    def hc(self,n_c1,n_c0):
        try:
            p_c1=n_c1*1.0/(n_c1+n_c0)	                #Probability of the first category
            p_c0=1-p_c1									#Probability of the second category
            hc=-p_c1*log(p_c1)/log(2)-p_c0*log(p_c0)/log(2)
        except:
            hc=0
        return hc

    def hct(self):
#        print "called hct"
        hist=histogram(self.X,self.BinArray)
#        print self.X
#        print hist
        P=hist[0]                  #probility of the values located into each bin
        C=hist[1]                  #bin bound values
        hct_sum=0
        total=len(self.X)

        for i in range(0,len(C)-1):
            # calculate conditional entropy
            ct1=0;ct0=0
            for j in range(0,total):
                if self.X[j]>=C[i] and self.X[j]<C[i+1]:
                    if self.Y[j]==1:
                        ct1=ct1+1
                    else:
                        ct0=ct0+1

            # sum each conditional entropy
#            print ct1,ct0,self.hc(ct1,ct0),P[i]/total
            hct_sum=hct_sum+P[i]*self.hc(ct1,ct0)/total

        return hct_sum


    def ig(self):
        return self.hc(self.n_c1,self.n_c0)-self.hct()					#Information Gain of feature t

###########ADDED Methods to return the T-test pvalue of the feature
    def pvalue(self):
        from scipy.stats import ttest_ind
        return ttest_ind(self.Xposi,self.Xnega)[1]


###########ADDED Methods to return the fold change of the feature
    def fc(self):
        from numpy import mean
        return abs(mean(self.Xposi)-mean(self.Xnega))


