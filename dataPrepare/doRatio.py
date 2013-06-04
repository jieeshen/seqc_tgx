#
#   Module: doRatio
#   Author: JShen
#   Date: 8/16/12
#   Time: 2:55 PM
#   Version: 
#
#
#	Description:  This module is processing the ensemble matrix with training, training controls, test, testing controls
#                   It will reads a matchfile which gives clues for sample:control relationships and do the ratio.
#                   It will return two matrix files (training and test)
#
#
#	Dependence:
#	
#	Usage: 
#
#!/bin/python
import math

Path="C:\\Documents and Settings\\jshen\\My Documents\\Research\\SEQC_TGx\\PredTGx\\"
InputFile=Path+"cleaned_RPM_ensemble.txt"
MatchList=Path+"matchtype.txt"
TrainingDataFile=Path+"training44_ratio_log2.txt"
TestingDataFile=Path+"testing36_ratio.txt"

mlistfile=open(MatchList,"r")
mlistdata=mlistfile.readlines()
matchdic={}
for line in mlistdata:
    datas=line.strip("\n").strip("\r").split("\t")
    matchdic[datas[0]]=datas[1]


inp=open(InputFile,"r")
dataMatrix=inp.readlines()
inp.close()
sampletitles=dataMatrix[0].strip("\n").strip("\r").split("\t")
newMatrix=[sampletitles[:-4]]

min2=1

for dataline in dataMatrix[1:]:
    datalist=dataline.strip("\n").strip("\r").split("\t")
    nline=[datalist[0]]
    for i in range(1,len(datalist)-4):
        controlname=matchdic[sampletitles[i]]
        controlid=sampletitles.index(controlname)
        dataratio=float(datalist[i])/float(datalist[controlid])
        nline.append(dataratio)
        if dataratio>0 and dataratio<min2:
            min2=dataratio
    newMatrix.append(nline)

print min2
min1=min2/1.2

outp=open(TrainingDataFile,"w")
for dataline in newMatrix:
    for data in dataline:
        try:
            value=float(data)
            if value<0.00000001:
                outp.write("%s\t" % math.log(min1,2))
            else:
                outp.write("%s\t" % math.log(data,2))
        except:
            outp.write("%s\t" % data)
    outp.write("\n")

outp.close()









