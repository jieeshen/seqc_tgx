#
#	Module: collectData.py
#	Author: Jie Shen
#	Date: 2012.06.20
#	Version 0.01
#
#
#
#	Description:
#		This module is used to collect data from different files and compose of an ensemble table.
#       Specifically, the data from NGS are separated by samples. There is also a fileName file which indicates
#       the samples identification.
#
#
#	Dependence:
#
#	Usage:
#
#!/bin/python

Path="C:\\Documents and Settings\\jshen\\My Documents\\Research\\SEQC_TGx\\PredTGx\\validation\\"
ListName=Path+"SampleList.txt"
DataName1="sum"
DataName3=".txt"
OutFileName=Path+"validation.txt"
OutControlFile=Path+"control.txt"

RepeatedList=["546","562","578","542","558","574"]

ControlList=["547","536","563","568","552","531"]

lfile=open(ListName, "r")
list=lfile.readlines()

DataMatrix=[]
nd=0
ControlMatrix=[]
nc=0
for i in range(0,len(list)/2):
    sampleid=list[i*2].strip("\n").strip("\r")
    if sampleid in ControlList:
        line=[sampleid]
        fileid="%03d" % (i+1)
        inp=open(Path+DataName1+fileid+DataName3,"r")
        datamatrix=inp.readlines()
        inp.close()
        genetitles=["RefSeqID"]
        for dataline in datamatrix:
            datas=dataline.split("\t");
            genetitles.append(datas[0])
            line.append(datas[2])
        ControlMatrix.append(line)
        nc=nc+1

    elif sampleid not in RepeatedList:
        line=[sampleid]
        fileid="%03d" % (i+1)
        inp=open(Path+DataName1+fileid+DataName3,"r")
        datamatrix=inp.readlines()
        inp.close()
        genetitles=["RefSeqID"]
        for dataline in datamatrix:
            datas=dataline.split("\t");
            genetitles.append(datas[0])
            line.append(datas[2])
        DataMatrix.append(line)
        nd=nd+1



outp=open(OutFileName,"w")

for i in range(0, len(genetitles)-1):
    outp.write("%s\t" % genetitles[i])
    for j in range(0,nd):
        outp.write("%s\t" % DataMatrix[j][i])
    outp.write("\n")

outp.close()

outc=open(OutControlFile,"w")

for i in range(0, len(genetitles)-1):
    outc.write("%s\t" % genetitles[i])
    for j in range(0,nc):
        outc.write("%s\t" % ControlMatrix[j][i])
    outc.write("\n")

outc.close()











