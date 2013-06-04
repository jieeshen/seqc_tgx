#
#	Module: fileRead.py
#	Author: Jie Shen @ NCTR, FDA
#	Date: 06/28/12
#	Version: 0.01
#		
#
#	Description:
#       This module define several useful methods to handle the input files
#		
#
#	Dependence:
#
#	
#	Usage: 
#
#
#
#!/usr/local/bin/python2.7

def readDic_int(labelfile):
    # read the label file and generate a dictionary
    ldic={}
    labelp=open(labelfile,"r")
    labelList=labelp.readlines()
    for line in labelList[1:]:
        items=line.strip("\n").strip("\r").split("\t")
        ldic[items[0]]=int(items[1])
    labelp.close()
    return ldic

def readfile_int(labelfile):
    # read the label file and generate a int matrix
    labelmatrix=[]
    labelp=open(labelfile,"r")
    labelList=labelp.readlines()
    for line in labelList[1:]:
        dataline=[]
        items=line.strip("\n").strip("\r").split("\t")
        dataline.append(items[0])
        dataline.append(int(items[1]))
        labelmatrix.append(dataline)
    labelp.close()
    return labelmatrix

def readfile_float(filename):
    # read a text file and return a matrix with all values are float except for some strings.
    matrix=[]
    readin=open(filename,"r")
    lines=readin.readlines()
    readin.close()
    for line in lines:
        datasStr=line.strip("\n").strip("\r").split("\t")
        datas=[]
        for data in datasStr:
            try:
                datas.append(float(data))
            except:
                datas.append(data)
        matrix.append(datas)
    return matrix
