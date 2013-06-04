#
#   Module: AffyDataProc.py
#   Author: JShen
#   Date: 8/20/12
#   Time: 9:16 PM
#   Version: 
#
#
#	Description:  This module process the Affy MicroArray Data. According to specified input, e.g. the repeat sample
#               list, the control sample list, the outlier list ..., it will generate a training/testing data matrix.
#               for the data, it will calculate the control mean and do the ratio.
#
#               Since the data has been transformed using log2, the control mean is calculated using
#               log2[(2^d1+2^d2+2^d3)/3] and the ratio is calculated using d1-(control mean)
#
#
#	Dependence:
#	
#	Usage: 
#
#!/bin/python
import getopt
import string
import sys
import math

def usage(name):
    print "\n"
    print "\tUSAGE\n"
    print "\t\t", name, "-i INPUTFILE -c CONTROLLIST -t OUTLIERLIST -r REPEATLIST -o OUTPUTFILE "
    print """
		PARAMETERS

			-h, --help	: Help
			-i,		: input text file, contain a matrix
			-c,     : match list file
			-t,     : outlier list, delimited using ","
			-r,     : repeat list, delimited using ","
			-o,		: output file
		"""
    sys.exit()

if __name__ == '__main__':
    is_win32 = (sys.platform == 'win32')
    if is_win32:
        outpath="C:\\Documents and Settings\\jshen\\My Documents\\Research\\SEQC_TGx\\AffyTGx\\Results\\"
    else:
        outpath="/home/hhong/seqc_tgx/AffyTGx/data/"

    #Define Control lists
    NNIPlist=["95865.CEL", "95836.CEL", "96606.CEL", "96211.CEL", "96483.CEL", "96496.CEL"]
    NNOGlist=["95716.CEL", "95841.CEL", "95893.CEL", "95793.CEL", "96006.CEL", "95375.CEL"]
    NUOGlist=["95686.CEL", "95746.CEL", "95915.CEL", "95797.CEL", "96706.CEL", "96729.CEL"]

    outlierlist=["96653.CEL"]

    repeatlist=[]

    #Define files
    infile=outpath+"SEQC_TGMX_Training_set-MAS5_Normalized.txt"
    outfile=outpath+"Training44.txt"
    matchlistfile=outpath+"Affymatchtype.txt"

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:c:t:r:o:",["help"])
    except getopt.GetoptError:
        usage(sys.argv[0])

    logfile=outpath+"log.txt"

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(sys.argv[0])
        elif opt in ("-i"):
            infile=arg
        elif opt in ("-c"):
            matchlistfile=arg
        elif opt in ("-t"):
            outlierlist=arg.split(",")
        elif opt in ("-r"):
            repeatlist=arg.split(",")
        elif opt in ("-o"):
            outfile=arg

    source="".join(args)

    mlistfile=open(matchlistfile,"r")
    mlistdata=mlistfile.readlines()
    matchdic={}
    for line in mlistdata:
        datas=line.strip("\n").strip("\r").split("\t")
        matchdic[datas[0]]=datas[1]

    inp=open(infile,"r")
    matrix=inp.readlines()
    inp.close()
    sampletitles=matrix[0].strip("\r").strip("\n").strip("\t").split("\t")
    outp=open(outfile,"w")
    datasampletitles=["Probe Set ID"]
    outp.write("%s" % "Probe Set ID")
    for title in sampletitles[1:]:
        if title in (NNIPlist+NNOGlist+NUOGlist):
            pass
        elif title in outlierlist:
            pass
        elif title in repeatlist:
            pass
        else:
            datasampletitles.append(title)
            outp.write("\t%s" % title)
    outp.write("\n")

    for line in matrix[58:]:
        dataline=line.strip("\r").strip("\n").strip("\t").split("\t")
        outp.write("%s" % dataline[0])
        nniptotal=0
        for datatitle in NNIPlist:
            nniptotal=nniptotal+2**float(dataline[sampletitles.index(datatitle)])
        nnipmean=math.log(nniptotal/len(NNIPlist),2)
        nnogtotal=0
        for datatitle in NNOGlist:
            nnogtotal=nnogtotal+2**float(dataline[sampletitles.index(datatitle)])
        nnogmean=math.log(nnogtotal/len(NNOGlist),2)
        nuogtotal=0
        for datatitle in NUOGlist:
            nuogtotal=nuogtotal+2**float(dataline[sampletitles.index(datatitle)])
        nuogmean=math.log(nuogtotal/len(NUOGlist),2)

        for datatitle in datasampletitles[1:]:
            if matchdic[datatitle]=="NNIP":
                newdata=float(dataline[sampletitles.index(datatitle)])-nnipmean
                outp.write("\t%s" % newdata)
            elif matchdic[datatitle]=="NNOG":
                newdata=float(dataline[sampletitles.index(datatitle)])-nnogmean
                outp.write("\t%s" % newdata)
            elif matchdic[datatitle]=="NUOG":
                newdata=float(dataline[sampletitles.index(datatitle)])-nuogmean
                outp.write("\t%s" % newdata)
            else:
                print("something is wrong.............")
        outp.write("\n")
    outp.close()