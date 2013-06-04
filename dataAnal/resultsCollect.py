#
#   Module: resultsCollect
#   Author: JShen
#   Date: 8/21/12
#   Time: 9:43 PM
#   Version: 0.01
#
#
#	Description: This module is used to collect the predict results of several models and generate a combined table
#
#
#	Dependence:
#	
#	Usage: resultsCollect.py -i file1 file2 file3 ... -o outfile
#
#!/bin/python
import getopt
import sys

def usage(name):
    print "\n"
    print "\tUSAGE\n"
    print "\t\t", name, "-i INPUTFILE -c CONTROLLIST -t OUTLIERLIST -r REPEATLIST -o OUTPUTFILE "
    print """
		PARAMETERS

			-h, --help	: Help
			-i,		: input text files, using ""
			-o,		: output file
		"""
    sys.exit()

if __name__ == '__main__':
    is_win32 = (sys.platform == 'win32')
    if is_win32:
        outpath="C:\\Documents and Settings\\jshen\\My Documents\\Research\\SEQC_TGx\\AffyTGx\\Results\\"
    else:
        outpath="/home/hhong/seqc_tgx/AffyTGx/data/"


    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["help"])
    except getopt.GetoptError:
        usage(sys.argv[0])

    logfile=outpath+"log.txt"

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(sys.argv[0])
        elif opt in ("-i"):
            infiles=arg.split(" ")
        elif opt in ("-o"):
            outfile=arg

    source="".join(args)

    inputdatalist=[]
    outp=open(outfile,"w")
    outp.write("sample")

    for infile in infiles:
        outp.write("\t%s" % infile[1:-12])
        inp=open(infile,"r")
        lines=inp.readlines()
        datamatrix=[]
        sampletitles=[]
        for line in lines[1:]:
            newline=[]
            datastrs=line.strip("\n").strip("\r").strip("\t").split("\t")
            sampletitles.append(datastrs[0])
            newline.append(int(datastrs[1]))
            if newline[0]*(float(datastrs[2])-0.5)<0:
                newline.append(1-float(datastrs[2]))
            else:
                newline.append(float(datastrs[2]))
            datamatrix.append(newline)
        inputdatalist.append(datamatrix)
        inp.close()

    resultmatrix=[]
    for i in range(0,len(sampletitles)):
        dataline=[]
        for j in range(0,len(infiles)):
            dataline.append(inputdatalist[j][i][1])
        sortedline=dataline[:]         #!!!!!!!!! if there is no ":", only the address were copied

        sortedline.sort(reverse=True)

        dataline.append(sortedline[0])  #Max
        dataline.append(sortedline[-1])  #Min
        datarange=sortedline[0]-sortedline[-1]
        dataline.append(datarange)  #Range
        datadiff=sortedline[0]-sortedline[1]
        dataline.append(datadiff)  #diff
        dataconf=datadiff/datarange
        dataline.append(dataconf) #conf
        if sortedline[0]>0.5:
            m1predict=infiles[dataline.index(sortedline[0])][1:-12]
            dataline.append(m1predict)
        else:
            dataline.append("")

        if dataconf>=0.5 and sortedline[0]>0.25:
            m2predict=infiles[dataline.index(sortedline[0])][1:-12]
            dataline.append(m2predict)
        else:
            dataline.append("")
        resultmatrix.append(dataline)

    outp.write("\tMax\tMin\tRange\tDiff\tConf\tModel1\tModel2\n")
    for i in range(0,len(sampletitles)):
        outp.write(sampletitles[i])
        for j in range(0,len(infiles)+7):
            outp.write("\t%s" % resultmatrix[i][j])
        outp.write("\n")
    outp.close()





