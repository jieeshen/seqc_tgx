#
#   Module: modeling_rbf_ig_v3.py
#	Author: Jie Shen  @ NCTR/FDA
#	Date: Sep. 02, 2012.
#
#   Description:
#       This update version updates several issues:
#       1. using "leave compound out" in the inner CV process.
#
#############################
#   Module: modeling_rbf_ig_v2.py
#	Author: Jie Shen  @ NCTR/FDA
#	Date: Aug. 24, 2012.
#
#   Description:
#       This update version updates several issues:
#       1. in information gain calculation, the bin number were fixed as 2
#       2. Information Gain ranking were calculated before grid search.#
#
############################
#   Module: modeling_p.py
#	Author: Jie Shen  @ NCTR/FDA
#	Date: Aug. 17, 2012.
#
#   Description:
#       This module is derived from cv_linear_p_v2.py. The difference is there are no cross validation performed on
#       the training set. It will reads the training and testing files and builds a modeling using training set. The
#       parameters including number of variables, C,(g) were determined in a cross validation process. it will output
#       the models and used feature sets and prediction results on test set.#
#
###########################
#
#   Module: cv_linear_p_v2.py
#	Author: Jie Shen  @ NCTR/FDA
#	Date: Jul. 14, 2012.
#
#   Description:
#       Previous version only run 1 time in the inner cv to determine the optimized parameter.
#       In this upgraded version v2, it runs RUN_INNER times in the inner cv process
#
###########################
#
#   Module: cv_linear_p.py
#	Author: Jie Shen  @ NCTR/FDA
#	Date: Jul. 14, 2012.
#
#	Description:
#		The previous v5 version is using inner k fold cv to obtain the svm parameters including "number of features"
#       and the C in SVM. This version is using the selected feature set replacing the "number of features"
#
#		This is totally different than previous versions. It will do svm training in the program by calling the
#       libSVM python module. It will combine the feature selection, grid search together.
#
#       it will run M times K-fold cross validation, in each fold in each run, the training set will do a inner K-fold
#       cross validation to determine the parameters, including feature sets, c (parameters in libSVM),
#       then it will build a model of the entire training set using the same parameters, to predict the testing set.
#
#       It uses information gain as the selected criteria
#       It will report the M*K models' features sets statistics and M predictive accuracy.
#
#	Dependence: infogain_v3.py
#
import sys
from numpy import *
from fileRead import *
from infogain_v3 import *
from svmutil import *


def k_fold_cross_validation(X, K, randomise=True):
    """
     Generates K (training, validation) pairs from the items in X.

     Each pair is a partition of X, where validation is an iterable
     of length len(X)/K. So each training iterable is of length (K-1)*len(X)/K.

     If randomise is true, a copy of X is shuffled before partitioning,
     otherwise its order is preserved in training and validation.
     """
    if randomise:
        from random import shuffle
        X=list(X)
        shuffle(X)
    for k in xrange(K):
        training = [x for i, x in enumerate(X) if i % K != k]
        validation = [x for i, x in enumerate(X) if i % K == k]
        yield training, validation

def process_options(argv=sys.argv):
    """

    """
    global maxruns, fold, maxinnerrun
    global n_begin, n_end, n_step
    global c_begin, c_end, c_step
    global g_begin, g_end, g_step
    global datafile, labelfile,outpath
    global labelMatrix, dataMatrix, sampletitles, featuretitles, testdatafile, testdataMatrix
    global param_n, param_c, param_g, param_fset
    global complistfile, complist, compmapmatrix

    usage = """\
Usage: modeling_rbf_ig.py [-n begin,end,step] [-log2c begin,end,step] [-log2g begin,end,step] [-v fold]
-i labelfile -j datafile -t testdatafile -d outputDIR"""

    fold=5
    maxinnerrun=10

    n_begin, n_end, n_step = 22, 20, -1
    c_begin, c_end, c_step = -15,  5, 2
    g_begin, g_end, g_step =  15, -5, -2
    is_win32 = (sys.platform == 'win32')
    if is_win32:
        Path="C:\\Documents and Settings\\jshen\\My Documents\\Research\\SEQC_TGx\\PredTGx\\"
        outpath="C:\\Documents and Settings\\jshen\\My Documents\\Research\\SEQC_TGx\\PredTGx\\"
    else:
        Path="/home/hhong/seqc_tgx/data/"
        outpath="/home/hhong/seqc_tgx/model/linear/y1/"

    datafile=Path+"train44_cleaned_ratio_log2.txt"
    testdatafile=Path+"test36_cleaned_ratio_log2.txt"
    labelfile=Path+"YAhR.txt"
    complistfile=outpath+"compmapping.txt"

    if len(argv) < 0:
        print(usage)
        sys.exit(1)

    pass_through_options=[]

    i = 1
    while i < len(argv) - 1:
        if argv[i] == "-m":
            i += 1
            maxruns = int(argv[i])
        elif argv[i] == "-n":
            i += 1
            (n_begin,n_end,n_step) = map(int,argv[i].split(","))
        elif argv[i] == "-log2c":
            i += 1
            (c_begin,c_end,c_step) = map(float,argv[i].split(","))
        elif argv[i] == "-log2g":
            i += 1
            (g_begin,g_end,g_step) = map(float,argv[i].split(","))
        elif argv[i] == "-v":
            i += 1
            fold = int(argv[i])
        elif argv[i] == '-i':
            i += 1
            labelfile = argv[i]
        elif argv[i] == '-j':
            i += 1
            datafile = argv[i]
        elif argv[i] == '-l':
            i += 1
            complistfile = argv[i]
        elif argv[i] == '-t':
            i += 1
            testdatafile = argv[i]
        elif argv[i] == '-d':
            i += 1
            outpath = argv[i]
        else:
            pass_through_options.append(argv[i])
        i += 1

    labelMatrix=readfile_int(labelfile)   #get label matrix without title line
    dataMatrix=readfile_float(datafile)     #get data matrix with title line
    sampletitles=dataMatrix[0]
    featuretitles=[]
    for line in dataMatrix[1:]:
        featuretitles.append(line[0])

    compmapmatrix=readmap_str(complistfile)
    templist=[]
    for dataline in compmapmatrix:
        templist.append(dataline[0])
    complist=list(set(templist))


def range_f(begin,end,step):
    # like range, but works on non-integer too
    seq = []
    while True:
        if step > 0 and begin > end: break
        if step < 0 and begin < end: break
        seq.append(begin)
        begin = begin + step
    return seq

def permute_sequence(seq):
    n = len(seq)
    if n <= 1: return seq

    mid = int(n/2)
    left = permute_sequence(seq[:mid])
    right = permute_sequence(seq[mid+1:])

    ret = [seq[mid]]
    while left or right:
        if left: ret.append(left.pop(0))
        if right: ret.append(right.pop(0))

    return ret


def getSVMinput(labelmatrix):              #for training set
    newtitles=["geneid"]
    newtitleindex=[0]
    Y=[]
    for line in labelmatrix:
        newtitles.append(line[0])
        Y.append(line[1])
        newtitleindex.append(sampletitles.index(line[0]))

    X=[]
    for i in newtitleindex[1:]:
        vectorX=[]
        for line in dataMatrix[1:]:
            vectorX.append(line[i])
        X.append(vectorX)

    return Y,X

def getTestSVMinput(testXfile): #this is specifically for test set data
    #it reads the file and return a Y and X list.
    testp=open(testXfile,"r")
    testX=testp.readlines()
    testsampletitles=testX[0].strip("\r").strip("\n").strip("\t").split("\t")
    Y=[]
    X=[]
    for i in range(1,len(testsampletitles)):
        Y.append(0)
        vectorX=[]
        for line in testX[1:]:
            datas=line.strip("\r").strip("\n").split("\t")
            vectorX.append(float(datas[i]))
        X.append(vectorX)
    return Y,X

def getTestTitles(testXfile):
    testp=open(testXfile,"r")
    testX=testp.readlines()
    testsampletitles=testX[0].strip("\r").strip("\n").strip("\t").split("\t")
    return testsampletitles




def genX(X,selectedfno):
    # this methods generate the X matrix  according to the feature selection result
    selectedX=[]
    for x in X:
        line=[]
        for i in selectedfno:
            line.append(x[i])
        selectedX.append(line)
    return selectedX

def selectX(Y,X,nf):
    datapool=[]
    for dataline in X:
        for data in dataline:
            datapool.append(data)
    #b=int(max(datapool)-min(datapool))+1
    b=2 #modified on 8/24/2012
    binArray=histogram(datapool,bins=b)[1]

    selected=[]
    for n in range(0,len(featuretitles)):
        featureX=[]
        for datas in X:
            featureX.append(datas[n])
        f=Feature(n, featuretitles[n], Y, featureX, binArray,"")
        featurepair=[f.ig(),f.no,f.name,f]
        selected.append(featurepair)
    selected.sort(reverse=True)

    selectedfno=[]
#    selectedfname=[]
    for fp in selected[:nf]:
        selectedfno.append(fp[1])
#        selectedfname.append(fp[2])

    selectedX=genX(X,selectedfno)
    return selectedX, selectedfno#, selectedfname





def main(argv):

# Give the default parameters
    """

    """

    process_options()

    # script begin

    out_feature=open(outpath+"features.txt","w")
    out_pred=open(outpath+"pred.txt","w")
    training=labelMatrix
    intrY_array=[]
    intrX_array=[]
    inteY_array=[]
    inteX_array=[]


    trainingY,trainingX=getSVMinput(training)
    selectedtrainingX,selectedfno=selectX(trainingY,trainingX,n_end if n_end>n_begin else n_begin)

    for comp in complist:
        trainingsamples=[]
        testingsamples=[]
        for line in compmapmatrix:
            if line[0]==comp:
                testingsamples.append(line[1])
            else:
                trainingsamples.append(line[1])
        innertraining=[]
        innertesting=[]
        for line in labelMatrix:
            if line[0] in testingsamples:
                innertesting.append(line)
            elif line[0] in trainingsamples:
                innertraining.append(line)
            else:
                print ("something is wrong!!!")
                exit(0)
        innertrainingY,innertrainingX=getSVMinput(innertraining)
        innertestingY,innertestingX=getSVMinput(innertesting)

        intrY_array.append(innertrainingY)
        intrX_array.append(innertrainingX)
        inteY_array.append(innertestingY)
        inteX_array.append(innertestingX)

    print "done partition...starting grid search..."
    #parameter grid search
    maxacc=0
    outacc=open(outpath+"acc.txt","w")
    for n in range_f(n_begin,n_end,n_step):
        for c in permute_sequence(range_f(c_begin, c_end, c_step)):
            for g in permute_sequence(range_f(g_begin, g_end, g_step)):
                param=svm_parameter("-t 2 -c "+str(2.0**c)+" -g "+str(2.0**g)+" -b 1")
                acc_list=[]
                print("n=%d, c=%d, g=%d\n" %(n,c,g))
                for innerrun in range(0,len(complist)):   #perform leave compound out cross validation...
                        innertrainingY=intrY_array[innerrun]
                        innertrainingX=intrX_array[innerrun]
                        innertestingY=inteY_array[innerrun]
                        innertestingX=inteX_array[innerrun]

                        selectedtrainingX=genX(innertrainingX,selectedfno[:n])
                        prob=svm_problem(innertrainingY,selectedtrainingX)
                        m=svm_train(prob,param)

                        selectedtestingX=genX(innertestingX,selectedfno[:n])

                        p_label, p_acc, p_val=svm_predict(innertestingY,selectedtestingX,m,'-b 1')
                        acc_list.append(p_acc[0])
                acc_ave=average(acc_list)
                outacc.write("%f\n" % acc_ave)
                if acc_ave>maxacc:
                    param_n=n
                    param_fset=selectedfno[:param_n]
                    param_c=c
                    param_g=g
                    maxacc=acc_ave
    outacc.close()

    print("Max Acc is %f, N of Features is %d" % (maxacc, len(param_fset)))
    # model building according to the parameters got from grid search

    testingY,testingX=getTestSVMinput(testdatafile)
    selectedtrainingX=genX(trainingX,param_fset)
    param=svm_parameter("-t 2 -c "+str(2.0**param_c)+" -g "+str(2.0**param_g)+" -b 1")
    prob=svm_problem(trainingY,selectedtrainingX)
    m=svm_train(prob,param)
    svm_save_model(outpath+"svm_model.model",m)

    selectedtestingX=genX(testingX,param_fset)


    testsampletitles=getTestTitles(testdatafile)
    p_label, p_acc, p_val=svm_predict(testingY,selectedtestingX,m,'-b 1')



    out_feature.write("(n=%d,c=%d,g=%d)\t" % (param_n,param_c,param_g)),
    for fno in param_fset:
        out_feature.write(featuretitles[fno]+" "),
    out_feature.write("\tMAX Average Accurate: %f\n" % maxacc)
    out_feature.close()

    out_pred.write("sample\tpredict_class\tprobability\n")
    for i in range(0,len(testingY)):
        out_pred.write("%s\t%d\t%f\n" % (testsampletitles[i+1],p_label[i],p_val[i][0]))

if __name__ == "__main__":
    main(sys.argv[1:])
