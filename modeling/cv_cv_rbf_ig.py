#
#   Module: cv_cv_rbf_ig.py
#	Author: Jie Shen  @ NCTR/FDA
#	Date: Sep. 4, 2012.
#
#   Description:
#       This module is based on the modeling_rbf_ig_v2.py. The purpose is to do "5-CV" cross validation
#       on the TGX Training data set. This is specifically for NGS data.
#############################
#   Module: ngs_cv_rbf_ig.py
#	Author: Jie Shen  @ NCTR/FDA
#	Date: Aug. 30, 2012.
#
#   Description:
#       This module is based on the modeling_rbf_ig_v2.py. The purpose is to do "leave compound out" cross validation
#       on the TGX data set. This is specifically for NGS data.
############################
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
    global labelMatrix, dataMatrix, sampletitles, featuretitles#, testdatafile, testdataMatrix
    global param_n, param_c, param_g, param_fset
    global complistfile, complist, compmapmatrix

    usage = """\
Usage: modeling_rbf_ig.py [-n begin,end,step] [-log2c begin,end,step] [-log2g begin,end,step] [-v fold]
-i labelfile -j datafile -d outputDIR"""

    fold=5
    maxinnerrun=10
    maxruns=10

    n_begin, n_end, n_step = 4, 100, 2
    c_begin, c_end, c_step = -5,  15, 2
    g_begin, g_end, g_step =  5, -15, -2
    is_win32 = (sys.platform == 'win32')
    if is_win32:
        Path="C:\\Documents and Settings\\jshen\\My Documents\\Research\\SEQC_TGx\\PredTGx\\"
        outpath="C:\\Documents and Settings\\jshen\\My Documents\\Research\\SEQC_TGx\\PredTGx\\"
    else:
        Path="/home/hhong/seqc_tgx/data/"
        outpath="/home/hhong/seqc_tgx/model/linear/y1/"

    datafile=Path+"train44_cleaned_ratio_log2.txt"
    labelfile=Path+"YAhR.txt"
#    complistfile=outpath+"compmapping.txt"
#    testdatafile=Path+"test36_cleaned_ratio_log2.txt"

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
#        elif argv[i] == '-t':
#            i += 1
#            testdatafile = argv[i]
#        elif argv[i] == '-l':
#            i += 1
#            complistfile = argv[i]
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


#    compmapmatrix=readmap_str(complistfile)
#    templist=[]
#    for dataline in compmapmatrix:
#        templist.append(dataline[0])
#    complist=list(set(templist))


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

def modelbuild(training, testing):
    # training and testing are the matrix of the data sets from label matrix without the title.

    #several parameters should be defined using this module:
    # n_end, n_begin, n_step
    # fold
    # maxinnerrun
    # c_begin, c_end, c_step
    # g_begin, g_end, g_step
    #
    returnlist=[]

    intrY_array=[]
    intrX_array=[]
    inteY_array=[]
    inteX_array=[]

    trainingY,trainingX=getSVMinput(training)
    testingY,testingX=getSVMinput(testing)
    selectedtrainingX,selectedfno=selectX(trainingY,trainingX,n_end if n_end>n_begin else n_begin)

    #generate the cross validation data sets:
    # Array[list1[],list2[],...,listn[]]#

    for innerrun in range(0,maxinnerrun):
        intrY_list=[]
        intrX_list=[]
        inteY_list=[]
        inteX_list=[]
        for innertraining,innertesting in k_fold_cross_validation(training,fold):
            innertrainingY,innertrainingX=getSVMinput(innertraining)
            intrY_list.append(innertrainingY)
            intrX_list.append(innertrainingX)
            innertestingY,innertestingX=getSVMinput(innertesting)
            inteY_list.append(innertestingY)
            inteX_list.append(innertestingX)
        intrY_array.append(intrY_list)
        intrX_array.append(intrX_list)
        inteY_array.append(inteY_list)
        inteX_array.append(inteX_list)

    #parameter grid search
    maxacc=0
    outacc=open(outpath+"acc.txt","w")
    for n in range_f(n_begin,n_end,n_step):
        for c in permute_sequence(range_f(c_begin, c_end, c_step)):
            for g in permute_sequence(range_f(g_begin, g_end, g_step)):
                param=svm_parameter("-t 2 -c "+str(2.0**c)+" -g "+str(2.0**g)+" -b 1")
                acc_list=[]
                print("n=%d, c=%d, g=%d\n" %(n,c,g))
                for innerrun in range(0,maxinnerrun):
                    for i in range(0,fold):
                        innertrainingY=intrY_array[innerrun][i]
                        innertrainingX=intrX_array[innerrun][i]
                        innertestingY=inteY_array[innerrun][i]
                        innertestingX=inteX_array[innerrun][i]

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

    returnlist.append(maxacc)
    returnlist.append([param_n,param_c,param_g])
    returnlist.append(param_fset)

    # model building according to the parameters got from grid search


    selectedtrainingX=genX(trainingX,param_fset)
    param=svm_parameter("-t 2 -c "+str(2.0**param_c)+" -g "+str(2.0**param_g)+" -b 1")
    prob=svm_problem(trainingY,selectedtrainingX)
    m=svm_train(prob,param)
    svm_save_model(outpath+"svm_model.model",m)

    selectedtestingX=genX(testingX,param_fset)


    p_label, p_acc, p_val=svm_predict(testingY,selectedtestingX,m,'-b 1')
    probset=[]
    for i in range(0,len(testingY)):
        if p_label[i]*(p_val[i][0]-0.5)<0:
            probset.append(1-p_val[i][0])
        else:
            probset.append(p_val[i][0])
    returnlist.append(p_label)
    returnlist.append(probset)
    returnlist.append(m)
    return returnlist
    # it will return a list containing
    # [maxACC in the CV, [n,c,g], [features], [predict labels], [predict probs], SVM model]



def main(argv):

# Give the default parameters
    """

    """

    process_options()

    # script begin

    out_feature=open(outpath+"features.txt","w")
    out_pred=open(outpath+"pred.txt","w")


    out_pred.write("sample")
    allresults=[]
    for run in range(0,maxruns):
        out_pred.write("\tCV%d" % run+1)
        #initialize cvresult
        cvresult=[]
        for i in range(0,len(sampletitles)):
            cvresult.append(0)
        ################
        for training,testing in k_fold_cross_validation(labelMatrix,fold):
            results=modelbuild(training,testing)
            params=results[1]
            out_feature.write("CV%d, (n=%d,c=%d,g=%d)\t" % (run, params[0],params[1],params[2])),
            for fno in results[2]:
                out_feature.write(str(featuretitles[fno])+" "),
            out_feature.write("\tMAX Average Accurate: %f\n" % results[0])
            for i in range(0,len(testing)):
                id=sampletitles.index(testing[i][0])
                cvresult[id]=results[4][i]
        allresults.append(cvresult)
    out_pred.write("\n")

#output the results
    for i in range(1,len(sampletitles)):
        out_pred.write("%s" % sampletitles[i])
        for j in range(0,maxruns):
            out_pred.write("\t%f" % allresults[j][i])
        out_pred.write("\n")

    out_feature.close()
    out_pred.close()



################
if __name__ == "__main__":
    main(sys.argv[1:])
