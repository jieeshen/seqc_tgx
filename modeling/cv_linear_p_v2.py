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
    global datafile, labelfile,outpath
    global labelMatrix, dataMatrix, sampletitles, featuretitles

    usage = """\
Usage: cv.py [-n begin,end,step] [-log2c begin,end,step] [-v fold]
-i labelfile -j datafile -d outputDIR"""

    fold=5
    maxruns=10
    maxinnerrun=5

    n_begin, n_end, n_step = 50, 500, 50
    c_begin, c_end, c_step = -15,  5, 2

    is_win32 = (sys.platform == 'win32')
    if not is_win32:
        datafile="/home/hhong/seqc_nb/processed-data/137_filtered.txt"
        labelfile="/home/hhong/seqc_nb/raw-data/137_label.txt"
        outpath="/home/hhong/seqc_nb/modeling/cv_v6_linear/"
    else:
        datafile="C:\Temp_data4\\137_filtered.txt"
        labelfile="C:\Temp_data4\\137_label.txt"
        outpath="C:\Temp_data4\\"

    if len(argv) < 0:
        print(usage)
        sys.exit(1)

    pass_through_options=[]

    i = 1
    while i < len(argv) - 1:
        if argv[i] == "-m":
            i += 1
            maxruns = argv[i]
        elif argv[i] == "-n":
            i += 1
            (n_begin,n_end,n_step) = map(int,argv[i].split(","))
        elif argv[i] == "-log2c":
            i += 1
            (c_begin,c_end,c_step) = map(float,argv[i].split(","))
        elif argv[i] == "-v":
            i += 1
            fold = argv[i]
        elif argv[i] == '-i':
            i += 1
            labelfile = argv[i]
        elif argv[i] == '-j':
            i += 1
            datafile = argv[i]
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

def range_f(begin,end,step):
    # like range, but works on non-integer too
    seq = []
    while True:
        if step > 0 and begin > end: break
        if step < 0 and begin < end: break
        seq.append(begin)
        begin = begin + step
    return seq

def getSVMinput(labelmatrix):
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
    b=int(max(datapool)-min(datapool))
    binArray=histogram(datapool,bins=b)[1]

    selected=[]
    for n in range(0,len(featuretitles)):
        featureX=[]
        for datas in X:
            featureX.append(datas[n])
        f=Feature(n, featuretitles[n], Y, featureX, binArray,"")
        featurepair=[f.pvalue(),f.no,f.name,f]
        selected.append(featurepair)
    selected.sort()

    selectedfno=[]
    for fp in selected[:nf]:
        selectedfno.append(fp[1])

    selectedX=genX(X,selectedfno)
    return selectedX, selectedfno





def main(argv):

# Give the default parameters
    """

    """
    global param_n, param_c, param_g, param_fset
    process_options()

    # script begin

    out_feature=open(outpath+"features.txt","w")
    out_acc=open(outpath+"acc.txt","w")
    acc_list=[]
    for runs in range(0,maxruns):
        #cross validation begins:
        count=0
        Y_list=[]
        predY_list=[]
        for training,testing in k_fold_cross_validation(labelMatrix,fold):
            count+=1
            maxacc=0
            #get inner cv datas
            intrY_array=[]
            intrX_array=[]
            inteY_array=[]
            inteX_array=[]
            fno_array=[]

            for innerrun in range(0,maxinnerrun):
                intrY_list=[]
                intrX_list=[]
                inteY_list=[]
                inteX_list=[]
                fno_list=[]
                for innertraining,innertesting in k_fold_cross_validation(training,fold):
                    innertrainingY,innertrainingX=getSVMinput(innertraining)
                    intrY_list.append(innertrainingY)
                    intrX_list.append(innertrainingX)
                    innertestingY,innertestingX=getSVMinput(innertesting)
                    inteY_list.append(innertestingY)
                    inteX_list.append(innertestingX)
                    selectedtrainingX,selectedfno=selectX(innertrainingY,innertrainingX,n_end)
                    fno_list.append(selectedfno)
                intrY_array.append(intrY_list)
                intrX_array.append(intrX_list)
                inteY_array.append(inteY_list)
                inteX_array.append(inteX_list)
                fno_array.append(fno_list)



            #parameter grid search
            for n in range_f(n_begin,n_end,n_step):
                for c in range_f(c_begin, c_end, c_step):
                    param=svm_parameter("-t 0 -c "+str(2.0**c))
                    acc_list=[]
                    for innerrun in range(0,maxinnerrun):
                        for i in range(0,fold):
                            innertrainingY=intrY_array[innerrun][i]
                            innertrainingX=intrX_array[innerrun][i]
                            innertestingY=inteY_array[innerrun][i]
                            innertestingX=inteX_array[innerrun][i]
                            selectedfno=fno_array[innerrun][i]

                            selectedtrainingX=genX(innertrainingX,selectedfno[:n])
                            prob=svm_problem(innertrainingY,selectedtrainingX)
                            m=svm_train(prob,param)

                            selectedtestingX=genX(innertestingX,selectedfno[:n])

                            p_label, p_acc, p_val=svm_predict(innertestingY,selectedtestingX,m)
                            acc_list.append(p_acc[0])
                    acc_ave=average(acc_list)
                    if acc_ave>maxacc:
                        param_n=n
                        param_fset=selectedfno[:param_n]
                        param_c=c
                        maxacc=acc_ave
            print("Max Acc is %f" % maxacc)
            # model building according to the parameters got from grid search
            trainingY,trainingX=getSVMinput(training)
            testingY,testingX=getSVMinput(testing)
            selectedtrainingX=genX(trainingX,param_fset)
            param=svm_parameter("-t 0 -c "+str(2.0**param_c))
            prob=svm_problem(trainingY,selectedtrainingX)
            m=svm_train(prob,param)

            selectedtestingX=genX(testingX,param_fset)
            p_label, p_acc, p_val=svm_predict(testingY,selectedtestingX,m)
            Y_list=Y_list+testingY
            predY_list=predY_list+p_label
            print("====================%d_%d(n=%d,c=%d):ACC=%f" % (runs,count,param_n,param_c,p_acc[0]))

            out_feature.write("%d_%d(n=%d,c=%d)\t" % (runs,count,param_n,param_c)),

            for fno in param_fset:
                out_feature.write(str(fno)+" "),
            out_feature.write("\t%f\n" % p_acc[0])

        acc,mse,scc=evaluations(Y_list,predY_list)
        acc_list.append(acc)
        out_acc.write("the %d run's accuracy:\t%f\n" % (runs,acc))
    out_acc.write("average ACC is %f" % average(acc_list))
    out_acc.close()
    out_feature.close()

if __name__ == "__main__":
    main(sys.argv[1:])
