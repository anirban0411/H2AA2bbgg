import numpy as np
import pandas as pd
import xgboost as xgb
import awkward as ak
import pickle
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import uproot as uproot
import operator
import ROOT
from optparse import OptionParser
from xgboost import *
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, roc_curve, auc
from copy import copy
from datetime import datetime
from sklearn import metrics

sys.path.append('tools/')
from xgboost2tmva import *

#############################################################
def exportTree(final_dataset,treeNameToExport,export_filename): 
    export_file = uproot.recreate(export_filename)
    print("\n\t  file exported  : ",export_filename)
    allBranchesToWrite={ ky : final_dataset[ky].dtype for ky in final_dataset}
    export_file.mktree( treeNameToExport ,allBranchesToWrite)
    export_file[treeNameToExport].extend(final_dataset)
    export_file.close()

#############################################################

#########################################################
def getTH1FromNumpHist(nph):
    hist = ROOT.TH1D("","",len(nph[1])-1,np.array(nph[1],dtype=float))
    for i,freq in enumerate(nph[0]):
        hist.SetBinContent(i+1,freq)
    hist.ResetStats() 	
    return hist
#########################################################

def loadSignal(inputVars):

  vars = inputVars
  if 'evt_wgt' not in vars:
    vars.insert(0,'evt_wgt')
    
  dataset_Sig = uproot.open("/home/abala/local_bkp/analysis/BDT/full_runII/sig_nom/WH_AllmA_bdt_v1_weight_nom.root")['tree'].arrays(inputVars,library='pd')
  df_sig = dataset_Sig[vars]
  return df_sig

def loadMCBackground(inputVars, inputVarsMC):
  
  vars = inputVarsMC
  if 'evt_wgt' not in vars:
    vars.insert(0,'evt_wgt')
   
  #dataset_Bkg = uproot.open('2018_files/training_root_files_wgt1/ALLbkg_2018_bdt_v1.root')['tree'].pandas.df()[vars]
  dataset_Bkg = uproot.open('/home/abala/local_bkp/analysis/BDT/2018/Data/data_bdt_v1_weight_nom.root')['tree'].arrays(inputVars,library='pd')
  df_bkg = dataset_Bkg

  renames = {}
  keys = inputVarsMC
  values = inputVars
  for i,key in enumerate(keys):
        renames[key] = values[i]
  print('Renames:',renames)   
      
  df_bkg.rename(columns = renames, inplace = True) 
  
  return df_bkg  

def drawScores(predY_train,predY_test,W_train,W_test,nBins,c_min,c_max,name):

  train_Sig = predY_train[Y_train==1]
  train_Bkg = predY_train[Y_train==0]
  test_Sig = predY_test[Y_test==1]
  test_Bkg = predY_test[Y_test==0]
  
  W_trainSig_sum = W_train[Y_train==1].sum(0) 
  W_trainBkg_sum = W_train[Y_train==0].sum(0)
  W_testSig_sum = W_test[Y_test==1].sum(0)
  W_testBkg_sum = W_test[Y_test==0].sum(0)
  
  #Get histograms of the classifiers
  Histo_train_Sig = np.histogram(train_Sig,bins=int(nBins),range=(c_min,c_max),weights=W_train[Y_train==1])
  Histo_train_Bkg = np.histogram(train_Bkg,bins=int(nBins),range=(c_min,c_max),weights=W_train[Y_train==0])
  Histo_test_Sig = np.histogram(test_Sig,bins=int(nBins),range=(c_min,c_max),weights=W_test[Y_test==1]*(W_trainSig_sum/W_testSig_sum))
  Histo_test_Bkg = np.histogram(test_Bkg,bins=int(nBins),range=(c_min,c_max),weights=W_test[Y_test==0]*(W_trainBkg_sum/W_testBkg_sum))
 
  #########################################################

  fOut=ROOT.TFile(name+'_hists.root','RECREATE')
  fOut.cd()
  th1f_Histo_train_Sig = getTH1FromNumpHist(Histo_train_Sig) ;  th1f_Histo_train_Sig.SetName(name+"_train_sig_mvaHist") ; th1f_Histo_train_Sig.Write()
  th1f_Histo_train_Bkg = getTH1FromNumpHist(Histo_train_Bkg) ;  th1f_Histo_train_Bkg.SetName(name+"_train_bkg_mvaHist") ; th1f_Histo_train_Bkg.Write()  
  th1f_Histo_test_Sig  = getTH1FromNumpHist(Histo_test_Sig ) ;  th1f_Histo_test_Sig.SetName(name+"_test_sig_mvaHist")   ; th1f_Histo_test_Sig.Write() 
  th1f_Histo_test_Bkg  = getTH1FromNumpHist(Histo_test_Bkg ) ;  th1f_Histo_test_Bkg.SetName(name+"_test_bkg_mvaHist")   ; th1f_Histo_test_Bkg.Write()  
  fOut.Close()
  
  final_dataset={
            'BDT' : train_Bkg
          }
  filename=name+'_trees_bkg_train.root'
  treeNameToExport='Tree'
  exportTree(final_dataset,treeNameToExport,export_filename=filename)

  final_dataset={
            'BDT' : test_Bkg
          }
  filename=name+'_trees_bkg_test.root'
  treeNameToExport='Tree'
  exportTree(final_dataset,treeNameToExport,export_filename=filename)

  final_dataset={
            'BDT' : train_Sig
          }
  filename=name+'_trees_sig_train.root'
  treeNameToExport='Tree'
  exportTree(final_dataset,treeNameToExport,export_filename=filename)

  final_dataset={
            'BDT' : test_Sig
          }
  filename=name+'_trees_sig_test.root'
  treeNameToExport='Tree'
  exportTree(final_dataset,treeNameToExport,export_filename=filename)

  #########################################################
  #Lets get the min/max of the Histograms
  AllHistos= [Histo_train_Sig,Histo_train_Bkg,Histo_test_Sig,Histo_test_Bkg]
  h_max = max([histo[0].max() for histo in AllHistos])*1.2
  h_min = max([histo[0].min() for histo in AllHistos])
 
  #Get the histogram properties (binning, widths, centers)
  bin_edges = Histo_train_Sig[1]
  bin_centers = ( bin_edges[:-1] + bin_edges[1:]  ) /2.
  bin_widths = (bin_edges[1:] - bin_edges[:-1])
 
  #To make error bar plots for the data, take the Poisson uncertainty sqrt(N)
  ErrorBar_test_Sig = np.sqrt(Histo_test_Sig[0])
  ErrorBar_test_Bkg = np.sqrt(Histo_test_Bkg[0])  
  
  #Plot  
  ax1 = plt.subplot(111)
  ax1.bar(bin_centers,Histo_train_Sig[0],facecolor='blue',linewidth=0,width=bin_widths,label='Sig (Train)',alpha=0.5)
  ax1.bar(bin_centers,Histo_train_Bkg[0],facecolor='red',linewidth=0,width=bin_widths,label='Bkg (Train)',alpha=0.5)
  ax1.errorbar(bin_centers, Histo_test_Sig[0], yerr=ErrorBar_test_Sig, xerr=None, ecolor='blue',c='blue',fmt='o',label='Sig (Test)')
  ax1.errorbar(bin_centers, Histo_test_Bkg[0], yerr=ErrorBar_test_Bkg, xerr=None, ecolor='red',c='red',fmt='o',label='Bkg (Test)')
  ax1.axvspan(0.0, c_max, color='blue',alpha=0.08)
  ax1.axvspan(c_min,0.0, color='red',alpha=0.08)
  ax1.axis([c_min, c_max, h_min, h_max])
  plt.title("BDT Signal--Bkg separation")
  plt.xlabel("BDT output")
  plt.ylabel("frequency")
  legend = ax1.legend(loc='upper center', shadow=True,ncol=2)
  for alabel in legend.get_texts():
    alabel.set_fontsize('small')
  plt.savefig(name+".png",dpi=300)
  plt.savefig(name+".pdf",dpi=300)  


def drawROC(Y_train,Y_test,predY_train,predY_test,W_train,W_test,name):

  f,ax = plt.subplots(figsize=(8,8))

  fpr, tpr, threshold = metrics.roc_curve(Y_train,  predY_train, sample_weight=W_train)
  train_auc = metrics.auc(fpr, tpr)
#  ax.plot(fpr, tpr,label=f'Train:  AUC = {round(train_auc,3)}')
  ax.plot(tpr, fpr,label=f'Train:  AUC = {round(train_auc,3)}')

  fpr, tpr, threshold = metrics.roc_curve(Y_test,  predY_test, sample_weight=W_test)
  test_auc = metrics.auc(fpr, tpr)
#  ax.plot(fpr, tpr,label=f'Test:  AUC = {round(test_auc,3)}')
  ax.plot(tpr, fpr,label=f'Test:  AUC = {round(test_auc,3)}')

  plt.legend(loc='lower right')
  plt.title("ROC")
  plt.xlabel('Signal Efficiency')
  plt.ylabel('1 - Background Efficiency')
#  plt.ylabel('True Positive Rate')
#  plt.xlabel('False Positive Rate')
  plt.savefig(name+".png",dpi=300)
  plt.savefig(name+".pdf",dpi=300)  
  
  print('AUC for training = %1.3f'%(train_auc) )
  print('AUC for test     = %1.3f'%(test_auc) )  



if __name__ == '__main__':

  parser = OptionParser()
  parser.add_option("", "--useMC",help="Use MC for training", action="store_true", dest="useMC")
  parser.add_option("", "--loadModel" , help="Load a previously trained model. [ provide the name of the model ]"  , default=None)
  parser.add_option("", "--export",help="Export the model to the file [ provide the output file path ]", default=None)

  (options, args) = parser.parse_args() 
  useMC = options.useMC  
  loadModel = options.loadModel
  export = options.export
  
  inputVars =   ["leppt","lepeta","b1eta","b2eta","pT1_by_mbb","pT2_by_mbb","Njets","pho1eta","pho2eta","pT1_by_mgg","pT2_by_mgg","leadbdeepjet","subleadbdeepjet","leadgMVA","subleadgMVA","delphi_gg","delphi_bb","delphi_bbgg","delphi_ggMET"]
  inputVarsMC = ["leppt","lepeta","b1eta","b2eta","pT1_by_mbb","pT2_by_mbb","Njets","pho1eta","pho2eta","pT1_by_mgg","pT2_by_mgg","leadbdeepjet","subleadbdeepjet","leadgMVA","subleadgMVA","delphi_gg","delphi_bb","delphi_bbgg","delphi_ggMET"]
  
  print("Loading Signal...")
  df_sig = loadSignal(inputVars)

  print("Loading MC Bkgs...")
  if useMC:
    df_bkg = loadMCBackground(inputVars,inputVarsMC)

  #Add target for signal (1)
  df_sig['target'] = np.ones(len(df_sig.index))
  df_sig['target'] = df_sig['target'].astype('int')

  #Add target for backgrounds (0)
  df_bkg['target'] = np.zeros(len(df_bkg.index))
  df_bkg['target'] = df_bkg['target'].astype('int')


  #Deal with negative weights (set to 1 and rescale accordingly)
  sig_sum = df_sig.sum(0)['evt_wgt']
  #df_sig['weight'] = df_sig['weight'].abs()
  df_sig['evt_wgt'] = df_sig['evt_wgt'].apply(lambda x: x if x>=0. else 1.)
  sig_sum_abs = df_sig.sum(0)['evt_wgt']
  df_sig['evt_wgt'] = df_sig['evt_wgt'].apply(lambda x: x*(sig_sum/sig_sum_abs))
  sig_sum = df_sig.sum(0)['evt_wgt']
  print(sig_sum)

  #Deal with negative weights (set to 1 and rescale accordingly)
  bkg_sum = df_bkg.sum(0)['evt_wgt']
  #df_bkg['weight'] = df_bkg['weight'].abs()
  df_bkg['evt_wgt'] = df_bkg['evt_wgt'].apply(lambda x: x if x>=0. else 1.)
  bkg_sum_abs = df_bkg.sum(0)['evt_wgt']
  df_bkg['evt_wgt'] = df_bkg['evt_wgt'].apply(lambda x: x*(bkg_sum/bkg_sum_abs))
  bkg_sum = df_bkg.sum(0)['evt_wgt']
#  df_bkg['evt_wgt'] = df_bkg['evt_wgt'].apply(lambda x: x*(sig_sum/bkg_sum))
  df_sig['evt_wgt'] = df_sig['evt_wgt'].apply(lambda x: x*(bkg_sum/sig_sum))
  print(bkg_sum)

  print("Loading done!...")
  
  if 'evt_wgt' in inputVars:
    inputVars.remove('evt_wgt')
  if 'Norm_SFs' in inputVars:
    inputVars.remove('Norm_SFs')
         
  seed = int(datetime.now().timestamp())
  seed = 12345
  np.random.seed(seed)
  print("Seed:",seed)

  X = pd.concat((df_sig.iloc[:,0:20] ,df_bkg.iloc[:,0:20] ))
  Y = pd.concat((df_sig['target'] ,df_bkg['target'] ))

  test_size = 0.0000001
  X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, shuffle=0)
  
  W_train = X_train['evt_wgt'].values
  X_train = X_train.drop('evt_wgt',axis=1).values
  Y_train = Y_train.values
  W_test  = X_test['evt_wgt'].values
  X_test = X_test.drop('evt_wgt',axis=1).values
  Y_test = Y_test.values

  params = {
    'booster': 'gbtree',
    'eval_metric': 'auc',
    'learning_rate': 0.5,
    'max_depth': 3,
    'gamma': 0.1,
    'objective': 'binary:logistic',
  }
  
  training = xgb.DMatrix(X_train,label=Y_train, weight=W_train,feature_names=inputVars)
  testing  = xgb.DMatrix(X_test, label=Y_test , weight=W_test, feature_names=inputVars) 
  

  if loadModel:
    with open(loadModel,'rb') as f:
        modelDict=pickle.load(f)
        model=modelDict['xgbModel']
        modelName=modelDict['modelName']
        print("Loaded model ",modelName)
  else:
    print("Training model! ")
    model = xgb.train(params,training) 
  if export:
    print("Model exported to : ",export)
    modelDict={
        'xgbModel' : model,
        'modelName' : 'HtoAATo4b2gVsTTBar',
    }
    with open(export,'wb') as f:
        pickle.dump(modelDict,f)
        

  predY_train = model.predict(training)
  predY_test  = model.predict(testing)
  
  label = '_UL18'
  if useMC:
    label = label + '_MC' 
  else:
    label = label + '_dataDriven' 
          
  drawScores(predY_train,predY_test,W_train,W_test,40,0.,1.,'BDT'+label)
  drawROC(Y_train,Y_test,predY_train,predY_test,W_train,W_test,'ROC'+label)
          
  mdl = model.get_dump()
  input_vars=[]
  for i,key in enumerate(inputVars):
    input_vars.append((key,'F'))
  convert_model(mdl,input_variables=input_vars,output_xml='H2AA2bbggMVA_xgboost_model'+label+'.xml')  
