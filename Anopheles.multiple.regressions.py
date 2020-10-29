
####Note that this script will be updated to a jupyter notebook format
import pandas as pd #loads the following packages to read data files, perform regressions, and make figures
from math import sqrt
from sklearn.linear_model import LinearRegression
from sklearn import metrics
from scipy.stats import linregress
import matplotlib.pyplot as plt
import numpy as np

final_null.to_csv("final.null.2L.GEA.frq.txt", index = None, sep = "\t") #loads neutral SNPs which serve as the null hypothesis
final_sig.to_csv("final.sig.2L.GEA.frq.txt", index = None, sep = "\t") # load significant SNPs from the RDA

null = final_null # change to new object to save original df 
sig = final_sig


null = null.set_index("Population/location") #setting and removing previous row/index names
null.index.names = [None]
sig = sig.set_index("Population/location")
sig.index.names = [None]


null_PreDry = null.drop(["N","Lats.","Lons.","srad","vapr","AnnTemp","Iso","TempSea","MaxTempWarm","MinTempCold","AnnPre","PreWet"],axis=1) #dropping the following columns for the multiple regression analysis. Note that this will change depending on the number and type of variable you are looking using. Here we lokk at regressions between allele frequencies and vprecipitation of the driest month
null_PreDry.head()
sig_PreDry = sig.drop(["N","Lats.","Lons.","srad","vapr","AnnTemp","Iso","TempSea","MaxTempWarm","MinTempCold","AnnPre","PreWet"],axis=1)
sig_PreDry.head()

null_PreDry = null_PreDry.astype(float) #change values to a float
sig_PreDry = sig_PreDry.astype(float)

null_var_SNPs = pd.DataFrame(columns=[0,1,2,3,4]) #makes an empty df with 5 columns (this will contain, "Slope","Intercept","Rvalue","Pvalue","Stderr")

for i in range(0,len(null_PreDry.columns)): #a loop that performs a series of linear regressions across all SNPs. (response = allele frequencies, predictor = precip. of the driest month) Depending on the number of SNPs, this may take some time.
    temp = (linregress(null_PreDry["PreDry"],null_PreDry.iloc[:,i]))
    temp = pd.DataFrame(temp)
    transposed = temp.transpose()
    null_var_SNPs = null_var_SNPs.append(transposed, ignore_index=True) #appends all iterations to a single df


null_var_SNPs.shape  #edits the final output and adds column names and r2 values
final_null_PreDry = null_var_SNPs
final_null_PreDry.columns = ["Slope","Intercept","Rvalue","Pvalue","Stderr"]
final_null_PreDry["R2"] = final_null_PreDry["Rvalue"]**2
final_null_PreDry = final_null_PreDry.drop(0)
final_null_PreDry.head()



sig_var_SNPs = pd.DataFrame(columns=[0,1,2,3,4])  #same loop as abov but  donw this significant SNPs

for i in range(0,len(sig_PreDry.columns)):
    temp = (linregress(sig_PreDry["PreDry"],sig_PreDry.iloc[:,i]))
    temp = pd.DataFrame(temp)
    transposed = temp.transpose()
    sig_var_SNPs = sig_var_SNPs.append(transposed, ignore_index=True)


sig_var_SNPs.shape
final_sig_PreDry = sig_var_SNPs
final_sig_PreDry.columns = ["Slope","Intercept","Rvalue","Pvalue","Stderr"]
final_sig_PreDry["R2"] = final_sig_PreDry["Rvalue"]**2
final_sig_PreDry = final_sig_PreDry.drop(0)
final_sig_PreDry.head()

final_null_PreDry.to_csv("final.null.3L.GEA.PreDry.txt", index = None, sep = "\t") #outputs final results for downstream analysis
final_sig_PreDry.to_csv("final.sig.3L.GEA.PreDry.txt", index = None, sep = "\t") #outputs final results for downstream analysis


