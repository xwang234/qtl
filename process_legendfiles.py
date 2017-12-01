#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 10:02:47 2017
Add filters of EUR and EAS/SAS to legend files
@author: xwang234
"""
import pandas as pd

infolder=outfolder="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3"

def makelegend(legendfile,outfile):
    data=pd.read_csv(legendfile,compression='gzip',sep=" ")
    data["filter"]=0
    idx=((data["EUR"]>0.001) & (data["EUR"]<0.999)) | ((data["SAS"]>0.001) & (data["SAS"]<0.999)) | ((data["EAS"]>0.001) & (data["EAS"]<0.999))
    data.loc[idx,"filter"]=1
    data.to_csv(outfile,compression='gzip',sep=" ",index=False)
 
#only consider EUR and EAS, used in iCOGS and ONCOarray
def makelegend1(legendfile,outfile):
    data=pd.read_csv(legendfile,compression='gzip',sep=" ")
    data["filter"]=0
    idx=((data["EUR"]>0.001) & (data["EUR"]<0.999)) | ((data["EAS"]>0.001) & (data["EAS"]<0.999))
    data.loc[idx,"filter"]=1
    data.to_csv(outfile,compression='gzip',sep=" ",index=False)
    
for i in range(1,23):
    print i
    legendfile=infolder+"/1000GP_Phase3_chr"+str(i)+".legend.gz"
    outfile=infolder+"/1000GP_Phase3_chr"+str(i)+".legend_filter.gz"
    makelegend(legendfile,outfile)

#for chrX    
legendfile=infolder+"/1000GP_Phase3_chrX_PAR1.legend.gz"
outfile=infolder+"/1000GP_Phase3_chrX_PAR1.legend_filter.gz"
makelegend(legendfile,outfile)

legendfile=infolder+"/1000GP_Phase3_chrX_PAR2.legend.gz"
outfile=infolder+"/1000GP_Phase3_chrX_PAR2.legend_filter.gz"
makelegend(legendfile,outfile)

legendfile=infolder+"/1000GP_Phase3_chrX_NONPAR.legend.gz"
outfile=infolder+"/1000GP_Phase3_chrX_NONPAR.legend_filter.gz"
makelegend(legendfile,outfile)

#makelegend1
for i in range(1,23):
    print i
    legendfile=infolder+"/1000GP_Phase3_chr"+str(i)+".legend.gz"
    outfile=infolder+"/1000GP_Phase3_chr"+str(i)+".legend_filter1.gz"
    makelegend1(legendfile,outfile)

#for chrX    
legendfile=infolder+"/1000GP_Phase3_chrX_PAR1.legend.gz"
outfile=infolder+"/1000GP_Phase3_chrX_PAR1.legend_filter1.gz"
makelegend1(legendfile,outfile)

legendfile=infolder+"/1000GP_Phase3_chrX_PAR2.legend.gz"
outfile=infolder+"/1000GP_Phase3_chrX_PAR2.legend_filter1.gz"
makelegend1(legendfile,outfile)

legendfile=infolder+"/1000GP_Phase3_chrX_NONPAR.legend.gz"
outfile=infolder+"/1000GP_Phase3_chrX_NONPAR.legend_filter1.gz"
makelegend1(legendfile,outfile)