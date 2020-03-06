# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 16:01:45 2019

@author: liuws
"""
import pandas as pd
import numpy as np

def devideSum(row,colsum):
	return pd.Series(np.round(((row.values)/colsum),4),index=row.index)

def tr_normalized(m):
	colsum=m.apply(sum)
	#如果这个人的所有transcript表达都为零，为避免除以零，让所有值为1
	colsum[colsum==0]=1 
	rs=m.apply(devideSum,colsum=colsum.values,axis=1)
	return rs

def get_normalized_Trdata(caseTr):
	caseTr.index=caseTr.trId.values
	caseTrInfo=caseTr.iloc[:,:2]
	caseTrData=caseTr.iloc[:,2:]
	rs=pd.DataFrame(columns=caseTrData.columns)
	for gene in list(set(caseTr.geneId.values)):
		df=caseTrData[caseTr.geneId==gene]
		df_normalized=tr_normalized(df)
		rs=rs.append(df_normalized)
	joinDf=caseTrInfo.join(other=rs)
	del(joinDf['trId'])
	del(joinDf['geneId'])
	return joinDf

def load_used_SNP_matrix(snp_file,snps):
	snp_df=pd.read_csv(snp_file,sep='\t')
	snp_df=snp_df[snp_df.snpId.isin(snps)]
	snp_df.index=snp_df.snpId.values
	new_cols=snp_df.columns.values[4:]
	return snp_df.reindex(columns=new_cols)


def prepareDiffReg(cwd,pairs_f_name='pairs.tsv',
	AD_transcript_f_name='AD_transcript.tsv',
	Control_transcript_f_name='Control_transcript.tsv',
	AD_snp_f_name='AD_snp.tsv',
	Control_snp_f_name='Control_snp.tsv'):

	if cwd[-1]!='/':
		cwd+='/'
	pairs_df=pd.read_csv(cwd+pairs_f_name,sep='\t')
	ad_tr=pd.read_csv(cwd+AD_transcript_f_name,sep='\t')
	ctrl_tr=pd.read_csv(cwd+Control_transcript_f_name,sep='\t')

	ad_tr.reindex(columns=['trId','geneId']).to_csv(cwd+'gene_tr.tsv',index=None)

	ad_tr_sub=ad_tr[ad_tr['geneId'].isin(pairs_df['gene'])]
	ctrl_tr_sub=ctrl_tr[ctrl_tr['geneId'].isin(pairs_df['gene'])]

	ad_tr_normalized=get_normalized_Trdata(ad_tr_sub)
	ctrl_tr_normalized=get_normalized_Trdata(ctrl_tr_sub)
	
	ad_tr_normalized.to_csv(cwd+'AD_transcript_normalized.tsv')
	ctrl_tr_normalized.to_csv(cwd+'Control_transcript_normalized.tsv')
	
	ad_snp_df=load_used_SNP_matrix(cwd+AD_snp_f_name,pairs_df['SNP'])
	ctrl_snp_df=load_used_SNP_matrix(cwd+Control_snp_f_name,pairs_df['SNP'])
	
	ad_snp_df.to_csv(cwd+'AD_snp_DiffReg_format.tsv')
	ctrl_snp_df.to_csv(cwd+'Control_snp_DiffReg_format.tsv')