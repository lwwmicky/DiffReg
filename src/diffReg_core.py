# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import multiprocessing
import os
def myAnd(x,y):
	if x>0 and y>0:
		return 1
	return 0
def diffReg_kernel(tr_case,snp_case,tr_control,snp_control,pairs,gene2tr,numOfpermutation):
	rs=pd.DataFrame(columns=['delta','pv'])
	for gene,snp in zip(pairs.gene.values,pairs.SNP.values):
		trs=gene2tr.trId[gene2tr.geneId==gene]
		#下面两处fillna(0)是因为有一些transcript在case和control中都是零，都在去除协变量方差中被过滤了，这里补回来
		tr_ingene_case=(tr_case.reindex(index=trs).fillna(0)).values
		tr_ingene_control=(tr_control.reindex(index=trs).fillna(0)).values
		case_snpfix=snp_case.loc[snp]
		control_snpfix=snp_control.loc[snp]
		deltaIndex=snp+'/'+trs.values
		deltacnt=np.zeros_like(trs.values)
		normal=True
		caselen0=sum(case_snpfix==0)
		caselen1=sum(case_snpfix==1)
		caselen2=sum(case_snpfix==2)
		controllen0=sum(control_snpfix==0)
		controllen1=sum(control_snpfix==1)
		controllen2=sum(control_snpfix==2)
		if myAnd(caselen0,controllen0)+myAnd(caselen1,controllen1)+myAnd(caselen2,controllen2)<2:
			continue
		for i in range(numOfpermutation+1):
			case_snpx=case_snpfix.copy().values
			control_snpx=control_snpfix.copy().values
			if(i!=0):
				np.random.shuffle(case_snpx)
				np.random.shuffle(control_snpx)	
			casem0=tr_ingene_case[:,case_snpx==0]
			casem1=tr_ingene_case[:,case_snpx==1]
			casem2=tr_ingene_case[:,case_snpx==2]
			controlm0=tr_ingene_control[:,control_snpx==0]
			controlm1=tr_ingene_control[:,control_snpx==1]
			controlm2=tr_ingene_control[:,control_snpx==2]
			casemean0=casemean1=casemean2=controlmean0=controlmean1=controlmean2=0
			if caselen0>0:
				casemean0=np.mean(casem0,axis=1)
			if caselen1>0:			
				casemean1=np.mean(casem1,axis=1)
			if caselen2>0:				
				casemean2=np.mean(casem2,axis=1)
			if controllen0>0:				
				controlmean0=np.mean(controlm0,axis=1)		 
			if controllen1>0:				
				controlmean1=np.mean(controlm1,axis=1)	  
			if controllen2>0:
				controlmean2=np.mean(controlm2,axis=1)
			delta01=delta12=delta02=0
			flag=0
			if caselen0!=0 and caselen1!=0 and controllen0!=0 and controllen1!=0:
				flag=1
				casek01=casemean0-casemean1
				controlk01=controlmean0-controlmean1
				delta01=abs(casek01-controlk01)
			if controllen1!=0 and controllen2!=0 and caselen1!=0 and caselen2!=0:
				flag=1
				casek12=casemean1-casemean2
				controlk12=controlmean1-controlmean2
				delta12=abs(casek12-controlk12)
			if caselen0!=0 and caselen2!=0 and controllen0!=0 and controllen2!=0:
				flag=1
				casek02=casemean0-casemean2
				controlk02=controlmean0-controlmean2
				delta02=abs(casek02-controlk02)  
			
			if flag==1:
				if(i==0):
					delta=delta01+delta12+delta02
				else:
					deltap=delta01+delta12+delta02
					deltacnt+=(deltap>delta).astype(int)
			else:
				normal=False
				break
		deltacnt=deltacnt+1
		deltacnt=deltacnt/numOfpermutation
		if(normal):
			rs=rs.append(pd.DataFrame(np.array([delta,deltacnt]).T,index=deltaIndex,columns=['delta','pv']))
	rs['delta']=rs['delta'].astype(float)
	return rs
	
def seprate_snp(snp_df,n):
	df_len=len(snp_df)
	sub_df_len=int(df_len/n)
	if(df_len%n!=0):
		n=n+1
	sub_dfs=[]
	for i in range(n):
		i_start=i*sub_df_len
		i_end=i_start+sub_df_len
		sub_dfs.append(snp_df.iloc[i_start:min(i_end,df_len)])
	return sub_dfs

def multiprocess_diffReg(ad_tr,ad_snp,ctrl_tr,
						ctrl_snp,pairs,gene2tr,
						timesOfPermutation,numOfParallel):
	real_processes=numOfParallel
	if len(ad_snp)%numOfParallel !=0:
		real_processes=numOfParallel+1
	ad_snps=seprate_snp(ad_snp,numOfParallel)
	ctrl_snps=seprate_snp(ctrl_snp,numOfParallel)
	pool=multiprocessing.Pool(processes=real_processes)
	rs=list(range(real_processes))
	idx_sub=0
	for (ad_sub_snp,ctrl_sub_snp) in zip(ad_snps,ctrl_snps):
		sub_pair=pairs[pairs.SNP.isin(ad_sub_snp.index.values)]
		rs[idx_sub]=pool.apply_async(diffReg_kernel,(ad_tr,ad_sub_snp,ctrl_tr,ctrl_sub_snp,sub_pair,gene2tr,timesOfPermutation))
		idx_sub+=1
	pool.close()
	pool.join()
	merge_rs=pd.DataFrame(columns=rs[0].get().columns)
	for sub_rs in rs:
		merge_rs=merge_rs.append(sub_rs.get())
	rs_order=merge_rs.sort_values(by='delta',ascending=False)
	rs_order['delta']=np.around(rs_order['delta'],4)
	return rs_order

def result_format(result_df,gene2tr):
	ids=pd.Series(result_df.index.values)
	result_df['snp']=ids.str.split('/').str.get(0).values
	result_df['transcript']=ids.str.split('/').str.get(1).values
	map_id=pd.Series(gene2tr.geneId.values,index=gene2tr.trId.values)
	result_df['gene']=result_df['transcript'].map(map_id)
	return result_df.reindex(columns=['snp','transcript','gene','delta','pv'])

def pv_adjust(p):
	"""Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
	p = np.asfarray(p)
	by_descend = p.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(p)) / np.arange(len(p), 0, -1)
	q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
	return q[by_orig]

def diffReg(cwd,n_runs,n_threads):
	print('[info] diffReg start...')
	if cwd[-1]!='/':
		cwd+='/'
	rs_wd=cwd+'DiffRegResult/'
	cwd+='tmp/'
	if not os.path.exists(rs_wd):
		print('[info] make directory '+rs_wd+' for DiffReg result.')
		os.mkdir(rs_wd)
	ad_tr=pd.read_csv(cwd+'AD_transcript_normalized.tsv',index_col=0)
	ctrl_tr=pd.read_csv(cwd+'Control_transcript_normalized.tsv',index_col=0)
	ad_snp=pd.read_csv(cwd+'AD_snp_DiffReg_format.tsv',index_col=0)
	ctrl_snp=pd.read_csv(cwd+'Control_snp_DiffReg_format.tsv',index_col=0)
	pairs=pd.read_csv(cwd+'pairs.tsv',sep='\t')
	gene2tr=pd.read_csv(cwd+'gene_tr.tsv')
	result_df=multiprocess_diffReg(ad_tr,ad_snp,ctrl_tr,ctrl_snp,pairs,gene2tr,n_runs,n_threads)
	result_df=result_format(result_df,gene2tr)
	result_df['p_adj']=pv_adjust(result_df['pv'])
	result_df.to_csv(rs_wd+'diffReg_result.tsv',sep='\t',index=None)
	sig_df=result_df[(result_df['delta']>0)&(result_df['p_adj']<0.05)]
	sig_df.drop_duplicates(subset='gene',keep='first').to_csv(rs_wd+'significant_predict_gene.tsv',sep='\t',index=None)
	print('[Done] DiffReg result stored in '+rs_wd)