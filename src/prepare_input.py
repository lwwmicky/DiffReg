import os
import pandas as pd
import shutil
def split_AD_Control_kerl(cwd,file,fileType,idx_first_sample,labels):
	df=pd.read_csv(file,sep='\t')
	info_df=df[df.columns[:idx_first_sample]]
	ad_df=make_group_dataframe(df,info_df,labels,1)
	ctrl_df=make_group_dataframe(df,info_df,labels,0)
	ad_df.to_csv(cwd+'tmp/AD_'+fileType+'.tsv',sep='\t',index=None)
	ctrl_df.to_csv(cwd+'tmp/Control_'+fileType+'.tsv',sep='\t',index=None)

def make_group_dataframe(df,info_df,labels,group):
	sid=labels.index[labels==group]
	sid_df=df.reindex(columns=sid)
	return pd.concat([info_df,sid_df],axis=1)

def load_samples_label_dict(file_label_dict):
	with open(file_label_dict,'r')as f:
		labs=eval(f.read())
	return labs

def split_AD_Control(cwd,tr_f_name,snp_f_name,sample_f_name):
	if cwd[-1]!='/':
		cwd=cwd+'/'
	labels=pd.Series(load_samples_label_dict(cwd+sample_f_name))
	print('[info] make directory \''+cwd+'tmp\' for temporary output during all the program')
	if(os.path.exists(cwd+'tmp/')):
		ch=input("[warning] "+cwd+"tmp directory has been made,files in \'"+cwd+"tmp\' will be deleted(yes/no)")
		if(ch in ['yes','y','Y']):
			shutil.rmtree(cwd+'tmp/')
			os.mkdir(cwd+'tmp/')
			print('[info] files in \'tmp\' removed.')
		else:
			print('[info] files in \'tmp\' reserved.')
	else:
		os.mkdir(cwd+'tmp/')
	split_AD_Control_kerl(cwd,cwd+tr_f_name,'transcript',2,labels)
	split_AD_Control_kerl(cwd,cwd+snp_f_name,'snp',4,labels)

