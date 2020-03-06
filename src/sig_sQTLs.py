import pandas as pd

def get_significant_sQTLs(cwd,rs_ad_f_name,rs_ctrl_f_name):
	if cwd[-1]!='/':
		cwd+='/'
	ad_pairs=significant_sqtls(cwd+rs_ad_f_name)
	ctrl_pairs=significant_sqtls(cwd+rs_ctrl_f_name)
	union_pairs=list(set(ad_pairs)|set(ctrl_pairs))
	print('[significant sQTLs summary]')
	show_sqtl_result(ad_pairs,'Disease')
	show_sqtl_result(ctrl_pairs,'Control')
	show_sqtl_result(union_pairs,'union Disease and Control')
	ser_pairs=pd.Series(union_pairs)
	genes=ser_pairs.str.split('/').str.get(0)
	snps=ser_pairs.str.split('/').str.get(1)
	pairs_df=pd.DataFrame([genes,snps],index=['gene','SNP']).T
	pairs_df.to_csv(cwd+'pairs.tsv',sep='\t',index=None)

def significant_sqtls(sqtl_result_file):
	sqtl_rs=pd.read_csv(sqtl_result_file,sep='\t')
	sig_keep=(sqtl_rs['qv']<=0.05)
	sqtl_rs['snp-gene']=sqtl_rs['geneId']+'/'+sqtl_rs['snpId']
	return list(sqtl_rs['snp-gene'][sig_keep].values)

def show_sqtl_result(pairs,grp):
	if len(pairs)==0:
		print('\t{0} : {1} significant sQTLs({2} genes,{3} SNPs).'.format(grp,0,0,0))
	else:
		ser_pairs=pd.Series(pairs)
		sig_genes=ser_pairs.str.split('/').str.get(0)
		num_genes=len(set(list(sig_genes)))
		sig_snps=ser_pairs.str.split('/').str.get(1)
		num_snps=len(set(list(sig_snps)))
		print('\t{0} : {1} significant sQTLs({2} genes,{3} SNPs).'.format(grp,len(pairs),num_genes,num_snps))
