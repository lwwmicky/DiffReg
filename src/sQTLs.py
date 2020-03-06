import os
import multiprocessing
from time import sleep
def Rsript_cmd(cwd,diffReg_dir,group,bed_f_name):
	if cwd[-1]!='/':
		cwd+='/'
	tr=cwd+'tmp/'+group+'_transcript.tsv'
	snp=cwd+'tmp/'+group+'_snp.tsv'
	out=cwd+'tmp/'+group+'_sQTL_result.tsv'
	bed=cwd+bed_f_name
	return run_sQTL_cmd(diffReg_dir,snp,tr,bed,out)

def run_sQTL_cmd(diffReg_dir,snp_f,tr_f,bed_f,out_f):
	cmd='Rscript '+diffReg_dir+'run_sQTLseekeR.R '+tr_f+' '+snp_f+' '+bed_f+' '+out_f
	return cmd

def run_sQTL(snp_f,tr_f,beb_f,out_f):
	os.system(run_sQTL_cmd(snp_f,tr_f,beb_f,out_f))

def multiprocessing_AD_Control_sQTL(cwd,diffReg_dir,bed_f_name):
	if diffReg_dir[-1]!='/':
		diffReg_dir+='/'
	pool = multiprocessing.Pool(processes = 2)
	ad_cmd=Rsript_cmd(cwd,diffReg_dir,'AD',bed_f_name)
	ctrl_cmd=Rsript_cmd(cwd,diffReg_dir,'Control',bed_f_name)
	pool.apply_async(os.system, (ad_cmd,))
	print('[info]Filtering non-sQTLs in disease group start...')
	pool.apply_async(os.system, (ctrl_cmd,))
	print('[info]Filtering non-sQTLs in control group start...')
	pool.close()
	pool.join()