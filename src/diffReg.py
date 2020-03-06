# -*- coding: utf-8 -*-
import click
import os
import prepare_input
import sQTLs
import prepare_diffReg
import sig_sQTLs
import diffReg_core

class InputDirectory():
	def __init__ (self,work_dir,bed_file_name,iso_file_name,snp_file_name,sample_file_name):
		if work_dir[-1]!='/':
			work_dir+='/'
		self.WORK_DIR=work_dir
		self.BED_FILE_NAME=bed_file_name
		self.TRANSCRIPT_FILE_NAME=iso_file_name
		self.SNP_FILE_NAME=snp_file_name
		self.SAMPLE_LABEL_DICT_FILE_NAME=sample_file_name
	def showDirectory(self):
		cwd=self.WORK_DIR
		print('work_directory:',cwd)
		print('bed_file_name:',cwd+self.BED_FILE_NAME)
		print('isoform_file_name:',cwd+self.TRANSCRIPT_FILE_NAME)
		print('snp_file_name:',cwd+self.SNP_FILE_NAME)
		print('sample_label_dict_file_name:',cwd+self.SAMPLE_LABEL_DICT_FILE_NAME)

@click.command()
@click.option('-w','--work-dir',required=True,type=str,help='all data for diffReg are stored in this directory.')
@click.option('-b','--bed-file-name',default='bed.tsv',show_default=True,type=str,help='gene bed information file\'s name(stored in work-dir)')
@click.option('-i','--iso-file-name',default='transcript.tsv',show_default=True,type=str,help='transcript isoform expression file\'s name(stored in work-dir)')
@click.option('-p','--snp-file-name',default='snp.tsv',show_default=True,type=str,help='SNP file\'s name(stored in work-dir)')
@click.option('-s','--sample-file-name',default='samples.txt',show_default=True,type=str,help='samples information file\'s name(stored in work-dir)')
@click.option('-t','--threads',default=20,show_default=True,help='numbers of process to run DiffReg')
@click.option('-n','--n-permutation',default=1000000,show_default=True,help='times of permutation test')
def pass_param(work_dir,bed_file_name,iso_file_name,snp_file_name,sample_file_name,threads,n_permutation):
	diffReg_dir = os.path.abspath(os.path.dirname(__file__))
	print('script_directory:',diffReg_dir)	
	indir=InputDirectory(work_dir,bed_file_name,iso_file_name,snp_file_name,sample_file_name)
	#indir.showDirectory()
	prepare_input.split_AD_Control(indir.WORK_DIR,
								indir.TRANSCRIPT_FILE_NAME,
								indir.SNP_FILE_NAME,
								indir.SAMPLE_LABEL_DICT_FILE_NAME)
	# temporary data are generate in directory ./tmp 
	sQTLs.multiprocessing_AD_Control_sQTL(indir.WORK_DIR,diffReg_dir,indir.BED_FILE_NAME)
	sig_sQTLs.get_significant_sQTLs(indir.WORK_DIR + 'tmp','AD_sQTL_result.tsv','Control_sQTL_result.tsv')
	prepare_diffReg.prepareDiffReg(indir.WORK_DIR+'tmp')
	diffReg_core.diffReg(indir.WORK_DIR,n_permutation,threads)
if __name__ == '__main__' :
	pass_param()