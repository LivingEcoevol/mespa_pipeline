#!/usr/bin/python
#-*- coding: UTF-8 -*-
import os,sys
import glob
from Bio import SeqIO

#Step 1: append amino acid sequences
#rename headers for the file '*_aa_scaffolds.fas'
def rename_header_aa(path, taxon):
	out=open(path+'/'+taxon+'_aa.fas','w')
	for record in SeqIO.parse(path+'/'+taxon+'_aa_scaffolds.fas','fasta'):
		record.id=record.id.split(';')[1]
		record.description=''
		SeqIO.write(record,out,'fasta')
	out.close()
	return

for line in open('/scratch1/li266/AHE/mespa/sample_names.txt'):
	taxon=line.strip()
	rename_header_aa('/scratch1/li266/AHE/mespa/mespa_results', taxon)
print 'rename_done'


#pool all sequences into folder /scratch1/li266/AHE/mespa/aa_taxa, which contains outgroup sequences
for line in open('/scratch1/li266/AHE/mespa/sample_names.txt'):
	taxon=line.strip()
	os.rename('/scratch1/li266/AHE/mespa/mespa_results/'+taxon+'_aa.fas','/scratch1/li266/AHE/mespa/aa_taxa/'+taxon+'.aa.fas')
print 'pool_done'


#sort by gene
path_in='/scratch1/li266/AHE/mespa/aa_taxa'
path_out='/scratch1/li266/AHE/mespa/aa_loci'
if os.path.exists(path_out):
	pass
else:
	os.makedirs(path_out)

gene_names=[]
filenames=glob.glob(path_in+'/'+'*.aa.fas')
for filename in filenames:
	for record in SeqIO.parse(filename,'fasta'):
		gene=record.id
		if gene not in gene_names:
			gene_names.append(gene)
print len(gene_names)

for gene in gene_names:
	outname=gene+'.aa.fas'
	out=open(path_out+'/'+outname,'w')
	for filename in filenames:
		taxon=filename.split('/')[6].split('.')[0]
		for r in SeqIO.parse(filename,'fasta'):
			if gene == r.id:
				r.id=taxon
				r.description=''
				SeqIO.write(r,out,'fasta')
	out.close()
print 'sort_by_gene_done'


#Step 2: append nucleotide sequences
def rename_header_nt(path, taxon):
	dic={}
	for line in open(path+'/'+taxon+'_order_scaff.txt'):
		scaff=line.split('\t')[0].strip()
		gene=line.split('\t')[1].strip()
		dic[scaff]=gene
	print dic

	out=open(path+'/'+taxon+'_nt.fas','w')
	for record in SeqIO.parse(path+'/'+taxon+'_onlygenemodels.fa','fasta'):
		if record.id in dic:
			record.id=dic[record.id]
			record.description=''
			SeqIO.write(record,out,'fasta')
		else:
			SeqIO.write(record,out,'fasta')
			print(record.id)
	out.close()
	return

for line in open('/scratch1/li266/AHE/mespa/sample_names.txt'):
	taxon=line.strip()
	rename_header_nt('/scratch1/li266/AHE/mespa/mespa_results', taxon)
print 'rename_done'


#pool all sequences into folder /scratch1/li266/AHE/mespa/nt_taxa, which contains outgroup sequences
for line in open('/scratch1/li266/AHE/mespa/sample_names.txt'):
	taxon=line.strip()
	os.rename('/scratch1/li266/AHE/mespa/mespa_results/'+taxon+'_nt.fas','/scratch1/li266/AHE/mespa/nt_taxa/'+taxon+'.nt.fas')
print 'pool_done'

#sort by gene
path_in='/scratch1/li266/AHE/mespa/nt_taxa'
path_out='/scratch1/li266/AHE/mespa/nt_loci'
if os.path.exists(path_out):
	pass
else:
	os.makedirs(path_out)

gene_names=[]
filenames=glob.glob(path_in+'/'+'*.nt.fas')
for filename in filenames:
	for record in SeqIO.parse(filename,'fasta'):
		gene=record.id
		if gene not in gene_names:
			gene_names.append(gene)
print len(gene_names)

for gene in gene_names:
	outname=gene+'.nt.fas'
	out=open(path_out+'/'+outname,'w')
	for filename in filenames:
		taxon=filename.split('/')[6].split('.')[0]
		for r in SeqIO.parse(filename,'fasta'):
			if gene == r.id:
				r.id=taxon
				r.description=''
				SeqIO.write(r,out,'fasta')
	out.close()
print 'sort_by_gene_done'