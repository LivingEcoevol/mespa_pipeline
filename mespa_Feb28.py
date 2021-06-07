#!/usr/bin/python
from Bio import SeqIO
from datetime import datetime
import argparse,os,subprocess,re,operator,time,shutil
from Bio.Blast.Applications import NcbiblastpCommandline
from os.path import basename
from collections import OrderedDict
def parse_config(filename):
    options = {}
    comment_char = '#' 
    option_char =  '=' 
    for line in filename:
        if comment_char in line:
            line, comment = line.split(comment_char, 1)
        if option_char in line:
            option, value = line.split(option_char, 1)
            option = option.strip()
            value = value.strip()
            options[option] = value
    filename.close()
    if args.choice == 1:
        try:
            with open(options['db']) as test:
	        pass
        except IOError:
            print "BLAST database not found"
            raise SystemExit
    if options['evalue'] == "":
        options['evalue'] = "1e-5"
    if options['numthreads'] == "":
        options['numthreads'] = "2"
    if options['spacer']== "":
	options['spacer']=100
#for duplication consider only alignments with 80 coverage and minimum protein length of 150
    if options['dup_protein_length'] =="":
	options['dup_protein_length'] = 300
    if options['dup_coverage'] == "":
	options['dup_coverage'] = 80
    if options['overlap_dist'] == "":
        options['overlap_dist'] = 25
    if options['disjoint_dist'] == "":
        options['disjoint_dist'] = 25
    return options

def spaln(protein,genome,iteration):
    spaln_pids = []
    timestr = ""
    spaln_calls = []
#    print genome
    base = basename(os.path.abspath(genome).split(os.extsep)[0])
    os.environ["ALN_DBS"]=options['path_to_spaln_seqdb']
    os.environ["PATH"]=os.environ["PATH"]+":"+options['path_to_spaln']
    os.environ["ALN_TAB"]=options['path_to_spaln_table']
#    print os.environ["ALN_DBS"],os.environ["ALN_TAB"],base
    time=datetime.now()
    time_tuple = time.timetuple()
    for i in time_tuple:
        timestr = timestr + str(i)
    link=options['path_to_spaln_seqdb']+"/"+base+"_"+timestr+".mfa"
#    print link
    os.symlink(os.getcwd()+"/"+genome, link)
    wd = os.getcwd()
    null = open(os.devnull,'w')
    os.chdir(options['path_to_spaln_seqdb'])
    index_call = "make "+base+"_"+timestr+".idx"
    subprocess.call(index_call, shell =True, stdout = null, stderr = null)
    index_call = "make "+base+"_"+timestr+".bkp"
    subprocess.call(index_call, shell =True, stdout = null, stderr = null)
    os.chdir(wd)
    output_basename = basename(os.path.abspath(genome).split(os.extsep)[0])
    spaln_calls.append(str(options['path_to_spaln']+"/spaln -Q7 -LS -O1 -S3 -M1 -d "+base+"_"+timestr+" -t "+options['numthreads']+" "+protein +" > spaln_"+ output_basename + ".aln"))
    if output_basename == "scaffolds":
        spaln_calls.append(str(options['path_to_spaln']+"/spaln -Q7 -LS -O0 -S3 -M1 -d "+base+"_"+timestr+" -t "+options['numthreads']+" "+protein +" > "+ output_basename + ".gff"))
    else:
        spaln_calls.append(str(options['path_to_spaln']+"/spaln -Q7 -LS -O0 -S3 -M1 -d "+base+"_"+timestr+" -t "+options['numthreads']+" "+protein +" > spaln_"+ output_basename + ".gff"))
#    spaln_calls.append(str(options['path_to_spaln']+"/spaln -Q7 -LS -O6 -S3 -M1 -d "+base+"_"+timestr+" -t "+options['numthreads']+" "+protein +" > spaln_"+ output_basename + ".fa"))
    spaln_calls.append(str(options['path_to_spaln']+"/spaln -Q7 -LS -O4 -S3 -M1 -d "+base+"_"+timestr+" -t "+options['numthreads']+" "+protein +" > spaln_"+ output_basename + ".mb"))
#    if iteration == 2:
    spaln_calls.append(str(options['path_to_spaln']+"/spaln -Q7 -LS -O7 -S3 -d "+base+"_"+timestr+" -t "+options['numthreads']+" "+protein +" > spaln_translated_proteins_"+ output_basename + ".fa"))
    for i in spaln_calls:
#	print i
        r = subprocess.Popen(i, shell = True, stdout = null, stderr =null)
        spaln_pids.append(r.pid)
    for j in spaln_pids:
        try:
            os.waitpid(j,0)
        except OSError:
            pass
    return

def spaln_format(spalnfa):
    out = open(basename(spalnfa.split(os.extsep)[0])+"_reformatted.fa",'w')
    for line in open(spalnfa):
        if line.startswith(">"):
            split = line.strip().split()
            out.write(str(split[0].split(".")[0]+";"+split[9]+"\n"))
#            header = line.strip().split(" ")[1]
#            out.write(str(">"+header+"\n"))
        elif line.startswith(";") or re.search("FrameShift:",line) is not None:
            pass 
        else:
            out.write(line)
    out.close()

def spaln_parse(aln_file,iteration):
    with open(aln_file) as aln:
	nested_list = []
	next(aln)
	xenof = {}
	global xeno_count
	coordinate_list = []
	first_protein_hit_check = {}
	for line in aln:
  	    split = line.split()
  	    if split[1] not in cdna_headers and args.choice != 4 and iteration == 1:
		xenof[split[1]]=""
		pass
	    elif line.startswith("@"):
		inside_list = []
		boundaries = re.findall("( [-]?\d+ [-]?\d+ )",line)
		score = str(int(round(float(re.findall("S: (\d+[\.]?\d+)",line)[0])*10)))
		cs,ce = boundaries[0].split()
		coordinate_list.sort(key=int)
		ps,pe = coordinate_list[0],coordinate_list[-1]
		frag_start,frag_end = boundaries[1].split()
                if int(cs)>int(ce):
                    orientation = "revcomp"
                else:
                    orientation = "forward"
  	        del coordinate_list[:]
		inside_list = [prot_id,frag_start,frag_end,ps,pe,contig_id,cs,ce,score,orientation]
		nested_list.append(inside_list)
	    else:
		contig_id = split[1]
		prot_id = split[0]
		coordinate_list.append(split[6])  
		coordinate_list.append(split[7])  
    if iteration == 1 and args.choice!=4:
        for i in xenof:
            xeno_count += 1 
            xenofil.write(i+"\n")
    return nested_list

def scaffold_coverage(list_exon,filehandle,text):
    cov55=0
    cov65=0
    cov75=0
    cov85=0
    cov95=0  
    temp = ""
    dic = {}
    for i in list_exon:  
        split = ('\t'.join(map(str,i))).split()
        if split[0] != temp:
            cov = (int(split[4])-int(split[3]))/float(prot_len[split[0].strip()])
            if (cov >0.9 and cov<=1):
                cov95 +=1
                dic[split[0]]=str("%.2f" % cov)
            elif (cov >0.8 and cov<=.9):
                cov85 +=1
                dic[split[0]]=str("%.2f" % cov)
            elif (cov >0.7 and cov<=.8):
                cov75 +=1
                dic[split[0]]=str("%.2f" % cov)
            elif (cov >0.6 and cov<=.7):
                cov65 +=1
                dic[split[0]]=str("%.2f" % cov)
            elif (cov >0.5 and cov<=.6):
                cov55 +=1
                dic[split[0]]=str("%.2f" % cov)
            else:
                dic[split[0]]=str("%.2f" % cov)
            temp = split[0]
    filehandle.write("\nCoverage statistics "+text+"-scaffolding\n")
    filehandle.write("Number of proteins covered at over 90% coverage:"+"\t"+str(cov95)+"\n")
    filehandle.write("Number of proteins covered at over 80% coverage:"+"\t"+str(cov85+cov95)+"\n")
    filehandle.write("Number of proteins covered at over 70% coverage:"+"\t"+str(cov75+cov85+cov95)+"\n")
    filehandle.write("Number of proteins covered at over 60% coverage:"+"\t"+str(cov65+cov75+cov85+cov95)+"\n")
    filehandle.write("Number of proteins covered at over 50% coverage:"+"\t"+str(cov55+cov65+cov75+cov85+cov95)+"\n")
    return dic


def move_files(f,tardir):
    if os.path.exists(tardir+"/"+f):
        os.remove(tardir+"/"+f)
        shutil.move(f,tardir)
    elif os.path.exists(tardir):
        shutil.move(f,tardir)
    else:
        os.makedirs(tardir)
        shutil.move(f,tardir)


parser=argparse.ArgumentParser(description="Generate gene models")
parser.add_argument('-a' , '--assembly' , type = file , help = 'path to genome assembly file in fasta format' , required = True)
parser.add_argument('-p' , '--pgs' , type = file , help = 'path to protein sequence file' , required = True)
parser.add_argument('-s' , '--config' , type = file , help = 'path to config file', required = True)
parser.add_argument('-c' , '--choice' , type = int , help = '1 if you would like to run the full pipeline \n 2 if you would to like run BLAST elsewhere \n 3 if you have run BLAST elsewhere and want to continue from scaffolding \n 4 if you would like to run without xenobiotic filtering', required = True, choices = range(1,5))
args = parser.parse_args()

prot_len = {}
align_score_spaln = []
best_align={}
pgs_group = []
anon_cont = {}
previous_pgs = ""
previous_end = 0
previous_start =0
scaffold_order ={}
contig_orientation = {}
contig_seq ={}
cdna_headers = {}
next_iteration_answer = {}
contig_scaffold_hash = {}

options = parse_config(args.config)
bufsize = 0 # to flush the output
tlog=open('timelog.txt','w+',bufsize)

protfile = open((os.path.splitext(os.path.basename(args.pgs.name))[0])+"_reformatted.txt",'w+')
# get length of protein seqs
with open(args.pgs.name) as prot:
    for line in prot:
	if line.startswith(">"):
	    original_header = line.replace(">","").strip()
	    reformatted_header = re.sub(r'[\W|\s]+', '_', original_header)[:240]
	    prot_len[reformatted_header] = 0 #protein length information to be later used while calculating coverage
    	    protfile.write(str(">"+reformatted_header+"\n"))
	else:
	    prot_len[reformatted_header]+=len(line.strip()) #protein length information to be later used while calculating coverage
    	    protfile.write(str(line))
protfile.seek(0)

edited_assem = open("assembly_edited.fa",'w+') 
cou = 1 
#get header info from assembly
for rec in SeqIO.parse(args.assembly.name,"fasta"):
    edited_assem.write(str(">contig_"+str(cou)+"\n"+rec.seq+"\n"))
    cdna_headers["contig_"+str(cou)]=[]
    anon_cont["contig_"+str(cou)] = []
    cou += 1
edited_assem.close()
if (args.choice ==1 or args.choice==2 or args.choice==4):
#run spaln
    localtime = time.asctime( time.localtime(time.time()) )
    print "Local current time (starting Spaln) :", localtime
    tlog.write(str("Local current time (starting Spaln) : "+localtime+"\n"))
    spaln(protfile.name, "assembly_edited.fa", 1)
    spaln_format("spaln_translated_proteins_assembly_edited.fa")    
    localtime = time.asctime( time.localtime(time.time()) )
    print "Local current time (finished Spaln):", localtime
    tlog.write(str("Local current time (finished Spaln) : "+localtime+"\n"))

if (args.choice == 2):
    raise SystemExit

#run blastx
if (args.choice==1):
    blast_call = NcbiblastpCommandline(cmd = options['path_to_blastp']+"/blastp", db = options['db'] , query =  "spaln_translated_proteins_assembly_edited_reformatted.fa" , out = "blastresults.tab" , evalue = options['evalue'] , max_target_seqs = "1", num_threads = options['numthreads'], outfmt ="6")
    print "Running blast"
    localtime = time.asctime( time.localtime(time.time()) )
    print "Local current time (starting blast):", localtime
    tlog.write(str("Local current time (starting blast) : "+localtime+"\n"))
    subprocess.call(str(blast_call), shell = True)
    localtime = time.asctime( time.localtime(time.time()) )
    print "Local current time (finished running blast):", localtime
    tlog.write(str("Local current time (finished blast) : "+localtime+"\n")) 
#list of expected taxonomy 

if (args.choice==1 or args.choice==3):
    taxonomy = {}
    with open(options['taxonomyinfo']) as taxfile:
        for line in taxfile:
	    taxcols = line.split("\t")
            taxonomy[taxcols[2].strip()]=""
	    lineage = taxcols[8].split(";")
	    for parent in lineage:
	        taxonomy[parent.strip()]="" 
    taxfile.close()
    fortest = {}
    temp = ""
    with open("blastresults.tab") as blastout:
        for line in blastout:
            cols = line.split("\t")
	    if cols[0].split(";")[0] != temp:
	        temp = cols[0].split(";")[0]
	        m = re.search('Tax=(.*);RepID',cols[1])
	        if float(cols[2]) > 80 and float(cols[11]) > 80 and m is not None:
	            org = m.group(1).replace(";"," ")
	            if org not in taxonomy:
		        cdna_headers.pop(cols[0].split(";")[0],"None")
#extract hits from exonerate_output that aren't xenobiotics
xeno_count = 0
xenofil = open("xeno_filtered.txt",'w')
align_score_spaln = spaln_parse("spaln_assembly_edited.mb",1)
xenofil.close()

all_align = open('spaln_alignment_summary.tab','w')
all_align.write("Protein\tFrag_start\tFrag_end\tProtein_start\tProtein_end\tContig\tContig_start\tContig_end\tScore\tOrientation\tProtein_length\tCoverage\n")
align_score_spaln.sort(key = lambda x : int(x[8]), reverse =True)
align_score_spaln.sort(key = operator.itemgetter(0)) # cluster contigs based on protein id's

for i in align_score_spaln:
    split = "\t".join(i).split()
    cov = ((int(split[4])-int(split[3])+1)*100)/float(prot_len[split[0].strip()])
    string = split[0]+"\t"+split[1]+"\t"+split[2]+"\t"+split[3]+"\t"+split[4]+"\t"+split[5]+"\t"+split[6]+"\t"+split[7]+"\t"+split[8]+"\t"+split[9]+"\t"+str(prot_len[split[0]])+"\t"+str("%.2f" %cov)+"\n"
    all_align.write(string) # spaln alignment summary with coverage information
all_align.close()

pgs_group = list(align_score_spaln)
pgs_group.sort(key = lambda x : int(x[8]), reverse =True)
pgs_group = sorted(pgs_group , key = operator.itemgetter(0))
next_iteration_answer = {}
iteration=0
scaffold_order_protein_start = {}
next_iteration = 0
scaffold_order_protein_end = {}
pgs_group_dup = []

# write out information all best alignments essentially a small subset of alignment summary file 
with open("list_before_scaff.txt",'w') as lis:
    lis.write("Protein\tFrag_start\tFrag_end\tProtein_start\tProtein_end\tContig\tContig_start\tContig_end\tScore\tOrientation\tProtein_length\tCoverage\n")
    for i in pgs_group:
        split = "\t".join(i).split()
        cov = ((int(split[4])-int(split[3])+1)*100)/float(prot_len[split[0].strip()])
        string = split[0]+"\t"+split[1]+"\t"+split[2]+"\t"+split[3]+"\t"+split[4]+"\t"+split[5]+"\t"+split[6]+"\t"+split[7]+"\t"+str(prot_len[split[0]])+"\t"+str("%.2f" %cov)+"\n"
        lis.write(str(string))
	contig_orientation[split[0]+split[5]]=split[9]

pgs_group.sort(key = lambda x : int(x[6]), reverse =True)
pgs_group = sorted(pgs_group , key = operator.itemgetter(0))

for j in list(pgs_group):
    if j[5] in contig_scaffold_hash:
        pass
    else:
        contig_scaffold_hash[j[5]] = []
    if j[0] == previous_pgs: #use bmor:101735811 for test
        if (previous_start < int(j[3]) and previous_end > int(j[4])): # completely included in the current alignment
	    pass
	elif (previous_end - int(j[3]) <= int(options['overlap_dist']) and previous_end - int(j[3]) >= -int(options['disjoint_dist'])): #(25AA overlap or disjoint by 25AA at the end)
	    scaffold_order[j[0]] = str(scaffold_order[j[0]])+"\t"+str(j[5])
	    contig_scaffold_hash[j[5]].append(j[0])
	    previous_end = int(j[4])
	    if next_iteration_question == 1:
	        next_iteration_answer[j[0]]=iteration
		next_iteration=1
	    pgs_group.remove(j)
            pgs_group_dup.append(j)
	elif (int(j[4]) - previous_start <= int(options['overlap_dist']) and int(j[4]) - previous_start >= -int(options['disjoint_dist'])): #(25AA overlap or disjoint by 25AA at the start)
	    scaffold_order[j[0]] = str(j[5]) + "\t" + str(scaffold_order[j[0]])
            contig_scaffold_hash[j[5]].append(j[0])
	    previous_start = int(j[3])
	    if next_iteration_question == 1:
	        next_iteration_answer[j[0]]=iteration	
		next_iteration=1
	    pgs_group.remove(j)
	    pgs_group_dup.append(j)
	elif (previous_start - int(j[4]) >int(options['disjoint_dist']) or int(j[3]) - previous_end >int(options['disjoint_dist'])): #(off by more than 25AA disjoint)
	    next_iteration_question =1
    else:
        previous_pgs = str(j[0])
	previous_start = int(j[3])
	previous_end = int(j[4])
	scaffold_order[j[0]] = j[5]
        contig_scaffold_hash[j[5]].append(j[0])
	next_iteration_question=0
	pgs_group.remove(j)
 	pgs_group_dup.append(j)
    scaffold_order_protein_start[j[0]] = previous_start
    scaffold_order_protein_end[j[0]] = previous_end
while next_iteration==1:
    next_iteration = 0
    iteration += 1 
    for j in list(pgs_group):
	if j[0] in next_iteration_answer and next_iteration_answer[j[0]]==iteration-1:
            if j[0] == previous_pgs: #use bmor:101735811 for test
	        if (previous_start < int(j[3]) and previous_end > int(j[4])): # completely included in the current alignment
	  	    pass
	        elif (previous_end - int(j[3]) <= int(options['overlap_dist']) and previous_end - int(j[3]) >= -int(options['disjoint_dist'])): #(25AA overlap or disjoint by 25AA at the end)
	            contig_scaffold_hash[j[5]] = tempo[j[5]].append(j[0])
	            scaffold_order[j[0]] = str(scaffold_order[j[0]])+"\t"+str(j[5])
	            previous_end = int(j[4])
	            if next_iteration_question == 1:
		        next_iteration_answer[j[0]]=iteration
		        next_iteration=1
		    pgs_group.remove(j)
                    pgs_group_dup.append(j)
		elif (int(j[4]) - previous_start <= int(options['overlap_dist']) and int(j[4]) - previous_start >= -int(options['disjoint_dist'])): #(25AA overlap or disjoint by 25AA at the start)
	            scaffold_order[j[0]] = str(j[5]) + "\t" + str(scaffold_order[j[0]])
	            contig_scaffold_hash[j[5]] = tempo[j[5]].append(j[0])
	            previous_start = int(j[3])
	            if next_iteration_question == 1:
		        next_iteration_answer[j[0]]=iteration	
		        next_iteration=1
		    pgs_group.remove(j)
                    pgs_group_dup.append(j)
	        elif (previous_start - int(j[4]) >int(options['disjoint_dist']) or int(j[3]) - previous_end >int(options['disjoint_dist'])): #(off by more than 25AA disjoint)
	            next_iteration_question = 1
            else:
	        previous_pgs = str(j[0])
	        previous_start = scaffold_order_protein_start[j[0]]
	        previous_end = scaffold_order_protein_end[j[0]]
	        next_iteration_question=0
		pgs_group.remove(j)
                pgs_group_dup.append(j)
 	else:
	    pass
        scaffold_order_protein_start[j[0]] = previous_start
        scaffold_order_protein_end[j[0]] = previous_end    

with open("list_after_scaff.txt",'w') as lis:
    lis.write("Protein\tFrag_start\tFrag_end\tProtein_start\tProtein_end\tContig\tContig_start\tContig_end\tScore\tOrientation\tProtein_length\tCoverage\n")
    for i in pgs_group:
        split = "\t".join(i).split()
        cov = ((int(split[4])-int(split[3]))*100)/float(prot_len[split[0].strip()])
        string = split[0]+"\t"+split[1]+"\t"+split[2]+"\t"+split[3]+"\t"+split[4]+"\t"+split[5]+"\t"+split[6]+"\t"+split[7]+"\t"+split[8]+"\t"+split[9]+"\t"+str(prot_len[split[0]])+"\t"+str("%.2f" % cov)+"\n"
        lis.write(str(string)) # write alignments that are not used for scaffolding

# remove duplicates from scaffold order in the same scaffold
for key in dict(scaffold_order):
    split = dict(scaffold_order)[key].strip().split()
    l = list(OrderedDict.fromkeys(split))
    scaffold_order[key] = '\t'.join(l)

t =  [(k, scaffold_order[k]) for k in scaffold_order]
#print t
t.sort()
scaffold_order_uniq = {}


# remove duplicates from scaffold order spread over different scaffolds
for k, v in t:
    if v in scaffold_order_uniq.values():
        continue
    scaffold_order_uniq[k] = v
 #   print scaffold_order_uniq[k] 

 #check for inverses
for i in scaffold_order_uniq:
    scaffold_order_uniq[i] = scaffold_order_uniq[i].strip().split()

t =  [(k, scaffold_order_uniq[k]) for k in scaffold_order_uniq]
t.sort()
scaffold_order_uniq = {}
for k, v in t:
    if v[::-1] in scaffold_order_uniq.values():
        continue
    scaffold_order_uniq[k] = v


for i in scaffold_order_uniq:
    scaffold_order_uniq[i] = "\t".join(scaffold_order_uniq[i])


#keep only longest scaffold
for key in contig_scaffold_hash:
    count_list = []
    q = []
    if len(contig_scaffold_hash[key])==1:
	pass
    elif len(contig_scaffold_hash[key])>1: # contains list of protein names
	for i in contig_scaffold_hash[key]:
	    if i in scaffold_order_uniq:
	        split = scaffold_order_uniq[i].strip().split()
	        count_list.append(len(split))
		c = 1
	if c == 1 and len(count_list)>1:
	    temp_clist = sorted(list(count_list), reverse = True)
            sorted(count_list, reverse=True)
  	    if temp_clist[0]!=1:
		q = filter(lambda x: x != 1, count_list)
		for i in contig_scaffold_hash[key]:
	            if i in dict(scaffold_order_uniq):
	                split = dict(scaffold_order_uniq)[i].strip().split()
		        if len(split)<max(q):
			    scaffold_order_uniq.pop(i)
	    q = []
	    c = 0

conflict = open("scaffolds_for_manual_curation.txt",'w+')
#combine scaffolds with same contigs
cont_in_ord = {}
for i in dict(scaffold_order_uniq):
    split = scaffold_order_uniq[i].strip().split()
#    split.pop(0)
    ts = list(split)
    for j in ts:#contigs in scaffold
        if j in cont_in_ord: #if exists
            if cont_in_ord[j][0][0]==ts[0]: # check the position of the contig in both scaffold if first contig
		if contig_orientation[cont_in_ord[j][1][0]+cont_in_ord[j][0][0]] == contig_orientation[i+split[0]]:
		    conflict.write(str(",".join(cont_in_ord[j][0])+"\n"+",".join(split)+"\n\n"))
		    split.pop(0)
		    scaffold_order_uniq[i]="\t".join(cont_in_ord[j][0]+split)
		    scaffold_order_uniq.pop(cont_in_ord[j][1][0],None)
		    e=[cont_in_ord[j][0]+split,[i]]
#		    print "1 2 1 3",cont_in_ord[j][0],ts,cont_in_ord[j][0][0],split[0],split
		elif contig_orientation[cont_in_ord[j][1][0]+cont_in_ord[j][0][0]] != contig_orientation[i+split[0]]:
		    split.pop(0)
		    for k in split:
		        if contig_orientation[i+j] == "revcomp":
			    contig_orientation[i+j] == "forward"
			else:
			    contig_orientation[i+j] == "revcomp"
		    scaffold_order_uniq[i]="\t".join(split[::-1]+cont_in_ord[j][0])
		    scaffold_order_uniq.pop(cont_in_ord[j][1][0],None)
	            e=[split[::-1]+cont_in_ord[j][0],[i]]
                for l in cont_in_ord[j][0]:
	            contig_orientation[i+l]=contig_orientation[cont_in_ord[j][1][0]+l]
	        cont_in_ord[j]=e
	    elif cont_in_ord[j][0][-1]==ts[0]:
		d = 0
		if contig_orientation[cont_in_ord[j][1][0]+cont_in_ord[j][0][-1]] == contig_orientation[i+split[0]]:
                    split.pop(0)
		    d = 1
                    scaffold_order_uniq[i]="\t".join(cont_in_ord[j][0]+split)
                    scaffold_order_uniq.pop(cont_in_ord[j][1][0],None)
		if d == 0 and contig_orientation[cont_in_ord[j][1][0]+cont_in_ord[j][0][-1]] != contig_orientation[i+split[0]]:
		    conflict.write(str(",".join(cont_in_ord[j][0])+"\n"+",".join(split)+"\n\n"))
		    split.pop(0)
		    scaffold_order_uniq[i]="\t".join(cont_in_ord[j][0]+split)
		    scaffold_order_uniq.pop(cont_in_ord[j][1][0],None)
                for l in cont_in_ord[j][0]:
  	            contig_orientation[i+l]=contig_orientation[cont_in_ord[j][1][0]+l]
		cont_in_ord[j]=[cont_in_ord[j][0]+split,[i]]
	    elif cont_in_ord[j][0][-1]==ts[-1]:
                if contig_orientation[cont_in_ord[j][1][0]+cont_in_ord[j][0][-1]] == contig_orientation[i+split[-1]] and (len(ts)!=1 and len(cont_in_ord[j][0])!=1) :
		    conflict.write(str(",".join(cont_in_ord[j][0])+"\n"+",".join(split)+"\n\n"))
		    split.pop(-1)
		    scaffold_order_uniq[i]="\t".join(cont_in_ord[j][0]+split[::-1])
		    scaffold_order_uniq.pop(cont_in_ord[j][1][0],None)
#                    print "1 2 3 2",cont_in_ord[j][0],ts,contig_orientation[cont_in_ord[j][1][0]+cont_in_ord[j][0][-1]],contig_orientation[i+split[-1]]
                else:
                    split.pop(-1)
                    for k in split:
                        if contig_orientation[i+j] == "revcomp":
                            contig_orientation[i+j] == "forward"
                        else:
                            contig_orientation[i+j] == "revcomp"
                    scaffold_order_uniq[i]="\t".join(cont_in_ord[j][0]+split[::-1])
                    scaffold_order_uniq.pop(cont_in_ord[j][1][0],None)
                for l in cont_in_ord[j][0]:
	            contig_orientation[i+l]=contig_orientation[cont_in_ord[j][1][0]+l]
	        cont_in_ord[j]=[cont_in_ord[j][0]+split[::-1],[i]]
	    elif cont_in_ord[j][0][0]==ts[-1]:
#		print "1 2 3 1",cont_in_ord[j][0],split
		d = 0
                if contig_orientation[cont_in_ord[j][1][0]+cont_in_ord[j][0][0]] == contig_orientation[i+split[-1]]:
                    split.pop(-1)
		    d = 1
                    scaffold_order_uniq[i]="\t".join(split+cont_in_ord[j][0])
                    scaffold_order_uniq.pop(cont_in_ord[j][1][0],None)
		    e=[split+cont_in_ord[j][0],[i]]
                if d== 0 and contig_orientation[cont_in_ord[j][1][0]+cont_in_ord[j][0][0]] != contig_orientation[i+split[-1]]:
                    conflict.write(str(",".join(cont_in_ord[j][0])+"\n"+",".join(split)+"\n\n"))
                    split.pop(-1)
                    scaffold_order_uniq[i]="\t".join(cont_in_ord[j][0]+split)
                    scaffold_order_uniq.pop(cont_in_ord[j][1][0],None)
		    e=[cont_in_ord[j][0]+split,[i]]
                for l in cont_in_ord[j][0]:
	            contig_orientation[i+l]=contig_orientation[cont_in_ord[j][1][0]+l]
		cont_in_ord[j]=e
	    elif split.index(j) == 0 or split.index(j)+1 == len(split):
		h = cont_in_ord[j][0].index(j)
		a,b = cont_in_ord[j][0][:h],cont_in_ord[j][0][h+1:]
                scaffold_order_uniq[i]="\t".join(a+split+b)
                scaffold_order_uniq.pop(cont_in_ord[j][1][0],None)
                for l in cont_in_ord[j][0]:
	            contig_orientation[i+l]=contig_orientation[cont_in_ord[j][1][0]+l]
		cont_in_ord[j]=[a+split+b,[i]]
	    elif cont_in_ord[j][0].index(j) == 0 or cont_in_ord[j][0].index(j)+1 == len(cont_in_ord[j][0]):
		h = split.index(j)
		a,b = split[:h],split[h+1:]
                scaffold_order_uniq[i]="\t".join(a+cont_in_ord[j][0]+b)
                #scaffold_order_uniq.pop(cont_in_ord[j][1][0],None)
                for l in cont_in_ord[j][0]:
	            contig_orientation[i+l]=contig_orientation[cont_in_ord[j][1][0]+l]
		cont_in_ord[j]=[a+cont_in_ord[j][0]+b,[i]]
	else:
	    cont_in_ord[j]=[split,[i]]

contigs_count = {}
#remove subsets
for key in scaffold_order_uniq:
    contigs = scaffold_order_uniq[key].strip().split()
    for i in contigs:
	if i in dict(contigs_count) and type(contigs_count[i]) is list:
	    z = sorted(contigs_count[i]+[key])
	    contigs_count[i] = '\t'.join(z)
	elif i in dict(contigs_count) and type(contigs_count[i]) is str:
	    print key+"\t"+contigs_count[i]
	else:
	    contigs_count[i]=[key]
	    #print contigs_count[i]

print "\n\n"
for key in dict(contigs_count):
#    print key,contigs_count[key],len(contigs_count[key])
    if len(contigs_count[key])==1:
	contigs_count.pop(key)

placehold = {}
prot_pop = {}
cont_pop = {}
for key in dict(contigs_count):
    if key in contigs_count and contigs_count[key] in placehold:
	split = contigs_count[key].split("\t")
#	print split
	s = scaffold_order_uniq[split[0]].strip().split()
#	print "s", s
	w = scaffold_order_uniq[split[1]].strip().split()
#	print "w", w
	intersect = [x for x in s if x in w]
#	print intersect
	if s.index(intersect[0]) == 0 and s.index(intersect[1])==1:
	    for y in intersect:
		s.pop(s.index(y))
		contigs_count.pop(y,None)
		cont_pop[y]=split
	    new_ord = w+s
	elif w.index(intersect[0]) == 0 and w.index(intersect[1])==1:
	    for y in intersect:
		w.pop(w.index(y))
		contigs_count.pop(y,None)
		cont_pop[y]=split
	    new_ord = s+w
        for l in w:
            contig_orientation[split[0]+l] = contig_orientation[split[1]+l]
        scaffold_order_uniq[split[0]]="\t".join(new_ord)
        scaffold_order_uniq.pop(split[1],None)
    elif key in contigs_count:
	placehold[contigs_count[key]]=key
conflict.seek(0)

count_dups =  {}
prot_for_pop_check = {}
for l in dict(scaffold_order_uniq):
    split = scaffold_order_uniq[l].strip().split()
    for h in split:
        if h in count_dups:
	    if count_dups[h]>len(split):
                scaffold_order_uniq.pop(l,None)
	    else:
		scaffold_order_uniq.pop(prot_for_pop_check[h],None)
        else:
            count_dups[h] = len(split)
	    prot_for_pop_check[h] = l

#build scaffolds
with open("assembly_edited.fa",'r') as assembly:
    contig_rec = SeqIO.parse(assembly, "fasta")
    for rec in contig_rec:
        contig_seq[rec.id]= rec.seq
localtime = time.asctime( time.localtime(time.time()) )
print "Local current time :", localtime
spacer = 'N'*int(options['spacer'])
num_scaffolds = 0
scaff_no = 0
m_scaff_map = {}
with open("onlygenemodels.fa" , 'w') as outfile:
    for key in scaffold_order_uniq:
	seq = ""
        contigs = scaffold_order_uniq[key].strip().split()
	temp_list = list(OrderedDict.fromkeys(contigs))
#	print temp_list
#	print key
	for contig in temp_list:
	    if key+contig not in contig_orientation:
		print key+contig
		a = 0
	    else:
 	        if contig_orientation[key+contig] == "forward":
		    if seq == "":
	                seq = contig_seq[contig]
		    else:
		        seq = seq + spacer + contig_seq[contig] 
	        elif contig_orientation[key+contig] == "revcomp":
		    if seq == "":
	                seq = contig_seq[contig].reverse_complement()
		    else:
		        seq = seq + spacer+ contig_seq[contig].reverse_complement()
	        anon_cont.pop(contig, None)
		a = 1
	if a ==1:
 	    scaff_no += 1
	    header = ">m_scaff_"+str(scaff_no)
	    m_scaff_map["m_scaff_"+str(scaff_no)] = key+"\t"+"\t".join(temp_list)
	    outfile.write(str(header+"\n"+seq+"\n"))
            num_scaffolds = num_scaffolds +1 
z = 0
with open("remaining_contigs.fa",'w+') as remain:
    for i in anon_cont:
	remain.write(str(">"+i+"\n"+contig_seq[i]+"\n"))
remain.close()
filenames = ['onlygenemodels.fa', 'remaining_contigs.fa']
with open('scaffolds.mfa', 'w+') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
with open("order_scaff.txt",'w') as order:
    for i in m_scaff_map:
        order.write(i +"\t"+m_scaff_map[i]+"\n")

print "scaffolds built"
if num_scaffolds==0:
    print "no scaffolds"
    raise SystemExit
#print "next run"
localtime = time.asctime( time.localtime(time.time()) )
#rerun_to_verify
localtime = time.asctime( time.localtime(time.time()) )
print "Local current time (starting Spaln):", localtime
tlog.write(str("Local current time (starting Spaln) : "+localtime+"\n"))
spaln(protfile.name,"scaffolds.mfa",2)
localtime = time.asctime( time.localtime(time.time()) )
print "Local current time (finished running Spaln):", localtime
tlog.write(str("Local current time (finished Spaln) : "+localtime+"\n"))

#presence
align_score_spaln.sort(key = lambda x : int(x[8]), reverse =True)
align_score_spaln = sorted(align_score_spaln, key = operator.itemgetter(0))
text = "Pre"
summary = open('validation_statistics.tab','w+',bufsize)
prefortest = scaffold_coverage(align_score_spaln,summary,text) # pre scaffolding statistics

#statistics post scaffolding
g_validate_post = spaln_parse("spaln_scaffolds.mb",2)
g_validate_post.sort(key = lambda x : int(x[8]), reverse =True)
g_validate_post = sorted(g_validate_post, key = operator.itemgetter(0))
text="Post"
postfortest = scaffold_coverage(g_validate_post,summary,text)
spaln_format("spaln_translated_proteins_scaffolds.fa")    
stop_codon_count = {}
with open("spaln_translated_proteins_scaffolds_reformatted.fa") as translate:
    for line in translate:
        if line.startswith(">"):
	    count = 0
	    split = line.split()
	    header = split[0].replace(">","")[:-2]
	else:
	    count += line.strip().count('X')
	    stop_codon_count[header]=count

if not os.path.exists('mespa_results'):
    os.makedirs('mespa_results')
covtab = open("cov_table.tab",'w+',bufsize)
covtab.write("Protein name\tProtein Length\tCov (Pre_scaffolding)\tCov (Post_scaffolding)\tNo of contigs\tNo of stop codons\n")
scaffold_contig = {} 
l = 0
for i in prefortest:
#    print "pre",i
    if i in postfortest:
#        print "post",i
        if i in scaffold_order_uniq:
            contigs = scaffold_order_uniq[i].split()
            l = len(contigs)
            if l in scaffold_contig:
                scaffold_contig[l]=scaffold_contig[l]+1
            else:
                scaffold_contig[l]=1
#            print i,prot_len[i],prefortest[i],str(l),postfortest[i]
        if l != 0:
	    if i in stop_codon_count:
                covtab.write(i+"\t"+str(prot_len[i])+"\t"+prefortest[i]+"\t"+postfortest[i]+"\t"+str(l)+"\t"+stop_codon_count[i]+"\n")
	    else:
		covtab.write(i+"\t"+str(prot_len[i])+"\t"+prefortest[i]+"\t"+postfortest[i]+"\t"+str(l)+"\tNA\n")
	if l == 0:
	    if i in stop_codon_count:
  	        covtab.write(i+"\t"+str(prot_len[i])+"\t"+prefortest[i]+"\t"+postfortest[i]+"\t1\t"+stop_codon_count[i]+"\n")
	    else:
		covtab.write(i+"\t"+str(prot_len[i])+"\t"+prefortest[i]+"\t"+postfortest[i]+"\t1\tNA\n")
        l = 0

summary.write("\nFragmentation statistics\n")
for key in scaffold_contig:
    if key == 1:
        summary.write("Number of unscaffolded proteins: "+str(scaffold_contig[key])+"\n")
    else:
        summary.write("Number of proteins scaffolded with "+str(key)+" contigs: "+str(scaffold_contig[key])+"\n")
    #duplication
temp = ""
count = 0
potdup = open("potential_duplicates.tab",'w+',bufsize)
potdup.write("Alignment\tGene_ID\tGene_start\tGene_stop\tContig_ID\tContig_start\tContig_stop\tScore\tOrientation\tProtein_length\tAlignment_coverage\n")
with open('spaln_alignment_summary.tab','r') as dup:
    next(dup)
    for line in dup:
	split = line.strip().split()
	if split[0]!=temp and int(split[10]) > options['dup_protein_length'] and float(split[11]) > options['dup_coverage']:
	    previous = line
	    t = 1
	    count = 1
	    temp = split[0]
	elif split[0]==temp and count == 1 and int(split[10]) > options['dup_protein_length'] and float(split[11]) > options['dup_coverage']:
	    if t == 1:
		potdup.write("best\t"+previous) # best alignment
		t = 0
#	    else:
	    potdup.write("potdup\t"+line)
	else:
	    temp = split[0]
	    count = 0
potdup.seek(0)
protfile.seek(0)
tlog.seek(0)
summary.seek(0)
covtab.seek(0)
spaln_format("spaln_translated_proteins_scaffolds.fa")    
spaln_format("spaln_translated_proteins_assembly_edited.fa")
outsum = open("output_summary.txt",'w+',bufsize)
outsum.write("Protein file name: "+args.pgs.name+"\n")
outsum.write("Number of proteins: "+str(len(prot_len))+"\n\nTime taken for analysis\n")
for line in tlog:
    outsum.write(line)
tlog.close()
outsum.write("\n\nValidation statistics\n")
for line in summary:
    outsum.write(line)    
if args.choice!=4:
    outsum.write("Number of contigs filtered during xenobiotic filtering stage: " +str(xeno_count)+"\n\n\n")
outsum.write("\n\nScaffolds that need manual curation\n\n\n")
for line in conflict:
    outsum.write(line)
outsum.close()
subprocess.call("perl " +options['path_to_AsmQC'] + "/AsmQC.pl " + args.assembly.name+" >>output_summary.txt",shell = True)
subprocess.call("perl " +options['path_to_AsmQC'] + "/AsmQC.pl scaffolds.mfa >>output_summary.txt",shell = True)
outsum = open("output_summary.txt",'a+',bufsize)
outsum.write("\n\nPotential duplicates\n")
for line in potdup:
    outsum.write(line)
outsum.write("\n\nCoverage info\n")
potdup.close()
for line in covtab:
    outsum.write(line)
covtab.close()
summary.close()
conflict.close()
scaff_prot_map = {}
prots_in_scaff = open("protein_scaff_relationship.tab","w+")
map_scaff = {}
gf = open("scaffolds.gff")
for line in gf:
    if not line.startswith("##") and c == 1:
        split = line.split()
        if split[2]=="cds":
            c = 0
            t = re.search("Target=(\w+)",split[8])
            s = t.groups()[0]
            if s in map_scaff:
                map_scaff[s]=map_scaff[t.groups()[0]]+"\t"+split[0]
            else:
                map_scaff[s]=split[0]
    elif line.startswith("##"):
        c = 1
for key in map_scaff:
    prots_in_scaff.write(str(key+"\t"+map_scaff[key]+"\n"))
    
prots_in_scaff.close()
move_files("spaln_assembly_edited.gff",'mespa_results')
move_files("spaln_assembly_edited.mb",'mespa_results')
move_files("spaln_scaffolds.mb",'mespa_results')
if args.choice!=4:
    move_files("blastresults.tab",'mespa_results')
move_files("spaln_assembly_edited.aln",'mespa_results')
move_files("spaln_scaffolds.aln",'mespa_results')
move_files((os.path.splitext(os.path.basename(args.pgs.name))[0])+"_reformatted.txt",'mespa_results')
move_files("scaffolds_for_manual_curation.txt",'mespa_results')
move_files("spaln_translated_proteins_assembly_edited.fa",'mespa_results')
move_files("spaln_translated_proteins_assembly_edited_reformatted.fa",'mespa_results')
move_files("spaln_translated_proteins_scaffolds_reformatted.fa",'mespa_results')
move_files("order_scaff.txt",'mespa_results')
move_files("protein_scaff_relationship.tab",'mespa_results')
move_files("list_before_scaff.txt",'mespa_results')
move_files("list_after_scaff.txt",'mespa_results')
move_files("spaln_alignment_summary.tab",'mespa_results')
move_files("potential_duplicates.tab",'mespa_results')
move_files("xeno_filtered.txt",'mespa_results')
move_files("validation_statistics.tab",'mespa_results')
move_files("timelog.txt",'mespa_results')
move_files("spaln_translated_proteins_scaffolds.fa",'mespa_results')
move_files("cov_table.tab",'mespa_results')
move_files("assembly_edited.fa",'mespa_results')
move_files("onlygenemodels.fa",'mespa_results')
move_files("output_summary.txt",'mespa_results')
move_files("remaining_contigs.fa",'mespa_results')
move_files("scaffolds.gff",'mespa_results')
move_files("scaffolds.mfa",'mespa_results')
outsum.close()
