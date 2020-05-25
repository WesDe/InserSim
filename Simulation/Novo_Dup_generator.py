import sys
import csv
from Bio.Seq import Seq
from Bio import SeqIO
import re
from collections import defaultdict
import string
import random
from random import randrange
import getopt
def main():
	print(sys.argv[1:])
	try:
		opts, args = getopt.getopt(sys.argv[1:], "g:b:n:s:p:o:c:t:a:e:i:r:x:v:", ["genome=","bed=","nb_ins=","size_ins=","snp=","out=","chr=","type=","pos=","encode=","insertion=","-random","-novo","vcf="])
	except getopt.GetoptError:
        # print help information and exit:
		#print ('error')  # will print something like "option -a not recognized"
		sys.exit(2)

    # Default parameters
	#print(opts)
	genome=""
	bed_file=""
	num_inser=0
	size_inser=int(0)
	snp=False
	output=""
	chromosome=""
	type_indel=""
	pos_snp=""
	encode_file=""
	insertion_file=""
	novo_genome=""
	randoms=True
	for opt, arg in opts:
		print(opt,arg)
		
		if opt in ('-g',"--genome"):
			genome= arg
			#print(i)
		elif opt in ('-b',"--bed"):
			bed_file= arg
			#print(r)
		elif opt in ('-n',"--nb_ins"):
			num_inser= int(arg)
			#print(i)
		elif opt in ('-s',"--size_ins"):
			size_inser=int(arg)
		elif opt in ('-p',"--snp"):
			snp=arg
		elif opt in ('-o',"--out"):
			output=arg
		elif opt in ('-c',"--chr"):
			chromosome=arg
		elif opt in ('-t',"--type"):
			type_indel=arg
		elif opt in ('-a',"--pos"):
			pos_snp=arg
		elif opt in('-e',"--encode"):
			encode_file=arg
		elif opt in('-i',"--insertion"):
			insertion_file=arg
		elif opt in ("-r","--random") :
			randoms=True
		elif opt in ("-x","--novo") :
			novo_genome=arg
		elif opt in ("-v","--vcf") :
			vcf=arg
		else:
	            assert False, "unhandled option"
	if snp=="True" or snp=="true" :
		snp=True
		if pos_snp!="up" and pos_snp!="down" :
			pos_snp="up"
	else : snp=False
	#print(pos_snp)
	#if genome=="" or num_inser==0 or size_inser==0 or output=="" or type_indel=="" :
	#	assert False,"ERROR mandatory parameters : genome file, number of insertion/deletion, size of insertion/deltion, type of INDEL, output name"
	if bed_file!="" and insertion_file=="":
		print("Exon insertion")
		dictionnary,dictionary2=extract_exon (bed_file,size_inser)
		make_insertion_exonic_inser (genome,dictionnary,type_indel,num_inser,size_inser,snp,output,pos_snp,dictionary2,randoms,novo_genome,vcf)
	if bed_file=="" and encode_file!="" and insertion_file=="":
		print("Gene insertion")
		dictionnary,dictionary2=extract_exon_encode (encode_file,size_inser)
		make_insertion_exonic_inser (genome,dictionnary,type_indel,num_inser,size_inser,snp,output,pos_snp,dictionary2,randoms,novo_genome,vcf)
	elif bed_file=="" :
		print("Whole genome insertion")
		make_insertion_WG (genome,type_indel,num_inser,size_inser,snp,output,pos_snp,randoms,novo_genome,vcf)
	
	#print(r,i)

def extract_exon (bed,size_inser):#,unmasked) :
	bed_parser = csv.reader(open(bed,"r"), delimiter="\t")
	full_position=[]
	dictionary=defaultdict(list)
	dictionary2=defaultdict(list)
	for elt in bed_parser :
		if "@" not in elt[0] and "#" not in elt[0] :
			#print(elt)
			couple=(int(elt[1])+50,int(elt[2])-50)
			key=elt[0].split(' ')[0]
			#print(key)
			exon_size=(int(elt[2])-100)-(int(elt[1])+100)
			if exon_size > 0 :
				dictionary[key].append(couple)	
			if exon_size >= size_inser :
				dictionary2[key].append(couple)
	#print (dictionary)
	return dictionary,dictionary2
	
def extract_exon_encode (encode_file,size_inser):#,unmasked) :
	encode_parser = csv.reader(open(encode_file,"r"), delimiter="\t")
	full_position=[]
	dictionary=defaultdict(list)
	dictionary2=defaultdict(list)
	for elt in encode_parser :
		if "@" not in elt[0] and "#" not in elt[0] :
			#print(elt)
			couple=(int(elt[3])+50,int(elt[4])-50)
			key=elt[1].split(' ')[0]
			#print(key)
			exon_size=(int(elt[4])-100)-(int(elt[3])+100)
			if exon_size > 0 :
				dictionary[key].append(couple)	
			if exon_size >= size_inser :
				dictionary2[key].append(couple)
	#print (dictionary)
	return dictionary,dictionary2

def check_INDEL_pos(size_inser,seen,pos,summ_fail,work_seq) :
	bad=0
	for i in range (size_inser) :
				if pos+i in seen :
					#check+=1
					summ_fail+=1
					bad=1
					break
				elif pos-i in seen :
					#check+=1
					summ_fail+=1
					bad=1
					break
				elif pos+i in seen or work_seq[pos+i]=='N'or work_seq[pos+i]=='Z' :
					summ_fail+=1
					bad=1
					break
	#print(summ_fail,bad)
	return summ_fail,bad			
				
def check_no_N (size_inser,seen,pos,summ_fail,work_seq) :	
	bad=0
	for i in range (50+size_inser) :
		if pos-i-size_inser in seen or pos+i+size_inser in seen :
			#check+=1
			summ_fail+=1
			bad=1
			break
		elif work_seq[pos-i-size_inser]=='N' or work_seq[pos+i+size_inser]=='N' :
			summ_fail+=1
			bad=1
			break
		elif work_seq[pos-i-size_inser]=='Z' or work_seq[pos+i+size_inser]=='Z' :
			summ_fail+=1
			bad=1
			break
	#print(summ_fail,bad)
	return summ_fail,bad
def make_snp (snp) :
	if snp=="A" :
		return "T"
	if snp=="G" :
		return "C"
	if snp=="C" :
		return "G"
	if snp=="T" :
		return "A"
	else :
		raise print( "Error SNP:",snp)
		
			
			
def altered_chrom2(seen,work_seq,out_m,size_inser,elts,old_sequence,head,duplicate_work_seq,out_m_no_ins_snp,seen2,randoms,stranger_sequence,liste_autor,del_writer) :
	count=0
	for position in seen2 :
		ensemble=[]
		ensemble.append(elts)
		Snp=[elts]
		treshold=0
		zygo=random.randint(0,1)
		zygos=str(zygo)+"/1"
		position_duplicated=random.choice(seen)
		size_inser=random.randint(50,1000)

		if count%4==0:
			inser_exon=random.choice(liste_autor)
			starts=inser_exon[0]
			ends=inser_exon[1]
			size_exon=ends-starts
			max_test=0
			while size_exon< size_inser and sequence_inser.count("N")!=0 :
				inser_exon=random.choice(liste_autor)
				starts=inser_exon[0]
				ends=inser_exon[1]
				max_test+=1
			inser=random.randint(starts,ends-(size_inser+1))
			insertion=stranger_sequence[inser:inser+size_inser]
			insertion_liste=list(str(insertion).upper())
			#print ("insertion",insertion)	
			if zygo==1 :
				duplicate_work_seq[position]+=str(stranger_sequence[inser:inser+size_inser])
				work_seq[position]+=str(stranger_sequence[inser:inser+size_inser])
			else :
				work_seq[position]+=str(stranger_sequence[inser:inser+size_inser])

		else :
			sequence_inser=list(old_sequence[position_duplicated:position_duplicated+size_inser])
			if sequence_inser.count("N")!=0 and randoms==True:
				while sequence_inser.count("N")!=0:
					size_inser=random.randint(50,1000)
					sequence_inser=list(old_sequence[position_duplicated:position_duplicated+size_inser])
			if zygo==1 :
				duplicate_work_seq[position]+=str(old_sequence[inser:inser+size_inser])
				work_seq[position]+=str(old_sequence[position_duplicated:position_duplicated+size_inser])
			else :
				work_seq[position]+=str(old_sequence[position_duplicated:position_duplicated+size_inser])


		seq_to_normalize="".join(work_seq[position-30:position+1])
		Normalized_seq,change_pos=left_normalization_insertion_ref(seq_to_normalize,size_inser,31)
		ensemble.append(int(position)+1-change_pos)
		ensemble.append(".")
		ensemble.append(Normalized_seq[0])
		ensemble.append(Normalized_seq)
		#ensemble.append(ensemble[2]-change_pos)
		#ensemble.append(size_inser)
		ensemble.append(".")
		ensemble.append("PASS")
		ensemble.append("TYPE=INS")
		ensemble.append("GT")

		ensemble.append(zygos)
		del_writer.writerow(ensemble)
		#del_writer.writerow(Snp)
		#print(incr)
		count+=1
	try :
		#print(work_seq)
		final=''.join(filter(lambda char: char != 'Z',work_seq))
		finals=''.join(filter(lambda char: char != 'Z',duplicate_work_seq))
	except AttributeError:
		final=string.join(work_seq,'')
		finals=string.join(duplicate_work_seq,'')
	print("final",len(final),"old_length",len(old_sequence))
	out_m_list=[head,final]
	out_m.writerow(out_m_list)
	out_m_lists=[head,finals]
	out_m_no_ins_snp.writerow(out_m_lists)
def left_normalization_insertion_ref(Sequence,Size_insertion,k_mer_size) :
	#print('Sequence : ',Sequence)
	k_mer_begin=Sequence[0:k_mer_size]
	Insertion=Sequence[k_mer_size:]
	#print (k_mer_begin,Insertion)
	k=k_mer_size
	l=len(Insertion)
	i=k_mer_size-1
	j=l-1
	repeatSize=0
	# to handle repeated fuzzy insertion
	while i>=0 and j>=0 :
		if k_mer_begin[i]==Insertion[j] :
			i-=1
			j-=1
			repeatSize+=1
			if j==-1 :
				j=l-1
		else :
			break
		
			
	begin = k_mer_begin[-(repeatSize+1):]
	
	insertion=begin+Insertion
	insertions=insertion[0:len(Insertion)+1]
	#print(begin)
	#print(insertion)
	return (insertions,repeatSize)	

def make_insertion_exonic_inser (genome,dictionnary,type_indel,num_inser,size_inser,snp,output,pos_snp,dictionary2,randoms,novo_genome,vcf):
	full_position=[]
	out_m=csv.writer(open(output+"_No_snp.fa","w"),delimiter="\n")
	if snp==True :
		out_m_snp=csv.writer(open(output+"_SNP.fa","w"),delimiter="\n")
		out_m_no_ins_snp=csv.writer(open(output+"_No_inser_Snp.fa","w"),delimiter="\n")


	else :
		out_m_snp=""
		out_m_no_ins_snp=""
	#ut_unm=csv.writer(unmasked,delimiter="\n")
	del_writer = csv.writer(open("Exon_Ref_"+type_indel+"_"+str(num_inser)+"n_"+str(size_inser)+"Size.vcf","w"), delimiter="\t")
	#head=["#num","CHR","old_pos",'new_pos',"Old_REF",'news_ref',"Left_normalized_ref","Normalized_pos","Size_INDel"]
	head=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","G1"]
	del_writer.writerow(head)
	genome_parser=SeqIO.parse(genome, "fasta")
	
	for record in genome_parser :
		elts=str(record.description)
		head=">" + str(record.description)
		old_sequence=str(record.seq).upper()
		old_length=len(old_sequence)
		sequence=""
		unmasked_seq=""
		num_exon=0
		total_length=0
		work_seq=list(str(record.seq).upper())
		duplicate_work_seq=list(str(record.seq).upper())
		work_seq_snp_no_inser=list(str(record.seq).upper())
		tot=0
		seen=[]
		seen2=[]
		seen_tot=[]
		summ_fail=0
		incr=0
		while tot<num_inser:
			#print(dictionnary)
			check=0
			a=random.choice(dictionary2[elts])	
			bad=0
			#print (a)
			pos=randrange(int(a[0]),int(a[1]))
			#while pos in seen :
			#	a=random.choice(dictionnary[elts])
			#	pos=randrange(int(a[0]),int(a[1]))
			summ_fail,bad=check_INDEL_pos(size_inser,seen,pos,summ_fail,work_seq)
			summ_fail,bad2=check_no_N (size_inser,seen,pos,summ_fail,work_seq)
			#print(summ_fail,bad)
			if bad==0 and bad2==0 :
				#print(a)
				seen.append(pos)
				seen_tot.append(pos)
				dictionary2[elts].remove((int(a[0]),int(a[1])))
					#print(i)
				tot+=1	
			if summ_fail >200 or not dictionary2[elts] :
					print("failed to make enough deletion, max number seq :", tot)
					break
					
		#for duplication insertion simulation
		seen.sort()
		tot=0
		while tot<num_inser :
			check=0
			a=random.choice(dictionnary[elts])	
			bad=0
			#print (a)
			pos=randrange(int(a[0]),int(a[1]))
			summ_fail,bad=check_INDEL_pos(size_inser,seen_tot,pos,summ_fail,work_seq)
			summ_fail,bad2=check_no_N (size_inser,seen_tot,pos,summ_fail,work_seq)
			#print(summ_fail,bad)
			if bad==0 and bad2==0 :
				#print(a)
				seen2.append(pos)
				seen_tot.append(pos)
					#print(i)
				tot+=1	
			if summ_fail >200 or not dictionnary[elts] :
					print("failed to make enough deletion, max number position :", tot)
					break
					

		seen2.sort()
		if len(seen2)!=0 and len(seen)!=0 :
			if type_indel.upper()=="INS":
				altered_chrom2(seen,work_seq,out_m,summ_fail,size_inser,snp,incr,elts,del_writer,old_sequence,head,out_m_snp,duplicate_work_seq,num_inser,out_m_no_ins_snp,work_seq_snp_no_inser,pos_snp,seen2,randoms,novo_genome,seen_tot,vcf)
			if  type_indel.upper() == "DEL" :
				print("make de novo insertion")
				altered_chrom2(seen,work_seq,out_m,summ_fail,size_inser,snp,incr,elts,del_writer,old_sequence,head,out_m_snp,duplicate_work_seq,num_inser,out_m_no_ins_snp,work_seq_snp_no_inser,pos_snp,seen2,randoms,novo_genome,seen_tot,vcf)
		
def make_insertion_WG (genome,type_indel,num_inser,size_inser,snp,output,pos_snp,randoms,novo_genome,vcf):#,unmasked) :
	full_position=[]
	out_m=csv.writer(open(output+"_No_snp.fa","w"),delimiter="\n")
	out_m_no_ins_snp=csv.writer(open(output+"_No_inser_Snp.fa","w"),delimiter="\n")

	del_writer = csv.writer(open("Exon_Ref_"+str(num_inser)+"n_"+str(size_inser)+"Size.vcf","w"), delimiter="\t")
	head=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","G1"]
	del_writer.writerow(head)
	genome_parser=SeqIO.parse(genome, "fasta")
	
	vcf_file =csv.reader(open(vcf,'r'),delimiter='\t')
	liste_autor=[]
	for elt in vcf_file :
		if "!" not in elt[0] and elt[2]=="ORF" :
			start=int(elt[3].split("-")[0])
			end=int(elt[3].split("-")[1])
			if (end-start) > size_inser :
				liste_autor.append((start,end))
	
	de_novo_genome=SeqIO.parse(open(novo_genome), "fasta")
	for records in de_novo_genome :
		stranger_sequence=str(records.seq).upper()
	for record in genome_parser :
		elts=str(record.description)
		head=">" + str(record.description)
		old_sequence=str(record.seq).upper()
		old_length=len(old_sequence)
		sequence=""
		unmasked_seq=""
		num_exon=0
		total_length=0
		work_seq=list(str(record.seq).upper())
		duplicate_work_seq=list(str(record.seq).upper())
		tot=0
		seen=[]
		seen2=[]
		seen_tot=[]
		#print(number)
		summ_fail=0
		while tot<=num_inser :
			check=0	
			#print (a)
			bad=0
			bad2=0
			pos=randrange(0,old_length)
			#print(type(size_inser),type(seen),type(pos),type(summ_fail),type(work_seq))
			
			summ_fail,bad=check_INDEL_pos(size_inser,seen,pos,summ_fail,work_seq)
			summ_fail,bad2=check_no_N (size_inser,seen,pos,summ_fail,work_seq)
			#print(summ_fail,bad)
			if bad==0 and bad2==0 :
				#print(a)
				seen.append(pos)
				seen_tot.append(pos)
					#print(i)
				tot+=1	
			if summ_fail >200 :
					break

		seen.sort()
		#print(seen)	
		#for duplication insertion simulation
		tot=0
		#print ('ok')
		while tot<=num_inser :
			check=0
			pos=randrange(0,old_length)	
			bad=0
			bad2=0
			#print(pos)
			#while pos in seen :
			#	a=random.choice(dictionnary[elts])
			#	pos=randrange(int(a[0]),int(a[1]))
			summ_fail,bad=check_INDEL_pos(size_inser,seen_tot,pos,summ_fail,work_seq)
			summ_fail,bad2=check_no_N (size_inser,seen_tot,pos,summ_fail,work_seq)
			#print(summ_fail,bad)
			if bad==0 and bad2==0 :
				seen2.append(pos)
				seen_tot.append(pos)
					#print(i)
				tot+=1	
			if summ_fail >200 :
					break
					

		seen2.sort()

		seen2.sort()
		print(seen2)
		print(seen)	
		if len(seen2)!=0 and len(seen)!=0 :
			altered_chrom2(seen,work_seq,out_m,size_inser,elts,old_sequence,head,duplicate_work_seq,out_m_no_ins_snp,seen2,randoms,stranger_sequence,liste_autor,del_writer)
	
if __name__ == "__main__":

    main()