	
import getopt
import sys
import csv
import statistics
import re
import gzip
csv.field_size_limit(sys.maxsize)


def main():
	print(sys.argv[1:])
	try:
		opts, args = getopt.getopt(sys.argv[1:], "t:v:m:f:c:", ["truth=","vcf=","microhomology","filter","compress"])
	except getopt.GetoptError:
        # print help information and exit:
		#print ('error')  # will print something like "option -a not recognized"
		sys.exit(2)

    # Default parameters
	#print(opts)
	truth_vcf=""
	vcf=""
	mh=0
	filtering="PASS"
	compress=False
	for opt, arg in opts:
		print(opt,arg)
		if opt in ('-t',"--truth"):
			truth_vcf= arg
			#print(i)
		elif opt in ('-v',"--vcf"):
			vcf= arg
		elif opt in ('-m',"--micrhomology"):
			mh= int(arg)
			#print(r)
		elif opt in ('-f',"--filter"):
			filtering=""
		elif opt in ('-c',"--compress"):
			compress= True
			#print(r)
		else:
	            assert False, "unhandled option"

	ref=dictionary_ref(truth_vcf)
	#dic_auto=identification(vcf)
	comparaisonRef_Vcf(ref,vcf,mh,filtering,compress)
	#return_found(liste_founds, lec_2)


def identification(vcf_file) :
	dic_autorized={}
	liste_valid=[]
	vcf_reader=csv.reader(open(vcf_file,"r"),delimiter='\t')
	for element in vcf_reader :
		if "#" not in element[0] :
			key=element[2]
			dic_autorized[key]=dic_autorized.setdefault(key,0)+1
	for a in dic_autorized :
		if dic_autorized[a]==2 or dic_autorized[a]==1:
			liste_valid.append(a)
	vcf_file.seek(0)
	return liste_valid

def is_valid(inser):
    allowed = "ATCGatcg"
    if all(c in allowed for c in inser):
        return True
    else:
        return False

def dictionary_ref(text_file) :
	vcf_reader=csv.reader(open(text_file,"r"),delimiter='\t')
	dic={}
	for element in vcf_reader :
		#print(element)
		if "#" not in element[0] and "@" not in element[0] and len(element[4])>50	:	
			key = (re.sub("chr", "", element[0]), int(element[1]))
			value =element[4]
			dic.setdefault(key,[]).append((value,element[9]))
	return dic

def check(chrom,pos,elements,elt,total_insertion,insertion_passed,bkpt_found,liste_found,liste_bkpt_unique,list_bkpt_seen_found,sum_seq,sum_site,mh,insertion_90,count_small_false,total_bkpt) :
	verif=False
	total_insertion+=1
	insertion_passed+=1
	key=(len(elt)-len(elements[3]))
	allowed="nNATCGatcg"
	if all(c in allowed for c in elt ) :
		insertion=elt
	else :
		insertion=""
	if (chrom,elements[1]) not in liste_bkpt_unique :
		liste_bkpt_unique.append((chrom,elements[1]))
		total_bkpt+=1
	for i in range(-mh,mh):
		newpos=pos+int(i)
		keys=(chrom,int(newpos))
		if keys in dicMills :
			verif=True
			if (chrom,newpos) not in list_bkpt_seen_found:
				result=dicMills[keys][0][0]
				z=dicMills[keys][0][1]
				bkpt_found+=1
				liste_found.append(keys)
				list_bkpt_seen_found.append((chrom,newpos))
				sum_site+=1
				if insertion =="" :
					break
				alignement_first=getAln(result,insertion)
				if alignement_first>=90 :
					sum_seq+=1
					insertion_90 += 1
					break
	if verif==False  :
		count_small_false+=1
	return total_insertion,insertion_passed,bkpt_found,liste_found,liste_bkpt_unique,list_bkpt_seen_found,sum_seq,sum_site,mh,insertion_90,count_small_false,total_bkpt

def check_valid(liste,allowed):
    for elt in liste :
        if all(c in allowed for c in elt) :
            return elt
def comparaisonRef_Vcf (dicMills,vcf_file,mh,filtering,compress) :
	if compress :
		vcf_reader=csv.reader(gzip.open(vcf_file,mode='rt'),delimiter='\t')
	else :
		vcf_reader=csv.reader(open(vcf_file,"r"),delimiter='\t')
	bkpt_found=0
	total_bkpt =0
	total_insertion=0
	list_bkpt_seen=[]
	list_bkpt_seen_found=[]
	total_length=[]
	insertion_90=0
	liste_bkpt_unique=[]
	liste_unique=[]
	insertion_passed=0
	count_small_false=0
	allowed="nNATCGatcg"
	liste_found=[]
	sum_seq=0
	sum_site=0
	for elements in vcf_reader :
		result=""
		if '#' not in elements[0] and "@" not in elements[0] : 
			chrom = re.sub("chr", "", elements[0])
			pos=int(elements[1])
			sep=elements[4].split(',')
			for elt in sep :
				verif=False
				insertion=check_valid(re.findall(r"[\w']+", elt),allowed)
				if (is_valid(elements[3]) and is_valid(insertion) and len(elements[3]) < len(insertion)) or elt == "<INS>" or "SVTYPE=INS" in elements[7] or "INSERTION=" in elements[7]:
					if filtering =="PASS" :
						if elements[6] == "PASS" or elements[6]=="LOWAS" :
							total_insertion+=1
							insertion_passed+=1
							key=(len(elt)-len(elements[3]))
							allowed="nNATCGatcg"
							if (chrom,elements[1]) not in liste_bkpt_unique :
								liste_bkpt_unique.append((chrom,elements[1]))
								total_bkpt+=1
							for i in range(-mh,mh):
								newpos=pos+int(i)
								keys=(chrom,int(newpos))
								if keys in dicMills :
									verif=True
									if (chrom,newpos) not in list_bkpt_seen_found:
										result=dicMills[keys][0][0]
										z=dicMills[keys][0][1]
										bkpt_found+=1
										liste_found.append(keys)
										list_bkpt_seen_found.append((chrom,newpos))
										sum_site+=1
										if insertion =="" :
											break
										alignement_first=getAln(result,insertion)
										if alignement_first>=90 :
											sum_seq+=1
											insertion_90 += 1
											break
							if verif==False  :
								count_small_false+=1
					else :
						total_insertion+=1
						insertion_passed+=1
						key=(len(elt)-len(elements[3]))
						allowed="nNATCGatcg"
						if (chrom,elements[1]) not in liste_bkpt_unique :
							liste_bkpt_unique.append((chrom,elements[1]))
							total_bkpt+=1
						for i in range(-mh,mh):
							newpos=pos+int(i)
							keys=(chrom,int(newpos))
							if keys in dicMills :
								verif=True
								if (chrom,newpos) not in list_bkpt_seen_found:
									result=dicMills[keys][0][0]
									z=dicMills[keys][0][1]
									bkpt_found+=1
									liste_found.append(keys)
									list_bkpt_seen_found.append((chrom,newpos))
									sum_site+=1
									if insertion =="" :
										break
									alignement_first=getAln(result,insertion)
									if alignement_first>=90 :
										sum_seq+=1
										insertion_90 += 1
										break
						if verif==False  :
							count_small_false+=1
	#print ('total_bkpt : ', total_bkpt, 'total_insertion : ', total_insertion, 'insertion > 95 % ', insertion_95, 'bkpt_found :', bkpt_found,'PASSED insertion_100: ', insertion_98,"FAILED insertion exact",good_failed,"Total Failed",failed_tot,"good_pos_good_length: ",good_pos_good_length) #,'insertion_50_80 :',insertion_90,'insertion_10_50 :',insertion_50,"good_pos_good_length: ",good_pos_good_length)
	print("Recall site : ", sum_site/len(dicMills), "Recall seq :",sum_seq/len(dicMills))#, "Sum_geno :", sum_genot/len(dicMills))
	print ("False positive discovery", count_small_false)

	liste_outp=[]
	#if insertion_90!=0 :
	#	precision =insertion_90/insertion_passed
	#	recall = insertion_90/len(dicMills)
	#	recall_bkpt=insertion_90/bkpt_found
	#	recall_soft=bkpt_found/len(dicMills)
	#	precision_soft=bkpt_found/insertion_passed


def return_found(liste_found, vcf_file):
	vcf_reader = csv.reader(vcf_file, delimiter='\t')
	otp_found = csv.writer(
		open("Found_Exon_250_size_final.vcf", "w"), delimiter="\t")
	for elts in vcf_reader:
		#print(liste_found)
		if "#" not in elts[0]:
			head = (re.sub("chr", "", elts[0]), int(elts[1]))
			print(head)
			if head in liste_found	:
				otp_found.writerow(elts)
def revcomp(sequence):
    baseComplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    rc=''.join([baseComplement[base] for base in reversed(list(sequence))])
    return rc

def initMatrix(S,T):
	matrix = [0]*(len(S)+1)
	for i in range(len(S)+1):
		matrix[i]=[0]*(len(T)+1)
	return matrix

def max3(a,b,c):
	if a>b:
		if a>c: return a
		else: return c
	else:
		if b>c:return b
		else: return c

def score(s,t,match,mismatch):
	if s==t:return match
	else: return mismatch

def initGlobal(matrix,indel):
	for i in range(1,len(matrix)):
		matrix[i][0]=matrix[i-1][0]+indel
	for j in range(1,len(matrix[0])):
		matrix[0][j]=matrix[0][j-1]+indel


def fillNW(matrix,indel,match,mismatch,S,T):
	for i in range(1,len(S)+1):
		for j in range(1,len(T)+1):
			## attention indices décalés pour matrix et S/T
			matrix[i][j]=max3(matrix[i-1][j-1]+score(S[i-1],T[j-1],match,mismatch),matrix[i-1][j]+indel,matrix[i][j-1]+indel)
	return matrix[len(S)][len(T)]


def whichMax3(a,b,c):
	if a>b:
		if a>c: return 1
		else: return 3
	else:
		if b>c:return 2
		else: return 3


def getAln(S,T):
	match=3
	mismatch=-1
	indel=-2
	#print('creation matrix')
	matrix=initMatrix(S,T)
	#print('init')
	initGlobal(matrix,indel)
	#print('fill')
	a=fillNW(matrix,indel,match,mismatch,S,T)
	alS=""
	alT=""
	matching=""
	i=len(S)
	j=len(T)
	ma=0
	oth=0
	#print ('start align')
	#print(len(matrix))
	while i>0 and j>0:
		coming=whichMax3(matrix[i-1][j]+indel,matrix[i][j-1]+indel,matrix[i-1][j-1]+score(S[i-1],T[j-1],match,mismatch))
		#print(coming)
		if coming==3: # match
			#alS=S[i-1]+alS
			#alT=T[j-1]+alT
			if S[i-1]==T[j-1]: 
				ma+=1
				#matching="|"+matching
			else: 
				#matching="."+matching
				oth+=1
			i=i-1
			j=j-1
			continue
		elif coming==1: # horizontal : gap dans T
			#alS=S[i-1]+alS
			#alT="-"+alT
			i=i-1
			#matching=" "+matching
			oth+=1
			continue
		else: # vertical : gap dans S
			#alT=T[j-1]+alT
			#alS="-"+alS
			j=j-1
			#matching=" "+matching
			oth+=1
	while i>0:
		#alS=S[i-1]+alS
		#alT="-"+alT
		i=i-1
		#matching=" "+matching
		oth+=1
	while j>0:
		#alT=T[j-1]+alT
		#alS="-"+alS
		j=j-1
		#matching=" "+matching
		oth+=1
	# print i,j	
	#return (alS,matching,alT)
	#print(ma,i,len(S))
	if ma!=0 :
		tot=(ma/max(len(T),len(S)))*100
	else :
		tot=0
	#print("pourcentage :",tot,"match:",ma,"others:",oth)
	return tot


	
if __name__ == "__main__":

    main()

