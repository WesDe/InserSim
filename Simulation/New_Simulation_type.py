import sys
import csv
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import re
import getopt
import random
from random import randrange
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:g:m:s:t:d:r:me:n:u:", ["vcf=","genome=","mobile=","size","TandemDup=","DisperDup=","TandemRep=","MobileElement=","Novel=","sizeRepeat="])
    except getopt.GetoptError:
        sys.exit(2)
    genome=""
    vcf_file=""
    mobile_element_file=""
    size=0
    td=False
    dd=False
    repeat=False
    mobile=False
    novel=False
    size_r=6
    for opt, arg in opts:
        print(opt,arg)

        if opt in ('-v',"--vcf"):
            vcf_file= arg
        #print(i)
        elif opt in ('-g',"--genome"):
            genome= arg
        elif opt in ('-m', "--mobile"):
            mobile_element_file = arg
        elif opt in ('-t', "--TandemDup"):
            td = True
        elif opt in ('-d', "--DisperDup"):
            dd = True
        elif opt in ('-r', "--TandemRep"):
            repeat = True
        elif opt in ('-me', "--MobileElement"):
            mobile = True
        elif opt in ('-n', "--Novel"):
            novel = True
        elif opt in ('-u', "--sizeR"):
            size_r = int(arg)
        elif opt in ('-s', "--size"):
            size = int(arg)
        else :
            assert False, "unhandled option"

    if genome!="" and vcf_file!="" and size!=0 :
        liste_position=get_position(vcf_file)
        if mobile_element_file !="" and mobile == True :
            chrom, work_seq = get_seq(genome)
            liste_me_valid=get_valid_ME(mobile_element_file,chrom,size)
            inser_types(work_seq, chrom, liste_position, "ME", liste_seq_me)
            liste_seq_me = mobile_element(work_seq, liste_me_valid, liste_position)
        if repeat==True :
            chrom, work_seq = get_seq(genome)
            liste_STR = STR_pattern(work_seq, size, liste_position,size_r)
            inser_types(work_seq, chrom, liste_position, "STR", liste_STR)
        if td==True :
            chrom, work_seq = get_seq(genome)
            inser_types(work_seq, chrom, liste_position, "TD", liste_TD)
    else :
        print ("Missing input")

def get_position(vcf_file) :
    inpt=csv.reader(open(vcf_file,"r"),delimiter="\t")
    liste=[]
    next(inpt)
    for elt in inpt :
        pos=int(elt[1])
        #seq=elt[4]
        liste.append(pos)
    return liste

def get_valid_ME(ME_file,chrom,size) :
    liste_valid=[]
    inpt = csv.reader(open(ME_file, "r"), delimiter="\t")
    for elt in inpt :
        if "#" not in elt[0] or "@" not in elt[0] and re.sub("chr", "", elt[0].description)==chrom and elt[-1]=="+" and "Alu" in elt[3]:
            #print(elt[1],elt[2])
            size_me=int(elt[2])-int(elt[1])
            #print(size_me,size)
            if size_me >int(size)-50 and size_me <int(size)+50 :
                liste_valid.append((int(elt[1]),int(elt[2])))
    
    return liste_valid

def get_seq(fasta_file):
    reader = SeqIO.parse(fasta_file, "fasta")
    for element in reader:
        chrom = re.sub("chr", "", element.description)
        work_seq = list(str(element.seq).upper())
    return (chrom,work_seq)

def STR_pattern(seq,size,liste,size_r):
    insertion_list=[]
    for elt in liste :
        pattern="".join(seq[elt:elt+size_r])
        motif=pattern
        while len(motif)<size :
            motif+=pattern
        insertion_list.append(motif)
    return insertion_list

def tandem_dup(seq, size,liste):
    insertion_list=[]
    for elt in liste :
        insertion_list.append("".join(seq[elt:elt+size]))
    return insertion_list

def mobile_element(seq,liste_me_valid,liste):
    insertion_list=[]
    for i in range(0,len(liste)):
        choice_alu = random.choice(liste_me_valid)
        insertion_list.append("".join(seq[choice_alu[0]:choice_alu[1]]))
    return insertion_list


def inser_types(work_seqs,chrom,liste,types,inser_list) :
    otp_hom=csv.writer(open("HOM_"+chrom+"_"+types+".fa","w"),delimiter="\n")
    otp_het = csv.writer(open("HET_"+chrom+"_"+types+".fa", "w"), delimiter="\n")
    otp_vcf = csv.writer(open("HOM_"+chrom+"_"+types+".vcf", "w"), delimiter="\t")
    old_seq="".join(work_seqs)
    work_seq=work_seqs
    header=">"+chrom
    head = ["#CHROM", "POS", "ID", "REF", "ALT","QUAL", "FILTER", "INFO", "FORMAT", "G1"]
    otp_vcf.writerow(head)
    for i in range(0,len(liste)) :
        inser = work_seq[liste[i]]+inser_list[i]
        work_seq[liste[i]]=inser
        liste_vcf_otp = [chrom, liste[i], ".", work_seq[liste[i]][0],inser, ".", "PASS", "TYPE=INS", "GT", "1/1"]
        otp_vcf.writerow(liste_vcf_otp)
    new_seq = "".join(work_seq)
    hom_list=[header,new_seq,header,new_seq]
    het_list=[header,new_seq,header,old_seq]
    otp_hom.writerow(hom_list)
    otp_het.writerow(het_list)
    print("Old size : ",len(old_seq), "New size :", len(new_seq))

		
if __name__ == "__main__":
    main()
