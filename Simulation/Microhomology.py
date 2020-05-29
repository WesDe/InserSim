#!/usr/bin/env python3

"""*******************************************************************************
	Name: InserSim 
	Description: InserSim aims to simulate insertion and benchmark insertion in vcf file
	Author: Wesley Delage
	Contact: wesley.delage@irisa.fr, IRISA/Univ Rennes/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
	
	Copyright (C) 2020 Inria
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Affero General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Affero General Public License for more details.
	You should have received a copy of the GNU Affero General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************"""
import sys
import csv
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import re
import getopt

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:g:m:t:", ["vcf=","genome=","micro_size=","inser_type="])
    except getopt.GetoptError:
        sys.exit(2)
    genome=""
    vcf_file=""
    homology_size=0
    inser_type=""
    for opt, arg in opts:
        print(opt,arg)

        if opt in ('-v',"--vcf"):
            vcf_file= arg
        #print(i)
        elif opt in ('-g',"--genome"):
            genome= arg
        elif opt in ('-m',"micro_size") :
            homology_size=int(arg)
        elif opt in ('-t',"inser_type") :
            inser_type=arg
        else:
            assert False, "unhandled option"

    if genome!="" and vcf_file!="" and inser_type!="" :
        recorded = get_position(vcf_file)
        chromo,work=get_seq(genome)
        inser_micro(work,chromo,recorded,homology_size,inser_type)
    else :
        print ("Missing input")
def get_position(vcf_file) :
    inpt=csv.reader(open(vcf_file,"r"),delimiter="\t")
    liste=[]
    next(inpt)
    for elt in inpt :
        pos=int(elt[1])
        seq=elt[4]
        liste.append((pos,seq))
    return liste


def get_seq(fasta_file):
    reader = SeqIO.parse(fasta_file, "fasta")
    for element in reader:
        chrom = re.sub("chr", "", element.description)
        work_seq = list(str(element.seq).upper())
    return (chrom,work_seq)

def inser_micro(work_seq,chrom,liste,micro_size,types) :
    otp_hom=csv.writer(open("HOM_"+chrom+"_"+types+"_"+str(micro_size)+"MH.fa","w"),delimiter="\n")
    #otp_het = csv.writer(open("HET_"+chrom+"_"+types+"_" +str(micro_size)+"MH.fa", "w"), delimiter="\n")
    old_seq="".join(work_seq)
    header=">"+chrom
    head = ["#CHROM", "POS", "ID", "REF", "ALT",
            "QUAL", "FILTER", "INFO", "FORMAT", "G1"]
    otp_vcf = csv.writer(
        open("HOM_"+chrom+"_"+types+"_"+str(micro_size)+ ".vcf", "w"), delimiter="\t")
    for elt in liste :
        
        inser = "".join(work_seq[elt[0]:elt[0]+micro_size])+elt[1][0:len(elt[1])-micro_size]
        work_seq[elt[0]]=inser
        liste_vcf_otp = [chrom, elt[0], ".", elt[1][0],
                         inser, ".", "PASS", "TYPE=INS", "GT", "1/1"]
        otp_vcf.writerow(liste_vcf_otp)
    new_seq = "".join(work_seq)
    hom_list=[header,new_seq,header,new_seq]
    #het_list=[header,new_seq,header,old_seq]
    otp_hom.writerow(hom_list)
    #otp_het.writerow(het_list)
    print("Old size : ",len(old_seq), "New size :", len(new_seq))

		
if __name__ == "__main__":
    main()
