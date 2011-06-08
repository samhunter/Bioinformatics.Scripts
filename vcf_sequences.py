#!/usr/bin/python
"""
This script reads in a vcf formatted result file plus a fasta/multiple fasta file. 
The m-fasta file is stored in a dictionary with the contig name as the key.
The VCF file is then parsed and each entry is used at the ">name" field for a new m-fasta file which contains 
details about the variant, as well as the +/-N bases that flank it (where N is set on the command line).

Output: 

   fasta input looks like:
      >GKF1FGD02EGOM3 0 0 0 1 79
      gactacacgtagtatCTCAATGATTCATACAGCTTCACATCGAGAGCAT
"""

import time
import sys
import string
from optparse import OptionParser #http://docs.python.org/library/optparse.html
#from Bio import pairwise2
#import thread

starttime = time.time()
parser = OptionParser()
parser.add_option('-c',  '--contig')
parser.add_option('-v',  '--vcf')
parser.add_option('-o',  '--output')
parser.add_option('-l',  '--length')

(options,  args) = parser.parse_args() #uncomment this line for command line support
#(options,  args) = parser.parse_args([ '-c', '454AllContigs.fna', '-v','R_filtered_variants_mod.vcf','-l','100','-o','variant_flanking_sequence.fa'])
if(options.contig is None or options.vcf is None or options.output is None):
   parser.print_help()
   sys.exit()

verbose = True

def read_fasta_to_dictionary(filename=""):
   try:
      fasta_inf = open(filename, mode='r')
   except IOError:
      print "Error opening input file: %s\n" % filename
      sys.exit()
   
   sequences = {}   
   rid = ''; seq=''
   for l in fasta_inf:
      if(l[0] == '>'):
         if(seq != ''):
            sequences[rid]['seq'] = seq
            rid=''
            seq=''
         #parse header
         l = l[1:].strip().split()
         rid = l[0]
         sequences[rid] = {} 
         for e in l[1:]:
            e = e.split("=")
            sequences[rid][e[0]] = e[1]
      else:
         seq = seq + l.strip()
   
   fasta_inf.close()
   return sequences

def read_vcf_to_list(filename=""):
   try:
      vcf_inf = open(filename,  mode='r')
   except IOError:
      print "Error opening input file: %s\n" % filename
      sys.exit()
   vcf = []
   vcf.append([])
   n_samps = 0
   for l in vcf_inf:
      if(l[0:2] == "#C"):
         fields = l[1:].strip().split("\t")
         vcf[0] = fields
      if(l[0] != "#"):
         data = l.strip().split("\t")
         rec = {}
         for i in range(len(data)):
            rec[fields[i]] = data[i]
         vcf.append(rec)
   return (vcf)


def main():
   #options.contig is None or options.vcf is None or options.output is None
   sequences =  read_fasta_to_dictionary(filename=options.contig)
   vcf = read_vcf_to_list(filename=options.vcf)

   try:
      outf = open(options.output,  mode='w')
   except:
      print "Error opening output file: %s\n" % options.output
      sys.exit()

   #Iterate over everything in the vcf list
   for i in range(1,  len(vcf)):
      chrom = vcf[i]['CHROM']
      pos = int(vcf[i]['POS'])
      alt = vcf[i]['ALT']
      ref = vcf[i]['REF']
      l = int(options.length)
      genotypes = ""
      score = 0
      for j in vcf[0][9:]: #genotypes of ALL
         genotypes += j + ":" + vcf[i][j][0:3] + "_"
      for j in vcf[0][10:]: #score is based on everything but the first (ancestral) entry
         score += sum(map(int, vcf[i][j][0:3].split("/")))
      genotypes = genotypes[:-1]
      if l > pos: 
         lc = pos 
      else: lc = l+1
      if l+pos > len(sequences[chrom]['seq']): 
         rc = len(sequences[chrom]['seq']) - pos 
      else: rc = l-1
      snp_seq = sequences[chrom]['seq'][pos-lc:pos+rc]
      snp_seq = snp_seq[0:lc-1] + " " + snp_seq[lc-1:]
      outf.write(">contig:%s_pos:%s_ref:%s_alt:%s_score:%s_genotypes:%s\n" % (chrom,  pos,  ref,  alt, score, genotypes))
      outf.write("%s\n" % snp_seq)
      
   outf.close()

main()
