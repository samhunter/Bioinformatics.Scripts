#!/usr/bin/python
"""
Implement a k-mer frequency calculator using dictionaries.
This is the extended version with the following improvements:
   -command line support for input fasta file
   -command line support for output file
   -count k-mer and its complement/reverse/reversecompliment as the same (i.e. ACGT = TGCA)

"""
import time
import sys
import string
from optparse import OptionParser #http://docs.python.org/library/optparse.html

parser = OptionParser()
parser.add_option('-i',  '--input')
parser.add_option('-o',  '--output')
parser.add_option('-k',  '--kmer',  action='store', type='int')

(options,  args) = parser.parse_args() #uncomment this line for command line support
#(options,  args) = parser.parse_args([ '-k', '6','-i', 'in.fa', '-o', 'out'])
if(options.input is None or options.output is None or options.kmer is None):
   parser.print_help()
   sys.exit()

k = options.kmer

try:
   inf = open(options.input, mode='r')
except IOError:
   print "Error opening input file: %s" % options.input
   sys.exit()

try:
   outf = open("%s-out-%smers.csv" % (options.output, k), mode="w")
except IOError:
   print "Error opening output file: %s" % options.output
   sys.exit()

#Iterate over lines in the file and count up k-mers
counts = {} #stores global counts for all k-mers
rcounts = {} #stores read-count to track copy number
j = 0
t = time.time()
startt = t
reads = 0

#total,kmer, revcomp, copynum
for line in inf:
  line = line.strip().upper()
  if(line[0] != ">"):
    reads += 1
    for i in range(len(line) - k +1): #add one to ensure analysis of full line
      #check global counts
      kmer = line[i:i+k]
      kmer_c = kmer.translate(string.maketrans('ATCG','TAGC'))#complement 
      #kmer_r = kmer[::-1]#reverse
      kmer_rc = kmer_c[::-1]#reverse complement
      if(kmer in counts):
         counts[kmer][0] += 1
         counts[kmer][1] += 1
         if(kmer not in rcounts):
            counts[kmer][3] += 1
            rcounts[kmer] = True
      elif(kmer_rc in counts):
         counts[kmer_rc][0] += 1
         counts[kmer_rc][2] += 1
         if(kmer_rc not in rcounts):
            counts[kmer_rc][3] += 1
            rcounts[kmer_rc] = True
      else: #total,kmer, revcomp, copynum
        counts[kmer] = [1, 1, 0, 1]
        rcounts[kmer] = True
    rcounts = {}
    j+=1
    if(j % 1000 == 0):
      print("record %s, %s records per second" % (j, 1000/(time.time() -t)))
      t = time.time()
    #if(j > 2): break #uncomment for testing

endt = time.time()

print("Processed %s reads in %s seconds" % (reads,  endt-startt))

l = []
#Generate a list of sorted keys
print("Sorting results....\n")
for kmer in counts.keys():
   l.append([kmer, counts[kmer][3]])

l = sorted(l,  key=lambda counts: counts[1],  reverse=True)

#Write out results
print("Writing results.....\n")
outf.write("kmer,occurrences,forward,revcomp,copynum\n")
for r in l:
   kmer = r[0]
   d = counts[kmer]
   outf.write("%s,%s,%s,%s,%s\n" %(kmer,d[0], d[1], d[2], d[3]))

outf.close()
