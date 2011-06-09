#!/usr/bin/python
"""
Parser for kmers2.py output:
   Functionality:
      -take 4 command line arguments, input1_name, input1 and input2_name, input2
      -read in results from both of these files
      -Generate two results:
         -list of shared repeats
            common_kmers.csv
            kmer,i1_occurrences,i1_forward,i1_revcomp,i1_copynum,  i2_occurrences,i2_forward,i2_revcomp,i2_copynum
         
         -list of unique repeats
            [input1]_unique.csv
               kmer, i1_occurrences, i1_forward, i1_revcomp, i1_copynum 
            [input2]_unique.csv
               kmer, i2_occurrences, i2_forward, i2_revcomp, i2_copynum
"""


import time
import sys
import string
from optparse import OptionParser #http://docs.python.org/library/optparse.html

parser = OptionParser()
parser.add_option('-a',  '--input1')
parser.add_option('-x',  '--input2')
parser.add_option('-b',  '--name1')
parser.add_option('-y',  '--name2')

(options,  args) = parser.parse_args() #uncomment this line for command line support
if(options.input1 is None or options.input2 is None):
   parser.print_help()
   sys.exit()

#if no names are given, use the input file names with extension stripped off
if(options.name1 is None):
   options.name1 = options.input1.strip().split(".")[0]
if(options.name2 is None):
   options.name2 = options.input2.strip().split(".")[0]


# Open input files
try:
   inf1 = open(options.input1, mode='r')
except IOError:
   print "Error opening input file: %s" % options.input1
   sys.exit()

try:
   inf2 = open(options.input2, mode='r')
except IOError:
   print "Error opening input file: %s" % options.input2
   sys.exit()

# Open output files
try:
   outf_common = open("%s_%s_common_kmers.csv" % (options.name1,  options.name2), mode='w')
except IOError:
   print "Error opening input file: %s_%s_common_kmers.csv" % (options.name1,  options.name2)
   sys.exit()

try:
   outf_unique1 = open("%s_unique.csv" % options.name1, mode='w')
except IOError:
   print "Error opening output file: %s_unique.csv" % options.name1
   sys.exit()

try:
   outf_unique2 = open("%s_unique.csv" % options.name2, mode='w')
except IOError:
   print "Error opening output file: %s_unique.csv" % options.name2
   sys.exit()


commonk = {}
unique1k = {}
unique2k = {}
recs1 = 0
recs2 = 0

#Read in first line (column names ) and discard
inf1.readline()
inf2.readline()

#Program flow:
   #Read in all of input1, store it in unique1k 
   #Read in input2, for each line:
      #if it is in unique1k
         #add both of commonk and delete from unique1k
      #else
         #add to unique2k
#Input looks like:
#kmer,occurrences,forward,revcomp,copynum
#GATAGATAGATAGATAGATAGATAG,2554,1472,1082,79
for line in inf1:
   recs1 += 1
   l = line.strip().split(",")
   unique1k[l[0]] = [l[1], l[2], l[3], l[4]]  #[kmer,] i1_occurrences, i1_forward, i1_revcomp, i1_copynum

print("%s records processed from %s\n" % (recs1,  options.input1))

for line in inf2:
   recs2 += 1
   l = line.strip().split(",")
   k = l[0]
   k_rc = k.translate(string.maketrans('ATCG', 'TAGC'))[::-1]
   if(k in unique1k): #should be in common
      #kmer,i1_occurrences,i1_forward,i1_revcomp,i1_copynum,  i2_occurrences,i2_forward,i2_revcomp,i2_copynum
      commonk[k] = [unique1k[k][0], unique1k[k][1], unique1k[k][2], unique1k[k][3], l[1], l[2], l[3], l[4]]
      del unique1k[k]
   elif(k_rc in unique1k): #should be in common
      #kmer,i1_occurrences,i1_forward,i1_revcomp,i1_copynum,  i2_occurrences,i2_forward,i2_revcomp,i2_copynum
      #note: forward and RC must be switched
      commonk[k_rc] = [unique1k[k_rc][0], unique1k[k_rc][1], unique1k[k_rc][2], unique1k[k_rc][3], l[1], l[3], l[2], l[4]]
      del unique1k[k_rc]
   else: #is unique
      unique2k[k] = [l[1], l[2], l[3], l[4]]  #[kmer,] i2_occurrences, i2_forward, i2_revcomp, i2_copynum


print("%s records processed from %s\n" % (recs2,  options.input2))


# Write out unique1k
l = []
#Generate a list of sorted keys
print("Sorting unique1 results....\n")
for kmer in unique1k.keys():
   l.append([kmer, float(unique1k[kmer][3])])

l = sorted(l,  key=lambda counts: counts[1],  reverse=True)

#Write out results
print("Writing unique1 results.....\n")
outf_unique1.write("kmer,occurrences,forward,revcomp,copynum\n")
for r in l:
   kmer = r[0]
   d = unique1k[kmer]
   outf_unique1.write("%s,%s,%s,%s,%s\n" %(kmer,d[0], d[1], d[2], d[3]))

outf_unique1.close()

# Write out unique2k
l = []
#Generate a list of sorted keys
print("Sorting unique2 results....\n")
for kmer in unique2k.keys():
   l.append([kmer, float(unique2k[kmer][3])])

l = sorted(l,  key=lambda counts: counts[1],  reverse=True)

#Write out results
print("Writing unique2 results.....\n")
outf_unique2.write("kmer,occurrences,forward,revcomp,copynum\n")
for r in l:
   kmer = r[0]
   d = unique2k[kmer]
   outf_unique2.write("%s,%s,%s,%s,%s\n" %(kmer,d[0], d[1], d[2], d[3]))

outf_unique2.close()

# Write out commonk: #kmer,i1_occurrences,i1_forward,i1_revcomp,i1_copynum,  i2_occurrences,i2_forward,i2_revcomp,i2_copynum
l = []
#Generate a list of sorted keys
print("Sorting common results....\n")
for kmer in commonk.keys():
   l.append([kmer, (float(commonk[kmer][3]) + float(commonk[kmer][7]))/2 ] )

l = sorted(l,  key=lambda counts: counts[1],  reverse=True)

#Write out results
print("Writing common results.....\n")
outf_common.write("kmer,%s_occurrences,%s_forward,%s_revcomp,%s_copynum,%s_occurrences,%s_forward,%s_revcomp,%s_copynum\n" % (
   options.name1,  options.name1,options.name1,options.name1,options.name2,options.name2,options.name2,options.name2))
for r in l:
   kmer = r[0]
   d = commonk[kmer]
   outf_common.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(kmer,d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]))

outf_common.close()

