#!/usr/bin/python
"""
This is a re-do and combination of kmers2.py and parse_kmers.py.  It combines the functionality of both, but rather than generate large
results files, it keeps everything in memory at each step, increasing speeds and efficiency.  It is also set up to handle a comparison for
only a single K, or for a range of them.  Additionally, it expects inputs as a list of fasta files, rather than a single fasta file.

Inputs:
   -file containing list of fasta files for group 1
   -file containing list of fasta files for group 2
   -a k-mer size to start at
   -a k-mer size to end at (if not supplied then generate output for single k-size only)
   -

Process: 
   -read command line input and parse
   -for each k in range(start, end+1)
      -for each file in group1 list
         -count kmers, store in unique1k
      -for each file in group2 list
         -count kmers, if kmer in unique1k, move to commonk, else store in unique2k
      -write to results.csv, [k, len(commonk), len(unique1k), len(unique2k)]
      -if options.verbose write options.groups[0]_unique_k.csv, options.groups[1]_unique_k.csv,  groups[1,2]_common_k.csv
"""

import time
import sys
import string
from optparse import OptionParser #http://docs.python.org/library/optparse.html

parser = OptionParser()
#parser.add_option('-a', '--input1')
#parser.add_option('-b', '--input2')
parser.add_option('-l', '--files', help="two file names for lists of fasta inputs", nargs=2)
parser.add_option('-s',  '--start', help="smallest kmer", action='store', type='int')
parser.add_option('-e',  '--end', help="largest kmer (optional)", action='store', type='int')
parser.add_option('-g', '--groups', help="group names (in order)", nargs=2, default=('g1','g2'))
parser.add_option("-v", action="store_true", dest="verbose", default=False)
#parser.add_option("-q", action="store_false", dest="verbose")

(options,  args) = parser.parse_args() #uncomment this line for command line support
#(options,  args) = parser.parse_args([ '-l', 'file1','file2','-s', '6', '-e', '10', '-g', 'male','female'])
#(options,  args) = parser.parse_args([ '-l', 'g1.lst','g2.lst','-s', '6', '-e', '10'])

if(options.files is None or options.start is None):
   parser.print_help()
   sys.exit()

if(options.end is None):
   options.end = options.start

if(options.end < options.start):
   print("Start and end values don't make sense.")
   sys.exit()
   
#Try opening list files:
try:
   g1_inf = open(options.files[0], mode='r')
   g1_files = g1_inf.readlines()
except IOError:
   print "Error opening input file: %s" % options.files[0]
   sys.exit()

try:
   g2_inf = open(options.files[1], mode='r')
   g2_files = g2_inf.readlines()
except IOError:
   print "Error opening input file: %s" % options.files[1]
   sys.exit()

try:
   results_outf = open("results.csv", mode='w')
   results_outf.write("k,%s_unique,%s_unique,common\n" % (options.groups[0], options.groups[1]))
except IOError:
   print "Error opening output file results.csv"
   sys.exit()
   
try:
   log_outf = open("log.txt", mode='w')
except IOError:
   print "Error opening output log.txt"
   sys.exit()

def log(msg):
   t = time.strftime("%m-%d-%Y %H:%M:%S", time.localtime())
   txt = "%s: %s" % (t, msg)
   print(txt)
   log_outf.write(txt+"\n")
   log_outf.flush() #flush file buffer
   #os.fsync() #and force writing to disk (so that tail -f works)
   
#global variables:
startt = time.time()
   
for k in range(options.start, options.end+1):
   uk_g1 = {} #stores counts for group1 kmers   [kmer] total, forward, revcomp, copynum
   uk_g2 = {} #stores counts for group2 kmers
   ck = {} #stores common kmers 
   reads = 0 #total reads processed 
   log("Starting k-mer computation for k=%s" % k)
   log("Processing files from group %s" % options.groups[0])
   for file in g1_files:
      file = file.strip()
      log("Processing %s" % file)
      try:
         inf = open(file, mode='r')
      except IOError:
         print "Error opening input: %s, exiting." % file
         sys.exit()
       
      freads = 0  #reads in file
      ftime = time.time()
      fkmers = 0
      for line in inf:
         if(line[0] != ">"):
            line = line.strip().upper()
            rcounts = {}
            reads += 1; freads +=1
            if(reads % 100000 == 0):
               txt = "File:%s, read: %s, total reads:%s, reads/sec:%s" %(file, freads, reads, reads/(time.time() - ftime) )
               log(txt)
            i = 0
            #for i in range(len(line) - k +1): #add one to ensure analysis of full line
            while i + k <= len(line):
               #check global counts
               kmer = line[i:i+k]
               if( 'N' in kmer):
                  i = i + kmer.rfind("N") + 1
               else:
                  i=i+1
                  #[kmer] tot,f, rc, copynum
                  kmer_rc = kmer.translate(string.maketrans('ATCG','TAGC'))[::-1]#reverse complement 
                  if kmer in uk_g1:
                     uk_g1[kmer][0] += 1
                     uk_g1[kmer][1] += 1
                     if(kmer not in rcounts):
                        uk_g1[kmer][3] += 1
                        rcounts[kmer] = True
                  elif(kmer_rc in uk_g1):
                     uk_g1[kmer_rc][0] += 1
                     uk_g1[kmer_rc][2] += 1
                     if(kmer_rc not in rcounts):
                        uk_g1[kmer_rc][3] += 1
                        rcounts[kmer_rc] = True
                  else: #[kmer] tot,f, rc, copynum
                     uk_g1[kmer] = [1, 1, 0, 1]
                     rcounts[kmer] = True
                     fkmers += 1
      inf.close()
      log("Finished processing %s, total lines: %s, total new kmers: %s" % (file, freads,  fkmers))
      log("Reads per second for %s: %s, total reads per second %s" % (file, freads/(time.time() - ftime), reads/(time.time() - startt)))
   log("Finished processing group %s, \n\ttotal reads: %s, \n\ttotal kmers: %s \n\treads per second:%s" % (options.groups[0], reads, len(uk_g1), reads/(time.time() - startt)))
   log("Processing files from group %s" % options.groups[1])
   for file in g2_files:
      file = file.strip()
      log("Processing %s" % file)
      try:
         inf = open(file.strip(), mode='r')
      except IOError:
         print "Error opening input: %s, exiting." % file
         sys.exit()
       
      freads = 0  #reads in file
      ftime = time.time()
      fkmers = 0
      for line in inf:
         if(line[0] != ">"):
            line = line.strip().upper()
            rcounts = {}
            reads += 1; freads +=1
            if(reads % 100000 == 0):
               txt = "File:%s, read: %s, total reads:%s, reads/sec:%s" % (file, freads, reads, reads/(time.time() - ftime))
               txt = txt + ("\n  %s unique: %s, %s unique: %s, common: %s" % (options.groups[0], len(uk_g1), options.groups[1], len(uk_g2), len(ck)))
               log(txt)
            i=0
            while i + k <= len(line):
               #check global counts
               kmer = line[i:i+k]
               if( 'N' in kmer):
                  i = i + kmer.rfind("N") + 1
               else:
                  i=i+1
                  #[kmer], 0:g1_total, 1:g1_f, 2:g1_rc, 3:g1_copynum,  4:g2_total, 5:g2_f, 6:g2_rc, 7:g2_copynum
                  kmer_rc = kmer.translate(string.maketrans('ATCG','TAGC'))[::-1]#reverse complement 
                  if kmer in ck:
                     ck[kmer][4] += 1
                     ck[kmer][5] += 1
                     if kmer not in rcounts:
                        ck[kmer][7] += 1
                        rcounts[kmer] = True
                  elif kmer_rc in ck:
                     ck[kmer_rc][4] += 1
                     ck[kmer_rc][6] += 1
                     if kmer_rc not in rcounts:
                        ck[kmer_rc][7] += 1
                        rcounts[kmer_rc] = True                     
                  elif kmer in uk_g1: #should be in ck
                     ck[kmer] = [uk_g1[kmer][0], uk_g1[kmer][1], uk_g1[kmer][2], uk_g1[kmer][3], 1, 1, 0, 1]
                     del uk_g1[kmer]
                     rcounts[kmer] = True
                  elif(kmer_rc in uk_g1): #should be in ck, but swap f and rc for g2
                     ck[kmer_rc] = [uk_g1[kmer_rc][0], uk_g1[kmer_rc][1], uk_g1[kmer_rc][2], uk_g1[kmer_rc][3], 1, 0, 1, 1]
                     del uk_g1[kmer_rc]
                     rcounts[kmer_rc] = True
                  elif(kmer in uk_g2): #kmer is already in g2 : [kmer] 0:tot, 1:f,  2:rc, 3:copynum
                     uk_g2[kmer][0] += 1
                     uk_g2[kmer][1] += 1
                     if(kmer not in rcounts):
                        uk_g2[kmer][3] += 1
                        rcounts[kmer] = True
                  elif(kmer_rc in uk_g2):
                     uk_g2[kmer_rc][0] += 1
                     uk_g2[kmer_rc][2] += 1
                     if(kmer_rc not in rcounts):
                        uk_g2[kmer_rc][3] += 1
                        rcounts[kmer_rc] = True
                  else: # new kmer : [kmer] tot,f, rc, copynum
                     uk_g2[kmer] = [1, 1, 0, 1]
                     rcounts[kmer] = True
                     fkmers += 1
      
      inf.close()
      log("Finished processing %s, total lines: %s, total new kmers: %s" % (file, freads,  fkmers))
      log("Reads per second for %s: %s, total reads per second %s" % (file, freads/(time.time() - ftime), reads/(time.time() - startt)))
      log("%s unique: %s, %s unique: %s, common: %s" % (options.groups[0], len(uk_g1), options.groups[1], len(uk_g2), len(ck)))
      
   log("Finished processing files for k=%s" % k)
   log("%s unique: %s, %s unique: %s, common: %s" % (options.groups[0], len(uk_g1), options.groups[1], len(uk_g2), len(ck)))
   results_outf.write("%s,%s,%s,%s\n" % (k,len(uk_g1), len(uk_g2), len(ck)))
   if(options.verbose):
      log("Vebose mode set, writing detailed results files for k=%s" % k)
      # ####### group 1 #########
      try:
         fname = "k%s_%s_unique.csv" % (k, options.groups[0])
         outf1 = open(fname, mode='w')
      except IOError:
         log(("Error opening output file %s" % fname))
         sys.exit()#write detailed output files
      l=[]
      log("Sorting %s results...." % options.groups[0])
      for kmer in uk_g1.keys():
         l.append([kmer, float(uk_g1[kmer][3])])

      l = sorted(l,  key=lambda counts: counts[1],  reverse=True)

      #Write out results
      log("Writing %s results....." % options.groups[0])
      outf1.write("kmer,total,forward,revcomp,copynum\n")
      for r in l:
         kmer = r[0]
         d = uk_g1[kmer]
         outf1.write("%s,%s,%s,%s,%s\n" %(kmer,d[0], d[1], d[2], d[3]))

      l = []
      uk_g1 = {}
      outf1.close()
      # ####### group 2 #########
      try:
         fname = "k%s_%s_unique.csv" % (k, options.groups[1])
         outf1 = open(fname, mode='w')
      except IOError:
         log(("Error opening output file %s" % fname))
         sys.exit()#write detailed output files
      l=[]
      log("Sorting %s results...." % options.groups[1])
      for kmer in uk_g2.keys():
         l.append([kmer, float(uk_g2[kmer][3])])

      l = sorted(l,  key=lambda counts: counts[1],  reverse=True)

      #Write out results
      log("Writing %s results....." % options.groups[1])
      outf1.write("kmer,total,forward,revcomp,copynum\n")
      for r in l:
         kmer = r[0]
         d = uk_g2[kmer]
         outf1.write("%s,%s,%s,%s,%s\n" %(kmer,d[0], d[1], d[2], d[3]))

      l=[]
      uk_g2 = {}
      outf1.close()
      
      # ####### ck #########
      try:
         fname = "k%s_common.csv" % (k)
         outf1 = open(fname, mode='w')
      except IOError:
         log(("Error opening output file %s" % fname))
         sys.exit()#write detailed output files
      l=[]
      log("Sorting common results....")
      # Write out ck: #kmer,i1_occurrences,i1_forward,i1_revcomp,i1_copynum,  i2_occurrences,i2_forward,i2_revcomp,i2_copynum
      l = []
      #Generate a list of sorted keys
      for kmer in ck.keys():
         l.append([kmer, (float(ck[kmer][3]) + float(ck[kmer][7]))/2 ] )

      l = sorted(l,  key=lambda counts: counts[1],  reverse=True)

      #Write out results
      log("Writing common results.....")
      outf1.write("kmer,%s_total,%s_forward,%s_revcomp,%s_copynum,%s_occurrences,%s_forward,%s_revcomp,%s_copynum\n" % (
         options.groups[0],  options.groups[0], options.groups[0], options.groups[0], options.groups[1], options.groups[1], options.groups[1], options.groups[1]))
      for r in l:
         kmer = r[0]
         d = ck[kmer]
         outf1.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(kmer,d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]))

      l = []
      ck = {}
      outf1.close()
      
   log("----------------------------------------\n\n")







   
# Close all files
g1_inf.close()
g2_inf.close()
results_outf.close()
log_outf.close()












"""

alp = ['A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','N']
sequences = []
for  _ in xrange(1000000):
   sequences.append("".join(random.choice(alp) for _ in xrange(random.uniform(20,500))))


t = time.time()
n = 0
for s in sequences:
   if "N" in s:
      n +=1

print time.time() - t


t = time.time()
n = 0
for s in sequences:
   if s.find("N") != -1:
      n +=1

print time.time() - t
"""
