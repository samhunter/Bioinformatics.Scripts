"""
  Given a big set of sequences, and a bunch of restriction enzyme cut site
  patterns, search the sequences for the number of occurances of cut-sites.
"""
import re

res = open("link_cutzymes", mode='r')

seq = open("/bio/Working/2010-07-23-Chinook-Salmon-Onch.tsh/filtering/lucy_clipped.fasta", mode='r')

patterns = {}
#Parse restriction patterns file
p = re.compile('.*[^AGCT]+.*') #match anything that isn't AGCT
for line in res:
  if(line[0] == 'C'):
    l = line.split()
    if not p.match(l[2]):
      patterns[l[1]] = l[2]
      
#Parse fasta
seqs = ""
for line in seq:
  if(line[0] != '>'):
    seqs += line.upper()
   
#for each pattern, count occurances
counts = {}
i=0
for k in patterns.keys():
  p = re.compile('.%s.'% patterns[k])
  counts[k] = len(p.findall(seqs))
  i+=1
  print("%s %s" % (k, i))
  
# Calculate frequencies and write them out
nchars = len(seqs)
outf = open("RE_freqs.csv", mode='w')
outf.write("RE,count,perbase\n")
for k in counts.keys():
  outf.write("%s,%s,%s,%s\n" % (k, patterns[k], counts[k], nchars/counts[k]))

outf.close()
  