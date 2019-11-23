#!/usr/bin/python
"""
#this script reads in a lucy result file plus a GLFLX xml file and produces a set of output files:
<base>_final_clips.csv
<base>_clipped.fasta

   fasta input looks like:
      >GKF1FGD02EGOM3 0 0 0 1 79
      gactacacgtagtatCTCAATGATTCATACAGCTTCACATCGAGAGCAT
"""

import time
import sys
import string
from optparse import OptionParser #http://docs.python.org/library/optparse.html
from Bio import pairwise2
import thread

starttime = time.time()
parser = OptionParser()
parser.add_option('-b',  '--base')

(options,  args) = parser.parse_args() #uncomment this line for command line support
#(options,  args) = parser.parse_args([ '-b', 'F58WHOC01_M7_Anc'])
if(options.base is None):
   parser.print_help()
   sys.exit()

base = options.base
verbose = True
#base = "F58WHOC01_M7_Anc"


#input files
try:
   fn = base+"_lucy.fasta"
   fasta_inf = open(fn, mode='r')
except IOError:
   print "Error opening input file: %s\n" % fn
   sys.exit()

try:
   fn = base+".xml"
   xml_inf = open(fn, mode='r')
except IOError:
   print "Error opening input file: %s\n" % fn
   sys.exit()

#output files
try:
   fn = base+"_errors.txt"
   errors_outf = open(fn, mode='w')
except IOError:
   print "Error opening output file: %s\n" % fn
   sys.exit()

try:
   fn = base+"_final_clips.txt"
   final_outf = open(fn, mode='w')
except IOError:
   print "Error opening output file: %s\n" % fn
   sys.exit()

try:
   fn = base+"_clipped.fasta"
   clipped_outf = open(fn, mode='w')
except IOError:
   print "Error opening output file: %s\n" % fn
   sys.exit()

try:
   fn = base+"_readlengths.csv"
   readl_outf = open(fn, mode='w')
   readl_outf.write("#rid,unclipped,roche,my,lucy,final\n")   
except IOError:
   print "Error opening output file: %s\n" % fn
   sys.exit()
   
try:
   fn = base+"_keepreads.txt"
   keep_outf = open(fn, mode='w')
except IOError:
   print "Error opening output file: %s\n" % fn
   sys.exit()   

def find_lrclip(seq,  key = "TCAG",  mids = ["ACGAGTGCGTAGATGTGTATAAGAGACAG", "ACGCTCGACAAGATGTGTATAAGAGACAG"],  B="CTGTCTCTTATACACATCTCTGA"):   
   lclip = 0
   rclip = len(seq)
   seq = seq.upper()
   
   #Look for tag/adaptor, typically these are from <tag><mid><tag><mid> situations that get missed by Roche 
   for m in mids:
      f1 = seq.rfind(key+m)
      if f1 != -1:
         lclip = f1 + len(key + m)
         break
         
   if lclip == 0:
      for m in mids:
         aln = pairwise2.align.localms(seq[0:200], (key+m), 1,-1,-1,-.5)
         for a in aln:
            if(a[2] > 12 and lclip < a[4]):
               lclip = a[4]
               if verbose:
                  print("Aligned better lclip.")
                  print(a[0][0:a[4]+10])
                  print(a[1][0:a[4]+10])
         if lclip>0: break #if we found a good lclip, stop
   
   #Attempt to find an R clip in the last 100 bp
   f = seq.find(B[0:13])
   if f != -1:
      rclip = f
   else:  
      aln = pairwise2.align.localms(seq[-100:], (B[0:18]), 1,-1,-1,-.5)
      for a in aln:
         if(a[2]>15 and rclip > len(seq) -100 + a[3]):
            rclip = len(seq) - 100 + a[3]
            if verbose:
               print("Aligned better rclip")
               print(a[0][a[3]-10:])
               print(a[1][a[3]-10:])

   return(lclip,  rclip)

def close_files():   
   fasta_inf.close()
   xml_inf.close()
   errors_outf.close()
   final_outf.close()
   clipped_outf.close()
   readl_outf.close()
   keep_outf.close()

def do_output(rid,  seq, lucy_lclip,  lucy_rclip):
   global lucylclip_counter; global lucyrclip_counter; global mylclip_counter; global myrclip_counter; global counter; global discarded_counter
   global errors;
   xml_lclip = 0; xml_rclip = len(seq)
   if rid in xml_clips:
      if xml_clips[rid][0] > xml_lclip:xml_lclip = xml_clips[rid][0]
      if xml_clips[rid][1] < xml_rclip:xml_rclip = xml_clips[rid][1]
   (my_lclip, my_rclip) = find_lrclip(seq)
   #grab the most conservative clip points
   #left clips
   lclip = max(lucy_lclip,  xml_lclip,  my_lclip)
   if(lclip == lucy_lclip and lucy_lclip > xml_lclip): lucylclip_counter += 1
   if(lclip == my_lclip and my_lclip > xml_lclip):  mylclip_counter+= 1
   #right clips, 
   if(lucy_rclip is None): lucy_rclip = len(seq)
   if(xml_rclip is None): xml_rclip = len(seq)
   if(my_rclip is None): my_rclip = len(seq)
   rclip = min(lucy_rclip, xml_rclip,  my_rclip)
   if(rclip == lucy_rclip and lucy_rclip < xml_rclip): lucyrclip_counter += 1
   if(rclip == my_rclip and my_rclip < xml_rclip): myrclip_counter += 1
   #write new clippoints
   if(rclip - lclip > 25):
      final_outf.write("%s\t%s\t%s\n" % (rid,lclip,rclip))
      clipped_outf.write(">%s\n%s\n" %(rid,seq[lclip-1:rclip-1]))
      #rid, unclipped, roche, my, lucy, final
      readl_outf.write("%s,%s,%s,%s,%s,%s\n" %(rid, len(seq),  xml_rclip-xml_lclip,  my_rclip-my_lclip,  lucy_rclip-lucy_lclip,  rclip-lclip))
      keep_outf.write("%s\n" % rid)
      counter += 1
      if(counter % 1000 == 0):
         print("Records processed: %s" % (counter))
         print("\tlucy lclips used: %s" % lucylclip_counter)
         print("\tlucy rclips used: %s" % lucyrclip_counter)
         print("\tmy lclips found: %s" % mylclip_counter)
         print("\tmy rclips found: %s" % myrclip_counter)
         print("\tDiscarded reads: %s" % discarded_counter)
   else:
      discarded_counter += 1
      print("Discard read %s, lclip=%s, rclip=%s, read length=%s" % (rid,  lclip,  rclip,  rclip-lclip))
      print("\txml_lclip:%s, xml_rclip:%s, lucy_lclip:%s, lucy_rclip:%s, my_lclip:%s, my_rclip:%s"% (xml_lclip,xml_rclip, lucy_lclip,lucy_rclip,my_lclip,  my_rclip) )
      errors_outf.write("Discard read %s, lclip=%s, rclip=%s, read length=%s\n" % (rid,  lclip,  rclip,  rclip-lclip))
      errors_outf.write("\txml_lclip:%s, xml_rclip:%s, lucy_lclip:%s, lucy_rclip:%s, my_lclip:%s, my_rclip:%s\n"% (xml_lclip,xml_rclip, lucy_lclip,lucy_rclip,my_lclip,  my_rclip))
      errors+=1

#Parse XML, I could use the python xml capabilities.. but they look annoying
xml_clips = {}
rid = None
lclip = None
rclip = None


print("Parsing XML file.\n")
for l in xml_inf:
   l = l.strip()
   if(l[0:12] == "<trace_name>"):
      if(rid is not None or lclip is not None or rclip is not None):
         errors_outf.write("Error reading XML, Read: %s was missing data.\n\tlclip=%s, rclip=%s" % (rid,  lclip,  rclip))
         print("Error reading XML, Read: %s was missing data.\n\tlclip=%s, rclip=%s" % (rid,  lclip,  rclip))
         xml_clips[rid] = [lclip,  rclip]
         lclip = None
         rclip = None
         rid = l[12:l.find("</trace_name>")]
      else:
         rid = l[12:l.find("</trace_name>")]
   if(l[0:18] == "<clip_vector_left>"):
      if(rid is None or lclip is not None):
         errors_outf.write("Error reading XML, Read: %s was missing data.\n\tlclip=%s, rclip=%s" % (rid,  lclip,  rclip))
         print("Error reading XML, Read: %s was missing data.\n\tlclip=%s, rclip=%s" % (rid,  lclip,  rclip))
         xml_clips[rid] = [lclip,  rclip]
         lclip = None
         rclip = None
         rid = l[12:l.find("</trace_name>")]
      else:
         lclip = int(l[18:l.find("</clip_vector_left>")])
   if(l[0:19] == "<clip_vector_right>"):
      if(rid is None or lclip is None or rclip is not None):
         errors_outf.write("Read: %s, XML contained weird right clip: %s" % (rid,  l))
         rclip = None
      else:
         rclip = int(l[19:l.find("</clip_vector_right>")])
   if(rid is not None and lclip is not None and rclip is not None):   
      xml_clips[rid] = [lclip,  rclip]
      rid = None
      lclip = None
      rclip = None

print("XML file parsed, %s total records found." % len(xml_clips.keys()))




#Parse fasta file
rid = None #read ID
seq = '' #sequence
lclip = None
rclip = None
errors = 0
counter = 0
min_lclip=0
mylclip_counter = 0
lucylclip_counter = 0
lucyrclip_counter = 0
myrclip_counter = 0
discarded_counter = 0

#The outer for loop iterates through all lines in the input file
print("\nParsing fasta file.")
for line in fasta_inf:
   if(line[0] == '>'):
      #first if is a special case for the first line
      if(seq == ''):
         pieces = line.strip().split(" ")
         rid = pieces[0][1:]
         lucy_lclip = int(pieces[4])
         lucy_rclip = int(pieces[5])
         if(rid == ''):
            errors += 1
      else:
         do_output(rid,  seq, lucy_lclip,  lucy_rclip)
         seq = ''; rid=None; lclip = None;  rclip = None;  xml_lclip = None; xml_rclip=None; lucy_lclip=None;  lucy_rclip=None
         
         pieces = line.strip().split(" ")
         rid = pieces[0][1:]
         lucy_lclip = int(pieces[4])
         lucy_rclip = int(pieces[5])
         if(id == ''):
            errors += 1  
   else:
       seq = seq + line.strip("\n")  


#take care of last line:
seq = line.strip("\n")
do_output(rid,  seq,  lucy_lclip,  lucy_rclip)

close_files()

endtime = time.time()

print("Total records processed %s.  Total errors identified %s" %(counter, errors))
print("Total lucy lclips used: %s" % lucylclip_counter)
print("Total lucy rclips used: %s" % lucyrclip_counter)
print("Total improved lclips found: %s" % mylclip_counter)
print("Total improved rclips found: %s" % myrclip_counter)
print("Total Discarded reads: %s" % discarded_counter)
print("Start Time: %s, End time: %s, total:%s, reads/time:%s" % (starttime,  endtime,  endtime-starttime, counter/(endtime-starttime) ))
