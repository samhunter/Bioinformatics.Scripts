#this script reads in a lucy result file and produces a set of output files:
# lucy_clips.csv
# lucy_clipped.fasta 
#>GKF1FGD02EGOM3 0 0 0 1 79
#gactacacgtagtatCTCAATGATTCATACAGCTTCACATCGAGAGCAT

inf = open("lucy.fasta", mode='r')

lclips =  open("lucy_clips.csv", mode='w')
lclips.write("ID,lclip,rclip\n")

lclipped =  open("lucy_clipped.fasta", mode='w')

id = "" #read ID
seq = "" #sequence
lclip = 0
rclip = 0
errors = 0
counter = 0
min_lclip=0

#The outer for loop iterates through all lines in the input file
for line in inf:
  if(line[0] == '>'):
    #first if is a special case for the first line
    if(id == ""):
      pieces = line.strip().split(" ")
      id = pieces[0]
      lclip = int(pieces[4])
      if(lclip < min_lclip): lclip = 16 #first 16 are tag and 454 bases
      rclip = int(pieces[5])
      if(id == ''):
        errors += 1
    else:
      lclips.write("%s,%s,%s\n" % (id[1:],lclip,rclip))
      lclipped.write("%s\n%s\n" %(id,seq[lclip-1:rclip-1]))
      counter += 1
      if(counter % 1000 == 0):
        print("Total records %s" % (counter))
      seq = ""
      pieces = line.strip().split(" ")
      id = pieces[0]
      lclip = int(pieces[4])
      if(lclip < 16): lclip = 16
      rclip = int(pieces[5])
      if(id == ''):
        errors += 1  
  else:
    seq = line.strip("\n")  


#take care of last line:
seq = line.strip("\n")

#write last record
lclips.write("%s,%s,%s\n" % (id[1:],lclip,rclip))
lclipped.write("%s\n%s\n" %(id,seq[lclip-1:rclip-1]))


lclips.close()   
lclipped.close()

print("Total records %s.  Total errors %s" %(counter, errors))
