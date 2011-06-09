def read_fasta_to_dictionary(filename=""):
   try:
      fasta_inf = open(filename, mode='r')
   except IOError:
      print "Error opening input file: %s\n" % filename
      sys.exit()
      
   

def reverse_complement(seq):
   from string import maketrans
   return(seq.translate(maketrans('ATCGatcg', 'TAGCtagc')))[::-1]


