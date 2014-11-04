#!/usr/bin/python
import sys, re

#pre: 1 the filename outputs from LSC Full reads (such as full_LR.fa), 
#     2 filename of corrected reads (such as corrected_LR.fa), 
#     3 a threshold (i.e. 0.90) for the fraction of length that corrected must be of full
#       in order to use the full rather than the corrected
#     4 an output fasta file where corrected reads have been switched in
#     5 an output of file to list the swapped out read names to
#post: writes to the output fasta to the file, and writes the list of swapped genes
#      it leaves read names the same as they were in the file they came from
#modifies: fileIO

def main():
  if len(sys.argv) != 6:
    print sys.argv[0] + ' <full_LR.fa (input)> <corrected_LR.fa (input)> <threshold (i.e. 0.9)> <output fasta> <output list>'
    sys.exit()
  full_filename = sys.argv[1]
  corrected_filename = sys.argv[2]
  threshold = float(sys.argv[3])
  output_fasta = sys.argv[4]
  output_list = sys.argv[5]

  sys.stderr.write("reading through full_LR to get lengths and names\n")
  fulllens = read_lengths(full_filename)
  sys.stderr.write("read "+str(len(fulllens))+" full_LR\n")

  sys.stderr.write("reading through corrected_LR length to get lengths and names\n")
  correctlens = read_lengths(corrected_filename)
  sys.stderr.write("read "+str(len(correctlens))+" correct_LR\n")

  sys.stderr.write("selecting reads to keep with threshold "+str(threshold)+"\n")
  keepers = set()
  for readname in correctlens:
    if readname in fulllens:
      if fulllens[readname] >= correctlens[readname] and fulllens[readname]*threshold <= correctlens[readname]:
        keepers.add(readname)
  sys.stderr.write("found "+str(len(keepers))+" sequences to swap\n")

  sys.stderr.write("reading through full_LR to get sequences and names\n")
  keeperseqs = get_keeperseqs(full_filename,keepers)
  sys.stderr.write("read "+str(len(keeperseqs))+" sequences to keep\n")

  sys.stderr.write("reading through corrected_LR and writing outputs\n")
  print_seqs(corrected_filename,keeperseqs,fulllens,correctlens,output_fasta,output_list)


def print_seqs(fasta,keeperseqs,fulllens,correctlens,outfasta,outlist):
  p = re.compile('^>([^\|]+)')
  ofasta = open(outfasta,'w')
  olist = open(outlist,'w')
  with open(fasta) as infile:
    name = ''
    header = ''
    keeper = 0
    seq = ''
    for line in infile:
      line = line.rstrip()
      m = p.match(line)
      if m:
        if seq != '' and keeper == 0: # if its not the first one and we haven't already printed it
          ofasta.write(header+"\n")
          ofasta.write(seq+"\n")
        header = line
        name = m.group(1)
        keeper = 0
        seq = '';
        if name in keeperseqs:
          keeper = 1
          ofasta.write(">"+name+"\n")
          ofasta.write(keeperseqs[name]+"\n")
          flen = fulllens[name]
          clen = correctlens[name]
          olist.write(name+"\t"+str(flen)+"\t"+str(clen)+"\t"+str(float(clen)/float(flen))+"\n")
      else:
        if keeper == 0:
          seq+=line
    if keeper == 0 and seq != '':
          ofasta.write(header+"\n")
          ofasta.write(seq+"\n")      
  ofasta.close()
  olist.close()

def get_keeperseqs(fasta,keepers):
  p = re.compile('^>([^\|]+)')
  seqs = {}
  with open(fasta) as infile:
    name = ''
    keeper = 0
    for line in infile:
      line = line.rstrip()
      m = p.match(line)
      if m:
        name = m.group(1)
        keeper = 0
        if name in keepers:
          keeper = 1
          seqs[name] = '';
      else:
        if keeper == 1:
          seqs[name]+=line
  return seqs

def read_lengths(fasta):
  lens = {}
  p = re.compile('^>([^\|]+)')
  with open(fasta) as infile:
    name = ''
    for line in infile:
      line = line.rstrip()
      m = p.match(line)
      if m:
        name = m.group(1)
        lens[name] = 0;
      else:
        lens[name]+=len(line)
  return lens

main()
