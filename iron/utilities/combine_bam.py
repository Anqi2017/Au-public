#!/usr/bin/python
import sys, argparse, re, subprocess

def main():
  parser = argparse.ArgumentParser(description = 'Combine sam files')
  parser.add_argument('bam_files',nargs='+',help='FILENAME for sam files')
  args = parser.parse_args()
  header = False
  seqs = set()
  tagorder = []
  tagseen = {}
  for file in args.bam_files:
    cmd = "samtools view -H "+ file
    process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
    while True:
        line = process.stdout.readline()
        if not line: break
        line = line.rstrip()
        f = line.split("\t")
        m = re.match('^(@\S\S)\s',line)
        if not m or len(f) > 10: break
        if m.group(1) == '@SQ':
          seqs.add(line)
        if m.group(1) not in tagseen:
          tagorder.append(m.group(1))
        tagseen[m.group(1)] = line
  for tag in tagorder:
    if tag != '@SQ':
      print tagseen[tag]
    else:
      for seq in sorted(seqs):
        print seq      
  #now go back through and do the sam data
  for file in args.bam_files:
    cmd = "samtools view "+file
    process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
    while True:
        line = process.stdout.readline()
        if not line: break
        line = line.rstrip()
        f = line.split("\t")
        if len(f) > 10:
          print line
main()
