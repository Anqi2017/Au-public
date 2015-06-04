#!/usr/bin/python
import sys
fi = sys.argv[1]
colcount = 0
with open(fi) as inf:
  for line in inf:
    f = line.rstrip().split("\t")
    if len(f) > colcount: colcount = len(f)
cols = []
with open(fi) as inf:
  for line in inf:
    col = line.rstrip().split("\t")
    for i in range(len(col),colcount):
      col.append("")
    cols.append(col)
for i in range(0,colcount):
  v = ""
  for j in range(0,len(cols)):
    v += cols[j][i] + "\t"      
  print v.rstrip("\t")
