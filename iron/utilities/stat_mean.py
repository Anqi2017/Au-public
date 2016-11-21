#!/usr/bin/python
import sys,re
from Bio.Statistics import average
values = []
for line in sys.stdin:
  m = re.match('(\S+)',line)
  values.append(float(m.group(1)))
print average(values)
