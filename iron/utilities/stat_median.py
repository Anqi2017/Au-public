#!/usr/bin/python
import sys,re
from Bio.Statistics import median
values = []
for line in sys.stdin:
  m = re.match('(\S+)',line)
  values.append(float(m.group(1)))
print median(values)
