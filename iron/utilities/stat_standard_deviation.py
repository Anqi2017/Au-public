#!/usr/bin/python
import sys,re
from Bio.Statistics import standard_deviation
values = []
for line in sys.stdin:
  m = re.match('(\S+)',line)
  values.append(float(m.group(1)))
print standard_deviation(values)
