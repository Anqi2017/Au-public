import sys
from math import log10, pow

def geometric_mean(arr):
  if len(arr) == 0:
    sys.stderr.write("ERROR: no content in array to take average\n")
    sys.exit()
  if len(arr) == 1: return arr[0]
  return pow(10,float(sum([log10(float(x)) for x in arr]))/float(len(arr)))

def average(arr):
  if len(arr) == 0:
    sys.stderr.write("ERROR: no content in array to take average\n")
    sys.exit()
  if len(arr) == 1:  return arr[0]
  return float(sum(arr))/float(len(arr))

def median(arr):
  if len(arr) == 0:
    sys.stderr.write("ERROR: no content in array to take average\n")
    sys.exit()
  if len(arr) == 1: return arr[0]
  quot = len(arr)/2
  rem = len(arr)%2
  if rem != 0:
    return sorted(arr)[quot]
  return float(sum(sorted(arr)[quot-1:quot+1]))/float(2)

