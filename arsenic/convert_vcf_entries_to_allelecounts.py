#!/usr/bin/python
import sys, argparse
#from __future__ import print_function


def main():
  for line in open("bcftools.merged.data.GJB4.cols10-205"):
    columns = line.split("\t")
    for c in columns:
      nums=[c[0], c[2]]
#      print(nums)
      if (nums[0] == "."):
        nums[0] = 0
      elif (nums[0] == "0"):
	nums[0] = 0
      else: 
        nums[0] = 1
      if (nums[1] == "."):
        nums[1] = 0
      elif (nums[1] == "0"):
	nums[1] = 0
      else:
        nums[1] = 1

#      print(nums[0])
      intnums=[int(nums[0]),int(nums[1])]
      print(sum(intnums)),
#      if (sum(intnums) > 0):
#        print(1),
#      elif (sum(intnums) == 0):
#        print(0),
#      else:
#        print("ERROR"),
      print("\t"),
    print
main()




