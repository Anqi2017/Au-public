#!/usr/bin/python
import argparse, re, sys

def main():
  parser = argparse.ArgumentParser(description="Concatonate files that contain an integer (and only one integer) in their name in that integer order.")
  parser.add_argument('file_names',nargs='+',help="File names containing integers that are their order")
  args = parser.parse_args()
  nums = {}
  for fname in args.file_names:
    m = re.search('(\d+)[^\d]*$',fname)
    if not m:
      sys.stderr.write("WARNING contains nonnumerical file name "+fname+"\n")
      sys.exit()
    num = int(m.group(1))
    if num in nums:
      sys.stderr.write("ERROR same number twice:\n"+nums[num]+"\n"+fname+"\n")
      sys.exit()
    nums[num] = fname
  for num in sorted(nums.keys()):
    with open(nums[num]) as inf:
      for line in inf:  print line.rstrip()

if __name__=="__main__":
  main()
