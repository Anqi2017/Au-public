#!/usr/bin/python
import sys, argparse

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  args = parser.parse_args()
  
  if args.input == '-':
    args.input = sys.stdin
  else: args.input = open(args.input)


if __name__=="__main__":
  main()
