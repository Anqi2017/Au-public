#!/usr/bin/python
import argparse, sys

# Import our main launchers for each mode
import analyze

version = 0.91

def main():
  operable_argv = [sys.argv[0]]+sys.argv[2:]
  sys.argv = sys.argv[:2]  
  #do our inputs
  args = do_inputs()
  if args.mode == 'analyze':
    analyze.external_cmd(" ".join(operable_argv),version=version)
  else:
    sys.stderr.write("Run mode not yet implemented\n")

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('mode',choices=['analyze','compare','combine'],help="MODE of program to run")
  args = parser.parse_args()
  return args


if __name__=="__main__":
  main()
