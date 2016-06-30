#!/usr/bin/python
import sys, argparse, os, inspect

#bring in the folder to the path for our utilities
pythonfolder_loc = "../utilities"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

import build_emitter
import emit
import train_errors

version = 0.9

def main():
  operable_argv = [sys.argv[0]]+sys.argv[2:]
  sys.argv = sys.argv[0:2]
  args = do_inputs()
  
  if args.mode == 'emit':
    #sys.stderr.write(" ".join(operable_argv)+"\n")
    emit.external_cmd(" ".join(operable_argv),version=version)
  elif args.mode == 'build_emitter':
    build_emitter.external_cmd(" ".join(operable_argv),version=version)
  if args.mode == 'train_errors':
    train_errors.external_cmd(" ".join(operable_argv),version=version)

def do_inputs():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('mode',choices=['emit','build_emitter','train_errors'],help="Run mode for simulator")
  args = parser.parse_args()
  return args

if __name__=="__main__":
  main()
