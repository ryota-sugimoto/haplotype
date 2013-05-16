#!/usr/bin/env python

#Copyright (c) 2013 All Right Reserved.
#
#author: Ryota Sugimoto
#institute: National Institute of Genetics
#email: ryota.sugimoto@gmail.com
#date: 2013-05-16

import re

def genotype(line):
  return line.split("\t")[9].split(":")[0]

def refReadNum(line):
  return float(line.split("\t")[9].split(":")[1].split(",")[0])

def altReadNum(line):
  return float(line.split("\t")[9].split(":")[1].split(",")[1])

def isAA(line, threshold):
  return altReadNum(line) / (refReadNum(line) + altReadNum(line)) < threshold

def isBB(line, threshold):
  return altReadNum(line) / (refReadNum(line) + altReadNum(line)) >= threshold

def rewrite2BB(line, bb_threshold):
  if re.match("^#", line):
    return line
  if isBB(line, bb_threshold):
    splited = line.split("\t")
    gen = splited[9].split(":")
    gen[0] = "1/1"
    splited[9] = ":".join(gen)
    return "\t".join(splited)
  else:
    return line

def rewriteVcf(file, aa_threshold, bb_threshold):
  file.seek(0)
  header = filter(lambda line: re.match("^#", line), file)
  file.seek(0)
  data = filter(lambda line: not re.match("^#", line), file)
  
  pf = lambda line: not isAA(line, aa_threshold)
  rewrited_data = map(lambda line: rewrite2BB(line, bb_threshold),
                      filter(pf,data))
  return "\n".join(map(lambda line: line.strip(), header + rewrited_data))
 
def parse():
  import argparse
  
  parser = argparse.ArgumentParser()
  parser.add_argument("vcf_file")
  parser.add_argument("-a", "--aa_threshold",
                      default=0.3,
                      type=float)
  parser.add_argument("-b", "--bb_threshold",
                      default=0.8,
                      type=float)
  return parser.parse_args()

if __name__ == "__main__":
  args = parse()
  vcf = open(args.vcf_file)
  print rewriteVcf(vcf, args.aa_threshold, args.bb_threshold)
