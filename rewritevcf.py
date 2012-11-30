#!/usr/bin/env python
import re

def genotype(line):
  return line.split("\t")[9].split(":")[0]

def refReadNum(line):
  return int(line.split("\t")[9].split(":")[1].split(",")[0])

def altReadNum(line):
  return int(line.split("\t")[9].split(":")[1].split(",")[1])

def isAA(line):
  return refReadNum(line) > altReadNum(line)

def isBB(line):
  return refReadNum(line) < altReadNum(line)

def rewrite2BB(line):
  if re.match("^#", line):
    return line
  else:
    splited = line.split("\t")
    gen = splited[9].split(":")
    gen[0] = "1/1"
    splited[9] = ":".join(gen)
    return "\t".join(splited)

def rewriteVcf(file):
  return "\n".join(map(lambda line: line.strip(),
                     map(rewrite2BB, 
                         filter(lambda line: re.match("^#",line) or isBB(line),
                                file))))
 
def parse():
  import argparse
  
  parser = argparse.ArgumentParser()
  parser.add_argument("vcf_file")
  return parser.parse_args()

if __name__ == "__main__":
  args = parse()
  vcf = open(args.vcf_file)
  print rewriteVcf(vcf)
