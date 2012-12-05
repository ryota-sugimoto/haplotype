#!/usr/bin/env python
import sys

def read_file(f):
  header = f.next().strip().split("\t")
  positions = map(int, header[2:])
  num_snp = len(positions)
  res = []
  for line in f:
    if not reduce(lambda a,b: a or b,
                  map(lambda c: c in line, [",","*","n"])):
      splited = line.strip().split("\t")
      read_name = splited[0]
      reference = splited[1]
      snps = "".join([ s if s != "" else "N" for s in splited[2:]]) + (num_snp-len(splited[2:])) * "N"
      d = {"read_name": read_name,
           "reference": reference}
      d["snp"] = dict(zip(positions, snps))
      res.append(d)
  return res

def snps_at_position(data, position):
  return map(lambda d: d["snp"][position], data)

def count_allele(alleles):
  count = {"A":0, "C":0, "G":0, "T":0, "N":0}
  for a in alleles:
    count[a] = count [a] + 1
  return count

def major_allele(count_, threshold=0.2):
  count = count_.copy()
  del count["N"]
  print count
  total = float(sum(count.values()))
  return filter(lambda a: count[a]/total > threshold,
                ["A","C","G","T"])

def count_pair(data, pos1, pos2):
  count = {}
  for read in data:
    pair = read["snp"][pos1] + read["snp"][pos2]
    if "N" not in pair:
      count[pair] = count.setdefault(pair,0) + 1
  return count

def major_pair(count, threshold=0.1):
  total = float(sum(count.values()))
  res = []
  for pair in count.keys():
    if count[pair]/total > threshold:
      res.append(pair)
  return res
  

def heterotype_position(data):
  positions = sorted(data[0]["snp"].keys())
  return filter(lambda p: len(major_allele(
                                count_allele(
                                  snps_at_position(data,p)))) == 2, positions)


def chain_haprotype(data, max_chain_distance=1300):
  het_pos = heterotype_position(data)
  
  heterotype_pairs = [ major_allele(
                         count_allele(
                           snps_at_position(data, p))) for p in het_pos ]

  chaining_positions = []
  for i1 in range(len(het_pos)-1):
    p1 = het_pos[i1]
#    for i2 in range(i1, len(het_pos))
    p2 = het_pos[i1+1]
    count = count_pair(data, p1, p2)
    chains = major_pair(count)
    if len(chains) != 2:
      msg1 = "WARN:Too many or little chains between %i and %i.\n" % (p1,p2)
      msg2 = ", ".join(["%s: %i" % (key, count[key]) for key in chains])
      sys.stderr.write(msg1 + "\t" +  msg2 + "\n")
 
    
  
    

import argparse
import time
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("snp_file")
  args = parser.parse_args()
  data = read_file(open(args.snp_file))
  positions = sorted(data[0]["snp"].keys())
  #for p in positions:
  #  count = count_allele(snps_at_position(data, p))
  #  print p,count, major_allele(count)
  #count = count_pair(data, 4492, 4568)
  #print count
  #print major_pair(count)
  #print heterotype_position(data)
  chain_haprotype(data)
