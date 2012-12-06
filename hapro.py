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
  total = float(sum(count.values()))
  try:
    return filter(lambda a: count[a]/total > threshold,
                  ["A","C","G","T"])
  except ZeroDivisionError:
    return ["N"]

def count_pair(data, pos1, pos2):
  count = {}
  for read in data:
    pair = read["snp"][pos1] + read["snp"][pos2]
    if "N" not in pair:
      count[pair] = count.setdefault(pair,0) + 1
  return count

def major_pair(count, threshold=0.2):
  total = float(sum(count.values()))
  res = []
  for pair in count.keys():
    if count[pair]/total > threshold:
      res.append(pair)
  return res
  
def rich_read_position(data, min_read_num = 1):
  positions = sorted(data[0]["snp"].keys())
  res = []
  for p in positions:
    count = count_allele(snps_at_position(data, p))
    del count["N"]
    if sum(count.values()) > min_read_num:
      res.append(p)
  return res

def heterotype_position(data):
  positions = rich_read_position(data)
  return filter(lambda p: len(major_allele(
                                count_allele(
                                  snps_at_position(data,p)))) == 2, positions)

def connect_pair(a, b, pair1, pair2):
  a_match_pair = filter(lambda p: a == p[0], [pair1, pair2])
  b_match_pair = filter(lambda p: b == p[0], [pair1, pair2])
  if len(a_match_pair) == len(b_match_pair) == 1:
    return [a_match_pair[0][1], b_match_pair[0][1]]
  else:
    return []
    
def chain_haprotype(data):
  het_pos = heterotype_position(data)
  
  res = []
  prev_pos = het_pos[0]
  prev_pair = major_allele(count_allele(snps_at_position(data, prev_pos)))
  current_haprotype = {prev_pos: prev_pair}
  for current_pos in het_pos[1:]:
    count = count_pair(data, prev_pos, current_pos)
    chains = major_pair(count)
    if len(chains) == 2:
      if (chains[0][0] != chains[1][0]) and (chains[0][1] != chains[1][1]):
        print "pos",current_pos
        print "prev_pair",prev_pair
        print "chains",chains
        current_pair = connect_pair(prev_pair[0],
                                    prev_pair[1],
                                    chains[0],
                                    chains[1])
        print "current_pair",current_pair
        prev_pair = current_pair
        current_haprotype[current_pos] = current_pair
      else:
        msg1 = "WARN: Chains have been converged at %i\n" % (current_pos)
        msg2 = "\t%s %s-%s\n\t%s %s-%s\n" % (prev_pair[0],
                                             chains[0][0],chains[0][1],
                                             prev_pair[1],
                                             chains[1][0],chains[1][1])
        sys.stderr.write(msg1 + msg2)
        prev_pair = major_allele(
                      count_allele(
                        snps_at_position(data, current_pos)))
        res.append(current_haprotype)
        current_haprotype = {prev_pos: prev_pair}
    else:
      msg1 = "WARN: Too many or little chains between %i and %i.\n"
      msg2 = ", ".join(["%s: %i" % (key, count[key]) for key in chains])
      sys.stderr.write(msg1 % (prev_pos, current_pos) + "\t" +  msg2 + "\n")
      prev_pair = major_allele(
                    count_allele(
                      snps_at_position(data, current_pos))) 
      res.append(current_haprotype)
      current_haprotype = {prev_pos: prev_pair}
    prev_pos = current_pos
  res.append(current_haprotype)
  return res

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
  hapros = chain_haprotype(data) 
  for hap in hapros:
    for key in sorted(hap.keys()):
      print key, hap[key]
    print
