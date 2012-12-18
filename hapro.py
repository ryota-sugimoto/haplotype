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
      snps = "".join([ s if s != "" else "N" for s in splited[2:]]) \
              + (num_snp-len(splited[2:])) * "N"
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

def major_allele(count_, threshold=0.3):
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
  
def rich_read_position(data, min_read_num = 100):
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
        current_pair = connect_pair(prev_pair[0],
                                    prev_pair[1],
                                    chains[0],
                                    chains[1])
        current_haprotype[current_pos] = current_pair
      else:
        msg1 = "WARN: Chains have been converged at %i\n" % (current_pos)
        msg2 = "\t%s %s-%s\n\t%s %s-%s\n" % (prev_pair[0],
                                             chains[0][0],chains[0][1],
                                             prev_pair[1],
                                             chains[1][0],chains[1][1])
        sys.stderr.write(msg1 + msg2)
        current_pair = major_allele(
                         count_allele(
                           snps_at_position(data, current_pos)))
        res.append(current_haprotype)
        current_haprotype = {current_pos: current_pair}
    else:
      msg1 = "WARN: Too many or little chains between %i and %i.\n"
      msg2 = ", ".join(["%s: %i" % (key, count[key]) for key in chains])
      sys.stderr.write(msg1 % (prev_pos, current_pos) + "\t" +  msg2 + "\n")
      current_pair = major_allele(
                       count_allele(
                         snps_at_position(data, current_pos))) 
      res.append(current_haprotype)
      current_haprotype = {current_pos: current_pair}
    prev_pair = current_pair
    prev_pos = current_pos
  res.append(current_haprotype)
  return res

def match_score(str1, str2):
  return reduce(lambda a,b: {"match": a["match"] + 1,
                             "mismatch":a["mismatch"]} if b > 0 \
                            else {"match": a["match"],
                                  "mismatch": a["mismatch"] + 1},
                map(lambda p: 1 if p[0] == p[1] else -1,
                    filter(lambda p: p[0] != "N" and p[1] != "N",
                           zip(str1, str2))),
                {"match": 0, "mismatch":0})

def align_reads(data, hapro, min_contig_length=10):
  hap = filter(lambda contig: len(contig.keys()) > min_contig_length, hapro) 
  res = []
  for contig in hap:
    positions = sorted(contig.keys())
    hapro_snp = map(lambda i: "".join([ contig[p][i] for p in positions ]),
                    [0,1])
    appending = {0: [], 1: []}
    for read in data:
      read_name = read["read_name"]
      read_snp = "".join([ read["snp"][p] for p in positions ])
      scores = map(lambda h: match_score(h, read_snp), hapro_snp)
      for hap_num in [0,1]:
        if (scores[hap_num]["match"] > 2) \
           and (scores[hap_num]["mismatch"] == 0):
          appending[hap_num].append(read_name)
    res.append(appending)
  return res

def filter_sam(readnames, in_sam, out_sam):
  in_sam.seek(0)
  set_readnames = set(readnames)
  for line in in_sam:
    if (line.split()[0] in set_readnames) or re.match("^@", line):
      print >> out_sam, line.strip()

def parse():
  import argparse
  
  parser = argparse.ArgumentParser()
  parser.add_argument("-m", "--min_haprotype_contig_size",
                      help="minimum size of haprotype contig",
                      default=10,
                      type=int)
  parser.add_argument("-o", "--out_dir",
                      default="./",
                      type=str)
  parser.add_argument("snp_file")
  parser.add_argument("sam_file")
  return parser.parse_args()

import os
import sys
def check_file(name, check_func=os.path.isfile):
  if not check_func(name):
    print >> sys.stderr, name, "doesn't exist."
    exit(1)

import re
def main():
  args = parse()
  check_file(args.out_dir, os.path.isdir)
  check_file(args.snp_file)
  check_file(args.sam_file)
  in_sam = open(args.sam_file)
  data = read_file(open(args.snp_file))
  hapro = chain_haprotype(data)
  for hap in hapro:
    for p in sorted(hap.keys()):
      print p, hap[p]
    print
  align = align_reads(data, hapro)
  contig_index = 0
  for contig in align:
    for hap in contig.keys():
      suffix = "_c%02i_h%02i.sam" % (contig_index, hap)
      out_sam_fn = re.sub("\.sam$",
                          suffix,
                          args.sam_file.split("/")[-1])
      out_sam_path = args.out_dir + "/" + out_sam_fn
      out_sam = open(out_sam_path,"w")
      filter_sam(contig[hap], in_sam, out_sam)
    contig_index = contig_index + 1
  
if __name__ == "__main__":
  main()
