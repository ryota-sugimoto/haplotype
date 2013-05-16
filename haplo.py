#!/usr/bin/env python

#Copyright (c) 2013 All Right Reserved.
#
#author: Ryota Sugimoto
#institute: National Institute of Genetics
#email: ryota.sugimoto@gmail.com
#date: 2013-05-16

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

def remove_range(data, bed_data):
  for read in data:
    for range in bed_data:
      if read["reference"] == range["reference"]:
        for pos in read["snp"].keys():
          if (pos > range["begin"]) and (pos <= range["end"]):
            del read["snp"][pos]

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

def heterotype_position(data, min_snp_rate, min_depth):
  positions = rich_read_position(data, min_depth)
  return filter(lambda p: len(major_allele(
                                count_allele(
                                  snps_at_position(data,p)),
                                          min_snp_rate)) == 2, positions)

def connect_pair(a, b, pair1, pair2):
  a_match_pair = filter(lambda p: a == p[0], [pair1, pair2])
  b_match_pair = filter(lambda p: b == p[0], [pair1, pair2])
  if len(a_match_pair) == len(b_match_pair) == 1:
    return [a_match_pair[0][1], b_match_pair[0][1]]
  else:
    return []
    
def no_converged_chain(chains_):
  chains = chains_[:]
  res = []
  while chains:
    c1 = chains.pop()
    for c2 in chains:
      if c1[0] != c2[0] and c1[1] != c2[1]:
        chains.remove(c2)
        res.append([c1,c2])
  if len(res) != 1:
    sys.stderr.write("I couldn't find a pair of unique chains\n")
    sys.stderr.write(str(res))
    exit(1)
  return res[0]

def chain_haplotype(data, snp_threshold,
                          haplo_chain_threshold,
                          min_snp_depth,
                          ignore_range,
                          force_connect=True):
  remove_range(data, ignore_range)
  het_pos = heterotype_position(data, snp_threshold, min_snp_depth)
  
  res = []
  try:
    prev_pos = het_pos[0]
  except IndexError:
    sys.stderr.write("ERROR: There is no snp found.\n")
    exit(1)
  prev_pair = major_allele(count_allele(snps_at_position(data, prev_pos)),
                           snp_threshold)
  current_haplotype = {prev_pos: prev_pair}
  for current_pos in het_pos[1:]:
    count = count_pair(data, prev_pos, current_pos)
    chains = major_pair(count,
                        haplo_chain_threshold)
    if len(chains) == 2 and len(prev_pair) == 2:
      if (chains[0][0] != chains[1][0]) and (chains[0][1] != chains[1][1]):
        current_pair = connect_pair(prev_pair[0],
                                    prev_pair[1],
                                    chains[0],
                                    chains[1])
        current_haplotype[current_pos] = current_pair
      else:
        msg1 = "WARN: Chains have been converged at %i\n" % (current_pos)
        msg2 = "\t%s %s-%s\n\t%s %s-%s\n" % (prev_pair[0],
                                             chains[0][0],chains[0][1],
                                             prev_pair[1],
                                             chains[1][0],chains[1][1])
        sys.stderr.write(msg1 + msg2)
        current_pair = major_allele(
                         count_allele(
                           snps_at_position(data, current_pos)),
                                   snp_threshold)
        res.append(current_haplotype)
        current_haplotype = {current_pos: current_pair}
    elif len(chains) >= 3 and len(prev_pair) == 2 and force_connect == True:
      msg1 = "WARN: There is %i possible chains between %i and %i."
      msg2 = ", ".join(["%s: %i" % (key, count[key]) for key in chains])
      msg3 = "I'm going to assume %s and %s are correct chains.\n"
      assumed_chain_1 = sorted(count.keys(),
                               reverse=True,
                               key=lambda key: count[key])[0]
      assumed_chains = no_converged_chain(chains)
      current_pair = connect_pair(prev_pair[0],
                                  prev_pair[1],
                                  assumed_chains[0],
                                  assumed_chains[1])
      current_haplotype[current_pos] = current_pair
      msg = msg1 % (len(chains),prev_pos,current_pos) + "\n" \
          + "\t" + msg2 + "\n" \
          + "\t" + msg3 % tuple(assumed_chains)
      sys.stderr.write(msg)
    else:
      msg1 = "WARN: Too many or small chains between %i and %i.\n"
      msg2 = ", ".join(["%s: %i" % (key, count[key]) for key in chains])
      sys.stderr.write(msg1 % (prev_pos, current_pos) + "\t" +  msg2 + "\n")
      current_pair = major_allele(
                       count_allele(
                         snps_at_position(data, current_pos)),
                                 snp_threshold) 
      res.append(current_haplotype)
      current_haplotype = {current_pos: current_pair}
    prev_pair = current_pair
    prev_pos = current_pos
  res.append(current_haplotype)
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

def align_reads(data, haplo, min_contig_length=10):
  hap = filter(lambda contig: len(contig.keys()) > min_contig_length, haplo) 
  res = []
  for contig in hap:
    positions = sorted(contig.keys())
    haplo_snp = map(lambda i: "".join([ contig[p][i] for p in positions ]),
                    [0,1])
    appending = {0: [], 1: []}
    for read in data:
      read_name = read["read_name"]
      read_snp = "".join([ read["snp"][p] for p in positions ])
      scores = map(lambda h: match_score(h, read_snp), haplo_snp)
      for hap_num in [0,1]:
        if (scores[hap_num]["match"] > 1) \
           and (scores[hap_num]["mismatch"] <= 0):
          appending[hap_num].append(read_name)
    res.append(appending)
  return res

def filter_sam(readnames, in_sam, out_sam):
  in_sam.seek(0)
  set_readnames = set(readnames)
  for line in in_sam:
    readname = line.split()[0]
    if (readname in set_readnames) or re.match("^@", line):
      print >> out_sam, line.strip()

def bed_parser(file):
  res = []
  for line in file:
    (ref, range) = line.strip().split(":")
    (begin, end) = map(int,range.split("-"))
    res.append({"reference": ref, "begin": begin, "end": end})
  return res

def parse():
  import argparse
  
  parser = argparse.ArgumentParser()
  parser.add_argument("-m", "--min_haplotype_contig_size",
                      help="minimum size of haplotype contig",
                      action="store",
                      default=10,
                      type=int)
  parser.add_argument("-o", "--out_dir",
                      help="specify output directory",
                      default="./",
                      type=str)
  parser.add_argument("-s", "--snp_rate_threshold",
                      help="minimum rate of snp number",
                      default=0.3,
                      type=float)
  parser.add_argument("-H", "--haplo_neighbor_rate_threshold",
                      help="minimum rate of chained snp number",
                      default=0.2,
                      type=float)
  parser.add_argument("-M", "--minimum_read_depth_at_snp_pos",
                      default=100,
                      type=int)
  parser.add_argument("-i", "--ignore_range_bed",
                      type=str)
  parser.add_argument("-f", "--force_connect",
                      help="force to connecting the haplo chains if possible.",
                      action="store_true")
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
  if args.ignore_range_bed:
    check_file(args.ignore_range_bed)
    ignore_range = bed_parser(open(args.ignore_range_bed))
  else:
    ignore_range = {}
  in_sam = open(args.sam_file)
  data = read_file(open(args.snp_file))
  haplo = chain_haplotype(data,
                          args.snp_rate_threshold,
                          args.haplo_neighbor_rate_threshold,
                          args.minimum_read_depth_at_snp_pos,
                          ignore_range,
                          args.force_connect)
  for hap in haplo:
    for p in sorted(hap.keys()):
      try:
        print >> sys.stderr, "%i: %s %s" % tuple([p] + hap[p])
      except TypeError:
        print >> sys.stderr, "%i: N N"  % tuple([p])
  align = align_reads(data, haplo, args.min_haplotype_contig_size)
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
