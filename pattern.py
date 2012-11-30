#!/usr/bin/env python
import re

class Call:
  def __init__(self, read_name, reference, positions, call):
    self.read_name = read_name
    self.reference = reference
    self.positions = positions
    num_snp = len(self.positions)
    last = len(call)
    self.call = "".join([ s if s != "" else "N" for s in call]) + (num_snp - last) * "N"

def readData(file):
  file.seek(0)
  header =  file.next().strip().split("\t")
  positions = map(int,header[2:])
  data = []
  for line in file:
    if not reduce(lambda a,b:a or b,
                  map(lambda c : c in line, [",","*","n"])): 
      splited = line.strip().split("\t")
      read_name = splited[0]
      reference = splited[1]
      call = splited[2:]
      data.append(Call(read_name, reference, positions, call))
  return data

def combine(snps):
  cols = zip(*snps)
  counts = map(countAllele, cols)
  res = ""
  for count in counts:
    del count["N"]
    m = max(count.items(), key=lambda i: i[1])
    if m[1] == 0:
       res = res + "N"
    else:
       res = res + m[0]
  return res
  
def countAllele(col):
  count = {"A":0,"C":0,"G":0,"T":0,"N":0}
  for a in col:
    count[a] = count[a] + 1
  return count

def isHetero(col,threshold=0.3):
  c = countAllele(col)
  del c["N"]
  v = sorted(c.values(),reverse=True)
  rate = v[1]/float(v[0])
  return rate > threshold

def countAlleleEachPosition(snps):
  return map(countAllele, zip(*snps))

def countNotN(snp):
  return reduce(lambda a,b: a if b == "N" else a + 1, snp, 0)

class Similarity:
  def __init__(self, a_length, b_length, overlap, score):
    self.a_length = a_length
    self.b_length = b_length
    self.overlap = overlap
    self.score = score

def similarity(snp1, snp2):
  overlap = filter(lambda p: not (p[0] == "N" or p[1] == "N"),
                   zip(snp1, snp2))
  score = reduce(lambda a,b: a if a < 0 else a+b ,
                 map(lambda p: 1 if p[0] == p[1] else -1, overlap),0)
  return Similarity(len(snp1), len(snp2), len(overlap), score)

def makePattern(size):
  base = ["A","C","G","T"]
  s = 1
  res = base
  while s < size:
    res = reduce(lambda a,b: a+b,
                 map(lambda r: map(lambda b: r+b,base),
                     res))
    s = s + 1
  return res

def filterHeterotype(calls_):
  readnames = [ c.read_name for c in calls_ ]
  calls = [ c.call for c in calls_ ]
  cols = zip(*calls)
  only_hapro = filter(isHetero, cols)
  hapro_calls = map(lambda t: reduce(lambda a,b:a+b,t), zip(*only_hapro))
  return zip(readnames, hapro_calls)
  
def alignSnps(calls_):
  calls = filter(lambda call: countNotN(call[1]) > 1,
                 sorted(calls_, key=lambda call: countNotN(call[1])))
  res = []
  while calls:
    sample = calls.pop()
    similar_readnames = [sample[0]]
    similar_calls = [sample[1]]
    NoSimilarCall = False
    while not NoSimilarCall and calls:
      s = filter(lambda c: similarity(c[1],sample[1]).score > 2, calls)
      if s:
        similar_calls.extend([ c[1] for c in s])
        similar_readnames.extend([ c[0] for c in s])
        for call in s:
          calls.remove(call)
      else:
        NoSimilarCall = True
      sample = combine(similar_calls)
    haprotype = sample
    res.append((haprotype, similar_readnames, similar_calls))
  return res

def alignHaprotype(hapros_):
  haprotypes = sorted(filter(lambda h:len(h[1]) > 100, hapros_),
                      key=lambda h:len(h[1]))
  res = []
  while haprotypes:
    init = haprotypes.pop()
    sample = init[0]
    grouped_readnames = init[1]
    grouped_calls = init[2]
    NoSimilarHaprotype = False
    while not NoSimilarHaprotype and haprotypes:
      hapros = filter(lambda h: similarity(sample, h[0]).score > 2, haprotypes)
      if hapros:
        for h in hapros:
          grouped_readnames.extend(h[1])
          grouped_calls.extend(h[2])
          haprotypes.remove(h)
      else:
        NoSimilarHaprotype = True
    new_haprotype = combine(grouped_calls)
    res.append((new_haprotype,grouped_readnames))
  return res

def filter_sam(readnames, in_sam, out_sam):
  for line in in_sam:
    if (line.split()[0] in readnames) or re.match("^@", line):
      print >> out_sam, line.strip()

def parse():
  import argparse
  
  parser = argparse.ArgumentParser()
  parser.add_argument("-n", "--num_haprotype", 
                      default=5,
                      type=int)
  parser.add_argument("-o", "--out_dir",
                      default="./",
                      type=str)
  parser.add_argument("snp_file")
  parser.add_argument("sam_file")
  return parser.parse_args()

debug = False
if debug:
  import sys
  args = parse()
  calls = readData(open(args.snp_file))
  count = countAlleleEachPosition([ c[1] for c in filterHeterotype(calls)])
  for c in count:
    print c

  __name__ = False

if __name__ == "__main__":
  args = parse()
  
  raw_calls = readData(open(args.snp_file))
  allreads = set([ c.read_name for c in raw_calls])
  haprotypes = alignHaprotype(alignSnps(filterHeterotype(raw_calls)))
  
  hapro_fn = re.sub(".txt$", ".hap", args.snp_file.split("/")[-1])
  hapro_file = open(args.out_dir + "/" + hapro_fn,"w")
  for hap in haprotypes:
    print >> hapro_file, "%5i %s" % (len(hap[1])," ".join(hap[0]))
  
  index=1
  for hap in haprotypes[:args.num_haprotype]:
    if index != 1:
      allreads.difference_update(set(hap[1]))
    sam_fn = re.sub(".sam$", 
                    ".hap%02i.sam"%index,
                    args.sam_file.split("/")[-1])
    hap_sam_fn = args.out_dir + "/" + sam_fn
    filter_sam(set(hap[1]),
               open(args.sam_file),
               open(hap_sam_fn,"w"))
    index = index + 1

  filter_sam(allreads,open(args.sam_file), open("test.sam","w"))
