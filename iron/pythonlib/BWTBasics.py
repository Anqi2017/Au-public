def bwt(inseq):
  seq = list(inseq)
  seq.append("|")
  strings = []
  for i in range(0,len(seq)):
    rot = seq[i:len(seq)]
    [rot.append(x) for x in seq[0:i]]
    strings.append(''.join(rot))
  return ''.join([x[-1] for x in sorted(strings)])
  
