import random, math
def poisson_random_number(lamb):
  x = 0
  p = math.exp(-1*float(lamb))
  s = p
  u = random.uniform(0,1)
  while u > s:
    x += 1
    p = p*lamb/x
    s = s + p
  return x

def probability_threshold(lamb,prob):
  vals = []
  for i in range(0,10000):
    vals.append(poisson_random_number(lamb))
  vals.sort()
  ind = int(prob*len(vals))
  if ind >= len(vals): ind = len(vals)-1
  return vals[int(prob*len(vals))]
