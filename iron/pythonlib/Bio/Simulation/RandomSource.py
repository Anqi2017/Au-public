import random

nts = ['A','C','G','T']

#You can asign it a seed or you can asign it another 
# rand = RandomSource
class RandomSource:
  def __init__(self,seed=None):
    self._random = random.Random()
    if seed: self._random.seed(seed)
  
  def random(self):
    return self._random.random()

  def gauss(self,mu,sigma):
    return self._random.gauss(mu,sigma)

  def randint(self,a,b):
    return self._random.randint(a,b)

  def different_random_nt(self,nt):
    global nts
    return self._random.choice([x for x in nts if x != nt.upper()])

  def random_nt(self):
    global nts
    return self._random.choice(nts)
  
  # weights is an array with floats ranging from (0,1]
  # if a random number between 0 and 1 is less than an index return the lowest index
  def get_weighted_random_index(self,weights):
    rnum = self._random.random()
    for i in range(len(weights)):
      if rnum < weights[i]: return i
