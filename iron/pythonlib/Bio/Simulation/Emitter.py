from Bio.Simulation.RandomSource import RandomSource

#Give it a transcriptome definition and a reference genome for it
#initialy give it uniform probability
class TranscriptomeEmitter:
  def __init__(self,transcriptome,seed=None,rand=None):
    if rand: self.random = rand
    elif seed: self.random = RandomSource(seed)
    else: self.random = RandomSource()

    self._transcriptome = transcriptome
    ######
    tcnt = len(self._transcriptome.get_transcripts())
    self._weights = [float(i+1)/float(tcnt) for i in range(0,tcnt)]
    ## _log stores what we are emitting ##
    self._log = []

  def emit_transcript(self):
    i = self.random.get_weighted_random_index(self._weights)
    return self._transcriptome.get_transcripts()[i]
