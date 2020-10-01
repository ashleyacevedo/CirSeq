"""
ConsensusModule converts tandem repeat reads to consensus sequences.

Author: Ashley Acevedo
"""

# Import required packages
import math
from scipy.stats import mode
from itertools import izip

def GetRepLength(Sequence):

  """
  GetRepLength identifies tandem repeats and returns their length.

  Repeats are identified by locating and saving the distance between substrings.
  The modal distance between substrings (excluding the distance between the last
  substring and the end of the sequence, which is unlikely to be the exact
  repeat length) is the expected repeat length.

  When no repeats are identified, the function returns 0.

  When repeats are either too short (often occurs with repetitive or low
  complexity sequence and with aberent reads) or too long the function returns 1.

  When appropriately sized repeats are identified, the function returns the
  repeat length.
  """

  cdef int i
  cdef int repLength

  SubstringLengths = []
  i = 1
  while i <= len(Sequence) - 8:
    Substrings = Sequence.split(Sequence[i:i+8])
    for String in Substrings[1:-1]:
      SubstringLengths.append(len(String))
    i += 1

  if len(SubstringLengths) == 0:
    repLength = 0
  else:
    repLength = int(mode(SubstringLengths)[0][0]) + 8
    if repLength  < 25 or repLength > math.floor(len(Sequence)/3):
      repLength = 1

  return repLength

def OptimizeReps(Sequence, Qscores, repLength):

  """
  OptimizeReps scores the identity of all possible repeat alignments and
  returns the optimal start site to slice the repeats from the input sequence.

  This procedure is applied because some reads are recombinant. That is, the
  beginning and ending of an input sequence do not originate from the same RNA
  molecule due to RT or PCR cross-overs. Thus, initiating slicing of the repeats
  from the first position in the input sequence may incorporate non-homologus
  sequence into the alignment and reduce the quality of the consensus sequence.

  Input sequence for which no repeat alignment results in an identity score
  above 0.85 are omitted from futher analysis and the function returns 0.

  For input sequences that pass the identity threshold, the function returns 1,
  the sliced repeats and sliced quality scores for those repeats.
  """

  cdef int i
  cdef int identity

  identityScores = []
  i = 0
  while i <= len(Sequence) - repLength * 3:
    identity = 0
    rep1 = Sequence[i:i+repLength]
    rep2 = Sequence[i+repLength:i+repLength*2]
    rep3 = Sequence[i+repLength*2:i+repLength*3]
    for repBase1, repBase2, repBase3 in zip(rep1, rep2, rep3):
      if repBase1 == repBase2 == repBase3 and repBase1 != "N":
        identity += 2
      elif repBase1 == repBase2 and repBase1 != "N":
        identity += 1
      elif repBase1 == repBase3 and repBase1 != "N":
        identity += 1
      elif repBase2 == repBase3 and repBase2 != "N":
        identity += 1
    identityScores.append(identity)
    i += 1

  if max(identityScores)/(repLength*2.0) < 0.85:
    return 0, "", "", "", "", "", ""
  else:
    start = identityScores.index(max(identityScores))
    rep1 = Sequence[start:start+repLength]
    rep2 = Sequence[start+repLength:start+repLength*2]
    rep3 = Sequence[start+repLength*2:start+repLength*3]
    q1 = Qscores[start:start+repLength]
    q2 = Qscores[start+repLength:start+repLength*2]
    q3 = Qscores[start+repLength*2:start+repLength*3]
    return 1, rep1, rep2, rep3, q1, q2, q3

def ConsensusProb(score):
  """Returns error probability for a base matching the consensus."""
  return 10**((ord(score)-33)/-10.0)

def NonConsensusProb(score):
  """Returns error probability for a base not matching the consensus."""
  return 1-(1/3.0)*(10**((ord(score)-33)/-10.0))

def CombineReps(rep1, rep2, rep3, q1, q2, q3):

  """
  CombineReps transforms aligned repeats to a consensus sequence and computes
  combined quality scores.

  The consensus sequence is the majority base in the repeat alignment. When no
  majority exists (when each repeat has a different base at the same position),
  the base with the highest quality score (lowest error probability) is chosen
  as the consensus base.

  Base-calling error probabilities for consensus bases are the product of the
  error probabilities for each repeat base. Error probabilities for each repeat
  base are derived from quality scores (stored as ascii characters) on the
  PHRED+33 scale. The error probabilities of Ns are set to 1 since they are, by
  nature, incorrectly called.
  """

  cdef float prob1
  cdef float prob2
  cdef float prob3
  cdef float probCombined

  consensus = ""
  consensusQ = ""
  for repBase1, repBase2, repBase3, score1, score2, score3 in zip(rep1, rep2, rep3, q1, q2, q3):
    if repBase1 == repBase2 == repBase3 and repBase1 != "N":
      consensus += repBase1
      consensusQ += chr(int(round(((ord(score1) + ord(score2) + ord(score3) - 99)/3.0))) + 33)
    else:
      if repBase1 == repBase2 and repBase1 != "N":
        consensus += repBase1
        prob1 = ConsensusProb(score1)
        prob2 = ConsensusProb(score2)
        prob3 = NonConsensusProb(score3)
      elif repBase1 == repBase3 and repBase1 != "N":
        consensus += repBase1
        prob1 = ConsensusProb(score1)
        prob2 = NonConsensusProb(score2)
        prob3 = ConsensusProb(score3)
      elif repBase2 == repBase3 and repBase2 != "N":
        consensus += repBase2
        prob1 = NonConsensusProb(score1)
        prob2 = ConsensusProb(score2)
        prob3 = ConsensusProb(score3)
      else:
        scores = [ord(score1), ord(score2), ord(score3)]
        repBases = [repBase1, repBase2, repBase3]
        consensus += repBases[scores.index(max(scores))]
        #Sort lists of Bases and Quality scores based on Quality scores
        sortedData = sorted(izip(repBases, scores), key=lambda x: x[1])
        repBases, scores = [[x[i] for x in sortedData] for i in range(2)]
        if repBases[0] == "N":
          prob1 = 1
        else:
          prob1 = NonConsensusProb(chr(scores[0]))
        if repBases[1] == "N":
          prob2 = 1
        else:
          prob2 = NonConsensusProb(chr(scores[1]))
        if repBases[2] == "N":
          prob3 = 1
        else:
          prob3 = ConsensusProb(chr(scores[2]))
      probCombined = prob1 * prob2 * prob3
      consensusQ += chr(int(round(-10*math.log10(probCombined)/3.0))+33)

  return consensus, consensusQ

def Consensus(infile, outfile):

  """
  Consensus reads in FASTQ formatted tandem-repeat sequence reads and writes
  processed consensus sequences and computed quality scores to a new FASTQ
  formatted file.

  This function processes sequences using other functions in the ConsensusModule.
  In addition to writing processed sequences, this function also tabulates the
  number of sequences with repeats of different lengths and keeps track of the
  number of reads that fail to form consensus sequences for a variety of reasons
  (noted in other functions in this module).
  """

  cdef int counter_PoorQuality
  cdef int counter_NoRepeats
  cdef int counter_AbnormalRepeatLength
  cdef int counter_LowIdentity
  cdef int counter_TotalReads
  cdef int counter_ConsensusSequences

  cdef int repLength
  cdef int flag

  counter_PoorQuality = 0
  counter_NoRepeats = 0
  counter_AbnormalRepeatLength = 0
  counter_LowIdentity = 0
  counter_TotalReads = 0
  counter_ConsensusSequences = 0

  RepeatLengths = [0]*115
	
  while True:

    # Read FASTQ formatted sequence entry
    SequenceID = infile.readline()
    Sequence = infile.readline()
    EmptyLine = infile.readline()
    QualityScores = infile.readline()

    if not SequenceID:
      break

    counter_TotalReads += 1

    #Remove the first base which is error prone
    Sequence = Sequence.rstrip("\n")
    Sequence = Sequence[1:]
    QualityScores = QualityScores.rstrip("\n")
    QualityScores = QualityScores[1:]
	
    #Remove/count reads with more than 30% Ns
    if Sequence.count("N") > len(Sequence)/3:
      counter_PoorQuality += 1
			
    else:
      repLength = GetRepLength(Sequence)

      if repLength == 0:
        counter_NoRepeats += 1
      elif repLength == 1:
        counter_AbnormalRepeatLength += 1
      else:
        RepeatLengths[repLength] += 1
        flag, rep1, rep2, rep3, q1, q2, q3 = OptimizeReps(Sequence, QualityScores, repLength)

        if flag == 0:
          counter_LowIdentity += 1
        else:
          counter_ConsensusSequences += 1
          consensus, consensusQ = CombineReps(rep1, rep2, rep3, q1, q2, q3)
								
          outfile.write(SequenceID)
          outfile.write(consensus + "\n")
          outfile.write(EmptyLine)
          outfile.write(consensusQ + "\n")

  return counter_PoorQuality, counter_NoRepeats, counter_AbnormalRepeatLength, \
         counter_LowIdentity, counter_ConsensusSequences, counter_TotalReads, \
         RepeatLengths
