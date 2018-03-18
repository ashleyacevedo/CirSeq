def Consensus(infile, outfile):
	
	import math
	from scipy.stats import mode
	from itertools import izip
	
	cdef int i
	cdef int counter_PoorQuality
	cdef int counter_NoRepeats
	cdef int counter_AbnormalRepeatLength
	cdef int counter_LowIdentity
	cdef int counter_TotalReads
	cdef int counter_ConsensusSequences
	
	cdef int StartPosition
	cdef int Identity
	
	counter_PoorQuality = 0
	counter_NoRepeats = 0
	counter_AbnormalRepeatLength = 0
	counter_LowIdentity = 0
	counter_TotalReads = 0
	counter_ConsensusSequences = 0

	RepeatLengths = [0]*100
	
	while True:
		
		SequenceID = infile.readline()		
		counter_TotalReads += 1
		
		if not SequenceID:					
			break						
		
		Sequence = infile.readline()
		Sequence = Sequence.rstrip("\n")
		EmptyLine = infile.readline()
		QualityScores = infile.readline()
		QualityScores = QualityScores.rstrip("\n")
	
		#Remove the first base which can be error prone
		Sequence = Sequence[1:]
		QualityScores = QualityScores[1:]
	
		#Remove/count reads with more than 30% Ns
		if Sequence.count("N") > 90:
			counter_PoorQuality += 1
			
		else:
		
			#Identify all possible sets of repeats and save their lengths
			SubstringLengths = []
			i = 1
			while i <= 291:
				Substrings = Sequence.split(Sequence[i:i+8])
				for String in Substrings[1:-1]: #[1:-1] avoids ends which will likely not be the exact repeat length
					SubstringLengths.append(len(String))
				i += 1
			
			#Remove/count reads with no detectable repeats
			if len(SubstringLengths) == 0:
				counter_NoRepeats += 1
				
			else:
				SubstringLengthsMode = mode(SubstringLengths) #mode is an array [0][0] = mode, [1][0] = counts of mode in list
				
				#Remove/count reads with abnormal/non-ideal repeat length
				if SubstringLengthsMode[0][0] < 17 or SubstringLengthsMode[0][0] > 91:
					counter_AbnormalRepeatLength += 1
				
				else:
					#Optimize sequence identity, all possible sets of repeats are considered
					RepeatLength = int(SubstringLengthsMode[0][0]) + 8
					RepeatLengths[RepeatLength] += 1
					RepeatSetIdentities = []
					i = 0
					while i < 299 - (RepeatLength*3):
						Identity = 0
						for RepeatBase1, RepeatBase2, RepeatBase3 in zip(Sequence[i:i+RepeatLength], Sequence[i+RepeatLength:i+RepeatLength*2], Sequence[i+RepeatLength*2:i+RepeatLength*3]):
							if RepeatBase1 == RepeatBase2 == RepeatBase3 and RepeatBase1 != "N":
								Identity += 2
							elif RepeatBase1 == RepeatBase2 and RepeatBase1 != "N":
								Identity += 1
							elif RepeatBase1 == RepeatBase3 and RepeatBase1 != "N":
								Identity += 1
							elif RepeatBase2 == RepeatBase3 and RepeatBase2 != "N":
								Identity += 1
						RepeatSetIdentities.append(Identity)
						i += 1
					
					#Select the highest identity set of repeats that is at least 85% identity
					if max(RepeatSetIdentities)/(RepeatLength*2.0) < 0.85:
						counter_LowIdentity += 1	
					else:
						counter_ConsensusSequences += 1
						StartPosition = RepeatSetIdentities.index(max(RepeatSetIdentities))
		
						#Form consensus sequences and recalculate quality scores
						ConsensusSequence = ""
						ReCalculatedQualityScores = ""
						for RepeatBase1, RepeatBase2, RepeatBase3, RepeatQualityScore1, RepeatQualityScore2, RepeatQualityScore3 in zip(Sequence[StartPosition:StartPosition+RepeatLength], Sequence[StartPosition+RepeatLength:StartPosition+RepeatLength*2], Sequence[StartPosition+RepeatLength*2:StartPosition+RepeatLength*3],QualityScores[StartPosition:StartPosition+RepeatLength], QualityScores[StartPosition+RepeatLength:StartPosition+RepeatLength*2], QualityScores[StartPosition+RepeatLength*2:StartPosition+RepeatLength*3]):
							if RepeatBase1 == RepeatBase2 == RepeatBase3 and RepeatBase1 != "N":
								ConsensusSequence += RepeatBase1
								ReCalculatedQualityScores += chr(int(round(((ord(RepeatQualityScore1) + ord(RepeatQualityScore2) + ord(RepeatQualityScore3) - 99)/3.0))) + 33)
							elif RepeatBase1 == RepeatBase2 and RepeatBase1 != "N":
								ConsensusSequence += RepeatBase1
								RawProbability = (10**((ord(RepeatQualityScore1)-33)/-10.0))*(10**((ord(RepeatQualityScore2)-33)/-10.0))*(1-(1/3.0)*(10**((ord(RepeatQualityScore3)-33)/-10.0)))
								ReCalculatedQualityScores += chr(int(round(-10*math.log10(RawProbability)/3.0))+33)
							elif RepeatBase1 == RepeatBase3 and RepeatBase1 != "N":
								ConsensusSequence += RepeatBase1
								RawProbability = (10**((ord(RepeatQualityScore1)-33)/-10.0))*(10**((ord(RepeatQualityScore3)-33)/-10.0))*(1-(1/3.0)*(10**((ord(RepeatQualityScore2)-33)/-10.0)))
								ReCalculatedQualityScores += chr(int(round(-10*math.log10(RawProbability)/3.0))+33)
							elif RepeatBase2 == RepeatBase3 and RepeatBase2 != "N":
								ConsensusSequence += RepeatBase2
								RawProbability = (10**((ord(RepeatQualityScore3)-33)/-10.0))*(10**((ord(RepeatQualityScore2)-33)/-10.0))*(1-(1/3.0)*(10**((ord(RepeatQualityScore1)-33)/-10.0)))
								ReCalculatedQualityScores += chr(int(round(-10*math.log10(RawProbability)/3.0))+33)
							else: #All bases are different or Ns
								RepeatQualityScores = [ord(RepeatQualityScore1), ord(RepeatQualityScore2), ord(RepeatQualityScore3)]
								RepeatBases = [RepeatBase1, RepeatBase2, RepeatBase3]
								ConsensusSequence += RepeatBases[RepeatQualityScores.index(max(RepeatQualityScores))]
								#Sort lists of Bases and Quality scores based on Quality scores
								SortedBasesAndQualityScores = sorted(izip(RepeatBases, RepeatQualityScores), key=lambda x: x[1])
								RepeatBases, RepeatQualityScores = [[x[i] for x in SortedBasesAndQualityScores] for i in range(2)]
								if RepeatBases[0] == "N":
									RawProbability1 = 1
								else:
									RawProbability1	= 1-(1/3.0)*(10**((RepeatQualityScores[0]-33)/-10.0))
								if RepeatBases[1] == "N":
									RawProbability2 = 1
								else:
									RawProbability2	= 1-(1/3.0)*(10**((RepeatQualityScores[1]-33)/-10.0))
								RawProbability = RawProbability1 * RawProbability2 * (10**((RepeatQualityScores[2]-33)/-10.0))
								ReCalculatedQualityScores += chr(int(round(-10*math.log10(RawProbability)/3.0))+33)
								
						outfile.write(SequenceID)
						outfile.write(ConsensusSequence + "\n")
						outfile.write(EmptyLine)
						outfile.write(ReCalculatedQualityScores + "\n")
							
	return counter_PoorQuality, counter_NoRepeats, counter_AbnormalRepeatLength, counter_LowIdentity, counter_ConsensusSequences, counter_TotalReads, RepeatLengths