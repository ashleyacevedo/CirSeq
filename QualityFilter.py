import sys,gzip
import numpy as np

workdir = sys.argv[1]
reffile = sys.argv[2]
QualityThreshold=20

infile = gzip.open(workdir + "data.sam.gz","rb")
parameterfile = open(workdir + "ProcessingStats.txt","r")
outfile1 = open(workdir + "Q%sthreshold.txt" % str(QualityThreshold),"w")
outfile2 = open(workdir + "QualityMetrics.txt","w")

reference = open(reffile,"r")
reference.readline()
ref = reference.readlines()
reference = "".join(ref)
reference = reference.replace("\n","")
reference = " " + reference

A = np.zeros((50,len(reference)+200),dtype=int)
C = np.zeros((50,len(reference)+200),dtype=int)
G = np.zeros((50,len(reference)+200),dtype=int)
T = np.zeros((50,len(reference)+200),dtype=int)

AlignedReadCount = 0
for line in infile:
	line = line.split()
	if line[5].count("I") > 0 or line[5].count("D") > 0 or line[5].count("S") > 0:
		pass
	else:
		start = int(line[3])
		i = 0
		while i < len(line[9]):
			Q = ord(line[10][i]) -33
			if line[9][i] == "A":
				A[Q][start + i] += 1
			elif line[9][i] == "C":
				C[Q][start + i] += 1
			elif line[9][i] == "G":
				G[Q][start + i] += 1
			elif line[9][i] == "T":
				T[Q][start + i] += 1	
			i += 1
		AlignedReadCount += 1

print "--------------------------------------"
print "Q", "\t", "Mismatches", "\t", "TotalBases", "\t", "Ts", "\t", "Tv"
outfile2.write("Q" + "\t" + "Mismatches" + "\t" + "TotalBases" + "\t" + "Ts" + "\t" + "Tv" + "\n")

j = 0
Athreshold = [0]*len(reference)
Cthreshold = [0]*len(reference)
Gthreshold = [0]*len(reference)
Tthreshold = [0]*len(reference)
for Q_A, Q_C, Q_G, Q_T in zip(A, C, G, T):
	match = 0
	mismatch = 0
	ts = 0
	tv = 0
	i = 0
	while i < len(reference):
		if reference[i] == "A":
			match += Q_A[i]
			mismatch += Q_C[i] + Q_G[i] + Q_T[i]
			tv += Q_T[i] + Q_C[i]
			ts += Q_G[i]
		elif reference[i] == "C":
			match += Q_C[i]
			mismatch += Q_A[i] + Q_G[i] + Q_T[i]
			tv += Q_A[i] + Q_G[i]
			ts += Q_T[i]
		elif reference[i] == "G":
			match += Q_G[i]
			mismatch += Q_A[i] + Q_C[i] + Q_T[i]
			tv += Q_T[i] + Q_C[i]
			ts += Q_A[i]
		elif reference[i] == "T":
			match += Q_T[i]
			mismatch += Q_A[i] + Q_G[i] + Q_C[i]
			tv += Q_A[i] + Q_G[i]
			ts += Q_C[i]
		i += 1
	if j >= QualityThreshold:
		i = 0
		while i < len(reference):
			Athreshold[i] += Q_A[i]
			Cthreshold[i] += Q_C[i]
			Gthreshold[i] += Q_G[i]
			Tthreshold[i] += Q_T[i]
			i += 1
			
	print j, "\t", mismatch, "\t", mismatch+match, "\t", ts, "\t", tv
	outfile2.write(str(j) + "\t" + str(mismatch) + "\t" + str(mismatch+match) + "\t" + str(ts) + "\t" + str(tv) + "\n")
	j += 1

position = range(0,len(reference))	
for a, c, g, t, base, p in zip(Athreshold, Cthreshold, Gthreshold, Tthreshold, reference, position):
	if p == 0:
		pass
	else:
		outfile1.write(str(p) + "\t" + base + "\t" + str(a) + "\t" + str(c) + "\t" + str(g) + "\t" + str(t) + "\n")

i = 1
for line in parameterfile:
	if i == 6:
		line = line.split()
		TotalConsensusSequences = float(line[0])
	i += 1

parameterfile.close()
parameterfile = open(workdir + "/ProcessingStats.txt", "a")
parameterfile.write("\n" + str(AlignedReadCount) + "\tAligned consensus sequences\n")
parameterfile.write(str(round((AlignedReadCount/TotalConsensusSequences)*100, 2)) + "\t% Consensus sequences aligned\n")
parameterfile.close()

infile.close()
outfile1.close()
outfile2.close()
