import sys,gzip

workdir = sys.argv[1]

infile   = gzip.open(workdir + "/2_alignment.sam.gz","rb")
outfile1 = gzip.open(workdir + "/3_alignment.sam.gz","wb") #perfectly matched, ungapped alignemnts
outfile2 = gzip.open(workdir + "/4_rearranged.fastq.gz","wb") #clipped reads with order of sequence blocks swapped
outfile3 = gzip.open(workdir + "/5_rotated.fastq.gz","wb") # every possible rotation of gapped and multi-clipped reads

def Rotate(line):
	i = 0
	while i < len(line[9]):
		outfile3.write("@" + line[0] + "\n")
		outfile3.write(line[9][-i:] + line[9][:-i] + "\n")
		outfile3.write("+" + "\n")
		outfile3.write(line[10][-i:] + line[10][:-i] + "\n")
		i += 1

for line in infile:
	line1 = line
	line = line.split()
	
	#make every possible rotation of reads with gaps and more than one clipped sequence
	if line[5].count("D") > 0 or line[5].count("I") > 0 or line[5].count("S") > 1 or line[5].count("*"): 
		Rotate(line)
	
	elif line[5].count("S") == 0:
		xm = line[13].split(":") #xm = number of mismatches
		
		#save perfectly matched, ungapped alignments
		if int(xm[2]) == 0:
			outfile1.write(line1)
		
		#make every possible rotation of reads with no gaps but have mismatches. Coincidentally matched bases near the
		#ends may suppress sequence clipping
		else:
			Rotate(line)
	
	#rearrange sequences with a single block of clipped reads
	elif line[5].count("S") == 1:
		
		#parse CIGAR
		numbers = ""
		letters = []
		for character in line[5]:
			if character == "M" or character == "S":
				letters.append(character)
				numbers += "X"
			else:
				numbers += character
		numbers = numbers.split("X")
		
		#rearrange seuences and quality scores
		Sequence = line[9][int(numbers[0]):] + line[9][0:int(numbers[0])]
		QualityScores = line[10][int(numbers[0]):] + line[10][0:int(numbers[0])]
		outfile2.write("@" + line[0] + "\n")
		outfile2.write(Sequence + "\n")
		outfile2.write("+" + "\n")
		outfile2.write(QualityScores + "\n")

infile.close()
outfile1.close()
outfile2.close()
outfile3.close()
