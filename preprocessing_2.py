import sys,gzip

workdir = sys.argv[1]

infile = gzip.open(workdir+"/6_alignment.sam.gz","rb")
outfile1 = gzip.open(workdir+"/7_alignment.sam.gz","wb")
outfile2 = gzip.open(workdir+"/8_rotated.fastq.gz","wb")

def Rotate(line):
	i = 0
	while i < len(line[9]):
		outfile2.write("@" + line[0] + "\n")
		outfile2.write(line[9][-i:] + line[9][:-i] + "\n")
		outfile2.write("+" + "\n")
		outfile2.write(line[10][-i:] + line[10][:-i] + "\n")
		i += 1

for line in infile:
	line1 = line
	line = line.split()
	
	#make every possible rotation of reads with gaps and more than one clipped sequence
	if line[5].count("D") > 0 or line[5].count("I") > 0 or line[5].count("S") > 1:
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
	
	#make every possible rotation of reads that still have clipped bases following rearrangement
	elif line[5].count("S") == 1:
		Rotate(line)

infile.close()
outfile1.close()
outfile2.close()
