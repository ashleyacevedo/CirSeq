import sys,gzip

workdir = sys.argv[1]

infiles = ["/9_alignment.sam.gz","/10_alignment.sam.gz"]
outfile = gzip.open(workdir + "/11_alignment.sam.gz","wb")

for file in infiles:
	infile = gzip.open(workdir + file,"rb")
	
	#running lists
	Alignment = []
	AlignmentScore = []
	
	#current read
	SequenceID = ""
	
	for line in infile:
		RawLine = line
		line = line.split()
		
		#add alignments of the same read to the running list
		if line[0] == SequenceID:
			#add only alignments lacking gaps and clipped bases
			if line[5].count("D") == 0 and line[5].count("I") == 0 and line[5].count("S") == 0:
				AS = line[11].split(":") 
				AlignmentScore.append(int(AS[2]))
				Alignment.append(RawLine)
				
		else:
			if len(Alignment) > 0:
				#write alignment with the best alignment score
				outfile.write(Alignment[AlignmentScore.index(max(AlignmentScore))])
			
			SequenceID = line[0]
			
			Alignment = []
			AlignmentScore = []
			
			if line[5].count("D") == 0 and line[5].count("I") == 0 and line[5].count("S") == 0:
				AS = line[11].split(":") 
				AlignmentScore.append(int(AS[2]))
				Alignment.append(RawLine)



	if len(Alignment) > 0:
		#write alignment with the best alignment score
		outfile.write(Alignment[AlignmentScore.index(max(AlignmentScore))])			