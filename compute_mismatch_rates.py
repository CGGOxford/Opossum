# Laura Oikkonen 7.12.2016

# Code requires version 0.8.1 of Pysam package
import pysam
import numpy as np


# Compute base mismatch rates for first and second strands separately

# Takes a bam file as input - note that file needs to contain md tags,
# they can be generated with samtools calmd.

# Code outputs two txt files, corresponding to first and second strands.
# The output consists of N 4x4 matrices (N=read length) with the number
# of A,C,G,T bases at each read position (more specifically, how many A,
# C,G,T bases when reference base is A, etc). Output can be plotted 
# using plot_mismatch_rates.py.

def main() :

#######################################################################
	# FILL THESE IN
	# Input bam file
	bamfile='input.bam'
	# Output dir (where firstbases.txt and secondbases.txt will be created)
	outputdir='~/output_folder/'
	# Read length (code expects all reads to be of same length)
	readlen=76
######################################################################

	try :
		samfile = pysam.AlignmentFile(bamfile, "rb")
	except :
		print 'Could not open file'



	firstbases = np.zeros((readlen,4,4)) # number of A, C, G, T bases in first strand at each position
	lastbases = np.zeros((readlen,4,4)) # bases in second strand

	discarded = 0
	firstforward = 0
	secondreverse = 0

	# Read in all reads from bam file
	cont = True
	while cont :

		try :
			tryread = samfile.next()

			# Discard reads with very poor quality
			if tryread.mapq < 20 :
				discarded += 1 
				continue

		except :
			cont = False
			break

		mismatches = BaseMismatches(tryread.opt('MD'), tryread.cigar)

		if (tryread.flag & 64) == 64 : # first strand in template

			# Go through each base in read
			readseq = list(tryread.seq)
	
			for m,pos in mismatches :
				alt = readseq[pos]
				if alt == 'N' : continue
				readseq[pos] = 'N'

				if (tryread.flag & 16) == 16 : pos = readlen-1-pos # reverse strand
				
				firstbases[pos][BaseToIndex(m)][BaseToIndex(alt)] += 1

			for n,base in enumerate(readseq) :
				if base == 'N' : continue

				if (tryread.flag & 16) == 16 : n = readlen-1-n # reverse strand

				# no base changes
				firstbases[n][BaseToIndex(base)][BaseToIndex(base)] += 1

			if (tryread.flag & 16) != 16 : firstforward += 1


		else : # second strand in template

			# Go through each base in read
			readseq = list(tryread.seq)

			for m,pos in mismatches :
				alt = readseq[pos]
				if alt == 'N' : continue

				readseq[pos] = 'N'

				if (tryread.flag & 16) == 16 : pos = readlen-1-pos # reverse strand
				
				lastbases[pos][BaseToIndex(m)][BaseToIndex(alt)] += 1

			for n,base in enumerate(readseq) :
				if base == 'N' : continue

				if (tryread.flag & 16) == 16 : n = readlen-1-n # reverse strand

				# no base changes
				lastbases[n][BaseToIndex(base)][BaseToIndex(base)] += 1

			if (tryread.flag & 16) == 16 : secondreverse += 1


	print 'Discarded reads ', discarded
	print 'First strand forward reads ', firstforward
	print 'Second strand reverse reads ', secondreverse


	try :
		outputfile = open(outputdir+'firstbases.txt', 'w')
	except :
		print 'Unable to create output file.'

	for box in firstbases :
		for row in box :
			outputfile.write('{:.0f} {:.0f} {:.0f} {:.0f}\n'.format(row[0], row[1], row[2], row[3]))

	outputfile.close()

	try :
		outputfile = open(outputdir+'secondbases.txt', 'w')
	except :
		print 'Unable to create output file.'

	for box in lastbases :
		for row in box :
			outputfile.write('{0:.0f} {1:.0f} {2:.0f} {3:.0f}\n'.format(row[0], row[1], row[2], row[3]))

	outputfile.close()

	

# Convert base into index
def BaseToIndex(base) :
	
	if base == 'A' :
		return 0
	elif base == 'C' :	
		return 1
	elif base == 'G' :
		return 2
	elif base == 'T' :
		return 3
	else : # base == 'N'
		print 'Base is N!'
		return -1



# Give MD string from Tags field and cigar as input, and return tuples of (base, position)
# of base mismatches - base here is the one appearing in reference sequence
def BaseMismatches(md, cigar) :

	inserts = [1 for cigartype, cigarlength in cigar if cigartype == 1]

	if md.isdigit() and not inserts :
		return [] 

	mismatches = []
	pos = 0 # sequence index
	number = ''
	i = 0 # md index

	while i < len(md) :
		if md[i] in '0123456789' :
			number += md[i]
			i += 1

		elif md[i] == '^' : # this means deletion - read in all following characters
			if number not in ['0',''] :
				pos += int(number)
				number = ''
			# no need to append as deletion does not appear in seq/qual strings
			i += 1
			while md[i].isalpha() : 
				i += 1			

		else : # symbol should be character, meaning it is base mismatch
			if number not in ['0',''] :
				pos += int(number)
				number = ''
			mismatches.append((md[i], pos))
			pos += 1
			i += 1
			# insertions must be checked from cigar 

	# Check if there are any insertions in Cigar, which will mess up indexing of md
	# and bring additional base-changes 
	cigarpos = 0

	if inserts :
		for cigartype, cigarlength in cigar :

			if cigartype == 1 :

				# Check if there are any items in mismatches coming after this - if
				# yes, then increase their indices by cigarlength
				for m,mpos in mismatches[:] :

					if mpos >= cigarpos :
						mismatches.remove((m,mpos)) 
						mismatches.append((m,mpos+cigarlength))

				cigarpos += cigarlength

			elif cigartype in [3,5] :
				pass	
			else :
				cigarpos += cigarlength

	return mismatches


main()




