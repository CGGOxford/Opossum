# Opossum version 0.2, released on 23.2.2017
# Laura Oikkonen


# Code requires version 0.10.0 of Pysam package
import pysam, itertools, argparse, os, sys

# Store reads temporarily in these dictionaries before processing them
forwardcache = {} # use qname as key for forward reads
reversecache = {} 
positioncache = {} # store end position of each pair, again use qname as key

# Store reads whose mate is unmapped
forwardsingle = {}
reversesingle = {}
positionsingle = {}

# Store number of discarded reads
secondary = 0 # secondary alignments
improper = 0 # improper reads
discardedsingledups = 0 # single reads (whose mates have not been mapped) that are merged with their non-duplicate counterpart

outwardreads = 0 # reads pointing outward and overlapping
duplicatereads = 0 # reads labelled as duplicates that are missing their non-duplicate counterpart
discardedduplicates = 0 # reads labelled as duplicates that are merged with their non-duplicate counterpart

hardclips = 0 # reads containing hard clips
basemismatch = 0 # read pairs having a mismatch in overlapping base sequence
exonmismatch = 0 # read pairs that do not align to same exons
poorquality = 0 # read pairs in which either read has mapping quality lower than set threshold value

# Store number of useful reads
mergedpairs = 0 # merged reads (without any splits)
splitmergedpairs = 0 # merged reads where at least 1 has also been split
separatepairs = 0 # reads pairs that do not overlap
individualreads = 0 # reads that do not have a pair (at least in same chromosome)


def main() :

	global forwardcache
	global reversecache
	global positioncache

	global reversesingle
	global forwardsingle
	global positionsingle

	global secondary
	global improper
	global discardedsingledups
	global outwardreads
	global duplicatereads
	global discardedduplicates
	global hardclips
	global basemismatch 
	global exonmismatch
	global poorquality 
	global mergedpairs 
	global splitmergedpairs 
	global separatepairs 
	global individualreads

	read_unmapped = 0 # unmapped read
	mate_unmapped = 0 # read whose mate is unmapped
	mate_diff_chr = 0 # read whose mate has been mapped to a different chromosome
	mate_junk = 0 # read whose mate has been mapped to junk
	samedir = 0 # both reads in pair have been aligned in the same direction

	parser = argparse.ArgumentParser()
	parser.add_argument('--BamFile', type=str, help='(Required) BAM file name, or entire path if the file is not located in the same directory as where the code is run.')
	parser.add_argument('--MinFlankEnd', default=0, type=int, help='Ignore base-changes closer than MinFlankEnd bases to the end of the read. The corresponding base qualities will be set to 0. Default = 0')
	parser.add_argument('--MinFlankStart', default=0, type=int, help='Ignore base-changes closer than MinFlankStart bases to the start of the read. The corresponding base qualities will be set to 0. Default = 0')
	parser.add_argument('--SoftClipsExist', default='False', choices=['True', 'False'], help='If set to False (default), no soft clips are expected in the reads, and MinFlank and MinFlankStart are applied to all reads. If set to True, then MinFlank and MinFlankStart are only applied to reads with soft clips. The value typically depends on the aligner used, e.g. for TopHat, it should be False, and for Star, it should be True.')
	parser.add_argument('--MapCutoff', default=40, type=int, help='Minimum mapping quality of a read pair. If either read in the pair has a lower mapping quality, both will be discarded. Default = 40')
	parser.add_argument('--ProperlyPaired', default='True', choices=['True', 'False'], help='If set to True (default), only properly paired reads are considered.')
	parser.add_argument('--KeepMismatches', default='False', choices=['True', 'False'], help='If set to False (default), overlapping paired-end reads having at least one base mismatch within the overlap region are discarded.')
	parser.add_argument('--OutFile', type=str, help='(Required) Output BAM file name, or entire path if the file will not be created in the same directory as where the code is run.')
	args = parser.parse_args()

	# Check that user provides at least some arguments
	if not len(sys.argv) > 1 : 
		print 'Please use the --help option to get usage information.'
		return

	elif not args.BamFile :
		print 'Please provide the name of the input BAM file.'
		return

	elif not args.OutFile :
		print 'Please provide the name of the output BAM file.'
		return

	if args.SoftClipsExist == 'True' : SoftClipsExist = True
	else : SoftClipsExist = False


	try :
		samfile = pysam.AlignmentFile(args.BamFile, "rb")
	except :
		print 'Could not open file', args.BamFile
		return

	filename = os.path.split(args.BamFile)

	# Store reads in this dictionary before their mate has also been read
	singlereads = {}

	# Store new read objects in this list
	newreads = []

	# Store single forward reads here until iterator arrives at next starting position
	forwardpos = 0 # position considered for single forward reads
	refpos = 0 # position of forward read that is properly paired

	outputfile = args.OutFile
	# Create output file
	try :
		newfile = pysam.Samfile(outputfile, "wb", template=samfile)
	except :
		print 'Could not create new file'

	contigname = ''

	cont = True
	# Read all entries, one chromosome at a time
	while cont or positioncache or positionsingle :

		# Read in reads until next position with read pairs terminating is reached
		position = 0

		try :
			min1 = min(positioncache.itervalues())
		except :
			min1 = []

		try :
			min2 = min(positionsingle.itervalues())
		except :
			min2 = []
	
		if min1 and min2 :
			minpos = min(min1, min2)
		elif min1 and not min2 :
			minpos = min1
		elif not min1 and min2 :
			minpos = min2
		else :
			minpos = 100000000

		# Read in all reads ending at minpos (or actually minpos-1 due to 0-based indexing)
		while position < minpos :

			try :
				tryread = samfile.next()
				position = tryread.reference_start
				readname = tryread.query_name # Store qname to speed up code
#				print readname # FOR DEBUGGING

			except :
				cont = False
				break

			# Check if chromosome has changed - if yes, write previous entries into file
			if contigname == '' :
				contigname = tryread.rname
				print samfile.get_reference_name(contigname)
			elif contigname == tryread.rname :
				pass
			else : # chromosome has changed
				contigname = tryread.rname

				if contigname == -1 : # read is unmapped
					pass
				else :
					try :
						print samfile.get_reference_name(contigname)
					except :
						print 'Could not print contig name'	

                                # all remaining reads corresponding to this chromosome have now
                                # been read into cache, process them before moving on to the 
				# next chromosome
				remainingpair = set(pos for pos in positioncache.values())
				remainingsingle = set(pos for pos in positionsingle.values())
				remainingpos = remainingpair.union(remainingsingle)

                                for rpos in remainingpos :
                                        addthesereads = ProcessReadsAtPosition(rpos, args.MinFlankEnd, args.MinFlankStart, SoftClipsExist, args.KeepMismatches)
                                        for r in addthesereads :
                                                newreads.append(r)

				# If single forward reads remain, process them :
				if forwardsingle : # Process forwardreads

					fnames = [key for key in forwardsingle.keys()]
					firstread = FindNonduplicate(fnames, False, True)

					# Process the non-duplicate single read 
					for n in CreateNewReads(forwardsingle[firstread], args.MinFlankEnd, args.MinFlankStart, SoftClipsExist) :
						n.query_name = n.query_name + '_single'
						newreads.append(n)

					individualreads += 1
					forwardsingle = {}	

				# Write all entries corresponding to this chromosome to output file
				for n in newreads :
					newfile.write(n)

				newreads = []
				forwardpos = 0 

				# Empty all caches (in theory they should all be empty, except if input bam file contains reads with hard clips)

#				print 'Number of leftover reads: ', len(forwardcache), len(reversecache), len(positioncache), len(singlereads), len(reversesingle), len(forwardsingle), len(positionsingle)
#				print 'Number of reads containing hard clips: ', hardclips

				if forwardcache : forwardcache = {}
				if reversecache : reversecache = {}
				if positioncache : positioncache = {}
				if singlereads : singlereads = {}
				if reversesingle :  reversesingle = {}
				if positionsingle : positionsingle = {}

			# Filter out unmapped reads
			if (tryread.flag & 4) == 4 :
				read_unmapped += 1
				continue

			# Filter out secondary reads
			elif not PrimaryAlignment(tryread.flag, args.ProperlyPaired) : # get rid of secondary alignments
				if (tryread.flag & 256) == 256 : # secondary alignment 
					secondary += 1
				else :
					improper += 1
				continue

			# Filter out reads containing hard clips (since these can involve read pairs with more than two reads, which the code cannot currently handle)
			elif CheckHardClips(tryread.cigar) : 
				hardclips += 1
				continue


			# If read mate is unmapped or doesn't exist, add it to separate cache
			elif (tryread.flag & 8) == 8 or (tryread.flag & 1) != 1 :
				
				if (tryread.flag & 8) == 8: 
					mate_unmapped += 1

				# Check that mapping quality is above threshold
				if tryread.mapping_quality < args.MapCutoff :
					poorquality += 1
					continue

				elif (tryread.flag & 16) == 16 : # check if it is a reverse read
					reversesingle[readname] = tryread
					newpos = tryread.reference_end
					positionsingle[readname] = newpos
					if newpos < minpos : minpos = newpos

				else : # forward read

					if forwardpos == position : 
						forwardsingle[readname] = tryread

					elif forwardpos < position :

						if refpos == forwardpos : 
							discardedsingledups += len(forwardsingle)
							
						elif forwardsingle : # Process forwardreads

							fnames = [key for key in forwardsingle.keys()]
							firstread = FindNonduplicate(fnames, False, True)
							# Process the non-duplicate single read 

							for n in CreateNewReads(forwardsingle[firstread], args.MinFlankEnd, args.MinFlankStart, SoftClipsExist) :
								n.query_name = n.query_name + '_single'
								newreads.append(n)

							individualreads += 1

						forwardsingle = {}
						forwardsingle[readname] = tryread
						forwardpos = position

				continue

			# get rid of reads whose mate has been mapped to a different chromosome
			elif tryread.rname != tryread.rnext :
				if samfile.get_reference_name(tryread.rnext) in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y'] :
					mate_diff_chr += 1
				else :
					mate_junk += 1
				continue

			# Check if mate of read has already been read
			elif readname in singlereads :

				mate = singlereads[readname]

				# Check that mapping quality for both reads is above threshold
				# otherwise discard all of them
				if mate.mapping_quality < args.MapCutoff or tryread.mapping_quality < args.MapCutoff :
					del singlereads[readname]
					poorquality += 2
					continue

				# Check that reads in pair are not pointing outwards (but allow
				# this if they are fully overlapping) - get rid of these cases
				elif ((tryread.flag & 16) == 16 and position <= mate.reference_start and not tryread.reference_end == mate.reference_end) or ((tryread.flag & 16) != 16 and mate.reference_start <= position and not tryread.reference_end == mate.reference_end) :
					outwardreads += 2
					del singlereads[readname]

					continue

				# Check that both reads have not been aligned in the same direction
				elif ((tryread.flag & 16) != 16 and (mate.flag & 16) != 16) or ((tryread.flag & 16) == 16 and (mate.flag & 16) == 16) :
					samedir += 2
					del singlereads[readname]
					continue

				# If reads are fine, then store them into forward / reverse cache
				elif (tryread.flag & 16) == 16 : # check if it is a reverse read 						
					reversecache[readname] = tryread
					
					# Calculate end position of read and store it
					# (If read has soft clips at the end of the read, then
					# adjust end position accordingly)
					softclip = AddClips(tryread.cigar, False)
					newpos = tryread.reference_end + softclip
					positioncache[readname] = newpos
					
					if newpos < minpos : minpos = newpos
						
					# Add read mate into corresponding dictionary
					forwardcache[readname] = mate
					del singlereads[readname]

				else : # it is a forward read
					forwardcache[readname] = tryread

					# Add read mate into corresponding dictionary
					reversecache[readname] = mate

					# Calculate end position of read and store it
					softclip = AddClips(tryread.cigar, False)
					newpos = mate.reference_end + softclip
					positioncache[readname] = newpos
					
					if newpos < minpos : minpos = newpos
						
					del singlereads[readname]

					# Also process single forward reads
					if refpos == position : 
						pass

					elif refpos < position :
						if refpos == forwardpos : 
							discardedsingledups += len(forwardsingle)
							forwardsingle = {}
							
						elif forwardpos == position :
							pass
								
						elif forwardsingle : # Process forwardreads
							fnames = [key for key in forwardsingle.keys()]
							firstread = FindNonduplicate(fnames, False, True)
							# Process the non-duplicate single read 

							for n in CreateNewReads(forwardsingle[firstread], args.MinFlankEnd, args.MinFlankStart, SoftClipsExist) :
								n.query_name = n.query_name + '_single'
								newreads.append(n)

							individualreads += 1

							forwardsingle = {}

						refpos = position


			# If mate has not been read, then store read until it will be
			else :
				singlereads[readname] = tryread

				# If read is forward read, also process single duplicate reads
				if not (tryread.flag & 16) == 16 :

					if refpos == position : 
						pass

					elif refpos < position :

						# Check whether properly paired read and single
						# read start at same position - if yes, you can
						# automatically discard single read
						if refpos == forwardpos : 
							discardedsingledups += len(forwardsingle)
							forwardsingle = {}
							
						elif forwardpos == position :
							pass

						elif forwardsingle : # Process forwardreads
							fnames = [key for key in forwardsingle.keys()]
							firstread = FindNonduplicate(fnames, False, True)
							# Process the non-duplicate single read 
							for n in CreateNewReads(forwardsingle[firstread], args.MinFlankEnd, args.MinFlankStart, SoftClipsExist) :
								n.query_name = n.query_name + '_single'
								newreads.append(n)

							individualreads += 1

							forwardsingle = {}

						refpos = position

				continue

                # Process all single reads and read pairs ending at minpos
                addthesereads = ProcessReadsAtPosition(minpos, args.MinFlankEnd, args.MinFlankStart, SoftClipsExist, args.KeepMismatches) 

                for r in addthesereads :
                        newreads.append(r)
			
	samfile.close()

			
	# If single forward reads remain, process them :
	if forwardsingle : # Process forwardreads

		fnames = [key for key in forwardsingle.keys()]
		firstread = FindNonduplicate(fnames, False, True)
		
		# Process the non-duplicate single read 
		for n in CreateNewReads(forwardsingle[firstread], args.MinFlankEnd, args.MinFlankStart, SoftClipsExist) :
			n.query_name = n.query_name + '_single' 
			newreads.append(n)

		individualreads += 1
		forwardsingle = {}	

	# Write entries in output file
	for n in newreads :
		newfile.write(n)

	newfile.close()

	print 'Number of discarded secondary reads: ', secondary
        if args.ProperlyPaired == 'True' :
                print 'Number of discarded improperly paired reads: ', improper
	print 'Number of discarded unmapped reads: ', read_unmapped
	print 'Number of reads whose mate is unmapped: ', mate_unmapped
	print 'Number of discarded reads whose mate has been mapped to a different chromosome: ', mate_diff_chr
	print 'Number of discarded reads whose mate has been mapped to junk: ', mate_junk
 
	print 'Number of discarded read pairs that are pointing outwards: ', outwardreads
	print 'Number of discarded reads where read and its mate have been mapped in the same direction: ', samedir
        if args.ProperlyPaired == 'True' :
                print 'Number of discarded duplicate reads that are missing non-duplicate counterpart: ', duplicatereads
	print 'Number of duplicate reads that have been merged into their non-duplicate counterpart: ', discardedduplicates
	print 'Number of duplicate single reads that have been merged to their non-duplicate counterpart: ', discardedsingledups

	print 'Number of discarded reads containing hard clips: ', hardclips
        if args.KeepMismatches == 'False' :
                print 'Number of discarded read pairs with base mismatch: ', basemismatch
	print 'Number of discarded read pairs that are not aligned to same exons: ', exonmismatch
	print 'Number of discarded reads / read pairs having too low mapping quality (threshold ' + str(args.MapCutoff) + '): ' + str(poorquality)
	print 'Number of merged read pairs: ', mergedpairs
	print 'Number of split merged read pairs: ', splitmergedpairs
	print 'Number of independently treated read pairs: ', separatepairs
	print 'Number of individual reads: ', individualreads

	print 'Number of leftover reads: ', len(forwardcache), len(reversecache), len(positioncache), len(singlereads), len(reversesingle), len(forwardsingle), len(positionsingle)


# Sort and index output file
	pysam.sort(outputfile, "-o", outputfile, "-T", outputfile[:-4]) # use outputfile[:-4] as prefix for temporary files
	pysam.index(outputfile)

	# Check that output file is not empty
#	try :
#		if os.path.getsize(outputfile) > 0 :

			# Sort and index output file
#			pysam.sort(outputfile, "-o", outputfile, "-T", outputfile[:-4]) # use outputfile[:-4] as prefix for temporary files
#			pysam.index(outputfile)
#		else :
#			print 'Output bam file ', outputfile, ' is empty.'
#	except OSError :
#		print 'Was not able to get size of file ', outputfile



# Go through those reads in cache that end at minpos and process them (de-duplicate, merge).
# Return single reads and/or read pairs that remain after processing.
def ProcessReadsAtPosition(minpos, MinFlankEnd, MinFlankStart, SoftClipsExist, KeepMismatches) :

        global positioncache
        global forwardcache
        global reversecache

        global positionsingle
        global reversesingle

        global individualreads
        global discardedsingledups
        global splitmergedpairs
        global mergedpairs
        global separatepairs

        newreads = []

        # These read pairs end at this position :	
        minpairs = [qname for qname, pos in positioncache.iteritems() if pos == minpos]
        # These single reverse reads end at this position :
        minsingles = [qname for qname, pos in positionsingle.iteritems() if pos == minpos]

        # Only single reads and not read pairs end at this position
        if minsingles and not minpairs :

                # Go through single reverse reads that end at this position and merge them
                firstread = FindNonduplicate(minsingles, False, False)

                # Process the non-duplicate single read
                r = reversesingle[firstread]

                for n in CreateNewReads(r, MinFlankEnd, MinFlankStart, SoftClipsExist) :
                        n.query_name = n.query_name + '_single' 
                        newreads.append(n)

		individualreads += 1
		del positionsingle[firstread], reversesingle[firstread]

                return newreads

        # If both single reads and read pairs exist : 
        # get rid of single reverse reads that end at this position since
        # they are considered to be duplicates			
        for m in minsingles :
                del positionsingle[m], reversesingle[m]
                discardedsingledups += 1

        # Go through these reads and merge each duplicate into corresponding read
        processedreads = FindDuplicatePairs(minpairs)

        # Then go through processedreads
        for p in processedreads :

                r1 = forwardcache[p]
                r2 = reversecache[p]

                # Find out if the reads have any overlap - if they do, then new
                # reads should be constructed base by base.
                if (r1.reference_start <= r2.reference_start and r1.reference_end > r2.reference_start) or (r2.reference_start <= r1.reference_start and r2.reference_end > r1.reference_start) :

                        possiblereads = MergeReads(r1, r2, MinFlankStart, KeepMismatches, SoftClipsExist)
                        tot = len(possiblereads)
                        ind = 1

                        if tot > 1 :
                                splitmergedpairs += 2
                        elif tot == 1 :
                                mergedpairs += 2
                        else :
                                del forwardcache[p], reversecache[p], positioncache[p]
                                continue # all reads have been discarded
				
                        # Additional tag containing number of base qualities
                        for read in possiblereads :
                                read.query_name = read.query_name + '_' + str(ind) + '/' + str(tot) + '_merged'
                                newreads.append(read)
                                ind += 1

                # If read1 and read2 do not overlap, they can be treated independently
                else :

                        for n in CreateNewReads(r1, MinFlankEnd, MinFlankStart, SoftClipsExist) :
                                n.query_name = n.query_name + '_pair' 
                                newreads.append(n)
					
                        for n in CreateNewReads(r2, MinFlankEnd, MinFlankStart, SoftClipsExist) :
                                n.query_name = n.query_name + '_pair' 
                                newreads.append(n)
		
                        separatepairs += 2

		del forwardcache[p], reversecache[p], positioncache[p]


        return newreads

	
# Checks if forward / reverse read has soft clips at the beginning / end and returns the number 
# of clips
def AddClips(cigar, Forward=True) :

	if not Forward : cigar = cigar[::-1]

	for cigartype, cigarlength in cigar :
		if cigartype == 4 :
			return cigarlength
		elif cigartype == 5 :
			print 'there are hard clips'
			return 0
		else :
			return 0


# Checks if read contains any hard clips. If yes, return True; if not, return False.
def CheckHardClips(cigar) :

	for cigartype, cigarlength in cigar :
		if cigartype == 5 :
			return True

	return False



# Returns actual read length, with insertions, soft and hard clippings removed but 
# inserts of type 'N' included
def ReadLength(cigar) :

	l = sum([c[1] for c in cigar if c[0] not in [1,4,5]])
	
	return l


# Identify read pair duplicates - we know that the rightmost base position is the same for 
# all pairs in readpairs, but now we have to check that leftmost position agrees as well
# Merge duplicates with primary reads and return all reads remaining from this processing step
def FindDuplicatePairs(readpairs) :

	global forwardcache
	global reversecache
	global discardedduplicates

	if len(readpairs) == 1 :
		return readpairs

	readgroups = {}
	newpairs = [] # Store read qnames that will be returned

	# Go through reads and store them based on the start position of forward read in
	# a dictionary
	for r in readpairs :

		startpos = forwardcache[r].reference_start - AddClips(forwardcache[r].cigar, True)
		if startpos in readgroups :
			readgroups[startpos].append(r)
		else :
			readgroups[startpos] = [r]

	# For each group in readgroups, decide which of the reads is non-duplicate and then merge 
	# the duplicates with it
	for key in readgroups :

		if len(readgroups[key]) == 1 :
			newpairs.append(readgroups[key][0])
			continue

#		print readgroups[key]

		firstread = FindNonduplicate(readgroups[key], True)
		newpairs.append(firstread)

	return newpairs


# Decide, out of a set of reads that have a common starting (forward reads) or end (reverse reads)
# position, which is the one with highest base quality, and merge the other reads to it.
# If a read has soft clips, don't let them affect the base qualities of the non-duplicate read.
# Return name of the resulting read
def FindNonduplicate(reads, PairedRead, Forward=True) :

	global forwardcache
	global reversecache

	global reversesingle
	global forwardsingle

	global discardedsingledups
	global discardedduplicates

	# Find read with highest base quality sum
	firstread = None
	qualsum = 0			
	
	# First try to find read pair that is properly paired
	if PairedRead :	
		for r in reads :
			
			# check which read has highest mapping quality while being properly paired
			if (forwardcache[r].flag & 2) == 2 :
				qualread = sum(forwardcache[r].query_qualities) + sum(reversecache[r].query_qualities)
				if qualread > qualsum : 
					qualsum = qualread
					firstread = r
	
	if firstread == None : # no reads are properly paired
		for r in reads :
			if PairedRead :
				qualread = sum(forwardcache[r].query_qualities) + sum(reversecache[r].query_qualities)
			elif Forward :
				qualread = sum(forwardsingle[r].query_qualities)
			else : # reverse single read
				qualread = sum(reversesingle[r].query_qualities)

			# check which read has highest mapping quality
			if qualread > qualsum : 
				qualsum = qualread
				firstread = r

	for r in reads :
		if r != firstread :
			if not PairedRead :
				MergeDuplicates(firstread, r, Forward, PairedRead)
				discardedsingledups += 1
			else :
				MergeDuplicates(firstread, r, True, True) # Merge forward reads
				MergeDuplicates(firstread, r, False, True) # Merge reverse reads
				discardedduplicates += 2

	return firstread


# Merge duplicate reads (d is duplicate) such that if they have a base mismatch, mark
# the corresponding base quality to zero (NOTE: if the other read has very poor quality, less
# than 10, at this position, then simply ignore it).
# If a read has soft clips, don't let them affect the base qualities of the non-duplicate read.
# If reads are forward reads, give True as input; for reverse reads, give False as input
# If read is paired, give True as fourth input; for single reads, give False as input
def MergeDuplicates(p, d, ForwardRead, PairedRead) :

	# In order to take into account deletions, soft clippings etc, we have to go through
	# the duplicates base by base

	global forwardcache
	global reversecache
	global positioncache

	global reversesingle
	global forwardsingle
	global positionsingle

	if ForwardRead and PairedRead :
		r1 = forwardcache[p]
		r2 = forwardcache[d]
	elif ForwardRead :
		r1 = forwardsingle[p]
		r2 = forwardsingle[d]
	elif not ForwardRead and PairedRead :
		r1 = reversecache[p]
		r2 = reversecache[d]
	else :
		r1 = reversesingle[p]
		r2 = reversesingle[d]

	firstseq = r1.query_sequence
	secondseq = r2.query_sequence

	# Check if sequences match : # they won't match if there are soft clips
	if firstseq == secondseq : 

		if ForwardRead and PairedRead : del forwardcache[d] # assume that if sequences are same, reads are same too (mapper should treat them in a similar way)
		elif ForwardRead : del forwardsingle[d]
		elif not ForwardRead and PairedRead : del reversecache[d], positioncache[d]
		else : del reversesingle[d], positionsingle[d]

		return 

	newqual = pysam.qualities_to_qualitystring(r1.query_qualities)
	secondqual = pysam.qualities_to_qualitystring(r2.query_qualities)

	if not ForwardRead : # reverse seq and qual
	
		firstseq = firstseq[::-1]
		secondseq = secondseq[::-1]
		newqual = newqual[::-1]
		secondqual = secondqual[::-1]
	# Check that the MD tag actually exists
	try:
		r1_MD = r1.get_tag('MD')
		r2_MD = r2.get_tag('MD')
	except KeyError:
		sys.stderr.write("Error: MD tag doesn't exist in input file\n")
		sys.exit(1)

	# Turn seq and qual into a list
	firstseq = list(firstseq)
	newqual = list(newqual)

	# Before comparing reads base by base, check if there are any misalignments (that is, other
	# read has cigar type 'N' while other does not) and set all following base qualities to 0

	newqual = CheckMisalignments(newqual, r1.cigar, r2.cigar) 

	# If any 'N's remain, they would occur at same positions in both reads so can be ignored
	firstcigar = CigarSequence(r1.cigar, [3]) 
	secondcigar = CigarSequence(r2.cigar, [3]) 

	if not ForwardRead :
		firstcigar = firstcigar[::-1]
		secondcigar = secondcigar[::-1]

	firstpos = 0 # index for first read sequence
	diffpos = 0 # relative index for second read sequence (in terms of first read)
	cigarpos = 0
	diffcigar = 0

	# Check how many bases match from the beginning of the sequences based on MD tag and cigar
	matchingbases = min(MatchingBases(r1_MD, r1.cigar, ForwardRead), MatchingBases(r2_MD, r2.cigar, ForwardRead))
	firstpos += matchingbases
	cigarpos += matchingbases

	# Store sequence lengths into variable in order to speed up code
	firstlen = len(firstseq)
	secondlen = len(secondseq)

	MoreBases = True
	while firstpos < firstlen and (firstpos+diffpos) < secondlen :

		cigartype1 = firstcigar[cigarpos]
		cigartype2 = secondcigar[cigarpos+diffcigar]

		# Read in seq from both reads and compare them
		if cigartype1 == cigartype2 :
			# Check for base mismatch
			if cigartype1 in '01' : # regular base or insertion
				if firstseq[firstpos] == secondseq[firstpos+diffpos] : 
					pass
				# base mismatch
				elif ord(secondqual[firstpos+diffpos]) < 42 :# compare ascii characters
					pass
				elif ord(newqual[firstpos]) < 42 :
					firstseq[firstpos] = secondseq[firstpos+diffpos]
					newqual[firstpos] = secondqual[firstpos+diffpos]
				else : # base mismatch - set quality of corresponding base to 0
					newqual[firstpos] = '!'
				firstpos += 1
				cigarpos += 1
			elif cigartype1 == '2' : # deletion
				cigarpos += 1
			elif cigartype1 == '4' : # soft clip
                                firstpos += 1
                                cigarpos += 1
			else :
				print 'Cigartype should not be this', cigartype1	

		elif cigartype1 in '35' or cigartype2 in '35' : # no useful bases left (since clippings from the beginning of sequences have been removed, these have to occur in the end)
			print 'There should not be clippings or N inserts remaining here'
			break

		elif cigartype1 == '0' :

			if cigartype2 == '1' : # insertion
				if ord(secondqual[firstpos+diffpos]) < 42 :
					pass
				else :	
					newqual[firstpos-1] = '!'
				diffpos += 1
				diffcigar += 1
				cigartype2 = secondcigar[cigarpos+diffcigar]
				while cigartype2 == '1' :
					diffpos += 1
					diffcigar += 1
					cigartype2 = secondcigar[cigarpos+diffcigar]

			elif cigartype2 == '2' : # deletion
				newqual[firstpos] = '!'
				firstpos += 1
				diffpos -= 1
				cigarpos += 1

			else : # cigartype2 == '4' soft clip
				firstpos += 1
				cigarpos += 1

		elif cigartype1 == '1' : # insertion
			
			if cigartype2 == '0' :
				if ord(secondqual[firstpos+diffpos]) < 42 :
					pass
				else :
					newqual[firstpos] = '!'
				firstpos += 1
				cigarpos += 1
				diffpos -= 1
				diffcigar -= 1

			elif cigartype2 == '2' :
				if ord(secondqual[firstpos+diffpos]) < 42 :
					pass
				else :
					newqual[firstpos] = '!'
				firstpos += 1
				diffpos -= 1
				cigarpos += 1
			
			else : # cigartype2 == '4' soft clip
				firstpos += 1
				cigarpos += 1
				diffpos -= 1
				diffcigar -= 1

		elif cigartype1 == '2' :# deletion

			if cigartype2 in ['0','4'] :
				# cannot change quality score of deletion
				cigarpos += 1
				diffpos += 1

			elif cigartype2 == '1' :
				# cannot change quality score of deletion
				diffcigar += 1
				diffpos += 1
 				cigartype2 = secondcigar[cigarpos+diffcigar]
				while cigartype2 == '1' :
					diffpos += 1
					diffcigar += 1
					cigartype2 = secondcigar[cigarpos+diffcigar]

		else : # cigartype1 == '4' soft clip

			if cigartype2 in ['0','4'] :
				firstpos += 1
				cigarpos += 1

			elif cigartype2 == '1' : # insertion
				diffcigar += 1
				diffpos += 1
				cigartype2 = secondcigar[cigarpos+diffcigar]
				while cigartype2 == '1' :
					diffpos += 1
					diffcigar += 1
					cigartype2 = secondcigar[cigarpos+diffcigar]

			else : # cigartype2 == '2' deletion
				firstpos += 1
				diffpos -= 1
				cigarpos += 1
							
	if ForwardRead and PairedRead :		
		forwardcache[p].query_sequence = ''.join(firstseq)
		forwardcache[p].query_qualities = pysam.qualitystring_to_array(''.join(newqual))
		del forwardcache[d]

	elif ForwardRead :
		forwardsingle[p].query_sequence = ''.join(firstseq)
		forwardsingle[p].query_qualities = pysam.qualitystring_to_array(''.join(newqual))
		del forwardsingle[d]

	elif not ForwardRead and PairedRead : # reverse read
		reversecache[p].query_sequence = ''.join(firstseq[::-1])
		reversecache[p].query_qualities = pysam.qualitystring_to_array(''.join(newqual[::-1]))
		del reversecache[d], positioncache[d]

	else : # reverse single read
		reversesingle[p].query_sequence = ''.join(firstseq[::-1])
		reversesingle[p].query_qualities = pysam.qualitystring_to_array(''.join(newqual[::-1]))
		del reversesingle[d], positionsingle[d]
	


# Compare two cigars and if other one has 'N' type while other does not, set the remaining
# base qualities in qual string to zero. Return updated qual string.
def CheckMisalignments(qual, firstcigar, secondcigar) :

	if firstcigar == secondcigar :
		return qual

#	print firstcigar, secondcigar

	# First check whether either cigar has any 'N' type parts
	n1 = sum([1 for cigartype, cigarlength in firstcigar if cigartype == 3])
	n2 = sum([1 for cigartype, cigarlength in secondcigar if cigartype == 3])

	if (n1 + n2) == 0 :
		return qual

	# In case of 'N' type parts, find start and end indices of each bit
	iq = 0 # base quality index (up to which the reads have possible N sections in same places)
	cigarindex1 = 0
	qualpos1 = 0
	cigarindex2 = 0
	qualpos2 = 0
	ntemp1 = 0 # number of N sections in read1
	ntemp2 = 0 # number of N sections in read2

	while (iq < len(qual) or qualpos2 < qualpos1) and (cigarindex2 < len(secondcigar) or qualpos1 < qualpos2 ) :

		if qualpos1 == qualpos2 : # if both have been read up to same point, read in both cigars
			cigartype1, cigarlength1 = firstcigar[cigarindex1]
			cigartype2, cigarlength2 = secondcigar[cigarindex2]

#			print cigartype1, cigarlength1, cigartype2, cigarlength2

			# Read in inserts
			if cigartype1 == 1 :
				iq += cigarlength1 # FIXED!
				cigarindex1 += 1
				continue
			elif cigartype2 == 1 :
				cigarindex2 += 1
				continue

			if cigartype1 == 3 and cigartype2 == 3 :
				ntemp1 += 1
				ntemp2 += 1
				if cigarlength1 != cigarlength2 : # alignment mismatch starts here
					break
				elif ntemp1 == n1 and ntemp2 == n2 : # all inserts have been checked
					iq = len(qual)
					break
				else :
					pass
			elif cigartype1 == 3 or cigartype2 == 3 :
				if cigartype1 == 3 and cigartype2 == 4 :
                                        ntemp1 += 1
                                        if ntemp1 == n1 and ntemp2 == n2 :
                                                iq = len(qual)
                                                break
					elif (cigarindex2+1) == len(secondcigar) :
						iq = len(qual)
						break
                                        else : print 'here is a problem with N inserts and soft clips', qname
                                elif cigartype1 == 4 and cigartype2 == 3 : 
                                        ntemp2 += 1
                                        if ntemp1 == n1 and ntemp2 == n2 :
                                                iq = len(qual)
                                                break
					elif (cigarindex1+1) == len(firstcigar) :
						iq = len(qual)
						break
                                        else : print 'here is a problem with N inserts and soft clips', qname
                                else :
                                        break
			if cigartype1 == 0 or cigartype1 == 4 : # regular base or soft clip
				qualpos1 += cigarlength1
				iq += cigarlength1
#				print iq, len(qual) # FOR DEBUGGING
			else : # deletion
				qualpos1 += cigarlength1
		
			if cigartype2 == 0 : # regular base
				qualpos2 += cigarlength2
			else : # deletion or soft clips
				qualpos2 += cigarlength2

#			if cigarindex1 < len(firstcigar) and cigarindex2 < len(secondcigar) :
			cigarindex1 += 1
			cigarindex2 += 1
#			else :
#				break


		elif qualpos1 < qualpos2 : 

			cigartype1, cigarlength1 = firstcigar[cigarindex1]
			if cigartype1 == 0 or cigartype1 == 4 : # regular base or soft clip
				qualpos1 += cigarlength1
				iq += cigarlength1
			elif cigartype1 == 1 : # insertion
				iq += cigarlength1
			elif cigartype1 == 3 : # N-type insertion; misalignment starts here
				break
			else : # deletion
				qualpos1 += cigarlength1
			cigarindex1 += 1	
		
		else :
			# If rest of first read are soft clips, ignore second read
			if cigarindex1 == len(firstcigar) - 1 : # last item 
				cigartype1, cigarlength1 = firstcigar[cigarindex1]
				if cigartype1 == 4 : # soft clips
					iq = len(qual)
					break

			cigartype2, cigarlength2 = secondcigar[cigarindex2]
			if cigartype2 == 0 : # regular base
				qualpos2 += cigarlength2
			elif cigartype2 == 1 : # insertion
				pass
			elif cigartype2 == 3 : # N-type insertion; misalignment starts here
				iq -= qualpos1 - qualpos2
				break
			else : # deletion
				qualpos2 += cigarlength2
			cigarindex2 += 1

	# Turn qual into a list
	qual = list(qual)

	# Set all base qualities coming after iq to 0 :
	for i in range(iq, len(qual)) :
		qual[i] = '!'

	return qual


# Returns cigar as a string of cigar symbol characters, one for each base etc
# ignoring cigartypes given as parameter (default is [])
def CigarSequence(cigar, ctypes = []) :
	
	seq = ''.join(str(cigartype) * cigarlength for cigartype, cigarlength in cigar if cigartype not in ctypes)

	return seq


# Give MD string from Tags field and cigar as input, and return the number of matching bases 
# at the beginning (FromBeginning = True) or from end (FromBeginning = False).
# Take soft clips into account as well.
def MatchingBases(md, cigar, FromBeginning) :

	if not FromBeginning : 
		md = md[::-1]
		cigar = cigar[::-1]

	i = 0 # md index
	number = ''
	mlen = len(md)

	while i < mlen and md[i] in '0123456789' :
		number += md[i]
		i += 1

	if not number : number = '0'
	if FromBeginning : number = int(number)
	else : number = int(number[::-1])

	# Go through cigar to find possible inserts
	clen = 0	
	for cigartype, cigarlength in cigar :
		if cigartype in [1,3] :
			break	
		elif cigartype == 5 :
			pass
		else :
			clen += cigarlength

	if number < clen :
		number += AddClips(cigar, True) # cigar has been reversed here!

	return min(number, clen)


# Takes flag as input and returns whether alignment is primary (True) or not (False).
# If additional input parameter is set to True, then also check whether read is properly
# paired (definition of properly paired depends on mapper).
def PrimaryAlignment(flag, properlypaired) :

	if properlypaired == 'True' : # parameter here is string and not boolean
		if (flag & 2) == 2 and not (flag & 256) == 256 : # check if properly paired and primary alignment
			return True
		else :
			return False

	elif (flag & 256) == 256 : # check if secondary alignment
		return False

	else :
		return True		



# Create new read(s) by splitting input row as necessary
def CreateNewReads(row, MinFlankEnd, MinFlankStart, SoftClipsExist) :

	# Trim cigar, sequence, and quality in case there are soft / hard clippings 
	newcigar, start, end = TrimCigar(row.cigar)

	#Split the reads if required (i.e. cigarType == 3) and create new object for each new read
	splits = NumberOfSplits(newcigar)

	# Check whether read is forward / reverse
	if (row.flag & 16) == 16 : ForwardRead = False
	else : ForwardRead = True

	# Check whether read originally had soft clips or not
	if SoftClipsExist :
		if newcigar == row.cigar :
			MinFlankEnd = 0
			MinFlankStart = 0
		else :
		# check if soft clips are at the beginning or end of the read
			cigarlist = [cigartype for cigartype, cigarlength in row.cigar]

			if ForwardRead :
				if cigarlist[0] != 4 :
					MinFlankStart = 0
				if cigarlist[-1] != 4 :
					MinFlankEnd = 0
			else : # reverse read
				if cigarlist[0] != 4 :
					MinFlankEnd = 0
				if cigarlist[-1] != 4 :
					MinFlankStart = 0

	# Check if the MD tag actually exists
	try:
		row_MD = row.get_tag('MD')
	except KeyError:
		sys.stderr.write("Error: MD tag doesn't exist in input file\n")
		sys.exit(1)

	mincut, maxcut = FindCutoffIndices(newcigar, MinFlankEnd, MinFlankStart, ForwardRead)

	# If read does not need to be split, just copy previous read object to new list
	if splits == 0 and newcigar == row.cigar :

		newqual = DiscardVariants(pysam.qualities_to_qualitystring(row.query_qualities), row_MD, newcigar, mincut, maxcut)
		newrow = CreateReadObject(row, row.query_sequence, newqual, newcigar, row.reference_start)
		yield newrow 

	elif splits == 0 : # cigar, seq, qual have been trimmed

		newseq = row.query_sequence[start:end]
		newqual = pysam.qualities_to_qualitystring(row.query_qualities)[start:end]

		newqual = DiscardVariants(newqual, row_MD, newcigar, mincut, maxcut)
		newrow = CreateReadObject(row, newseq, newqual, newcigar, row.reference_start)
		yield newrow 
		
	else :
		if newcigar != row.cigar : 
			row.cigar = newcigar
			q = pysam.qualities_to_qualitystring(row.query_qualities)[start:end]
			row.query_sequence = row.query_sequence[start:end]
			row.query_qualities = pysam.qualitystring_to_array(q)

		newreads = []
		i = 0 #index for sequence and quality
		k = 0 #elements in cigar
		position = row.reference_start #position of read

		for part in range(splits+1) :

			readname = row.query_name + "_" + str(part+1) + '/' + str(splits+1)
			newrow, newpos, newk, newi = CreateSplitRead(row, position, k, i)
			position = newpos

			# Check that new read does not consist only of deletions - this would not 
			# make sense and would result in empty quality string
			if len(newrow.cigar) == 1 :
				cigartype, cigarlength = newrow.cigar[0]
				if cigartype == 2 : # deletion
					k = newk
					i = newi
					continue

			newrow.query_qualities = pysam.qualitystring_to_array(DiscardVariants(newrow.qual, row_MD, row.cigar, mincut, maxcut, i))
			k = newk
			i = newi 
			newrow.query_name = readname
			yield newrow 


# Create a new split read as part X out of Y splits from original read
def CreateSplitRead(row, position, k, i) :

	startread = position
	newcigar = []
	newstring = ''
	oldseq = row.query_sequence
	oldqual = pysam.qualities_to_qualitystring(row.query_qualities)
	newqual = ''

	while k < len(row.cigar) :
	
		cigartype, cigarlength = row.cigar[k]
						
		# cigartype 1 stands for insertion
		# in this case, position should be kept as it is
		if cigartype == 1 :
			newstring +=  oldseq[i:(cigarlength+i)]
			newqual += oldqual[i:(cigarlength+i)]
			i += cigarlength
			newcigar.append((cigartype, cigarlength))	
			k += 1	

		# cigartype 2 stands for deletion
		# in this case, nothing should be added to the sequence
		elif cigartype == 2 :
			position += cigarlength
			newcigar.append((cigartype, cigarlength))
			k += 1

		# cigartype 4 (soft clipping) should be skipped - bases present in original seq.
		# Position does not take this into account
		elif cigartype == 4 :
			print 'Cigartype should not be 4!'
			i += cigarlength
			k += 1

		# cigartype 5 (hard clipping) does not include bases in original seq
		elif cigartype == 5 :
			print 'Cigartype should not be 5!'
			k += 1

		# cigartype 3 should be skipped 
		elif cigartype != 3 :
			newstring +=  oldseq[i:(cigarlength+i)]
			newqual += oldqual[i:(cigarlength+i)]
			i += cigarlength
			position += cigarlength 
			newcigar.append((cigartype, cigarlength))
			k += 1	

		else :
			position += cigarlength
			k += 1
			break	

	a = CreateReadObject(row, newstring, newqual, newcigar, startread)

	return a, position, k, i


# Create read object based on information from original read and new sequence, quality, cigar
# and position.
# Return new read object.
def CreateReadObject(read, newseq, newqual, newcigar, startread, basetag=[]) :

	a = pysam.AlignedSegment()
	a.query_name = read.query_name
	a.query_sequence = newseq
	a.query_qualities = pysam.qualitystring_to_array(newqual)
	a.cigar = newcigar
	a.reference_start = startread

	# If (Star) mapper has assigned a value of 255 to mapping quality,
	# change it to 50
	mapqual = read.mapping_quality 
	if mapqual == 255 :
		a.mapping_quality = 50
	else :
		a.mapping_quality = mapqual

	a.reference_id = read.reference_id

	# If read has RG read group tag, keep it
	try :
		r_RG = read.get_tag('RG')
		a.tags = ()
		a.set_tag('RG', r_RG)
	except :
		a.tags = ()

	a.next_reference_id = -1
	a.next_reference_start = -1
	a.template_length = 0
	a.flag = UpdateFlag(read.flag)

	return a


# Goes through the given cigar list and reports how many times cigarType is equal to 3
# Optional field is finalpos, which is cutoff position for counting the splits (relative to start pos)
# Default value for finalpos is something very large
def NumberOfSplits(cigarlist, finalpos=100000) :

	number = 0
	pos = 0

	if type(cigarlist) is not list:
		print 'Strange CIGAR:', str(cigarlist)
		return 0

	for (cigartype, cigarlength) in cigarlist :

		if pos >= finalpos :
			break

		if cigartype == 3 :
			number += 1

		if cigartype not in [1,4,5] :
			pos += cigarlength

	return number


# Checks if read has soft or hard clippings and removes them from cigar. Also returns
# new indices for sequence / quality string so that it can also be trimmed.
def TrimCigar(cigar) :

	clippings = [1 for cigartype, cigarlength in cigar if cigartype in [4,5]]

	if not clippings :
		return cigar, 0, None

	start = 0
	end = 0

	# soft and hard clippings can either be in the beginning or end of a read
	for cigartype, cigarlength in cigar[:] :

		if cigartype == 4 : # soft clipping
			start = cigarlength
			cigar.remove((cigartype, cigarlength))

		elif cigartype == 5 : # hard clipping
			cigar.remove((cigartype, cigarlength))

		else :
			break

	for cigartype, cigarlength in reversed(cigar[:]) :

		if cigartype == 4 : # soft clipping
			end -= cigarlength
			cigar.remove((cigartype, cigarlength))

		elif cigartype == 5 : # hard clipping
			cigar.remove((cigartype, cigarlength))

		else :
			break

	return cigar, start, end if not end == 0 else None


			

# Give as input the quality string, md from Tags field, cigar, and mincut and maxcut,
# and ipos (starting index of current read in terms of quality sequence)
# and adjustbase (only needed if two reads
# are being merged, and second read starts later than first read).
# Md tells whether there are any base-changes lying closer or equal to mincut from
# read start edge, or at least maxcut away from same edge. These base-changes will be 
# discarded, that is, their corresponding base
# quality will be set to 0 (or ! in ascii format).
# Returns updated quality string.
def DiscardVariants(qual, md, cigar, mincut, maxcut, ipos=0, adjustbase=0) :

	inserts = [1 for cigartype, cigarlength in cigar if cigartype == 1]

	if md.isdigit() and not inserts :
		return qual 

	# Turn qual into a list
	qual = list(qual)

	#First go through MD and mark all variant positions to a list
	variantpos = [] # corresponding to indices in qual string
	pos = 0 # qual string index
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
			i += 1

			while md[i].isalpha() : 
				i += 1	

		else : # symbol should be character, meaning it is variant
			if number not in ['0',''] :
				pos += int(number)
				number = ''

			variantpos.append(pos)
			pos += 1
			i += 1
			# character may also be insertion, this must be checked from cigar 

	# Check if there are any insertions in Cigar, which will mess up indexing of md
	# and bring additional base-changes 
	cigarpos = 0
	if inserts :
		for cigartype, cigarlength in cigar :

			if cigartype == 1 :
				# Check if there are any items in variantpos coming after this - if there
				# are, then increase their indices by cigarlength
				for var in variantpos[:] :
					if var >= cigarpos :
						variantpos.remove(var) 
						variantpos.append(var+cigarlength)
						

				for i in range(cigarlength) : 
					variantpos.append(cigarpos) 
					cigarpos += 1

			elif cigartype in [2,3,4,5] :
				pass	
			else :
				cigarpos += cigarlength

	# Finally replace base quality of these variants to 0 (which is !)
	for var in variantpos :
		if var >= ipos and var <= (ipos+len(qual)-1) and (var < mincut or var > maxcut) : # qual may be shorter than original read
			qual[var-ipos+adjustbase] = '!'
		
	qual = ''.join(qual)

	return qual


# Find out what are the min and max cutoff indices for quality string if we want
# to discard variants that are closer than MinFlank bases away from read edges.
# They are defined with respect to start edge of the read.
# Indices are 0-based.
# Assume that soft and hard clippings have been discarded or
# otherwise taken into account. Deletions will not affect cutoff value.
def FindCutoffIndices(cigar, MinFlankEnd, MinFlankStart, ForwardRead) :

	cigarstring = CigarSequence(cigar, [3,5])

	if ForwardRead :
		firstcut = MinFlankStart
		lastcut = MinFlankEnd
	else :
		firstcut = MinFlankEnd
		lastcut = MinFlankStart

	mincut = 0
	maxcut = 0

	i = 1	
	for c in cigarstring :

		if i < firstcut :
			if c == '1' : # indicates insert
				mincut += 1
			elif c == '2' : # indicates deletion
			#	i += 1
				pass
			elif c in '345' : # indicates soft or hard clipping / N-type insert
				print 'There should not be a clip or N insert here anymore'
			else :
				mincut += 1
				i += 1
		else : break

	j = 1
	for d in cigarstring[::-1] :

		if j < lastcut :
			if d == '1' : # indicates insert
				maxcut += 1
			elif d == '2' : # indicates deletion
				pass
			#	j += 1
			elif d in '345' : # indicates soft or hard clipping / N-type insert
				print 'There should not be a clip or N insert here anymore'
			else :
				maxcut += 1
				j += 1

		else : break

	l = sum([cigarlength for cigartype, cigarlength in cigar if cigartype not in [2,3,5]])

	return mincut, l-maxcut-1


#Takes flag as input and returns whether read is forward or reverse (0 / 16) and/or duplicate (1024)
def UpdateFlag(flag) :

	if (flag & 16) == 16 : return 16
	else : return 0


# Returns the cigar length when cigar is expressed base by base,
# up to the index k of cigarlist (default is a very large index)
def ReadCigarLength(cigar, k=100) :

	if k < len(cigar) :
		cigarlist = cigar[:(k+1)]
	else :
		cigarlist = cigar

	l = sum([c[1] for c in cigarlist])
		
	return l


# Input cigar and pos denoting index within a cigar transformed into a string sequence.
# Return what is the index in the original cigar and how many extra positions from the beginning
# of this index is the current position.
def CigarIndex(cigar, pos) :

	position = 0
	k = 0 # cigar index
	addpos = pos # extra positions

	for cigartype, cigarlength in cigar :

		position += cigarlength

		if position < pos : 
			k += 1
			addpos -= cigarlength
		else : return k, addpos



# Converts string of cigar symbol characters into the format used by Pysam
def CigarConcatenate(cigarseq) :

	if not cigarseq :
		return []

	return [(int(k), len(list(g))) for k, g in itertools.groupby(cigarseq)]


# Merge and split (if necessary) the two reads 
# Return new reads
def MergeReads(r1, r2, MinFlankStart, keepit, SoftClipsExist) :

	global basemismatch
	global exonmismatch

	if keepit == 'True' : # keepit is a string here and not boolean
		KeepMismatches = True	
	else :
		KeepMismatches = False

	startbase = min(r1.reference_start, r2.reference_start)
	laterbase = max(r1.reference_start, r2.reference_start)
	endbase = min(r1.reference_end, r2.reference_end)
	finalbase = max(r1.reference_end, r2.reference_end)	
	
	if startbase == r1.reference_start : 
		firstread = r1
		secondread = r2
	else : 
		firstread = r2
		secondread = r1

	# Trim cigar, sequence, and quality in case there are soft / hard clippings 
	newcigar1, start1, end1 = TrimCigar(firstread.cigar)
	newcigar2, start2, end2 = TrimCigar(secondread.cigar)

	# Check if there were any soft clips in original read
	if SoftClipsExist :
		if newcigar1 == firstread.cigar and newcigar2 == secondread.cigar :
			MinFlankStart = 0

	mincut, maxcutno = FindCutoffIndices(newcigar1, MinFlankStart, MinFlankStart, True)
	mincutno, maxcut = FindCutoffIndices(newcigar2, MinFlankStart, MinFlankStart, False)

	firstcigar = CigarSequence(newcigar1) 
	secondcigar = CigarSequence(newcigar2)

	firstseq = firstread.query_sequence[start1:end1]
	firstqual = pysam.qualities_to_qualitystring(firstread.query_qualities)[start1:end1] # retrieve quality strings - this will speed up code

	secondseq = secondread.query_sequence[start2:end2]
	secondqual = pysam.qualities_to_qualitystring(secondread.query_qualities)[start2:end2]

	try:
		firstmd = firstread.get_tag('MD') # retrieve MD tags - this also speeds up code
		secondmd = secondread.get_tag('MD')
	except KeyError:
		sys.stderr.write("Error: MD tag doesn't exist in input file\n")
		sys.exit(1)

	firstpos = 0
	cigarpos = 0
	possiblereads = [] # Store here reads before adding them to newreads

	# While there is no overlap, check if first read needs to be split
	# If not, then just add the entire first part to the sequence
	moresplits = NumberOfSplits(newcigar1, laterbase - startbase)
	i = 0 # index of first quality base of new read
	k = 0 # elements in cigar
	position = startbase # position of read

	if moresplits > 0 and newcigar1 != firstread.cigar :
		firstread.query_sequence = firstseq
		firstread.query_qualities = pysam.qualitystring_to_array(firstqual)
		firstread.cigar = newcigar1
			
	for n in range(moresplits) :
		newrow, newpos, newk, newi = CreateSplitRead(firstread, position, k, i)
		position = newpos

                # Check that new read does not consist only of deletions - this would not  
	        # make sense and would result in empty quality string
	        if len(newrow.cigar) == 1 :
		        cigartype, cigarlength = newrow.cigar[0]
		        if cigartype == 2 : # deletion
			       k = newk
			       i = newi
			       continue

		# Check if mincut extends to this new read
		newrow.query_qualities = pysam.qualitystring_to_array(DiscardVariants(newrow.qual, firstmd, newcigar1, mincut, 100000, i))

		k = newk
		i = newi		
		possiblereads.append(newrow)

	firstpos += i
				
	# use k-1 here because we want to know cigarlength up to the 
	# k-th cigar term
	cigarpos += ReadCigarLength(newcigar1, k-1)

	# Start new sequence to combine the two reads
	newseq = ''
	newqual = ''
	newcigar = ''
	startread = position # base position indicating start of a read

	# If position < laterbase then add this bit to the new sequence.
	# This bit cannot contain any 'N' inserts, otherwise position would 
	# have jumped to at least laterbase.
				
	while position < laterbase :

		cigartype, cigarlength = newcigar1[k]

		if position + cigarlength < laterbase or cigartype in [1,4,5] : # read in whole cigar block - position will not increase in case of insertions or clippings
			lastpos = firstpos + cigarlength

		else : # read only up to laterbase - 1
			lastpos = firstpos + laterbase - position 
			cigarlength = laterbase - position 

		if cigartype == 1 : # insertion
			newseq += firstseq[firstpos:lastpos]
			newqual += firstqual[firstpos:lastpos]
			newcigar += str(cigartype) * cigarlength	
			firstpos += cigarlength	

		elif cigartype == 2 : # deletion
			position += cigarlength
			newcigar += str(cigartype) * cigarlength

		elif cigartype in [3,4,5] :
			print 'Cigartype should not be 3 / 4 / 5 here.'
			position += cigarlength

		else :
			newseq +=  firstseq[firstpos:lastpos]
			newqual += firstqual[firstpos:lastpos]
			position += cigarlength 
			newcigar += str(cigartype) * cigarlength
			firstpos += cigarlength

		cigarpos += cigarlength
		k += 1

	# If position was actually > laterbase, it means that the overlapping
	# part of the first base starts with 'N' sequence. We cannot have  
	# a base for one read and 'N' for the other, therefore, remove
	# both reads
	if position > laterbase :
		exonmismatch += 2
		return []
					
	# Now go through bases one by one.

	basediff = firstpos  # define seq index for secondread with respect to firstread
	cigardiff = cigarpos  # same with cigar string

	# This is the part where the two reads are overlapping
	while position < endbase :

		nextc1 = firstcigar[cigarpos]
		nextc2 = secondcigar[cigarpos-cigardiff]
		
		if nextc1 == nextc2 :

			# Read in sequence of bases.
			if nextc1 == '0' : # normal base

				if firstseq[firstpos] == secondseq[firstpos-basediff] : # bases match
					if secondqual[firstpos-basediff] != '!' : # base quality for second read has not been set to 0
						newseq += firstseq[firstpos]
						newqual += max(firstqual[firstpos], secondqual[firstpos-basediff]) # choose the highest base quality score
					else :
						newseq += firstseq[firstpos]
						newqual += '!'
				elif KeepMismatches : # bases do not match
					if ord(secondqual[firstpos-basediff]) < 42 :
						newseq += firstseq[firstpos]
						newqual += firstqual[firstpos]
					elif ord(firstqual[firstpos]) < 42 :
						newseq += secondseq[firstpos-basediff]
						newqual += secondqual[firstpos-basediff]
					else :
						newseq += firstseq[firstpos]
						newqual += '!'
				else :
					basemismatch += 2
					return []

				firstpos += 1
				position += 1
				newcigar += nextc1
				cigarpos += 1

			elif nextc1 == '1' : # insertion
				if firstseq[firstpos] == secondseq[firstpos-basediff] : # bases match
					if secondqual[firstpos-basediff] != '!' :
						newseq += firstseq[firstpos]
						newqual += max(firstqual[firstpos], secondqual[firstpos-basediff]) # choose the highest base quality score
					else :
						newseq += firstseq[firstpos]
						newqual += '!'
				elif KeepMismatches : # bases do not match
					if ord(secondqual[firstpos-basediff]) < 42 :
						newseq += firstseq[firstpos]
						newqual += firstqual[firstpos]
					elif ord(firstqual[firstpos]) < 42 :
						newseq += secondseq[firstpos-basediff]
						newqual += secondqual[firstpos-basediff]
					else :
						newseq += firstseq[firstpos]
						newqual += '!'
				else :
					basemismatch += 2
					return []

				firstpos += 1
				newcigar += nextc1
				cigarpos += 1

			elif nextc1 == '2' : # deletion
				position += 1
				newcigar += nextc1
				cigarpos += 1	

			# Both reads are cut here. Make previous sequence into a 
			# new read and start next sequence
			elif nextc1 == '3':

				# bases at the beginning of the read should be checked 
				if i <= mincut : 
					newqual = DiscardVariants(newqual, firstmd, newcigar1, mincut, 100000, i) 
				# new read starts before secondread
				if i >= (maxcut + basediff) and startread < secondread.reference_start : 
					newqual = DiscardVariants(newqual, secondmd, newcigar2, 0, maxcut, 0, basediff+morepos-i) 
				elif i >= (maxcut + basediff) :
					newqual = DiscardVariants(newqual, secondmd, newcigar2, 0, maxcut, i-basediff)

				a = CreateReadObject(firstread, newseq, newqual, CigarConcatenate(newcigar), startread)
				possiblereads.append(a)

				newseq = ''
				newqual = ''
				newcigar = ''
				i = firstpos

				# Move both reads forward until at least the other 
				# has sequence again
				# THIS STEP PROBABLY TAKES A LONG TIME			
				while firstcigar[cigarpos] == '3' and secondcigar[cigarpos-cigardiff] == '3' :	
					cigarpos += 1
					position += 1
					
				startread = position

			else :
				print 'Strange cigartype'

		# Reads have different cigartypes
		elif nextc1 == '3' or nextc2 == '3': # inconsistent alignment between the reads. 							     # Discard them.
			exonmismatch += 2
			return []

		elif nextc1 in '45' : 
			print 'Soft or hard clipping here'

		elif KeepMismatches and nextc1 == '0' : 

			if nextc2 == '1' : # insertion	
				if ord(secondqual[firstpos-basediff]) < 42 :
					pass
				elif ord(firstqual[firstpos-1]) < 42 : # add insertion
					newseq += secondseq[firstpos-basediff]
					newqual += secondqual[firstpos-basediff]
					newcigar += nextc2
				else :
					newqual = newqual[:-1] + '!'

				basediff -= 1
				cigardiff -= 1

				while secondcigar[cigarpos-cigardiff] == '1': # read in all insertions
					basediff -= 1
					cigardiff -= 1

			else : # nextc2 == '2' : # deletion
				newseq += firstseq[firstpos]
				newqual += '!' # deletion does not have base quality
				newcigar += nextc1
				firstpos += 1
				cigarpos += 1
				position += 1
				basediff += 1

		elif KeepMismatches and nextc1 == '1' :

			while firstcigar[cigarpos] == '1' : # read in all insertions
				newseq += firstseq[firstpos]
				newqual += '!'
				newcigar += firstcigar[cigarpos]
				firstpos += 1
				cigarpos += 1
				basediff += 1
				cigardiff += 1

		elif KeepMismatches and nextc1 == '2' : # deletion

			if nextc2 == '0' : # normal base

				newseq += secondseq[firstpos-basediff] #new line
				newqual += '!' #new line
				newcigar += nextc2 #modified line
				cigarpos += 1
				basediff -= 1
				position += 1

			elif nextc2 == '1' : # insertion
				newqual = newqual[:-1] + '!'
				basediff -= 1
				cigardiff -= 1

				while secondcigar[cigarpos-cigardiff] == '1': # read in all insertions
					basediff -= 1
					cigardiff -= 1

		else : # KeepMismatches is False, so get rid of reads

			basemismatch += 2
			return []

	# Check if both reads have reached their end :
	if cigarpos == ReadCigarLength(newcigar1) and cigarpos-cigardiff == ReadCigarLength(newcigar2) :
		pass

	# Check here which read has now reached its end:
	elif cigarpos == ReadCigarLength(newcigar1) : # First read has reached its end (this equation would not otherwise hold but cigarpos is +1 than it actually should as the result of the last iteration step)

		k, addpos = CigarIndex(newcigar2, cigarpos-cigardiff)
	
		while position < finalbase :
			
			cigartype, cigarlength = newcigar2[k]
			cigarlength -= addpos

			if cigartype == 1 : # insertion
				newseq += secondseq[(firstpos-basediff):(firstpos-basediff+cigarlength)]
				newqual += secondqual[(firstpos-basediff):(firstpos-basediff+cigarlength)]
				newcigar += str(cigartype) * cigarlength	
				firstpos += cigarlength	

			elif cigartype == 2 : # deletion
				position += cigarlength
				newcigar += str(cigartype) * cigarlength

			elif cigartype in [4,5] : # soft or hard clipping
				pass

			elif cigartype == 3 : # Create new read object
				if i <= mincut : # bases at the beginning of the read should be checked 
					newqual = DiscardVariants(newqual, firstmd, newcigar1, mincut, 100000, i) 

				if i >= (maxcut + basediff) and startread < secondread.reference_start : # new read starts before secondread
					newqual = DiscardVariants(newqual, secondmd, newcigar2, 0, maxcut, 0, basediff+morepos-i)
				elif i >= (maxcut + basediff) :
					newqual = DiscardVariants(newqual, secondmd, newcigar2, 0, maxcut, i-basediff) 
				a = CreateReadObject(secondread, newseq, newqual, CigarConcatenate(newcigar), startread)		
				possiblereads.append(a)

				newseq = ''
				newqual = ''
				newcigar = ''
				i = firstpos

				position += cigarlength
				startread = position

			else :
				newseq += secondseq[(firstpos-basediff):(firstpos-basediff+cigarlength)]
				newqual += secondqual[(firstpos-basediff):(firstpos-basediff+cigarlength)]
				position += cigarlength 
				newcigar += str(cigartype) * cigarlength
				firstpos += cigarlength

			cigarpos += cigarlength
			k += 1
			addpos = 0

	else : # Second read has reached its end
	
		k, addpos = CigarIndex(newcigar1, cigarpos)

		while position < finalbase :

			cigartype, cigarlength = newcigar1[k]
			cigarlength -= addpos

			if cigartype == 1 : # insertion
				newseq += firstseq[firstpos:(firstpos+cigarlength)]
				newqual += firstqual[firstpos:(firstpos+cigarlength)]
				newcigar += str(cigartype) * cigarlength	
				firstpos += cigarlength	

			elif cigartype == 2 : # deletion
				position += cigarlength
				newcigar += str(cigartype) * cigarlength

			elif cigartype in [4,5] : # soft or hard clipping
				pass

			elif cigartype == 3 : # Create new read object
				
				newqual = DiscardVariants(newqual, firstmd, newcigar1, 0, maxcutno, i) 
				if i <= mincut : # bases at the beginning of the read should be checked
					newqual = DiscardVariants(newqual, firstmd, newcigar1, mincut, 100000, i) 
				a = CreateReadObject(firstread, newseq, newqual, CigarConcatenate(newcigar), startread)		
				possiblereads.append(a)

				newseq = ''
				newqual = ''
				newcigar = ''
				i = firstpos

				position += cigarlength
				startread = position

			else :
				newseq += firstseq[firstpos:(firstpos+cigarlength)]
				newqual += firstqual[firstpos:(firstpos+cigarlength)]
				position += cigarlength 
				newcigar += str(cigartype) * cigarlength
				firstpos += cigarlength

			cigarpos += cigarlength
			k += 1
			addpos = 0


	# Write down last read
	if newseq != '' and newcigar:
			
		if i <= mincut : # bases at the beginning of the read should be checked 
			newqual = DiscardVariants(newqual, firstmd, newcigar1, mincut, 100000, i) 
	
		if startread < secondread.reference_start : # new read starts before secondread	
			newqual = DiscardVariants(newqual, secondmd, newcigar2, 0, maxcut, 0, basediff-i) 
		else :	
			newqual = DiscardVariants(newqual, secondmd, newcigar2, 0, maxcut, i-basediff) 

		a= CreateReadObject(firstread, newseq, newqual, CigarConcatenate(newcigar), startread)
		
		possiblereads.append(a)

	else :
		print 'empty read', firstread.query_name

	return possiblereads


if __name__ == "__main__":
	main()
