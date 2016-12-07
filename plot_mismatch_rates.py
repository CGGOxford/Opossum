# Laura Oikkonen 7.12.2016

import matplotlib.pyplot as plt
import numpy as np

# Plot mismatch rates for the four first base positions starting 
# from input txt file that has been generated with 
# compute_mismatch_rates.py. 


def main() :

#################################################################
	# FILL THESE IN
	# Input txt file
	inputname='firstbases.txt'
	# Read length
	readlen=76
	# Plot title name
	ptitle='First strand nucleotides'
	# Output eps file name
	outputname='firststrand.eps'
#################################################################

	try :
		inputfile = open(inputname, 'r')
	except :
		print 'Could not open file'

	firstbases = np.zeros((readlen,4,4))
	a = 0 # index of 'box'
	b = 0 # index of reference base

	for row in inputfile :
		row = row.rstrip()
		row = row.split()
		row = [int(x) for x in row]
		basesum = sum(row)
		firstbases[a][b][0:4] = [1.0 * i / basesum if n!=b else 0 for n,i in enumerate(row)]
		
		if b == 3 : 
			b = 0
			a += 1	
		else : 
			b += 1

	inputfile.close()

	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, sharex='col', sharey='row')

	ind = np.arange(4) # x locations for the groups
	width = 0.35 # width of the bars

	abases0 = firstbases[0,:,0]
	cbases0 = firstbases[0,:,1]
	gbases0 = firstbases[0,:,2]
	tbases0 = firstbases[0,:,3]

	ax1.bar(ind, abases0, width, color='r')
	ax1.bar(ind, cbases0, width, color='y', bottom=abases0)
	ax1.bar(ind, gbases0, width, color='g', bottom=abases0+cbases0)
	ax1.bar(ind, tbases0, width, color='b', bottom=abases0+cbases0+gbases0)
	ax1.set_ylabel('% of error nucleotides')
	ax1.set_title('First position')
	ax1.set_xlim(-width/2., len(ind)-width)


	abases1 = firstbases[1,:,0]
	cbases1 = firstbases[1,:,1]
	gbases1 = firstbases[1,:,2]
	tbases1 = firstbases[1,:,3]

	ax2.bar(ind, abases1, width, color='r')
	ax2.bar(ind, cbases1, width, color='y', bottom=abases1)
	ax2.bar(ind, gbases1, width, color='g', bottom=abases1+cbases1)
	ax2.bar(ind, tbases1, width, color='b', bottom=abases1+cbases1+gbases1)
	ax2.set_title('Second position')
	ax2.set_xlim(-width/2., len(ind)-width)

	abases2 = firstbases[2,:,0]
	cbases2 = firstbases[2,:,1]
	gbases2 = firstbases[2,:,2]
	tbases2 = firstbases[2,:,3]

	ax3.bar(ind, abases2, width, color='r')
	ax3.bar(ind, cbases2, width, color='y', bottom=abases2)
	ax3.bar(ind, gbases2, width, color='g', bottom=abases2+cbases2)
	ax3.bar(ind, tbases2, width, color='b', bottom=abases2+cbases2+gbases2)
	ax3.set_xlabel('Reference base')
	ax3.set_ylabel('% of error nucleotides')
	ax3.set_title('Third position')
	ax3.set_ylim((0,0.1))
	ax3.set_xlim(-width/2., len(ind)-width)
	ax3.set_xticks(ind+width/2.)
	ax3.set_xticklabels(['A','C','G','T'])
	
	abases3 = firstbases[3,:,0]
	cbases3 = firstbases[3,:,1]
	gbases3 = firstbases[3,:,2]
	tbases3 = firstbases[3,:,3]

	ax4.bar(ind, abases3, width, color='r', label='A')
	ax4.bar(ind, cbases3, width, color='y', bottom=abases3, label='C')
	ax4.bar(ind, gbases3, width, color='g', bottom=abases3+cbases3, label='G')
	ax4.bar(ind, tbases3, width, color='b', bottom=abases3+cbases3+gbases3, label='T')
	ax4.set_title('Fourth position')
	ax4.set_xlim(-width/2., len(ind)-width)

	ax1.set_ylim((0,0.1))

	plt.xticks(ind+width/2., ('A','C','G','T'))
	plt.xlabel('Reference base')
	plt.legend()
	
	plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['FreeSans']	
	plt.rcParams.update({'font.size': 16})

	plt.suptitle(ptitle, size=20)

	plt.savefig(outputname, dpi=96)



main()
