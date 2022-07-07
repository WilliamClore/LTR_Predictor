#LTR_Predictor requires Biopython and Matplotlib to be loaaded to run

from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
import pylab
from Bio.SeqUtils import lcc
import math

#user promted to enter FASTA file target by LTR_Predictor
file = input("Enter Target File Name: ")

#user prompted to enter yes or no to create dotplots for LTR candidates
dt = input("LTR Candidate Dotplots (Yes or No): ")

#Creates Headings and adds headings to open text file
with open("OutputFile", "a") as f:

	print("Sequence\t", "Local_Score\t", "Sequence_Length\t", "LTR_Length\t", "Global_Score\t", "Percent_Similarity\t", "PBS\t","Sequence_Complexity\t","Pass_Fail_Condition",file=f)

	blosum62 = substitution_matrices.load("BLOSUM62")


	#List of nonexhaustive list of primer binding sites associated with LTR retrotransposons
	PBSs = [] 
	with open ("PrimerBindingSites.txt", "rt") as PBS_File:

		for PBS in PBS_File:

			PBSs.append(PBS.rstrip('\n, \t'))

	#reads sequences as strings and formats those strings
	for seq_record in SeqIO.parse( file , "fasta"):
		sequence_string = seq_record.seq

		sequence = sequence_string.upper()

		Sequence_Complexity = lcc.lcc_simp(sequence)

		#this sets up the dotplots produced
		#the window is the number of bases being compared at once to create a point on the dotplot and is adjustable depending on the accuracy of the data
		#the dictionaries are the sequences converted into python strings and then the strings are searched according to window size. A match will create a point on the dotplot
		
		if dt == 'Yes' :

			window = 10
							
			dict_one = {}
			dict_two = {}

			for(seq, section_dict) in [
				(str(seq_record.seq).upper(), dict_one),
				(str(seq_record.seq).upper(), dict_two),
			]:
				for i in range(len(seq) - window):
					section = seq[i : i + window]
					try:
						section_dict[section].append(i)
					except KeyError:
						section_dict[section] = [i]

			matches = set(dict_one).intersection(dict_two)
		

			x = []
			y = []
			for section in matches:
				for i in dict_one[section]:
					for j in dict_two[section]:
						x.append(i)
						y.append(j)

		#python's find function spits out a "-1" if there is no match or the index of the first occurence of the substring if there is a match
		#the way the code as written works is that python searches for a match in the PBSs list outputting values based on whether or not there is a match
		#these values are put into the list names "results". Then results is searched for unique values and those unique values are put into the list "results2"
		#if the number of values in results2 is greater than 1, that means that the sequence has a primer binding site present in the PBSs list
		#if the number of values in results2 is equal to 1, that means the find function only ever delivered a "-1" meaning no primer binding site
		results = []

		def get_unique_numbers(numbers):
			unique = []

			for number in numbers:
				if number in unique:
					continue
				else:
					unique.append(number)
			return unique

		for i in range(len(PBSs)):
			results.append(sequence.find(PBSs[i]))
			results2 = get_unique_numbers(results)


		#this is the first of the if_else loops, the minimum base limit is 300 bases. If it is equal or greater than 300 bases, it moves on to the next step
		if len(sequence) >= 300:

			#Creating variables for local alignment
			midpoint = round((len(seq_record.seq))/2) 					
			endpoint = len(seq_record.seq)
			a_sequence = (seq_record.seq[0:99])
			b_sequence = (seq_record.seq[midpoint:endpoint])

			#local alignment comparing the first 100 bases to the second half of the sequence.
			#match = 4, mismatch = -1, opening gap = -4, extending gap = -2, perfect local score = 400
			alignments = pairwise2.align.localms(a_sequence.upper(), b_sequence.upper(), 4, -1, -4, -2)
			for alignment in alignments:
				local_score = alignment[2] #local alignment score
				local_begins = alignment[3] #coordinates for where local alignment starts in second half of sequence
				local_ends = alignment[4] #coordinates for where local alignment ends in second half of sequence

			#creating variable for global alignment
			d_sequence = b_sequence[local_begins:len(b_sequence)]
			c_sequence = seq_record.seq[0:len(d_sequence)]

			#if local alignment is equal to or greater than 200, the sequence moves on the next step.
			if local_score >= 200:
				
				#if estimated LTR length (from where the local alignment starts to the end of the sequence) is equal or greater than 90 bases, moves onto next step.
				if len(c_sequence) >= 90:

					#global alignment comparing two parts of the sequence equal in length, that length being determined by where the local alignment began to the end of the sequence
					#match =4, mismatch = -1, opening gap = -3, extending gap = -3
					new_alignments = pairwise2.align.globalms(d_sequence.upper(), c_sequence.upper(), 4, -1, -3, -3)
					for new_alignment in new_alignments:
						global_score = new_alignment[2] #global alignment score

						percent_sim = (global_score/(len(d_sequence)*4))

					#the percent similarity is determined by dividing the global alignment score by what the perfect score would have been
					#if the percent similarity is .75 or greater, it moves on to the next step
					if percent_sim >= .75:

						#this step removes sequence made up of simple repeats from passing
						#the sequence complexity takes into account the sequence length and the number of occurences of the four nucleotides (ACTG)
						#the sequence complexity is greater than or equal to 1.5
						if Sequence_Complexity >= 1.5:

							#if results2 has more than 1 value, there is a PBS present in the sequence. 
							if len(results2) > 1:
								print(f'{seq_record.id}\t{local_score}\t{len(seq_record)}\t{len(c_sequence)}\t{global_score}\t{percent_sim}\tYes PBS\t{Sequence_Complexity}\tLTR_Candidate', file = f)

								with open("LTR_Candidates.fasta", "a") as g:

									print(f'>{round(percent_sim*100)}_{seq_record.id}',file=g)
									print(sequence,file=g)
								
								if dt == 'Yes' :

									pylab.cla()
									pylab.gray()
									pylab.scatter(x, y, s=.20)
									pylab.xlim(0, len(seq_record.seq) - window)
									pylab.ylim(0, len(seq_record.seq) - window)
									pylab.xlabel("%s (length %i bp)" % (seq_record.id, len(seq_record.seq)))
									pylab.ylabel("%s (length %i bp)" % (seq_record.id, len(seq_record.seq)))
									pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
									savefile = f'{round(percent_sim*100)}_{seq_record.id.replace(".","_").replace(":","_").replace("-","_")}.png'
									pylab.savefig( savefile )
							
							#No Primer Binding Site
							else:
								print(f'{seq_record.id}\t{local_score}\t{len(seq_record)}\t{len(c_sequence)}\t{global_score}\t{percent_sim}\tNo PBS\t{Sequence_Complexity}\tLTR_Candidate', file = f)


								with open("LTR_Candidates.fasta", "a") as g:
									print(f'>{round(percent_sim*100)}_{seq_record.id}',file=g)
									print(sequence, file=g)
							
								if dt == 'Yes' :

									pylab.cla()
									pylab.gray()
									pylab.scatter(x, y, s=.20)
									pylab.xlim(0, len(seq_record.seq) - window)
									pylab.ylim(0, len(seq_record.seq) - window)
									pylab.xlabel("%s (length %i bp)" % (seq_record.id, len(seq_record.seq)))
									pylab.ylabel("%s (length %i bp)" % (seq_record.id, len(seq_record.seq)))
									pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
									savefile = f'{round(percent_sim*100)}_{seq_record.id.replace(".","_").replace(":","_").replace("-","_")}.png' 
									pylab.savefig( savefile )

						elif Sequence_Complexity < 1.5:
							#Simple Sequence w/ PBS
							if len(results2) > 1:
								print(f'{seq_record.id}\t{local_score}\t{len(seq_record)}\t{len(c_sequence)}\t{global_score}\t{percent_sim}\tYes PBS\t{Sequence_Complexity}\tSimple_Sequence', file = f)
							#Simple Sequence w/o PBS
							else:
								print(f'{seq_record.id}\t{local_score}\t{len(seq_record)}\t{len(c_sequence)}\t{global_score}\t{percent_sim}\tNo PBS\t{Sequence_Complexity}\tSimple_Sequence', file = f)




					elif percent_sim < .75:
						#Global Failure w/ PBS
						if len(results2) > 1:
							print(f'{seq_record.id}\t{local_score}\t{len(seq_record)}\t{len(c_sequence)}\t{global_score}\t{percent_sim}\tYes PBS\t{Sequence_Complexity}\tGlobal_Failure', file = f)
						#Global Failure w/o PBS
						else:
							print(f'{seq_record.id}\t{local_score}\t{len(seq_record)}\t{len(c_sequence)}\t{global_score}\t{percent_sim}\tNo PBS\t{Sequence_Complexity}\tGlobal_Failure', file = f)

				elif len(c_sequence) < 90:
					#LTR too small w/ PBS
					if len(results2) > 1:
						print(f'{seq_record.id}\t{local_score}\t{len(seq_record)}\t{len(c_sequence)}\t{percent_sim}\tN/A\tYes PBS\t{Sequence_Complexity}\tLTR_base_min', file = f)
					#LTR too small w/o PBS
					else:
						print(f'{seq_record.id}\t{local_score}\t{len(seq_record)}\t{len(c_sequence)}\t{percent_sim}\tN/A\tNo PBS\t{Sequence_Complexity}\tLTR_base_min', file = f)

			elif local_score < 200:
				#local score failure w/ PBS
				if len(results2) > 1:
					print(f'{seq_record.id}\t{local_score}\t{len(seq_record)}\tN/A\tN/A\tN/A\tYes PBS\t{Sequence_Complexity}\tLocal_Failure', file = f)
				#local score failure w/o PBS
				else:
					print(f'{seq_record.id}\t{local_score}\t{len(seq_record)}\tN/A\tN/A\tN/A\tNo PBS\t{Sequence_Complexity}\tLocal_Failure', file = f)


		elif len(sequence) < 300:
			#sequence too small w/ PBS
			if len(results2) > 1:
				print(f'{seq_record.id}\tN/A\t{len(seq_record)}\tN/A\tN/A\tN/A\tYes PBS\t{Sequence_Complexity}\tToo_Short', file = f)
			#sequence too small w/o PBS
			else:
				print(f'{seq_record.id}\tN/A\t{len(seq_record)}\tN/A\tN/A\tN/A\tNo PBS\t{Sequence_Complexity}\tToo_Short', file = f)
