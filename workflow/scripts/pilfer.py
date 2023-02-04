#!/usr/bin/python
import numpy as np
from operator import itemgetter
import sys, csv, os, argparse
import pandas as pd


parser = argparse.ArgumentParser(usage="pilfer.py [-h] [--factor <value>] -i input.BED", description="Predicts piRNA clusters from an input BED file with the reads and the raw read counts in the 5th column.\nBy default output is printed in stardard output")
parser.add_argument("-i", metavar="<input filename>", help="Input BED file. - if input is stdin", required=True)
parser.add_argument("-f", default=3.0, metavar="<value>", type=float, help="The factor by which the read should be away from standard deviation to be called a peak")
args = parser.parse_args()

#Function for extracting columns
def column(matrix, i):
	return [row[i] for row in matrix]

def read_data_pandas(infile_pt):
	map = pd.read_csv(infile_pt, header=0, sep='\t', names = ['chr', 'start', 'end', 'pir', 'rpm', 'strand', 'count'])

	rpm_fact = map['count'].sum()/1000000
	map['rpm'] = map['count']*rpm_fact
	#print(map.head())
	mean = map['rpm'].mean()
	sd = map['rpm'].std()
	#print(rpm_fact,mean,sd)
	map = map.values.tolist()
	return mean, sd, map

def create_dict(map):
	for row in map:
		if row[0] in chrom_dict:
			chrom_dict[row[0]].append(row)
		else:
			chrom_dict[row[0]] = [row]
	return chrom_dict
#PILFER algorith

def pilfer2(chrom_dict, mean, sd, columnidx):
	pilfer_result = []
	for key in chrom_dict:
		j = 0
		while j < len(chrom_dict[key]):
			row = chrom_dict[key][j]
			if (row[columnidx] - mean)/sd >= sd_factor:
				start_bp = row[2] - 100000
				if start_bp < 0:
					start_bp = 0
				end_bp = row[2]
				score = 0
				new_score = 0
				start_read_index = -1
				end_read_index = 0
				
				cur_read_index = chrom_dict[key].index(row)

				#Calculate the initial score and start index
				for read in chrom_dict[key]:
					if read[1] >= start_bp and read[1] <= start_bp + 100000:
						score += read[columnidx]
						if start_read_index == -1:
							start_read_index = chrom_dict[key].index(read)
							start_bp = read[1]
						end_read_index = chrom_dict[key].index(read)

				new_score = score
				max_start = start_read_index
				max_end = end_read_index
				#calculating optimum 100KB window
				for i in range(start_read_index+1,cur_read_index+1):
					new_score = score - chrom_dict[key][i-1][columnidx]
					while  end_read_index+1 < len(chrom_dict[key]) and chrom_dict[key][end_read_index +1][1] <= chrom_dict[key][i][1] + 100000 :
						new_score += chrom_dict[key][end_read_index+1][columnidx]
						end_read_index += 1
					
					if new_score > score:
						score = new_score
						max_start = i
						max_end = end_read_index

				#print key + "\t" + str(chrom_dict[key][max_start][1]) + "\t" + str(chrom_dict[key][max_end][2]) + "\t" + str(score)
				#print(key + ":" + str(chrom_dict[key][max_start][1]) + "-" + str(chrom_dict[key][max_end][2]) + "\t" + str(score))
				pilfer_result.append([key,chrom_dict[key][max_start][1],chrom_dict[key][max_end][2],score])
				j = max_end
			j += 1
	return pilfer_result


#Variables
map = []
chrom_dict = {}
columnidx = 4

infile = args.i
sd_factor = args.f

if infile == "-":
	infile_pt = sys.stdin
else:
	if not os.path.isfile(infile):
		sys.stderr.write("Not a valid input file\n")
		sys.exit()
	else:
		infile_pt = open(infile,"r")


#Reading the BED records
mean, sd, map = read_data_pandas(infile_pt)

#Creating the dictionary
chrom_dict = create_dict(map)

#Sorting the dictionary
for key in chrom_dict:
	chrom_dict[key] = sorted(chrom_dict[key],key=itemgetter(1))
	

#Rul PILFER2
pilfer_result = pilfer2(chrom_dict, mean, sd, columnidx)

for row in pilfer_result:
	print(str(row[0]) + ":" + str(row[1]) + "-" + str(row[2]) + "\t" + str(row[3]))
