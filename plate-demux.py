#!/usr/bin/env python

"""
Demultiplex S3-ATAC-Seq data using plates.
Input: paired-end FASTQ files
Input: config file of well positions, tn5 indices, and biological samples used
Output: forward and reverse reads at the biological sample level in FASTQ format.
Code was inspired by James Adler
This program uses a nested dictionary to keep track of reads per sample.
To prevent excessive RAM usage (~ 1GB RAM/1M reads), this program stores reads in a buffer
that gets written to file periodically (user-defined).
"""

import os
import argparse
import sys
import logging
import pandas as pd
import pysam
import gzip
from collections import defaultdict

def parse_args():

	parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter, description = """
	
plate-demux.py demultiplexes ATAC-Seq/S3-ATAC-Seq reads performed in 96-well plates. This script assumes that a Tn5 index (e.g. an 8bp DNA barcode) tags for a specific biological sample, and that index appears in a specific position in a 96-well plate (e.g. index ACTAAGTAA in A12). By relating the index to the coordinates of a plate, this sript will demultiplex a FASTQ file of mixed samples into separate files.\n

This script accepts paired FASTQ files processed with unidex () and a configuration file. The FASTQ file should contain a mix of samples, while the configuration is a table that specifies an index and its coordinates in a 96-well plate. The output is a folder that contains foward and reverse reads per individual sample.""")

	parser.add_argument(
		"-R1",
		help = "Forward read (R1) of a paired sequencing run.",
		required = True
	)

	parser.add_argument(
		"-R2",
		help = "Reverse read (R2) of a paired sequencing run.",
		required = True
	)

	parser.add_argument(
		"-o",
		"--outdir",
		help = "Output directory",
		required = True
	)

	parser.add_argument(
		"-c",
		"--config",
		type = str,
		help = "File containing tn5 indices, well positions, and biological samples.",
		required = True
	)

	parser.add_argument(
		"-b",
		"--buffer",
		type = int,
		help = "Maximum number of reads to load into memory. The program generally requires approximately 1GB of memory per 1M reads.",
		default = None
	)

	parser.add_argument(
		"-v",
		"--verbose",
		help = "Adjust the verbosity of the program",
		default = False,
		action='store_true'
	)

	args = parser.parse_args()
	return args


def read_indices(index_file):
	logging.info("Extracting indices from {}".format(index_file))
	df = pd.read_table(index_file)
	return df


def parse_indices(df):

	""" Associate each biological sample with a Tn5 index. """

	tn5_dict = defaultdict(list)

	for i in range(0, df.shape[0]):
		sample_name = df.iloc[i,].loc["sample"]
		tn5_index = df.iloc[i,].loc["Tn5_index"]

		if tn5_index not in tn5_dict[sample_name]:
			tn5_dict[sample_name].append(tn5_index)

	return tn5_dict

def index_length(df):

	""" Get the length of Tn5 indices """
	barcode_lengths = df["Tn5_index"].apply(len)
	if len( barcode_lengths.unique() ) == 1:
		return int( df["Tn5_index"].apply(len).unique() )
	else:
		logging.info("Length of Tn5 indices not uniform")
		os.exit()

def write_reads(sample_info, write_counter):

	"""
	Write reads to designated output file defined by the nested dict `sample_info`.
	sample_info:
		sample_1:
			R1: outdir/sample_1_R1.fastq.gz
			R2: outdir/sample_1_R2.fastq.gz
			R1_reads: [read_1, read_2, ..., read_n]
			R2_reads: [read_1, read_2, ..., read_n]
		sample_2:
			R1: outdir/sample_2_R1.fastq.gz
			R2: outdir/sample_2_R2.fastq.gz
			R1_reads: [read_1, read_2, ..., read_n]
			R2_reads: [read_1, read_2, ..., read_n]
	"""

	all_samples = list(sample_info.keys())

	for s in all_samples:

		# create a new file, or re-open a file to write into.

		# write-mode = bulk or initiate buffered file
		# write_counter counts from 0
		if write_counter == 0:
			r1_out = gzip.open(sample_info[s]["R1"], "wb+")
			r2_out = gzip.open(sample_info[s]["R2"], "wb+")

		# write-mode = append compressed results to existing file
		elif write_counter >= 1:
			r1_out = gzip.open(sample_info[s]["R1"], "ab")
			r2_out = gzip.open(sample_info[s]["R2"], "ab")

		# define and write the reads per sample.
		all_fw = sample_info[s]["R1_reads"]
		all_rv = sample_info[s]["R2_reads"]
		if len(all_fw) == len(all_rv):
			for fw,rv in zip(all_fw, all_rv):
				r1_out.write(  str(str(fw) + '\n').encode()  )
				r2_out.write(  str(str(rv) + '\n').encode()  )
		else:
			logging.info("ERROR: Number of read pairs do not match in sample {}".format(s))
			os.exit(1)

		r1_out.close()
		r2_out.close()

def leftover_reads(sample_info):
	# check if there are reads left in the buffer.
	# if any sample contains > 0 read, return True
	all_samples = list(sample_info.keys())
	leftover_results = []
	for s in all_samples:
		r1_len = len(sample_info[s]["R1_reads"])
		r2_len = len(sample_info[s]["R2_reads"])
		res = any([ r1_len > 0, r2_len > 0 ])
		leftover_results.append(res)
	if any(leftover_results):
		return(True)
	else:
		return(False)

def wipe_buffer(sample_info):
	# remove all read contents in sample_info
	all_samples = list(sample_info.keys())
	for s in all_samples:
		sample_info[s]["R1_reads"] = list()
		sample_info[s]["R2_reads"] = list()
	return sample_info

def output_message(sample_info, verbose):

	all_samples = sample_info.keys()
	
	if verbose:
		for s in all_samples:
			logging.info(  "Sample {a} R1 output file: {b}".format(a = s, b = sample_info[s]["R1"])  )
			logging.info(  "Sample {a} R2 output file: {b}".format(a = s, b = sample_info[s]["R2"])  )
			logging.info(  "Sample {a} R1 number of reads: {b}".format(a = s, b = str(len(sample_info[s]["R1_reads"]))  ))
			logging.info(  "Sample {a} R2 number of reads: {b}".format(a = s, b = str(len(sample_info[s]["R2_reads"]))  ))
			logging.info(  "Sample {a} total number of reads: {b}".format(a = s, b = str(sample_info[s]["number_of_reads"])  ))

	logging.info("Completed barcode splitting.")
	logging.info("Total reads processed: {}".format(read_counter))
	logging.info("Total reads with no matching tn5 index: {}".format(dropped_reads))

def parse_reads(R1, R2, outdir, buffer_size, tn5_dict, tn5_len, df, verbose):

	all_samples = list( tn5_dict.keys() )

	# initialize dict of sample information. Example below
	sample_info = defaultdict(dict)
	for s in all_samples:
		r1 = os.path.join(outdir, s) + "_R1.fastq.gz"
		r2 = os.path.join(outdir, s) + "_R2.fastq.gz"
		# r1 = os.path.join(outdir, s) + "_R1.fastq"
		# r2 = os.path.join(outdir, s) + "_R2.fastq"
		sample_info[s]["R1"] = r1
		sample_info[s]["R2"] = r2
		sample_info[s]["R1_reads"] = list()
		sample_info[s]["R2_reads"] = list()
		sample_info[s]["number_of_reads"] = 0

	"""
	sample_info
		sample_1:
			R1: outdir/sample_1_R1.fastq.gz
			R2: outdir/sample_1_R2.fastq.gz
			R1_reads: [forward_1, forward_2, forward_3]
			R2_reads: [reverse_1, reverse_2, reverse_3]
		sample_2:
			R1: outdir/sample_2_R1.fastq.gz
			R2: outdir/sample_2_R2.fastq.gz
			R1_reads: [forward_1, forward_2, forward_3]
			R2_reads: [reverse_1, reverse_2, reverse_3]
	"""

	# demultiplex reads

	R1 = pysam.FastxFile(R1)
	R2 = pysam.FastxFile(R2)

	# define constants
	read_counter = 0
	buffer_counter = 0
	write_counter = 0
	dropped_reads = 0

	for fw, rv in zip(R1, R2):

		read_counter += 1
		buffer_counter += 1

		# test dataset
		# if read_counter == 1_000_000 + 1:
		# 	break

		# logging info
		if verbose:
			if read_counter % 100_000 == 0:
				logging.info("Processing {} reads".format(read_counter))

		# identify forward and reverse indices
		fw_index = fw.name.split(":")[0][-tn5_len:]
		rv_index = rv.name.split(":")[0][-tn5_len:]

		# test fw and rv indices match
		if fw_index == rv_index:
			sample = "".join(df.loc[df["Tn5_index"] == fw_index, "sample"])
			sample_info[sample]["R1_reads"].append(fw)
			sample_info[sample]["R2_reads"].append(rv)
			sample_info[sample]["number_of_reads"] += 1
		else:
			dropped_reads += 1

		# write and reset reads buffer
		if buffer_size != None:
			if buffer_counter == buffer_size:

				# write the reads to each sample
				write_reads(sample_info, write_counter)

				# wipe the reads in the nested dictionary sample_info
				sample_info = wipe_buffer(sample_info)

				logging.info("Writing to output this many times: {}".format(write_counter))

				write_counter += 1
				buffer_counter = 0

	# EOF

	# if in buffered mode and if there are reads left in the buffer,
	# append the reads to the sample outfile one last time.
	if buffer_size != None:
		if leftover_reads(sample_info):
			logging.info("Writing results to file for the last time.")
			write_reads(sample_info, write_counter = write_counter)

	# else if in "bulk" mode, write the whole buffer to each sample.
	else:
		logging.info("Writing results to file.")
		write_reads(sample_info, write_counter = 0) # 0 means writing in 'bulk' mode

	# print summary information
	logging.info("Exporting read split summary")
	output_message(sample_info, verbose)
	return None

def main():

	args = parse_args()

	df = read_indices(index_file = args.config)
	tn5_dict = parse_indices(df)

	if not os.path.exists(args.outdir):
		logging.info("Exporting reads in " + args.outdir)
		os.mkdir(args.outdir)

	tn5_len = index_length(df)

	# what if buffer > number of reads?

	parse_reads(
		R1 = args.R1,
		R2 = args.R2,
		outdir = args.outdir,
		buffer_size = args.buffer,
		tn5_dict = tn5_dict,
		tn5_len = tn5_len,
		df = df,
		verbose = args.verbose
	)

	logging.info("PCR plate split complete.")


if __name__ == "__main__":
	logging.basicConfig(
		format = '%(asctime)s: %(levelname)s: %(message)s',
		level = logging.INFO
	)
	main()


