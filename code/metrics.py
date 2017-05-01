#!/usr/bin/python

import os, sys, re, statistics
from itertools import groupby
from collections import defaultdict, OrderedDict

def computeMetrics(dir):
	problems = defaultdict(list)
 	for filename in os.listdir(dir):
		if filename.endswith(".txt"):	
 			basename, extension = os.path.splitext(filename)
 			print "#Debug basename: {}".format(basename)
 			genopt, problem, seed = basename.split('-')
 			problems[problem].append(filename) # Agrupamos todos los nombres de ficheros del mismo problema
 	for problem_id, problem_files in problems.iteritems(): # problem_files es un array con todos los ficheros
 		problem_files.sort()
 		metric_median = dir + "/median/median_" + str(problem_id) + ".metric"
 		metric_boxplot = dir + "/boxplots/boxplots_" + str(problem_id) + ".metric"
 		print "#Debug files to open {}, {}".format(metric_boxplot, metric_median)
 		boxplot_file = open(metric_boxplot, "a")
 		median_file = open(metric_median, "a")
 		index = 1
 		for filename in problem_files:
 			file = open(filename, "r")
			median_results = defaultdict(list)
			index = 3
			while index < 49:
				line = file.readlines()[index].decode()
				file.seek(0, 0)
				splitted = re.split(r'\t+', line.rstrip('\t'))
				check_point = splitted[1]
				value = splitted[-1]
				print "#Debug index: {}, line: {}, check_point {}, value: {}".format(index, line, check_point, value)
				median_results[check_point].append(value)
				if index == 48:
					boxplot_file.write(value)
					print "#Debug line_to_write in boxplot: {}  ".format(line)
				index = index + 1
			 #median_results = OrderedDict(sorted(median_results.items()))
			for check_point, median in median_results.iteritems():
				line = str(check_point) + "  {}  ".format(statistics.median(median_results[check_point]))
				print "#Debug line_to_write: {}  ".format(line)
				median_file.write(line)

def createDirsAndFiles(dir):
		boxplots = dir + "/boxplots"
		median = dir + "/median"
		print "#Debug new directories \n {}, {}".format(boxplots, median)
		if not os.path.exists(boxplots):
			os.makedirs(boxplots)
		if not os.path.exists(median):
			os.makedirs(median)
		for i in xrange(0,18):
			boxplots_i = boxplots + "/boxplots_" + str(i) + ".metric"
			median_i = median + "/median_" + str(i) + ".metric"
			print "#Debug files to create: {}, {}".format(boxplots_i, median_i)
			open(boxplots_i, 'a').close()
			open(median_i, 'a').close()

def main(argv):
	if(len(argv) != 1):
		print "Error al iniciar el script. \n Usage: python metrics.py dir"
	else:
		print "#Debug: Dir : {}".format(argv[0])
		createDirsAndFiles(argv[0])
		computeMetrics(argv[0])

if __name__ == '__main__':
	main(sys.argv[1:])