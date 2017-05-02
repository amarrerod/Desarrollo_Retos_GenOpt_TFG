#!/usr/bin/python

import os, sys, re, statistics
from itertools import groupby, chain
from collections import defaultdict, OrderedDict

def computeMetrics(source_dir, destiny_dir):
	problems = defaultdict(list)
 	for filename in os.listdir(source_dir):
		if filename.endswith(".txt"):	
 			basename, extension = os.path.splitext(filename)
 			print "#Debug basename: {}".format(basename)
 			genopt, problem, seed = basename.split('-')
 			problems[problem].append(filename) # Agrupamos todos los nombres de ficheros del mismo problema
 	for problem_id, problem_files in problems.iteritems(): # problem_files es un array con todos los ficheros
 		problem_files.sort()
 		metric_median = destiny_dir + "/median/median_" + str(problem_id) + ".metric"
 		metric_boxplot = destiny_dir + "/boxplots/boxplots_" + str(problem_id) + ".metric"
 		print "#Debug files to open {}, {}".format(metric_boxplot, metric_median)
 		boxplot_file = open(metric_boxplot, "a")
 		median_file = open(metric_median, "a")
 		index = 1
 		median_results = defaultdict(list)
 		# Recorremos todos los ficheros del problema 1-100
 		for filename in problem_files:
 			file = open(filename, "r")
			index = 3
			# Leemos todas las lineas de cada fichero
			while index < 49:
				line = file.readlines()[index]
				file.seek(0, 0)
				splitted = re.split(r'\t+', line.rstrip('\t'))
				check_point = splitted[1]
				value = float(splitted[-1])
				print "#Debug index: {}, check_point {}, value: {}".format(index, check_point, value)
				median_results[check_point].append(value)
				print "#Debug nuevo valor! checkpoint: {} value: {}".format(check_point, median_results[check_point])
				if index == 48: ##Ultima linea del fichero --> Valor final
					boxplot_file.write(str(value))
					print "#Debug line_to_write in boxplot: {}  ".format(line)
				index = index + 1
				print "\n"
		# Rellenamos los ficheros de metricas para cada problema
		for check_point in sorted(median_results.keys(), key = lambda x: int(x)):
			print "#Debug check_point {}".format(check_point)
			print "#Debug valores: {}".format(median_results[check_point]) 
			line = " " + check_point + "  {}  \n".format(statistics.median(median_results[check_point]))
			print "#Debug line_to_write: {}  ".format(line)
			median_file.write(line)
		median_results.clear() # Vaciamos el diccionario

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
	if(len(argv) != 2):
		print "Error al iniciar el script. \n Usage: python metrics.py source-dir destiny-dir"
	else:
		print "#Debug: source: {}, destiny: {}".format(argv[0], argv[1])
		createDirsAndFiles(argv[1])
		computeMetrics(argv[0], argv[1])

if __name__ == '__main__':
	main(sys.argv[1:])