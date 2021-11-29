import csv
import numpy as np


def csvImport(filename = "test.csv"):
	csvArray = np.empty((0,3))

	with open(filename, mode='r') as csv_file:
		csv_reader = csv.DictReader(csv_file)
		line_count = 0
		for row in csv_reader:
		    if line_count == 0:
		        line_count += 1
		    csvArray = np.append(csvArray,np.array([[row["Points:0"],row["Points:1"],row["Points:2"]]]),axis = 0)
		    line_count += 1
		csvArray = csvArray[np.where(csvArray[:,2] == '0')].astype(np.float)
		return(csvArray)
	return("File not found!")
