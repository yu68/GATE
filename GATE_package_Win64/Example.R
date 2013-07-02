source("FMM-HMM-functions.R")

#input and preprocess the data
data=Inputdata("data.txt",6,3)

#run the core program to get the clustering and hidden state information, set ncluster=20 here.
results=FMM.HMM.program(data$observation,20,1500,30,0.5,1)

#group clusters with similar patterns (lambda), set ngroup=5 here.
results.data.cl.gr=grouping.clusters(data$location,results,5)
#a file named "result_cl-gr-hid.txt" will also be generated with cluster, group and hidden state info for each region.

#generate colorful BED for visualization in UCSC genome browser
color.BED(results.data.cl.gr,results$hidden)
#a file named "data_col_BED.bed" will be generated in the same folder
