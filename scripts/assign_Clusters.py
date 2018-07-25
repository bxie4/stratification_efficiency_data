import sys, os
import argparse as args
import numpy as np
from scipy.cluster import hierarchy

parser = args.ArgumentParser()
parser.add_argument('--distances',default='distances/fingerprint_jaccard.dat',
  help = 'A netCDF file with a distance matrix')
parser.add_argument('--split_array', default = 2, type=int,
  help = 'An option to split the arrays to get statistics for estimates.\
    If the value is greater than 1, the output will have a set of \
    assignments and medoids for each distance cutoff.')

parser.add_argument('--criterion', default='fraction',
  choices = ['distance','maxclust','fraction'],
                    
  help = 'The criterion to use in forming flat clusters.')
parser.add_argument('--t_range', nargs=3, default = [1.0, 6.0, 26],
  help = 'Range of thresholds used for hierarchical clustering.' + \
    'Format for input is *start_value end_value number_of_steps*.' + \
    '(Ignored for fraction)')
args = parser.parse_args()

distance_name = 'fingerprint_jaccard'
distance_matrix = np.loadtxt(args.distances)
#from scipy.spatial.distance import squareform
#distance_matrix = squareform(distance_matrix)# for CA.dat, siteCA.dat,siteHA.dat

if args.criterion=='fraction':
  args.criterion = 'maxclust'
  sz = distance_matrix.shape[0]/args.split_array
  t_range = np.array(np.hstack((np.linspace(1,14,14),\
    np.floor(sz*np.logspace(np.log(15./sz)/np.log(10),0,20)))),dtype=int)
else:
  start_t, end_t, number_of_steps = args.t_range
  t_range = np.linspace(float(start_t), float(end_t), int(number_of_steps))

cluster_assignments = []
for c in range(args.split_array):
  distance_matrix_c = distance_matrix[c::args.split_array,c::args.split_array]
  
  for linkage_method in ['complete','single','average','weighted']:
    cluster = hierarchy.linkage(squareform(distance_matrix_c), \
      method=linkage_method)

    nclusters_o = 0
    for t in t_range:
      fcluster = hierarchy.fcluster(cluster,t,criterion=args.criterion)
      if max(fcluster)!=nclusters_o:
        medoids = []
        for n in range(1,max(fcluster)+1):
          inds = [i for i in range(len(fcluster)) if fcluster[i]==n]
          distance_matrix_n = distance_matrix_c[inds][:,inds]
          medoids.append([inds[i] \
            for i in np.argsort(np.mean(distance_matrix_n,0))])
    
      else:
        print 'fcluster',fcluster,'t_range value', t,'nclusters_o',nclusters_o
      nclusters_o = max(fcluster)
      cluster_assignments.append(\
        {'method':'hierarchical-%s-%s'%(linkage_method,args.criterion),
         'threshold':t,
         'subset':c,
         'assignments':fcluster,
         'medoids':medoids})

if not os.path.isdir('assignments_split_%d'%args.split_array):
  os.mkdir('assignments_split_%d'%args.split_array)

import pickle, gzip
methods = set([a['method'] for a in cluster_assignments])
for method in methods:
  cluster_assignments_c = [a for a in cluster_assignments if a['method']==method]
  F = gzip.open('assignments_split_%d/%s_%s.pkl.gz'%(args.split_array,distance_name,method),'w')
  pickle.dump(cluster_assignments_c,F)
  F.close()
