import numpy as np
def proportional_efficiency(assignments,vals,boostrap):
  N = len(vals)#400
  if N!= len(assignments):
    raise Exception('assignments array length is not the same length as value array')
  unique_assignments = list(set(assignments))
  H = len(unique_assignments)
  vals_h = [[vals[ind] for ind in range(len(assignments)) \
   if assignments[ind]==unique_assignments[h] ] \
  for h in range(H)]
  N_h = [len(nz) for nz in vals_h ]
  nsamples_h = [int(round(H*(float(nh)/np.sum(N_h)),0)) for nh in N_h]
  prop_se_per_boostrap=[]
  for i in range(boostrap):
    numerator_variance_per_boostrap = []
    denominator_sample_variance_tem =[]
    for n_index in range(len(nsamples_h)):
      if nsamples_h[n_index] >0:
        new_nz = np.random.choice(vals_h[n_index],size=nsamples_h[n_index])
        min_h_boostrap= min(vals_h[n_index])
        diff_square =np.square(np.min(new_nz)-min_h_boostrap)
        if new_nz[0]!=min_h_boostrap:
          print new_nz, min_h_boostrap
        random_square=[np.square(n-np.min(vals)) for n in new_nz]
        numerator_variance_per_boostrap.append(diff_square)
        for r in random_square:
          denominator_sample_variance_tem.append(r)
      else:
        numerator_variance_per_boostrap.append(0.0)

    random_variance= np.sum(denominator_sample_variance_tem)/len(denominator_sample_variance_tem)
    var_prop = np.sum(np.array(N_h)/float(N)*numerator_variance_per_boostrap)
    print var_prop
    val_random=(float(N-H+1)/float(N))*random_variance
    prop_se_per_boostrap.append(var_prop/val_random)

  prop_se=np.mean(prop_se_per_boostrap)
  prop_se_err=np.std(prop_se_per_boostrap)
  return prop_se,prop_se_err

def optimal_efficiency(assignments,vals,boostrap):
  N = len(vals)
  if N!= len(assignments):
    raise Exception('assignments array length is not the same length as value array')
  unique_assignments = list(set(assignments))
  H = len(unique_assignments)
  vals_h = [[vals[ind] for ind in range(len(assignments)) \
    if assignments[ind]==unique_assignments[h] ] \
    for h in range(H)]
  N_h = [len(nz) for nz in vals_h]
  nsamples_h = [int(round(H*(float(nh)/np.sum(N_h)),0)) for nh in N_h]
  opt_se_per_boostrap=[]
  for i in range(boostrap):
    numerator_variance_per_boostrap = []
    denominator_sample_variance_tem =[]
    for n_index in range(len(nsamples_h)):
      if nsamples_h[n_index] >0:
        new_nz = np.random.choice(vals_h[n_index],size=nsamples_h[n_index])
        min_h_boostrap= min(vals_h[n_index])
        diff_square =np.square(np.min(new_nz)-min_h_boostrap)
        random_square=[np.square(n-np.min(vals)) for n in new_nz]
        numerator_variance_per_boostrap.append(diff_square)
        for r in random_square:
          denominator_sample_variance_tem.append(r)
      else:
        numerator_variance_per_boostrap.append(0.0)

    random_variance= np.sum(denominator_sample_variance_tem)/len(denominator_sample_variance_tem)
    opti_v = np.sum(np.array(N_h)/float(N)*np.sqrt(random_variance))**2
    opti_random=(float(N-H+1)/float(N))**2*random_variance
    opt_se_per_boostrap.append(opti_v/opti_random)
      
  opti_se=np.mean(opt_se_per_boostrap)
  opti_se_err=np.std(opt_se_per_boostrap)
  return opti_se,opti_se_err


def senatorial_efficeincy_simple(assignments,vals,boostrap):
  N = len(vals)
  if N!= len(assignments):
    raise Exception('assignments array length is not the same length as value array')
  unique_assignments = list(set(assignments))
  H = len(unique_assignments)

  vals_h = [[vals[ind] for ind in range(len(assignments)) \
    if assignments[ind]==unique_assignments[h] ] \
    for h in range(H)]
  N_h = [len(nz) for nz in vals_h]
  
  sen_se_per_boostrap =[]
  for i in range(boostrap):
    numerator_variance_per_boostrap = []
    denominator_sample_variance_tem = []
    
    for n_index in range(len(vals_h)):
      if len(vals_h[n_index])>0:
        new_nz = np.random.choice(vals_h[n_index], size=1)
        min_h_boostrap= min(vals_h[n_index])
        diff_square =np.square(new_nz-min_h_boostrap)
        numerator_variance_per_boostrap.append(diff_square)
        random_square=[np.square(x-np.min(vals)) for x in vals_h[n_index]]
        for r in random_square:
          denominator_sample_variance_tem.append(r)
      else:
        numerator_variance_per_boostrap.append(0.0)
    random_variance= np.sum(denominator_sample_variance_tem)/len(denominator_sample_variance_tem)
    sen_v =H*np.sum((np.array(N_h)/float(N))**2 *numerator_variance_per_boostrap)
    sen_se_per_boostrap.append(sen_v/random_variance)
  sen_se=np.mean(sen_se_per_boostrap)
  sen_se_err=np.std(sen_se_per_boostrap)
  return sen_se,sen_se_err


def archetypal_senatorial_efficiency(assignments,vals,medoids):
  N = len(vals)
  if N!= len(assignments):
    raise Exception('assignments array length is not the same length as value array')
  unique_assignments = list(set(assignments))
  H = len(unique_assignments)
  vals_h = [[vals[ind] for ind in range(len(assignments)) \
    if assignments[ind]==unique_assignments[h] ] \
    for h in range(H)]
  N_h = [len(nz) for nz in vals_h if len(nz)>0]

  min_h_ave = [np.min(nz) for nz in vals_h if len(nz)>0]
  nsamples = 0
  reprs_ave = []
  total_val_random = np.mean([np.square( v - min(vals)) for v in vals])
  random_var_normal = total_val_random/n

  for medoid in medoids:
    for ind in medoid:
      nsamples += 1
      try:
        reprs_ave.append(vals[ind])
        break
      except:
        raise Exception(vals.shape, ind)
  reprs_ave = np.array(reprs_ave)
  if len(medoids) == len(min_h_ave):
    clusters_ave_variance = []
    for i in range(len(medoids)):
      clusters_ave_variance.append( (reprs_ave[i] - min_h_ave[i])**2)
    archetypal_sen_v = np.sum((np.array(N_h)/float(N))**2*clusters_ave_variance)/(float(n)/float(H))
    archetypal_sen_se = archetypal_sen_v/random_var_normal
    return archetypal_sen_se
  else:
    raise Exception('check the medoid')

#-----starting assignment---------#
import pickle, gzip
import glob
import os.path

import numpy as np
import sys, os
import fileLoader
import argparse as args
parser = args.ArgumentParser()
parser.add_argument('--d_matrix', default='PCA')
parser.add_argument('--linkage', default='average')
parser.add_argument('--split_num', default='2')
parser.add_argument('--observable_directory_0', default='observables_with_0.0')
parser.add_argument('--plot', action='store_true')
args = parser.parse_args()

pkl_name ='assignments_split_'+args.split_num+'/'+args.d_matrix+'_hierarchical-'+args.linkage+'-maxclust.pkl.gz'
F = gzip.open(pkl_name,'r')
cluster_assignments = pickle.load(F)
F.close()
split = max([a['subset'] for a in cluster_assignments])+1
  
for dirN in ['Stratification_stats_min_split_%s'%split]:
  if not os.path.isdir(dirN):
    os.mkdir(dirN)

# Load observables into vals and compute true_stats
FNs = glob.glob(os.path.join(args.observable_directory_0,'*.dat'))
vals_0 = {}
for FN in FNs:
  label = os.path.basename(FN)[:-3]
  vals_0[label] = np.loadtxt(FN)
labels = [os.path.basename(FN)[:-3] for FN in FNs]
del FNs

stats = []
boostrap = 100
for assignment in cluster_assignments:
  prop_SE=[]
  prop_SE_err=[]
  optimal_SE=[]
  optimal_SE_err = []
  for label in labels:
      
    prop_se,prop_se_err = proportional_efficiency(assignment['assignments'],vals_0[label][assignment['subset']::split],boostrap)
    prop_SE.append(prop_se)
    prop_SE_err.append(prop_se_err)
      
    opti_se,opti_se_err = optimal_efficiency(assignment['assignments'],vals_0[label][assignment['subset']::split],boostrap)
    optimal_SE.append(opti_se)
    optimal_SE_err.append(opti_se_err)
  
    sen_se,sen_se_err = senatorial_efficeincy_simple(assignment['assignments'],vals_0[label][assignment['subset']::split],boostrap)
    senatorial_SE.append(sen_se)
    senatorial_SE_err.append(sen_se_err)
    
    archetypal_sen_se = archetypal_senatorial_efficiency(assignment['assignments'],vals_0[label][assignment['subset']::split],assignment['medoids'])
    archetypal_sen_SE.append(archetypal_sen_se)


  prop_SE = np.array(prop_SE)
  prop_SE_err = np.array(prop_SE_err)
  optimal_SE = np.array(optimal_SE)
  optimal_SE_err = np.array(optimal_SE_err)
  senatorial_SE = np.array(senatorial_SE)
  senatorial_SE_err = np.array(senatorial_SE_err)
  archetypal_sen_SE = np.array(archetypal_sen_SE)

  stats.append({'threshold':assignment['threshold'],\
                'nclusters':max(assignment['assignments']), \
                'eta_prop':prop_SE,'eta_prop_std':prop_SE_err,\
                'eta_opt':optimal_SE, 'eta_opt_std':optimal_SE_err,\
                'eta_sen':senatorial_SE, 'eta_sen_finite':senatorial_SE_err,\
                'eta_archetypal_sen':archetypal_sen_SE})
  print assignment['threshold']
# Store stats
F = gzip.open('Stratification_stats_min_split_%s/%s'%(split, os.path.basename(pkl_name)),'w')
pickle.dump(stats,F)
F.close()
