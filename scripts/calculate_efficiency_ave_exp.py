import numpy as np
TEMPERATURE = 300.
KB = 8.3144621E-3/4.184  # kcal/mol/K
BETA = 1. / TEMPERATURE / KB

def proportional_efficiency(assignments,vals,EXPONENTIAL_bool,n):
  N = len(vals)
  CORRECTION =float(N-n)/float(N-1)
  if EXPONENTIAL_bool:
    Exp_vals = []
    for v in vals:
      Exp_vals.append(np.exp(-BETA * v))
    Exp_vals = np.array(Exp_vals)
    vals = Exp_vals
  else:
    vals = vals

  if N!= len(assignments):
    raise Exception('assignments array length is not the same length as value array')
  unique_assignments = list(set(assignments))
  H = len(unique_assignments)
  vals_h = [[vals[ind] for ind in range(len(assignments)) \
           if assignments[ind]==unique_assignments[h] ] \
          for h in range(H)]
  N_h = [len(nz) for nz in vals_h if len(nz)>0]
  var_h = [np.var(nz) for nz in vals_h if len(nz)>0]
  prop_v = np.sum(np.array(N_h)/float(N)*var_h)/n
  #prop_v_f = np.sum(np.array(N_h)/float(N)*CORRECTION*var_h)/n         #**finite****

  random_var_normal = np.var(vals)/n
  #random_var_finite = np.var(vals)*CORRECTION/n
  prop_se = prop_v/random_var_normal
  #prop_se_f = prop_v_f/random_var_finite
  return prop_se #,prop_se_f


def optimal_efficiency(assignments,vals,EXPONENTIAL_bool,n):
  N = len(vals)
  CORRECTION =float(N-n)/float(N-1)
  if EXPONENTIAL_bool:
    Exp_vals = []
    for v in vals:
      Exp_vals.append(np.exp(-BETA * v))
    Exp_vals = np.array(Exp_vals)
    vals = Exp_vals
  else:
    vals = vals

  if N!= len(assignments):
    raise Exception('assignments array length is not the same length as value array')
  unique_assignments = list(set(assignments))
  H = len(unique_assignments)
  vals_h = [[vals[ind] for ind in range(len(assignments)) \
           if assignments[ind]==unique_assignments[h] ] \
          for h in range(H)]
  N_h = [len(nz) for nz in vals_h if len(nz)>0]
  var_h = [np.var(nz) for nz in vals_h if len(nz)>0]

  opti_v = np.sum(np.array(N_h)/float(N)*np.sqrt(var_h))**2/n
  #opti_v_f = np.sum((np.array(N_h)/float(N))*np.sqrt(CORRECTION)*np.sqrt(var_h))**2/n #**finite****
  
  random_var_normal = np.var(vals)/n
  #random_var_finite = np.var(vals)*CORRECTION/n

  opti_se = opti_v/random_var_normal
  #opti_se_f = opti_v_f/random_var_finite
  return opti_se#,opti_se_f


def senatorial_efficeincy_simple(assignments,vals,EXPONENTIAL_bool,n):
  N = len(vals)
  CORRECTION =float(N-n)/float(N-1)
  if EXPONENTIAL_bool:
    Exp_vals = []
    for v in vals:
      Exp_vals.append(np.exp(-BETA * v))
    Exp_vals = np.array(Exp_vals)
    vals = Exp_vals
  else:
    vals = vals

  if N!= len(assignments):
    raise Exception('assignments array length is not the same length as value array')
  unique_assignments = list(set(assignments))
  H = len(unique_assignments)
  vals_h = [[vals[ind] for ind in range(len(assignments)) \
             if assignments[ind]==unique_assignments[h] ] \
            for h in range(H)]
  N_h = [len(nz) for nz in vals_h if len(nz)>0]
  var_h = [np.var(nz) for nz in vals_h if len(nz)>0]
  sen_v = H*np.sum((np.array(N_h)/float(N))**2 *var_h)/n
  #sen_v_f = H*np.sum((np.array(N_h)/float(N))**2 * CORRECTION *var_h)/n#**finite****
  random_var_normal = np.var(vals)/n
  #random_var_finite = np.var(vals)*CORRECTION/n
  sen_se = sen_v/random_var_normal
  #sen_se_f = sen_v_f/random_var_finite
  return sen_se#,sen_se_f


def archetypal_senatorial_efficiency(assignments,vals,medoids,EXPONENTIAL_bool,n):
  N = len(vals)
  CORRECTION =float(N-n)/float(N-1)

  if EXPONENTIAL_bool:
    Exp_vals = []
    for v in vals:
      Exp_vals.append(np.exp(-BETA * v))
    Exp_vals = np.array(Exp_vals)
    vals = Exp_vals
  else:
    vals = vals

  if N!= len(assignments):
    raise Exception('assignments array length is not the same length as value array')
  unique_assignments = list(set(assignments))
  H = len(unique_assignments)
  vals_h = [[vals[ind] for ind in range(len(assignments)) \
             if assignments[ind]==unique_assignments[h] ] \
            for h in range(H)]
  N_h = [len(nz) for nz in vals_h if len(nz)>0]
  mean_h_ave = [np.mean(nz) for nz in vals_h if len(nz)>0]
  nsamples = 0
  reprs_ave = []
  
  random_var_normal = np.var(vals)/n
  random_var_finite = np.var(vals)*CORRECTION/n


  for medoid in medoids:
    for ind in medoid:
      nsamples += 1
      try:
        reprs_ave.append(vals[ind])
        break
      except:
        raise Exception(vals.shape, ind)
  reprs_ave = np.array(reprs_ave)
  if len(medoids) == len(mean_h_ave):
    clusters_ave_variance = []
    for i in range(len(medoids)):
      clusters_ave_variance.append( (reprs_ave[i] - mean_h_ave[i])**2)
    
    archetypal_sen_v = H*np.sum((np.array(N_h)/float(N))**2*clusters_ave_variance)/n
    #archetypal_sen_v_f = H*np.sum((np.array(N_h)/float(N))**2*CORRECTION*clusters_ave_variance)/n #**finite****

    archetypal_sen_se = archetypal_sen_v/random_var_normal
    #archetypal_sen_se_f = archetypal_sen_v_f/random_var_finite
  else:
    raise Exception('check the medoid')
  return archetypal_sen_se#,archetypal_sen_se_f

#----------------- function end -------------------------#
#-----starting assignment---------#
import pickle, gzip
import glob
import os.path

pkl_names = glob.glob("assignments_split_2/*A._hierarchical-*-maxclust.pkl.gz")+glob.glob("assignments_split_2/fingerprint*_hierarchical-*-maxclust.pkl.gz")

for pkl_name in pkl_names:
  import numpy as np
  import sys, os
  import fileLoader
  import argparse as args
  parser = args.ArgumentParser()
  parser.add_argument('--assignments', default=pkl_name)
  parser.add_argument('--observable_directory_inf', default='observables_with_inf')
  parser.add_argument('--observable_directory_0', default='observables_with_0.0') 
  parser.add_argument('--plot', action='store_true')
  args = parser.parse_args()

  F = gzip.open(args.assignments,'r')
  cluster_assignments = pickle.load(F)
  F.close()
  split = max([a['subset'] for a in cluster_assignments])+1

  for dirN in ['Stratification_stats_ave_split_%s'%split,'Stratification_plots_ave_split_%s'%split]:
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

  FNs = glob.glob(os.path.join(args.observable_directory_inf,'*.dat'))
  vals_inf = {}
  for FN in FNs:
    label = os.path.basename(FN)[:-3]
    vals_inf[label] = np.loadtxt(FN)

  stats = []
  EXPONENTIAL = 'False'      ########## average way or exponential average way
  n = 200                    ########## for finite
  for assignment in cluster_assignments:
    prop_SE=[]
    #prop_SE_finite=[]
    optimal_SE=[]
    #optimal_SE_finite = []
    senatorial_SE =[]
    #senatorial_SE_finite = []
    archetypal_sen_SE = []
    #archetypal_sen_SE_finite = []
    
    for label in labels:
      
      prop_se = proportional_efficiency(assignment['assignments'],vals_0[label][assignment['subset']::split],EXPONENTIAL,n)
      prop_SE.append(prop_se)
      #prop_SE_finite.append(prop_se_f)
      
      opti_se = optimal_efficiency(assignment['assignments'],vals_0[label][assignment['subset']::split],EXPONENTIAL,n)
      optimal_SE.append(opti_se)
      #optimal_SE_finite.append(opti_se_f)
      
      sen_se = senatorial_efficeincy_simple(assignment['assignments'],vals_0[label][assignment['subset']::split],EXPONENTIAL,n)
      senatorial_SE.append(sen_se)
      #senatorial_SE_finite.append(sen_se_f)
      
      archetypal_sen_se = archetypal_senatorial_efficiency(assignment['assignments'],vals_0[label][assignment['subset']::split],assignment['medoids'],EXPONENTIAL,n)
      archetypal_sen_SE.append(archetypal_sen_se)
      #archetypal_sen_SE_finite.append(archetypal_sen_se_f)
    
    prop_SE = np.array(prop_SE)
    #prop_SE_finite = np.array(prop_SE_finite)
    optimal_SE = np.array(optimal_SE)
    #optimal_SE_finite = np.array(optimal_SE_finite)
    senatorial_SE = np.array(senatorial_SE)
    #senatorial_SE_finite = np.array(senatorial_SE_finite)
    archetypal_sen_SE = np.array(archetypal_sen_SE)
    #archetypal_sen_SE_finite = np.array(archetypal_sen_SE_finite)
    
    stats.append({'threshold':assignment['threshold'],\
                 'nclusters':max(assignment['assignments']), \
                 'eta_prop':prop_SE,\
                 'eta_opt':optimal_SE,\
                 'eta_sen':senatorial_SE,\
                 'eta_archetypal_sen':archetypal_sen_SE})
    print assignment['threshold']
  # Store stats
  F = gzip.open('Stratification_stats_ave_split_%s/%s'%(split, os.path.basename(args.assignments)),'w')
  pickle.dump(stats,F)
  F.close()













