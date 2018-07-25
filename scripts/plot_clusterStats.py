##############
# START HERE #
##############
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.ion()
figure_size=(3.2, 3.2*6/8)
dpi=300
fontsize=8
markersize=5
lw=1
text_pos = [ 0.55, 0.1 ]
font = { "fontname": "Arial"}
fmt_list = ['-o','-s','-^','-H','-v','-8','-d']
###########################################
import numpy as np
import glob
import os.path
import math
import pickle, gzip

whole_x_axis_new_mean_nclusters = {}
whole_y_axis_medoid_var = {}
whole_x_std_nclusters = {}
whole_y_std_variance = {}

for linkage in ['single','average','weighted','complete']:
  whole_x_axis_new_mean_nclusters[linkage] = {}
  whole_y_axis_medoid_var[linkage] = {}
  whole_x_std_nclusters[linkage]  = {}
  whole_y_std_variance[linkage]  = {}

  for eta_type in ['eta_prop','eta_opt','eta_sen','eta_archetypal_sen']:
    whole_x_axis_new_mean_nclusters[linkage][eta_type] = {}
    whole_y_axis_medoid_var[linkage][eta_type] = {}
    whole_x_std_nclusters[linkage][eta_type] = {}
    whole_y_std_variance[linkage][eta_type] = {}
    stats_directory = 'Stratification_stats_ave'+ '_split_2/'
  
    FNs = glob.glob(os.path.join(stats_directory,'*.pkl.gz'))
    FNs = sorted(FNs)
    all_stats = {}
    for FN in FNs:
      method = os.path.basename(FN)[:-7]
      F = gzip.open(FN,'r')
      all_stats[method] = pickle.load(F)
      F.close()
    methods = sorted(all_stats.keys())
      
    criterion ='maxclust'
    linkage_stats = []
    methods_order=[]
    for key in methods:
      if (key.find(linkage)>-1) and (key.find(criterion)>-1):
        linkage_stats.append(all_stats[key])
        methods_order.append(key)
    for i in range(len(linkage_stats)):
      whole_x_axis_new_mean_nclusters[linkage][eta_type][i] =[]
      whole_y_axis_medoid_var[linkage][eta_type][i] =[]
      whole_x_std_nclusters[linkage][eta_type][i] =[]
      whole_y_std_variance[linkage][eta_type][i] =[]
        
      thresholds = sorted(list(set([a['threshold'] for a in linkage_stats[i]])))
      mean_nclusters = [np.mean([np.mean(s['nclusters']) for s in linkage_stats[i] if s['threshold']==threshold]) for threshold in thresholds]
        
      std_nclusters = []
      for threshold in thresholds:
        tem_nclusters = []
        for s in linkage_stats[i]:
          if s['threshold']==threshold:
            tem_nclusters.append(s['nclusters'])
        std_nclusters.append(np.std(tem_nclusters))
      std_eta=[]
      mean_eta =[]
      for threshold in thresholds:
        tem_list = []
        tem_std_list = []
        for s in linkage_stats[i]:
          if s['threshold']==threshold:
            new_s_eta_type_total_times = []
            for times in range(100):
              new_s_eta_type_per_time=[]
              resamples_list = np.random.choice(len(s[eta_type]), len(s[eta_type]))
              for num_index in resamples_list :
                value = s[eta_type][num_index]
                if not math.isnan(value) and not math.isinf(value):
                  new_s_eta_type_per_time.append(np.sqrt(value))
              this_times_mean_value = np.mean(new_s_eta_type_per_time)
              new_s_eta_type_total_times.append(this_times_mean_value)
            bootstrap_std = np.std(new_s_eta_type_total_times)
            bootstrap_mean =  np.mean(new_s_eta_type_total_times)
            tem_list.append(bootstrap_mean)
            tem_std_list.append(bootstrap_std)

        mean_value= np.mean(tem_list)
        mean_std_value = np.mean(tem_std_list)
        std_eta.append(mean_std_value)
        mean_eta.append(mean_value)
  
      new_mean_nclusters=[]
      for m in mean_nclusters:
        new_mean_nclusters.append(np.sqrt(m))

      whole_x_axis_new_mean_nclusters[linkage][eta_type][i]  =new_mean_nclusters
      whole_y_axis_medoid_var[linkage][eta_type][i]  =mean_eta
      whole_x_std_nclusters[linkage][eta_type][i]  =std_nclusters
      whole_y_std_variance[linkage][eta_type][i]  =std_eta

#-----------------------------------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fignum = 0
#plt.ion()
figure_size=(6.66,4.8)
dpi=300
fontsize=8
markersize=5
lw=1
text_pos = [ 0.55, 0.1 ]
font = { "fontname": "Arial"}
fmt_list = ['-o','-s','-^','-H','-v','-d','*','D']
picked_list = [0,1,2,3,5,11,12,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33]

for linkage in ['single','average','weighted','complete']:
  for eta_type in ['eta_prop','eta_opt','eta_sen','eta_archetypal_sen']: # SE
    plt.figure(fignum,figsize=figure_size)
    colors = ['red','green','blue','purple','cyan','magenta','orange']

    shift = 0.0
      
    total_mean_eta = []
    
    for i in range(len(linkage_stats)):
      new_mean_nclusters_org = whole_x_axis_new_mean_nclusters[linkage][eta_type][i]
      std_nclusters_org = whole_x_std_nclusters[linkage][eta_type][i]
      mean_eta_org = whole_y_axis_medoid_var[linkage][eta_type][i]
      std_eta_org = whole_y_std_variance[linkage][eta_type][i]
        
      new_mean_nclusters = []
      mean_eta = []
      std_eta = []
      std_nclusters=[]
        
      for l in picked_list:
        new_mean_nclusters.append(new_mean_nclusters_org[l])
        mean_eta.append(mean_eta_org[l])
        std_eta.append(std_eta_org[l])
        std_nclusters.append(std_nclusters_org[l])
        
        total_mean_eta.append(mean_eta_org[l])
        
      new_mean_nclusters = np.array(new_mean_nclusters)
      mean_eta = np.array(mean_eta)
      std_eta = np.array(std_eta)
      std_nclusters = np.array(std_nclusters)
        
      shift+=0.1
      plt.errorbar(new_mean_nclusters+shift,mean_eta,xerr=std_nclusters,yerr = std_eta, color=colors[i], label =methods_order[i].split('.')[0],fmt=fmt_list[i],markersize=5,linewidth = 2, elinewidth=1,capthick = 1,capsize =2)
      Area = np.trapz(new_mean_nclusters,mean_eta)
      print Area
    xl  = np.array([ num for num in range(1,21)])
    yl = xl/xl
      
    max_value = max(total_mean_eta)
    #plt.axis([0, 20, 0, 3])
    #plt.ylim(0, round(max_value*1.25,1))
    #plt.xlim(1, 20)

    
    x_label = r'$\sqrt[2]{\mathcal{H}}$'
    if not eta_type in ['eta_archetypal_sen','eta_archetypal_sen_finite']:
      y_label = r'$\eta$'
    else:
      y_label = r'$\eta^*$'
    
    plt.xlabel( x_label, fontsize=10,**font)
    plt.ylabel( y_label, fontsize=10,**font)
    plt.tight_layout()
    plt.show()
    plt.savefig("ave_%s-%s.png"%(linkage,eta_type))
    fignum += 1
    plt.clf()

#-------------------------------------------------------------------------------------------------------------#



