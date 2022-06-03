# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 07:19:58 2019

@author: amyk0
"""
#%%

import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np

def read_db_nondups(filename):
    df = pd.read_excel(filename, sheet_name = 'NonDups', skiprows = [0])
    return df

# def get_lengths(df):
#     lengths = df['Length'].values
#     return lengths
    
# def get_net_charge(df):
#     net_charge = df['Net Charge'].values
#     return net_charge

#%%

directory = os.getcwd()
subfolder = directory + '\\input'

db = []
p_title = 'POSH Lipobead'
p_title2 = 'Sample Comparison'
d = True # True returns probability distn, area = 1

with os.scandir(subfolder) as it:
    for entry in it:
        print('Processing file', entry.name)
        df = read_db_nondups(entry)
        db.append( (df, entry.name[:-4]) )
            
ls = []
cs = []
hs = []
leg = []
for d in db:
    df = d[0]
    l = df['Length'].values
    l = np.array(l)
    ls.append(l)
    
    c = df['Net Charge'].values
    c = np.array(c)
    cs.append(c)
    
    h = df['Hydropathy'].values
    h = np.array(h)
    hs.append(h)
    
    leg.append(d[1])

#%%
# https://matplotlib.org/3.1.1/gallery/statistics/histogram_multihist.html    
fig, axes = plt.subplots(nrows=3, ncols=1)
ax0, ax1, ax2 = axes.flatten()

ax0.hist(ls, range=(0,66), bins=33, histtype='bar', density=d)
ax0.tick_params(which='both', labelsize='xx-small')
ax0.set_title(p_title + '\n' + p_title2 + '\n\nLength Distribution', fontsize='xx-small')
# ax0.legend(labels=leg, bbox_to_anchor=(1.01,1), loc='upper left', fontsize='xx-small')
ax0.legend(labels=leg, fontsize='xx-small')

ax1.hist(cs, range=(-8,16), bins=25, histtype='bar', density=d)
ax1.tick_params(which='both', labelsize='xx-small')
ax1.set_title('Net Charge Distribution', fontsize='xx-small')
# ax1.legend(labels=leg, bbox_to_anchor=(1.01,1), loc='upper left', fontsize='xx-small')            
ax1.legend(labels=leg, fontsize='xx-small')

ax2.hist(hs, histtype='bar', density=d)
ax2.set_title('Hydropathy Distribution', fontsize='xx-small')
ax2.tick_params(which='both', labelsize='xx-small')
# ax2.legend(labels=leg, bbox_to_anchor=(1.01,1), loc='upper left', fontsize='xx-small')
ax2.legend(labels=leg, fontsize='xx-small')

fig.tight_layout()
fig.savefig(p_title + ' ' + p_title2 + '.png', dpi=300)
plt.close()

  
            
