import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import numpy as np
import seaborn as sns
import pandas as pd

sns.set()
sns.set(style="whitegrid")
sns.set_context("paper")

for p in ["host", "mic"]:
    for m in ["ff", "pt"]:
        compl_seq = pd.DataFrame()
        compl_m = pd.DataFrame()

        print("Plotting ", p, m)
        for i in [5000, 10000, 15000, 30000]:
            df_seq = pd.read_csv("results/res_seq_" + p + "_" + str(i) + ".csv", delimiter=' *, *', engine='python')
            df_m = pd.read_csv("results/res_"+ m + "_" + p + "_" + str(i) +".csv", delimiter=' *, *', engine='python')
       
            seq_time = df_seq.loc[df_m['nw'] == 1]['time'][0]
            m_time_0 = df_m.loc[df_m['nw'] == 1]['time'][0]
                
            df_seq['N'] = df_seq.apply(lambda row: i, axis=1)
        
            df_m['N'] = df_m.apply(lambda row: i, axis=1)
            df_m['s'] = df_m.apply(lambda row: seq_time/row['time'], axis=1)
            df_m['scalab'] = df_m.apply(lambda row: m_time_0/row['time'], axis=1)
            df_m['eff'] = df_m.apply(lambda row: row['s']/row['nw'], axis=1)

            compl_seq = compl_seq.append(df_seq, ignore_index=True)
            compl_m = compl_m.append(df_m, ignore_index=True)

        # Add ideal values for efficiency, scalability and speedup to the dataframe
        for nw in compl_m['nw'].unique():
            nRow = pd.DataFrame([["ideal", nw, nw, nw, 1.0]],columns=['N', 'nw','s','scalab', 'eff'])
            compl_m = compl_m.append(nRow)    

        maxW = compl_m['nw'].max()
    
        fig, axs = plt.subplots(ncols=3, figsize=(11.7*4, 8.27))
        sns.factorplot(x='nw', y='s', data=compl_m, hue='N',  ax=axs[0])
        sns.factorplot(x='nw', y='scalab', data=compl_m, hue='N', ax=axs[1])
        sns.factorplot(x='nw', y='eff', data=compl_m, hue='N',  ax=axs[2])
    
        axs[0].set(ylabel="Speedup", xlabel="# Workers")
        axs[1].set(ylabel="Scalability", xlabel="# Workers")
        axs[2].set(ylabel="Efficiency", xlabel="# Workers")

        
        fig.savefig("graphs/graph_" + m + "_" + p + ".png")
