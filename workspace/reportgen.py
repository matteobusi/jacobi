import matplotlib.pyplot as plt
import matplotlib.ticker as plticker


import numpy as np
import seaborn as sns
import pandas as pd

sns.set()
sns.set(style="whitegrid")
sns.set(font_scale=2) 
plt.rc('text', usetex=True)
font = {'family' : 'Times-Roman',
        'size'   : 22}

plt.rc('font', **font)
color = sns.color_palette("Set2", 10)

for p in ["host", "mic"]:
    for m in ["ff", "th"]:
        compl_seq = pd.DataFrame()
        compl_m = pd.DataFrame()

        print("Plotting", p, m)
        for i in [5000, 10000, 15000, 30000]:
            df_seq = pd.read_csv("results/res_seq_" + p + "_" + str(i) + ".csv", delimiter=' *, *', engine='python')
            df_m = pd.read_csv("results/res_"+ m + "_" + p + "_" + str(i) +".csv", delimiter=' *, *', engine='python')
       
            seq_time = df_seq.loc[df_m['nw'] == 1]['latency'][0]
            m_time_0 = df_m.loc[df_m['nw'] == 1]['latency'][0]
                
            df_seq['N'] = df_seq.apply(lambda row: i, axis=1)
        
            df_m['N'] = df_m.apply(lambda row: i, axis=1)
            df_m['s'] = df_m.apply(lambda row: seq_time/row['latency'], axis=1)
            df_m['scalab'] = df_m.apply(lambda row: m_time_0/row['latency'], axis=1)
            df_m['eff'] = df_m.apply(lambda row: row['s']/row['nw'], axis=1)

            compl_seq = compl_seq.append(df_seq, ignore_index=True)
            compl_m = compl_m.append(df_m, ignore_index=True)

        # Add ideal values for efficiency, scalability and speedup to the dataframe
        for nw in compl_m['nw'].unique():
            nRow = pd.DataFrame([["ideal", nw, nw, nw, 1.0]],columns=['N', 'nw','s','scalab', 'eff'])
            compl_m = compl_m.append(nRow, ignore_index = True)

        maxW = compl_m['nw'].max()

        for (em, lbl) in [('s', 'Speedup'), ('scalab', 'Scalability'), ('eff', 'Efficiency')]:
            g = sns.FacetGrid(compl_m, hue='N', aspect=1.6, size=11)
            g.map(plt.scatter, 'nw', em)
            g.map(plt.plot, 'nw', em)
            g.add_legend()
            g.set_xlabels("# Workers")
            g.set_ylabels(lbl)
            plt.savefig("graphs/graph_" + m + "_" + p + "_" + em + ".png")

        # Compute the table and put into a file.
