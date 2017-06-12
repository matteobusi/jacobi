import matplotlib.pyplot as plt
import matplotlib.ticker as plticker


import numpy as np
import seaborn as sns
import pandas as pd

sns.set()
sns.set(style="whitegrid")
sns.set(font_scale=2.5) 
plt.rc('text', usetex=True)
font = {'family' : 'Times-Roman'}

plt.rc('font', **font)
color = sns.color_palette("Set2", 10)

for p in ["host", "mic"]:
    for m in ["ff", "th"]:
        compl_seq = pd.DataFrame()
        compl_m = pd.DataFrame()
        df_min = pd.DataFrame()

        print("Plotting", p, m)
        print("========== proc: {} = impl: {} ==========".format(p, m))
        
        for i in [5000, 10000, 15000, 30000]:
            df_seq = pd.read_csv("results/res_seq_" + p + "_" + str(i) + ".csv", delimiter=' *, *', engine='python')
            df_m = pd.read_csv("results/res_"+ m + "_" + p + "_" + str(i) +".csv", delimiter=' *, *', engine='python')
       
            seq_time = df_seq.loc[df_seq['nw'] == 1]['latency'][0]
            m_time_0 = df_m.loc[df_m['nw'] == 1]['latency'][0]
                
            df_seq['N'] = df_seq.apply(lambda row: i, axis=1)
        
            df_m['N'] = df_m.apply(lambda row: i, axis=1)
            df_m['s'] = df_m.apply(lambda row: seq_time/row['latency'], axis=1)
            df_m['scalab'] = df_m.apply(lambda row: m_time_0/row['latency'], axis=1)
            df_m['eff'] = df_m.apply(lambda row: row['s']/row['nw'], axis=1)

            compl_seq = compl_seq.append(df_seq, ignore_index=True)
            compl_m = compl_m.append(df_m, ignore_index=True)

            df_1 = df_m[ df_m['nw'] == 1 ][['N', 'latency']]
            df_best = df_m[ df_m['latency'] == df_m['latency'].min() ][['N', 'nw', 'latency']]
            df_best['latency_1'] = df_best.apply(lambda row: df_1[df_1['N'] == row['N']]['latency'], axis=1)
            df_best['latency_seq'] = df_best.apply(lambda row: seq_time, axis=1)
            df_min = df_min.append(df_best, ignore_index=True)

        # Add ideal values for efficiency, scalability and speedup to the dataframe
        for nw in compl_m['nw'].unique():
            nRow = pd.DataFrame([["ideal", nw, nw, nw, 1.0]],columns=['N', 'nw','s','scalab', 'eff'])
            compl_m = compl_m.append(nRow, ignore_index = True)

        maxW = compl_m['nw'].max()

        for (em, lbl) in [('s', 'Speedup'), ('scalab', 'Scalability'), ('eff', 'Efficiency')]:
            g = sns.FacetGrid(compl_m, hue='N', aspect=1.3, size=8)
            g.map(plt.scatter, 'nw', em)
            g.map(plt.plot, 'nw', em)
            g.add_legend()
            g.set_xlabels("Number of workers")
            g.set_ylabels(lbl)
            plt.savefig("graphs/graph_" + m + "_" + p + "_" + em + ".png")

        # Produce tables. 
        # For each N and target processor compute 
        #   - FF : best nw, latency (best), latency (1)
        #   - TH : best nw, latency (best), latency(1), ratio
        #   - Seq : latency
        print(df_min)
        print("==============================================")