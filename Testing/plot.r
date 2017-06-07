df_seq <- read.csv("results/res_seq.csv")
df_ff <- read.csv("results/res_ff.csv")
df_pt <- read.csv("results/res_pt.csv")

df_ff <- df_ff[order(df_ff$nw),]
df_pt <- df_pt[order(df_pt$nw),]

time_seq <- df_seq$time
time_ff_1 <- df_ff$time[1]
time_pt_1 <- df_pt$time[1]


scalab_ff <- time_ff_1/df_ff$time
speedup_ff <- time_seq/df_ff$time
eff_ff <- speedup_ff/df_ff$nw

scalab_pt <- time_pt_1/df_pt$time
speedup_pt <- time_seq/df_pt$time
eff_pt <- speedup_pt/df_pt$nw

plot(df_ff$nw, scalab_ff)
plot(df_ff$nw, speedup_ff)
plot(df_pt$nw, eff_ff)

plot(df_pt$nw, scalab_pt)
plot(df_pt$nw, speedup_pt)
plot(df_pt$nw, eff_pt)
