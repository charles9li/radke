import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set style
sns.set_context('talk')
plt.rcParams['figure.figsize'] = 6.25, 4.5
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})

# Produce figure
plt.figure("Pashley fig 13")
filepath = "data/pashley_fig13_data.csv"
df_exp = pd.read_csv(filepath_or_buffer=filepath,
                     names=('D', 'F/R'))
sns.scatterplot(x='D', y='F/R', data=df_exp)

filepath = "data/pashley_fig13_model.csv"
df_model = pd.read_csv(filepath_or_buffer=filepath)
sns.lineplot(x='D', y='F/R', data=df_model)

plt.xlabel('D [nm]')
plt.ylabel('F/R [$\mu$N m$^{-1}$]')
plt.xlim((-1, 60))
plt.ylim((5e1, 1e5))
plt.tight_layout()
plt.yscale('log')
plt.show()
