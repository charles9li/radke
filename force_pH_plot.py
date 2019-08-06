import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set style
sns.set_context('talk')
plt.rcParams['figure.figsize'] = 6.25, 4.5
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})

# Parameters
c = 1e-3
cation = 'K'
pKa = 5.3

# Lists
pH_list = [4, 5, 6, 8, 10]

# Product figures
for pH in pH_list:
    filepath = 'data/force_' + cation + 'Cl_c' + str(c) + '_pH' + str(pH) + '_model.csv'
    df = pd.read_csv(filepath_or_buffer=filepath)
    sns.lineplot(x='D', y='F/R', data=df, label=str(pH))
plt.legend(title='pH', frameon=False)
plt.yscale('log')
plt.tight_layout()
plt.xlim((-1, 50))
plt.ylim((1e1, 1e4))
plt.xlabel('D [nm]')
plt.ylabel('F/R [$\mu$N m$^{-1}$]')
plt.text(17, 2e3, str(c) + r' $M$ ' + cation + 'Cl')
plt.show()
