from line_marker import line_cycle, marker_cycle, color_cycle
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set style
sns.set_context('talk')
plt.rcParams['figure.figsize'] = 6.25, 4.5
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})

# Lists
cation_list = ['Li', 'Na', 'K', 'Cs']

# Plot data
plt.figure('grahame')
i = 0
for cation in cation_list:
    df = pd.read_csv(filepath_or_buffer='data/grahame_' + cation + 'Cl.csv')
    plt.plot(df.get('D'), df.get('sigma_d')/0.32,
             color='C'+str(i),
             label=cation+'Cl')
    i += 1

i = 0
for cation in cation_list:
    df = pd.read_csv(filepath_or_buffer='data/grahame_' + cation + 'Cl.csv')
    plt.plot(df.get('D'), df.get('sigma_d_grahame')/0.32,
             color='C'+str(i),
             linestyle='--',
             label='(Grahame eqn) (model)')
    i += 1

plt.legend(frameon=False,
           ncol=2)
plt.xlabel('D [nm]')
plt.ylabel('$q_d$ / |$q_0$|')
plt.text(10, 0.04, 'c = 1e-3 M\npH = 5.8')
plt.tight_layout()
plt.show()
