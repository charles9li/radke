from line_marker import line_cycle, marker_cycle, color_cycle
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set style
sns.set_context('talk')
plt.rcParams['figure.figsize'] = 6.25, 4.5
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})

# Parameters
pH = 5.8
del_pH = 0.3
pKa = 5.3

# Lists
c_list = ['1e-4', '1e-3', '1e-2', '1e-1']

# Plot data
colors = color_cycle()
markers = marker_cycle()
lines = line_cycle()
for c in c_list:
    color = next(colors)

    filepath = 'data/israelachvili_fig3_data_'+str(c)+'.csv'
    df_exp = pd.read_csv(filepath_or_buffer=filepath,
                         names=('D', 'F/R'))
    sns.scatterplot(x='D',
                    y='F/R',
                    data=df_exp,
                    marker=next(markers),
                    color=color,
                    label=c)

    filepath = 'data/israelachvili_' + c + '_model.csv'
    df_model = pd.read_csv(filepath_or_buffer=filepath)
    plt.plot(df_model.get('D'),
             df_model.get('F/R_pH_mid'),
             linestyle=next(lines),
             label='c')
    plt.fill_between(df_model.get('D'),
                     df_model.get('F/R_pH_low'),
                     df_model.get('F/R_pH_high'),
                     alpha=0.3,
                     color=color)

plt.legend(frameon=False,
           ncol=2,
           title='conc [M]')
plt.yscale('log')
plt.xlabel('D [nm]')
plt.ylabel('F/R [$\mu$N/m]')
plt.text(60, 1e3, 'model')
plt.text(60, 2e3, 'exp')
plt.text(60, 1e2, 'pH = 5.8$\pm$0.3')
plt.xlim((-1, 125))
plt.ylim((1e1, 1e4))
plt.tight_layout()
plt.show()
