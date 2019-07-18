from line_marker import *
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set style
sns.set_context('talk')
plt.rcParams['figure.figsize'] = 6.25, 4.5
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})

# Lists
pH_list = [3, 10]
c_list = [1, 5, 10]

# Plot data
for pH in pH_list:
    colors = color_cycle()
    markers = marker_cycle()
    lines = line_cycle()
    plt.figure('donaldson_pH' + str(pH))
    for c in c_list:
        color = next(colors)
        filepath = 'data/donaldson_pH' + str(pH) + '_' + str(c) + 'mM_data.csv'
        df_exp = pd.read_csv(filepath_or_buffer=filepath,
                             names=('D', 'F/R'))
        sns.scatterplot(x='D',
                        y='F/R',
                        data=df_exp,
                        color=color,
                        marker=next(markers),
                        label=str(c))

        filepath = 'data/donaldson_pH' + str(pH) + '_' + str(c) + 'mM_model.csv'
        df_model = pd.read_csv(filepath_or_buffer=filepath)
        plt.plot(df_model.get('D'),
                 df_model.get('F/R') * 1e-3,
                 color=color,
                 linestyle=next(lines),
                 label=str(c))

    if pH == 3:
        plt.xlim((0, 25))
    else:
        plt.xlim((0, 40))
    plt.yscale('log')
    plt.tight_layout()

plt.show()
