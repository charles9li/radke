from line_marker import line_cycle, marker_cycle, color_cycle
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set style
sns.set_context('talk')
plt.rcParams['figure.figsize'] = 6.25, 4.5
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})

# Lists
fig_index_list = [1, 2, 3, 4]
cation_list = ['Li', 'Na', 'K', 'Cs']
c_list_list = [['6e-2', '1e-2', '1e-3', '1e-4'],
               ['1e-2', '1e-3', '1e-4', '4e-5'],
               ['1e-3', '3e-4', '4e-5'],
               ['1e-3', '1e-4', '4e-5']]
pH_list = [5.4, 5.7, 5.7, 5.8]
xlim_list = [(-1, 40), (-1, 40), (-1, 40), (-1, 40)]
ylim_list = [(1e2, 1e5), (1e2, 1e5), (1e2, 1e5), (1e2, 1e5)]

# Produce figures
for i in range(len(fig_index_list)):
    fig_index = str(fig_index_list[i])
    cation = cation_list[i]
    c_list = c_list_list[i]
    pH = pH_list[i]
    xlim = xlim_list[i]
    ylim = ylim_list[i]
    lines = line_cycle()
    markers = marker_cycle()
    colors = color_cycle()
    plt.figure(cation + 'Cl')
    for j in range(len(c_list)):
        c = c_list[j]
        color = next(colors)
        filepath = "data/pashley_fig" + fig_index + "_data_" + c + ".csv"
        df_exp = pd.read_csv(filepath_or_buffer=filepath,
                             names=('D', 'F/R'))
        plt.scatter(df_exp.get('D'),
                    df_exp.get('F/R'),
                    label=c,
                    color=color,
                    marker=next(markers))

        filepath = "data/pashley_fig" + fig_index + "_model_" + c + ".csv"
        df_model = pd.read_csv(filepath_or_buffer=filepath)
        plt.plot(df_model.get('D'),
                 df_model.get('F/R_pH_mid'),
                 color=color,
                 linestyle=next(lines),
                 label=c)
        plt.fill_between(df_model.get('D'),
                         df_model.get('F/R_pH_low'),
                         df_model.get('F/R_pH_high'),
                         alpha=0.3,
                         color=color)

    plt.yscale('log')
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [$\mu$N m$^{-1}$]')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.legend(frameon=False,
               title='Conc [M]',
               ncol=2)
    plt.text(12, 3e3, cation + 'Cl\n' + r'pH$=%.1f \pm0.3$' % pH)
    plt.text(0, 2e2, 'model')
    plt.text(0, 5e2, 'exp')
    plt.tight_layout()
plt.show()
