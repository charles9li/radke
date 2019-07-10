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
xlim_list = [(-1, 40), (-1, 40), (-1, 40), (-1, 40)]
ylim_list = [(1e2, 1e5), (1e2, 1e4), (1e2, 1e5), (1e2, 2e4)]

# Produce figures
for i in range(len(fig_index_list)):
    fig_index = str(fig_index_list[i])
    cation = cation_list[i]
    c_list = c_list_list[i]
    xlim = xlim_list[i]
    ylim = ylim_list[i]
    plt.figure(cation + 'Cl')
    for j in range(len(c_list)):
        c = c_list[j]
        color = 'C' + str(j)
        filepath = "data/pashley_fig" + fig_index + "_data_" + c + ".csv"
        df_exp = pd.read_csv(filepath_or_buffer=filepath,
                             names=('D', 'F/R'))
        sns.scatterplot(x='D',
                        y='F/R',
                        data=df_exp,
                        label=c,
                        color=color)

        filepath = "data/pashley_fig" + fig_index + "_model_" + c + ".csv"
        df_model = pd.read_csv(filepath_or_buffer=filepath)
        sns.lineplot(x='D',
                     y='F/R_pH_mid',
                     data=df_model,
                     color=color)
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
               title='Conc [M]')
    plt.text(15, 7e3, cation + 'Cl')
    plt.tight_layout()
plt.show()
