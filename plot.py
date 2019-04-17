from line_marker import *
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
sns.set_context("paper")

##########
# SCALES #
##########

# Figure 1
plot_scales_fig1 = False
if plot_scales_fig1:
    salts = ['LiCl', 'NaCl', 'KCl', 'CsCl']
    plt.figure('Scales_fig1')
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    for salt in salts:
        color = next(colors)
        df_exp = pd.read_csv(filepath_or_buffer='data/scales_fig1_data_'+salt+'.csv',
                             names=('c', 'zeta'))
        g = sns.scatterplot(x='c', y='zeta',
                        data=df_exp,
                        marker=next(markers),
                        label=salt,
                        color=color)
        g.set(xscale='log',
              xlabel='c [mol dm$^{-3}$]',
              ylabel='$\zeta$ [mV]')
        df_model = pd.read_csv(filepath_or_buffer='data_figures/scales_fig1_model_'+salt+'.csv')
        plt.plot(df_model.get('c'), df_model.get('zeta_pH5.8'),
                 linestyle=next(lines),
                 label=" ",
                 color=color)
        plt.fill_between(df_model.get('c'), df_model.get('zeta_pH5.5'), df_model.get('zeta_pH6.1'),
                         label='_nolegend_',
                         alpha=0.2,
                         color=color)
    plt.legend(frameon=False,
               ncol=2)


#########
# OSMAN #
#########

plot_osman = True
if plot_osman:
    plt.figure('Osman')
    for fig_index in [1, 2, 4]:
        df = pd.read_csv(filepath_or_buffer='data_figures/osman_fig'+str(fig_index)+'_exp.csv')
        sns.scatterplot(x='frac_A', y='frac_ads',
                        data=df)

    plt.xlim((0, 1))
    plt.ylim((0, 1))


#################
# ISRAELACHVILI #
#################

plt.show()
