from line_marker import *
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
context = "talk"
sns.set_context(context)

##############
# PARAMETERS #
##############

# plot all
plot_all = False


##########
# SCALES #
##########

# Figure 1
plot_scales_fig1 = False
if plot_scales_fig1 or plot_all:
    salts = ['LiCl', 'NaCl', 'KCl', 'CsCl']
    plt.rcParams['figure.figsize'] = 9.5, 7
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

# Figure 2
plot_scales_fig2 = False
if plot_scales_fig2 or plot_all:
    plt.rcParams['figure.figsize'] = 8, 6.5
    plt.figure('Scales_fig2')
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    df_exp = pd.read_csv(filepath_or_buffer='data/scales_fig2_data.csv',
                         names=('pH', 'zeta'))
    sns.scatterplot(x='pH', y='zeta',
                    data=df_exp)
    df_model = pd.read_csv(filepath_or_buffer='data_figures/scales_fig2_model.csv')
    g = sns.lineplot(x='pH', y='zeta',
                 data=df_model)
    g.set(xlabel='pH',
          ylabel='$\zeta$ [mV]')
    plt.text(6, -60, 'conc = 1e-3')


#########
# OSMAN #
#########

plot_osman = False
if plot_osman or plot_all:
    plt.rcParams['figure.figsize'] = 8, 6.5
    plt.figure('Osman')
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    label_list = ['Na-Li', 'K-Li', 'Cs-Li']
    i = 0
    for fig_index in [1, 2, 4]:
        label = label_list[i]
        color = next(colors)
        df_exp = pd.read_csv(filepath_or_buffer='data_figures/osman_fig'+str(fig_index)+'_exp.csv')
        sns.scatterplot(x='frac_A', y='frac_ads',
                        data=df_exp,
                        marker=next(markers),
                        color=color,
                        label=label)
        df_model = pd.read_csv(filepath_or_buffer='data_figures/osman_fig'+str(fig_index)+'_model.csv')
        plt.plot(df_model.get('frac_A'), df_model.get('frac_ads'),
                 color=color,
                 linestyle=next(lines),
                 label=" ")
        i += 1
    plt.legend(frameon=False,
               ncol=2)
    overhang = 0.04
    plt.xlim((0 - overhang, 1 + overhang))
    plt.ylim((0 - overhang, 1 + overhang))
    plt.xlabel('[M$^+$] / [Cl$^-$]')
    plt.ylabel('$\Gamma_M$ / CEC')


#################
# ISRAELACHVILI #
#################

plot_israelachvili = True
if plot_israelachvili or plot_all:
    plt.figure('israelachvili')
    c_list = ['1e-4', '1e-3', '1e-2', '1e-1']
    for c in c_list:
        df_exp = pd.read_csv(filepath_or_buffer='data/israelachvili_fig3_data_'+str(c)+'.csv',
                         names=('D', 'F/R'))
        sns.scatterplot(x='D', y='F/R',
                        data=df_exp,
                        marker=next(markers))
        df_model = pd.read_csv(filepath_or_buffer='data_figures/israelachvili_'+str(c)+'_model.csv')
        sns.lineplot(x='D', y='F/R',
                     data=df_model)
    plt.yscale('log')


plt.show()
