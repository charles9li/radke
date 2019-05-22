from line_marker import *
from radii import *
from solver import *
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
context = "talk"
sns.set_context(context)
plt.rcParams['figure.figsize'] = 6.25, 4.5
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})

##############
# PARAMETERS #
##############

# K_ads
K_Li = 10**0.7
K_Na = K_Li*2.5
K_K = 10**2.8
K_Cs = 10**3
K_ads_list = [K_Li, K_Na, K_K, K_Cs]

# pH
pH = 5.8

# pKa
pKa = 5.3

# Permittivity
eps1 = 6*epsilon_0
eps2 = 30*epsilon_0
eps_bulk = 79*epsilon_0

# C1
C1_Li = eps1/R_Li
C1_Na = eps1/R_Na
C1_K = eps1/R_K
C1_Cs = eps1/R_Cs
C1_list = [C1_Li, C1_Na, C1_K, C1_Cs]

# C2
C2_Li = eps2/(2*R_Li_hyd)
C2_Na = eps2/(2*R_Na_hyd)
C2_K = eps2/(2*R_K_hyd)
C2_Cs = eps2/(2*R_Cs_hyd)
C2_list = [C2_Li, C2_Na, C2_K, C2_Cs]

# plot all
plot_all = False


##########
# SCALES #
##########

# Figure 1
plot_scales_fig1 = False
if plot_scales_fig1 or plot_all:
    salts = ['LiCl', 'NaCl', 'KCl', 'CsCl']
    # textlocations = [(1e-3, -120), (1e-3, -120), (1e-3, -120), (1e-3, -120)]
    # xlims = [(1e-4, 2e-2), (9e-5, 2e-2), (2e-5, 2e-2), (5e-5, 1e-2)]
    i = 0
    for salt in salts:
        plt.figure('Scales_fig1_'+salt)
        plt.rcParams['xtick.top'] = True
        plt.rcParams['ytick.right'] = True
        color = next(colors)
        df_exp = pd.read_csv(filepath_or_buffer='data/scales_fig1_data_'+salt+'.csv',
                             names=('c', 'zeta'))
        g = sns.scatterplot(x='c', y='zeta',
                        data=df_exp,
                        marker=next(markers),
                        label='_nolegend_',
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
        plt.text(1e-3, -120, 'pH = 5.8$\pm$0.3')
        parameter_label = salt+'\n'+\
                          '$C_1$=%i $\mu$F/m$^2$\n' +\
                          '$C_2$=%i $\mu$F/m$^2$\n' +\
                          'p$K_M$=%.2f\n' +\
                          'p$K_a$=%.2f'
        plt.text(3e-5, -90, parameter_label % (C1_list[i]*100, C2_list[i]*100, np.log10(K_ads_list[i]), pKa))
        plt.tight_layout()
        plt.xlim((2e-5, 2e-2))
        plt.ylim((-150, -25))
        i += 1

# Figure 2
plot_scales_fig2 = False
if plot_scales_fig2 or plot_all:
    # plt.rcParams['figure.figsize'] = 8, 6.5
    plt.figure('Scales_fig2')
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    df_exp = pd.read_csv(filepath_or_buffer='data/scales_fig2_data.csv',
                         names=('pH', 'zeta'))
    sns.scatterplot(x='pH', y='zeta',
                    data=df_exp,
                    marker='v',
                    color='C1')
    df_model = pd.read_csv(filepath_or_buffer='data_figures/scales_fig2_model.csv')
    g = sns.lineplot(x='pH', y='zeta',
                 data=df_model,
                 color='C1')
    g.set(xlabel='pH',
          ylabel='$\zeta$ [mV]')
    label = 'KCl\n' + \
            'conc=10$^{-3}$ M\n' + \
            'p$K_M$=%.2f\n' + \
            'p$K_a$=%.2f\n' + \
            '$C_1$=%i $\mu$F/m$^2$\n' + \
            '$C_2$=%i $\mu$F/m$^2$\n'
    plt.text(7, -75, label % (np.log10(K_K), pKa, C1_K*100, C2_K*100))
    plt.tight_layout()


#########
# OSMAN #
#########

plot_osman = False
if plot_osman or plot_all:
    # plt.rcParams['figure.figsize'] = 8, 6.5
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
                        color='C'+str(i),
                        label=label)
        df_model = pd.read_csv(filepath_or_buffer='data_figures/osman_fig'+str(fig_index)+'_model.csv')
        plt.plot(df_model.get('frac_A'), df_model.get('frac_ads'),
                 color='C'+str(i),
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
    plt.tight_layout()


#################
# ISRAELACHVILI #
#################

plot_israelachvili = False
if plot_israelachvili or plot_all:
    plt.figure('israelachvili')
    c_list = ['1e-4', '1e-3', '1e-2', '1e-1']
    i = 0
    for c in c_list:
        df_exp = pd.read_csv(filepath_or_buffer='data/israelachvili_fig3_data_'+str(c)+'.csv',
                         names=('D', 'F/R'))
        sns.scatterplot(x='D', y='F/R',
                        data=df_exp,
                        marker=next(markers),
                        color='C'+str(i),
                        label=c)
        df_model = pd.read_csv(filepath_or_buffer='data_figures/israelachvili_'+str(c)+'_model.csv')
        plt.plot(df_model.get('D'), df_model.get('F/R'),
                 linestyle=next(lines),
                 color='C'+str(i),
                 label='_nolegend_')
        i += 1
    plt.legend(frameon=False,
               ncol=2,
               title='conc [mol L$^{-1}$]')
    plt.yscale('log')
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [$\mu$N/m]')
    plt.tight_layout()


#############
# DONALDSON #
#############

plot_donaldson = False
if plot_donaldson or plot_all:
    for pH in [3, 10]:
        plt.figure('donaldson_pH' + str(pH))
        c_list = [1, 5, 10]
        i = 0
        for c in c_list:
            df_exp = pd.read_csv(filepath_or_buffer='data/donaldson_data_pH'+str(pH)+'_'+str(c)+'mM.csv',
                                 names=('D', 'F/R'))
            sns.scatterplot(x='D', y='F/R',
                            data=df_exp,
                            marker=next(markers),
                            color='C' + str(i),
                            label=c)
            df_model = pd.read_csv(filepath_or_buffer='data_figures/donaldson_pH'+str(pH)+'_'+str(c)+'mM_model.csv')
            plt.plot(df_model.get('D'), df_model.get('F/R'),
                     linestyle=next(lines),
                     color='C' + str(i),
                     label='_nolegend_')
            i += 1
        plt.legend(frameon=False,
                   title='conc [mM]')
        plt.xlabel('D [nm]')
        plt.ylabel('F/R [mN/m]')
        plt.tight_layout()
        plt.yscale('log')


######################
# POTENTIAL PROFILES #
######################

plot_LiCl_conc = False
if plot_LiCl_conc or plot_all:
    pH = 5.8
    sns.set_palette('viridis')
    plt.figure('LiCl potential profile vary concentration')
    c_list = np.array([1e-1, 1e-2, 1e-3, 1e-4])
    sigma_d_LiCl_conc = np.zeros(len(c_list))
    i = 0
    for c in c_list:
        c_list1 = np.array([c+10**-pH, c, 10**-pH])
        K_list = np.array([K_Li, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1_Li, C_2=C2_Li)
        sol.solver_PB()
        sigma_d_LiCl_conc[i] = sol.sigma_d
        x = np.concatenate((np.array([-(R_Li+2*R_Li_hyd), -2*R_Li_hyd]), sol.x))
        psi = np.concatenate((np.array([sol.psi_0, sol.psi_beta]), sol.psi))
        plt.plot(x*1e9, psi*1000, label='10$^{%i}$' % np.log10(c))
        i += 1
    plt.plot(np.array([0, 0]), np.array([-700, 0]), color='k', linestyle='--', linewidth=1)
    plt.plot(np.array([-2*R_Li_hyd, -2*R_Li_hyd])*1e9, np.array([-700, 0]), color='k', linestyle='--', linewidth=1)
    plt.plot(np.array([-2*R_Li_hyd-R_Li, -2*R_Li_hyd-R_Li])*1e9, np.array([-700, 0]), color='k', linestyle='--', linewidth=1)
    plt.xlabel('x [nm]')
    plt.ylabel('$\psi$ [mV]')
    plt.text(0, -500, '10$^{-1}$ M')
    plt.text(0, -600, '10$^{-4}$ M')
    plt.legend(title='Conc [M]',
               frameon=False)
    plt.xlim(((-2*R_Li_hyd-R_Li)*1e9*1.2, 1))
    plt.tight_layout()

# LiCl, vary K
plot_LiCl_K = False
if plot_LiCl_K:
    c = 1e-1
    pH = 5.8
    sns.set_palette('viridis')
    plt.figure('LiCl potential profile vary K')
    K_list = np.array([1e2, 1e1, 1e0, 1e-1])
    sigma_d_LiCl_K = np.zeros(len(K_list))
    i = 0
    for K in K_list:
        c_list1 = np.array([c+10**-pH, c, 10**-pH])
        K_list = np.array([K, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1_Li, C_2=C2_Li)
        sol.solver_PB()
        sigma_d_LiCl_K[i] = sol.sigma_d
        x = np.concatenate((np.array([-(R_Li+2*R_Li_hyd), -2*R_Li_hyd]), sol.x))
        psi = np.concatenate((np.array([sol.psi_0, sol.psi_beta]), sol.psi))
        plt.plot(x*1e9, psi*1000, label='%i' % np.log10(K))
        i += 1
    plt.plot(np.array([0, 0]), np.array([-700, 0]), color='k', linestyle='--', linewidth=1)
    plt.plot(np.array([-2*R_Li_hyd, -2*R_Li_hyd])*1e9, np.array([-700, 0]), color='k', linestyle='--', linewidth=1)
    plt.plot(np.array([-2*R_Li_hyd-R_Li, -2*R_Li_hyd-R_Li])*1e9, np.array([-700, 0]), color='k', linestyle='--', linewidth=1)
    plt.xlabel('x [nm]')
    plt.ylabel('$\psi$ [mV]')
    plt.legend(title='p$K_M$',
               frameon=False)
    plt.xlim(((-2*R_Li_hyd-R_Li)*1e9*1.2, 1))
    plt.tight_layout()


#########
# SIGMA #
#########

plot_sigma_d_conc = False
if plot_sigma_d_conc:
    plt.figure('sigma_d_conc')
    df = pd.read_csv(filepath_or_buffer='data_figures/sigma_d_conc.csv')
    sns.lineplot(x='c', y='beta_frac',
                 data=df)
    plt.tight_layout()
    plt.xscale('log')
    plt.xlabel('c [mol dm$^{-3}$]')
    plt.ylabel('$\sigma_d$ / |$\sigma_0$|')
    plt.text(3e-5, 0.04, 'LiCl\npH=5.8')

plot_sigma_beta_conc_H = False
if plot_sigma_beta_conc_H:
    plt.figure('sigma_beta_conc_H')
    df = pd.read_csv(filepath_or_buffer='data_figures/sigma_beta_conc_H.csv')
    sns.lineplot(x='c', y='beta_frac_H',
                 data=df,
                 label='H')
    ax = sns.lineplot(x='c', y='beta_frac_M',
                 data=df,
                 label='Li')
    ax.lines[1].set_linestyle("--")
    plt.tight_layout()
    plt.xscale('log')
    plt.xlabel('c [mol dm$^{-3}$]')
    plt.ylabel('$\sigma_{β}$ / |$\sigma_0$|')
    plt.text(1e-3, 0.5, 'LiCl\npH=5.8')
    plt.legend(frameon=False)

plot_grahame = True
if plot_grahame or plot_all:
    plt.figure('grahame')
    salts = ['LiCl', 'NaCl', 'KCl', 'CsCl']
    i = 0
    for salt in salts:
        df = pd.read_csv(filepath_or_buffer='data_figures/grahame_'+salt+'.csv')
        plt.plot(df.get('c'), df.get('sigma_d')/0.32,
                 color='C'+str(i),
                 label=salt)
        i += 1

    i = 0
    for salt in salts:
        df = pd.read_csv(filepath_or_buffer='data_figures/grahame_'+salt+'.csv')
        plt.plot(df.get('c'), df.get('sigma_d_grahame')/0.32,
                 color='C'+str(i),
                 linestyle='--',
                 label='(Grahame eqn) (model)')
        i += 1

    plt.legend(frameon=False,
               ncol=2)
    plt.xscale('log')
    plt.xlabel('c [mol dm$^{-3}$]')
    plt.ylabel('$\sigma_d$ / |$\sigma_0$|')
    plt.tight_layout()



plt.show()
