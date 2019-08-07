from line_marker import color_cycle, line_cycle
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Parameters
c = 0.6
pH = 5.8

# Set style
context = 'talk'
sns.set_context(context)
plt.rcParams['figure.figsize'] = 6.25, 4.5
sns.set_style({'xtick.direction': 'in', 'ytick.direction': 'in'})

# Lists
cation_list = ["Na"]

# Plot
for cation in cation_list:
    colors = color_cycle()
    lines = line_cycle()
    plt.figure(cation)
    filepath = "data_figures/isotherm_" + cation + "_model.csv"
    df = pd.read_csv(filepath_or_buffer=filepath)
    plt.plot(df.get("fraction"), df.get("gamma_M"),
             color=next(colors),
             linestyle=next(lines),
             label=r"$\Gamma$ (S$^-$" + cation + "$^+$)")
    plt.plot(df.get("fraction"), df.get("gamma_Li"),
             color=next(colors),
             linestyle=next(lines),
             label=r"$\Gamma$ (S$^-$Li$^+$)")
    plt.plot(df.get("fraction"), df.get("gamma_H"),
             color=next(colors),
             linestyle=next(lines),
             label=r"$\Gamma$ (S$^-$H$^+$)")
    plt.xlabel(r"[%s$^+$] / [Cl$^-$]" % cation)
    plt.ylabel(r"$\Gamma$ / CEC")
    plt.text(0.6, 0.4,
             "[Cl$^-$] = {0} M\npH = 5.8".format(c, pH))
    plt.legend(frameon=False)
    plt.ylim((-0.05, 1.05))
    plt.tight_layout()

plt.show()
