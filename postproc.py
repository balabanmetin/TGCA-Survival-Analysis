import pandas
from matplotlib import pyplot as plt
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

import sys

if __name__ == "__main__":

    kmf = KaplanMeierFitter()
    all_clin = pandas.read_csv(sys.argv[1], delimiter = "\t")
    z_rna = pandas.read_csv(sys.argv[2], delimiter = "\t")
    genes = pandas.read_csv("genes.csv", delimiter = "\t")
    epig = list(genes.iloc[1:]["HGNC approved symbol"])


    all_clin = all_clin[all_clin['new_death'].notnull()]
    """
    T = all_clin["new_death"]
    E = all_clin["death_event"]

    kmf.fit(T, event_observed = E)
    s = kmf.survival_function_
    kmf.median_
    kmf.plot()
    plt.title('Survival function of political regimes');
    plt.show()"""

    r, k = z_rna.shape
    print(k)
    for index, row in z_rna.iterrows():
        if index not in epig:
            continue
        upper = row[lambda x: x > np.median(row)]
        lower = row[lambda x: x <= np.median(row)]
        CU = all_clin.filter(items = upper.index, axis = 0)
        TU = CU["new_death"]
        EU = CU["death_event"]
        CL = all_clin.filter(items = lower.index, axis = 0)
        TL = CL["new_death"]
        EL = CL["death_event"]

        results = logrank_test(TU, TL, EU, EL, alpha=.95)
        if results.p_value < 0.001:
            print(results.p_value)
        else:
            continue
        fig = plt.figure()
        ax = fig.add_subplot(111)
        t = np.linspace(0, 4000, 1001)
        kmf.fit(TU, event_observed = EU, timeline = t, label = "High Expression")
        ax = kmf.plot(ax=ax)
        kmf.fit(TL, event_observed = EL, timeline = t,  label = "Low Expression")
        ax = kmf.plot(ax=ax)
        plt.ylim(0,1.05)
        plt.title("Kaplan-Meier plot of gene " + str(index))
        plt.text(0.9, 0.1, "p < {0:.4g}".format(results.p_value),
        horizontalalignment='center',
        verticalalignment='center',
        transform = ax.transAxes)
        plt.savefig('plots/'+ sys.argv[2] + index + '.png')
        plt.close()
        
        
