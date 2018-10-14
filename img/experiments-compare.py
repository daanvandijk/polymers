#!/usr/bin/python3
import pypolymer, numpy, json
import matplotlib.pyplot as plt
def cm_to_inch(cm):
    return cm / 2.54
plt.rc('text', usetex=True)

# pypolymer.get_experiments("../experiments.json")
# import the parameters of the experiments
comparison = []
with open('../comparison.json') as handler:
    comparison = json.load(handler)

# only read data once
data = {} 
for c in comparison:
    for name in c:
        if (name not in data.keys()):
            print("Loading experiment %s:" % name)
            data[name] = pypolymer.read_experiment('../data/%s.dat'%(name))

# now do the comparisons
for c in comparison:
    comparison_name = ','.join(c)
    print("Comparison: %s" % comparison_name)
    plt.figure(figsize=(cm_to_inch(8), cm_to_inch(6)), dpi=300)
    for name in c:
        plt.loglog(
                data[name].time[0, 1:],
                data[name].cm_sigma[0, 1:],
                label = name)
    plt.xlabel("t")
    plt.ylabel("\sigma_{\\text{cm}}")
    plt.legend()
    plt.savefig('comparisons/cm %s.png' % comparison_name, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(cm_to_inch(8), cm_to_inch(6)), dpi=300)
    for name in c:
        plt.loglog(
                data[name].time[0, 1:], 
                data[name].mm_sigma[0, 1:],
                label = name)
    plt.xlabel("t")
    plt.ylabel("\sigma_{\\text{mm}}")
    plt.legend()
    plt.savefig('comparisons/mm %s.png' % comparison_name, bbox_inches='tight')
    plt.close()
