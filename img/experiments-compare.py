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

uniques = []
data = []

for c in comparison:
    for name in c:
        if (name not in uniques):
            uniques.append(name)

for name in uniques:
    print(name)
    pypolymer.read_experiment('../data/%s.dat'%(name))
