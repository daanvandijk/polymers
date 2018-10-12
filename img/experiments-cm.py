#!/usr/bin/python3
import pypolymer, numpy, json
import matplotlib.pyplot as plt
def cm_to_inch(cm):
    return cm / 2.54
plt.rc('text', usetex=True)

# pypolymer.get_experiments("../experiments.json")
# import the parameters of the experiments
parameters = []
with open('../experiments.json') as handler:
    obj = json.load(handler)
    for p in obj['vec']:
        parameters.append(p)

# for each experiment, make a plot of cm
for p in parameters:
    path = '../data/%s.dat' % (p['title'])
    img_path = 'experiments-cm-%s.png' % (p['title'])
    data = pypolymer.read_experiment(path)
    print(path)

    plt.figure(figsize=(cm_to_inch(8), cm_to_inch(6)), dpi=300)
    plt.loglog(data.time[0, 1:], data.cm_sigma[0, 1:])
    T = [data.time[0,1], data.time[0,-1]]
    Y = numpy.power(T, 0.5) * (data.cm_sigma[0,1]/numpy.sqrt(data.time[0,1]))
    plt.loglog(T, Y, 'k--') 
    plt.xlabel("t")
    plt.ylabel("\sigma_{\\text{cm}}")
    plt.title("Rouse model")
    # plt.legend()
    plt.savefig(img_path, bbox_inches='tight')
    plt.close()
