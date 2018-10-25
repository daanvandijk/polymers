#!/usr/bin/python3
import json

comparison = []
experiments = {}
with open('../comparison.json') as handler:
    comparison = json.load(handler)

with open('../experiments.json') as handler:
    experiments_raw = json.load(handler)
    for exp in experiments_raw['vec']:
        experiments[exp['title']] = exp

with open('comparisons.tex', 'w') as out:
    out.write('\\section*{Experiment comparisons}\n')

    for c in comparison:
        comparison_name = ','.join(c)
        l = ', '.join(c)

        out.write('\\begin{figure}[h]\n')
        out.write('\\centering\n')
        out.write('\\begin{minipage}{0.45\\textwidth}\n')
        out.write('\\includegraphics[width=8cm]{../img/comparisons/cm %s}\n' % comparison_name)
        out.write('\\caption{Center of mass comparison between experiments %s}\n' % l)
        out.write('\\label{fig:cm %s}\n' % comparison_name)
        out.write('\\end{minipage} \\hfil\n')
        out.write('\\begin{minipage}{0.45\\textwidth}\n')
        out.write('\\includegraphics[width=8cm]{../img/comparisons/mm %s}\n' % comparison_name)
        out.write('\\caption{Middle monomer comparison between %s}\n' % l)
        out.write('\\label{fig:mm %s}\n' % comparison_name)
        out.write('\\end{minipage}\n')
        out.write('\\end{figure}\n')

    out.write('\\begin{table}[h]\n')
    out.write('\\centering\n')
    out.write('\\footnotesize\n')
    out.write('\\begin{tabular}{l | l | l | l | l | l | l | l | l}\n')
    out.write('Title & $N_p$ & $N_t$ & $dt$ & $k_{BT}$ & $k$ & $a$ & $C$ & $l$ \\\\ \\hline \n')
    for title, exp in experiments.items():
        out.write('%s & %i & %i & %.1e & %.1e & %.1e & %.1e & %.1e & %.1e \\\\\n' % (
            exp['title'], 
            exp['Np'],
            exp['Nt'],
            exp['dt'],
            exp['k_BT'],
            exp['k'],
            exp['a'],
            exp['C'],
            exp['l']))

    out.write('\\end{tabular}\n')
    out.write('\\caption{Parameters of the experiments}\n')
    out.write('\\label{tab:experiment parameters}\n')
    out.write('\\end{table}\n')

for c in comparison:
    print(c)

