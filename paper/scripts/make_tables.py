from experiments import *
def read(name):return list(csv.DictReader((OUT/'results'/name).open(encoding='utf-8')))
sim=read('simulation_summary.csv');real=read('archive_metrics.csv')
text=r'''\begin{table}[htbp]\centering
\caption{Mean scalar RMSE across independent runs.}\label{tab:sim}
{\small\renewcommand{\baselinestretch}{1}\selectfont
\begin{tabular}{lrrrrrr}\hline
Method & Gaussian & Alt. & Outlier & Gaps & Combined & Maneuver\\\hline
'''
for m in METHODS:
    vals=[next(x for x in sim if x['case']==case and x['method']==m) for case in ['Gaussian','Alternating','Outliers','Gaps','Combined','Maneuver']]
    text+=m+' & '+' & '.join(f"{float(v['mean']):.3f}" for v in vals)+r'\\'+'\n'
text+=r'''\hline\end{tabular}}
\par\small Normalized observation units; 60 runs per case. Confidence intervals and all individual outcomes are provided in the supplementary CSV tables. Alt. denotes the alternating-disturbance case.
\end{table}
\begin{table}[htbp]\centering
\caption{Ship-data predictive disagreement and output roughness.}\label{tab:real}
{\small\renewcommand{\baselinestretch}{1}\selectfont
\begin{tabular}{lrrrr}\hline
&\multicolumn{2}{c}{Range (m)}&\multicolumn{2}{c}{Bearing (rad)}\\
Method & Next RMSE & Roughness & Next RMSE & Roughness\\\hline
'''
for m in METHODS:
    vs=[next(x for x in real if x['channel']==kind and x['method']==m) for kind in ['range','bearing']]
    text+=m+' & '+' & '.join(f"{float(v[k]):.5f}" for v in vs for k in ['next_rmse','roughness'])+r'\\'+'\n'
text+=r'''\hline\end{tabular}}
\par\small Next RMSE compares current output with the next adjacent observed epoch; roughness is the median absolute output increment. Neither is an independent truth error.
\end{table}
'''
(OUT/'sections/results_tables.tex').write_text(text,encoding='utf-8')
print('Tables generated from experiment CSV files')
