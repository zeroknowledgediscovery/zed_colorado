from pathlib import Path
import re

src=Path('/mnt/data/full_theory_manuscript/original_theory.tex').read_text()

# Preamble: add PGFPlots support.
src=src.replace('\\usepackage{booktabs}\n', '\\usepackage{booktabs}\n\\usepackage{array}\n\\usepackage{graphicx}\n\\usepackage{pgfplots,pgfplotstable}\n\\usepgfplotslibrary{groupplots}\n\\usepackage{tikz}\n\\pgfplotsset{compat=1.18}\n')

abstract=r'''\begin{abstract}
Genomic risk and longitudinal clinical history are distinct information channels for forecasting future disease. We study when a static germline genomic channel can add predictive value once a longitudinal clinical-history score has become strongly informative. The central object is the residual genomic likelihood ratio conditional on clinical history. If this residual genomic evidence is uniformly bounded while sequential clinical evidence can develop an extended likelihood-ratio tail, then two consequences follow. First, genomic information can change a likelihood-ratio decision only inside a bounded neighborhood of the clinical decision threshold; outside that region the clinical evidence is already too strong for bounded genomic evidence to reverse the decision. Second, at sufficiently extreme specificity, the clinical likelihood-ratio test strictly dominates every genomic-only rule whose likelihood ratio is bounded. We also show that any deterministic genomic model inherits the bound and that decision-local benefit can coexist with negligible global AUC change. We evaluate these predictions in a Colorado fibrotic interstitial lung disease (ILD) analysis using ZeBRA, a longitudinal EHR-derived score, and a curated genomic panel containing \textit{MUC5B} rs35705950. Across repeated held-out splits, ZeBRA dominates global discrimination, broad genomic augmentation provides essentially no global FILD/FILA gain after conditioning on ZeBRA, and ZeBRA predicts \textit{MUC5B} carrier status only at chance level. Nevertheless, \textit{MUC5B} retains significant conditional disease information in restricted ZeBRA operating bands, and a two-threshold genomic rescue rule increases held-out sensitivity by up to approximately 4.5 percentage points at essentially matched false-positive rate. The empirical pattern is therefore not global genomic redundancy but regime-dependent complementarity localized near the clinical decision boundary.
\end{abstract}'''
src=re.sub(r'\\begin\{abstract\}.*?\\end\{abstract\}', lambda m: abstract, src, count=1, flags=re.S)

# Update motivating empirical paragraph but keep biological context.
old_start='Our recent fibrotic ILD analyses ask whether genomic information improves prediction once ZeBRA is known.'
old_end='These observations are exploratory and the particular percentile bands were selected post hoc; they motivate the theory below rather than constitute independent confirmation of it.'
start=src.index(old_start)
end=src.index(old_end,start)+len(old_end)
new_para=r'''Our recent fibrotic ILD analyses ask whether genomic information improves prediction once ZeBRA is known. In 30 repeated held-out splits of the current operational FILD/FILA endpoint, mean AUC is approximately 0.814 for ZeBRA, 0.566 for the selected genomic panel, and 0.799 for a nonlinear combined model. Because that comparison alone does not isolate incremental information, we also performed a paired regularized analysis in which ZeBRA-only and ZeBRA-plus-genomics models were evaluated on the same held-out subjects. For FILD/FILA, the mean paired AUC increment from the complete genomic panel is only $+0.00043$, the median increment is 0, and the split-wise 2.5th--97.5th percentile range spans approximately $-0.0051$ to $+0.0041$; average precision decreases and both Brier score and log loss worsen after broad genomic augmentation.

The absence of global gain is not explained by ZeBRA reconstructing the major genomic signal. Among subjects with called rs35705950 genotype, ZeBRA predicts \textit{MUC5B} carrier status at chance level (AUC approximately 0.509). A modest pooled enrichment of carriers in the extreme ZeBRA tail disappears after stratification by FILD/FILA status, showing that the pooled association is consistent with both channels tracking disease rather than with clinical history encoding genotype. Yet residual \textit{MUC5B} disease information is strongly heterogeneous across ZeBRA score space. In the ZeBRA band corresponding approximately to 5--10\% control FPR, adding \textit{MUC5B} after ZeBRA yields a likelihood-ratio statistic of 16.55 ($p\approx4.7\times10^{-5}$). A two-threshold rule that consults \textit{MUC5B} only between a stringent upper and lower ZeBRA threshold increases held-out sensitivity by up to approximately 4.5 percentage points while preserving essentially the same realized overall FPR. These observations motivate and correspond to the theory below; they are not taken as proofs of its regularity assumptions.'''
src=src[:start]+new_para+src[end:]

# A1: make boundedness an explicit assumption, not an equivalence implied by non-Mendelian architecture.
a1_old=re.compile(r'\\textbf\{A1 \(Non-Mendelian genomic evidence\)\.\}.*?is bounded above by \$M_G<\\infty\$\.', re.S)
a1_new=r'''\textbf{A1 (Uniformly bounded residual genomic evidence).} There exist constants $0<m\le 1\le M<\infty$ such that
\begin{equation}
0<m\le \LR_{G\mid C}(g;c)\le M<\infty
\label{eq:boundedconditional}
\end{equation}
for the relevant support of $(g,c)$, and the marginal genomic likelihood ratio
\begin{equation}
\LR_G(g)=\frac{p(g\mid Y=1)}{p(g\mid Y=0)}
\end{equation}
is bounded above by $M_G<\infty$. Complex, non-Mendelian disease architecture motivates finite residual genomic Bayes factors, but does not by itself imply a uniform bound over all conditioned clinical states; boundedness is therefore an explicit regularity assumption.'''
src=a1_old.sub(lambda m: a1_new,src,count=1)

# Insert an AUC-localization proposition after the second corollary proof.
marker='''A conditional expectation of a random variable bounded above by $M_G$ is itself bounded above by $M_G$. The analogous lower bound follows if $\\LR_G$ is bounded away from zero.\n\\end{proof}\n\n\\section{Why the clinical channel can have an extended tail}'''
prop=r'''A conditional expectation of a random variable bounded above by $M_G$ is itself bounded above by $M_G$. The analogous lower bound follows if $\LR_G$ is bounded away from zero.
\end{proof}

\begin{proposition}[Local decision gain can coexist with negligible global AUC change]
\label{prop:auc-local}
Let $S(C)$ be a clinical score and $S'(C,G)$ an augmented score. Suppose $S'=S$ except on a measurable region $B$ of the clinical--genomic state space. Then
\begin{equation}
\left|\operatorname{AUC}(S')-\operatorname{AUC}(S)\right|
\le \Pone(B)+\Pzero(B)-\Pone(B)\Pzero(B)
\le \Pone(B)+\Pzero(B).
\label{eq:auc-local-bound}
\end{equation}
Hence an augmented rule can produce a substantial change in decisions within a small boundary region while producing an arbitrarily small change in global AUC as the case and control mass of that region tends to zero.
\end{proposition}

\begin{proof}
AUC is the probability that a randomly drawn case receives a higher score than a randomly drawn control, with the usual half-credit convention for ties. If neither member of a case--control pair lies in $B$, then $S'=S$ for both members and the pairwise ranking contribution is unchanged. Therefore the AUC contribution can differ only when at least one member of the pair lies in $B$. The probability of that event is
\begin{align}
1-(1-\Pone(B))(1-\Pzero(B))
&=\Pone(B)+\Pzero(B)-\Pone(B)\Pzero(B),
\end{align}
which bounds the maximum possible AUC change.
\end{proof}

\section{Why the clinical channel can have an extended tail}'''
if marker not in src:
    raise RuntimeError('marker for proposition not found')
src=src.replace(marker,prop,1)

# Replace Interpretation through Scope with a stronger interpretation + complete empirical section.
sec_start=src.index('\\section{Interpretation of the ILD/IPF observations}')
sec_end=src.index('\\section{Scope and limitations}',sec_start)
interp_emp=r'''\section{Interpretation of the theoretical result}

The theorem and Proposition~\ref{prop:auc-local} make four distinct claims that should not be conflated.

\subsection{Global discrimination is not the quantity controlled by the theorem}
Theorem~\ref{thm:main}(ii) is about whether bounded residual genomic evidence can change a threshold decision, not about global ranking. If genomic information matters only in the boundary band \eqref{eq:boundaryband}, then only a restricted fraction of case--control pairs can be reordered. Proposition~\ref{prop:auc-local} formalizes the consequence: a meaningful local operating-point improvement can coexist with a negligible global AUC change. Thus ``no AUC gain'' and ``local genomic utility'' are mathematically compatible rather than contradictory.

\subsection{Weak ZeBRA--genotype association does not imply zero conditional genomic disease information}
The relevant residual quantity is $\LR_{G\mid C}$, not the marginal dependence between the clinical score and genotype. ZeBRA can be nearly uninformative about rs35705950 carrier status while \textit{MUC5B} still changes disease odds conditional on a particular clinical state. The channels can therefore be statistically distinct yet jointly useful near selected boundaries.

\subsection{Decision utility is expected to move with the operating threshold}
Corollary~\ref{cor:zebra} shows that the genomic decision-relevant region is defined relative to the chosen ZeBRA threshold. There is no unique genomic ``useful band'' independent of the decision problem. Changing the tolerated false-positive rate changes the clinical boundary and therefore changes which patients can be reclassified by bounded genomic evidence.

\subsection{Clinical-tail dominance is a high-specificity statement, not an AUC statement}
Theorem~\ref{thm:main}(iii)--(v) concerns the likelihood-ratio tail. A genomic-only rule has LR$+$ bounded by $M_G$, whereas a sufficiently extreme clinical likelihood-ratio threshold has LR$+$ at least $K$ with FPR at most $1/K$. The empirical quantities most relevant to this part of the theory are therefore stringent-specificity likelihood ratios and low-FPR operating behavior, not global AUC by itself.

\section{Empirical correspondence in the Colorado ILD analysis}
\label{sec:empirical}

\subsection{Cohort and analysis hierarchy}
The current analysis uses the Colorado genomic--clinical cohort assembled in \texttt{zebra\_comp}. Three stored phenotype fields are evaluated: \texttt{FILD or FILA ADJUDICATED}, \texttt{NONFIBROTIC ILD ADJUDICATED}, and \texttt{NON FIBROTIC ILA ADJUDICATED}. The historical notebook construction maps \texttt{Y/y} to 1 and \texttt{N/n} to 0, converts the field to numeric, assigns missing target values to 0, and excludes observations without an available 104-week ZeBRA score. The reported endpoint is therefore the notebook-defined operational target; not every target-zero observation is an independently adjudicated negative. This distinction matters especially when interpreting FPR.

The empirical analyses are intentionally hierarchical. Notebook~32 asks whether genomics alone or broad nonlinear ZeBRA+genomics improves global held-out discrimination. Notebook~33 asks the more specific paired question of incremental genomic value conditional on ZeBRA. Separate analyses then test whether ZeBRA reconstructs \textit{MUC5B}, where residual \textit{MUC5B} information lies in ZeBRA score space, and whether a two-threshold genomic rule improves boundary decisions at matched FPR.

\begin{table*}[!t]
\caption{Theory--evidence correspondence. The empirical analyses illustrate predicted consequences; they do not prove A1 or A2.}
\label{tab:theory-evidence-map}
\centering\footnotesize\setlength{\tabcolsep}{5pt}
\begin{tabular}{p{0.23\textwidth}p{0.31\textwidth}p{0.37\textwidth}}
\toprule
Theory & Analysis & Empirical meaning\\
\midrule
Distinct channels and residual genomic evidence & ZeBRA-to-\textit{MUC5B} pooled and FILD-stratified analyses & Weak global genomic increment is not explained by ZeBRA reconstructing genotype.\\
Boundary localization, Theorem~\ref{thm:main}(ii) & Local ZeBRA--\textit{MUC5B} information curves and FPR-defined bands & Conditional genomic disease information is heterogeneous across clinical-score space.\\
Local gain with little global AUC change, Proposition~\ref{prop:auc-local} & Notebook~32 and notebook~33 plus matched-FPR rescue & Global augmentation yields essentially no FILD AUC gain even though selected boundary decisions improve.\\
Clinical-tail dominance, Theorem~\ref{thm:main}(iii)--(v) & Low-FPR sensitivity and LR$+$ behavior & Stringent clinical operating points can dominate genomic-only evidence even when genomics remains locally useful.\\
\bottomrule
\end{tabular}
\end{table*}

\pgfplotstableread[col sep=comma]{data/global_auc_summary.csv}\globalauc
\pgfplotstableread[col sep=comma]{data/incremental_summary.csv}\incremental
\pgfplotstableread[col sep=comma]{data/local_info_bands.csv}\localbands
\pgfplotstableread[col sep=comma]{data/boundary_rescue_summary_ordered_ext.csv}\rescuedata
\pgfplotstableread[col sep=comma]{data/low_fpr_operating_points.csv}\lowfprdata
\pgfplotstableread[col sep=comma]{data/low_fpr_upper_envelope.csv}\lowenv

\subsection{Global repeated-split discrimination}
Across 30 repeated held-out outer splits, ZeBRA is the strongest predictor for all three endpoints. For FILD/FILA, mean AUC is 0.8142 for ZeBRA, 0.5662 for genomics alone, and 0.7988 for the nonlinear combined model. The same qualitative ordering is stronger for the two nonfibrotic endpoints. These results establish the global pattern, but do not by themselves test whether genomics has a small conditional increment after ZeBRA is known.

\begin{figure*}[!t]
\centering
\begin{tikzpicture}
\begin{axis}[width=0.91\textwidth,height=0.31\textwidth,
ylabel={Held-out ROC AUC},ymin=0.35,ymax=0.90,
xtick={0,1,2},xticklabels={FILD/FILA,Nonfibrotic ILD,Nonfibrotic ILA},
xticklabel style={font=\footnotesize},yticklabel style={font=\footnotesize},label style={font=\footnotesize},
legend style={at={(0.5,1.03)},anchor=south,legend columns=3,draw=none,font=\footnotesize},
ymajorgrids=true,major grid style={dotted},enlarge x limits=0.18]
\addplot+[only marks,mark=*,mark size=2.5pt,error bars/.cd,y dir=both,y explicit] coordinates {(0,0.8142)+-(0,0.0089) (1,0.8202)+-(0,0.0257) (2,0.8281)+-(0,0.0162)};\addlegendentry{ZeBRA}
\addplot+[only marks,mark=square*,mark size=2.5pt,error bars/.cd,y dir=both,y explicit] coordinates {(0,0.5662)+-(0,0.0245) (1,0.5141)+-(0,0.0719) (2,0.4717)+-(0,0.0443)};\addlegendentry{Genomic}
\addplot+[only marks,mark=triangle*,mark size=2.8pt,error bars/.cd,y dir=both,y explicit] coordinates {(0,0.7988)+-(0,0.0150) (1,0.6145)+-(0,0.0953) (2,0.7112)+-(0,0.0541)};\addlegendentry{Combined}
\end{axis}\end{tikzpicture}
\caption{Global repeated-split discrimination from notebook~32. Points are mean held-out AUC over 30 outer splits; bars indicate one standard deviation. The clinical-history channel dominates globally, and broad nonlinear combination does not improve ZeBRA.}
\label{fig:global-auc}
\end{figure*}

\begin{table*}[!t]
\caption{Held-out AUC over 30 repeated outer splits from notebook~32. Values are mean $\pm$ SD.}
\label{tab:auc-repeated}\centering\footnotesize
\begin{tabular}{lccc}\toprule
Target & ZeBRA & Genomic & Combined\\\midrule
FILD/FILA & $0.8142\pm0.0089$ & $0.5662\pm0.0245$ & $0.7988\pm0.0150$\\
Nonfibrotic ILD & $0.8202\pm0.0257$ & $0.5141\pm0.0719$ & $0.6145\pm0.0953$\\
Nonfibrotic ILA & $0.8281\pm0.0162$ & $0.4717\pm0.0443$ & $0.7112\pm0.0541$\\\bottomrule
\end{tabular}
\end{table*}

\subsection{Paired incremental value conditional on ZeBRA}
Notebook~33 evaluates ZeBRA-only and ZeBRA-plus-genomics models on the same held-out subjects in 30 repeated splits. For FILD/FILA, the mean paired AUC difference is only $+0.00043$, the median is 0, and the split-wise 2.5th--97.5th percentile interval is $[-0.00507,0.00408]$. Average precision decreases, while Brier score and log loss worsen. The conclusion is therefore not merely that one nonlinear combined model underperformed, but that broad genomic augmentation supplies no reproducible global FILD/FILA discrimination gain conditional on ZeBRA.

\begin{figure*}[!t]
\centering
\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.45cm},width=0.43\textwidth,height=0.28\textwidth,
tick label style={font=\footnotesize},label style={font=\footnotesize},ymajorgrids=true,major grid style={dotted}]
\nextgroupplot[ylabel={Paired $\Delta$AUC},ymin=-0.24,ymax=0.03,xtick={1,2,3},xticklabels={FILD/FILA,NF-ILD,NF-ILA},xticklabel style={rotate=20,anchor=east}]
\addplot+[only marks,mark=*,mark size=2.5pt] table[col sep=comma,x expr=\coordindex+1,y=delta_auc_mean]{data/incremental_summary.csv};
\addplot+[dashed,no marks] coordinates {(0.5,0) (3.5,0)};
\draw (axis cs:1,-0.0050703)--(axis cs:1,0.0040843);\draw (axis cs:0.94,-0.0050703)--(axis cs:1.06,-0.0050703);\draw (axis cs:0.94,0.0040843)--(axis cs:1.06,0.0040843);
\draw (axis cs:2,-0.2257717)--(axis cs:2,-0.0050563);\draw (axis cs:1.94,-0.2257717)--(axis cs:2.06,-0.2257717);\draw (axis cs:1.94,-0.0050563)--(axis cs:2.06,-0.0050563);
\draw (axis cs:3,-0.0738433)--(axis cs:3,-0.0070462);\draw (axis cs:2.94,-0.0738433)--(axis cs:3.06,-0.0738433);\draw (axis cs:2.94,-0.0070462)--(axis cs:3.06,-0.0070462);
\nextgroupplot[ylabel={Mean paired metric change (FILD/FILA)},ymin=-0.04,ymax=0.01,xtick={1,2,3},xticklabels={$\Delta$AP,$\Delta$Brier,$\Delta$LogLoss}]
\addplot+[ybar,bar width=15pt] coordinates {(1,-0.0311) (2,0.000312) (3,0.00468)};\addplot+[dashed,no marks] coordinates {(0.5,0) (3.5,0)};
\end{groupplot}\end{tikzpicture}
\caption{Paired incremental-value analysis from notebook~33. Left: mean paired $\Delta$AUC from adding the complete genomic panel to ZeBRA, with split-wise 2.5th--97.5th percentile intervals. FILD/FILA is centered essentially at zero. Right: FILD/FILA average precision, Brier score, and log loss changes all move in the unfavorable direction.}
\label{fig:incremental}
\end{figure*}

\subsection{ZeBRA does not reconstruct \textit{MUC5B}}
Among subjects with called rs35705950 genotype, ZeBRA predicts carrier status at essentially chance level: pooled AUC 0.5088, explicit FILD/FILA-negative AUC 0.4847, and FILD/FILA-positive AUC 0.4817. Yet carrier prevalence is approximately 39.3\% in FILD/FILA-positive subjects versus 19.2\% among explicit negatives, and the pooled top 1\% ZeBRA tail shows 1.47-fold enrichment. The correct interpretation is therefore that both ZeBRA and \textit{MUC5B} are associated with disease, not that ZeBRA has reconstructed genotype.

\begin{figure*}[!t]\centering\begin{tikzpicture}
\begin{groupplot}[group style={group size=3 by 1,horizontal sep=1.1cm},width=0.295\textwidth,height=0.27\textwidth,tick label style={font=\footnotesize},label style={font=\footnotesize},ymajorgrids=true,major grid style={dotted}]
\nextgroupplot[ylabel={AUC for MUC5B carrier},ymin=0.45,ymax=0.53,xtick={0,1,2},xticklabels={Pooled,Explicit neg.,FILD+},xticklabel style={rotate=22,anchor=east}]
\addplot+[ybar,bar width=16pt] coordinates {(0,0.5088) (1,0.4847) (2,0.4817)};\addplot+[dashed,no marks] coordinates {(-0.5,0.5) (2.5,0.5)};
\nextgroupplot[ylabel={MUC5B T-carrier prevalence (\%)},ymin=0,ymax=45,xtick={0,1,2},xticklabels={Pooled,Explicit neg.,FILD+},xticklabel style={rotate=22,anchor=east}]
\addplot+[ybar,bar width=16pt] coordinates {(0,19.58) (1,19.17) (2,39.29)};
\nextgroupplot[ylabel={Fold enrichment},ymin=0,ymax=1.8,xtick={0},xticklabels={Top 1\% ZeBRA tail},xticklabel style={rotate=15,anchor=east}]
\addplot+[ybar,bar width=22pt] coordinates {(0,1.47)};\node[font=\scriptsize] at (axis cs:0,1.62) {$p\approx0.009$};
\end{groupplot}\end{tikzpicture}
\caption{The weak pooled ZeBRA--\textit{MUC5B} relationship is disease-mediated rather than genotype reconstruction. Left: ZeBRA-to-carrier AUC is at chance in pooled and FILD-stratified analyses. Middle: carrier prevalence is much higher among FILD/FILA-positive subjects. Right: the pooled extreme ZeBRA tail is modestly enriched for carriers.}
\label{fig:muc5b}
\end{figure*}

\subsection{Residual genomic information is localized}
The local-information analysis directly addresses the quantity suggested by Theorem~\ref{thm:main}(ii). Conditional \textit{MUC5B} evidence varies strongly across ZeBRA-defined control-FPR bands. In the 5--10\% band, the likelihood-ratio statistic for adding \textit{MUC5B} after ZeBRA is 16.55 ($p=4.74\times10^{-5}$); in the 20--40\% band it is 13.69 ($p=2.15\times10^{-4}$). This is not evidence that \textit{MUC5B} globally outperforms ZeBRA. It is evidence that residual genomic disease information is spatially heterogeneous in clinical-score space.

\begin{figure*}[!t]\centering\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.45cm},width=0.43\textwidth,height=0.28\textwidth,tick label style={font=\footnotesize},label style={font=\footnotesize},ymajorgrids=true,major grid style={dotted}]
\nextgroupplot[ylabel={Conditional MUC5B LRT after ZeBRA},ymin=0,ymax=19,xtick={1,2,3,4},xticklabels={2--5,5--10,10--20,20--40},xlabel={Control-FPR band (\%)}]
\addplot+[ybar,bar width=15pt] table[col sep=comma,x=band_index,y=lrt]{data/local_info_bands.csv};\addplot+[dashed,no marks] coordinates {(0.5,3.84) (4.5,3.84)};
\nextgroupplot[ylabel={Local ROC AUC},ymin=0.3,ymax=0.75,xtick={1,2,3,4},xticklabels={2--5,5--10,10--20,20--40},xlabel={Control-FPR band (\%)},legend style={at={(0.5,1.03)},anchor=south,legend columns=2,draw=none,font=\scriptsize}]
\addplot+[mark=*,thick] table[col sep=comma,x=band_index,y=auc_zebra]{data/local_info_bands.csv};\addlegendentry{ZeBRA}
\addplot+[mark=square*,thick] table[col sep=comma,x=band_index,y=auc_muc5b]{data/local_info_bands.csv};\addlegendentry{MUC5B}
\end{groupplot}\end{tikzpicture}
\caption{Localized genomic disease information. Left: conditional \textit{MUC5B} likelihood-ratio statistics after ZeBRA; the dashed line is the $\chi^2_1=3.84$ reference. Right: local discrimination varies by operating region, with \textit{MUC5B} stronger than ZeBRA in selected bands.}
\label{fig:local-info}
\end{figure*}

\begin{table*}[!t]\caption{Local FILD/FILA information in ZeBRA-defined control-FPR bands.}\label{tab:local-info}\centering\footnotesize
\begin{tabular}{lrrrrrr}\toprule
Band & $N$ & Cases & LRT $G\mid Z$ & $p$ & AUC ZeBRA & AUC \textit{MUC5B}\\\midrule
2--5\% & 401 & 26 & 4.44 & 0.0352 & 0.577 & 0.595\\
5--10\% & 652 & 26 & 16.55 & $4.74\times10^{-5}$ & 0.367 & 0.677\\
10--20\% & 1,281 & 30 & 4.77 & 0.0289 & 0.558 & 0.586\\
20--40\% & 2,540 & 37 & 13.69 & $2.15\times10^{-4}$ & 0.516 & 0.635\\\bottomrule
\end{tabular}\end{table*}

\subsection{Boundary-targeted genomics improves the operating decision}
A two-threshold rescue rule consults \textit{MUC5B} only for subjects lying between an upper and a lower ZeBRA threshold. Across all six evaluated threshold pairs, held-out sensitivity increases while realized overall FPR remains essentially unchanged. For the 0.5\%/10\% pair, sensitivity rises from 0.3587 to 0.4033 ($+0.0446$), FPR changes from 0.02335 to 0.02310, and LR$+$ increases from 15.58 to 17.54. This is the most direct empirical counterpart of boundary-localized utility.

\begin{figure*}[!t]\centering\begin{tikzpicture}
\begin{groupplot}[group style={group size=3 by 1,horizontal sep=1.1cm},width=0.30\textwidth,height=0.27\textwidth,tick label style={font=\footnotesize},label style={font=\footnotesize},ymajorgrids=true,major grid style={dotted}]
\nextgroupplot[ylabel={Held-out sensitivity},ymin=0.28,ymax=0.42,xtick={0,1,2,3,4,5},xticklabels={0.5/2,1/2,0.5/5,1/5,0.5/10,1/10},xticklabel style={rotate=34,anchor=east},xlabel={Upper/lower target FPR (\%)},legend style={at={(0.5,1.03)},anchor=south,legend columns=2,draw=none,font=\scriptsize}]
\addplot+[mark=*,thick] table[col sep=comma,x expr=\coordindex,y=zebra_sens]{data/boundary_rescue_summary_ordered_ext.csv};\addlegendentry{Matched ZeBRA}
\addplot+[mark=square*,thick] table[col sep=comma,x expr=\coordindex,y=rescue_sens]{data/boundary_rescue_summary_ordered_ext.csv};\addlegendentry{MUC5B rescue}
\nextgroupplot[xlabel={$\Delta$FPR},ylabel={$\Delta$sensitivity},xmin=-0.00032,xmax=0.00006,ymin=0,ymax=0.05,scaled x ticks=false,xtick={-0.0003,-0.0002,-0.0001,0},xticklabel style={/pgf/number format/fixed,/pgf/number format/precision=4,rotate=28,anchor=east}]
\addplot+[only marks,mark=*,mark size=2.4pt] table[col sep=comma,x=delta_fpr,y=delta_sens]{data/boundary_rescue_summary_ordered_ext.csv};\addplot+[dashed,no marks] coordinates {(-0.00032,0) (0.00006,0)};\addplot+[dashed,no marks] coordinates {(0,0) (0,0.05)};
\nextgroupplot[ylabel={LR$+$},ymin=12,ymax=42,xtick={0,1,2,3,4,5},xticklabels={0.5/2,1/2,0.5/5,1/5,0.5/10,1/10},xticklabel style={rotate=34,anchor=east},xlabel={Upper/lower target FPR (\%)}]
\addplot+[mark=*,thick] table[col sep=comma,x expr=\coordindex,y=zebra_lrplus]{data/boundary_rescue_summary_ordered_ext.csv};\addplot+[mark=square*,thick] table[col sep=comma,x expr=\coordindex,y=rescue_lrplus]{data/boundary_rescue_summary_ordered_ext.csv};
\end{groupplot}\end{tikzpicture}
\caption{Matched-FPR boundary rescue. Left: sensitivity of matched ZeBRA and boundary-targeted \textit{MUC5B} rescue. Center: every evaluated pair has positive sensitivity gain while $\Delta$FPR is near zero. Right: LR$+$ also increases under rescue.}
\label{fig:boundary-rescue}
\end{figure*}

\begin{table*}[!t]\caption{Selected matched-FPR boundary-rescue results over 100 repeated splits.}\label{tab:boundary-selected}\centering\footnotesize
\begin{tabular}{lrrrrrrrr}\toprule
Upper/lower & ZeBRA sens. & Rescue sens. & $\Delta$ sens. & ZeBRA FPR & Rescue FPR & ZeBRA LR$+$ & Rescue LR$+$ & Rescue LR$-$\\\midrule
1\%/5\% & 0.3303 & 0.3528 & +0.0225 & 0.01844 & 0.01826 & 18.19 & 19.52 & 0.659\\
0.5\%/5\% & 0.3121 & 0.3487 & +0.0366 & 0.01404 & 0.01400 & 22.74 & 25.21 & 0.661\\
1\%/10\% & 0.3740 & 0.4075 & +0.0335 & 0.02753 & 0.02735 & 13.75 & 14.95 & 0.609\\
0.5\%/10\% & 0.3587 & 0.4033 & +0.0446 & 0.02335 & 0.02310 & 15.58 & 17.54 & 0.611\\\bottomrule
\end{tabular}\end{table*}

\subsection{Low-FPR operating geometry}
The low-FPR view is most relevant to the clinical-tail part of the theorem. Across the evaluated matched operating points, rescue shifts sensitivity upward with very small FPR changes and systematically increases LR$+$. The pointwise envelope in Fig.~\ref{fig:low-fpr} is descriptive rather than a third classifier; it summarizes the best sensitivity achieved among the evaluated operating rules.

\begin{figure*}[!t]\centering\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.35cm},width=0.43\textwidth,height=0.29\textwidth,tick label style={font=\footnotesize},label style={font=\footnotesize},xmajorgrids=true,ymajorgrids=true,major grid style={dotted}]
\nextgroupplot[xlabel={False-positive rate},ylabel={Sensitivity},xmin=0,xmax=0.03,ymin=0.27,ymax=0.43,legend style={at={(0.5,1.03)},anchor=south,legend columns=3,draw=none,font=\scriptsize}]
\addplot+[only marks,mark=*,mark size=2.4pt] coordinates {(0.00801,0.2919) (0.01227,0.3046) (0.01404,0.3121) (0.01844,0.3303) (0.02335,0.3587) (0.02753,0.3740)};\addlegendentry{Matched ZeBRA}
\addplot+[only marks,mark=square*,mark size=2.4pt] coordinates {(0.00802,0.3101) (0.01228,0.3142) (0.01400,0.3487) (0.01826,0.3528) (0.02310,0.4033) (0.02735,0.4075)};\addlegendentry{MUC5B rescue}
\addplot+[dashed,thick,no marks] table[col sep=comma,x=fpr,y=tpr]{data/low_fpr_upper_envelope.csv};\addlegendentry{Pointwise envelope}
\nextgroupplot[xlabel={False-positive rate},ylabel={LR$+$},xmin=0.007,xmax=0.0285,ymin=12,ymax=42]
\addplot+[mark=*,thick] coordinates {(0.00801,38.16) (0.01227,25.43) (0.01404,22.74) (0.01844,18.19) (0.02335,15.58) (0.02753,13.75)};
\addplot+[mark=square*,thick] coordinates {(0.00802,39.67) (0.01228,26.03) (0.01400,25.21) (0.01826,19.52) (0.02310,17.54) (0.02735,14.95)};
\end{groupplot}\end{tikzpicture}
\caption{Low-FPR operating geometry. Boundary-targeted rescue occupies a favorable low-FPR region and increases LR$+$ without requiring improved global AUC.}
\label{fig:low-fpr}
\end{figure*}

'''
src=src[:sec_start]+interp_emp+src[sec_end:]

# Replace scope and conclusion.
scope_start=src.index('\\section{Scope and limitations}')
concl_start=src.index('\\section{Conclusion}',scope_start)
new_scope=r'''\section{Scope and limitations}

The theory is a statement about information channels, not about a particular machine-learning architecture. It does not require independence between genomics and clinical history; dependence is retained explicitly through $\LR_{G\mid C}$. It also applies to any deterministic genomic score derived from the measured genomic state through Corollary~2.

Several qualifications are essential. First, A1 is an explicit uniform boundedness assumption; non-Mendelian disease architecture motivates finite genomic evidence but does not prove a uniform bound after conditioning on every clinical state. Second, A2 is also an assumption. Sequential observations make likelihood-ratio accumulation possible, but high AUC alone does not establish an unbounded or sufficiently extended clinical tail. The empirical quantities most relevant to A2 are calibration, score-stratum likelihood ratios, stability of extreme-tail LR$+$, and fixed-specificity or partial-AUC behavior. Third, Theorem~\ref{thm:main}(ii) concerns threshold changes, not posterior-probability changes: genomics may still shift posterior risk away from the boundary even when it cannot reverse the decision. Fourth, the operational FILD/FILA target in the historical notebook construction maps missing adjudication values to zero; consequently, the target-zero group should not be described as a fully adjudicated control cohort without further confirmation. Finally, the particular genomic rescue rule remains a modeling result that requires independent or prospectively prespecified validation. The empirical analyses support the predicted geometry but do not validate A1 or A2 by themselves.

A stronger future validation of the theory should therefore test ZeBRA calibration, convert score strata into empirical likelihood ratios, quantify uncertainty in extreme-tail LR$+$, evaluate prespecified fixed-specificity thresholds, and replicate boundary-localized genomic rescue in an external cohort.'''
src=src[:scope_start]+new_scope+'\n\n'+src[concl_start:]

concl_start=src.index('\\section{Conclusion}')
bib_start=src.index('\\begin{thebibliography}',concl_start)
new_concl=r'''\section{Conclusion}

The central result is not that clinical history ``contains'' genomics, nor that genomics is globally useless. The two modalities are distinct evidence channels. Bayes' rule factorizes their joint evidence into a clinical likelihood ratio and a residual genomic likelihood ratio conditional on clinical state. If the residual genomic factor is bounded while sequential clinical evidence can extend into a strong likelihood-ratio tail, then genomic information can reverse a clinical threshold decision only inside a bounded neighborhood of that threshold, and sufficiently extreme clinical evidence must eventually dominate every bounded genomic-only rule. Any deterministic genomic predictor inherits the same finite-evidence constraint, while Proposition~\ref{prop:auc-local} explains why clinically meaningful boundary reclassification can coexist with almost no change in global AUC.

The current Colorado results match this qualitative geometry. ZeBRA dominates global discrimination; broad genomic augmentation provides essentially zero incremental FILD/FILA AUC after conditioning on ZeBRA; ZeBRA does not reconstruct \textit{MUC5B}; residual \textit{MUC5B} disease information is concentrated in selected ZeBRA operating bands; and a boundary-targeted rescue rule improves sensitivity and LR$+$ at essentially matched overall FPR. The resulting picture is one of \emph{regime-dependent complementarity}: genomics is most useful where the clinical evidence is close enough to the action threshold that a finite residual Bayes factor can still change the decision, whereas longitudinal clinical evidence dominates in the extreme tail.'''
src=src[:concl_start]+new_concl+'\n\n'+src[bib_start:]

Path('/mnt/data/full_theory_manuscript/zebra_genomics_boundary_full_theory_pgf.tex').write_text(src)
print('wrote', len(src), 'bytes')
