# PHYSICELL RETINAL DEVELOPMENT AND LAMINTION MODEL
This PhysiCell model seeks to recapitulate the proliferative and migratory processes in the developing human retina that give rise to characteristic adult lamina. The model is built in the PhysiCell ABM framework:

A. Ghaffarizadeh, R. Heiland, S.H. Friedman, S.M. Mumenthaler, and P. Macklin, PhysiCell: an open source physics-based cell simulator for 3-D multicellular systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: 10.1371/journal.pcbi.1005991.

# MODEL OVERVIEW
The model introduces an initial RPC cell population that evolves over time into adult tissue by neurogenic production of Retinal Ganglion, Amacrine, Müller Glial, Bipolar, Horizontal, Rod and Cone cells chosen via intrinsic developmental competence and proliferative cell divisions across all resulting cell types.

# RPC competence 
Given an atlas $A$ of cell type fates, the characteristic competence of the $i$th RPC in-simulation is realized as a discrete pdf $\left(\lambda_{\alpha\in A}^{i}\right)$ of the RPC toward individual cell fates $\alpha\in A$; i.e., $$0\leq \lambda_{\alpha}^{i}\leq 1,\quad \forall \alpha\in A$$ and $$\sum_{\alpha\in A}\lambda_{\alpha}^{i}=1,$$ such that $\lambda_{\alpha'}^{i}$ gives the probability of the $i$th RPC producing cell type $\alpha'$ during a neurogenic cell division.

The individual probabilities change according to 7 catalogued TIFs $\tau\in T=\{\}$ that evolve determinsticallyas a sum $\sum_{\i=1}^{n_{t}}{G^{t}_{i}}$ of 7-dimensional Gaussian curves with Expectation and Variance vectors $\{\mu_{t}\}$ $\{(\sigma^{2})_{t}\}$. 

Prolonged exposure to sufficiently homogeneous cell neighborhoods also change RPC cell competency according to the following ODE:

$$\frac{d\lambda_{\alpha}^{i}}{dt} = \sum_{\} + \[\]$$


 
