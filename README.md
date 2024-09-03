# VOAM

The VoltOhmAmpereMaxwell.f95 (VOAM.f95) FORTRAN code is based on the Field Dependent Barrier (FDB) method presented by Pau Besalú-Sala _et al._ 

$$ \Delta E^\ddagger (**F**) = \Delta E^\ddagger (**F**=0) - \sum_{i=1}^{x,y,z} \Delta \mu_i F_i $$ - \frac{1}{2} \sum_{i,j}^{x,y,z} \Delta \alpha_{ij} F_i F_j - \frac{1}{6} \sum_{i,j,k}^{x,y,z} \Delta \beta_{ijk} F_i F_j F_k 

together with the analytic solution of the cubic equation published in the 16th Century by Cardano inspired both by del Ferro and Tartaglia's work. 
$$

