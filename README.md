# VOAM

The VoltOhmAmpereMaxwell.f95 (VOAM.f95) FORTRAN code is based on the Field Dependent Barrier (FDB) method presented by Pau Besalú-Sala _et al._ 

$$ \Delta E^\ddagger (**F**) = \Delta E^\ddagger (**F**=0) - \sum_{i=1}^{x,y,z} \Delta \mu_i F_i - \frac{1}{2} \sum_{i,j}^{x,y,z} \Delta \alpha_{ij} F_i F_j - \frac{1}{6} \sum_{i,j,k}^{x,y,z} \Delta \beta_{ijk} F_i F_j F_k \dots $$ 

together with the analytic solution of the cubic equation published in the 16th Century by Cardano inspired both by del Ferro and Tartaglia's work

$$ ax^3 + bx^2 + cx + d = 0 $$

$$ x= \Gamma^n \cdot \sqrt[3]{-\frac{q}{2} \sqrt{\frac{q^2}{4}+\frac{p^3}{27}}} + \Gamma^{n+3} \cdot \sqrt[3]{-\frac{q}{2} \sqrt{\frac{q^2}{4}-\frac{p^3}{27}}} \quad \text{for}\; n=0,1,2$$



