# VOAM

The VoltOhmAmpereMaxwell.f95 (VOAM.f95) FORTRAN code is based on the Field Dependent Barrier (FDB) method presented by Pau Besalú-Sala _et al._ 

$$ \Delta E^\ddagger (**F**) = \Delta E^\ddagger (**F**=0) - \sum_{i=1}^{x,y,z} \Delta \mu_i F_i - \frac{1}{2} \sum_{i,j}^{x,y,z} \Delta \alpha_{ij} F_i F_j - \frac{1}{6} \sum_{i,j,k}^{x,y,z} \Delta \beta_{ijk} F_i F_j F_k \dots $$ 

together with the analytic solution of the cubic equation published in the 16th Century by Cardano inspired both by del Ferro and Tartaglia's work

$$ ax^3 + bx^2 + cx + d = 0 $$

$$ x= \Gamma^n \cdot \sqrt[3]{-\frac{q}{2} \sqrt{\frac{q^2}{4}+\frac{p^3}{27}}} + \Gamma^{n+3} \cdot \sqrt[3]{-\frac{q}{2} \sqrt{\frac{q^2}{4}-\frac{p^3}{27}}} - \frac{b}{3a} \quad \text{for} \quad n=0,1,2$$

where $\Gamma$ is the primitive cube root of unity defined as

$$ \Gamma^1 = -\frac{1}{2} + i\frac{\sqrt{3}}{2} \quad \Gamma^2 = \Gamma^* = -\frac{1}{2} - i\frac{\sqrt{3}}{2} \quad \Gamma^3 = 1$$

and the p and q coefficients of the general solutions read as

$$ p = \frac{3ac-b^2}{3a^2} \quad q=\frac{2b^3-9abc+27a^2d}{27a^3} $$

and are the consequence of making a change of variable $ z = x-\frac{b}{3a} $ such as the second order of the cubic equation vanishes.
