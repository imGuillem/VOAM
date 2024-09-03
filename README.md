# VOAM

The VoltOhmAmpereMaxwell.f95 (VOAM.f95) FORTRAN code is based on the Field Dependent Barrier (FDB) method presented by Pau Besalú-Sala _et al._ 

$$ \Delta E^\ddagger (**F**) = \Delta E^\ddagger (**F**=0) - \sum_{i=1}^{x,y,z} \Delta \mu_i F_i - \frac{1}{2} \sum_{i,j}^{x,y,z} \Delta \alpha_{ij} F_i F_j - \frac{1}{6} \sum_{i,j,k}^{x,y,z} \Delta \beta_{ijk} F_i F_j F_k \dots $$ 

where μ, α and β are the expansion coefficients defined as the electric dipole moment, the (static) polarisability matrix and the (static) hyperpolarisability first order tensor. This equation is combined with the analytic solution of the cubic equation published in the 16th Century by Cardano inspired both by del Ferro and Tartaglia's work

$$ ax^3 + bx^2 + cx + d = 0 $$

$$ x= \Gamma^n \cdot \sqrt[3]{-\frac{q}{2} \sqrt{\frac{q^2}{4}+\frac{p^3}{27}}} + \Gamma^{n+3} \cdot \sqrt[3]{-\frac{q}{2} \sqrt{\frac{q^2}{4}-\frac{p^3}{27}}} - \frac{b}{3a} \quad \text{for} \quad n=0,1,2$$

where $\Gamma$ is the primitive cube root of unity defined as

$$ \Gamma^1 = -\frac{1}{2} + i\frac{\sqrt{3}}{2} \quad \Gamma^2 = \Gamma^* = -\frac{1}{2} - i\frac{\sqrt{3}}{2} \quad \Gamma^3 = 1$$

and the p and q coefficients of the general solutions read as

$$ p = \frac{3ac-b^2}{3a^2} \quad q=\frac{2b^3-9abc+27a^2d}{27a^3} $$

and are the consequence of making a change of variable z = x-b/3a such as the second order of the cubic equation vanishes and the coefficient of the third order term is 1.

In addition, VOAM does also consider the reorientation of small molecules along their dipole moment when an electric field is applied. Hence, the methodology applied in VOAM.f95 begins with a quick scan of the FDB method equation while considering that the effect of the electric field in small molecules is independent of the orientation of the electric field, _i.e._ acting as an absolute value. Consequently, the FDB equation is going to be different according to the orientation (sign) of the electric field:

$$ \Delta E^\ddagger_{reaction} (**F**) = \Delta E^\ddagger_{main} (**F**) + \Delta E^\ddagger_{small molecules} $$

where the label "main" refers to those chemical species that are not easily reoriented with the electric field, opposed to the small molecules. In this regard, the "main" Maclaurin series is going to have the same shape as the first equation (_vide supra_) but the "small molecules" expansion in one dimension reads as

$$ \Delta E^\ddagger (**F**) = \Delta E^\ddagger (**F**=0) - \Delta \mu |F| - \frac{1}{2} \Delta \alpha |F^2| - \frac{1}{6} \Delta \beta |F^3| $$

and knowing that the $|x^2|$ is still $x^2$ it can be alternatively be written as

$$ \Delta E^\ddagger (**F**) = \Delta E^\ddagger (**F**=0) - \Delta \mu \abs(F) - \frac{1}{2} \alpha F^2 - \frac{1}{6} \beta |F^3| $$

In consquence, the total expression can be gathered as

$$ \Delta E^\dagger (**F**) = [ \Delta E^\ddagger (**F**=0) - \Delta \mu F - \frac{1}{2} \alpha F^2 - \frac{1}{6} \beta F^3 ] $$

