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

$$ \Delta E^\ddagger (**F**) = \Delta E^\ddagger (**F**=0) - \Delta \mu^\ddagger |F| - \frac{1}{2} \Delta \alpha^\ddagger F^2 - \frac{1}{6} \Delta \beta^\ddagger |F^3| $$

In consquence, the total expression can be gathered as

$$ \Delta E^\dagger (**F**) = \Delta E^\ddagger_{main} (**F**=0) - \Delta \mu^\ddagger_{main} F - \frac{1}{2} \Delta \alpha^\ddagger_{main} F^2 - \frac{1}{6} \Delta \beta^\ddagger_{main} F^3 + \Delta E^\ddagger_{smallmolecules} (**F**=0) - \Delta \mu^\ddagger_{smallmolecules} |F| - \frac{1}{2} \Delta \alpha^\ddagger_{smallmolecules} F^2 - \frac{1}{6} \beta^\ddagger_{smallmolecules} |F^3| $$

or simply as

$$ \Delta E^\dagger (**F**) = \Delta E^\ddagger (**F**=0) - (\Delta \mu^\ddagger_{main} F + \Delta \mu^\ddagger_{smallmolecules} |F|) - \frac{1}{2} \Delta \alpha F^2 - \frac{1}{6} (\Delta \beta^\ddagger_{main} F^3 + \beta^\ddagger_{smallmolecules} |F^3|) $$

This equation is hard to solve analytically and there are no simple approximations to solve this equation, which forces us to separate the expression for negative and positive electric fields:

$$ \Delta E^\dagger (**F^-**) = \Delta E^\ddagger (**F**=0) - (\Delta \mu^\ddagger_{main} + \Delta \mu^\ddagger_{smallmolecules})F^{(-),1} - \frac{1}{2} \Delta \alpha F^2 - \frac{1}{6} (\Delta \beta^\ddagger_{main} + \beta^\ddagger_{smallmolecules}) F^{(-),3} $$

$$ \Delta E^\dagger (**F^+**) = \Delta E^\ddagger (**F**=0) - (\Delta \mu^\ddagger_{main} + \Delta \mu^\ddagger_{smallmolecules})F^{(+),1} - \frac{1}{2} \Delta \alpha F^2 - \frac{1}{6} (\Delta \beta^\ddagger_{main} + \beta^\ddagger_{smallmolecules}) F^{(+),3} $$

So the quick scan are performed according to the just presented equations. Then, the cubic equation to be solved is tunned according to where the input target barrier is closest (and lowest field strength) to _i.e._ in the negative or positive field domain of fields.

When trying to compute 2- and 3-D this approach is still applied but using the very modulus of the electric field instead of a single component of the electric field.

The input file necessary for the execution of VOAM, though it may seem complex, contains enough information to compute the extra necessary data for the computation of the optimal electric field. It is separated in two sections: the keyword and the NLOPs section.

The keyword section is then split into 4 lines of keywords: Method, Thermochemistry, Computation and Extra.
- **Method**: in this line we only define which kind of scan wants to be performed (X, YZ or XYZ) and the level of approximation of the FDB method expansion.
  - Axis scan: here we define the dimensionality of the scan to be performed. The following set of keywords will trigger the different resolutions
    - 0D: Scan, scan, SCAN
    - 1D: Can analitically compute the solution for the X, Y and Z axis.
      
        X axis: ("00X","0X0","X00","0x0","x00","00x","X","x","XXX","xxx","xx","XX")
      
        Y axis: ("00Y","0Y0","Y00","0y0","y00","00y","Y","y","YYY","yyy","yy","YY")
      
        Z axis: ("00Z","0Z0","Z00","0z0","z00","00z","Z","z","ZZZ","zzz","zz","ZZ")
      
    - 2D: Numerically compute the solution in the XY, YZ and XZ planes.
      
        XY plane: ("XY0","xy0","xy","XY","0XY","0xy","x0y","X0Y")
      
        YZ plane: ("0YZ","0yz","YZ0","yz0","ZY0","zy0","YZ","yz","y0z","Y0Z")
      
        XZ plane: ("X0Z","x0z","0xz","xz0","0XZ","XZ0","xz","XZ")
      
    - 3D: Compute the optimal strength for a three-dimensiona electric field vector. It is activated with "XYZ" or "xyz"
      
  - Taylor: level approximation of the FDB method expansion.
    
    - First order: ("Dipole","dipole","mu","Mu")
      
    - Second order: ("Alpha","alpha")
      
    - Third order: ("Beta","beta")
      
  Any value different than the proposed ones will automatically stoped the execution of the code.
  
- **Thermochemistry**: Line regarding the computation of the $\Delta$G of the input reaction.
  
  - Target barrier: It allows float values with up to two decimal points.
    
  - Redox: Number of electrons that are being exchanged in the reaction. Only integers.
    
      - Input negative values are defined as reductions _i.e._ electrons are in the reactants reducing the main species.
        
      - Input positive values are defined as oxidations _i.e._ electrons leave the system oxidising the main species.
        
  - Potential: Redox potential (in V) vs SHE. It then can be converted to the ferrocene (0/+) reference electrode. More electrodes are yet to be implemented.
    
- **Computation**:
  
  - Initial point: Defines the constant electric field to be considered for the computation of the energy and the NLOPs. VOAM.f95 can also detect the possible incompatibilities _i.e._ performing a 3D with constant electric fields.
    
  - Modulus: Defines the maximum modulus to be scanned and be allowed as a solution. Recommended value of 100 ($\cdot 10^{-4}$ a.u.)
    
  - Grid: The number of partitions to be applied in the $\theta$ angle of the spherical polar coordinates of the $\phi$ in the cylindrical coordinate system. Recommended value of 1000 or 2000.
  
- **Extra**: at the moment is only applied for the numerical tolerance. The default value is 1E-8.
  
  - Tolerance: tolerance

 After the keyword section, as previously mentioned, it is followed the NLOP section. The different chemical species participating in the chemical reaction are separated by delimiters such "-#-#-#- Chemical species/Small molecules: Reactant/Product -#-#-#-". It is recommended to look into the "Examples" folder, which contains 4 different examples, to have a better grasp of the structure of the input file.
