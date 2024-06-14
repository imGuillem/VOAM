!-- Program written by Pau Besalú and Guillem Pey* - TFM Guillem Pey
!-- VOAM.f95 (VoltOhmAmpereMaxwell.f95) 
!-- OpEEF searcher within the FDB_beta method using the third degree equation and polar spherical coordinates
Module VOAM
implicit none
! Non-linear optical properties (NLOPs) and energies
    double precision :: E_0,E_r,E_p,E_f,target_barrier
    double precision :: E_mu,E_alpha
    integer, dimension (3) :: axis = 0
    double precision, dimension(3) :: F,mu,mu_r,mu_p,mu_f
    double precision, dimension(3,3) :: alpha,alpha_r,alpha_p,alpha_f,alpha_tmp
    double precision, dimension(3,3,3) :: beta,beta_r,beta_p,beta_f
! VOAM keywords
    double precision :: redox_potential,c_i,radius,max_radius,tol
    integer order_nlop,n_dim,x_axis,y_axis,z_axis
    integer stoi,redox_coef,grid 
    integer gpc
! Guillem keywords
    logical, dimension(10) :: Guillem = .FALSE.
    logical :: sudo = .FALSE.
    logical :: print_NLOP_react,print_NLOP_prod,print_NLOP_lig ! Guillem(1-3)
    logical :: print_coeffs = .FALSE. ! Guillem(4)
    logical :: print_EQ = .FALSE. ! Guillem(5)
    logical :: print_sols = .FALSE. ! Guillem(6)
    logical :: print_checks = .FALSE. ! Guillem(7)
! Others
    character*2,dimension(3) :: axis_name
    double precision :: PI = 4.0d0*datan(1.0d0)
    double precision :: eps = cmplx(-0.5d0,dsqrt(3.0d0/4.0d0))
    integer i,j,k
contains
    Function E_FDB(F,order_nlop)
    implicit none
    double precision :: E_FDB,Fx,Fy,Fz,theta,phi
    !double precision :: E_0,E_mu,E_alpha,redox_potential
    double precision, dimension(3) :: mu,mu_f,F
    double precision, dimension(3,3) :: alpha,alpha_f
    double precision, dimension(3,3,3) :: beta,beta_f
    integer i,j,k,order_nlop,redox_coef,stoi

    E_FDB=E_0-abs(redox_coef)*redox_potential
    !-- Just in case
    !if (mu_f(1).lt.0) then
    !    write(*,*) "WARNING! LIGAND NOT PROPERL ALIGNED! STOP!"
    !    stop
    !end if
    !-- Linear term
    do i=1,3
        E_FDB=E_FDB-mu(i)*F(i)
    end do
        !-- Adding solvent correction
    E_FDB=E_FDB+stoi*mu_f(1)*sqrt(F(1)**2+F(2)**2+F(3)**2)
    E_mu=E_FDB

    !-- Quadratic term
    do i=1,3
        do j=1,3
            E_FDB=E_FDB-0.5d0*alpha_f(i,j)*F(i)*F(j)
        end do
    end do
        !-- Adding solvent correction
    E_FDB=E_FDB+stoi*0.5d0*alpha_f(1,1)*(sqrt(F(1)**2+F(2)**2+F(3)**2))**2
    E_alpha=E_FDB

    !-- Cubic term
    do i=1,3
        do j=1,3
            do k=1,3
                E_FDB=E_FDB+(1.0d0/6.0d0)*beta(i,j,k)*F(i)*F(j)*F(k)
            end do
        end do
    end do
        !-- Adding solvent correction
    E_FDB=E_FDB-stoi*(1.0d0/6.0d0)*beta(1,1,1)*(sqrt(F(1)**2+F(2)**2+F(3)**2))**3
    !-- Choosing the approximation of the run
    if (order_nlop.eq.0) then
        write(*,*) "WARNING: You must choose an approximation to optimise the system. STOP!"
        stop
    else if (order_nlop.eq.1) then
        E_FDB=E_mu
    else if (order_nlop.eq.2) then
        E_FDB=E_alpha
    else if (order_nlop.eq.3) then
        E_FDB=E_FDB
    else
        write(*,*) "WARNING: NOT A VALID APPROXIMATION KEYWORD. STOP"
        stop
    end if
    return
    End

!###################################################################################################
    !=========================================================================!
    !============================== 1D REGIME ================================!
    !=========================================================================!
!###################################################################################################

    Subroutine VOAM_1D(Guillem,F,gpc,tol,radius,axis,axis_name,order_nlop)
    implicit none
    logical,dimension(10) :: Guillem
    double precision, dimension (3) :: F
    character*2,dimension(3) :: axis_name
    integer, dimension(3) :: axis
    double precision :: E_FF,tmp_E,tmp_mu,tmp_alpha,radius
    integer step,gpc,order_nlop
        !---Specific parameters for the 3rd degree solver---!
    integer i,n
    double complex, dimension(3) :: root

    double complex :: eenergy,mmu,aalpha,bbeta
    double complex :: p,q,mark,root1,root2,root3

    double precision :: tol,discriminant
    double precision :: root_1,root_2,root_3,half_root1,half_root2,min_root
    double precision :: omega,angle,pp,qq,markk,markkk

    !--         Reference of the 3rd degree solution:
    !-- http://olewitthansen.dk/Mathematics/The_cubic_equation.pdf
    E_FF=E_0-abs(redox_coef)*redox_potential
    tmp_E=E_FF;tmp_E=tmp_E-target_barrier/6.2751d2                  !-- Energy 
    mu(gpc)=mu(gpc)+stoi*mu_f(gpc)                                  !-- Dipole moment
    alpha(gpc,gpc)=alpha(gpc,gpc)+stoi*alpha_f(gpc,gpc)             !-- Polarizability matrix
    beta(gpc,gpc,gpc)=beta(gpc,gpc,gpc)+stoi*beta_f(gpc,gpc,gpc)    !-- Hyperpolarizability tensor
    tmp_mu=mu(gpc);tmp_alpha=0.5d0*alpha(gpc,gpc)                   !-- Tmp values for variations in case of constant fields

                !--                 Step of the do loop according to the gpc value:
                !-- Generated through as a two-valued function that can only yield values of 1 or 2
    step=nint(1+abs(cos(0.5d0*PI*gpc)))    
                            !-- Value of the dummy indices for the nested loops
            ! X_axis(step=1) = 2,3 // 2,3; Y_axis(step=2) = 1,3 // 1,3; Z_axis(step=1) = 1,2 // 1,2
            ! X_axis--> gpc=1            ; Y_axis--> gpc=2            ; Z_axis--> gpc=3
            ! X_axis--> axis(1)=1        ; Y-axis--> axis(2)=2        ; Z_axis--> axis(3)=3

                !-- "General" formula for the computation of the parameters of the equations within the FDB beta method
                !-- Beta tensor components are changed in sign because of G16  
    do i=(1+axis(1)),(3-axis(3)),step
        do j=(1+axis(1)),(3-axis(3)),step
            tmp_E=tmp_E-mu(i)*radius*F(i)-0.5d0*alpha(i,j)*(radius**2)*F(i)*F(j)+(1.0d0/6.0d0)*(radius**3)*beta(i,i,i)*F(i)**3+0.5d0*(radius**3)*beta(i,i,j)*F(j)*F(i)**2 
    
            tmp_mu=tmp_mu+alpha(gpc,i)*radius*F(i)-0.5d0*beta(gpc,i,j)*radius*F(i)*F(j)

            tmp_alpha=tmp_alpha-0.5d0*beta(gpc,gpc,i)*radius*F(i) 
        end do
    end do

        !-- Declaration of each parameter as a complex number for its proper handling
    eenergy=cmplx(tmp_E,0.0d0)
    mmu=-cmplx(tmp_mu,0.0d0) !-- Change of sign because of the series expansion definition
    aalpha=-cmplx(tmp_alpha,0.0d0) !-- Change of sign because of the series expansion definition
    bbeta=(1.0d0/6.0d0)*cmplx(beta(gpc,gpc,gpc),0.0d0) !-- Gaussian already changes the sign of the beta tensor, and then later changed by the FDB expansion.
                                                       !   So, it is a positive value. 

    write(*,'(" Target barrier to be solved:",xF6.2," kcal/mol")') target_barrier
    if(Guillem(4).eqv..TRUE.) then
        write(*,*)
        write(*,'(" Beta coefficient:           ",xF18.6)') real(bbeta)
        write(*,'(" Alpha coefficient:          ",xF18.6)') real(aalpha)
        write(*,'(" Dipole moment coefficient:  ",xF18.6)') real(mmu)
        write(*,'(" Energy coefficient:         ",xF18.6)') real(eenergy)
        write(*,*)
    end if

        !-- Considering a change of variable to remove the second degree term of the third order equation in order to simplify 
    p=(3.0d0*bbeta*mmu-aalpha**2.0d0)/3.0d0/bbeta**2.0d0
    q=(2.0d0*aalpha**3.0d0-9.0d0*bbeta*aalpha*mmu+27.0d0*eenergy*bbeta**2.0d0)/27.0d0/bbeta**3.0d0
    omega=real(0.25d0*q**2.0d0+(1.0d0/27.0d0)*p**3.0d0) 

        !-- The impossible case --> ΔΔG=0
    if (abs(bbeta).lt.abs(tol).and.abs(aalpha).lt.abs(tol).and.abs(mmu).lt.abs(tol)) then
        write(*,*) "Input error. Target barrier is the same as the desired thermochemistry"
        write(*,*) "Either consider changing the tolerance or another thermochemistry value"
        write(*,*)
        stop
    end if

!===============================================================================
!###############################################################################
!===============================================================================

        !===== μ-FDB-1D solution ------ General solution for 1st degree equation
    if (order_nlop.eq.1.or.(order_nlop.gt.1.and.abs(aalpha).lt.real(tol).and.abs(bbeta).lt.real(tol)**2)) then
        if(Guillem(5).eqv..TRUE.) then ! Print the equation to be solved
            write(*,'(" (μ-FDB-",A2,") - Equation to solve":,xF12.6"x^1 +",xF12.6,"x^0 = 0")') axis_name(gpc),real(mmu),real(eenergy)
        end if
        if (order_nlop.gt.1.and.abs(aalpha).lt.real(tol).and.abs(bbeta).lt.real(tol)**2) then
            write(*,*) "---------------------------------------------"
            write(*,*) "    Alternative first order solution"
            write(*,*) "---------------------------------------------"
        end if
        root1=-real(eenergy)/real(mmu)
        if (real(root1).gt.max_radius) then
            write(*,*) "The solution is outside of the proposed radius"
        else
            write(*,'(" The minimum field for the FDB-1D,μ approximation is ",A2,"=",xF18.2,"(*10^-4 a.u)")') axis_name(gpc),real(root1)*1.0d4
        end if
        if (Guillem(6).eqv..TRUE.) then
            write(*,'(" (μ-FDB-",A2,") - Specific solution:",xF20.6," a.u")') axis_name(gpc),real(root1)
            write(*,*)
        end if
    return
    end if  
        !===== μ-FDB-1D solution -----------------------------------------------

!===============================================================================
!###############################################################################
!===============================================================================

        !===== α-FDB-1D solution(s) --- General solution for 2nd degree equation
    if (order_nlop.eq.2.or.(order_nlop.gt.2.and.abs(bbeta).lt.real(tol)**2)) then

            !-- Print the equation to be solved
        if(Guillem(5).eqv..TRUE.) then 
            write(*,'(" (α-FDB-",A2,") - Equation to solve:",xF12.6,"x^2 +",xF12.6,"x^1 +",xF12.6,"x^0 = 0")') axis_name(gpc),real(aalpha),real(mmu),real(eenergy)
            write(*,*)
        end if

            !-- Consideration of two alternative cases
        if (order_nlop.gt.2.and.abs(bbeta).lt.real(tol)**2) then !-- 1st hyperpolarizabilty is neglegible
            write(*,*) "--------------------------------------------------------------------------------"
            write(*,*) "     Alternative second order solution - Neglegible 1st hyperpolarizability     "
            write(*,*) "--------------------------------------------------------------------------------"
        end if
        if (abs(bbeta).lt.real(tol)**2.and.abs(aalpha).gt.abs(mmu)*1.0d5) then ! Dipole moment is neglegible
            write(*,*) "--------------------------------------------------------------------------------"
            write(*,*) "         Alternative second order solution - Neglegible dipole moment           "
            write(*,*) "--------------------------------------------------------------------------------"
            root1=+sqrt(-eenergy/aalpha);root2=-sqrt(-eenergy/aalpha)
            if (abs(root1)*1.0d4.gt.max_radius.or.abs(root2)*1.0d4.gt.max_radius) then ! Out of the proposed radius
                write(*,'("The solutions are outside of the proposed radius (",xF20.2," +i",xF20.2," a.u)")') root1*1.0d4
            else
                write(*,'(" (α-FDB-",A2,") - First solution:  ",xF15.6," +i",xF15.6," a.u")') axis_name(gpc),root1*1.0d4
                write(*,'(" (α-FDB-",A2,") - Second solution: ",xF15.6," +i",xF15.6," a.u")') axis_name(gpc),root2*1.0d4
                write(*,*)
                return
            end if
        end if

                !-- Regular 2nd degree equation solution --!
        discriminant=mmu**2.0d0-4.0d0*aalpha*eenergy
        if(Guillem(6).eqv..TRUE.) then ! Specific solutions do not have the 1.0d4 factor
            write(*,'(" Discriminant of the polynomial:",xF10.4)') discriminant
            write(*,'(" (α-FDB-",A2,") - Specific solution 1:",xF15.6," +i*",xF15.6," a.u")') axis_name(gpc),(-mmu+sqrt(discriminant))/2.0d0/aalpha
            write(*,'(" (α-FDB-",A2,") - Specific solution 2:",xF15.6," +i*",xF15.6," a.u")') axis_name(gpc),(-mmu-sqrt(discriminant))/2.0d0/aalpha
            write(*,*)
        end if
        if (discriminant.ge.0.0d0) then
            root1=(-mmu+sqrt(discriminant))/2.0d0/aalpha
            root2=(-mmu-sqrt(discriminant))/2.0d0/aalpha
        else
            write(*,*) "There does not exist a *REAL* solution for the proposed system"
        end if

            !-- Get the closest-to-zero solution 
        min_root=root1
        if (abs(root2).lt.abs(min_root)) min_root=root2

            !-- Check whether the solutions are inside the resolution sphere or not
        if(abs(root1)*1.0d4.gt.max_radius.and.abs(root2)*1.0d4.gt.max_radius) then
            write(*,'(" The solution is outside of the proposed radius (",xF18.6,"*10^-4 a.u)")') min_root*1.0d4
            write(*,*)
        else 
            write(*,*) "The proposed solutions for the system are:" 
            write(*,'(xA2,"=",xF10.2,"(*10^-4) a.u and",xA2,"=",xF10.2,"(*10^-4) a.u")') axis_name(gpc),real(root1)*1.0d4,axis_name(gpc),real(root2)*1.0d4
            write(*,*)
            write(*,'(" (α-FDB-",A2,") - The minimum field strength for the proposed system is:",xA2,"=",xF12.2,"(*10^-4) a.u")') axis_name(gpc),axis_name(gpc),min_root*1.0d4
            write(*,*)
        end if
    
            !-- Print the checks of the solutions
        if (Guillem(7).eqv..TRUE.) then
            write(*,'(" (α-FDB-",A2,") - First root: ",xF18.6 )') axis_name(gpc),root1
            write(*,*) abs(aalpha*root1**2+mmu*root1+eenergy)
            write(*,'(" (α-FDB-",A2,") - Second root:",xF18.6 )') axis_name(gpc),root2
            write(*,*) abs(aalpha*root2**2+mmu*root2+eenergy)
        end if
    end if
        !===== α-FDB-1D solution(s) --------------------------------------------

!===============================================================================
!###############################################################################
!===============================================================================

        !===== β-FDB-1D solution(s) --- General solution for 3rd degree equation 
    if(order_nlop.eq.3) then
            !-- Print the specific parameters for the (depressed) cubic equation
        if(Guillem(4).eqv..TRUE.) then
            write(*,*) "Omega parameter:    ",omega 
            write(*,*) "P parameter:        ",real(p)
            write(*,*) "Q parameter:        ",real(q)
            write(*,*)
        end if
            
            !-- Print the equation that is being solved
        if(Guillem(5).eqv..TRUE.) then
            write(*,'(" (β-FDB-",A2,") - Equation to solve:",xF12.6,"x^3 +",xF12.6,"x^2 +",xF12.6,"x^1 +",xF12.6,"x^0 = 0")') axis_name(gpc),real(bbeta),real(aalpha),real(mmu),real(eenergy)
            write(*,*)
        end if

            !-- Three different scenarios
        !================== CASE 1 CASE 1 CASE 1 CASE 1 CASE 1 CASE 1 CASE ===============!
        if(abs(omega).lt.real(tol)) then ! We have three real roots, two of them being equal // omega.eq.0.0d0
            write(*,*) "            Case 1 - Three real solutions, two of them equal" 
            write(*,*) "--------------------------------------------------------------------------------------------------------------------  "
            root_1=+2.0d0*sign(abs(-0.5d0*real(q))**(1.0d0/3.0d0),-0.5d0*real(q))-aalpha/3.0d0/bbeta
            root_2=-1.0d0*sign(abs(-0.5d0*real(q))**(1.0d0/3.0d0),-0.5d0*real(q))-aalpha/3.0d0/bbeta
            root_3=root_2
            root(1)=root_1;root(2)=root_2;root(3)=root_3

                !-- Print the raw specific roots
            if(Guillem(6).eqv..TRUE.) then 
                write(*,*) "        Specific solutions (a.u):"
                write(*,'(" (β-FDB-",A2,") - Case 1 - First solution: ",xF18.6)') axis_name(gpc),root_1
                write(*,'(" (β-FDB-",A2,") - Case 1 - Second solution:",xF18.6)') axis_name(gpc),root_2
                write(*,'(" (β-FDB-",A2,") - Case 1 - Third solution: ",xF18.6)') axis_name(gpc),root_3
                write(*,*)
            end if

                !-- Obtaining the closest-to-zero solution
            min_root=real(root(1))
            do i=2,3
                if (abs(root(i)).lt.abs(min_root)) min_root=real(root(i))
            end do

                !-- Check whether the solutions are inside the resolution sphere or not
            if (min_root*1.0d4.gt.max_radius) then
                write(*,*) "The minimum solution is out of the proposed radius"
                write(*,*) "Please consider either a bigger radius or an other thermochemistry value"
                write(*,*) 
            else
                write(*,'(" (β-FDB-",A2,") - Case 1 - The minimum field strength is:",xF12.2,"(*10^-4 a.u)")') axis_name(gpc),min_root*1.0d4
                write(*,*)
            end if

            if (Guillem(7).eqv..TRUE.) then    
                do i=1,3
                    write(*,'(" Root",xI1,":",xF18.6," +i",xF18.6)') i,root(i)
                    check(i)=bbeta*real(root(i))**3.0d0+aalpha*real(root(i))**2.0d0+mmu*real(root(i))+eenergy
                    write(*,'("      Check:",xF10.8)') abs(check(i))
                end do
            end if
            write(*,*) "--------------------------------------------------------------------------------------------------------------------  "
            return

            !================== CASE 2 CASE 2 CASE 2 CASE 2 CASE 2 CASE 2 CASE ===============!
        else if (omega.gt.real(tol)) then ! We have one real root and two complex conjugated roots. 
            write(*,*) "            Case 2 - One real root and two complex conjugated ones"
            write(*,*) "--------------------------------------------------------------------------------------------------------------------  "
            do i=1,3
                n=i-1 
                    !-- Generating the first([u_0]) and ([v_0])second partial solution ([u_0])
                half_root1=sign(abs(-0.5d0*q-sqrt(0.25d0*q**2+p**3/27.0d0))**(1.0d0/3.0d0),real(-0.5d0*q-sqrt(0.25d0*q**2+p**3/27.0d0)))
                half_root2=sign(abs(-0.5d0*q+sqrt(0.25d0*q**2+p**3/27.0d0))**(1.0d0/3.0d0),real(-0.5d0*q+sqrt(0.25d0*q**2+p**3/27.0d0)))

                    !-- Applying the unitary complex root for the answers and undoing the variable change to recover the original answer
                root(i)=half_root1*(eps**n)+half_root2*(eps**(2*n))-aalpha/3.0d0/bbeta
            end do

                !-- Check whether the solutions are inside the resolution sphere or not
            if (abs(root(1)).gt.max_radius) then ! The solution is out of the scan
                write(*,*) " The real root is out of the predefined scanning radius"
                write(*,*)
            else 
                write(*,'(" (β-FDB-",A2,") - Case 2 - The minimum field strength is:",xF12.2,"(*10^-4 a.u)")') axis_name(gpc),real(root(1))*1.0d4
                write(*,*)
            end if
                
                !-- Print the raw specific roots
            if (Guillem(6).eqv..TRUE.) then
                write(*,'(" (β-FDB-",A2,") - Case 2 - Unique specific real solution:   ",xF18.6," a.u")') axis_name(gpc),root(1)
                write(*,'(" (β-FDB-",A2,") - Case 2 - First specific complex solution: ",xF18.6," +i",xF18.6" a.u")') axis_name(gpc),root(2)
                write(*,'(" (β-FDB-",A2,") - Case 2 - Second specific complex solution:",xF18.6," +i",xF18.6" a.u")') axis_name(gpc),root(3)
                write(*,*)
            end if
        
                !-- Print the checks of the solutions
            if (Guillem(7).eqv..TRUE.) then
                do i=1,3
                    write(*,'(" Root",xI1,":",xF18.6," +i",xF18.6)') i,root(i)
                    check(i)=bbeta*real(root(i))**3.0d0+aalpha*real(root(i))**2.0d0+mmu*real(root(i))+eenergy
                    write(*,'("      Check:",xF10.8)') abs(check(i))
                end do
            end if 
            write(*,*) "--------------------------------------------------------------------------------------------------------------------  "
            return

            !================== CASE 3 CASE 3 CASE 3 CASE 3 CASE 3 CASE 3 CASE ===============!
        else if (omega.lt.0.0d0) then ! We have three real solutions
            write(*,*) "            Case 3 - Three real different solutions"
            write(*,*) "--------------------------------------------------------------------------------------------------------------------  "
                !-- Computation of the three real solutions --!
            angle=dacos(-0.5d0*real(q)/sqrt(-1.0d0*real(p)**3/27.0d0))

                        !-- Formulae for each root --!
            root_1=2*sqrt(-pp/3.0d0)*dcos(angle/3.0d0)-aalpha/3.0d0/bbeta
            root_2=2*sqrt(-pp/3.0d0)*dcos((angle+2*PI)/3.0d0)-aalpha/3.0d0/bbeta
            root_3=2*sqrt(-pp/3.0d0)*dcos((angle+4*PI)/3.0d0)-aalpha/3.0d0/bbeta
            root(1)=root_1;root(2)=root_2;root(3)=root_3

            if(Guillem(6).eqv..TRUE.) then !-- Print raw roots
                write(*,'(" The solutions given by VOAM_1D are:")')
                write(*,'(" (β-FDB-",A2,") - Case 3 - First specific solutions (a.u):  ",xF15.6)') axis_name(gpc),root_1
                write(*,'(" (β-FDB-",A2,") - Case 3 - Second specific solutions (a.u): ",xF15.6)') axis_name(gpc),root_2
                write(*,'(" (β-FDB-",A2,") - Case 3 - Third specific solutions (a.u):  ",xF15.6)') axis_name(gpc),root_3
                write(*,*)
            end if
    
                !-- Finding the closest-to-zero root --!
            min_root=real(root(1))
            do i=2,3
                if (abs(root(i)).lt.abs(min_root)) min_root=real(root(i))
            end do

                !-- Print of the minimum solution --!
            if (min_root*1.0d4.gt.max_radius) then
                write(*,'(" The minimum solution is out of the proposed radius (",F10.2," *10^-4 a.u)")') min_root*1.0d4 
            else
                write(*,'(" (β-FDB-",A2,") - Case 3 - The minimum field strength is:",xF20.6,"(*10^-4) a.u")') axis_name(gpc),min_root*1.0d4
            end if

                !-- Print the checks of the solutions
            if (Guillem(7).eqv..TRUE.) then 
                do i=1,3
                   write(*,'(" Root",xI1,":",xF18.6," +i",xF18.6)') i,root(i)
                   check(i)=bbeta*real(root(i))**3.0d0+aalpha*real(root(i))**2.0d0+mmu*real(root(i))+eenergy
                   write(*,'("      Check:",xF10.8)') abs(check(i))
                end do
            end if
            write(*,*) "--------------------------------------------------------------------------------------------------------------------  "
            return
        end if
    end if

!===============================================================================
!###############################################################################
!===============================================================================

    if (order_nlop.eq.4) then
        write(*,*) " Gamma approach is not implemented yet!"
        write(*,*) " Please consider lower approximations."
    end if
    End subroutine VOAM_1D

!###################################################################################################
    !=========================================================================!
    !============================== 2D REGIME ================================!
    !=========================================================================!
!###################################################################################################

    Subroutine VOAM_2D(Guillem,gpc,radius,axis,axis_name,order_nlop,grid) 
    implicit none
    logical, dimension(10) :: Guillem
    integer, dimension(3) :: axis
    character*2, dimension(3) :: axis_name
    double precision, dimension (3) :: F
    double complex, dimension (grid,2) :: tmp_root ! Stores only the modulus
    double precision, dimension (grid,3) :: root ! Modulus, theta and phi
    double precision :: radius
        !-- Only for FDB_μ
    double precision :: E_FF,numerator,denominator 
        !-- For FDB_α/β 
    double complex :: eenergy,mmu,aalpha,bbeta
    double complex :: p,q
    double precision :: omega
    integer :: sign_bbeta,sign_aalpha,sign_mmu,sign_eenergy
    double precision :: tmp_E,tmp_mu,tmp_alpha,tmp_beta
    double precision :: theta,phi,angle_step
    double precision :: discriminant,angle
        !-- Integers for iteration
    integer gpc,step,order_nlop,upper,lower
    integer grid,n,i,j

!================================ μ-VOAM-2D ==========================================    
    
    if (order_nlop.eq.1) then
        upper=nint(0.5d0*axis(1)+1.5d0*axis(2)+0.5d0*axis(3))
        lower=nint(0.5d0*axis(1)+0.5d0*axis(2)+2.5d0*axis(3))
            !-- axis_scan=XY --> (upper,lower)==(2,1);  axis_scan=XZ --> (upper,lower)==(1,3);  axis_scan=YZ --> (upper,lower)==(2,3)
        mu(upper)=mu(upper)+stoi*mu_f(3);mu(lower)=mu(lower)+stoi*mu_f(3)
        E_FF=E_0-abs(redox_coef)*redox_potential-(target_barrier/6.2751d2)-mu(gpc)*F(gpc)
        if (Guillem(5).eqv..TRUE.) then
            write(*,'(" (μ-FDB-(",A2,",",A2,")) Equation to solve: 0 = ",F12.6," +",xF12.6,xA2," +",xF12.6,xA2)') &
            & axis_name(upper),axis_name(lower),E_FF,mu(upper),axis_name(upper),mu(lower),axis_name(lower)
        end if
        numerator=mu(upper)*E_FF;denominator=mu(upper)**2+mu(lower)**2
        F(upper)=numerator/denominator; F(lower)=(mu(lower)/mu(upper))*F(upper)
        write(*,'(x"The minimum field strength is (",A2,",",A2,";",A2,")=(",F7.2,","F7.2,","F7.2,")(*10^-4) a.u)")') & 
        & axis_name(upper),axis_name(lower),axis_name(gpc),F(upper)*1.0d4,F(lower)*1.0d4,F(gpc)*1.0d4
    end if

!================================= μ-VOAM-2D ==========================================    
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!================================= α-VOAM-2D ==========================================    

    n=0
    if (order_nlop.eq.2) then
            ! XY=(1,1,0) gpc=3 // XZ=(1,0,1) gpc=2 // YZ=(0,1,1) gpc=1
        step=nint(1+abs(cos(0.5d0*PI*gpc)))
        tmp_E=E_0-abs(redox_coef)*redox_potential-target_barrier/6.2751d2     !-- Energy associated to thermochemistry
        tmp_mu=0.0d0;tmp_alpha=0.0d0
           !##############################################################
               !-- Scanning the XZ or YZ plane
           !##############################################################
        if (axis(1).eq.0.or.axis(2).eq.0) then
            do i=1,grid
                    ! Theta: 0.le.θ.le.π
                theta=((i-1)/(1.0d0*grid-1))*PI
                if(axis(1).eq.1) then
                    phi=0.0d0*PI;root(i,3)=phi !-- Scanning the XZ plane
                else if (axis(2).eq.1) then
                    phi=0.5d0*PI;root(i,3)=phi !-- Scanning the YZ plane
                end if

                    ! Field vector computation
                F(1)=dsin(theta)*dcos(phi);F(2)=dsin(theta)*dsin(phi);F(3)=dcos(theta)
                tmp_E=tmp_E-mu(gpc)*radius*F(gpc)-0.5d0*(radius**2)*alpha(gpc,gpc)*F(gpc)**2 !-- Energy associated to constant fields

                    ! Initialization of angular alpha calculation
                tmp_alpha=tmp_alpha+2.0d0*alpha(2-axis(1),2+axis(3))*F(2-axis(1))*F(2+axis(3))
                do j=(2-axis(1)),(2+axis(3)),step
                    tmp_mu=tmp_mu+F(j)*(mu(j)+alpha(j,gpc)*radius*F(gpc))
                    tmp_alpha=tmp_alpha+alpha(j,j)*F(j)**2
                end do
                eenergy=cmplx(tmp_E,0.0d0); mmu=-cmplx(tmp_mu,0.0d0)
                aalpha=-0.5d0*cmplx(tmp_alpha,0.0d0)
                if (Guillem(4).eqv..TRUE.) then
                    write(*,'(" Angles: Theta=",xF6.4," and Phi=",xF6.4)') theta,phi
                    write(*,'(" Energy coefficient:        ",xF18.6)') real(eenergy)
                    write(*,'(" Dipole moment coefficient: ",xF18.6)') real(mmu)
                    write(*,'(" Polarizability coefficient:",xF18.6)') real(aalpha)
                    write(*,*)
                end if

                if (abs(mmu).lt.tol) then ! Angular-mmu is neglegible --> ax**2+b=0
                    sign_aalpha=sign(1.0d0,real(aalpha));sign_eenergy=sign(1.0d0,real(eenergy))
                    if(sign_aalpha.eq.sign_eenergy) then ! The solution is complex. Skipped
                        if(Guillem(6).eqv..TRUE.) then ! Print the specific solutions
                            !-- Coded in this way so it is not stored in memory and is computed on-the-fly
                            write(*,'(" (α-FDB-",A2,A2," - Positive complex root:",xF12.6," +i",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),+1.0d0*sqrt(-eenergy/aalpha)
                            write(*,'(" (α-FDB-",A2,A2," - Negative complex root:",xF12.6," +i",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),-1.0d0*sqrt(-eenergy/aalpha)
                        end if
                        continue
                    else !-- Double signed root --> ax**2 -b=0
                            !-- It is defined a complex matrix to store both solutions. If they are they share the absolute value, we are going to apply a factor of
                            !   1/sqrt(2). If one of the solutions is deemed to be zero, we will conserve the module as is considered to be the lowest one
                        tmp_root(i,1)=+1.0d0*sqrt(-real(eenergy)/real(aalpha));tmp_root(i,2)=-1.0d0*sqrt(-real(eenergy)/real(aalpha))
                        root(i,3)=theta
                        if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                            write(*,'(" (α-FDB-",A2,A2,") - Neglegible dipole moment positive root:",xF18.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),real(tmp_root(i,1))
                            write(*,'(" (α-FDB-",A2,A2,") - Neglegible dipole moment negative root:",xF18.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),real(tmp_root(i,2))
                            write(*,*)
                        end if
                    end if
                else if (abs(aalpha).lt.tol) then ! Angular-aalpha is neglegible --> ax**1 +b=0
                    tmp_root(i,1)=-real(eenergy)/real(mmu);tmp_root(i,2)=0.0d0
                    root(i,3)=theta
                    if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                        write(*,'(" (α-FDB-",A2,A2,") - Neglegible alpha root:",xF12.6," a.u")') real(tmp_root(i,1))
                    end if
                else !-- Complete second degree solution
                    discriminant=mmu**2-4.0d0*aalpha*eenergy
                    if (discriminant.ge.0.0d0) then
                        root(i,3)=theta
                        tmp_root(i,1)=-mmu+sqrt(discriminant)/2.0d0/aalpha
                        tmp_root(i,2)=-mmu-sqrt(discriminant)/2.0d0/aalpha
                        if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                            write(*,'(" (α-FDB-",A2,A2,") - Positive pure second degree root:",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),real(tmp_root(i,1))
                            write(*,'(" (α-FDB-",A2,A2,") - Negative pure second degree root:",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),real(tmp_root(i,2))
                            write(*,*)
                        end if
                    end if
                end if
                if (abs(tmp_root(i,1)).eq.abs(tmp_root(i,2))) then ! The root is double signed
                    root(i,1)=abs(tmp_root(i,1)) !-- Absolute value because there is no need of negative modulus
                                                 !   the direction of the electric field can be tuned with the theta and phi
                end if
                    !-- Re-declare tmp_mu and tmp_alpha to restart the scan    
                tmp_mu=0.0d0;tmp_alpha=0.0d0
            end do
            !##################################################################
                !-- Scanning the XY plane
            !##################################################################
        else if (axis(3).eq.0) then
            do i=1,2*grid
                    !-- Phi: 0.le.φ.le.2π
                phi=((i-1)/(2.0d0*grid-1))*2.0d0*PI;theta=0.5d0*PI;root(i,2)=theta

                    !-- Declare F as a variable vector
                F(1)=dsin(theta)*dcos(phi);F(2)=dsin(theta)*dsin(phi);F(3)=dcos(theta)
                tmp_E=tmp_E-mu(gpc)*radius*F(gpc)-0.5d0*(radius**2)*alpha(gpc,gpc)*F(gpc)**2 !-- Energy associated to constant fields

                    !-- Initialization of angular alpha's calculation
                tmp_alpha=tmp_alpha+2.0d0*alpha(2-axis(1),2+axis(3))*F(2-axis(1))*F(2+axis(3))
                do j=(2-axis(1)),(2+axis(3)),step
                    tmp_mu=tmp_mu+F(j)*(mu(j)+alpha(j,gpc)*radius*F(gpc)) 
                    tmp_alpha=tmp_alpha+alpha(j,j)*F(j)**2
                end do
                
                    !-- Declaration of the parameters for the second degree solution
                eenergy=cmplx(tmp_E,0.0d0); mmu=-cmplx(tmp_mu,0.0d0)
                aalpha=-0.5d0*cmplx(tmp_alpha,0.0d0)

                    !-- Print the coefficients
                if (Guillem(4).eqv..TRUE.) then
                    write(*,'(" Angles: Theta=",xF6.4,"& Phi=",xF6.4)') theta,phi
                    write(*,'(" Energy coefficient:          ",xF18.6)') real(eenergy)
                    write(*,'(" Dipole moment coefficient:   ",xF18.6)') real(mmu)
                    write(*,'(" Polarizability coefficient:  ",xF18.6)') real(aalpha)
                    write(*,*)
                end if

                !###############################################################
                    !-- Computation of the solutions
                !###############################################################
                if (abs(mmu).lt.tol) then 
                        !-- Angular-mmu is neglegible --> ax**2+-b=0
                    sign_aalpha=sign(1.0d0,real(aalpha));sign_eenergy=sign(1.0d0,real(eenergy))
                    if(sign_aalpha.eq.sign_eenergy) then 
                            !-- The solution belongs to the complex numbers. Skip solution
                        if(Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                            write(*,'(" (α-FDB-FxFy) - Complex positive root:",xF12.6," +i",xF12.6," a.u")') +1.0d0*sqrt(-eenergy/aalpha)
                            write(*,'(" (α-FDB-FxFy) - Complex negative root:",xF12.6," +i",xF12.6," a.u")') -1.0d0*sqrt(-eenergy/aalpha)
                            write(*,*)
                        end if
                        continue
                    else 
                        !-- Double signed root --> ax**2 -b=0
                                !-- We define a complex number for each of the signs and whether the number has a value different or
                                !   equal to zero, it is going to be applied a factor of sqrt(2) so we can still conserve the modulus
                                !   of the solution
                        tmp_root(i,1)=sqrt(-real(eenergy)/real(aalpha));tmp_root(i,2)=-sqrt(-real(eenergy)/real(aalpha))
                        root(i,3)=phi
                        if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                            write(*,'(" (α-FDB-FxFy) - Neglegible dipole moment positive root:",xF18.6)') tmp_root(i,1)
                            write(*,'(" (α-FDB-FxFy) - Neglegible dipole moment negative root:",xF18.6)') tmp_root(i,1)
                            write(*,*)
                        end if
                    end if
                else if (abs(aalpha).lt.tol) then 
                        ! Angular-aalpha is neglegible --> ax**1+b=0
                    tmp_root(i,1)=-real(eenergy)/real(mmu);tmp_root(i,2)=0.0d0
                    root(i,3)=phi
                    if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                        write(*,'(" (α-FDB-FxFy) - Negelgible alpha solution:",xF18.6)') tmp_root(i,1)
                    end if
                else 
                        ! Actual second degree equation --> ax**2+bx**1+c=0
                    discriminant=real(mmu)**2-4.0d0*real(aalpha)*real(eenergy)
                    if (discriminant.ge.0.0d0) then !-- The solution is not going to be complex
                        tmp_root(i,1)=-mmu+sqrt(discriminant)/2.0d0/aalpha
                        tmp_root(i,2)=-mmu-sqrt(discriminant)/2.0d0/aalpha
                        root(i,3)=phi
                    end if
                    if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                        write(*,'(" (α-FDB-FxFy) - Positive pure second degree root:",xF18.6)') -mmu+sqrt(discriminant)/2.0d0/aalpha
                        write(*,'(" (α-FDB-FxFy) - Negative pure second degree root:",xF18.6)') -mmu-sqrt(discriminant)/2.0d0/aalpha
                        write(*,*)
                    end if
                end if
                if (abs(tmp_root(i,1)).eq.abs(tmp_root(i,2))) then ! The root is double signed
                    root(i,1)=abs(tmp_root(i,1)) !-- Absolute value because there is no need of negative modulus
                                                 !   the direction of the electric field can be tuned with the theta and phi
                end if  
                    ! Declare tmp_mu and tmp_alpha as 0 to restart the computation of the roots
                tmp_mu=0.0d0;tmp_alpha=0.0d0
            end do
        end if 
    end if

!================================= α-VOAM-2D ==========================================
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!================================= β-VOAM-2D ==========================================

    if (order_nlop.eq.3) then
            ! XY=(1,1,0) gpc=3 // XZ=(1,0,1) gpc=2 // YZ=(0,1,1) gpc=1
        step=nint(1+abs(cos(0.5d0*PI*gpc)))
        tmp_E=E_0-abs(redox_coef)*redox_potential-target_barrier/6.2751d2     !-- Energy associated to thermochemistry
        tmp_mu=0.0d0;tmp_alpha=0.0d0;tmp_beta=0.0d0
            !##############################################################
                !-- Scanning the XZ or YZ plane
            !##############################################################
        if (axis(3).eq.1) then ! Scanning either XZ or YZ
            do i=1,grid
                    !-- Theta: 0.le.θ.le.π
                theta=((i-1)/(1.0d0*grid-1))*PI
                if (axis(1).eq.1) then
                    phi=0.0d0*PI;root(i,3)=phi !-- Scanning the XZ plane
                else if (axis(2).eq.1) then
                    phi=0.5d0*PI;root(i,3)=phi !-- Scanning the YZ plane
                end if

                    !-- Declare F as a variable vector
                F(1)=dsin(theta)*cos(phi);F(2)=dsin(theta)*dsin(phi);F(3)=dcos(theta)

                    !-- Initialization of angular alpha's calculation
                tmp_alpha=tmp_alpha+2.0d0*F(2-axis(1))*F(2+axis(3))*alpha(2-axis(1),2+axis(3))-beta(2-axis(1),2+axis(3),gpc)*radius*F(gpc)
                do j=(2-axis(1)),(2+axis(3)),step
                    tmp_mu=tmp_mu+F(j)*(mu(j)+alpha(j,gpc)*radius*F(gpc)-beta(j,gpc,gpc)*(radius**2)*F(gpc)**2)
                    tmp_alpha=tmp_alpha+F(j)**2*(alpha(j,j)+beta(j,j,gpc)*radius*F(gpc))
                    tmp_beta=tmp_beta+beta(j,j,j)*F(j)**3+3.0d0*beta(2-axis(1),j,2+axis(3))*F(2-axis(1))*F(2+axis(3))*F(j)
                end do
        
                    !-- Declaration of the cubic equation parameters
                eenergy=cmplx(tmp_E,0.0d0); mmu=-cmplx(tmp_mu,0.0d0)
                aalpha=-0.5d0*cmplx(tmp_alpha,0.0d0)
                bbeta=(1.0d0/6.0d0)*cmplx(tmp_beta,0.0d0)

                    !-- Computation of the specific parameters for the third degree solution
                p=(3.0d0*bbeta*mmu-aalpha**2)/3.0d0/bbeta**2
                q=(2.0d0*aalpha**3-9.0d0*bbeta*aalpha*mmu+27.0d0*eenergy*bbeta**2)/27.0d0/bbeta**3
                omega=real(0.25d0*q**2+(1.0d0/27.0d0)*p**3)
    
                    !-- Print the coefficients for the third degree equation
                if (Guillem(4).eqv..TRUE.) then
                    write(*,'(" Angles: Theta=",xF6.4," and Phi=",xF6.4)') theta,phi
                    write(*,'(" Energy coefficient:             ",xF18.6)') real(eenergy)
                    write(*,'(" Dipole moment coefficient:      ",xF18.6)') real(mmu)
                    write(*,'(" Polarizability coefficient:     ",xF18.6)') real(aalpha)
                    write(*,'(" Hyperpolarizability coefficient:",xF18.6)') real(aalpha)
                    write(*,'(" P=",xF18.6,", Q=",xF18.6," and ω=",xF18.6," parameters")') p,q,omega
                    write(*,*)
                end if

                !###############################################################
                    !-- There is a single real root
                !###############################################################
                if(omega.gt.0.0d0) then
                    tmp_root(1,1)=sign(abs(-0.5d0*q-sqrt(0.25d0*q**2+p**3/27.0d0))**(1.0d0/3.0d0),real(-0.5d0*q-sqrt(0.25d0*q**2+p**3/27.0d0)));tmp_root(1,2)=0.0d0
                    tmp_root(1,1)=tmp_root(1,1)+sign(abs(-0.5d0*q+sqrt(0.25d0*q**2+p**3/27.0d0))**(1.0d0/3.0d0),real(-0.5d0*q+sqrt(0.25d0*q**2+p**3/27.0d0)))
                    tmp_root(1,1)=tmp_root(1,1)-aalpha/3.0d0/bbeta
                
                        !-- Print the specific raw solutions
                    if (Guillem(6).eqv..TRUE.) then
                        write(*,'(" Angles: Theta=",xF6.4," and Phi=",xF6.4)') theta,phi
                        write(*,'(" (β-FDB-"A2,A2") - Unique specific solution:",xF18.6)') axis_name(2-axis(1)),axis_name(2+axis(3)),real(root(i,1))
                        write(*,*) "Complex solutions skipped."
                        write(*,*)
                    end if

                        !-- Checking whether the solution is inside the radius of the sphere
                    if (abs(root(1,1)).gt.max_radius) then
                        continue
                    else !it is inside
                        root(i,1)=real(tmp_root(1,1))
                    end if

                !###############################################################
                    !-- There are three real roots, two of them equal
                !###############################################################
                else if (abs(omega).lt.tol) then !-- Three real roots, two of them are equal
                    tmp_root(1,1)=2.0d0*sign(abs(-0.5d0*q)**(1.0d0/3.0d0),real(-0.5d0*q))-aalpha/3.0d0/bbeta;tmp_root(2,1)=0.0d0
                    tmp_root(2,1)=1.0d0*sign(abs(-0.5d0*q)**(1.0d0/3.0d0),real(-0.5d0*q))-aalpha/3.0d0/bbeta;tmp_root(2,2)=0.0d0
                    tmp_root(3,1)=tmp_root(2,1);tmp_root(3,2)=0.0d0

                        !-- Print the specific raw solutions
                    if (Guillem(6).eqv..TRUE.) then
                        write(*,'(" Angles: Theta=",xF6.4," and Phi=",xF6.4)') theta,phi
                        write(*,'(" (β-FDB-"A2,A2") - First specific solution: ",xF18.6)') axis_name(2-axis(1)),axis_name(2+axis(3)),real(tmp_root(1,1))
                        write(*,'(" (β-FDB-"A2,A2") - Second specific solution:",xF18.6)') axis_name(2-axis(1)),axis_name(2+axis(3)),real(tmp_root(2,1))
                        write(*,'(" (β-FDB-"A2,A2") - Third specific solution: ",xF18.6)') axis_name(2-axis(1)),axis_name(2+axis(3)),real(tmp_root(3,1))
                        write(*,*)
                    end if

                        !-- Checking whether the solutions are inside the sphere
                    if (abs(tmp_root(1,1)).gt.max_radius.and.abs(tmp_root(2,1)).gt.max_radius) then !-- Both solutions are outside the sphere radius
                        continue
                    else
                        root(i,1)=real(tmp_root(1,1))
                        if (abs(tmp_root(2,1)).lt.root(i,1)) root(i,1)=abs(tmp_root(2,1))
                    end if

                !###############################################################
                    !-- There three different real roots
                !###############################################################
                else if (omega.lt.0.0d0) then
                    angle=dacos(-0.5d0*real(q)/sqrt(-1.0d0*real(p)**3/27.0d0))
                    tmp_root(1,1)=2.0d0*sqrt(-p/3.0d0)*cos(angle/3.0d0);tmp_root(1,2)=0.0d0
                    tmp_root(2,1)=2.0d0*sqrt(-p/3.0d0)*cos((angle+2.0d0*PI)/3.0d0);tmp_root(2,2)=0.0d0
                    tmp_root(3,1)=2.0d0*sqrt(-p/3.0d0)*cos((angle+4.0d0*PI)/3.0d0);tmp_root(3,2)=0.0d0

                        !-- Print the specific raw solutions
                    if (Guillem(6).eqv..TRUE.) then
                        write(*,'(" Angles: Theta=",xF6.4," and Phi=",xF6.4)') theta,phi
                        write(*,'(" (β-FDB-FxFy) - First specific solution: ",xF18.6)') real(tmp_root(1,1))
                        write(*,'(" (β-FDB-FxFy) - Second specific solution:",xF18.6)') real(tmp_root(2,1))
                        write(*,'(" (β-FDB-FxFy) - Third specific solution: ",xF18.6)') real(tmp_root(3,1))
                        write(*,*)
                    end if

                        !-- Check whether there are solutions inside the scanning sphere or not
                    if (abs(tmp_root(1,1)).gt.max_radius.and.abs(tmp_root(2,1)).gt.max_radius.and.abs(tmp_root(3,1)).gt.max_radius) then ! All solutions are outside the sphere radius. Skip solution
                        continue
                    else !at least one of them is inside the sphere
                        root(i,1)=real(tmp_root(1,1))
                        do j=2,3 ! Get the minimum solution for each iteration
                            if (abs(tmp_root(j,1)).lt.abs(root(i,1))) root(i,1)=abs(tmp_root(j,1))
                        end do
                    end if
                end if
    
                    !-- Setting tmp_beta, tmp_alpha and tmp_mu to zero for the following iteration
                tmp_mu=0.0d0;tmp_alpha=0.0d0;tmp_beta=0.0d0
            end do
        else if (axis(3).eq.0) then ! Scanning XY
            do i=1,2*grid
                    ! Phi: 0.le.φ.le.2π
                phi=((i-1)/(2.0d0*grid-1))*2.0d0*PI;theta=0.5d0*PI; root(i,2)=theta; root(i,3)=phi

                    !-- Declare F as a variable vector
                F(1)=dsin(theta)*dcos(phi);F(2)=dsin(theta)*dsin(phi);F(3)=dcos(theta)
                tmp_E=tmp_E-mu(gpc)*radius*F(gpc)-0.5d0*(radius**2)*alpha(gpc,gpc)*F(gpc)**2+(1.0d0/6.0d0)*(radius**3)*beta(gpc,gpc,gpc)*F(gpc)**3 !-- Energy associated to constant fields

                    !-- Initialization of angular alpha's calculation
                tmp_alpha=tmp_alpha+2.0d0*F(2-axis(1))*F(2+axis(3))*(alpha(2-axis(1),2+axis(3))-beta(2-axis(1),2+axis(3),gpc)*radius*F(gpc))
                do j=(2-axis(1)),(2+axis(3)),step
                    tmp_mu=tmp_mu+F(j)*(mu(j)+alpha(j,gpc)*radius*F(gpc)-beta(j,gpc,gpc)*(radius**2)*F(gpc)**2)
                    tmp_alpha=tmp_alpha+F(j)**2*(alpha(j,j)+beta(j,j,gpc)*radius*F(gpc))
                    tmp_beta=tmp_beta+beta(j,j,j)*F(j)**3+3.0d0*beta(2-axis(1),j,2+axis(3))*F(2-axis(1))*F(2+axis(3))*F(j)
                end do 

                    !-- Declaration of the cubic equation parameters
                eenergy=cmplx(tmp_E,0.0d0); mmu=-cmplx(tmp_mu,0.0d0)
                aalpha=-0.5d0*cmplx(tmp_alpha,0.0d0)
                bbeta=(1.0d0/6.0d0)*cmplx(tmp_beta,0.0d0)

                    !-- Computation of the specific parameters for the third degree solution
                p=(3.0d0*bbeta*mmu-aalpha**2)/3.0d0/bbeta**2
                q=(2.0d0*aalpha**3-9.0d0*bbeta*aalpha*mmu+27.0d0*eenergy*bbeta**2)/27.0d0/bbeta**3
                omega=real(0.25d0*q**2+(1.0d0/27.0d0)*p**3)

                    !-- Print the coefficients for the third degree equation
                if (Guillem(4).eqv..TRUE.) then
                    write(*,'(" Angles: Theta=",xF6.4," and Phi=",xF6.4)')  theta,phi
                    write(*,'(" Energy coefficient:             ",xF18.6)') real(eenergy)
                    write(*,'(" Dipole moment coefficient:      ",xF18.6)') real(mmu)
                    write(*,'(" Polarizability coefficient:     ",xF18.6)') real(aalpha)
                    write(*,'(" Hyperpolarizability coefficient:",xF18.6)') real(aalpha)
                    write(*,'(" P=",xF18.6,", Q=",xF18.6," and ω=",xF18.6," parameters")') p,q,omega
                    write(*,*)
                end if

                !###############################################################
                    !-- There is a single real root
                !###############################################################
                if(omega.gt.0.0d0) then
                    tmp_root(1,1)=sign(abs(-0.5d0*q-sqrt(0.25d0*q**2+p**3/27.0d0))**(1.0d0/3.0d0),real(-0.5d0*q-sqrt(0.25d0*q**2+p**3/27.0d0)));tmp_root(1,2)=0.0d0
                    tmp_root(1,1)=tmp_root(1,1)+sign(abs(-0.5d0*q+sqrt(0.25d0*q**2+p**3/27.0d0))**(1.0d0/3.0d0),real(-0.5d0*q+sqrt(0.25d0*q**2+p**3/27.0d0)))
                    tmp_root(1,1)=tmp_root(1,1)-aalpha/3.0d0/bbeta
                    if (Guillem(6).eqv..TRUE.) then
                        write(*,'(" Angles: Theta=",xF6.4," and Phi=",xF6.4)')  theta,phi
                        write(*,'(" (β-FDB-FxFy) - Unique specific solution: ",xF18.6)') real(root(i,1))
                        write(*,*) "Complex solutions skipped."
                        write(*,*)
                    end if

                        !-- Checking whether the solution is inside the radius of the sphere
                    if (abs(root(1,1)).gt.max_radius) then
                        continue
                    else !it is inside
                        root(i,1)=real(tmp_root(1,1))
                    end if

                !###############################################################
                    !-- There are three real roots, two of them equal
                !###############################################################
                else if (abs(omega).lt.tol) then !-- Three real roots, two of them are equal
                    tmp_root(1,1)=2.0d0*sign(abs(-0.5d0*q)**(1.0d0/3.0d0),real(-0.5d0*q))-aalpha/3.0d0/bbeta;tmp_root(2,1)=0.0d0
                    tmp_root(2,1)=1.0d0*sign(abs(-0.5d0*q)**(1.0d0/3.0d0),real(-0.5d0*q))-aalpha/3.0d0/bbeta;tmp_root(2,2)=0.0d0
                    tmp_root(3,1)=tmp_root(2,1);tmp_root(3,2)=0.0d0

                        !-- Print all the specific raw roots
                    if (Guillem(6).eqv..TRUE.) then
                        write(*,'(" Angles: Theta=",xF6.4," and Phi=",xF6.4)')  theta,phi
                        write(*,'(" (β-FDB-FxFy) - First specific solution: ",xF18.6)') real(tmp_root(1,1))
                        write(*,'(" (β-FDB-FxFy) - Second specific solution:",xF18.6)') real(tmp_root(2,1))
                        write(*,'(" (β-FDB-FxFy) - Third specific solution: ",xF18.6)') real(tmp_root(3,1))
                        write(*,*)
                    end if

                        !-- Checking whether the solutions are inside the sphere
                    if (abs(tmp_root(1,1)).gt.max_radius.and.abs(tmp_root(2,1)).gt.max_radius) then ! -- Both solutions are outside the sphere radius
                        continue
                    else 
                        root(i,1)=real(tmp_root(1,1))
                        if (abs(tmp_root(2,1)).lt.root(i,1)) root(i,1)=abs(tmp_root(2,1))
                    end if
                !###############################################################
                    !-- There three different real roots
                !###############################################################
                else if (omega.lt.0.0d0) then !-- Three real non-equal roots
                    angle=dacos(-0.5d0*real(q)/sqrt(-1.0d0*real(p)**3/27.0d0))
                    tmp_root(1,1)=2.0d0*sqrt(-p/3.0d0)*cos(angle/3.0d0);tmp_root(1,2)=0.0d0
                    tmp_root(2,1)=2.0d0*sqrt(-p/3.0d0)*cos((angle+2.0d0*PI)/3.0d0);tmp_root(2,2)=0.0d0
                    tmp_root(3,1)=2.0d0*sqrt(-p/3.0d0)*cos((angle+4.0d0*PI)/3.0d0);tmp_root(3,2)=0.0d0
                    if (Guillem(6).eqv..TRUE.) then !-- Print all the specific solutions
                        write(*,'(" Angles: Theta=",xF6.4," and Phi=",xF6.4)')  theta,phi
                        write(*,'(" (β-FDB-FxFy) - First specific solution: ",xF18.6)') real(tmp_root(1,1))
                        write(*,'(" (β-FDB-FxFy) - Second specific solution:",xF18.6)') real(tmp_root(2,1))
                        write(*,'(" (β-FDB-FxFy) - Third specific solution: ",xF18.6)') real(tmp_root(3,1))
                        write(*,*)
                    end if

                        !-- Check whether there are solutions inside the scanning sphere or not
                    if (abs(tmp_root(1,1)).gt.max_radius.and.abs(tmp_root(2,1)).gt.max_radius.and.abs(tmp_root(3,1)).gt.max_radius) then ! All solutions are outside the sphere radius. Skip solution
                        continue
                    else !at least one of them is inside the sphere 
                        root(i,1)=real(tmp_root(1,1))
                        do j=2,3 ! Get the minimum solution for each iteration
                            if (abs(tmp_root(j,1)).lt.abs(root(i,1))) root(i,1)=abs(tmp_root(j,1))
                        end do
                    end if
                end if

                    !-- Setting tmp_beta, tmp_alpha and tmp_mu to zero for the following iteration
                tmp_mu=0.0d0;tmp_alpha=0.0d0;tmp_beta=0.0d0
            end do
        end if    
        ! Fer el sort de la llista de N lements
    end if

!================================= β-VOAM-2D ==========================================

    End subroutine VOAM_2D

!###################################################################################################
    !=========================================================================!
    !============================== 3D REGIME ================================!
    !=========================================================================!
!###################################################################################################

    Subroutine VOAM_3D(Guillem,F,radius,axis,axis_name,order_nlop)
    implicit none
    double precision, dimension(3) :: F
    logical, dimension(10) :: Guillem
    integer, dimension(3) :: axis
    character*2,dimension(3) :: axis_name
    double precision :: radius
        !-- Only for FDB_μ
    double precision :: E_FF,numerator,denominator !-- Only for FDB_μ
        !-- For FDB_α/β 
    integer order_nlop
    if (order_nlop.eq.1) then
        mu=mu+stoi*mu_f(3)
        E_FF=E_0-abs(redox_coef)*redox_potential-(target_barrier/6.2751d2)
        if (Guillem(5).eqv..TRUE.) then
            write(*,'(" (μ-FDB-XYZ) The equation to be solved is: 0 =",xF7.3," +",xF7.3,xA2," +",xF7.3,xA2," +",xF7.3,xA2)') &
            & E_FF,mu(1),axis_name(1),mu(2),axis_name(2),mu(3),axis_name(3)
        end if
        numerator=mu(3)*E_FF; denominator=mu(1)**2+mu(2)**2+mu(3)**2
        F(3)=numerator/denominator;F(1)=(mu(1)/mu(3))*F(3);F(2)=(mu(2)/mu(3))*F(3)
        write(*,'(x"(μ-FDB-XYZ) The minimum field strength is (",A2,",",A2,",",A2,") = (",F7.2,","F7.2,","F7.2,")(*10^-4) a.u")') &
        & axis_name(1),axis_name(2),axis_name(3),F(1)*1.0d4,F(2)*1.0d4,F(3)*1.0d4
    end if
    if (order_nlop.eq.2) then
    end if
    if (order_nlop.eq.3) then
    end if
    End subroutine VOAM_3D
End module VOAM
Module fieldvector
implicit none
double precision :: x0,y0,z0
double precision :: theta,phi
contains
    Function Fx(theta,phi)
    implicit none
    double precision :: Fx,theta,phi,x0
    if (x0.eq.0.0d0) then
        Fx=sin(theta)*cos(phi)
    else
        Fx=x0
    end if
    !write(*,*) "Fx",Fx,theta,phi
    End
!------------------------------------------------------------------------------!
    Function Fy(theta,phi)
    implicit none
    double precision :: Fy,theta,phi,y0
    if (y0.eq.0.0d0) then
        Fy=sin(theta)*sin(phi)
    else
        Fy=y0
    end if
    !write(*,*) "Fy",Fy,theta,phi
    End
!------------------------------------------------------------------------------!
    Function Fz(theta)
    implicit none
    double precision :: Fz,theta,z0
    if (z0.eq.0.0d0) then
        Fz=cos(theta)
    else
        Fz=z0
    end if
    !write(*,*) "Fz",Fz,theta
    End
End module fieldvector
!==============================================================================!
!=======================  MAIN VOAM STARTS HERE  ==============================!
!==============================================================================!
Program mainVOAM
use VOAM
use fieldvector
implicit none
character*10 ET,axis_scan,approximation,stoichiometry
character*80 file ! Name of the file
character*100 line,title ! Characters to read the input file
character*80 name_reactant,name_product,name_ligand ! Name of the involved chemical species
character*10 mister
!---- Integers for the input reading
integer :: index_NLOP,index_scan !-- Indices  for the method line
integer :: index_barrier,index_stoi,index_redox,index_potential !-- Indices for the thermochemistry line
integer :: index_initial,index_modulus,index_trust,index_grid
integer :: index_tol
!---- Integers for Guillem
integer :: index_Guillem
!---- Other
!----------------Comments for future implementations----------------!
    !1)     !--- Programar els minimitzadors amb móduls
                !--- En principi, el fet de tenir tots els elements en memòria, no hauria de suposar un problema per cridar-ho tot dins d'un módul
    !2)     !--- Refer tota la lectura d'input perquè llegeixi qualsevol tipus
            !    de reacció
    !4)     !--- Implementar una nova keyword amb dues opcions: opt i scan (?) --> Calcular el valor de la dG per un rang de camps (-modulus,+modulus) en grid de 500. 
            !    Pels casos 2D es força que el grid sigui una fracció del grid original per disminuir l'emmagatzematge en memòria    
!-------------------------------------------------------------------!
call getarg(1,file)
open (2,file=TRIM(file),status="old")
!open (3,file="VOAM.csv",status="new") --> Only for the 2D scans!!!!!!!

!=============================================================================!
!              Conditional statements for Guillem's purposes                  !
    Guillem(1)=print_NLOP_react;Guillem(2)=print_NLOP_prod
    Guillem(3)=print_NLOP_lig;  Guillem(4)=print_coeffs;
    Guillem(5)=print_EQ;        Guillem(6)=print_sols; 
    Guillem(7)=print_checks;    !Guillem(8)=
    !Guillem(9)=                !Guillem(10)=
!                                                                             !
!=============================================================================!


!================================KEYWORDS BLOCK================================!
!-----------------------Basic/Reaction settings of the program-----------------!
read(2,'(A80)') title
read(2,'(A80)') line !-- Line before the "route section"
read(2,'(A100)') line !-- Line corresponding to the 'Method' line
!------------------------------------------------------------------------------!
        !============================!
        !   Reading the method line  !
        !============================!
if(index(line,"Method:").ne.0.or.index(line,"method:").ne.0) then ! Try reading the axis scan and the NLOPs
    !--- Indices for the method-related keywords
    index_scan=index(line,"scan") !-- Detecting the first coincidence for the 'scan' word to detect the scanning axis
    index_NLOP=index(line,"NLOP") !-- Same as index_scan
    index_Guillem=index(line,"Guillem")             ! It's me, hi! I'm the developer, it's me
    index_Guillem=index_Guillem+index(line,"SUDO")
    index_Guillem=index_Guillem+index(line,"sudo")
        !-- Reading axis scan for the minimization
    if(index(line,"Axis scan=").ne.0.or.index(line,"axis scan=").ne.0) then
        read(line(index_scan+5:index_scan+8),*) axis_scan
        select case(axis_scan)
                !====================== 1D CASES ==========================!
            case ("00X","0X0","X00","0x0","x00","00x","X","x","XXX","xxx","xx","XX")
                    !-- Solve for R but with X-oriented NLOP
                !-- Theta is PI halves and phi is zero, unless stated in the initial point
                n_dim=1;gpc=1;axis(gpc)=1;theta=0.5d0*PI; phi=0.0d0
            case ("00Y","0Y0","Y00","0y0","y00","00y","Y","y","YYY","yyy","yy","YY")
                    !-- Solve for R but with Y-oriented NLOP
                !-- Theta and phi are PI halves, unless stated in the initial point
                n_dim=1;gpc=2;axis(gpc)=1;theta=0.5d0*PI; phi=0.5d0*PI
            case ("00Z","0Z0","Z00","0z0","z00","00z","Z","z","ZZZ","zzz","zz","ZZ")
                    !-- Solve for R but with Z-oriented NLOP
                !-- Theta is zero, unless stated in the initial point
                n_dim=1;gpc=3;axis(gpc)=1;theta=0.0d0
                !===================== 2D CASES ===========================!
            case ("XY0","xy0","xy","XY","0XY","0xy","x0y","X0Y")
                n_dim=2;axis=1;gpc=3;axis(gpc)=0;theta=0.5d0*PI
            case ("X0Z","x0z","0xz","xz0","0XZ","XZ0","xz","XZ")
                n_dim=2;axis=1;gpc=2;axis(gpc)=0;phi=2*PI
            case ("0YZ","0yz","YZ0","yz0","ZY0","zy0","YZ","yz","y0z","Y0Z")
                n_dim=2;axis=1;gpc=1;axis(gpc)=0;phi=0.5d0*PI
                !==================== 3D CASE  ============================!
            case ("XYZ","xyz")
                n_dim=3;axis=1
            case default
                write(*,*) "Bad input axis scan! STOP!"
                stop 
        end select
    else
        ! Potser generar un cas default?
        write(*,*) "Bad input scan keyword! STOP!"
        stop
    end if
    !
        !-- Reading the order of NLOP
    if(index(line,"NLOP=").ne.0.or.index(line,"nlop=").ne.0) then
        read(line(index_NLOP+5:index_NLOP+10),*) approximation
        select case(approximation)
            case ("Dipole","dipole","mu","Mu")
                order_nlop=1
            case ("Alpha","alpha")
                order_nlop=2
            case ("Beta","beta")
                order_nlop=3
            case ("Gamma","gamma")
                order_nlop=4
            case ("Delta","delta")
                order_nlop=5
            case default
                write(*,*) "Bad input approximation! STOP!"
        end select
    else
        write(*,*) "Bad input NLOP keyword! STOP!"
        stop
    end if
    !
    if (index(line,"SUDO=").ne.0.or.index(line,"sudo=").ne.0) then ! Guillem is here
        read(line(index_Guillem+5:index_Guillem+15),*) mister
        select case (mister)
            case ("Guillem")
                Guillem=.TRUE.
            case ("Reactant")
                Guillem(1)=.TRUE.
            case ("Product")
                Guillem(2)=.TRUE.
            case ("Ligand")
                Guillem(3)=.TRUE.
            case ("Coefficient","coefficient","coeff")
                Guillem(4)=.TRUE.
            case ("Equations","equation","equations","Equation")
                Guillem(5)=.TRUE.
            case ("Solutions","solutions","solution","Solution","sols","Sols","SOLS")
                Guillem(6)=.TRUE.
            case ("Check","check")
                Guillem(7)=.TRUE.
            case ("NLOP","nlop")
                Guillem(1)=.TRUE.;Guillem(2)=.TRUE.;Guillem(3)=.TRUE.
            case ("full check","Full check","fullcheck")
                Guillem(4)=.TRUE.;Guillem(6)=.TRUE.
        end select
    else ! No sudo
        continue
    end if
    !
else
    write(*,*) "Error reading the method line! Stop!"
    stop
end if
!------------------------------------------------------------------------------!
read(2,'(A100)') line !-- Line corresponding to the 'Thermochemistry' line
            !=====================================!
            !   Reading the thermochemistry line  !
            !=====================================!
if(index(line,"Thermochemistry:").ne.0.or.index(line,"thermochemistry:").ne.0) then ! Try reading thermodynamics 
    !--- Indices for the thermochemistry-related keywords
    index_barrier=index(line,"barrier")
    index_stoi=index(line,"chiometry")
    index_redox=index(line,"edox")
    index_potential=index(line,"tential")
        ! -- Reading the thermodynamic barrier (in kcal/mol) that is going to be scanned
    if(index(line,"Target barrier=").ne.0.or.index(line,"target barrier=").ne.0) then
        read(line(index_barrier+8:index_barrier+13),*) target_barrier
        if(target_barrier.le.0.0d0) then ! The program stops because such scenario cannot exist
            write(*,*) "Non suitable barrier! Stop!"
            stop
        end if
    else
        continue
    end if
    !
        ! -- Reading the stoichiometry of the ligand: whether it acts as a reactant, as a product or there is none involved
    if(index(line,"Stoichiometry=").ne.0.or.index(line,"stoichiometry=").ne.0) then
        read(line(index_stoi+10:index_stoi+12),*) stoi
        select case(stoi)
        !-- Value of the stoichiometry parameter abiding to the regular thermochemistry conventions
            case (1) ! Ligand acting as a product
                stoi=1
            case (0) ! Reaction without ligand
                stoi=0
            case (-1) ! Ligand acting as a reactant
                stoi=-1
            case default
                write(*,*) "Not a valid ligand stoichiometry! STOP!"
                stop
        end select
    else
        !stoi=0 --> Default?
        stop
    end if
    !
        ! -- Reading whether there is an ET process, or not, and how it behaves: oxidation or reduction
    if(index(line,"Redox=").ne.0.or.index(line,"redox=").ne.0) then
        read(line(index_redox+5:index_redox+6),*) redox_coef
        select case(redox_coef) !-- Reformular perquè representi el nombre d'electrons en el procés redox
            ! -- The following lines are only redundant statements
            case (1)
                redox_coef=+1 ! +1 in oxidation state = Oxidation
            case (0)
                redox_coef=0 ! No redox process involved
            case (-1)
                redox_coef=-1 ! -1 in oxidation state = Reduction
            case default
                write(*,*) "Not a valid redox reaction! STOP!"
                stop
        end select
    else
        !-- Default value if it is not found
        redox_coef=0
    end if
    !
    if(index(line,"Potential=").ne.0.or.index(line,"potential=").ne.0) then
        read(line(index_potential+8:index_potential+13),*) redox_potential ! vs SHE !-- Reformular segons redox_coeff
                    !-- S'ha de mirar si s'ha de corregir la correcció de la fórmula
            if(redox_potential.ne.0.0) then
                    ! Potential: Converting the SHE red_potential into Gibbs energy
                redox_potential=(-redox_potential*9.6485d4)/4184.0d0-98.6991
                    ! Conversion from Volts to kcal/mol
                    ! 98.6991 --> Conversion of -4.28eV to kcal/mol
                redox_potential=redox_potential*redox_coef/6.2751d2 ! SHE red_potential into a.u adapted for the redox transoformation 
            end if
    else
        write(*,*) "Bad input potential keyword! STOP!"
        stop
    end if
    !
else
    write(*,*) "Error in reading the Thermochemistry line! Stop!"
    stop
end if
!------------------------------------------------------------------------------!
read(2,'(A100)') line !-- Line corresponding to the 'Computation' line
                !=================================!
                !   Reading the computation line  !
                !=================================!
if(index(line,"Computation:").ne.0.or.index(line,"computation:").ne.0) then
    index_initial=index(line,"oint")
    index_modulus=index(line,"dulus")
    index_trust=index(line,"dius")
    index_grid=index(line,"rid") 
    !
    if(index(line,"Initial point=").ne.0.or.index(line,"initial point=").ne.0) then
        read(line(index_initial+5:index_initial+21),*) x0,y0,z0
        x0=x0*1.0d-4;y0=y0*1.0d-4;z0=z0*1.0d-4
    else
        write(*,*) "Bad input of initial point! STOP!"
        stop 
    end if
    !
    if(index(line,"Modulus=").ne.0.or.index(line,"modulus=").ne.0) then
        read(line(index_modulus+6:index_modulus+8),*) radius; max_radius=radius; radius=radius*1.0d-4
    else
        write(*,*) "Bad input of maximum scanning modulus! STOP!"
        stop
    end if
    !
    if(index(line,"Trust radius=").ne.0.or.index(line,"trust radius=").ne.0) then
        read(line(index_trust+5:index_trust+7),*) c_i
        if(c_i.lt.0) then ! Negative trust radius? Do you trust this program so much?
            write(*,*) "Invalid trust radius value. Stop!"
            stop
        end if
    else
        write(*,*) "Bad input of trust radius! STOP!"
        stop
    end if
    !
    if(index(line,"Grid=").ne.0.or.index(line,"grid=").ne.0) then
        read(line(index_grid+4:index_grid+8),*) grid
        if(grid.ge.1E3) write(*,*) "Warning: grid too dense for 3D"
        if(grid.le.1E2) write(*,*) "Warning: grid too thin for scanning"
    else
        write(*,*) "Bad input of grid! STOP!"
        stop
    end if
    !
else
    write(*,*) "Bad input line! Stop!"
    stop
end if
!-----------------------------------------------------------------------------!
read(2,'(A100)') line !-- Line corresponding to the 'Extra' line
                !=================================!
                !       Reading the extra line    !
                !=================================!
if(index(line,"Extra:").ne.0.or.index(line,"extra:").ne.0) then ! Tolerance line exits
    index_tol=index(line,"lerance")
    read(line(index_tol+8:index_tol+14),*) tol
    if(tol.lt.0) then ! Tolerance cannot be lower than zero
        write(*,*) "Tolerance cannot be negative! Stop!"
        stop
    end if
else
    tol=1E-4
end if
    !===================== Specific Guillem keywords
!read(2,'(A100)') line !-- Line corresponding to MY keywords
!if(index(line,"Guillem:").ne.0.or.index(line,"guillem").ne.0) then ! Specific keywords for my purpose
!    index_pNLOP=index(line,"pNLOP") 
!    if(index(line,"SUDO").ne.0.or.index(line,"sudo").ne.0) then ! Activate SUDO access
!        print_NLOP_react=.TRUE.
!        print_NLOP_prod=.TRUE.
!        print_NLOP_lig=.TRUE.
!    end if
!    if(index(line,"pNLOP=").ne.0) then !-- Print NLOPs
!        read(line(index_pNLOP+6:index_pNLOP+7) pNLOP
!        select case(pNLOP)
!            case(1)
!                print_NLOP_react=.true.
!            case(2)
!                print_NLOP_prod=.true.
!            case(3)
!                print_NLOP_lig=.true.
!            case(4)
!                print_NLOP_react=.true.;print_NLOP_prod=.true.
!            case default
!                continue
!        end select         
!else
!    continue
!end if

!-- Defining the electric fields vectors in spherical polar coordinates for the following scan
!-----------------------------------------------------------------------------!
        F(1)=Fx(theta,phi);F(2)=Fy(theta,phi);F(3)=Fz(theta)
        axis_name(1)="Fx";axis_name(2)="Fy";axis_name(3)="Fz"    
                if(sudo.eqv..TRUE.) Guillem = .TRUE.
!-----------------------------------------------------------------------------!

read(2,'(A80)') line ! DO NOT TOUCH THIS LINE! SOMEHOW IT MAKES IT WORK
read(2,'(A80)') line ! Line after the "route section": corresponds to ">> FDB Parameters"
!===========================END OF ROUTE SECTION===============================!


            !-----------------Echoing the input-----------------!
write(*,*) "-------------------------------------------------------------------"
write(*,*) 
write(*,*) "            The input for this run is:"
write(*,'(" Central point ((x,y,z)·10^-4 a.u)",xF8.3,xF8.3,xF8.3)') x0,y0,z0
write(*,'(" Radius of ",F7.5," a.u and confidence radius of",xF4.1"%")') radius,c_i
write(*,'(" Number of points ",I5," (1D) ",I10," (2D) ",I15," (3D) ")') grid,grid**2,grid**3
if (redox_potential.ne.0.0) then
    write(*,'(" Potential (a.u)",xF8.3)') redox_potential
end if
write(*,'(" Chosen order of the approximation for this run:",I2)') order_nlop
if (order_nlop.eq.1) then
    write(*,*) "Linear approximation -- Dipole moment (μ)"
else if (order_nlop.eq.2) then
    write(*,*) "Quadratic approximation -- Polarizability (α)"
else if (order_nlop.eq.3) then
    write(*,*) "Cubic approximation -- Hyperpolarizability (β)"
end if
write(*,*) 
write(*,*) "------------------------------------------------------------------"
read(2,'(A80)') name_reactant; read(2,*) line
!do while(index(line,"-----Products-----") 
!end do
!-----------------Chemical species-----------------!
call readvalues(E_r,mu_r,alpha_r,beta_r,alpha_tmp)
    ! -- Just to check --
if(Guillem(1).eqv..TRUE.) then
    write(*,*) "=========================================="
    write(*,*) "Properties for reactants (R)"
    write(*,*) "------------------------------------------"
    write(*,*) "Dipole moment"
    write(*,*) (mu_r(i),i=1,3)
    write(*,*) "Polarizability matrix"
    do i=1,3
            write(*,*) (alpha_r(i,j),j=1,3)
    end do
    write(*,*) "First hyperpolarizability tensor"
    write(*,*) "X _ _"
    do j=1,3
            write(*,*) (beta_r(1,j,k),k=1,3)
    end do
    write(*,*) "Y _ _"
    do j=1,3
            write(*,*) (beta_r(2,j,k),k=1,3)
    end do
    write(*,*) "Z _ _"
    do j=1,3
            write(*,*) (beta_r(3,j,k),k=1,3)
    end do
    write(*,*) "=========================================="
end if
!-----------------------------------------------------------------------------!
read(2,'(A80)') name_product; read(2,*) line
!do while(index(line,"----Ligands products-----")
!end do
call readvalues(E_p,mu_p,alpha_p,beta_p,alpha_tmp)
if(Guillem(2).eqv..TRUE.) then
    write(*,*) "=========================================="
    write(*,*)
    write(*,*) "Properties for products (P)"
    write(*,*) "------------------------------------------"
    write(*,*) "Dipole moment"
    write(*,*) (mu_p(i),i=1,3)
    write(*,*) "Polarizability matrix"
    do i=1,3
            write(*,*) (alpha_p(i,j),j=1,3)
    end do
    write(*,*) "First hyperpolarizability tensor"
    write(*,*) "X _ _"
    do j=1,3
            write(*,*) (beta_p(1,j,k),k=1,3)
    end do
    write(*,*) "Y _ _"
    do j=1,3
            write(*,*) (beta_p(2,j,k),k=1,3)
    end do
    write(*,*) "Z _ _"
    do j=1,3
            write(*,*) (beta_p(3,j,k),k=1,3)
    end do
    write(*,*) "=========================================="
end if

!--- For ligands the subroutine does not apply:
!--- We are considering ligands orient along their dipole moment and always stabilising the chemical system
!--- Therefore, it has to be done "manually"
    !--- Generalitzar la lectura dels lligands i diferenciar si son reactius o productes
read(2,'(A80)') name_ligand; read(2,*) line
read(2,*) E_f
read(2,*) line
read(2,*) mu_f(1),mu_f(2),mu_f(3)

if (mu_f(3).ge.0.0.or.mu_f(2).ge.0.0.or.mu_f(1).ge.0.0) then ! Check whether the solvent is properly oriented
    continue
    if(mu_f(3).eq.0.0.or.mu_f(2).eq.0.0.or.mu_f(1).eq.0.0) then !
        continue
        if(mu_f(3).eq.0.0.or.mu_f(2).eq.0.0.or.mu_f(1).eq.0.0) then
            ! The molecule is properly oriented!
            continue
        else
            write(*,*) "WARNING! The molecule seems that is not properly oriented"
        end if
    else
        write(*,*) "WARNING! The molecule seems that is not properly oriented"
    end if
else
    write(*,*) "WARNING! The ligand molecule must be oriented towards your reaction axis!"
    write(*,*) "Please, insert a positive value of dipole moment for this molecule"
end if

        ! Unique conversion for ligands
    mu_f(1)=mu_f(3);mu_f(2)=mu_f(3)
read(2,*) line
read(2,*) alpha_f(1,1),alpha_f(1,2),alpha_f(2,2),alpha_f(1,3),alpha_f(2,3)
read(2,*) alpha_f(3,3)
        ! Unique conversion for ligands
    alpha_f(1,1)=alpha_f(3,3);alpha_f(2,2)=alpha_f(3,3)
        ! Conversion of alpha to its transposed values
    alpha_f(2,1)=alpha_f(1,2);alpha_f(3,1)=alpha_f(1,3);alpha_f(3,2)=alpha_f(2,3)
read(2,*) line
read(2,*) beta_f(1,1,1),beta_f(1,1,2),beta_f(1,2,2),beta_f(2,2,2),beta_f(1,1,3)
read(2,*) beta_f(1,2,3),beta_f(2,2,3),beta_f(1,3,3),beta_f(2,3,3),beta_f(3,3,3)
        ! Unique conversion for ligands
    beta_f(1,1,1)=beta_f(3,3,3);beta_f(2,2,2)=beta_f(3,3,3)
        ! Conversion of beta to its transposed values
                        !xxy
    beta_f(2,1,1)=beta_f(1,1,2);beta_f(1,2,1)=beta_f(1,1,2)
                        !xyy
    beta_f(2,1,2)=beta_f(1,2,2);beta_f(2,2,1)=beta_f(1,2,2)
                        !xxz
    beta_f(3,1,1)=beta_f(1,1,3);beta_f(1,3,1)=beta_f(1,1,3)
                        !yyz
    beta_f(3,2,2)=beta_f(2,2,3);beta_f(2,3,2)=beta_f(2,2,3)
                        !xzz
    beta_f(3,1,3)=beta_f(1,3,3);beta_f(3,3,1)=beta_f(1,3,3)
                        !yzz
    beta_f(3,3,2)=beta_f(2,3,3);beta_f(3,2,3)=beta_f(2,3,3)
                        !xyz
    beta_f(1,3,2)=beta_f(1,2,3);beta_f(2,1,3)=beta_f(1,2,3)
    beta_f(2,3,1)=beta_f(1,2,3);beta_f(3,1,2)=beta_f(1,2,3)
    beta_f(3,2,1)=beta_f(1,2,3)
read(2,*) line
read(2,*) alpha_tmp(1,1),alpha_tmp(1,2),alpha_tmp(1,3)
read(2,*) alpha_tmp(2,1),alpha_tmp(2,2),alpha_tmp(2,3)
read(2,*) alpha_tmp(3,1),alpha_tmp(3,2),alpha_tmp(3,3)
do i=1,3
        do j=1,3
                alpha_f(i,j)=alpha_f(i,j)+alpha_tmp(i,j)
        end do
end do

if( Guillem(3).eqv..TRUE.) then
    write(*,*) "=========================================="
    write(*,*)
    write(*,*) "Properties for fakes (F)"
    write(*,*) "------------------------------------------"
    write(*,*) "Dipole moment"
    write(*,*) (mu_f(i),i=1,3)
    write(*,*) "Polarizability matrix"
    do i=1,3
            write(*,*) (alpha_f(i,j),j=1,3)
    end do
    write(*,*) "First hyperpolarizability tensor"
    write(*,*) "X _ _"
    do j=1,3
            write(*,*) (beta_f(1,j,k),k=1,3)
    end do
    write(*,*) "Y _ _"
    do j=1,3
            write(*,*) (beta_f(2,j,k),k=1,3)
    end do
    write(*,*) "Z _ _"
    do j=1,3
            write(*,*) (beta_f(3,j,k),k=1,3)
    end do
    write(*,*) "=========================================="
end if

!-----------------System properties-----------------!
mu=0.0; alpha=0.0; beta=0.0
E_0=E_p-E_r-stoi*E_f ! Energy in atomic units
do i=1,3
    mu(i)=mu_p(i)-mu_r(i)
    do j=1,3
        alpha(i,j)=alpha_p(i,j)-alpha_r(i,j)
        do k=1,3
            beta(i,j,k)=beta_p(i,j,k)-beta_r(i,j,k)
        end do
    end do
end do
    !-- Echoing the properties of the run
write(*,*)
write(*,*) "=============================================================================="
write(*,'("     Title of the job: ",xA80)') title
write(*,*) "------------------------------------------------------------------------------"
write(*,*)
write(*,'(" Reaction studied",xA18,"+",xI2,"*",xA9,"--->",xA18)') name_reactant,-1*stoi,name_ligand,name_product
write(*,*)
write(*,*) "------------------------------------------------------------------------------"
write(*,*)
write(*,'(" Stoichiometry of the ligand:",xxI2)') stoi
write(*,'(" Target barrier for the study:",xF6.2," kcal/mol")') target_barrier
write(*,*) 
write(*,*) "=============================================================================="
write(*,*) "     Thermodynamics and non-linear optical properties of the system"
write(*,*) "------------------------------------------------------------------------------"
write(*,'("    Gibbs free energy (ΔG)",xF15.4,x"kcal/mol")') E_0*6.2751d2
write(*,*) "            Dipole moment vector (μ)"
write(*,'(xxxxxE11.4,xxxxxE11.4,xxxxxE11.4)') (mu(i)+stoi*mu_f(i),i=1,3)
write(*,*)
write(*,*) "          Polarizability matrix (α)"
do i=1,3
        write(*,'(xxxxxE11.4,xxxxxE11.4,xxxxxE11.4)') (alpha(i,j)+stoi*alpha_f(i,j),j=1,3)
end do
write(*,*)
write(*,*) "      First hyperpolarizability tensor (β)"
write(*,*) "X _ _"
do j=1,3
        write(*,'(xxxxxE11.4,xxxxxE11.4,xxxxxE11.4)') (beta(1,j,k)+stoi*beta_f(1,j,k),k=1,3)
end do
write(*,*) "Y _ _"
do j=1,3
        write(*,'(xxxxxE11.4,xxxxxE11.4,xxxxxE11.4)') (beta(2,j,k)+stoi*beta_f(2,j,k),k=1,3)
end do
write(*,*) "Z _ _"
do j=1,3
        write(*,'(xxxxxE11.4,xxxxxE11.4,xxxxxE11.4)') (beta(3,j,k)+stoi*beta_f(3,j,k),k=1,3)
end do
write(*,*)
write(*,*) "------------------------------------------------------------------------------"
write(*,*) "=============================================================================="
write(*,*)
write(*,*) "            Please check everything is correct before running."
write(*,'(" If you wish to make any change at the input file",A20,",")') trim(file)
write(*,*) "                        please press Ctrl+C"
write(*,*) 
!call sleep(3)
!--- S'ha de fer un if segons el n_dim
write(*,'("     Searching the ",A2,x"minimum field strength for the input-imposed conditions")') axis_name(gpc)
write(*,*) "                           Please wait..."
!--- S'ha de fer un if segons el n_dim
write(*,*) 
write(*,*)
write(*,*) "======================================================================================================================"
!=============================================================================!
!----------------------------VOAM STARTS HERE---------------------------------!
!=============================================================================!
if (n_dim.eq.1) then !-- 1D resolution is performed
    call VOAM_1D(Guillem,F,gpc,tol,radius,axis,axis_name,order_nlop)
end if
if (n_dim.eq.2) then !-- 2D resolution is performed
    call VOAM_2D(Guillem,gpc,radius,axis,axis_name,order_nlop,grid)  
end if
if (n_dim.eq.3) then !-- 3D resolution is performed
    call VOAM_3D(Guillem,F,radius,axis,axis_name,order_nlop)
end if
!=============================================================================!
!----------------------------VOAM ENDS HERE ----------------------------------!
!=============================================================================!
write(*,*) "======================================================================================================================"
write(*,*)
write(*,*)
write(*,*) " Program written by Pau Besalú and Guillem Pey* for his Master Thesis"
write(*,*) "            supervised by Dr. Josep Maria Luis Luis in the"
write(*,*) "    Master in Advanced Catalysis and Molecular Modeling (MACMoM)"
write(*,*) "                at Universitat de Girona (UdG) on 2024" 
write(*,*) 
write(*,*) 
write(*,*) 
End
!------------------------------------------------------------------------------!
Subroutine readvalues(E_i,mu_i,alpha_i,beta_i,alpha_tmp)
implicit none
character*100 line
double precision :: E_i
double precision, dimension (3) :: mu_i
double precision, dimension (3,3) :: alpha_i,alpha_tmp
double precision, dimension (3,3,3) :: beta_i
integer i,j

read(2,*) E_i
read(2,*) line
read(2,*) mu_i(1),mu_i(2),mu_i(3)
read(2,*) line
read(2,*) alpha_i(1,1),alpha_i(1,2),alpha_i(2,2),alpha_i(1,3),alpha_i(2,3)
read(2,*) alpha_i(3,3)
        ! Conversion of alpha to its transposed values
alpha_i(2,1)=alpha_i(1,2);alpha_i(3,1)=alpha_i(1,3);alpha_i(3,2)=alpha_i(2,3)

read(2,*) line
read(2,*) beta_i(1,1,1),beta_i(1,1,2),beta_i(1,2,2),beta_i(2,2,2),beta_i(1,1,3)
read(2,*) beta_i(1,2,3),beta_i(2,2,3),beta_i(1,3,3),beta_i(2,3,3),beta_i(3,3,3)
        ! Conversion of beta to its transposed values
                                !xxy
        beta_i(2,1,1)=beta_i(1,1,2);beta_i(1,2,1)=beta_i(1,1,2)
                                !xyy
        beta_i(2,1,2)=beta_i(1,2,2);beta_i(2,2,1)=beta_i(1,2,2)
                                !xxz
        beta_i(3,1,1)=beta_i(1,1,3);beta_i(1,3,1)=beta_i(1,1,3)
                                !yyz
        beta_i(3,2,2)=beta_i(2,2,3);beta_i(2,3,2)=beta_i(2,2,3)
                                !xzz
        beta_i(3,1,3)=beta_i(1,3,3);beta_i(3,3,1)=beta_i(1,3,3)
                                !yzz
        beta_i(3,3,2)=beta_i(2,3,3);beta_i(3,2,3)=beta_i(2,3,3)
                                !xyz
        beta_i(1,3,2)=beta_i(1,2,3);beta_i(2,1,3)=beta_i(1,2,3)
        beta_i(2,3,1)=beta_i(1,2,3);beta_i(3,1,2)=beta_i(1,2,3)
        beta_i(3,2,1)=beta_i(1,2,3)
read(2,*) line
read(2,*) alpha_tmp(1,1),alpha_tmp(1,2),alpha_tmp(1,3)
read(2,*) alpha_tmp(2,1),alpha_tmp(2,2),alpha_tmp(2,3)
read(2,*) alpha_tmp(3,1),alpha_tmp(3,2),alpha_tmp(3,3)
do i=1,3
    do j=1,3
        alpha_i(i,j)=alpha_i(i,j)+alpha_tmp(i,j)
    end do
end do
End
