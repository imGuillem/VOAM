!-- Program written by Pau Besalú and Guillem Pey* - TFM Guillem Pey
!-- VOAM.f95 (VoltOhmAmpereMaxwell.f95) 
!-- OpEEF searcher within the FDB_beta method using the third degree equation, cylindrical and polar spherical coordinates
Module VOAM
implicit none
! Non-linear optical properties (NLOPs) and energies
    double precision :: domain,E_0,E_r,E_p,E_lig_r,E_lig_p,target_barrier
    double precision, dimension(3) :: F,initial_position
    double precision, dimension(3) :: mu,mu_r,mu_p,mu_lig,mu_lig_r,mu_lig_p
    double precision, dimension(3,3) :: alpha,alpha_r,alpha_p,alpha_lig,alpha_lig_r,alpha_lig_p,alpha_tmp
    double precision, dimension(3,3,3) :: beta,beta_r,beta_p,beta_lig,beta_lig_r,beta_lig_p
! VOAM keywords
    integer, dimension (3) :: axis = 0
    double precision :: redox_potential,radius,max_radius,tol
    double precision :: x0,y0,z0,theta,phi
    integer order_nlop,n_dim,x_axis,y_axis,z_axis
    integer stoi,redox_coef,grid 
    integer gpc
! Guillem keywords
    logical, dimension(10) :: Guillem = .FALSE.
    logical :: sudo = .FALSE.
    logical :: print_NLOP_react,print_NLOP_prod,print_NLOP_lig ! Guillem(1-3)
    logical :: print_coeffs = .FALSE.   ! Guillem(4)
    logical :: print_EQ = .FALSE.       ! Guillem(5)
    logical :: print_sols = .FALSE.     ! Guillem(6)
    logical :: print_checks = .FALSE.   ! Guillem(7)
    logical :: print_mins = .FALSE.     ! Guillem(8)
    logical :: print_bubble = .FALSE.   ! Guillem(9)
! Others
    character*2,dimension(3) :: axis_name
    double precision :: PI = 4.0d0*datan(1.0d0)
    double precision :: eps = cmplx(-0.5d0,dsqrt(3.0d0/4.0d0))
    double precision :: Faraday = 9.64853321233100184d4 ! C/mol
    double precision :: hartree2kcal = 6.2751d2
    double precision :: G_SHE = 4.281d0 ! V
    integer i,j,k
contains

    Subroutine VOAM_0D(axis_name,grid,radius) !-- Scanning subroutine
    implicit none
    character*2,dimension(3) :: axis_name
    double precision, allocatable, dimension (:,:) :: dG_field
    double precision :: radius,step,field,G_field,G_main,G_side
    integer i,j,k,grid
    
        !-- Field scanning
    allocate(dG_field(int(grid/10)+1,4))
    step=(2.0d0*radius/(grid*1.0d-1))
    E_0=E_0-abs(redox_coef)*redox_potential
    do i=1,3
        do j=1,(grid/10)+1
            field=-radius+(j-1)*step
            G_field=E_0-(mu(i)+sign(1.0d0,field)*mu_lig(3))*field-0.5d0*(alpha(i,i)+alpha_lig(3,3))*field**2-(1.0d0/6.0d0)*(beta(i,i,i)+sign(1.0d0,field)*beta_lig(3,3,3))*field**3
                !G_side=(E_lig_p-E_lig_r)-mu_lig(3)*abs(field)-0.5d0*alpha_lig(3,3)*field**2-(1.0d0/6.0d0)*beta_lig(3,3,3)*abs(field)**3
                !G_main=(E_p-E_r)-mu(1)*field-0.5*alpha(1,1)*field**2-(1.0d0/6.0d0)*beta(1,1,1)*field**3
                !G_field=G_main_G_side
            dG_field(j,1)=field;dG_field(j,i+1)=G_field*hartree2kcal
        end do
    end do
    write(*,*) "        VOAM - Scan of the input file"
    write(*,*) "=================================================================="
    write(*,'("      F(a.u.)",XXXX"G-"A2,"(kcal/mol)",XXXX"G-"A2,"(kcal/mol)",XXXX"G-"A2,"(kcal/mol)")') axis_name(1),axis_name(2),axis_name(3)
    do i=1,(grid/10)+1
        write(*,'(xxF12.6,xxxxF12.6,xxxxxxxxF12.6,xxxxF12.6)') (dG_field(i,j),j=1,4)
    end do

    deallocate(dG_field)

    End subroutine VOAM_0D
!###################################################################################################
    !=========================================================================!
    !============================== 1D REGIME ================================!
    !=========================================================================!
!###################################################################################################

    Subroutine VOAM_1D(Guillem,order_nlop,gpc,radius,axis,initial_position,axis_name,tol)
    implicit none
    logical,dimension(10) :: Guillem
    double precision, dimension (3) :: F,initial_position
    double precision, allocatable, dimension (:,:) :: dG_field
    character*2,dimension(3) :: axis_name
    integer, dimension(3) :: axis
    double precision :: tmp_energy,tmp_mu,tmp_alpha
    double precision :: radius,domain,field
    integer step,gpc,order_nlop
        !---Specific parameters for the 3rd degree solver---!
    integer i,n,sign_eenergy,sign_aalpha
    double complex, dimension(3) :: root
    double complex :: eenergy,mmu,aalpha,bbeta
    double complex :: p,q,mark,root1,root2,root3
    double precision :: tol,discriminant
    double precision :: root_1,root_2,root_3,half_root1,half_root2,min_root
    double precision :: omega,angle,pp,qq,markk,markkk

            !--         Reference of the 3rd degree solution:
            !-- http://olewitthansen.dk/Mathematics/The_cubic_equation.pdf

        !-- Defining whether it is going to be used the negative- or positive-field definition of the 1D-ΔG(F) expression
    allocate(dG_field(21,2))
    E_0=E_0-abs(redox_coef)*redox_potential
    do i=1,21
            !-- Computing an internal VOAM_0D to know the shape of the 1D-ΔG(F) function
        field=-radius+(i-1)*1.0d-3
        dG_field(i,2)=E_0-(mu(gpc)+sign(1.0d0,field)*mu_lig(3))*field-0.5d0*(alpha(gpc,gpc)+alpha_lig(3,3))*field**2-(1.0d0/6.0d0)*(beta(gpc,gpc,gpc)+sign(1.0d0,field)*beta_lig(3,3,3))*field**3
        dG_field(i,1)=field
    end do
    do i=1,21
            !-- Defining the region where the solution should be with a 3 kcal/mol threshold
        if (abs(dG_field(i,2)-target_barrier/hartree2kcal)*hartree2kcal.le.3.0d0) then
            field=dG_field(i,1)
            if (abs(field).le.abs(dG_field(i,1))) then
                field=dG_field(i,1); domain=1.0d0*sign(1.0d0,field) 
                exit
            end if
        end if
    end do               
    deallocate(dG_field)
 
    tmp_energy=E_0-target_barrier/hartree2kcal                                      !-- Energy
    mu(gpc)=mu(gpc)+domain*mu_lig(3)                                                !-- Dipole moment
    alpha(gpc,gpc)=alpha(gpc,gpc)+alpha_lig(3,3)                                    !-- Polarizability matrix
    beta(gpc,gpc,gpc)=beta(gpc,gpc,gpc)+domain*beta_lig(3,3,3)                      !-- Hyperpolarizability tensor
    
        !-- Initialization of the iterable NLOPs
    tmp_mu=mu(gpc);tmp_alpha=0.5d0*alpha(gpc,gpc)                           !-- Tmp values for variations in case of constant fields

                !--                 Step of the do loop according to the gpc value:
                !-- Generated through as a two-valued function that can only yield values of 1 or 2
    step=nint(1+abs(cos(0.5d0*PI*gpc)))    
                            !-- Value of the dummy indices for the nested loops
            ! X_axis(step=1) = 2,3 // 2,3; Y_axis(step=2) = 1,3 // 1,3; Z_axis(step=1) = 1,2 // 1,2
            ! X_axis--> gpc=1            ; Y_axis--> gpc=2            ; Z_axis--> gpc=3
            ! X_axis--> axis(1)=1        ; Y-axis--> axis(2)=2        ; Z_axis--> axis(3)=3

        !-- Declaration of the constant fields
    F(1+axis(1))=initial_position(1+axis(1));F(3-axis(3))=initial_position(3-axis(3))
    max_radius=dsqrt(radius**2-F(1+axis(1))**2-F(3-axis(3))**2)
     
        !-- "General" formula for the computation of the parameters of the equations within the FDB beta method
    do i=(1+axis(1)),(3-axis(3)),step
        do j=(1+axis(1)),(3-axis(3)),step
                !-- Constant-field energy coefficient
            tmp_energy=tmp_energy-mu(i)*F(i)-0.5d0*alpha(i,j)*F(i)*F(j)-(1.0d0/6.0d0)*beta(i,i,i)*F(i)**3-0.5d0*beta(i,i,j)*F(j)*F(i)**2 
    
                !-- Constant-field dipole moment coefficient
            tmp_mu=tmp_mu+alpha(gpc,i)*F(i)+0.5d0*beta(gpc,i,j)*F(i)*F(j)

                !-- Constant-field polarizability matrix coefficient
            tmp_alpha=tmp_alpha+0.5d0*beta(gpc,gpc,i)*F(i) 
        end do
    end do

        !-- Declaration of each parameter as a complex number for its proper handling: defined as negative because of the expansion coefficients
    eenergy=cmplx(tmp_energy,0.0d0)
    mmu=-cmplx(tmp_mu,0.0d0)
    aalpha=-cmplx(tmp_alpha,0.0d0)
    bbeta=-(1.0d0/6.0d0)*cmplx(beta(gpc,gpc,gpc),0.0d0) 
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

!================================ μ-VOAM-1D ==========================================    

        ! --- General solution for 1st degree equation ---
    if (order_nlop.eq.1.) then
            !-- Print the equation to be solved
        if(Guillem(5).eqv..TRUE.) then
            write(*,'(" (μ-FDB-",A2,") - Equation to solve":,xF12.6"x^1 +",xF12.6,"x^0 = 0")') axis_name(gpc),real(mmu),real(eenergy)
        end if

            !-- Check whether the solution is lower than the maximum radius
        root_1=-real(eenergy)/real(mmu)
        if (root_1.le.max_radius) then
            write(*,'(" (μ-FDB-",A2,") - The minimum field strength is ",A2,"=",xF18.2,"(*10^-4) a.u")') axis_name(gpc),axis_name(gpc),real(root1)*1.0d4
            write(*,*)
        else
            write(*,'(" (μ-FDB-",A2,") - Maximum radius:",xF10.6)') axis_name(gpc),max_radius
            write(*,'(" (μ-FDB-",A2,") - The solution is outside of the proposed radius: ",A2," = ",xF18.6,"(*10^-4) a.u)")') axis_name(gpc),axis_name(gpc),min_root*1.0d4
            write(*,*)
        end if
    end if  

!================================ μ-VOAM-1D ==========================================    
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!================================ α-VOAM-1D ==========================================    

        ! --- General solution for 2nd degree equation ---
    if (order_nlop.eq.2) then
            !-- Print the equation to be solved
        if(Guillem(5).eqv..TRUE.) then 
            write(*,'(" (α-FDB-",A2,") - Equation to solve:",xF12.6,"x^2 +",xF12.6,"x^1 +",xF12.6,"x^0 = 0")') axis_name(gpc),real(aalpha),real(mmu),real(eenergy)
            write(*,*)
        end if

            !-- Neglegible dipole moment coefficient
        if (abs(mmu).lt.tol) then 
            write(*,*) "--------------------------------------------------------------------------------"
            write(*,*) "         Alternative second order solution - Neglegible dipole moment           "
            write(*,*) "--------------------------------------------------------------------------------"
            sign_aalpha=sign(1.0d0,real(aalpha));sign_eenergy=sign(1.0d0,real(eenergy))

            if(sign_aalpha.eq.sign_eenergy) then
                write(*,'(" (α-FDB-",A2,") - The system does not have real solutions")') axis_name(gpc)
                write(*,'(" (α-FDB-",A2,") - First pure complex solution :",xF15.6," +i",xF15.6," a.u")') axis_name(gpc),+1.0d4*sqrt(-eenergy/aalpha)
                write(*,'(" (α-FDB-",A2,") - Second pure complex solution:",xF15.6," +i",xF15.6," a.u")') axis_name(gpc),-1.0d4*sqrt(-eenergy/aalpha)
                write(*,*)

            else
                root1=sqrt(-eenergy/aalpha);root2=-root1
                if (Guillem(6).eqv..TRUE.) then
                    write(*,'(" (α-FDB-",A2,") - First pure real solution :",xF15.6," a.u")') axis_name(gpc),real(root1) 
                    write(*,'(" (α-FDB-",A2,") - Second pure real solution:",xF15.6," a.u")') axis_name(gpc),real(root2)
                    write(*,*)
                end if
            end if            

            !-- Neglegible polarizability coefficient
        else if (abs(aalpha).lt.tol) then
            write(*,*) "--------------------------------------------------------------------------------"
            write(*,*) "        Alternative second order solution - Neglegible polarizability           "
            write(*,*) "--------------------------------------------------------------------------------"
            root1=-real(eenergy)/real(mmu); root2=1E+9

            !-- Second degree equation
        else
            discriminant=mmu**2.0d0-4.0d0*aalpha*eenergy
            if (discriminant.ge.0.0d0) then
                root1=(-mmu+sqrt(discriminant))/2.0d0/aalpha
                root2=(-mmu-sqrt(discriminant))/2.0d0/aalpha
                if(Guillem(6).eqv..TRUE.) then ! Specific solutions do not have the 1.0d4 factor
                    write(*,'(" (α-FDB-",A2,") - Discriminant of the polynomial:",xF10.4)') axis_name(gpc),discriminant
                    write(*,'(" (α-FDB-",A2,") - Specific solution 1:",xF18.6," +i*",xF18.6," a.u")') axis_name(gpc),real(root1)
                    write(*,'(" (α-FDB-",A2,") - Specific solution 2:",xF18.6," +i*",xF18.6," a.u")') axis_name(gpc),real(root2)
                    write(*,*)
                end if
            else
                write(*,'(" (α-FDB-",A2,") - There does not exist a real solution for the propose system and input approximation") ') axis_name(gpc)
                write(*,'(" (α-FDB-",A2,") - Discriminant of the polynomial:",xF18.6)') axis_name(gpc),discriminant
                write(*,'(" (α-FDB-",A2,") - Specific solution 1:",xF15.6," +i*",xF15.6,"(*10^-4) a.u")') axis_name(gpc),1.0d4*(-mmu+sqrt(mmu**2.0d0-4.0d0*aalpha*eenergy))/2.0d0/aalpha
                write(*,'(" (α-FDB-",A2,") - Specific solution 2:",xF15.6," +i*",xF15.6,"(*10^-4) a.u")') axis_name(gpc),1.0d4*(-mmu-sqrt(mmu**2.0d0-4.0d0*aalpha*eenergy))/2.0d0/aalpha
                write(*,*)
            end if
        end if

            !-- Get the closest-to-zero solution 
        min_root=abs(root1)
        if (abs(root2).le.min_root) min_root=real(root2)

            !-- Check whether the solutions are inside the resolution sphere or not
        if(min_root.le.max_radius) then
            write(*,'(" (α-FDB-",A2,") - The minimum field strength for the proposed system is:",xA2,"=",xF12.2,"(*10^-4) a.u")') axis_name(gpc),axis_name(gpc),min_root*1.0d4
            write(*,*)
        else 
            write(*,'(" (α-FDB-",A2,") - Maximum radius:",xF10.6)') max_radius
            write(*,'(" (α-FDB-",A2,") - The solution is outside of the proposed radius: ",A2," = ",xF18.6,"(*10^-4) a.u)")') axis_name(gpc),axis_name(gpc),min_root*1.0d4
            write(*,*)
        end if
    
            !-- Print the checks of the solutions
        if (Guillem(7).eqv..TRUE.) then
            write(*,'(" (α-FDB-",A2,") - First root: ",xF18.6," +i",xF18.6," a.u ")') axis_name(gpc),root1
            write(*,*) abs(aalpha*root1**2+mmu*root1+eenergy)
            write(*,'(" (α-FDB-",A2,") - Second root:",xF18.6," +i",xF18.6," a.u ")') axis_name(gpc),root2
            write(*,*) abs(aalpha*root2**2+mmu*root2+eenergy)
        end if

    end if

!================================ α-VOAM-1D ==========================================    
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!================================ β-VOAM-1D ==========================================    

        ! --- General solution for 3rd degree equation ---
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

            !-- Neglegible hyperpolarizability
        if (abs(bbeta).lt.tol**2) then
            
                !-- Neglegible hyperpolarizability and dipole moment coefficients
            if (abs(mmu).lt.tol**2) then 
                write(*,*) "-------------------------------------------------------------------------------------------------------"
                write(*,*) "         Alternative third order solution - Neglegible hyperpolarizability and dipole moment           "
                write(*,*) "-------------------------------------------------------------------------------------------------------"
                sign_aalpha=sign(1.0d0,real(aalpha));sign_eenergy=sign(1.0d0,real(eenergy))

                if(sign_aalpha.eq.sign_eenergy) then
                    write(*,'(" (β-FDB-",A2,") - The system does not have real solutions")') axis_name(gpc)
                    write(*,'(" (β-FDB-",A2,") - First pure complex solution :",xF15.6," +i",xF15.6," a.u")') axis_name(gpc),+1.0d4*sqrt(-eenergy/aalpha)
                    write(*,'(" (β-FDB-",A2,") - Second pure complex solution:",xF15.6," +i",xF15.6," a.u")') axis_name(gpc),-1.0d4*sqrt(-eenergy/aalpha)
                    write(*,*)

                else
                    root1=sqrt(-eenergy/aalpha);root2=-root1
                    if (Guillem(6).eqv..TRUE.) then
                        write(*,'(" (β-FDB-",A2,") - First pure real solution :",xF15.6," a.u")') axis_name(gpc),real(root1) 
                        write(*,'(" (β-FDB-",A2,") - Second pure real solution:",xF15.6," a.u")') axis_name(gpc),real(root2)
                        write(*,*)
                    end if
                end if            

                !-- Neglegible hyperpolarizability and polarizability coefficients
            else if (abs(aalpha).lt.tol**2) then
                write(*,*) "------------------------------------------------------------------------------------------------------"
                write(*,*) "        Alternative thid order solution - Neglegible hyperpolarizability and polarizability           "
                write(*,*) "------------------------------------------------------------------------------------------------------"
                root1=-real(eenergy)/real(mmu); root2=1E+9

                !-- Only neglegible hyperpolarizability
            else
                write(*,*) "------------------------------------------------------------------------------------------------------"
                write(*,*) "                Alternative thid order solution - Neglegible hyperpolarizability                      "
                write(*,*) "------------------------------------------------------------------------------------------------------"
                discriminant=mmu**2.0d0-4.0d0*aalpha*eenergy

                if (discriminant.ge.0.0d0) then
                    root1=(-mmu+sqrt(discriminant))/2.0d0/aalpha
                    root2=(-mmu-sqrt(discriminant))/2.0d0/aalpha
                    if(Guillem(6).eqv..TRUE.) then ! Specific solutions do not have the 1.0d4 factor
                        write(*,'(" (β-FDB-",A2,") - Discriminant of the polynomial:",xF10.4)') axis_name(gpc),discriminant
                        write(*,'(" (β-FDB-",A2,") - Specific solution 1:",xF18.6," +i*",xF18.6," a.u")') axis_name(gpc),real(root1)
                        write(*,'(" (β-FDB-",A2,") - Specific solution 2:",xF18.6," +i*",xF18.6," a.u")') axis_name(gpc),real(root2)
                        write(*,*)
                    end if
                else
                    write(*,'(" (β-FDB-",A2,") - There does not exist a real solution for the propose system and input approximation") ') axis_name(gpc)
                    write(*,'(" (β-FDB-",A2,") - Discriminant of the polynomial:",xF18.6)') axis_name(gpc),discriminant
                    write(*,'(" (β-FDB-",A2,") - Specific solution 1:",xF15.6," +i*",xF15.6,"(*10^-4) a.u")') axis_name(gpc),1.0d4*(-mmu+sqrt(mmu**2.0d0-4.0d0*aalpha*eenergy))/2.0d0/aalpha
                    write(*,'(" (β-FDB-",A2,") - Specific solution 2:",xF15.6," +i*",xF15.6,"(*10^-4) a.u")') axis_name(gpc),1.0d4*(-mmu-sqrt(mmu**2.0d0-4.0d0*aalpha*eenergy))/2.0d0/aalpha
                    write(*,*)
                end if
            end if

                !-- Get the closest-to-zero solution 
            min_root=abs(root1)
            if (abs(root2).le.min_root) min_root=real(root2)

                !-- Check whether the solutions are inside the resolution sphere or not
            if(min_root.le.max_radius) then
                write(*,'(" (β-FDB-",A2,") - The minimum field strength for the proposed system is:",xA2,"=",xF12.2,"(*10^-4) a.u")') axis_name(gpc),axis_name(gpc),min_root*1.0d4
                write(*,*)
            else 
                write(*,'(" (β-FDB-",A2,") - Maximum radius:",xF10.6)') max_radius
                write(*,'(" (β-FDB-",A2,") - The solution is outside of the proposed radius: ",A2," = ",xF18.6,"(*10^-4) a.u)")') axis_name(gpc),axis_name(gpc),min_root*1.0d4
                write(*,*)
            end if

            if (Guillem(7).eqv..TRUE.) then
                write(*,'(" (β-FDB-",A2,") - First root: ",xF18.6," +i",xF18.6," a.u ")') axis_name(gpc),root1
                write(*,*) abs(aalpha*root1**2+mmu*root1+eenergy)
                write(*,'(" (β-FDB-",A2,") - Second root:",xF18.6," +i",xF18.6," a.u ")') axis_name(gpc),root2
                write(*,*) abs(aalpha*root2**2+mmu*root2+eenergy)
            end if

            !-- Hyperpolarizability coefficient is not neglegible
        else

                !-- Three different scenarios
            !================== CASE 1 CASE 1 CASE 1 CASE 1 CASE 1 CASE 1 CASE ===============!
            if(abs(omega).lt.tol) then ! We have three real roots, two of them being equal // omega.eq.0.0d0
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
                if (abs(min_root).le.max_radius) then
                    write(*,'(" (β-FDB-",A2,") - Case 1 - The minimum field strength is:",xF12.2,"(*10^-4 a.u)")') axis_name(gpc),min_root*1.0d4
                    write(*,*)
                else
                    write(*,'(" (β-FDB-",A2,") - Maximum radius:",xF10.6)') axis_name(gpc),max_radius
                    write(*,'(" (β-FDB-",A2,") - The solution is outside of the proposed radius: ",A2," = ",xF18.6,"(*10^-4) a.u)")') axis_name(gpc),axis_name(gpc),min_root*1.0d4
                    write(*,'(" (β-FDB-",A2,") - Please conside an other stoichiometry")') axis_name(gpc)
                    write(*,*) 
                end if

                if (Guillem(7).eqv..TRUE.) then    
                    do i=1,3
                        write(*,'(" (β-FDB-",A2,") - Root",xI1,":",xF18.6," +i",xF18.6)') axis_name(gpc),i,root(i)
                        write(*,'("      Check:",xF10.8)') abs(bbeta*real(root(i))**3.0d0+aalpha*real(root(i))**2.0d0+mmu*real(root(i))+eenergy)
                    end do
                end if
                write(*,*) "--------------------------------------------------------------------------------------------------------------------  "
                return

                !================== CASE 2 CASE 2 CASE 2 CASE 2 CASE 2 CASE 2 CASE ===============!
            else if (omega.gt.tol) then ! We have one real root and two complex conjugated roots. 
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
                min_root=real(root(1))
                if (abs(min_root).le.max_radius) then ! The solution is out of the scan
                    write(*,'(" (β-FDB-",A2,") - Case 2 - The minimum field strength is:",xA2," = ",xF12.2,"(*10^-4) a.u")') axis_name(gpc),axis_name(gpc),real(root(1))*1.0d4
                    write(*,*)
                else 
                    write(*,'(" (β-FDB-",A2,") - Maximum radius:",xF10.6)') axis_name(gpc),max_radius
                    write(*,'(" (β-FDB-",A2,") - The solution is outside of the proposed radius: ",A2," = ",xF18.6,"(*10^-4) a.u)")') axis_name(gpc),axis_name(gpc),min_root*1.0d4
                    write(*,'(" (β-FDB-",A2,") - Please conside an other stoichiometry")') axis_name(gpc)
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
                        write(*,'("      Check:",xF10.8)') abs(bbeta*real(root(i))**3.0d0+aalpha*real(root(i))**2.0d0+mmu*real(root(i))+eenergy)
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
                root_1=2*sqrt(-real(p)/3.0d0)*dcos(angle/3.0d0)-aalpha/3.0d0/bbeta
                root_2=2*sqrt(-real(p)/3.0d0)*dcos((angle+2*PI)/3.0d0)-aalpha/3.0d0/bbeta
                root_3=2*sqrt(-real(p)/3.0d0)*dcos((angle+4*PI)/3.0d0)-aalpha/3.0d0/bbeta
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
                if (abs(min_root).le.max_radius) then
                    write(*,'(" (β-FDB-",A2,") - Case 3 - The minimum field strength is:",xA2," = ",xF12.2,"(*10^-4) a.u")') axis_name(gpc),axis_name(gpc),min_root*1.0d4
                    write(*,*)
                else
                    write(*,'(" (β-FDB-",A2,") - Maximum radius:",xF10.6)') axis_name(gpc),max_radius
                    write(*,'(" (β-FDB-",A2,") - The solution is outside of the proposed radius: ",A2," = ",xF18.6,"(*10^-4) a.u)")') axis_name(gpc),axis_name(gpc),min_root*1.0d4
                    write(*,'(" (β-FDB-",A2,") - Please conside an other stoichiometry")') axis_name(gpc)
                    write(*,*)
                end if

                    !-- Print the checks of the solutions
                if (Guillem(7).eqv..TRUE.) then 
                    do i=1,3
                       write(*,'(" Root",xI1,":",xF18.6," +i",xF18.6)') i,root(i)
                       write(*,'("      Check:",xF10.8)') abs(bbeta*real(root(i))**3.0d0+aalpha*real(root(i))**2.0d0+mmu*real(root(i))+eenergy)
                    end do
                end if
                write(*,*) "--------------------------------------------------------------------------------------------------------------------  "
                return

            end if
        end if
    end if

!================================ β-VOAM-1D ==========================================    
!#####################################################################################
!#####################################################################################
!#####################################################################################
!#####################################################################################
!================================ γ-VOAM-1D ==========================================    

    if (order_nlop.eq.4) then
            !-- Galois is proud to tell this can still be analytical
        write(*,*) " Gamma approach is not implemented yet!"
        write(*,*) " Please consider lower approximations."
        return
    end if

!================================ γ-VOAM-1D ==========================================    

    End subroutine VOAM_1D

!###################################################################################################
    !=========================================================================!
    !============================== 2D REGIME ================================!
    !=========================================================================!
!###################################################################################################

    Subroutine VOAM_2D(Guillem,order_nlop,gpc,radius,axis,initial_position,axis_name,grid) 
    implicit none
    logical, dimension(10) :: Guillem
    integer, dimension(3) :: axis
    character*2, dimension(3) :: axis_name
    double precision, dimension (3) :: F,initial_position
    double precision :: radius,x0,y0,z0
    integer grid
        !-- Only for FDB_μ
    double precision :: E_FF,numerator,denominator,constant 
        !-- For FDB_α/β 
    double precision, allocatable, dimension (:,:) :: sort ! Sorting matrix for the Bubble Sort algorithm
    double precision, dimension (2) :: tmp_sort ! Temporal matrix for the sorting
    double precision, dimension (3) :: tmp_root ! Stores only the modulus
    double precision, dimension (grid,2) :: root ! Modulus and phi
    double complex :: eenergy,mmu,aalpha,bbeta
    double complex :: p,q
    double precision :: min_root,tmp_energy,tmp_mu,tmp_alpha,tmp_beta
    double precision :: theta,phi,angle,discriminant,omega
    integer :: sign_aalpha,sign_eenergy
        !-- Integers for iteration
    integer gpc,step,order_nlop,upper,lower
    integer i,j,n,dimsort

!================================ μ-VOAM-2D ==========================================    
    
    if (order_nlop.eq.1) then
        upper=nint(0.5d0*axis(1)+1.5d0*axis(2)+0.5d0*axis(3))
        lower=nint(0.5d0*axis(1)+0.5d0*axis(2)+2.5d0*axis(3))
        !-- Falta trobar com F(GPC) es relaciona amb initial_position()
            !-- axis_scan=XY --> (upper,lower)==(2,1);  axis_scan=XZ --> (upper,lower)==(1,3);  axis_scan=YZ --> (upper,lower)==(2,3)
        mu(upper)=mu(upper)+mu_lig(3);mu(lower)=mu(lower)+mu_lig(3)
        E_FF=E_0-abs(redox_coef)*redox_potential-(target_barrier/hartree2kcal)-mu(gpc)*F(gpc)
        if (Guillem(5).eqv..TRUE.) then
            write(*,'(" (μ-FDB-(",A2,",",A2,")) Equation to solve: 0 = ",F12.6," +",xF12.6,xA2," +",xF12.6,xA2)') &
            & axis_name(upper),axis_name(lower),E_FF,mu(upper),axis_name(upper),mu(lower),axis_name(lower)
        end if
        numerator=mu(upper)*E_FF;denominator=mu(upper)**2+mu(lower)**2
        F(upper)=numerator/denominator; F(lower)=(mu(lower)/mu(upper))*F(upper)
        !if (radius**2-F(gpc)**2
        write(*,'(x"The minimum field strength is (",A2,",",A2,";",A2,")=(",F7.2,","F7.2,","F7.2,")(*10^-4) a.u)")') & 
        & axis_name(upper),axis_name(lower),axis_name(gpc),F(upper)*1.0d4,F(lower)*1.0d4,F(gpc)*1.0d4
    end if

!================================= μ-VOAM-2D ==========================================    
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!================================= α-VOAM-2D ==========================================    

    if (order_nlop.eq.2) then

            ! XY=(1,1,0) gpc=3,axis(gpc)=0 --> Step=1 // XZ=(1,0,1) gpc=2,axis(gpc)=0 --> Step=2 // YZ=(0,1,1) gpc=1,axis(gpc)=0 --> Step=1
        step=nint(1+abs(cos(0.5d0*PI*gpc)))

            !-- Initialization of the tmp_NLOP values
        tmp_energy=E_0-abs(redox_coef)*redox_potential-target_barrier/hartree2kcal     !-- Energy associated to thermochemistry

        do i=1,grid
            tmp_mu=0.0d0; tmp_alpha=0.0d0; tmp_root=0.0d0

                !-- 0.le.φ.le.2π 
            phi=((i-1)/(1.0d0*grid-1))*2.0d0*PI

                !-- Declaration of the field vector in cylindrical coordinates
            F(gpc)=initial_position(gpc);F(2-axis(1))=dcos(phi);F(2+axis(3))=dsin(phi)
                    !-- XY plane: F(1)=cos(phi); F(2)=sin(phi); F(3)=Fz
                    !-- XZ plane: F(1)=cos(phi); F(3)=sin(phi); F(2)=Fy
                    !-- YZ plane: F(2)=cos(phi); F(3)=sin(phi); F(1)=Fx
            
                !-- The maximum radius of scan is going to decrease with the increase of the constant field. The general expression if its decrement is the following
            max_radius=dsqrt(radius**2-F(gpc)**2)

                !-- Initialization of the calculation of the (constant-field dependent) energy and angular alpha
            tmp_energy=tmp_energy-mu(gpc)*F(gpc)-0.5d0*alpha(gpc,gpc)*F(gpc)**2
            tmp_alpha=tmp_alpha+2.0d0*alpha(2-axis(1),2+axis(3))*dsin(2.0d0*phi)
                
                !-- Computation of the angular NLOP and declaration of the parameters of the second degree equation
            do j=(2-axis(1)),(2+axis(3)),step
                    !-- Angular dipole moment
                tmp_mu=tmp_mu+F(j)*(mu(j)+alpha(j,gpc)*F(gpc))

                    !-- Angular polarizability matrix
                tmp_alpha=tmp_alpha+alpha(j,j)*F(j)**2

            end do
            eenergy=cmplx(tmp_energy,0.0d0); mmu=-cmplx(tmp_mu,0.0d0); aalpha=-0.5d0*cmplx(tmp_alpha,0.0d0)
            if (Guillem(4).eqv..TRUE.) then
                write(*,'(" Initial point:(",F6.2,",",F6.2,",",F6.2,")" )') initial_position(1),initial_position(2),initial_position(3)
                write(*,'(" Phi=                           ",xF6.4)') phi
                write(*,'(" Energy coefficient:        ",xF18.12)') real(eenergy)
                write(*,'(" Dipole moment coefficient: ",xF18.12)') real(mmu)
                write(*,'(" Polarizability coefficient:",xF18.12)') real(aalpha)
                write(*,*)
            end if 

            !##################################################################!
            !           Computing the second degree solutions                  !
            !##################################################################!

            if (abs(mmu).lt.tol) then ! angular-mmu is neglegible
                sign_aalpha=sign(1.0d0,real(aalpha));sign_eenergy=sign(1.0d0,real(eenergy))

                if(sign_aalpha.eq.sign_eenergy) then 
                        ! The sign is the same ---> ax**2+b=0 ---> The solution is proven to be complex. Skipped
                    if(Guillem(6).eqv..TRUE.) then ! Print the specific solutions
                            !-- Coded in this way so it is not stored in memory and is computed on-the-fly
                        write(*,'(" (α-FDB-",A2,A2,") - Positive complex root:",xF12.6," +i",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),+1.0d0*sqrt(-eenergy/aalpha)
                        write(*,'(" (α-FDB-",A2,A2,") - Negative complex root:",xF12.6," +i",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),-1.0d0*sqrt(-eenergy/aalpha)
                    end if
                    continue
                else
                    !    ax**2-b=0 ---> The solution is proven to be real. Double signed root == tmp_root(i)=positive_root tmp_root(i+1)=negative_root
                    tmp_root(1)=+1.0d0*dsqrt(-real(eenergy)/real(aalpha)); tmp_root(2)=-tmp_root(1)
                    if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                        write(*,'(" (α-FDB-",A2,A2,") - Neglegible dipole moment positive root:",xF18.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),tmp_root(1)
                        write(*,'(" (α-FDB-",A2,A2,") - Neglegible dipole moment negative root:",xF18.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),tmp_root(2)
                        write(*,*)
                    end if
                       
                        !-- Check whether the solution is inside the maximum scanning radius 
                    min_root=abs(tmp_root(1))
                    if (min_root.le.max_radius) then
                        root(i,1)=min_root; root(i,2)=phi
                    end if
                end if

            else if (abs(aalpha).lt.tol) then ! angular-aalpha is neglegible --> ax**1+b=0
                tmp_root(1)=-real(eenergy)/real(mmu)
                min_root=abs(tmp_root(1))
                if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                    write(*,'(" (α-FDB-",A2,A2,") - Neglegible alpha root:",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),min_root
                end if

                    !-- Check whether the solution is inside the maximum scanning radius
                if (min_root.le.max_radius) then
                    root(i,1)=min_root; root(i,2)=phi
                end if

            else ! complete second degree equation
                discriminant=mmu**2.0d0-4.0d0*aalpha*eenergy

                if (discriminant.ge.0.0d0) then
                    tmp_root(1)=-mmu+dsqrt(discriminant)/2.0d0/aalpha
                    tmp_root(2)=-mmu-dsqrt(discriminant)/2.0d0/aalpha
                    if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                        write(*,'(" (α-FDB-",A2,A2,") - Discriminant of the second degree equation:",xF18.12)') axis_name(2-axis(1)),axis_name(2+axis(3)),discriminant
                        write(*,'(" (α-FDB-",A2,A2,") - Positive pure second degree root:",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),tmp_root(1)
                        write(*,'(" (α-FDB-",A2,A2,") - Negative pure second degree root:",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),tmp_root(2)
                        write(*,*)
                    end if

                        !-- Get the closest-to-zero solution; tmp_root(1) and tmp_root(2) are not necessary equal; only when the discriminant is zero
                    min_root=abs(tmp_root(1))
                    if (abs(tmp_root(2)).le.min_root) min_root=abs(tmp_root(2))
                    
                        !-- Check whether the solution is inside the maximum scanning radius
                    if (min_root.le.max_radius) then
                        root(i,1)=min_root; root(i,2)=phi
                    end if

                else ! complex solutions
                    if (Guillem(6).eqv..TRUE.) then
                            !-- They are computed as such so they are not stored in memory and only affect the time because of the printing
                        write(*,'(" (α-FDB-",A2,A2,") - Discriminant of the second degree equation:",xF18.12)') axis_name(2-axis(1)),axis_name(2+axis(3)),discriminant
                        write(*,'(" (α-FDB-",A2,A2,") - Positive pure complex second degree root:",xF12.6," +i",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),-mmu+sqrt(mmu**2-4.0d0*aalpha*eenergy)/2.0d0/aalpha
                        write(*,'(" (α-FDB-",A2,A2,") - Negative pure complex second degree root:",xF12.6," +i",xF12.6," a.u")') axis_name(2-axis(1)),axis_name(2+axis(3)),-mmu-sqrt(mmu**2-4.0d0*aalpha*eenergy)/2.0d0/aalpha
                        write(*,*)
                    end if
                    continue
                end if

            end if
            if (Guillem(8).eqv..TRUE.) then
                write(*,'(" (α-FDB-",A2,A2,") - Cylindrical solutions (ρ,φ;",A2,") = (",F15.10,",",xF6.4,",",xF10.6,") a.u ")') axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(gpc),min_root,phi,F(gpc)
                write(*,'(" (α-FDB-",A2,A2,") - Cartesian solutions (",A2,",",A2,";",A2,") = (",xF12.6,",",xF12.6,",",xF12.6,") a.u ")') &
                & axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(gpc),min_root*F(2-axis(1)),min_root*F(2+axis(3)),F(gpc)
                write(*,*)
            end if
        end do

        !#################################################################################!
        !                  Bubble Sort algorithm for matrix root(i,j)                     !
        !#################################################################################!

            !-- Get the number of solutions (dimsort)
        dimsort=0
        do i=1,grid
            if (root(i,1).gt.0.0d0) then
                dimsort=dimsort+1
            end if 
        end do
        allocate(sort(dimsort,2))

            !-- Pass the solutions to the sorting matrix
        n=1
        do i=1,grid
            if (root(i,1).gt.0.0d0) then
                sort(n,1)=root(i,1); sort(n,2)=root(i,2)
                n=n+1
            end if
        end do

            !-- Bubble Sort the sorting matrix
        do i=dimsort,2,-1
            do j=1,i-1
                if(sort(j,1).gt.sort(j+1,1)) then
                    tmp_sort(1)=sort(j+1,1); tmp_sort(2)=sort(j+1,2)
                    sort(j+1,1)=sort(j,1); sort(j+1,2)=sort(j,2)
                    sort(j,1)=tmp_sort(1); sort(j,2)=tmp_sort(2)
                end if
            end do
        end do

            !-- Print the Bubble Sort matrix
        if (Guillem(9).eqv..TRUE.) then
            write(*,*)
            write(*,'(" (α-FDB-",A2,A2,") - Number of solutions:",xI2)') axis_name(2-axis(1)),axis_name(2+axis(3)),dimsort
            write(*,*) "----------------------------------------------------------------------------------------------" 
            do i=1,dimsort
                write(*,'(" (α-FDB-",A2,A2,") - Bubble sort matrix - Row ",I3," == (ρ,φ;",A2,") = (",xF12.8,",",xF6.4,",",xF6.4")")') axis_name(2-axis(1)),axis_name(2+axis(3)),i,axis_name(gpc),sort(i,1),sort(i,2),F(gpc)
            end do
            write(*,*) "----------------------------------------------------------------------------------------------" 
        end if

            !-- Print of the absolute minimum solution
        if (dimsort.ne.0) then
            write(*,*)
            write(*,'(" (α-FDB-",A2,A2,") - The minimum solution in cylindrical cooridnates: (ρ,φ;",A2,")   = (",xF12.8,",",xF6.4,",",xF6.4,") (*10^-4) a.u ")') &
            & axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(gpc),sort(1,1)*1.0d4,sort(1,2),F(gpc)*1.0d4
            write(*,'(" (α-FDB-",A2,A2,") - The minimum solution in Cartesian cooridnates:   (",A2,",",A2,";",A2,") = (",F12.6,",",F12.6,",",F12.6,") (*10^-4) a.u ")') &
            & axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(gpc),sort(1,1)*dcos(sort(1,2))*1.0d4,sort(1,1)*dsin(sort(1,2))*1.0d4,F(gpc)*1.0d4
            write(*,*)
        else
            write(*,'(" (α-FDB-",A2,A2,") - No solution within the maximum radius was found. Consider using a lower approximation. ")') axis_name(2-axis(1)),axis_name(2+axis(3))
            write(*,*)
        end if
        deallocate(sort)
        
    end if

!================================= α-VOAM-2D ==========================================
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!================================= β-VOAM-2D ==========================================

    if (order_nlop.eq.3) then

            ! XY=(1,1,0) gpc=3,axis(gpc)=0 --> Step=1 // XZ=(1,0,1) gpc=2,axis(gpc)=0 --> Step=2 // YZ=(0,1,1) gpc=1,axis(gpc)=0 --> Step=1
        step=nint(1+abs(cos(0.5d0*PI*gpc)))

            !-- Initialization of the tmp_NLOP values
        tmp_energy=E_0-abs(redox_coef)*redox_potential-target_barrier/hartree2kcal !-- Energy associated to thermochemistry

        do i=1,grid
            tmp_mu=0.0d0; tmp_alpha=0.0d0; tmp_beta=0.0d0; tmp_root=0.0d0

                !-- 0.le.φ.le.2π
            phi=((i-1)/(1.0d0*grid-1))*2.0d0*PI

                !-- Declaration of the field vector in cylindrical coordinates
            F(gpc)=initial_position(gpc);F(2-axis(1))=dcos(phi);F(2+axis(3))=dsin(phi)
                    !-- XY plane: F(1)=cos(phi); F(2)=sin(phi); F(3)=Fz
                    !-- XZ plane: F(1)=cos(phi); F(3)=sin(phi); F(2)=Fy
                    !-- YZ plane: F(2)=cos(phi); F(3)=sin(phi); F(1)=Fx

                !-- The maximum radius of scan is going to decrease with the increase of the constant field. The general expression if its decrement is the following
            max_radius=dsqrt(radius**2-F(gpc)**2)

                !-- Initialization of the calculation of the (constant-field dependent) energy and angular alpha
            tmp_energy=tmp_energy-mu(gpc)*F(gpc)-0.5d0*alpha(gpc,gpc)*F(gpc)**2-(1.0d0/6.0d0)*beta(gpc,gpc,gpc)*F(gpc)**3
            tmp_alpha=tmp_alpha+2.0d0*(alpha(2-axis(1),2+axis(3))+beta(2-axis(1),2+axis(3),gpc)*F(gpc))*sin(2.0d0*phi)

                !-- Computation of the angular NLOP and declaration of the parameters for the third degree equation
            do j=(2-axis(1)),(2+axis(3)),step
                    !-- Angular dipole moment
                tmp_mu=tmp_mu+F(j)*(mu(j)+alpha(j,gpc)*F(gpc)+beta(j,gpc,gpc)*F(gpc)**2)

                    !-- Angular polarizabilty matrix
                tmp_alpha=tmp_alpha+(alpha(j,j)+beta(j,j,gpc)*F(gpc))*F(gpc)**2

                    !-- Angular hyperpolarizability tensor
                tmp_beta=tmp_beta+beta(j,j,j)*F(j)**3+3.0d0*beta(2-axis(1),j,2+axis(1))*F(j)*dsin(2.0d0*phi)
            end do
            eenergy=cmplx(tmp_energy,0.0d0) 
            mmu=-cmplx(tmp_mu,0.0d0)
            aalpha=-0.5d0*cmplx(tmp_alpha,0.0d0)
            bbeta=-(1.0d0/6.0d0)*cmplx(tmp_beta,0.0d0)

                !-- Computation of the specific parameters for the third degree solution
            p=(3.0d0*bbeta*mmu-aalpha**2)/3.0d0/bbeta**2
            q=(2.0d0*aalpha**3-9.0d0*bbeta*aalpha*mmu+27.0d0*eenergy*bbeta**2)/27.0d0/bbeta**3
            omega=real(0.25d0*q**2+(1.0d0/27.0d0)*p**3)
            if (Guillem(4).eqv..TRUE.) then
                write(*,'(" Phi=    ",xF6.4)') phi
                write(*,'(" Energy coefficient:             ",xF18.12)') real(eenergy)
                write(*,'(" Dipole moment coefficient:      ",xF18.12)') real(mmu)
                write(*,'(" Polarizability coefficient:     ",xF18.12)') real(aalpha)
                write(*,'(" Hyperpolarizability coefficient:",xF18.12)') real(bbeta)
                write(*,'(" P=",xF18.12,", Q=",xF18.12," and ω=",xF18.12)') real(p),real(q),omega
                write(*,*)
            end if
           
            !##################################################################!
            !           Computation of the third degree solutions              !
            !##################################################################!
 
            if (omega.gt.tol) then ! ω > 0 --> There is a single real root
                tmp_root(1)=sign(abs(-0.5d0*q-sqrt(0.25d0*q**2+p**3/27.0d0))**(1.0d0/3.0d0),real(-0.5d0*q-sqrt(0.25d0*q**2+p**3/27.0d0)))
                tmp_root(1)=tmp_root(1)+sign(abs(-0.5d0*q+sqrt(0.25d0*q**2+p**3/27.0d0))**(1.0d0/3.0d0),real(-0.5d0*q+sqrt(0.25d0*q**2+p**3/27.0d0)))
                tmp_root(1)=tmp_root(1)-aalpha/3.0d0/bbeta
             
                    !-- Print the specific raw solutions
                if (Guillem(6).eqv..TRUE.) then
                    write(*,'(" Phi=    ",xF6.4)') phi
                    write(*,'(" (β-FDB-"A2,A2") - Unique specific solution:",xF18.6)') axis_name(2-axis(1)),axis_name(2+axis(3)),tmp_root(1)
                    write(*,*) "Complex solutions skipped."
                    write(*,*)
                end if

                    !-- Checking whether the solution is inside the radius of the sphere
                min_root=abs(tmp_root(1))
                if (min_root.le.max_radius) then
                    root(i,1)=tmp_root(1); root(i,2)=phi
                end if

            else if (abs(omega).lt.tol) then ! ω = 0 --> Three real solutions, two of them equal 
                tmp_root(1)=+2.0d0*sign(abs(-0.5d0*q)**(1.0d0/3.0d0),real(-0.5d0*q))-aalpha/3.0d0/bbeta
                tmp_root(2)=-0.5d0*tmp_root(1); tmp_root(3)=tmp_root(2)
                !tmp_root(2)=-1.0d0*sign(abs(-0.5d0*q)**(1.0d0/3.0d0),real(-0.5d0*q))-aalpha/3.0d0/bbeta

                    !-- Print the specific raw solutions
                if (Guillem(6).eqv..TRUE.) then
                    write(*,'(" Phi=    ",xF6.4)') phi
                    write(*,'(" (β-FDB-"A2,A2") - First specific solution: ",xF18.6)') axis_name(2-axis(1)),axis_name(2+axis(3)),tmp_root(1)
                    write(*,'(" (β-FDB-"A2,A2") - Second specific solution:",xF18.6)') axis_name(2-axis(1)),axis_name(2+axis(3)),tmp_root(2) 
                    write(*,'(" (β-FDB-"A2,A2") - Third specific solution: ",xF18.6)') axis_name(2-axis(1)),axis_name(2+axis(3)),tmp_root(3)
                    write(*,*)
                end if

                    !-- Get the closest-to-zero solution
                min_root=tmp_root(1)
                if (abs(tmp_root(2)).lt.abs(min_root)) min_root=abs(tmp_root(2))

                    !-- Checking whether the solutions are inside the sphere
                if (min_root.le.max_radius) then
                    root(i,1)=min_root; root(i,2)=phi
                end if

            else if (omega.lt.0.0d0) then ! ω < 0 --> Three different real solutions
                angle=dacos(-0.5d0*real(q)/sqrt(-1.0d0*real(p)**3/27.0d0))
                tmp_root(1)=2.0d0*sqrt(-p/3.0d0)*cos(angle/3.0d0)-aalpha/3.0d0/bbeta           
                tmp_root(2)=2.0d0*sqrt(-p/3.0d0)*cos((angle+2.0d0*PI)/3.0d0)-aalpha/3.0d0/bbeta
                tmp_root(3)=2.0d0*sqrt(-p/3.0d0)*cos((angle+4.0d0*PI)/3.0d0)-aalpha/3.0d0/bbeta

                    !-- Print the specific raw solutions
                if (Guillem(6).eqv..TRUE.) then
                    write(*,'(" Phi=    ",xF6.4)') phi
                    write(*,'(" (β-FDB-",A2,A2,") - First specific solution: ",xF18.6)') axis_name(2-axis(1)),axis_name(2+axis(3)),tmp_root(1)
                    write(*,'(" (β-FDB-",A2,A2,") - Second specific solution:",xF18.6)') axis_name(2-axis(2)),axis_name(2+axis(3)),tmp_root(2)
                    write(*,'(" (β-FDB-",A2,A2,") - Third specific solution: ",xF18.6)') axis_name(2-axis(3)),axis_name(2+axis(3)),tmp_root(3)
                    write(*,*)
                end if
                    
                    !-- Get the closest-to-zero solution
                min_root=tmp_root(1)
                do j=2,3
                    if (abs(tmp_root(j)).le.abs(min_root)) min_root=abs(tmp_root(j))
                end do
                
                    !-- Check whether there are solutions inside the scanning sphere or not
                if (min_root.lt.max_radius) then
                    root(i,1)=min_root; root(i,2)=phi
                end if 

            end if
            if (Guillem(8).eqv..TRUE.) then
                write(*,'(" (β-FDB-",A2,A2,") - Cylindrical solutions (ρ,φ;",A2,") = (",xF13.10,",",xF6.4,",",xF10.6,") a.u ")') axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(gpc),min_root,phi,F(gpc)
                write(*,'(" (β-FDB-",A2,A2,") - Cartesian solutions (",A2,",",A2,";",A2,") = (",F18.12,",",F18.12,",",F18.12,") a.u ")') &
                & axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(gpc),min_root*F(2-axis(1)),min_root*F(2+axis(3)),F(gpc)
                write(*,*)
            end if

        end do
    
        !#################################################################################!
        !                  Bubble Sort algorithm for matrix root(i,j)                     !
        !#################################################################################!

            !-- Get the number of solutions (dimsort)
        dimsort=0
        do i=1,grid
            if (root(i,1).gt.0.0d0) then
                dimsort=dimsort+1
            end if 
        end do
        allocate(sort(dimsort,2))

            !-- Pass the solutions to the sorting matrix
        n=1
        do i=1,grid
            if (root(i,1).gt.0.0d0) then
                sort(n,1)=root(i,1); sort(n,2)=root(i,2)
                n=n+1
            end if
        end do

            !-- Bubble Sort the sorting matrix
        do i=dimsort,2,-1
            do j=1,i-1
                if(sort(j,1).gt.sort(j+1,1)) then
                    tmp_sort(1)=sort(j+1,1); tmp_sort(2)=sort(j+1,2)
                    sort(j+1,1)=sort(j,1); sort(j+1,2)=sort(j,2)
                    sort(j,1)=tmp_sort(1); sort(j,2)=tmp_sort(2)
                end if
            end do
        end do

            !-- Print the Bubble Sort matrix
        if (Guillem(9).eqv..TRUE.) then
            write(*,*)
            write(*,'(" (β-FDB-",A2,A2,") - Number of solutions:",xI2)') axis_name(2-axis(1)),axis_name(2+axis(3)),dimsort
            write(*,*) "----------------------------------------------------------------------------------------------" 
            do i=1,dimsort
                write(*,'(" (β-FDB-",A2,A2,") - Bubble sort matrix - Row ",I2," == (ρ,φ;",A2,") = (",xF12.10,",",xF6.4,",",xF6.4")")') axis_name(2-axis(1)),axis_name(2+axis(3)),i,axis_name(gpc),sort(i,1),sort(i,2),F(gpc)
            end do
            write(*,*) "----------------------------------------------------------------------------------------------" 
        end if

            !-- Print of the absolute minimum solution
        if (dimsort.ne.0) then
            write(*,*)
            write(*,'(" (β-FDB-",A2,A2,") - The minimum solution in cylindrical cooridnates: (ρ,φ;",A2,")   = (",xF12.10,",",xF6.4,",",xF6.4,") (*10^-4) a.u ")') &
            & axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(gpc),sort(1,1)*1.0d4,sort(1,2),F(gpc)*1.0d4
            write(*,'(" (β-FDB-",A2,A2,") - The minimum solution in Cartesian cooridnates:   (",A2,",",A2,";",A2,") = (",F12.6,",",F12.6,",",F12.6,") (*10^-4) a.u ")') &
            & axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(2-axis(1)),axis_name(2+axis(3)),axis_name(gpc),sort(1,1)*dcos(sort(1,2))*1.0d4,sort(1,1)*dsin(sort(1,2))*1.0d4,F(gpc)*1.0d4
            write(*,*)
        else
            write(*,'(" (β-FDB-",A2,A2,") - No solution within the maximum radius was found. Consider using a lower approximation. ")') axis_name(2-axis(1)),axis_name(2+axis(3))
            write(*,*)
        end if
        deallocate(sort)

    end if
    
!================================= β-VOAM-2D ==========================================

    End subroutine VOAM_2D

!###################################################################################################
    !=========================================================================!
    !============================== 3D REGIME ================================!
    !=========================================================================!
!###################################################################################################

    Subroutine VOAM_3D(Guillem,order_nlop,gpc,radius,axis,initial_position,axis_name,grid)
    implicit none
    logical, dimension(10) :: Guillem
    integer, dimension(3) :: axis
    character*2,dimension(3) :: axis_name
    double precision, dimension(3) :: F,initial_position
    double precision, allocatable, dimension (:,:) :: dG_field
    integer order_nlop,grid
    double precision :: max_radius,radius
        !-- Only for μ-FDB
    double precision :: E_FF,numerator,denominator
    double precision :: modulus
        !-- For α/β-FDB 
    double precision, dimension (3) :: tmp_root,tmp_sort    ! Stores only the modulus
    double precision, dimension (grid,3) :: root            ! Modulus, theta and phi
    double precision, allocatable, dimension (:,:) :: sort  ! Bubble Sort matrix
    double complex :: eenergy,mmu,aalpha,bbeta
    double complex :: p,q
    double precision :: omega
    double precision :: tmp_energy,tmp_mu,tmp_alpha,tmp_beta,min_root
    double precision :: tmp_energy_lig,tmp_mu_lig,tmp_alpha_lig,tmp_beta_lig
    double precision :: theta,phi,angle,discriminant,field,orientation
    integer :: sign_aalpha,sign_eenergy
    integer :: dimsort,n
        !-- Integers for iteration
    integer gpc,step,upper,lower,tracker
    integer i,j,k

!================================ μ-VOAM-3D ==========================================    

    if (order_nlop.eq.1) then
        mu=mu+mu_lig(3)
        E_FF=E_0-abs(redox_coef)*redox_potential-(target_barrier/hartree2kcal)
        if (Guillem(5).eqv..TRUE.) then
            write(*,'(" (μ-FDB-XYZ) - The equation to be solved is: 0 =",xF7.3," +",xF7.3,xA2," +",xF7.3,xA2," +",xF7.3,xA2)') &
            & E_FF,mu(1),axis_name(1),mu(2),axis_name(2),mu(3),axis_name(3)
        end if

            !-- Lagrange formulae to compute the Fz component
        numerator=mu(3)*E_FF; denominator=mu(1)**2+mu(2)**2+mu(3)**2
        F(3)=numerator/denominator;F(1)=(mu(1)/mu(3))*F(3);F(2)=(mu(2)/mu(3))*F(3)
        modulus=dsqrt(F(1)**2+F(2)**2+F(3)**2)

            !-- Check whether the solution is inside the solving sphere
        if (modulus.lt.radius) then 
            write(*,'(" (μ-FDB-XYZ) - The minimum field strength is (",A2,",",A2,",",A2,") = (",F7.2,","F7.2,","F7.2,")(*10^-4) a.u")') &
            & axis_name(1),axis_name(2),axis_name(3),F(1)*1.0d4,F(2)*1.0d4,F(3)*1.0d4
        else
            if (Guillem(6).eqv..TRUE.) then
                write(*,'(" (μ-FDB-XYZ) - Specific solution is (",A2",",A2,",",A2,") = (",F18.12,F18.12,F18.12,") a.u")') axis_name(1),axis_name(2),axis_name(3),F(1),F(2),F(3)
                write(*,*)
            end if
            write(*,'(" (μ-FDB-XYZ) - The modulus of the vector is larger than the maximum input radius.")')
            write(*,'(" (μ-FDB-XYZ) - Consider using a larger radius or check the input file for any mistakes.")')
            write(*,*)
        end if

    end if

!================================ μ-VOAM-3D ==========================================    
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!================================= α-VOAM-3D ==========================================    

    if (order_nlop.eq.2) then
        
            !-- Initialization of the tmp_NLOP values
        tmp_energy=E_0-abs(redox_coef)*redox_potential-target_barrier/hartree2kcal !-- Energy associated to thermochemistry
        do i=1,3
            mu(i)=mu(i)+mu_lig(3)
            do j=1,3
                alpha(i,j)=alpha(i,j)+alpha_lig(i,j)
            end do
        end do
        
        do i=1,grid         !-- Loop for \theta  == π/{grid}
            do j=1,2*grid   !-- Loop for \varphi == 2π/{2*grid}
    
                    !-- Initialization of the tmp_NLOP values
                tmp_mu=0.0d0;tmp_alpha=0.0d0

                    !-- Declaration of the scanning angles: 0.le.φ.le.2π && 0.le.θ.le.π
                theta=((i-1)/(1.0d0*grid-1))*PI; phi=((j-1)/(2.0d0*grid-1))*2.0d0*PI

                    !-- Declaration of the field vectors
                F(1)=dsin(theta)*dcos(phi);F(2)=dsin(theta)*dsin(phi);F(3)=dcos(theta)
            
                    !-- Computation of the angular NLOP
                do k=1,3
                        !-- Angular dipole moment
                    tmp_mu=tmp_mu+mu(k)*F(k)
    
                        !-- Angular polarizability matrix
                    tmp_alpha=tmp_alpha+alpha(k,k)*F(k)**2+2.0d0*alpha(nint(1+0.2d0*k),nint(2+k/3.0d0))*F(nint(1+0.2d0*k))*F(nint(2+k/3.0d0))
                end do
                eenergy=cmplx(tmp_energy,0.0d0); mmu=-cmplx(tmp_mu,0.0d0); aalpha=-0.5d0*cmplx(tmp_alpha,0.0d0)
                if (Guillem(4).eqv..TRUE.) then
                    write(*,'(" Theta=    ",xF6.4x,"Phi=    ",xF6.4)') theta,phi
                    write(*,'(" Energy coefficient:        ",xF18.12)') real(eenergy)
                    write(*,'(" Dipole moment coefficient: ",xF18.12)') real(mmu)
                    write(*,'(" Polarizability coefficient:",xF18.12)') real(aalpha)
                    write(*,*)
                end if 
                 
                !##################################################################!
                !           Computing the second degree solutions                  !
                !##################################################################!

                if (abs(mmu).lt.tol) then ! angular-mmu is neglegible
                    sign_aalpha=sign(1.0d0,real(aalpha));sign_eenergy=sign(1.0d0,real(eenergy))

                    if(sign_aalpha.eq.sign_eenergy) then
                            ! The sign is the same ---> ax**2+b=0 ---> The solution is proven to be complex. Skipped
                        if(Guillem(6).eqv..TRUE.) then ! Print the specificsolutions
                                !-- Coded in this way so it is not stored in memory and is computed on-the-fly
                            write(*,'(" (α-FDB-XYZ) - Positive complex root:",xF12.6," +i",xF12.6," a.u")') +1.0d0*sqrt(-eenergy/aalpha)
                            write(*,'(" (α-FDB-XYZ) - Negative complex root:",xF12.6," +i",xF12.6," a.u")') -1.0d0*sqrt(-eenergy/aalpha)
                        end if
                        continue
                    else
                                !    ax**2-b=0 ---> The solution is proven to be real.
                            !    Double signed root == tmp_root(i)=positive_root  tmp_root(i+1)=negative_root
                        tmp_root(1)=+1.0d0*dsqrt(-real(eenergy)/real(aalpha)); tmp_root(2)=-tmp_root(1)
                        if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                            write(*,'(" (α-FDB-XYZ) - Neglegible dipole moment positive root:",xF18.6," a.u")') tmp_root(1)
                            write(*,'(" (α-FDB-XYZ) - Neglegible dipole moment negative root:",xF18.6," a.u")') tmp_root(2)
                            write(*,*)
                        end if

                            !-- Check whether the solution is inside the maximum scanning radius
                        min_root=abs(tmp_root(1))
                        if (min_root.lt.max_radius) then
                            root(i,1)=min_root; root(i,2)=theta; root(i,3)=phi
                        end if
                    end if

                else if (abs(aalpha).lt.tol) then ! angular-aalpha is neglegible --> ax**1+b=0
                    tmp_root(1)=-real(eenergy)/real(mmu)

                    if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                        write(*,'(" (α-FDB-XYZ) - Neglegible alpha root:",xF12.6," a.u")') tmp_root(1)
                    end if
                        
                        !-- Check whether the solution is inside the maximum scanning radius
                    min_root=abs(tmp_root(1))
                    if (max_radius.gt.abs(tmp_root(1))) then
                        root(i,1)=min_root; root(i,2)=theta; root(i,3)=phi
                    end if

                else ! complete second degree equation
                    discriminant=mmu**2.0d0-4.0d0*aalpha*eenergy

                    if (discriminant.ge.0.0d0) then
                        tmp_root(1)=-mmu+dsqrt(discriminant)/2.0d0/aalpha
                        tmp_root(2)=-mmu-dsqrt(discriminant)/2.0d0/aalpha
                        if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                            write(*,'(" (α-FDB-XYZ) - Discriminant of the second degree equation:",xF18.12)') discriminant
                            write(*,'(" (α-FDB-XYZ) - Positive pure second degree root:",xF12.6," a.u")') tmp_root(1)
                            write(*,'(" (α-FDB-XYZ) - Negative pure second degree root:",xF12.6," a.u")') tmp_root(2)
                            write(*,*)
                        end if

                            !-- Get the closest-to-zero solution; tmp_root(1) and tmp_root(2) are not necessary equal; only when the discriminant is zero
                        min_root=abs(tmp_root(1))
                        if (abs(tmp_root(2)).lt.min_root) min_root=abs(tmp_root(2))

                            !-- Check whether the solution is inside the maximum scanning radius
                        if (max_radius.gt.min_root) then
                            root(i,1)=min_root; root(i,2)=theta; root(i,3)=phi
                        end if

                    else ! complex solutions
                        if (Guillem(6).eqv..TRUE.) then
                                !-- They are computed as such so they are not stored in memory and only affect the time because of the printing
                            write(*,'(" (α-FDB-XYZ) - Discriminant of the second degree equation:",xF18.12)') discriminant
                            write(*,'(" (α-FDB-XYZ) - Positive pure complex second degree root:",xF12.6," +i",xF12.6," a.u")') -mmu+sqrt(mmu**2-4.0d0*aalpha*eenergy)/2.0d0/aalpha
                            write(*,'(" (α-FDB-XYZ) - Negative pure complex second degree root:",xF12.6," +i",xF12.6," a.u")') -mmu-sqrt(mmu**2-4.0d0*aalpha*eenergy)/2.0d0/aalpha
                            write(*,*)
                        end if
                        continue
                    end if

                end if
                
                    !-- Print all min_root's in both spherical and Cartesian coordinates
                if (Guillem(8).eqv..TRUE.) then
                    write(*,'(" (α-FDB-XYZ) - Polar spherical solutions (R,θ,φ)  = (",xF12.10,",",xF6.4,",",xF6.4,") a.u ")') min_root,theta,phi
                    write(*,'(" (α-FDB-XYZ) - Cartesian solutions (Fx,Fy,Fz) = (",F12.6,",",F12.6,",",F12.6,") a.u ")') min_root*F(1),min_root*F(2),min_root*F(3)
                    write(*,*)
                end if

            end do
        end do

        !#################################################################################!
        !                  Bubble Sort algorithm for matrix root(i,j)                     !
        !#################################################################################!

            !-- Get the number of solutions (dimsort)
        dimsort=0
        do i=1,grid
            if (root(i,1).gt.0.0d0) then
                dimsort=dimsort+1
            end if 
        end do
        allocate(sort(dimsort,3))

            !-- Pass the solutions to the sorting matrix
        n=1
        do i=1,grid
            if (root(i,1).gt.0.0d0) then
                sort(n,1)=root(i,1); sort(n,2)=root(i,2); sort(n,3)=root(i,3)
                n=n+1
            end if
        end do

            !-- Bubble Sort the sorting matrix
        do i=dimsort,2,-1
            do j=1,i-1
                if(sort(j,1).gt.sort(j+1,1)) then
                    tmp_sort(1)=sort(j+1,1); tmp_sort(2)=sort(j+1,2); tmp_sort(3)=sort(j+1,3)
                    sort(j+1,1)=sort(j,1);   sort(j+1,2)=sort(j,2);   sort(j+1,3)=sort(j,3)
                    sort(j,1)=tmp_sort(1);   sort(j,2)=tmp_sort(2);   sort(j,3)=tmp_sort(3)
                end if
            end do
        end do

            !-- Print the Bubble Sort matrix
        if (Guillem(9).eqv..TRUE.) then
            write(*,*)
            write(*,'(" (α-FDB-XYZ) - Number of solutions:",xI2)') dimsort
            write(*,*) "----------------------------------------------------------------------------------------------" 
            do i=1,dimsort
                write(*,'(" (α-FDB-XYZ) - Bubble sort matrix - Row ",I2," == (R,θ,φ) = (",xF12.10,",",xF6.4,",",xF6.4")")') i,sort(i,1),sort(i,2),sort(i,3)
            end do
            write(*,*) "----------------------------------------------------------------------------------------------" 
        end if

            !-- Print of the absolute minimum solution
        if (dimsort.ne.0) then
            write(*,*)
            write(*,'(" (α-FDB-XYZ) - The minimum solution in polar spherical cooridnates: (R,θ,φ)    = (",xF12.10,",",xF6.4,",",xF6.4,") (*10^-4) a.u ")') &
            & sort(1,1)*1.0d4,sort(1,2),sort(1,3)
            write(*,'(" (α-FDB-XYZ) - The minimum solution in Cartesian cooridnates: (Fx,Fy,Fz) = (",F12.6,",",F12.6,",",F12.6,") (*10^-4) a.u ")') &
            & sort(1,1)*dsin(sort(1,2))*dcos(sort(1,3))*1.0d4,sort(1,1)*dsin(sort(1,2))*dsin(sort(1,3))*1.0d4,sort(1,1)*dcos(sort(1,2))*1.0d4
            write(*,*)
        else
            write(*,'(" (α-FDB-XYZ) - No solution within the maximum radius was found. Consider using a lower approximation. ")') 
            write(*,*)
        end if
        deallocate(sort)

    end if

!================================= α-VOAM-3D ==========================================    
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!################################################################################################################################
!================================= β-VOAM-3D ==========================================    

    if (order_nlop.eq.3) then
    
            !-- Initialization of the tmp_NLOP values
        tmp_energy=E_0-abs(redox_coef)*redox_potential !-- Energy associated to thermochemistry

        do i=1,grid         !-- Loop for \theta  == π/{grid}
             do j=1,2*grid   !-- Loop for \varphi == 2π/{2*grid}

                    !-- Declaration of the scanning angles: 0.le.φ.le.2π && 0.le.θ.le.π
                theta=((i-1)/(1.0d0*grid-1))*PI; phi=((j-1)/(2.0d0*grid-1))*2.0d0*PI

                    !-- Declaration of the field vectors:
                F(1)=dsin(theta)*dcos(phi); F(2)=dsin(theta)*dsin(phi); F(3)=dcos(theta)

                    !-- Initialization of the tmp_NLOP values
                tmp_mu=0.0d0; tmp_alpha=0.0d0; tmp_beta=3.0d0*beta(1,2,3)*F(1)*F(2)*F(3)
                !tmp_mu_lig=0.0d0; tmp_alpha_lig=0.0d0; Tmp_beta_lig=3.0d0*beta_lig(1,2,3)*F(1)*F(2)*F(3)
                tmp_mu_lig=mu_lig(3);tmp_alpha_lig=alpha_lig(3,3);tmp_beta_lig=beta_lig(3,3,3)
   
                    !-- Computation of the angular NLOP
                do k=1,3
                        !-- Angular dipole moment coefficient
                    tmp_mu=tmp_mu+mu(k)*F(k)
                    !tmp_mu_lig=tmp_mu+mu(3)*F(k)

                        !-- Angular polarizability coefficient
                    tmp_alpha=tmp_alpha+alpha(k,k)*F(k)**2+alpha(nint(1+0.2d0*k),nint(2+k/3.0d0))*F(nint(1+0.2d0*k))*F(nint(2+k/3.0d0))
                    !tmp_alpha_lig=tmp_alpha+alpha_lig(k,k)*F(k)**2+alpha_lig(nint(1+0.2d0*k),nint(2+k/3.0d0))*F(nint(1+0.2d0*k))*F(nint(2+k/3.0d0))

                        !-- Angular hyperpolarizability coefficient - Part 1
                    tmp_beta=tmp_beta+beta(k,k,k)*F(k)**3+3.0d0*beta(1,k,3)*F(1)*F(k)*F(3)
                    !tmp_beta_lig=tmp_beta_lig+beta_lig(k,k,k)*F(k)**3+3.0d0*beta_lig(1,k,3)*F(1)*F(k)*F(3)
                end do
                do k=1,2
                        !-- Angular hyperpolarizability coefficient - Part 2
                    tmp_beta=tmp_beta+3.0d0*(beta(k,k,k+1)*F(k)+beta(k,k+1,k+1)*F(k+1))*F(k)*F(k+1)
                    !tmp_beta_lig=tmp_beta_lig+3.0d0*(beta_lig(k,k,k+1)*F(k)+beta_lig(k,k+1,k+1)*F(k+1))*F(k)*F(k+1)
                end do

                if(allocated(dG_field)) deallocate (dG_field)
                allocate(dG_field(21,2))
                do k=1,21
                        !-- Computing an internal VOAM_0D to know the shape of the radial 3D-ΔG(F) function
                    field=-radius+(k-1)*1.0d-3
                    dG_field(k,2)=tmp_energy-(tmp_mu+sign(1.0d0,field)*tmp_mu_lig)*field-0.5d0*(tmp_alpha+tmp_alpha_lig)*field**2-(1.0d0/6.0d0)*(tmp_beta+sign(1.0d0,field)*tmp_beta_lig)*field**3
                    dG_field(k,1)=field
                end do
                tracker=0
                do k=1,21
                        !-- Defining the region where the solution should be with a 3 kcal/mol threshold
                    if (abs(dG_field(k,2)-target_barrier/hartree2kcal)*hartree2kcal.le.3.0d0) then
                        field=dG_field(k,1)
                        if (abs(field).le.abs(dG_field(k,1))) then
                            field=dG_field(k,1); orientation=1.0d0*sign(1.0d0,field)
                            exit
                        end if
                    else
                        tracker=tracker+1
                        if (tracker.eq.21) goto 9999
                    end if
                end do
                deallocate(dG_field)
                
                eenergy=cmplx(tmp_energy-target_barrier/hartree2kcal,0.0d0)
                mmu=-cmplx(tmp_mu+orientation*tmp_mu_lig,0.0d0)
                aalpha=-0.5d0*cmplx(tmp_alpha+tmp_alpha_lig,0.0d0)
                bbeta=-(1.0d0/6.0d0)*cmplx(tmp_beta+orientation*tmp_beta_lig,0.0d0)

                    !-- Computation of the specific parameters for the third degree solution
                p=(3.0d0*bbeta*mmu-aalpha**2)/3.0d0/bbeta**2
                q=(2.0d0*aalpha**3-9.0d0*bbeta*aalpha*mmu+27.0d0*eenergy*bbeta**2)/27.0d0/bbeta**3
                omega=real(0.25d0*q**2+(1.0d0/27.0d0)*p**3)
                if (Guillem(4).eqv..TRUE.) then
                    write(*,'(" Theta=    ",xF6.4," Phi=    ",xF6.4)') theta,phi
                    write(*,'(" Energy coefficient:             ",xF18.12)') real(eenergy)
                    write(*,'(" Dipole moment coefficient:      ",xF18.12)') real(mmu)
                    write(*,'(" Polarizability coefficient:     ",xF18.12)') real(aalpha)
                    write(*,'(" Hyperpolarizability coefficient:",xF18.12)') real(bbeta)
                    write(*,'(" P=",xF18.12,", Q=",xF18.12," and ω=",xF18.12)') real(p),real(q),omega
                    write(*,*)
                end if

                    !##################################################################!
                    !              Computing the third degree solutions                !
                    !##################################################################!
                
                if (abs(bbeta).lt.tol) then !angular-bbbeta is neglegible
                    if (abs(mmu).lt.tol) then ! angular-mmu and -beta are neglegible --> ax**2+-b=0
                        sign_aalpha=sign(1.0d0,real(aalpha));sign_eenergy=sign(1.0d0,real(eenergy))

                        if(sign_aalpha.eq.sign_eenergy) then
                                ! The sign is the same ---> ax**2+b=0 ---> The solution is proven to be complex. Skipped
                            if(Guillem(6).eqv..TRUE.) then ! Print the specificsolutions
                                    !-- Coded in this way so it is not stored in memory and is computed on-the-fly
                                write(*,'(" (β-FDB-XYZ) - Positive complex root:",xF12.6," +i",xF12.6," a.u")') +1.0d0*sqrt(-eenergy/aalpha)
                                write(*,'(" (β-FDB-XYZ) - Negative complex root:",xF12.6," +i",xF12.6," a.u")') -1.0d0*sqrt(-eenergy/aalpha)
                                write(*,*)
                            end if
                        else
                                    !    ax**2-b=0 ---> The solution is proven to be real.
                                !    Double signed root == tmp_root(i)=positive_root  tmp_root(i+1)=negative_root
                            tmp_root(1)=+1.0d0*dsqrt(-real(eenergy)/real(aalpha)); tmp_root(2)=-tmp_root(1)
                            if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                                write(*,'(" (β-FDB-XYZ) - Neglegible dipole moment positive root:",xF18.6," a.u")') tmp_root(1)
                                write(*,'(" (β-FDB-XYZ) - Neglegible dipole moment negative root:",xF18.6," a.u")') tmp_root(2)
                                write(*,*)
                            end if

                                !-- Check whether the solution is inside the maximum scanning radius
                            min_root=abs(tmp_root(1))
                            if (min_root.le.radius) then
                                root(i,1)=min_root; root(i,2)=theta; root(i,3)=phi
                            end if
                        end if

                    else if (abs(aalpha).lt.tol) then ! angular-aalpha and -bbeta are neglegible --> ax**1+b=0
                        tmp_root(1)=-real(eenergy)/real(mmu)

                        if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                            write(*,'(" (β-FDB-XYZ) - Neglegible alpha root:",xF12.6," a.u")') tmp_root(1)
                        end if
                            
                            !-- Check whether the solution is inside the maximum scanning radius
                        min_root=abs(tmp_root(1))
                        if (min_root.le.radius) then
                            root(i,1)=min_root; root(i,2)=theta; root(i,3)=phi
                        end if

                    else ! complete second degree equation
                        discriminant=mmu**2.0d0-4.0d0*aalpha*eenergy

                        if (discriminant.ge.0.0d0) then
                            tmp_root(1)=-mmu+dsqrt(discriminant)/2.0d0/aalpha
                            tmp_root(2)=-mmu-dsqrt(discriminant)/2.0d0/aalpha
                            if (Guillem(6).eqv..TRUE.) then !-- Print the specific solutions
                                write(*,'(" (β-FDB-XYZ) - Discriminant of the second degree equation:",xF18.12)') discriminant
                                write(*,'(" (β-FDB-XYZ) - Positive pure second degree root:",xF12.6," a.u")') tmp_root(1)
                                write(*,'(" (β-FDB-XYZ) - Negative pure second degree root:",xF12.6," a.u")') tmp_root(2)
                                write(*,*)
                            end if

                                !-- Get the closest-to-zero solution; tmp_root(1) and tmp_root(2) are not necessary equal; only when the discriminant is zero
                            min_root=abs(tmp_root(1))
                            if (abs(tmp_root(2)).lt.min_root) min_root=abs(tmp_root(2))

                                !-- Check whether the solution is inside the maximum scanning radius
                            if (min_root.le.radius) then
                                root(i,1)=min_root; root(i,2)=theta; root(i,3)=phi
                            end if

                        else ! complex solutions
                            if (Guillem(6).eqv..TRUE.) then
                                    !-- They are computed as such so they are not stored in memory and only affect the time because of the printing
                                write(*,'(" (β-FDB-XYZ) - Discriminant of the second degree equation:",xF18.12)') discriminant
                                write(*,'(" (β-FDB-XYZ) - Positive pure complex second degree root:",xF12.6," +i",xF12.6," a.u")') -mmu+sqrt(mmu**2-4.0d0*aalpha*eenergy)/2.0d0/aalpha
                                write(*,'(" (β-FDB-XYZ) - Negative pure complex second degree root:",xF12.6," +i",xF12.6," a.u")') -mmu-sqrt(mmu**2-4.0d0*aalpha*eenergy)/2.0d0/aalpha
                                write(*,*)
                            end if
                        end if
                    end if
                    
                        !-- Print all min_root's in both spherical and Cartesian coordinates
                    if (Guillem(8).eqv..TRUE.) then
                        write(*,'(" (β-FDB-XYZ) - Spherical polar solutions (R,θ,φ)  = (",xF12.10,",",xF6.4,",",xF6.4,") a.u ")') min_root,theta,phi
                        write(*,'(" (β-FDB-XYZ) - Cartesian solutions (Fx,Fy,Fz) = (",F12.6,",",F12.6,",",F12.6,") a.u ")') min_root*F(1),min_root*F(2),min_root*F(3)
                        write(*,*)
                    end if

                else !-- Angular-bbeta is not neglegible
                    if (omega.gt.tol) then ! ω > 0 --> There is a single real root
                        tmp_root(1)=sign(abs(-0.5d0*q-sqrt(0.25d0*q**2+p**3/27.0d0))**(1.0d0/3.0d0),real(-0.5d0*q-sqrt(0.25d0*q**2+p**3/27.0d0)))
                        tmp_root(1)=tmp_root(1)+sign(abs(-0.5d0*q+sqrt(0.25d0*q**2+p**3/27.0d0))**(1.0d0/3.0d0),real(-0.5d0*q+sqrt(0.25d0*q**2+p**3/27.0d0)))
                        tmp_root(1)=tmp_root(1)-aalpha/3.0d0/bbeta
                     
                            !-- Print the specific raw solutions
                        if (Guillem(6).eqv..TRUE.) then
                            write(*,'(" Theta=    ",xF6.4," Phi=    ",xF6.4)') theta,phi
                            write(*,'(" (β-FDB-XYZ) - Unique specific solution:",xF18.6)') tmp_root(1)
                            write(*,*) "Complex solutions skipped."
                            write(*,*)
                        end if

                            !-- Checking whether the solution is inside the radius of the sphere
                        min_root=abs(tmp_root(1))
                        if (min_root.le.radius) then
                            root(i,1)=min_root; root(i,2)=theta; root(i,3)=phi
                        end if

                    else if (abs(omega).lt.tol) then ! ω = 0 --> Three real solutions, two of them equal 
                        tmp_root(1)=+2.0d0*sign(abs(-0.5d0*q)**(1.0d0/3.0d0),real(-0.5d0*q))-aalpha/3.0d0/bbeta
                        tmp_root(2)=-0.5d0*tmp_root(1); tmp_root(3)=tmp_root(2)
                        !tmp_root(2)=-1.0d0*sign(abs(-0.5d0*q)**(1.0d0/3.0d0),real(-0.5d0*q))-aalpha/3.0d0/bbeta

                            !-- Print the specific raw solutions
                        if (Guillem(6).eqv..TRUE.) then
                            write(*,'(" Theta=    ",xF6.4," Phi=    ")') theta,phi
                            write(*,'(" (β-FDB-XYZ) - First specific solution: ",xF18.6)') tmp_root(1)
                            write(*,'(" (β-FDB-XYZ) - Second specific solution:",xF18.6)') tmp_root(2) 
                            write(*,'(" (β-FDB-XYZ) - Third specific solution: ",xF18.6)') tmp_root(3)
                            write(*,*)
                        end if

                            !-- Get the closest-to-zero solution
                        min_root=abs(tmp_root(1))
                        if (abs(tmp_root(2)).lt.abs(min_root)) min_root=abs(tmp_root(2))

                            !-- Checking whether the solutions are inside the sphere
                        if (min_root.le.radius) then
                            root(i,1)=min_root; root(i,2)=theta; root(i,3)=phi
                        end if

                    else if (omega.lt.0.0d0) then ! ω < 0 --> Three different real solutions
                        angle=dacos(-0.5d0*real(q)/sqrt(-1.0d0*real(p)**3/27.0d0))
                        tmp_root(1)=2.0d0*sqrt(-p/3.0d0)*cos(angle/3.0d0)-aalpha/3.0d0/bbeta     
                        tmp_root(2)=2.0d0*sqrt(-p/3.0d0)*cos((angle+2.0d0*PI)/3.0d0)-aalpha/3.0d0/bbeta
                        tmp_root(3)=2.0d0*sqrt(-p/3.0d0)*cos((angle+4.0d0*PI)/3.0d0)-aalpha/3.0d0/bbeta

                            !-- Print the specific raw solutions
                        if (Guillem(6).eqv..TRUE.) then
                            write(*,'(" Theta=    ",xF6.4," Phi=    ",xF6.4)') phi
                            write(*,'(" (β-FDB-XYZ) - First specific solution: ",xF18.6)') tmp_root(1)
                            write(*,'(" (β-FDB-XYZ) - Second specific solution:",xF18.6)') tmp_root(2)
                            write(*,'(" (β-FDB-XYZ) - Third specific solution: ",xF18.6)') tmp_root(3)
                            write(*,*)
                        end if
                            
                            !-- Get the closest-to-zero solution
                        min_root=tmp_root(1)
                        do k=2,3
                            if (abs(tmp_root(k)).le.abs(min_root)) min_root=abs(tmp_root(k))
                        end do
                        
                            !-- Check whether there are solutions inside the scanning sphere or not
                        if (min_root.le.radius) then
                            root(i,1)=min_root; root(i,2)=theta; root(i,3)=phi
                        end if 

                    end if
                        
                        !-- Print the minimum solutions, though they might not be inside the solving sphere
                    if (Guillem(8).eqv..TRUE.) then
                        write(*,'(" (β-FDB-XYZ) - Polar spherical solutions (R,θ,φ) = (",xF13.10,",",xF6.4,",",xF10.6,") a.u ")') min_root,theta,phi
                        write(*,'(" (β-FDB-XYZ) - Cartesian solutions (Fx,Fy,Fz) = (",F18.12,",",F18.12,",",F18.12,") a.u ")') &
                        & min_root*F(1),min_root*F(2),min_root*F(3)
                        write(*,*)
                    end if

                end if
                9999 continue
                tmp_mu=0.0d0;tmp_alpha=0.0d0;tmp_beta=0.0d0;tmp_mu_lig=0.0d0;tmp_alpha_lig=0.0d0;tmp_beta_lig=0.0d0
            end do
        end do

        !#################################################################################!
        !                  Bubble Sort algorithm for matrix root(i,j)                     !
        !#################################################################################!

            !-- Get the number of solutions (dimsort)
        dimsort=0
        do i=1,grid
            if (root(i,1).gt.0.0d0) then
                dimsort=dimsort+1
            end if 
        end do
        allocate(sort(dimsort,3))

            !-- Pass the solutions to the sorting matrix
        n=1
        do i=1,grid
            if (root(i,1).gt.0.0d0) then
                sort(n,1)=root(i,1); sort(n,2)=root(i,2); sort(n,3)=root(i,3)
                n=n+1
            end if
        end do

            !-- Bubble Sort the sorting matrix
        do i=dimsort,2,-1
            do j=1,i-1
                if(sort(j,1).gt.sort(j+1,1)) then
                    tmp_sort(1)=sort(j+1,1); tmp_sort(2)=sort(j+1,2); tmp_sort(3)=sort(j+1,3)
                    sort(j+1,1)=sort(j,1);   sort(j+1,2)=sort(j,2);   sort(j+1,3)=sort(j,3)
                    sort(j,1)=tmp_sort(1);   sort(j,2)=tmp_sort(2);   sort(j,3)=tmp_sort(3)
                end if
            end do
        end do

            !-- Print the Bubble Sort matrix
        if (Guillem(9).eqv..TRUE.) then
            write(*,*)
            write(*,'(" (β-FDB-XYZ) - Number of solutions:",xI3)') dimsort
            write(*,*) "----------------------------------------------------------------------------------------------" 
            do i=1,dimsort
                write(*,'(" (β-FDB-XYZ) - Bubble sort matrix - Row ",I3," == (R,θ,φ) = (",xF12.8,",",xF6.4,",",xF6.4")")') i,sort(i,1)*1.0d4,sort(i,2),sort(i,3)
            end do
            write(*,*) "----------------------------------------------------------------------------------------------" 
        end if

            !-- Print of the absolute minimum solution
        if (dimsort.ne.0) then
            write(*,*)
            write(*,'(" (β-FDB-XYZ) - The minimum solution in polar spherical cooridnates: (R,θ,φ)    = (",xF12.8,",",xF6.4,",",xF6.4,") (*10^-4) a.u ")') &
            & sort(1,1)*1.0d4,sort(1,2),sort(1,3)
            write(*,'(" (β-FDB-XYZ) - The minimum solution in Cartesian cooridnates: (Fx,Fy,Fz) = (",F12.6,",",F12.6,",",F12.6,") (*10^-4) a.u ")') &
            & sort(1,1)*dsin(sort(1,2))*dcos(sort(1,3))*1.0d4,sort(1,1)*dsin(sort(1,2))*dsin(sort(1,3))*1.0d4,sort(1,1)*dcos(sort(1,2))*1.0d4
            write(*,*)
        else if(dimsort.eq.0.or.sort(1,1).gt.radius) then
            write(*,'(" (β-FDB-XYZ) - No solution within the maximum radius was found. Consider using a lower approximation. ")') 
            write(*,*)
        end if
        deallocate(sort)

    end if

!================================= β-VOAM-3D ==========================================    
    End subroutine VOAM_3D

End module VOAM

!##############################################################################!
!==============================================================================!
!=======================  MAIN VOAM STARTS HERE  ==============================!
!==============================================================================!
!##############################################################################!

Program mainVOAM
use VOAM
!use fieldvector
implicit none
character*10 ET,axis_scan,approximation,stoichiometry,mister,field_domain
character*80 file ! Name of the file
character*100 line,title ! Characters to read the input file
character*80, dimension(10) :: name_reactant_reactant
character*80, dimension(10) :: name_product_product
character*80, dimension(10) :: name_reactant_ligand
character*80, dimension(10) :: name_product_ligand
!---- Algebraic elements for the NLOP reading
double precision :: E_iter
double precision, dimension(3) :: mu_iter
double precision, dimension(3,3) :: alpha_iter
double precision, dimension(3,3,3) :: beta_iter
!---- Integers for the input reading
integer :: index_NLOP,index_scan,index_field !-- Indices  for the method line
integer :: index_barrier,index_stoi,index_redox,index_potential !-- Indices for the thermochemistry line
integer :: index_initial,index_initial_end,index_modulus,index_trust,index_grid
integer :: index_tol
!---- Integers for the NLOP reading
!integer :: index_reactant_reactant,index_product_product
!integer :: index_ligand_reactant,index_ligand_product
integer :: n,m,step
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

!==============================================================================!
!               Conditional statements for Guillem's purposes                  !
    Guillem(1)=print_NLOP_react;Guillem(2)=print_NLOP_prod
    Guillem(3)=print_NLOP_lig;  Guillem(4)=print_coeffs;
    Guillem(5)=print_EQ;        Guillem(6)=print_sols; 
    Guillem(7)=print_checks;    Guillem(8)=print_mins;
    Guillem(9)=print_bubble;    !Guillem(10)=
!                                                                              !
!==============================================================================!

    !-- Declaration of a character vector for the later output handling
        axis_name(1)="Fx";axis_name(2)="Fy";axis_name(3)="Fz"    

!================================KEYWORDS BLOCK================================!
!                                                                              !
!-----------------------Basic/Reaction settings of the program-----------------!
read(2,'(A80)') title
read(2,'(A80)') line !-- Line before the "route section"

    !======================================================================!
                        !   Reading the method line   !
    !======================================================================!

read(2,'(A100)') line !-- Line corresponding to the 'Method' line
if(index(line,"Method:").ne.0.or.index(line,"method:").ne.0) then ! Try reading the axis scan and the NLOPs

        !--- Indices for the method-related keywords
    index_scan=index(line,"scan")                   !-- Detecting the first coincidence for the 'scan' word to detect the scanning axis
    index_NLOP=index(line,"aylor")                  !-- Same as index_scan
    index_NLOP=index_NLOP+index(line,"AYLOR")

    index_Guillem=index(line,"Guillem")             ! It's me, hi! I'm the developer, it's me
    index_Guillem=index_Guillem+index(line,"SUDO")
    index_Guillem=index_Guillem+index(line,"sudo")

        !-- Reading axis scan for the minimization
    if (index(line,"Axis scan=").ne.0.or.index(line,"axis scan=").ne.0) then
        read(line(index_scan+5:index_scan+8),*) axis_scan
        select case(axis_scan)
            case ("Scan","scan")
                n_dim=0
                !====================== 1D CASES ==========================!
            case ("00X","0X0","X00","0x0","x00","00x","X","x","XXX","xxx","xx","XX")
                        !-- Solve for R but with X-oriented NLOP
                    !-- Theta is PI halves and phi is zero, unless stated in the initial point
                n_dim=1;gpc=1;axis(gpc)=1
            case ("00Y","0Y0","Y00","0y0","y00","00y","Y","y","YYY","yyy","yy","YY")
                        !-- Solve for R but with Y-oriented NLOP
                    !-- Theta and phi are PI halves, unless stated in the initial point
                n_dim=1;gpc=2;axis(gpc)=1
            case ("00Z","0Z0","Z00","0z0","z00","00z","Z","z","ZZZ","zzz","zz","ZZ")
                        !-- Solve for R but with Z-oriented NLOP
                    !-- Theta is zero, unless stated in the initial point
                n_dim=1;gpc=3;axis(gpc)=1
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
        write(*,*) "Bad input scan keyword! STOP!"
        stop
    end if
        
        !-- Reading the order of NLOP
    if (index(line,"Taylor=").ne.0.or.index(line,"taylor=").ne.0.or.index(line,"TAYLOR=").ne.0) then
        read(line(index_NLOP+6:index_NLOP+12),*) approximation
        select case(approximation)
            case ("Dipole","dipole","mu","Mu")
                order_nlop=1
            case ("Alpha","alpha")
                order_nlop=2
            case ("Beta","beta")
                order_nlop=3
            case ("Gamma","gamma")
                order_nlop=4
            case default
                write(*,*) "Bad input approximation! STOP!"
        end select
    else
        write(*,*) "Bad input NLOP keyword! STOP!"
        stop
    end if

        !-- Reading the SUDO keywords
    if (index(line,"SUDO=").ne.0.or.index(line,"sudo=").ne.0) then ! Guillem is here
        read(line(index_Guillem+5:index_Guillem+25),*) mister
        select case (mister)
            case ("Guillem")
                Guillem=.TRUE.
            case ("Reactant","reactant","Reactants","reactants","Reacts","reacts","react","React")
                Guillem(1)=.TRUE.
            case ("Product","product","Products","products","prod","prods","Prod","Prods")
                Guillem(2)=.TRUE.
            case ("Ligand","ligand","Ligands","ligands","Lig","lig","Ligs","ligs")
                Guillem(3)=.TRUE.
            case ("Coefficient","coefficient","coeff","coeffs")
                Guillem(4)=.TRUE.
            case ("Equations","equation","equations","Equation")
                Guillem(5)=.TRUE.
            case ("Solutions","solutions","solution","Solution","sols","Sols","SOLS")
                Guillem(6)=.TRUE.
            case ("Check","check")
                Guillem(7)=.TRUE.
            case ("min_sols","Min sols","min sols","Min_sols")
                Guillem(8)=.TRUE.
            case ("sorting","bubble sort","bubblesort","Bubblesort","Bubble Sort","Bubble","bubble")
                Guillem(9)=.TRUE.
            case ("NLOP","nlop")
                Guillem(1)=.TRUE.; Guillem(2)=.TRUE.; Guillem(3)=.TRUE.
            case ("full check","Full check","fullcheck","All check","allcheck","all check")
                Guillem(4)=.TRUE.; Guillem(6)=.TRUE.
            case ("multicheck","Multicheck")
                Guillem(8)=.TRUE.; Guillem(9)=.TRUE.
        end select
    else ! No sudo
        continue
    end if

else
    write(*,*) "Error reading the method line! Stop!"
    stop
end if

!------------------------------------------------------------------------------!
!##############################################################################!
!------------------------------------------------------------------------------!

    !=====================================================================!
                    !   Reading the thermochemistry line   !
    !=====================================================================!

read(2,'(A100)') line !-- Line corresponding to the 'Thermochemistry' line
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
        write(*,*) "Redox coeficient not found. Set to 0"
        redox_coef=0
    end if
    
    if(index(line,"Potential=").ne.0.or.index(line,"potential=").ne.0) then
        read(line(index_potential+8:index_potential+13),*) redox_potential ! vs SHE !-- Reformular segons redox_coeff
                    !-- S'ha de mirar si s'ha de corregir la correcció de la fórmula
            if(redox_potential.ne.0.0) then
                    ! Potential: Converting the SHE red_potential into Gibbs energy
                redox_potential=(-redox_potential*Faraday)/4184.0d0-98.6991
                    ! Conversion from Volts to kcal/mol
                    ! 98.6991 --> Conversion of -4.28eV to kcal/mol
                redox_potential=redox_potential*redox_coef/hartree2kcal ! SHE red_potential into a.u adapted for the redox transoformation 
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
!##############################################################################!
!------------------------------------------------------------------------------!

    !======================================================================!
                        !   Reading the computation line   !
    !======================================================================!

read(2,'(A100)') line !-- Line corresponding to the 'Computation' line
if(index(line,"Computation:").ne.0.or.index(line,"computation:").ne.0) then
    index_initial=index(line,"oint")
    index_modulus=index(line,"dulus")
    index_trust=index(line,"dius")
    index_grid=index(line,"rid") 

    if(index(line,"Initial point=").ne.0.or.index(line,"initial point=").ne.0) then
            !-- Check the initial point of the scan
        index_initial_end=index(line,")")
        read(line(index_initial+6:index_initial_end-1),*) x0,y0,z0
        x0=x0*1.0d-4;y0=y0*1.0d-4;z0=z0*1.0d-4
        initial_position(1)=x0;initial_position(2)=y0;initial_position(3)=z0
    
            !-- Check there are no "missunderstandings" in the scans, i.e. computing a 3D with a constant field
        if(n_dim.eq.1) then
            if(axis(gpc).eq.1.and.initial_position(gpc).ne.0.0d0) then
                write(*,'(" Cannot perform the scan of ",A2," with a non-zero position! Stop!")') axis_name(gpc)
                write(*,*)
                stop
            end if
        end if
        if (n_dim.eq.2) then
            step=nint(1+abs(cos(0.5d0*PI*gpc)))
            do i=1,3,step 
                if (axis(gpc).eq.1.and.initial_position(i).ne.0.0d0) then
                        !-- Axis(gpc) --> Constant axis
                    write(*,'(" Cannot perform the 2D scan with a non-zero",A2," starting point")') axis_name(i)
                    write(*,*) "Stop!"
                    write(*,*)
                    stop 
                end if
            end do
        end if
        if (n_dim.eq.3.and.all(initial_position.ne.0.0d0)) then
            write(*,*) "The 3D scan cannot be performed with a non-zero initial position! Stop!"
            write(*,*)
            stop
        end if
    else
        write(*,*) "Bad input of initial point! STOP!"
        write(*,*)
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
    if(index(line,"Grid=").ne.0.or.index(line,"grid=").ne.0) then
        read(line(index_grid+4:index_grid+8),*) grid
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
!#############################################################################!
!-----------------------------------------------------------------------------!

    !=====================================================================!
                    !       Reading the extra line    !
    !=====================================================================!

read(2,'(A100)') line !-- Line corresponding to the 'Extra' line
if(index(line,"Extra:").ne.0.or.index(line,"extra:").ne.0) then ! Tolerance line exits
    index_tol=index(line,"lerance")
    read(line(index_tol+8:index_tol+14),*) tol
    if(tol.lt.0) then ! Tolerance cannot be lower than zero
        write(*,*) "Tolerance cannot be negative! Stop!"
        stop
    end if
else
    tol=1E-8
end if

!==============================================================================!
            !-----------------Echoing the input-----------------!
!==============================================================================!

write(*,*) "-------------------------------------------------------------------"
write(*,*) 
write(*,*) "            The input for this run is:"
write(*,'(" Central point ((x,y,z)·10^-4 a.u)",xF8.3,xF8.3,xF8.3)') x0,y0,z0
write(*,'(" Maximum radius of ",F7.5," a.u for the iterative resolutions")') radius
write(*,'(" Number of points computed: 0 (1D) ",I9," (2D) ",I10," (3D) ")') ,nint(PI*(radius**2-initial_position(gpc)**2)*1.0d8),nint(4.0d0*PI*radius**2*1.0d8)
if(grid.ge.1E+4) write(*,*) "       Warning! Grid too dense for scanning"
if(grid.le.1E+2) write(*,*) "       Warning! Grid too  thin for scanning"
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

read(2,'(A80)') line ! DO NOT TOUCH THIS LINE! SOMEHOW IT MAKES IT WORK
!===========================END OF ROUTE SECTION===============================!

!==============================================================================!
!##############################################################################!
!==============================================================================!

!==============================================================================!
                 !   Reading the reactant NLOPs   !
!==============================================================================!

read(2,'(A80)') line ! Line after the blank line right after the keywords section --> Begin reading reactant NLOPs
                     ! Stands for:  -#-#-#- Chemical species: Reactants -#-#-#-
if(index(line,"#- Chemical species: Reactants -#").ne.0.or.index(line,"#- chemical species: reactants -#").ne.0 &
& .or.index(line,"#- Chemical species: reactants -#").ne.0.or.index(line,"#- chemical species: Reactants -#").ne.0) then

    n=1     !-- At least there is going to be one reactant species
    E_iter=0.0d0; mu_iter=0.0d0; alpha_iter=0.0d0; beta_iter=0.0d0

    read(2,'(A80)') name_reactant_reactant(n)
    read(2,*) line
    call readvalues(E_r,mu_r,alpha_r,beta_r)
    
    do!while (index(line,"#- Chemical species: Products -#").eq.0.or.index(line,"#- chemical species: products -#").eq.0 &
    !& .or.index(line,"#- Chemical species: products -#").eq.0.or.index(line,"#- chemical species: Products -#").eq.0)
        E_iter=E_iter+E_r
        do i=1,3
            mu_iter(i)=mu_iter(i)+mu_r(i)
            do j=1,3
                alpha_iter(i,j)=alpha_iter(i,j)+alpha_r(i,j)
                do k=1,3
                    beta_iter(i,j,k)=beta_iter(i,j,k)+beta_r(i,j,k)
                end do
            end do
        end do  
        read(2,'(A80)') line
        if (index(line,"#- Chemical species: Products -#").eq.0.or.index(line,"#- chemical species: products -#").eq.0 &
            & .or.index(line,"#- Chemical species: products -#").eq.0.or.index(line,"#- chemical species: Products -#").eq.0) then
            exit
        else
            n=n+1
            name_reactant_reactant=trim(line)
            read(2,*) line
            call readvalues(E_r,mu_r,alpha_r,beta_r)
        end if
    end do
else
    write(*,*) "No input reactant chemical species! Stop!"
    call sleep(2)
    write(*,*) "If you consider this might not be true, remember adding:"
    write(*,*) " ' -#-#-#- Chemical species: Reactants -#-#-#- ' "
    write(*,*) "Right after the blank line after the keywords section"
    write(*,*)
    stop
end if

    !--- Redefinition of the original names for the later data handling
E_r=E_iter; mu_r=mu_iter; alpha_r=alpha_iter; beta_r=beta_iter

    !--- Print the (joint) properties of the reactants
if(Guillem(1).eqv..TRUE.) then
    write(*,*) "=============================================================================="
    write(*,*) "                     Joint properties for the reactants (R)"
    write(*,'(" Reactants: ",xA)') (name_reactant_reactant(i),i=1,n)
    write(*,*) "  --------------------------------------------------------------------------  "
    write(*,*) "                            Dipole moment"
    write(*,*) (mu_r(i),i=1,3)
    write(*,*)
    write(*,*) "                        Polarizability matrix"
    do i=1,3
            write(*,*) (alpha_r(i,j),j=1,3)
    end do
    write(*,*)
    write(*,*) "                     First hyperpolarizability tensor"
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
    write(*,*) "=============================================================================="
    write(*,*)
end if

!==============================================================================!
!##############################################################################!
!==============================================================================!

!==============================================================================!
                 !   Reading the product NLOPs   !
!==============================================================================!

if (index(line,"#- Chemical species: Products -#").ne.0.or.index(line,"#- chemical species: products -#").ne.0 &
& .or.index(line,"#- Chemical species: products -#").ne.0.or.index(line,"#- chemical species: Products -#").ne.0) then
    
    n=1 !-- At least there is going to be one reactant speciest
    E_iter=0.0d0; mu_iter=0.0d0; alpha_iter=0.0d0; beta_iter=0.0d0
    
    read(2,'(A80)') name_product_product(n)
    read(2,*) line
    call readvalues(E_p,mu_p,alpha_p,beta_p)

    do !while (index(line,"#- Small molecules: Reactants -#").eq.0.or.index(line,"#- small molecules: reactants -#").eq.0 &
    !& .or.index(line,"#- Small molecules: reactants -#").eq.0.or.index(line,"#- small molecules: Reactants -#").eq.0) 
        E_iter=E_iter+E_p
        do i=1,3
            mu_iter(i)=mu_iter(i)+mu_p(i)
            do j=1,3
                alpha_iter(i,j)=alpha_iter(i,j)+alpha_p(i,j)
                do k=1,3
                    beta_iter(i,j,k)=beta_iter(i,j,k)+beta_p(i,j,k)
                end do
            end do
        end do
        read(2,'(A80)') line
        if (index(line,"#- END OF FILE -#").ne.0.or.index(line,"#- End of file -#").ne.0.or.index(line,"#- end of file -#").ne.0) then
            write(*,*) "Warning! No small molecules at the input!"
            E_lig_p=0.0d0; mu_lig_p=0.0d0; alpha_lig_p=0.0d0; beta_lig_p=0.0d0
            E_lig_r=0.0d0; mu_lig_r=0.0d0; alpha_lig_r=0.0d0; beta_lig_r=0.0d0
            goto 99999
        else if (index(line,"#- Small molecules: Reactants -#").ne.0.or.index(line,"#- small molecules: reactants -#").ne.0 &
        & .or.index(line,"#- Small molecules: reactants -#").ne.0.or.index(line,"#- small molecules: Reactants -#").ne.0) then
            exit
        else
            n=n+1
            name_product_product(n)=trim(line)
            read(2,*)
            call readvalues(E_p,mu_p,alpha_p,beta_p) 
        end if
    end do
else
    write(*,*) "No input product chemical species! Stop!"
    call sleep(2)
    write(*,*) "If you consider this might not be true, remember adding:"
    write(*,*) " ' -#-#-#- Chemical species: Products -#-#-#- ' "
    write(*,*) "Right after the blank line of the last line of the nuclear relaxation polarizability matrix reactant"
    write(*,*)
    stop
end if

    !--- Redefinition of the original names for the later data handling
E_p=E_iter; mu_p=mu_iter; alpha_p=alpha_iter; beta_p=beta_iter

    !--- Print the (joint) properties of the products
if(Guillem(2).eqv..TRUE.) then
    write(*,*) "=============================================================================="
    write(*,*)
    write(*,*) "                     Joint properties for the products (P)"
    write(*,*) (name_product_product(i),i=1,n)
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) "                            Dipole moment"
    write(*,*) (mu_p(i),i=1,3)
    write(*,*)
    write(*,*) "                        Polarizability matrix"
    do i=1,3
            write(*,*) (alpha_p(i,j),j=1,3)
    end do
    write(*,*)
    write(*,*) "                    First hyperpolarizability tensor"
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
    write(*,*) "=============================================================================="
    write(*,*)
end if

!==============================================================================!
!##############################################################################!
!==============================================================================!

!==============================================================================!
                !   Reading the small molecules NLOPs   !
!==============================================================================!

        !--- Reading the NLOPs for the reactants small molecules
n=0;m=0
if (index(line,"#- Small molecules: Reactants -#").ne.0.or.index(line,"#- small molecules: reactants -#").ne.0 &
& .or.index(line,"#- small molecules: Reactants -#").ne.0.or.index(line,"#- Small molecules: reactants -#").ne.0) then

    n=1 !-- At least there is going to be one small molecule if this option is considered
    E_iter=0.0d0; mu_iter=0.0d0; alpha_iter=0.0d0; beta_iter=0.0d0

    read(2,'(A10)') name_reactant_ligand(n)
    read(2,*) line
    call readvalues_ligands(tol,E_lig_r,mu_lig_r,alpha_lig_r,beta_lig_r)
    
    do  !while (index(line,"#- Small molecules: Products -#").eq.0.or.index(line,"#- small molecules: products").eq.0 &
        !& .or.index(line,"#- small molecules: Products -#").eq.0.or.index(line,"#- Small molecules: products -#").eq.0)
        E_iter=E_iter+E_lig_r
        do i=1,3
            mu_iter(i)=mu_iter(i)+mu_lig_r(i)
            do j=1,3
                alpha_iter(i,j)=alpha_iter(i,j)+alpha_lig_r(i,j)
                do k=1,3
                    beta_iter(i,j,k)=beta_iter(i,j,k)+beta_lig_r(i,j,k)
                end do
            end do
        end do
        read(2,'(A80)') line
        if (index(line,"#- END OF FILE -#").ne.0.or.index(line,"#- End of file -#").ne.0.or.index(line,"#- end of file -#").ne.0) then
            write(*,*) "Warning! No small molecules detected at the products."
            E_lig_r=E_iter; mu_lig_r=mu_iter; alpha_lig_r=alpha_iter; beta_lig_r=beta_iter
            E_lig_p=0.0d0; mu_lig_p=0.0d0; alpha_lig_p=0.d0; beta_lig_p=0.0d0
            goto 99999
        else if (index(line,"#- Small molecules: Products -#").ne.0.or.index(line,"#- small molecules: products -#").ne.0 &
        & .or.index(line,"#- Small molecules: products -#").ne.0.or.index(line,"#- small molecules: Products -#").ne.0) then
            !-- The reaction contains small molecules at the reactants 
            exit
        else
            n=n+1
            name_reactant_ligand(n)=trim(line)
            read(2,*) line
            call readvalues_ligands(tol,E_lig_r,mu_lig_r,alpha_lig_r,beta_lig_r)
        end if
    end do
else
    ! No reactant ligand was found
    read(2,*) line
end if
     
    !--- Redefinition of the original names for the later data handling
E_lig_r=E_iter; mu_lig_r=mu_iter; alpha_lig_r=alpha_iter; beta_lig_r=beta_iter

        !--- Reading the NLOPs for the products small molecules
if (index(line,"#- Small molecules: Products -#").ne.0.or.index(line,"#- small molecules: products -#").ne.0 &
& .or.index(line,"#- small molecules: Products -#").ne.0.or.index(line,"#- Small molecules: products -#").ne.0) then


    m=1 !-- At least there is going to be one small molecule if this option is considered
    E_iter=0.0d0; mu_iter=0.0d0; alpha_iter=0.0d0; beta_iter=0.0d0
    
    read(2,'(A10)') name_product_ligand(m)
    read(2,*) line
    call readvalues_ligands(tol,E_lig_p,mu_lig_p,alpha_lig_p,beta_lig_p)
     
    do!while (index(line,"#- END OF FILE -#").eq.0.or.index(line,"#- End of file -#").eq.0.or.index(line,"#- end of file -#").eq.0)
        E_iter=E_iter+E_lig_p
        do i=1,3
            mu_iter(i)=mu_iter(i)+mu_lig_p(i)
            do j=1,3
                alpha_iter(i,j)=alpha_iter(i,j)+alpha_lig_p(i,j)
                do k=1,3
                    beta_iter(i,j,k)=beta_iter(i,j,k)+beta_lig_p(i,j,k)
                end do
            end do
        end do
        read(2,'(A80)') line
        if (index(line,"#- END OF FILE -#").ne.0.or.index(line,"#- End of file -#").ne.0.or.index(line,"#- end of file -#").ne.0) then 
            !-- End of the file. Reading complete
            exit
        else
            m=m+1
            name_product_ligand(m)=trim(line)
            read(2,*) line
            call readvalues_ligands(tol,E_lig_p,mu_lig_p,alpha_lig_p,beta_lig_p) 
        end if
    end do
else
    !--- End of file
    read(2,*) line
end if

    !--- Redefinition of the original names for the later data handling
E_lig_p=E_iter; mu_lig_p=mu_iter; alpha_lig_p=alpha_iter; beta_lig_p=beta_iter
99999 continue
do i=1,3
    mu_lig(i)=mu_lig_p(i)-mu_lig_r(i)
    do j=1,3
        alpha_lig(i,j)=alpha_lig_p(i,j)-alpha_lig_r(i,j)
        do k=1,3
            beta_lig(i,j,k)=-(beta_lig_p(i,j,k)-beta_lig_r(i,j,k))
        end do
    end do
end do

    !--- Print the joint properties of the small molecules
if(Guillem(3).eqv..TRUE.) then
    write(*,*) "=============================================================================="
    write(*,*)
    write(*,*) "                Joint properties for the small molecules (Lig)"
    if (n.ge.1) write(*,'(" Reactants: ",Ax)') (trim(name_reactant_ligand(i)),i=1,n)
    if (m.ge.1) write(*,'(" Products:  ",Ax)') (trim(name_product_ligand(i)),i=1,m)
    write(*,*) "------------------------------------------------------------------------------"
    write(*,*) "                            Dipole moment"
    write(*,*) (mu_lig(i),i=1,3)
    write(*,*)
    write(*,*) "                        Polarizability matrix"
    do i=1,3
            write(*,*) (alpha_lig(i,j),j=1,3)
    end do
    write(*,*)
    write(*,*) "                    First hyperpolarizability tensor"
    write(*,*) "X _ _"
    do j=1,3
            write(*,*) (beta_lig(1,j,k),k=1,3)
    end do
    write(*,*) "Y _ _"
    do j=1,3
            write(*,*) (beta_lig(2,j,k),k=1,3)
    end do
    write(*,*) "Z _ _"
    do j=1,3
            write(*,*) (beta_lig(3,j,k),k=1,3)
    end do
    write(*,*) "=============================================================================="
    write(*,*)
end if

!==============================================================================!
!##############################################################################!
!==============================================================================!
    
!==============================================================================!
          !   Computing the molecular properties of the system   !
!==============================================================================!

E_0=E_p-E_r+E_lig_p-E_lig_r ! Energy in atomic units
do i=1,3
    mu(i)=mu_p(i)-mu_r(i)!+mu_lig_p(i)-mu_lig_r(i)
    do j=1,3
        alpha(i,j)=alpha_p(i,j)-alpha_r(i,j)!+alpha_lig_p(i,j)-alpha_lig_r(i,j)
        do k=1,3
            beta(i,j,k)=-(beta_p(i,j,k)-beta_r(i,j,k))!+beta_lig_p(i,j,k)-beta_lig_r(i,j,k))
        end do
    end do
end do

    !-- Echoing the properties of the run
write(*,*)
write(*,*) "=============================================================================="
write(*,'("     Title of the job: ",xA80)') title
write(*,*) "------------------------------------------------------------------------------"
write(*,*)
write(*,'(" Target barrier for the study:",xF6.2," kcal/mol")') target_barrier
write(*,*) 
write(*,*) "=============================================================================="
write(*,*) "     Thermodynamics and non-linear optical properties of the system"
write(*,*) "------------------------------------------------------------------------------"
write(*,'("    Gibbs free energy (ΔG)",xF15.4,x"kcal/mol")') E_0*hartree2kcal
write(*,*) "            Dipole moment vector (μ)"
write(*,'(xxxxxE11.4,xxxxxE11.4,xxxxxE11.4)') (mu(i)+mu_lig_p(i)-mu_lig_r(i),i=1,3)
write(*,*)
write(*,*) "          Polarizability matrix (α)"
do i=1,3
        write(*,'(xxxxxE11.4,xxxxxE11.4,xxxxxE11.4)') (alpha(i,j)+alpha_lig_p(i,j)-alpha_lig_r(i,j),j=1,3)
end do
write(*,*)
write(*,*) "      First hyperpolarizability tensor (β)"
write(*,*) "X _ _"
do j=1,3
        write(*,'(xxxxxE11.4,xxxxxE11.4,xxxxxE11.4)') (beta(1,j,k)+(-beta_lig_p(1,j,k)+beta_lig_r(1,j,k)),k=1,3)
end do
write(*,*) "Y _ _"
do j=1,3
        write(*,'(xxxxxE11.4,xxxxxE11.4,xxxxxE11.4)') (beta(2,j,k)+(-beta_lig_p(2,j,k)+beta_lig_r(2,j,k)),k=1,3)
end do
write(*,*) "Z _ _"
do j=1,3
        write(*,'(xxxxxE11.4,xxxxxE11.4,xxxxxE11.4)') (beta(3,j,k)+(-beta_lig_p(3,j,k)+beta_lig_r(3,j,k)),k=1,3)
end do
write(*,*)
write(*,*) "------------------------------------------------------------------------------"
write(*,*) "=============================================================================="
write(*,*)
write(*,*) "            Please check everything is correct before running."
write(*,'(" If you wish to make any change at the input file",A20,",")') trim(file)
write(*,*) "                        please press Ctrl+C"
write(*,*) 
!--- S'ha de fer un if segons el n_dim
!write(*,'("     Searching the ",A2,x"minimum field strength for the input-imposed conditions")') axis_name(gpc)
!write(*,*) "                           Please wait..."
!--- S'ha de fer un if segons el n_dim
write(*,*) 
write(*,*)
write(*,*) "===================================================================================================================================================================================="
!=============================================================================!
!----------------------------VOAM STARTS HERE---------------------------------!
!=============================================================================!
if (n_dim.eq.0) then !-- Scanning VOAM is performed
    call VOAM_0D(axis_name,grid,radius)
end if
if (n_dim.eq.1) then !-- 1D resolution is performed
    call VOAM_1D(Guillem,order_nlop,gpc,radius,axis,initial_position,axis_name,tol)
end if
if (n_dim.eq.2) then !-- 2D resolution is performed
    call VOAM_2D(Guillem,order_nlop,gpc,radius,axis,initial_position,axis_name,grid)  
end if
if (n_dim.eq.3) then !-- 3D resolution is performed
    call VOAM_3D(Guillem,order_nlop,gpc,radius,axis,initial_position,axis_name,grid)
end if
!=============================================================================!
!----------------------------VOAM ENDS HERE ----------------------------------!
!=============================================================================!
write(*,*) "===================================================================================================================================================================================="
write(*,*)
write(*,*)
write(*,*) "  Program written by Dr. Pau Besalú and Guillem Pey* for his Master Thesis"
write(*,*) "supervised by Dr. Josep Maria Luis Luis and Prof. Julio Lloret-Fillol in the"
write(*,*) "    Master in Advanced Catalysis and Molecular Modeling (MACMoM)"
write(*,*) "                at Universitat de Girona (UdG) on 2024" 
write(*,*) 
write(*,*) 
write(*,*) 
End
!------------------------------------------------------------------------------!
Subroutine readvalues(E_i,mu_i,alpha_i,beta_i)
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
!##############################################################################!
Subroutine readvalues_ligands(tol,E_i,mu_i,alpha_i,beta_i)
implicit none
character*100 line
double precision :: E_i,tol
double precision, dimension (3) :: mu_i
double precision, dimension (3,3) :: alpha_i,alpha_tmp
double precision, dimension (3,3,3) :: beta_i
integer i,j

read(2,*) E_i
read(2,*) line
read(2,*) mu_i(1),mu_i(2),mu_i(3)

if (mu_i(3).ge.0.0d0) then !and.abs(mu_i(3)).gt.tol) then
    continue
    if (abs(mu_i(1)).lt.tol*1E+6.and.abs(mu_i(2)).lt.tol*1E+6) then
        continue
    else
        write(*,*) "Warning! The X and Y components of the dipole moment are not neglegible!"
    end if
else !if (mu_i(3).lt.0.0d0.and.abs(mu_i(3)).gt.tol) then
    write(*,*)
    write(*,*) "Warning! One of the input small molecules is not properly oriented!"
    write(*,*) "The Z component of the dipole moment should be positive."
end if

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
