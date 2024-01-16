
SUBROUTINE newton_0D( &
    & component_guess, component_total, &
    & tableau, logK,                    &
    & N_components, N_species,          &
    & atol, maxiter,                    &
    & jacobian_template,                &
    & iter, info,                       &
    & difference, species_conc          &
    & )

    IMPLICIT NONE

    ! SUBROUTINE ARGUMENTS
    DOUBLE PRECISION, DIMENSION (N_components)              :: component_guess   !current guess for component concentrations
    DOUBLE PRECISION, DIMENSION (N_components)              :: component_total   !known total amount of components
    DOUBLE PRECISION, DIMENSION (N_species,N_components)    :: tableau           !the tableau ...
    DOUBLE PRECISION, DIMENSION (N_species)                 :: logK              !equilibrium constants
    INTEGER                                                 :: N_components      !number of components in tableau
    INTEGER                                                 :: N_species         !number of species in tableau
    DOUBLE PRECISION                                        :: atol              !absolute tolerance for deviating from component_total
    INTEGER                                                 :: maxiter           !maximum iterations to conduct
    DOUBLE PRECISION, DIMENSION (N_species,N_components**2) :: jacobian_template !part of jacobian that stays constant for a given tableau
    INTEGER                                                 :: iter              !counts iterations conducted (result)
    INTEGER                                                 :: info              !state (success/error) of dgesv-solver (result)
    DOUBLE PRECISION, DIMENSION (N_components,1)            :: difference        !remaining difference between guessed and known component_totals (result)
    DOUBLE PRECISION, DIMENSION (N_species)                 :: species_conc      !species concentrations (result)
    
    ! INTRA SUBROUTINE VARIABLES
    DOUBLE PRECISION, DIMENSION (N_components,N_components) :: jacobian          !the jacobian matrix for the current component-concentration guess
    INTEGER, DIMENSION (N_components)                       :: pivot             !pivot indices that define the permutation matrix for dgesv-solver
    DOUBLE PRECISION, DIMENSION (N_components)              :: relax             !relaxation factor for calculating new guess
    
    ! convert component_guess to logarithmic form
    WHERE (component_guess < 1e-44)
        component_guess = -44
    ELSEWHERE
        component_guess = log10(component_guess)
    END WHERE

    DO iter = 1, maxiter

        species_conc = 10**(matmul(tableau, component_guess) + logK)
             
        difference = reshape((-(matmul(transpose(tableau), species_conc) - component_total)), [N_components,1])

        IF (all(abs(difference) < atol)) THEN
            exit
        END IF

        jacobian = reshape( &
            & log(10.0) * matmul(transpose(jacobian_template), species_conc), &
            & [N_components,N_components], order = [2,1]                      &
        & )

        call dgesv(N_components, 1, jacobian, N_components, pivot, difference, N_components, info)

        IF (info .ne. 0) THEN
            exit
        END IF

        relax = -1 * [difference] / (0.5 * component_guess)
        WHERE (relax < 1) relax = 1
        relax = relax**(-1)
        component_guess = component_guess + relax * reshape(difference, [N_components])

        !WHERE (component_guess < -44) component_guess = -44
        !WHERE (component_guess < 1) component_guess = -1
            
    END DO

    ! recalculate difference
    difference = reshape((-(matmul(transpose(tableau), species_conc) - component_total)), [N_components,1])
    ! convert component_guess back to non-logarithmic form
    component_guess = 10**component_guess

END SUBROUTINE newton_0D


SUBROUTINE pcfm_0D( &
    & component_guess, component_total, &
    & tableau, logK,                    &
    & N_components, N_species,          &
    & epsilon, maxiter,                 &
    & iter,                             &
    & difference, species_conc          &
    & )

    IMPLICIT NONE

    ! SUBROUTINE ARGUMENTS
    DOUBLE PRECISION, DIMENSION (N_components)              :: component_guess   !current guess for component concentrations
    DOUBLE PRECISION, DIMENSION (N_components)              :: component_total   !known total amount of components
    DOUBLE PRECISION, DIMENSION (N_species,N_components)    :: tableau           !the tableau ...
    DOUBLE PRECISION, DIMENSION (N_species)                 :: logK              !equilibrium constants
    INTEGER                                                 :: N_components      !number of components in tableau
    INTEGER                                                 :: N_species         !number of species in tableau
    DOUBLE PRECISION                                        :: epsilon           !stop criterion
    INTEGER                                                 :: maxiter           !maximum iterations to conduct
    INTEGER                                                 :: iter              !counts iterations conducted (result)
    DOUBLE PRECISION, DIMENSION (N_components,1)            :: difference        !remaining difference between guessed and known component_totals (result)
    DOUBLE PRECISION, DIMENSION (N_species)                 :: species_conc      !species concentrations (result)
    
    ! INTRA SUBROUTINE VARIABLES
    DOUBLE PRECISION, DIMENSION (N_components)              :: p0j               !smallest positive coefficient in tableau for a component 
    DOUBLE PRECISION, DIMENSION (N_components,2)            :: sum_quantities    !sum_j_reac (column1) and sum_j_prod (column2)
    DOUBLE PRECISION, DIMENSION (N_components)              :: attenuation_factors !attenuation factors
    DOUBLE PRECISION, DIMENSION (N_components)              :: component_guess_log !log10 of component_guess
    INTEGER, ALLOCATABLE                                    :: index_reac(:), index_prod(:)
    INTEGER, ALLOCATABLE                                    :: index_tot0(:)
    DOUBLE PRECISION                                        :: attenuation_factor_new
    INTEGER                                                 :: i, j, ix          !iterators


    ! fill p0j
    DO j = 1, N_components
        p0j(j) = minval(tableau(:,j), MASK = tableau(:,j) > 0)
    END DO

    ! where TOT_X is 0
    index_tot0 = pack([(ix, ix = 1, N_components)], component_total == 0)
    
    DO iter = 1, maxiter

        ! initialize attenuation_factors
        attenuation_factors = 0.3!0.01

        !calculate species concentrations
        component_guess_log = log10(component_guess)
        component_guess_log(index_tot0) = -100
        component_guess(index_tot0) = 0
        species_conc = 10**(matmul(tableau, component_guess_log) + logK)

        difference = reshape((-(matmul(transpose(tableau), species_conc) - component_total)), [N_components,1])
        IF (all(abs([difference]) < epsilon)) exit

        !for each component ...
        DO j = 1, N_components

            !... calculate sum_j_reac, sum_j_prod
            index_reac = pack([(ix,ix=1,size(tableau(:,j)))], tableau(:,j) > 0)
            index_prod = pack([(ix,ix=1,size(tableau(:,j)))], tableau(:,j) < 0)

            IF (component_total(j) >= 0) THEN
                sum_quantities(j,1) = sum(tableau(index_reac,j) * species_conc(index_reac))
                sum_quantities(j,2) = component_total(j) + sum(abs(tableau(index_prod,j)) * species_conc(index_prod))
            ELSE
                sum_quantities(j,1) = abs(component_total(j)) + sum(tableau(index_reac,j) * species_conc(index_reac))
                sum_quantities(j,2) = sum(abs(tableau(index_prod,j)) * species_conc(index_prod))
            END IF

            ! ... calculate attenuation factor
            !attenuation_factor_new = attenuation_factors(j)
            !IF (sum_quantities(j,1) > sum_quantities(j,2)) THEN
            !    attenuation_factor_new = 0.9 - 0.8 * (sum_quantities(j,2) / sum_quantities(j,1))
            !ELSE IF (sum_quantities(j,2) > sum_quantities(j,1)) THEN
            !    attenuation_factor_new = 0.9 - 0.8 * (sum_quantities(j,1) / sum_quantities(j,2))
            !END IF

            !IF (attenuation_factor_new > attenuation_factors(j)) THEN
            !    attenuation_factors(j) = attenuation_factor_new
            !END IF

            ! ... calculate updated concentration
            component_guess(j) = attenuation_factors(j) * component_guess(j) *           &
            &                    (sum_quantities(j,2) / sum_quantities(j,1))**(1/p0j(j)) &
            &                    + (1 - attenuation_factors(j)) * component_guess(j)
            
        END DO

    END DO

END SUBROUTINE pcfm_0D


SUBROUTINE solve_tableau(               &  
    & component_total,                  &
    & tableau, logK,                    &
    & N_components, N_species, N_grid,  &
    & iter_pcfm,                        &
    & iter_newton, info_newton,         &
    & difference, species_conc,         &
    & success                           &
    & )
    
    IMPLICIT NONE

    ! SUBROUTINE ARGUMENTS
    DOUBLE PRECISION, DIMENSION (N_components * N_grid)     :: component_total     !known total amounts of components
    DOUBLE PRECISION, DIMENSION (N_species,N_components)    :: tableau             !the tableau ...
    DOUBLE PRECISION, DIMENSION (N_species)                 :: logK                !equilibrium constants
    INTEGER                                                 :: N_components        !number of components in tableau
    INTEGER                                                 :: N_species           !number of species in tableau
    INTEGER                                                 :: N_grid              !number of grid layers
    INTEGER, DIMENSION (N_grid)                             :: iter_pcfm           !iterations conducted for pcfm
    INTEGER, DIMENSION (N_grid)                             :: iter_newton         !iterations conducted for each layer
    INTEGER, DIMENSION (N_grid)                             :: info_newton         !exit status of solving soutine for each layer
    DOUBLE PRECISION, DIMENSION (N_grid,N_components)       :: difference          !remaining difference between guessed and known total components
    DOUBLE PRECISION, DIMENSION (N_grid,N_species)          :: species_conc        !result: species concentrations
    INTEGER, DIMENSION (N_grid)                             :: success             !convergence reached? 0: no, 1: first round-pcfm, 2: first round-newton, 3: second round-pcfm, 4: second round-newton
        
    ! INTRA SUBROUTINE VARIABLES
    DOUBLE PRECISION, DIMENSION (N_components * N_grid)     :: component_guess     !guess for component concentrations
    DOUBLE PRECISION, DIMENSION (N_species,N_components**2) :: jacobian_template   !part of jacobian matrix that remains the same for a given tableau
    DOUBLE PRECISION, DIMENSION (N_components)              :: component_guess_ly  !guess for component concentrations (current layer)
    DOUBLE PRECISION, DIMENSION (N_components)              :: component_total_ly  !known total components (current layer)
    INTEGER                                                 :: iter_pcfm_ly        !number of iterations of last pcfm run
    INTEGER                                                 :: iter_newton_ly      !number of iterations conducted (for current layer)
    INTEGER                                                 :: info_newton_ly      !exit status of solving routine (for current layer)
    DOUBLE PRECISION, DIMENSION (N_components,1)            :: difference_ly       !remaining difference (current layer)
    DOUBLE PRECISION, DIMENSION (N_species)                 :: species_conc_ly     !species concentrations (current layer)
    INTEGER                                                 :: success_ly          !success status for current layer
    INTEGER, DIMENSION (N_components)                       :: index               !index for getting components of one layer from input vectors
    INTEGER                                                 :: i, j, column, layer !counters

    ! SOLVER PARAMETERS
    DOUBLE PRECISION, PARAMETER                             :: atol = 1e-12        !absolute tolerance for deviating from component_total
    INTEGER, PARAMETER                                      :: maxiter_newton = 500!maximum iterations for newton method
    INTEGER, PARAMETER                                      :: maxiter_pcfm = 1000 !maximum iterations for pcfm method
    INTEGER, PARAMETER                                      :: maxiter_pcfm_final = 2000 ! maximum pcfm iterations of problem has tp be solved with it ...
    DOUBLE PRECISION, PARAMETER                             :: epsilon1 = 0.5      !tolerated atol in pcfm round 1
    DOUBLE PRECISION, PARAMETER                             :: epsilon2 = 1e-3     !tolerated atol in pcfm round 2

    !construct jacobian_template
    column = 1
    DO i = 1, N_components
        DO j = 1, N_components
            jacobian_template(:,column) = tableau(:,i) * tableau(:,j)
            column = column + 1
        END DO
    END DO

    !initial guess for component concentrations
    component_guess = abs(component_total / 1000)
    
    DO layer = 1, N_grid

        !initialiation
        success_ly     = 0
        info_newton_ly = 999
        iter_newton_ly = 999
        iter_pcfm_ly   = 999
        difference_ly  = 999

        ! extract component_total and component_guess for this layer
        index = [(i, i = layer, N_grid * N_components, N_grid)]

        IF (layer == 1) THEN
            component_guess_ly = component_guess(index)
        ELSE IF (success(layer-1) == 0) THEN
            component_guess_ly = component_guess(index)
        ELSE
            component_guess_ly = component_guess(index - 1)
        END IF
       
        component_total_ly = component_total(index)

        ! pcfm-preconditioning 1
        call pcfm_0D ( &
            & component_guess_ly, component_total_ly, &
            & tableau, logK,                          &
            & N_components, N_species,                &
            & epsilon1, maxiter_pcfm,                 &
            & iter_pcfm_ly,                           &
            & difference_ly, species_conc_ly          &
            & )

        ! test if is guess is already good enough
        IF (all(abs(difference_ly) < atol)) THEN
            success_ly = 1
            GOTO 10
        END IF

        ! newton 1
        call newton_0D( &
            & component_guess_ly, component_total_ly, &
            & tableau, logK,                          &
            & N_components, N_species,                &
            & atol, maxiter_newton,                   &
            & jacobian_template,                      &
            & iter_newton_ly, info_newton_ly,         &
            & difference_ly, species_conc_ly          &
            & )

        ! test if is guess is already good enough
        IF (all(abs(difference_ly) < atol)) THEN
            success_ly = 2
            GOTO 10
        END IF

        ! new try if failed before with longer pcfm
        IF (layer == 1) THEN
            component_guess_ly = component_guess(index)
        ELSE IF (success(layer-1) == 0) THEN
            component_guess_ly = component_guess(index)
        ELSE
            component_guess_ly = component_guess(index - 1)
        END IF
        
        ! pcfm-preconditioning 2
        call pcfm_0D ( &
            & component_guess_ly, component_total_ly, &
            & tableau, logK,                          &
            & N_components, N_species,                &
            & epsilon2, maxiter_pcfm,                 &
            & iter_pcfm_ly,                           &
            & difference_ly, species_conc_ly          &
            & )

        ! test if is guess is already good enough
        IF (all(abs(difference_ly) < atol)) THEN
            success_ly = 3
            GOTO 10
        END IF

        ! newton 2
        call newton_0D( &
            & component_guess_ly, component_total_ly, &
            & tableau, logK,                          &
            & N_components, N_species,                &
            & atol, maxiter_newton,                   &
            & jacobian_template,                      &
            & iter_newton_ly, info_newton_ly,         &
            & difference_ly, species_conc_ly          &
            & )

        ! test if is guess is already good enough
        IF (all(abs(difference_ly) < atol)) THEN
            success_ly = 4
            GOTO 10
        END IF

        ! new try if failed before -> solve solely with pcfm
        IF (layer == 1) THEN
            component_guess_ly = component_guess(index)
        ELSE IF (success(layer-1) == 0) THEN
            component_guess_ly = component_guess(index)
        ELSE
            component_guess_ly = component_guess(index - 1)
        END IF
        
        call pcfm_0D ( &
            & component_guess_ly, component_total_ly, &
            & tableau, logK,                          &
            & N_components, N_species,                &
            & atol, maxiter_pcfm_final,               &
            & iter_pcfm_ly,                           &
            & difference_ly, species_conc_ly          &
            & )

        
        ! test if is finally good enough
        IF (all(abs(difference_ly) < atol)) success_ly = 5
        

10      component_guess(index) = component_guess_ly
        iter_newton(layer)     = iter_newton_ly
        info_newton(layer)     = info_newton_ly
        iter_pcfm(layer)       = iter_pcfm_ly
        difference(layer,:)    = [difference_ly]
        species_conc(layer,:)  = species_conc_ly
        success(layer)         = success_ly

    END DO

END SUBROUTINE solve_tableau

