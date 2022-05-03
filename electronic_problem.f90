!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: electronic_problem                                                   !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: On-the-fly electronic-structure calculations                    !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief On-the-fly electronic-structure calculations.
MODULE electronic_problem
  USE variables
  USE kinds
  USE analytical_potentials

  IMPLICIT NONE

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Calculation of electronic properties on-the-fly              !
  !---------------------------------------------------------------------------!
  !> Electronic energies (adiabatic or spin-(a)diabatic), forces and
  !! non-adiabatic couplings are compueted at the trajectory position.
  !> @param[in] trajlabel label of the trajectory
  !> @param[in] Q position of the trajectory
  !> @param istate integer index running over the electronic states
  !> @param V array of (a)diabatic Hamiltonian
  !> @param G array of gradients of the (a)diabatic Hamiltonian
  !> @param NAC array of non-adiabatic couplings
  !> @return Energies, forces and non-adiabatic couplings are stored in the
  !! arrays BOenergy, BOforce, coup.
  SUBROUTINE BOproblem(Q,trajlabel)

    INTEGER      ,INTENT(IN) :: trajlabel
    REAL(KIND=DP),INTENT(IN) :: Q(n_dof)
    INTEGER                  :: i_state,i,j
    REAL(KIND=DP)            :: V(nstates,nstates),G(nstates,nstates,n_dof),&
                                NAC(nstates,nstates,n_dof)   

    IF(new_potential)       &
       CALL new_model_potentials(V,G,NAC,Q)
    IF(.NOT. new_potential) &
       CALL sub_model1_VG_NAC(V,G,NAC,Q,n_dof,nstates,model_potential,option)

    DO i_state=1,nstates
       BOenergy(trajlabel,i_state)  =   V(i_state,i_state)
       BOforce(trajlabel,i_state,:) = - G(i_state,i_state,:)
    END DO
    coup(trajlabel,:,:,:) =  NAC(:,:,:)

    IF(spin_dia) THEN
      coup_so(trajlabel,:,:) = CMPLX(0.0_dp,0.0_dp)
      DO i=1,nstates
        DO j=1,nstates
          IF(i/=j) coup_so(trajlabel,i,j)=V(i,j)
        ENDDO
      ENDDO
    ENDIF


  END SUBROUTINE BOproblem

END MODULE electronic_problem





