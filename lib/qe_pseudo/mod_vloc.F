!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------
MODULE vlocal
  !
  !! The variables needed for the local potential in reciprocal space.
  !
  USE kinds,       ONLY : DP
  !USE parameters,  ONLY : ntypx
  !
  SAVE
  !
  ! [FP] strf already obtained by calling struc_fact()
  !      e.g. by qe_pseudo_init in p2y or PP_uspp_init.
  !      We comment this one.
  !COMPLEX(DP), ALLOCATABLE :: strf(:,:)
  !! the structure factor
  REAL(DP), ALLOCATABLE :: vloc(:,:)
  !! the local potential for each atom type
  !REAL(DP) :: starting_charge(ntypx)
  !! the atomic charge used to start with
  !
END MODULE vlocal

! This module contains the variables and routines necessary to the
! implementation
! of the two-dimensional Coulomb cutoff. Details of the implementation can be
! found in:
!
! Sohier, T., Calandra, M., & Mauri, F. (2017),
! "Density functional perturbation theory for gated two-dimensional
! heterostructures:
! Theoretical developments and application to flexural phonons in graphene."
! Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
!
!----------------------------------------------------------------------------
MODULE Coul_cut_2D
  !----------------------------------------------------------------------------
  !! This module contains the variables and subroutines needed for the
  !! 2D Coulomb cutoff.
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : tpi, pi
  !
  SAVE
  !
  LOGICAL :: do_cutoff_2D = .FALSE.
  !! flag for 2D cutoff. If true, the cutoff is active.
  REAL(DP) :: lz
  !! The distance in the out-plne direction after which potential are cut off.
  !
  REAL(DP), ALLOCATABLE :: cutoff_2D(:)
  !! The factor appended to the Coulomb Kernel to cut off potentials
  REAL(DP), ALLOCATABLE :: lr_Vloc(:,:)
  !! The long-range part of the local part of the ionic potential
  !
  ! [FP] This module contains several subroutines 
  !       in the original QE implementation
END MODULE Coul_cut_2D

MODULE esm
  !--------------------------------------------------------------------------
  !! This module contains the variables and subroutines needed for the
  !! EFFECTIVE SCREENING MEDIUM (ESM) METHOD.
  !
  !! Original version by Minoru Otani (AIST), Yoshio Miura (Tohoku U.),
  !! Nicephore Bonet (MIT), Nicola Marzari (MIT), Brandon Wood (LLNL),
  !! Tadashi Ogitsu (LLNL).
  !! Constant bias potential (constant-mu) method by Minoru Otani (AIST) and
  !! Nicephore Bonnet (AIST).
  !
  !! Contains subroutines for implementation of:
  !! 1) ESM (Effective Screening Medium Method) developed by M. Otani and
  !!    O. Sugino (see PRB 73, 115407 [2006]);
  !! 2) Constant-mu method developed by N. Bonnet, T. Morishita, O. Sugino,
  !!    and M. Otani (see PRL 109, 266101 [2012]).
  !
  !! ESM enables description of a surface slab sandwiched between two
  !! semi-infinite media, making it possible to deal with polarized surfaces
  !! without using dipole corrections. It is useful for simulating interfaces
  !! with vacuum, one or more electrodes, or an electrolyte.
  !
  !! Constant-mu scheme with the boundary condition 'bc2' and 'bc3' enables
  !! description of the system is connected to a potentiostat which preserves
  !! the Fermi energy of the system as the target Fermi energy (mu).
  !
  !! Modified subroutines for calculating the Hartree potential, the local
  !! potential, and the Ewald sum are contained here, along with subroutines for
  !! calculating force contributions based on the modified local potential and
  !! Ewald term. Constant-mu parts are contained in the fcp.f90.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: do_comp_esm, esm_nfit, esm_efield, esm_w, esm_a, esm_bc,   &
            mill_2d, imill_2d, ngm_2d,                                 &
            esm_init, esm_hartree, esm_local, esm_ewald, esm_force_lc, &
            esm_force_ew, esm_printpot, esm_summary
  PUBLIC :: esm_stres_har, esm_stres_ewa, esm_stres_loclong
  !
  !
  LOGICAL :: do_comp_esm=.FALSE.
  !! if TRUE the ESM is active
  INTEGER :: esm_nfit
  !! number of grid points for fit at edges
  REAL(DP) :: esm_efield
  !! field strength
  REAL(DP) :: esm_w
  !! ESM offset from cell edge
  REAL(DP) :: esm_a
  !! smoothness parameter
  CHARACTER(LEN=3) :: esm_bc
  !! pbc, bc1, ..., bc4
  INTEGER, ALLOCATABLE :: mill_2d(:,:)
  !! miller index on vectors
  INTEGER, ALLOCATABLE :: imill_2d(:,:)
  INTEGER :: ngm_2d = 0
  !! ngm_2d = total number of vectors (h,k) on this proc, excluding
  !! duplicates with different l values
  !
  REAL(DP), EXTERNAL   :: qe_erf, qe_erfc
  !
  ! [FP] This module contains several subroutines
  !      in the original QE implementation 
END MODULE esm
