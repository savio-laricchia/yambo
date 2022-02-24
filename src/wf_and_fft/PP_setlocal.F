!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! TB
! setup of the gate, search for 'TB'
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE setlocal
  !----------------------------------------------------------------------
  !! This routine computes the local potential in real space vltot(ir).
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : eps8
  USE ions_base,         ONLY : zv, ntyp => nsp
  USE cell_base,         ONLY : omega
  USE gvect,             ONLY : igtongl, gg
  !USE scf,               ONLY : v_of_0
  USE vlocal,            ONLY : vloc
  USE fft_base,          ONLY : dfftp
  USE gvect,             ONLY : ngm
  USE control_flags,     ONLY : gamma_only
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_sum
  USE esm,               ONLY : esm_local, esm_bc, do_comp_esm
  USE Coul_cut_2D,       ONLY : do_cutoff_2D, cutoff_local 
  !
  ! [FP] Remember to add call to struc_fact() to obtain strf!!!
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:), v_corr(:)
  ! auxiliary variable
  INTEGER :: nt, ng
  ! counter on atom types
  ! counter on g vectors
  !
  ALLOCATE( aux(dfftp%nnr) )
  aux(:) = (0.d0,0.d0)
  !
  DO nt = 1, ntyp
     DO ng = 1, ngm
        aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vloc(igtongl(ng),nt) &
                            * strf(ng,nt)
     ENDDO
  ENDDO
  !
  IF (gamma_only) THEN
     DO ng = 1, ngm
        aux(dfftp%nlm(ng)) = CONJG( aux(dfftp%nl(ng)) )
     ENDDO
  ENDIF
  !
  ! [FP] Effective Screening Medium: DISABLED
  IF ( do_comp_esm .AND. ( esm_bc .NE. 'pbc' ) ) THEN
     !
     ! ... Perform ESM correction to local potential
     !
     CALL esm_local( aux )
     !
  ENDIF
  !
  ! [FP] QE 2D cutoff by TB: DISABLED
  IF ( do_cutoff_2D ) THEN
     !
     ! ... re-add the CUTOFF fourier transform of erf function
     !
     CALL cutoff_local( aux )
     !
  ENDIF 
  !
  ! ... v_of_0 is (Vloc)(G=0)
  !
  !v_of_0 = 0.0_DP
  !IF (gg(1) < eps8) v_of_0 = DBLE( aux(dfftp%nl(1)) )
  !
  !CALL mp_sum( v_of_0, intra_bgrp_comm )
  !
  ! ... aux = potential in G-space.
  !
  DEALLOCATE( aux )
  !
  !
  RETURN
  !
END SUBROUTINE setlocal

