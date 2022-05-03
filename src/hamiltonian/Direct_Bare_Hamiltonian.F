!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): FP, DS
!
! This file is distributed under the terms of the GNU
! General Public License. You can redistribute it and/or
! modify it under the terms of the GNU General Public
! License as published by the Free Software Foundation;
! either version 2, or (at your option) any later version.
!
! This program is distributed in the hope that it will
! be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU General Public
! License along with this program; if not, write to the Free
! Software Foundation, Inc., 59 Temple Place - Suite 330,Boston,
! MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
!
subroutine Direct_Bare_Hamiltonian(E,Xk,k)
 !
 use pars,           ONLY:SP,cZERO
 use electrons,      ONLY:levels,n_sp_pol
 use QP_m,           ONLY:QP_nk
 use interfaces,     ONLY:WF_load
 use hamiltonian,    ONLY:Hzero,H_ref_bands,WF_G_max,WF_Go_indx,Hext_nl,Hext_loc!,E_reference
 use parallel_m,     ONLY:PAR_IND_WF_k
 use parallel_int,   ONLY:PP_redux_wait
 use wave_func,      ONLY:WF,wf_ng
 !use gvect,          ONLY:ngm
 use R_lattice,      ONLY:bz_samp,nkibz
 use timing_m,       ONLY:timing
 use pseudo,         ONLY:becp,pp_is_uspp,qe_atoms_map
 use uspp,           ONLY:nkb,ikb_i=>indv_ijkb0,deeq!,deeq_nc
 use uspp_param,     ONLY:nh
 use ions_base,      ONLY:nat,ityp,ntyp=>nsp
 use vlocal,         ONLY:vloc_Gfull
 use wrapper,        ONLY:Vstar_dot_VV
#if defined _RT
 use drivers,        ONLY:l_real_time
#endif
#include<memory.h>
 !
 type(levels)             :: E ! This is so far not needed by the subroutine 
 type(bz_samp)            :: Xk,k
 !
 ! Work space
 !
 integer     :: ik,ikb,ib1,ib2,i_sp_pol,WFbands(2),na,nt,i1,i2,ia_qe,ia_per_type
 integer     :: ib1_wf,ib2_wf,ig
 !
 Hzero   =cZERO
 Hext_nl =cZERO
 Hext_loc=cZERO
 !
 call timing('Direct_Bare_Hamiltonian',OPR='start')
 !
 if (pp_is_uspp) call error("[PPs] Currently, USPP not implemented for ED")
 !
 ! WF Loading
 !
 ! (NB: this turns pp_is_uspp to .true., later on it is restored to .false.)
 !=====================================================================
 WFbands=(/H_ref_bands(1),max(H_ref_bands(2),maxval(E%nbm))/)
 call WF_load(WF,WF_G_max,WF_Go_indx,WFbands,(/1,nkibz/),space='G',title='-SC')
 !write (*,*) shape(WF%c),ngm,wf_ng
 !
 ! Calculation of Hzero
 !
 ! Hzero = T + Hext = -P**2./2. + V_ext^loc + V_ext^nl
 !
 ! Hzero_n,m = P**2/2. + sum_G [WF%c(G,n)]* Vloc(G) * WF%c(G,m)  + sum_i [becp%k(i,n)]* D_ij * becp%k(j,m)
 !
 !=======================================================================================================
 !
 do i_sp_pol=1,n_sp_pol
   !
   do ik=1,QP_nk
     !
     write (*,*) "H_pseudo (ik,i_sp_pol) ",ik,i_sp_pol
     !
     Hext_loc=cZERO
     Hext_nl=cZERO
     do ib1=H_ref_bands(1),H_ref_bands(2)
       do ib2=H_ref_bands(1),H_ref_bands(2)
         !
         !
         ! Calculate the local part of Hext ...
         !
         ! This is a sandwich in G-space of vloc(G) with the wave functions
         ! ==================================================================
         ib1_wf = WF%index(ib1,ik,i_sp_pol)
         ib2_wf = WF%index(ib2,ik,i_sp_pol)
         !
         ! [FP] This loop should be made more efficient with a function/wrapper
         do ig=1,wf_ng
           Hext_loc(ib1,ib2) = Hext_loc(ib1,ib2) + &
                        conjg(WF%c(ig,i_sp_pol,ib1_wf))*Vloc_Gfull(ig)*WF%c(ig,i_sp_pol,ib2_wf)
         enddo
         !
         ! Calculate the nonlocal part of Hext ...
         !
         ! (Spin-orbit / Non-collinear cases thus far not implemented.
         !  Use deeq_nc insted of deeq for those)
         !
         ! NB: Unlike for (i) the dipoles or (ii) the Pseudo_Hamiltonian cases,
         !     here we have to *explicitly* use the nondiagonal deeq matrix to
         !     compute the pseudo Hamiltonian matrix elements since deeq is *not*
         !     already diagonalised and incorporated into the KB projectors.
         !
         ! The external loops over ntyp, nat are adapted from add_vuspsi_k() 
         ! in PW/src/add_vuspsi.f90
         ! =================================================================
         ikb = 0
         do nt=1,ntyp 
           !
           if ( nh(nt) == 0 ) cycle
           !
           ia_per_type=0
           do na=1,nat
             !
             if ( ityp(na) == nt ) then
               ia_per_type=ia_per_type+1
               !
               ia_qe=qe_atoms_map(ia_per_type,nt)
               !
               !write(*,*) nt,na,ia_qe,ikb_i(na)
               !
               do i1=1,nh(nt) ! number of projectors for atom na which is of species nt
                 do i2=1,nh(nt) ! deeq is nondiagonal if more than one projector per atom
                   
                   !write(*,*) i1,i2,deeq(i1,i2,ia_qe,i_sp_pol)
                   !
                   Hext_nl(ib1,ib2) = Hext_nl(ib1,ib2) + &
                                      conjg(becp(ik,i_sp_pol)%k(ikb_i(na)+i1,ib1)) * &
                                                  deeq(i1,i2,ia_qe,i_sp_pol)/2._SP * &
                                            becp(ik,i_sp_pol)%k(ikb_i(na)+i2,ib2)
                 enddo
                 ikb = ikb +1  
               enddo
               !
             endif
             !
           enddo
           !
         enddo
         !
         ! ikb and nkb are here just for debugging, may be removed
         if (ikb/=nkb) call error("[Pseudo Ham.] Something wrong in counting KB projectors and deeq")
         !do ikb=1,nkb
           !
           !Hext_nl(ib1,ib2) = Hext_nl(ib1,ib2) + conjg(becp(ik,i_sp_pol)%k(ikb,ib1)) * becp(ik,i_sp_pol)%k(ikb,ib2)
           !
         !enddo
       enddo
     enddo
     !
     do ib1=H_ref_bands(1),H_ref_bands(2)
       write(*,*) Hext_nl(ib1,:)
       !write(*,*) Hext_loc(ib1,:)
     enddo
     !
     ! Fill Hzero
     ! ==============
     !
     Hzero(:,:,ik,i_sp_pol) = Hzero(:,:,ik,i_sp_pol) + Hext_loc(:,:) + Hext_nl(:,:)
     !
#if defined _RT
     if(l_real_time) then
       if (.not.PAR_IND_WF_k%element_1D(ik)) cycle
     endif
#endif
     !
   enddo
   !
 enddo
 !
 call PP_redux_wait(Hzero)
 pp_is_uspp=.false. 
 ! Here we call again WF_load in default (real) space for compatibility with SC driver 
 call WF_load(WF,WF_G_max,WF_Go_indx,WFbands,(/1,nkibz/),title='-SC')
 !
 write(*,*) "End of Direct hamiltonian "
 !
 call timing('Direct_Bare_Hamiltonian',OPR='stop')
 !
end subroutine Direct_Bare_Hamiltonian
