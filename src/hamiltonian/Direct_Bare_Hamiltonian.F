!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): DV DS
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
 use pars,           ONLY:cZERO
 use electrons,      ONLY:levels,n_sp_pol,spin,n_spin
 use QP_m,           ONLY:QP_nk,QP_Vxc,QP_n_states
 use hamiltonian,    ONLY:Hzero,E_reference,V_hartree_sc,rho_n,magn_n,V_xc_sc,&
&                         H_ref_bands,WF_G_max,WF_Go_indx,Hext_nl
 use parallel_m,     ONLY:PAR_IND_WF_k,PAR_IND_WF_linear
 use parallel_int,   ONLY:PP_redux_wait
 use wave_func,      ONLY:WF
 use R_lattice,      ONLY:bz_samp,nkibz
 use xc_functionals, ONLY:V_xc
 use interfaces,     ONLY:el_density_and_current,el_magnetization,WF_load,WF_free
 use H_interfaces,   ONLY:V_real_space_to_H
 use timing_m,       ONLY:timing
 use DIPOLES,        ONLY:DIPOLE_t
 use pseudo,         ONLY:becp
 use uspp,           ONLY:nkb
#if defined _SC
 use drivers,        ONLY:l_sc_run
#endif
#if defined _RT
 use drivers,        ONLY:l_real_time
 use real_time,      ONLY:REF_V_xc_sc,REF_V_hartree_sc,rho_reference,magn_reference
#endif
#if defined _SC
 use drivers,        ONLY:l_sc_magnetic
#endif
#include<memory.h>
 !
 type(levels)             :: E
 type(bz_samp)            :: Xk,k
 !
 ! Work space
 !
 integer     :: ik,ikb,ib1,ib2,i_sp_pol,WFbands(2)
 !
 Hzero=cZERO
 Hext_nl=cZERO
 !
 call timing('Direct_Bare_Hamiltonian',OPR='start')
 !
 ! WF Loading, density and magnetization
 !=======================================
 WFbands=(/H_ref_bands(1),max(H_ref_bands(2),maxval(E%nbm))/)
 !call WF_load(WF,WF_G_max,WF_Go_indx,WFbands,(/1,nkibz/),title='-SC')
 !
 ! The Hzero = -nabla**2. + V_ext is computed directly:
 !
 ! Hzero_n,m = Dip**2/2. + sum_i [becp%k(i,n)]* * becp%k(i,m)
 !
 !========================================
 !
 call WF_load(WF,WF_G_max,WF_Go_indx,WFbands,(/1,nkibz/),title='-SC')
 !call get_Hext_nl() 
 if (allocated(becp)) then
   write (*,*) "BECP ALLOCATED!",shape(becp)
 else
   write (*,*) "BECP NOOOOT ALLOCATED!"
 endif
 !write (*,*) "DIPOLES:", shape(Dip)
 do i_sp_pol=1,n_sp_pol
   !
   do ik=1,QP_nk
     !
     write (*,*) "becp", allocated(becp(ik,i_sp_pol)%k),shape(becp(ik,i_sp_pol)%k)
     !
     ! Calculate the nonlocal part of Hext ...
     ! =========================================
     !
     Hext_nl=cZERO
     do ib1=H_ref_bands(1),H_ref_bands(2)
       do ib2=H_ref_bands(1),H_ref_bands(2)
         Hext_nl = dot_product( conjg(becp(ik,i_sp_pol)%k(:,ib1)), becp(ik,i_sp_pol)%k(:,ib2) )
         write(*,*) Hext_nl(H_ref_bands(1),H_ref_bands(1))
         !do ikb=1,nkb
         !  !
         !  Hext_nl = Hext_nl + conjg(becp(ik,i_sp_pol)%k(ikb,ib1)) * becp(ik,i_sp_pol)%k(ikb,ib2)
         !  !
         !enddo
       enddo
     enddo
     !
     ! Fill Hzero
     ! ==============
     !
     Hzero(:,:,ik,i_sp_pol) = Hzero(:,:,ik,i_sp_pol) + Hext_nl(:,:)
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
 !
 !call WF_free(WF)
 !
 write(*,*) "End of Direct hamiltonina "
 !
 call timing('Direct_Bare_Hamiltonian',OPR='stop')
 !
end subroutine Direct_Bare_Hamiltonian

subroutine get_Hext_nl()
 !
 ! Calculate non-local projectors ("becp") for current pseudopotential.
 ! This will be used to construct nonlocal part of Hext.
 !
 ! Flow is taken from WF_load using ultrasoft-specific parts:
 !   PP_uspp_init -> init_us_2 -> PP_compute_becp
 !
 use pars,          ONLY:SP
 use vec_operate,   ONLY:c2a
 use R_lattice,     ONLY:nkibz,k_pt,b
 use qe_pseudo_m,   ONLY:vkb,tpiba
 use wave_func,     ONLY:WF,wf_nc_k,wf_igk,wf_nb_to_load,wf_b_indx
 use pseudo,        ONLY:becp
 use electrons,     ONLY:n_sp_pol
 use hamiltonian,   ONLY:Hext_nl
 !
 ! Work space
 ! 
 integer                  :: i_sp_pol,ikibz
 real(SP)                 :: xk(3)
 complex(SP), allocatable :: wf_disk(:,:,:)
 !
 call PP_uspp_init()
 !
 do i_sp_pol=1,n_sp_pol
   !
   if (any( (/i_sp_pol<WF%sp_pol(1), i_sp_pol>WF%sp_pol(2)/) )) cycle
   !
   do ikibz=1,nkibz
     !
     if (any( (/ikibz<WF%k(1), ikibz>WF%k(2)/) )) cycle
     !
     call c2a(b,k_pt(ikibz,:),xk,mode="ki2c")
     xk=xk/real(tpiba,SP)
     !
     npwk=wf_nc_k(ikibz)
     call init_us_2(npwk,wf_igk(1:npwk,ikibz),xk,vkb)
     !
     call PP_compute_becp( becp(ikibz,i_sp_pol), wf_nc_k(ikibz), wf_nb_to_load,wf_disk, wf_b_indx)
     !
   enddo
   !
 enddo

end subroutine get_Hext_nl
