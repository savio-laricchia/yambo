!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AF, IM, FP
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
subroutine PP_uspp_init()
 !
 ! checks and reads USPP DB, and, in case, initialize
 ! general USPP vars
 !
 use drivers,       ONLY:l_HF_and_locXC,l_col_cut,&
&                        list_dbs,l_setup,l_rim
 use pars,          ONLY:SP,DP
 use vec_operate,   ONLY:c2a
 use D_lattice,     ONLY:nsym,dl_sop,atom_mapper
 use R_lattice,     ONLY:b,g_vec,ng_vec
 use pseudo,        ONLY:pp_is_uspp,qe_pseudo_alloc,PP_uspp_free,qe_atoms_map,&
&                        l_use_p2_becp,PP_uspp_alloc
 use IO_int,        ONLY:io_control
 use IO_m,          ONLY:OP_RD_CL,REP
 !
 use qe_pseudo_m,   ONLY:qe_alat=>alat,qe_tpiba=>tpiba,qe_nat=>nat,qe_tau=>tau,qe_nsp=>nsp,qe_ityp=>ityp
 use qe_pseudo_m,   ONLY:nkb,qe_nsym=>nsym,d1,d2,d3
 use gvect,         ONLY:qe_eigts1=>eigts1, qe_eigts2=>eigts2, qe_eigts3=>eigts3,&
&                        qe_g=>g, qe_gg=>gg, qe_ngm=>ngm, qe_mill=>mill
 use timing_m,      ONLY:timing
 use parser_m,       ONLY:parser
 !
#include<memory.h>
 !
 ! Work Space 
 !
 integer  :: io_err,ID
 integer  :: gv_min,gv_max,i,ig,nfft(3)
 real(SP) :: v1(3)
 real(DP) :: qe_b(3,3)
 complex(DP), allocatable :: strf(:,:)
 logical :: l_will_need_becp
 !
 integer, external :: io_USPP_pwscf

 call parser('BareHfromBecp', l_will_need_becp)
 !
 ! in case, assume qe_pseudo already alloc and init
 !
 if (qe_pseudo_alloc) return
 call timing("PP_uspp_init","start")
 ! [FP]
 if (l_use_p2_becp) then
   call PP_uspp_alloc()
 endif
 ! [FP] 
 !
 ! perform main data IO
 !
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1/),ID=ID)
 io_err=io_USPP_pwscf(ID)
 !
 pp_is_uspp=(io_err==0)
 !
 ! [FP] 
#ifndef _USPP
 if (pp_is_uspp.and..not.l_will_need_becp) then
   call error("[PPs] Ultrasoft PP not supported. [Is ns.uspp_pp_pwscf present but should not?]")
 endif
#endif
 !
 if ((pp_is_uspp.and.l_will_need_becp.and..not.l_use_p2_becp).or.&
    (.not.pp_is_uspp.or.nkb<=0)) then
 !if (.not.pp_is_uspp.or.nkb<=0) then
   ! [FP]
   pp_is_uspp=.false.
   call PP_uspp_free()
   call timing("PP_uspp_init","stop")
   return
 endif
 


 !
 ! checks
 ! [FP]
 if (.not.( l_setup.or.l_col_cut.or.l_rim.or.list_dbs.or. &
            l_HF_and_locXC.or.l_use_p2_becp) ) then
    ! [FP]
    call error("[PPs] USPP not implemented for current runlevel")
 endif

 !
 ! symmetries
 !
 qe_nsym=nsym
 call d_matrix(nsym,real(dl_sop,DP),d1,d2,d3)

 !
 ! atom mapping
 !
 call atom_mapper(qe_nat,real(qe_tau*qe_alat,SP),"cc",qe_atoms_map)
 
 !
 ! init g vectors and structure factors
 !
 qe_ngm=ng_vec
 !

 YAMBO_ALLOC(qe_g,(3,qe_ngm))
 YAMBO_ALLOC(qe_gg,(qe_ngm))
 YAMBO_ALLOC(qe_mill,(3,qe_ngm))

 !
 do ig = 1, ng_vec
    !
    !qe_g(:,ig)=g_vec(ig,:)
    call c2a(b,g_vec(ig,:),v1,'ki2c')
    qe_g(:,ig)=v1/qe_tpiba
    !
    qe_gg(ig)=sum(qe_g(:,ig)**2)
    !
    call c2a(b,g_vec(ig,:),v1,'ki2a')
    qe_mill(:,ig)=nint(v1)
    !
 enddo
 !

 do i = 1, 3
   gv_min=minval(qe_mill(i,:))
   gv_max=maxval(qe_mill(i,:))
   nfft(i)=gv_max-gv_min+1
 enddo
 !
 YAMBO_ALLOC(qe_eigts1,(-nfft(1):nfft(1),qe_nat))
 YAMBO_ALLOC(qe_eigts2,(-nfft(2):nfft(2),qe_nat))
 YAMBO_ALLOC(qe_eigts3,(-nfft(3):nfft(3),qe_nat))
 YAMBO_ALLOC(strf,(qe_ngm,qe_nsp))
 !
 qe_b=transpose(b)/qe_tpiba
 !
 call struc_fact(qe_nat,qe_tau,qe_nsp,qe_ityp,qe_ngm,qe_g,qe_b,&
&                nfft(1),nfft(2),nfft(3), &
&                .false.,strf,.true.,qe_eigts1,qe_eigts2,qe_eigts3)
 YAMBO_FREE(strf)
 !
 ! [FP] Compute D_ij coefficients  
 if (l_use_p2_becp) call PP_compute_deeq()
 ! [FP]
 !
 call timing("PP_uspp_init","stop")
 return
 !
end subroutine PP_uspp_init

subroutine PP_compute_deeq()
 !
 ! This subroutine gets the D-term needed to compute the
 ! nonlocal part of the Hamiltonian together with the KB
 ! projectors.
 !
 ! Adapted from subroutine newd() in PW/src/newd.f90 in QE (ovkan=.false. case).
 !
 ! USPP part (Q function) NOT implemented.
 !
 use pars,                   ONLY:SP,DP
 use pseudo,                 ONLY:l_use_p2_becp
 ! data already stored by p2y
 use uspp,                   ONLY:deeq,deeq_nc,dvan,dvan_so
 use uspp_param,             ONLY:nh
 use ions_base,              ONLY:nat,ityp
 use lsda_mod,               ONLY:nspin
 use noncollin_module,       ONLY: noncolin
 ! qe_pseudo internal modules
 use spin_orb,               ONLY:qe_lspinorb=>lspinorb
 !
 ! Work Space
 !
 integer :: nt, na, is, nht
 ! counters atom type, atoms, spin, aux
 !
 do na = 1, nat
   !
   nt  = ityp(na)
   nht = nh(nt)
   !
   if ( qe_lspinorb ) then
     !
     deeq_nc(1:nht,1:nht,na,1:nspin) = dvan_so(1:nht,1:nht,1:nspin,nt)
     !
   elseif ( noncolin ) then
     !
     deeq_nc(1:nht,1:nht,na,1) = dvan(1:nht,1:nht,nt)
     deeq_nc(1:nht,1:nht,na,2) = ( 0.D0, 0.D0 )
     deeq_nc(1:nht,1:nht,na,3) = ( 0.D0, 0.D0 )
     deeq_nc(1:nht,1:nht,na,4) = dvan(1:nht,1:nht,nt)
     !
   else
     !
     do is = 1, nspin
       !
       deeq(1:nht,1:nht,na,is) = dvan(1:nht,1:nht,nt)
       !
     enddo
     !
   endif
   !
 enddo
 !
end subroutine
