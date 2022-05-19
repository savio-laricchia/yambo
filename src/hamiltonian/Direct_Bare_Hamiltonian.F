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
subroutine Direct_Bare_Hamiltonian(E,Xk,k,Dip)
 !
 use pars,           ONLY:SP,cZERO
 use electrons,      ONLY:levels,n_sp_pol,l_spin_orbit
 use dipoles,        ONLY:DIPOLE_t,compute_P2_dipoles,P_square
 use QP_m,           ONLY:QP_nk
 use interfaces,     ONLY:WF_load
 use hamiltonian,    ONLY:Hzero,H_ref_bands,WF_G_max,WF_Go_indx,Hext_nl,Hext_loc!,E_reference
 use parallel_m,     ONLY:PAR_IND_WF_k
 use parallel_int,   ONLY:PP_redux_wait
 use wave_func,      ONLY:WF
 use R_lattice,      ONLY:bz_samp,nkibz
 use timing_m,       ONLY:timing
 use pseudo,         ONLY:becp,pp_is_uspp,qe_atoms_map
 use uspp,           ONLY:nkb,ikb_i=>indv_ijkb0,deeq,deeq_nc
 use uspp_param,     ONLY:nh
 use ions_base,      ONLY:nat,ityp,ntyp=>nsp
 use vlocal,         ONLY:vloc_dealloc
 use wrapper,        ONLY:Vstar_dot_VV
 use fft_m,          ONLY:fft_size
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:NONE,OP_RD_CL,REP,DUMP
 use gvect,          ONLY:ngm
 use stderr,         ONLY:intc
#if defined _RT
 use drivers,        ONLY:l_real_time
#endif
#include<dev_defs.h>
#include<memory.h>
 !
 type(levels)             :: E ! This is so far not needed by the subroutine 
 type(bz_samp)            :: Xk,k
 type(DIPOLE_t)           :: Dip
 !
 ! Work space
 !
 integer     :: ik,ikb,ib1,ib2,i_sp_pol,WFbands(2),ia,it,i1,i2,ia_qe,ia_per_type
 integer     :: ib1_wf,ib2_wf,io_DIP_err,ID_DIP
 complex(SP), allocatable :: vloc_full_r(:)
 !
 integer, external :: io_DIPOLES
 !
 ! Debug
 !
 open(1, file='TEST_Hbare_real.dat')
 open(2, file='TEST_Hbare_imag.dat')
 !
 Hzero   =cZERO
 Hext_nl =cZERO
 Hext_loc=cZERO
 allocate(vloc_full_r(fft_size))
 vloc_full_r=cZERO
 !
 call timing('Direct_Bare_Hamiltonian',OPR='start')
 !
 if (pp_is_uspp) call error("[PPs] Currently, USPP not implemented for ED")
 !
 ! Read dipoles
 !
 ! NB: this is a workaround. In the final version this should be done externally
 !=====================================================================
 if (.not.compute_P2_dipoles) call error("[DIP] P_square not activated")
 call DIPOLE_dimensions(E,Dip,H_ref_bands,(/0._SP,0._SP,0._SP/))
 ! I use this double call to io_DIPOLES since DIPOLE_IO does not work
 !
 call io_control(ACTION=OP_RD_CL,COM=REP,MODE=DUMP,SEC=(/0/),ID=ID_DIP)
 io_DIP_err = io_DIPOLES(Dip,ID_DIP)
 call io_control(ACTION=OP_RD_CL,COM=NONE,MODE=DUMP,SEC=(/2/),ID=ID_DIP)
 io_DIP_err = io_DIPOLES(Dip,ID_DIP)
 ! The four lines above could be replaced with this one but it gives
 ! an error because PAR_K_scheme%ibz_index ends up not being allocated 
 !
 !call DIPOLE_IO(k,E,Dip,'read',io_DIP_err,"SC")
 if (.not.allocated(P_square)) call error("[DIP] P_square not found")
 !
 ! Load WF, compute Vloc and becp and setup FFT
 !
 ! (NB:  this turns pp_is_uspp to .true., later on it is restored to .false.)
 ! (NB2: this calls vloc_alloc(), it is deallocated below) 
 !=====================================================================
 WFbands=(/H_ref_bands(1),max(H_ref_bands(2),maxval(E%nbm))/)
 call WF_load(WF,WF_G_max,WF_Go_indx,WFbands,(/1,nkibz/),title='-SC')
 if (WF_G_max/=ngm) call error('[FFT] FFT needs ng='//trim(intc(ngm))//' but has ng='//trim(intc(WF_G_max))//'')
 !
 ! FFT of vloc into real space
 !
 ! (NB: grid dimensions are determined by the FFT_setup previously run by WF_load)
 !=====================================================================
 call FFT_vloc_G_to_R()
 !
 ! Calculation of Hzero
 !
 ! Hzero = T + Hext = -P**2/2. + V_ext^loc + V_ext^nl
 !
 ! Hzero_n,m = -P_square(n,m)/2. + int_r [WF%c(r,n)]* Vloc(r) * WF%c(r,m)  + sum_ij [becp%k(i,n)]* D_ij * becp%k(j,m)
 !
 !=======================================================================================================
 !
 do i_sp_pol=1,n_sp_pol
   !
   do ik=1,QP_nk
     !
     write (*,*) "H_zero (ik,i_sp_pol) ",ik,i_sp_pol
     !
     Hext_loc=cZERO
     Hext_nl=cZERO
     do ib1=H_ref_bands(1),H_ref_bands(2)
       do ib2=H_ref_bands(1),H_ref_bands(2)
         !
         !
         ! Calculate the local part of Hext ...
         !
         ! This is a sandwich in R-space of vloc(r) with the wave functions
         ! ==================================================================
         ib1_wf = WF%index(ib1,ik,i_sp_pol)
         ib2_wf = WF%index(ib2,ik,i_sp_pol)
         !
         Hext_loc(ib1,ib2) = Hext_loc(ib1,ib2) + &
                             Vstar_dot_VV(fft_size,WF%c(:,i_sp_pol,ib1_wf),vloc_full_r,WF%c(:,i_sp_pol,ib2_wf))
         !
         ! Calculate the nonlocal part of Hext ...
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
         do it=1,ntyp 
           !
           if ( nh(it) == 0 ) cycle
           !
           ia_per_type=0
           do ia=1,nat
             !
             if ( ityp(ia) == it ) then
               ia_per_type=ia_per_type+1
               !
               ia_qe=qe_atoms_map(ia_per_type,it)
               !
               !write(*,*) it,na,ia_qe,ikb_i(ia)
               !
               do i1=1,nh(it) ! number of projectors for atom ia which is of species it
                 do i2=1,nh(it) ! deeq is nondiagonal if more than one projector per atom
                   
                   !write(*,*) i1,i2,deeq(i1,i2,ia_qe,i_sp_pol)
                   !
                   if (l_spin_orbit) then
                     !
                     Hext_nl(ib1,ib2) = Hext_nl(ib1,ib2) + &
                                       conjg(becp(ik,i_sp_pol)%k(ikb_i(ia)+i1,ib1)) * &
                                                 deeq_nc(i1,i2,ia_qe,i_sp_pol)/2._SP * &
                                              becp(ik,i_sp_pol)%k(ikb_i(ia)+i2,ib2)
                     !
                   else
                     !
                     Hext_nl(ib1,ib2) = Hext_nl(ib1,ib2) + &
                                        conjg(becp(ik,i_sp_pol)%k(ikb_i(ia)+i1,ib1))* &
                                                    deeq(i1,i2,ia_qe,i_sp_pol)/2._SP* &
                                              becp(ik,i_sp_pol)%k(ikb_i(ia)+i2,ib2)
                     !
                   endif
                   !  
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
     !do ib1=H_ref_bands(1),H_ref_bands(2)
       !write(*,*) Hext_nl(ib1,:)
       !write(*,*) Hext_loc(ib1,:)
     !enddo
     !
     ! Fill Hzero
     ! ==============
     !
     Hzero(:,:,ik,i_sp_pol) = Hzero(:,:,ik,i_sp_pol) & 
                              - P_square(:,:,ik,i_sp_pol)/2._SP + Hext_loc(:,:) + Hext_nl(:,:)
     !
     write(1,*) "#Re H_bare (ik,i_sp_pol) ",ik,i_sp_pol
     write(2,*) "#Im H_bare (ik,i_sp_pol) ",ik,i_sp_pol
     do ib1=H_ref_bands(1),H_ref_bands(2)
       write(*,*) Hzero(ib1,:,ik,i_sp_pol)
       write(1,*) real(Hzero(ib1,:,ik,i_sp_pol))
       write(2,*) aimag(Hzero(ib1,:,ik,i_sp_pol))
     enddo
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
 call vloc_dealloc() ! Dealloc here since we couldn't do it in PP_uspp_init 
 deallocate(vloc_full_r)
 !
 close(1)
 close(2)
 write(*,*) "End of Direct hamiltonian "
 !
 call timing('Direct_Bare_Hamiltonian',OPR='stop')
 !
 contains
   !
   subroutine FFT_vloc_G_to_R()
    !
    use pars,       ONLY:DP  
    use vlocal,     ONLY:vloc_full
    use fft_m,      ONLY:fft_dim,fftw_plan,fft_g_table
    use gvect,      ONLY:ngm
    !
    ! Work space
    !
    integer                   :: ig
    complex(DP), allocatable  :: vloc_aux(:)
    !
    allocate(vloc_aux(fft_size))
    !
    !$omp parallel do default(shared), private(ig)
    do ig=1,ngm
      vloc_aux(fft_g_table(ig,1))=vloc_full(ig)
    enddo
    !$omp end parallel do
    !
#if defined _FFTW
    call dfftw_destroy_plan(fftw_plan)
    fftw_plan = 0
    call fft_3d(vloc_aux,fft_dim,1,fftw_plan)
#else
    call fft_3d(vloc_aux,fft_dim,1)
#endif
    !
    ! [FP] BUG: vloc(r) actually has a non-negligible imaginary part!!!!
    !
    !write (*,*) "Vloc(r):"
    !do ig=1,fft_size
    !  write (*,*) vloc_aux(ig)
    !enddo
    !write (*,*) "End Vloc(r)"
    !
    ! Convert to single precision and real but keep complex type
    ! even though vloc(r) is real, as required by Vstar_dot_VV
    vloc_full_r = real(cmplx(vloc_aux))
    !
    deallocate(vloc_aux)
    !
    end subroutine FFT_vloc_G_to_R
    !
end subroutine Direct_Bare_Hamiltonian
