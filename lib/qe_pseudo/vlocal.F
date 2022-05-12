!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AF, IM
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
module vlocal
   use kinds, only : DP
   real(DP), allocatable    :: vloc(:,:)
   complex(DP), allocatable :: vloc_full(:)
   complex(DP), allocatable :: strf(:,:)

   contains

     subroutine vloc_alloc()
       !
       ! Allocate Vlocal-related quantities
       !
       ! NB: strf and vloc also allocated independently in p2y step
       !
       use gvect,      ONLY: ngl,ngm,gshells
       use ions_base,  ONLY: nsp
       use pseudo,     ONLY: l_use_p2_becp
       !
       ! Allocate vloc if reconstructing bare Hamiltonian...
       !
       if (l_use_p2_becp) then
         !
         if (.not.allocated(vloc)) then
           call gshells()
           allocate(vloc(ngl,nsp))
         endif
         !
         if (.not.allocated(vloc_full)) allocate(vloc_full(ngm)) 
         !
       endif
       !
       ! ... otherwise just allocate structure factors
       !
       if (.not.allocated(strf)) allocate(strf(ngm,nsp))
       !
     end subroutine vloc_alloc
     !
     subroutine vloc_dealloc()
       !
       ! Dellocate Vlocal-related quantities
       !
       if (allocated(vloc))      deallocate(vloc) 
       if (allocated(vloc_full)) deallocate(vloc_full)
       if (allocated(strf))      deallocate(strf)
       !
     end subroutine vloc_dealloc
     !
end module
