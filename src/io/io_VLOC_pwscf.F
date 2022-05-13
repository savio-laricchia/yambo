!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): FP
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
integer function io_VLOC_pwscf(ID)
 !
 use pars,             ONLY: SP,LP
 use stderr,           ONLY: intc
 use IO_m,             ONLY: io_sec,read_is_on,write_is_on
 use IO_int,           ONLY: io_connect,io_disconnect, &
&                            io_elemental,io_bulk,io_header
 use qe_pseudo_m,      ONLY: nsp
 use gvect,            ONLY: ngm,ngl
 use vlocal,           ONLY: vloc_yambo
 !
 implicit none
 !
 integer,intent(in) :: ID
 !
 ! Work Space
 !
 io_VLOC_pwscf=io_connect(desc="vloc_pp_pwscf",type=0,ID=ID)
 if (io_VLOC_pwscf/=0) goto 1
 !
 ! Dimensions
 !
 if (any((/io_sec(ID,:)==1/))) then
   !
   call io_elemental(ID,VAR="PARS",VAR_SZ=3,MENU=0)
   call io_elemental(ID,VAR="Number of G-vectors   :",I0=ngm,CHECK=.true.,OP=(/"=="/))
   call io_elemental(ID,VAR="Number of G-shells    :",I0=ngl,CHECK=.true.,OP=(/"=="/))
   call io_elemental(ID,VAR="Number of at. species :",I0=nsp,CHECK=.true.,OP=(/"=="/))
   call io_elemental(ID,VAR="",VAR_SZ=0,MENU=1)
   !
 endif
 !
 ! g-vector variables
 call io_bulk(ID,VAR="ngm",VAR_SZ=(/1/))
 call io_bulk(ID,I0=ngm)
 call io_bulk(ID,VAR="ngl",VAR_SZ=(/1/))
 call io_bulk(ID,I0=ngl)
 call io_bulk(ID,VAR="nsp",VAR_SZ=(/1/))
 call io_bulk(ID,I0=nsp)
 !
 ! vloc
 call io_bulk(ID,VAR="vloc",VAR_SZ=shape(vloc_yambo))
 call io_bulk(ID,D2=vloc_yambo)
 !
1 call io_disconnect(ID=ID)
 !
end function io_VLOC_pwscf
