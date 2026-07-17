# AMPS Lookup-Table Data: Extraction for Python Conversion

## 1. Data files on disk

### `/Users/jcanton/projects/scale_amps/scale-rm/test/case/cloudlab/AMPS_DATA/`

Three subdirectories, matching the three namelist directory variables `DRCETB` (collision_data), `DRAPTB` (apact), `DRSTTB` (statpack) declared in `com_amps.F90:74-75` and read via namelist in `mod_amps_utility.F90:1384-1386`. Readers are called from `mod_amps_lib.F90:626` (RDCETB), `:632` (RDAPTB), `:638` (RDSTTB).

```
AMPS_DATA/apact/
  ap_act.bin.dat            7987272 bytes   (Fortran unformatted binary)

AMPS_DATA/collision_data/
  bbc_drop_mat.dat            10182
  bbc_drop_Nre_Ec.dat         30957   * read
  bbc_drop.dat                 4089
  col_drop_mat.dat            10178
  col_drop_Nre_Ec.dat         33818   * read
  col_drop.dat                 6212
  col_frq.dat                 77326   * read
  drop_col_mat.dat              516
  drop_col.dat                  754
  drop_drop_Rey3.dat          82990
  drop_drop_Rey4.dat         565774   * read
  drop_drop.dat               89238
  drop_pla_mat.dat              596
  drop_pla.dat                  862
  grap01_drop.dat              8585
  grap04_drop.dat              4935
  grap08_drop.dat              3758
  grp01_ratNre_Ec.dat         26868   * read
  grp04_ratNre_Ec.dat         19674   * read
  grp08_ratNre_Ec.dat         14522   * read
  hex_drop_mat.dat            10169
  hex_drop_Nre_Ec.dat         30839   * read
  hex_drop.dat                 5749
  lmt_mass_col.dat             2020   * read
  lmt_mass_pla.dat             2072   * read
  lmt_mass.dat                 1188
  mta_map.dat                 52723
  mtc_map.dat                 42196
  pla_frq.dat                 77324   * read
  pol_frq.dat                 56717   * read
  ppo_frq.dat                 77328   * read
  ros_frq.dat                 77325   * read
  seed_input.org                448
  tmd_map_col_r.dat           68382
  tmd_map_col.dat             54732   * read
  tmd_map_pla_r.dat           63830
  tmd_map_pla.dat             54730   * read
  tmd_map.dat                 38632
  tmp_map_col.dat             68380   * read
  tmp_map_pla.dat             68378   * read
  tmp_map.dat                 52671

AMPS_DATA/statpack/
  stdnorm.dat                 44183   * read  (451 lines)
```

Files marked `* read` are the only ones RDCETB/RDSTTB actually open (18 files). The rest (`*_mat.dat`, `drop_drop.dat`, `drop_drop_Rey3.dat`, `grap*_drop.dat`, `*_map.dat` without `_col`/`_pla`, `_r` variants, `lmt_mass.dat`, `seed_input.org`) are legacy/unused by the current code. `ap_act.bin.dat` is opened by RDAPTB but that subroutine is **dead code** (see below).

`/Users/jcanton/projects/scale_amps/data/mp/AMPS/` contains only `AMPSTASK.F` (8409 bytes, namelist sample) and `seedini_0m` (70 bytes) — **no lookup tables there**.

## 2. Reader subroutines (verbatim)

### 2.1 RDCETB — `mod_amps_utility.F90:96-294` (complete)

```fortran
      subroutine RDCETB(nbr,ncr,nbi,nci)
      use com_amps
      implicit none
      integer, intent(in) :: nbr,nbi,ncr,nci
!      character*100   DRCETB
      character(len=110) :: ifname
      character(len=300) :: header
      real(PS) :: dum1,dum2,dum3
      integer :: i,j,k,l,nrow,ncol,npar,nt
      integer :: io_rdcetb
!     +++++++++++ read drop_drop.dat +++++++++++++++++++++++++++++
      io_rdcetb = IO_get_available_fid()
      open(io_rdcetb, file=trim(DRCETB)//"/drop_drop_Rey4.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_drpdrp,nc_drpdrp,xs_drpdrp,dx_drpdrp &
           ,ys_drpdrp,dy_drpdrp
      do i = 1, nr_drpdrp
         read(io_rdcetb,*) (drpdrp(i,j), j = 1, nc_drpdrp)
      end do
      close(io_rdcetb)

!     +++++++++++ read hex_drop_Nre_Ec.dat +++++++++++++++++++++++++++++
      open(io_rdcetb, file=trim(DRCETB)//"/hex_drop_Nre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_hexdrp,nc_hexdrp,xs_hexdrp,dx_hexdrp &
           ,ys_hexdrp,dy_hexdrp
      do i = 1, nr_hexdrp
         read(io_rdcetb,*) (hexdrp(i,j), j = 1, nc_hexdrp)
      end do
      close(io_rdcetb)

!     +++++++++++ read bbc_drop_Nre_Ec.dat +++++++++++++++++++++++++++++
      open(io_rdcetb, file=trim(DRCETB)//"/bbc_drop_Nre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_bbcdrp,nc_bbcdrp,xs_bbcdrp,dx_bbcdrp &
           ,ys_bbcdrp,dy_bbcdrp
      do i = 1, nr_bbcdrp
         read(io_rdcetb,*) (bbcdrp(i,j), j = 1, nc_bbcdrp)
      end do
      close(io_rdcetb)

!     +++++++++++ read col_drop_Nre_Ec.dat +++++++++++++++++++++++++++++
      open(io_rdcetb, file=trim(DRCETB)//"/col_drop_Nre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_coldrp,nc_coldrp,xs_coldrp,dx_coldrp &
           ,ys_coldrp,dy_coldrp
      do i = 1, nr_coldrp
         read(io_rdcetb,*) (coldrp(i,j), j = 1, nc_coldrp)
      end do
      close(io_rdcetb)

!     +++++++++++ read grp01_ratNre_Ec.dat +++++++++++++++++++++++++++++
!ccc      ifname=DRCETB(1:LEN(DRCETB))//"/grap01_drop.dat"
      open(io_rdcetb, file=trim(DRCETB)//"/grp01_ratNre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_gp1drp,nc_gp1drp,xs_gp1drp,dx_gp1drp &
           ,ys_gp1drp,dy_gp1drp
      do i = 1, nr_gp1drp
         read(io_rdcetb,*) (gp1drp(i,j), j = 1, nc_gp1drp)
      end do
      close(io_rdcetb)

!     +++++++++++ read grp04_ratNre_Ec.dat +++++++++++++++++++++++++++++
!ccc      ifname=DRCETB(1:LEN(DRCETB))//"/grap04_drop.dat"
      open(io_rdcetb, file=trim(DRCETB)//"/grp04_ratNre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_gp4drp,nc_gp4drp,xs_gp4drp,dx_gp4drp &
           ,ys_gp4drp,dy_gp4drp
      do i = 1, nr_gp4drp
         read(io_rdcetb,*) (gp4drp(i,j), j = 1, nc_gp4drp)
      end do
      close(io_rdcetb)

!     +++++++++++ read grp08_ratNre_Ec.dat +++++++++++++++++++++++++++++
!ccc      ifname=DRCETB(1:LEN(DRCETB))//"/grap08_drop.dat"
      open(io_rdcetb, file=trim(DRCETB)//"/grp08_ratNre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_gp8drp,nc_gp8drp,xs_gp8drp,dx_gp8drp &
           ,ys_gp8drp,dy_gp8drp
      do i = 1, nr_gp8drp
         read(io_rdcetb,*) (gp8drp(i,j), j = 1, nc_gp8drp)
      end do
      close(io_rdcetb)


!     +++++++++++++++++ read frequency of habit (T<-20) ++++++
      open(io_rdcetb, file=trim(DRCETB)//"/pol_frq.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (pol_frq(i,j), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/pla_frq.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (pla_frq(i,j), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/col_frq.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (col_frq(i,j), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/ros_frq.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (ros_frq(i,j), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/ppo_frq.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (ppo_frq(i,j), j = 1, ncol)
      end do
      close(io_rdcetb)

!     +++ normalize the frequency +++
      do i=1,nrow
         do j=1,ncol
            pol_frq(i,j)=max(0.0_RP,pol_frq(i,j))
            pla_frq(i,j)=max(0.0_RP,pla_frq(i,j))
            col_frq(i,j)=max(0.0_RP,col_frq(i,j))
            ros_frq(i,j)=max(0.0_RP,ros_frq(i,j))
            ppo_frq(i,j)=max(0.0_RP,ppo_frq(i,j))

            dum1=pol_frq(i,j)+pla_frq(i,j)+col_frq(i,j)
            pol_frq(i,j)=pol_frq(i,j)/dum1
            pla_frq(i,j)=pla_frq(i,j)/dum1
            col_frq(i,j)=col_frq(i,j)/dum1

!            dum2=pol_frq(i,j)+pla_frq(i,j)+col_frq(i,j)
!            write(fid_alog,*) "i,j,bf sum of frq, af",i,j,dum1,dum2
         end do
      end do

!     +++++++++++++++++ read map data for diagnose a and c ++++++
      open(io_rdcetb, file=trim(DRCETB)//"/tmp_map_col.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (mtac_map_col(i,j,1), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/tmd_map_col.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (mtac_map_col(i,j,2), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/lmt_mass_col.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) dum1,lmt_mass_col(i),dum2,dum3
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/tmp_map_pla.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (mtac_map_pla(i,j,1), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/tmd_map_pla.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (mtac_map_pla(i,j,2), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/lmt_mass_pla.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) dum1,lmt_mass_pla(i),dum2,dum3
      end do
      close(io_rdcetb)
      call cal_indx_lmtmass()

      if ( IsMaster .and. debug ) then
        write(fid_alog,*) "bottom of RDCETB"
      end if

      return
      end subroutine RDCETB
```

All reads are list-directed (`read(...,*)`), so a Python parser can just whitespace-split; note `read(io_rdcetb,*) header` on a `#`-comment line consumes exactly one line (list-directed read into a single character variable stops at the first blank-delimited token but still advances the record).

### 2.2 RDAPTB — `mod_amps_utility.F90:296-330` (complete) — **DEAD CODE**

Note the bare `return` at line 302, *before* any I/O: `ap_act.bin.dat` is never read in the current code.

```fortran
      subroutine RDAPTB()
      use com_amps
      implicit none
      character(len=300) :: header
      integer :: i, j, k,l,id
      integer :: io_rdaptb
      return
!     +++++++++++ read ap_act.bin.dat +++++++++++++++++++++++++++++
      io_rdaptb = IO_get_available_fid()
      open(io_rdaptb, file=trim(DRAPTB)//"/ap_act.bin.dat",form='unformatted')
      read(io_rdaptb) (apt_st(i),i=1,4),(apt_add(i),i=1,4),(apt_max(i),i=1,4) &
           ,(napt(i),i=1,4)
      do i=1,napt(1)
         do j=1,napt(2)
            do k=1,napt(3)
               do l=1,napt(4)
                  id=(j-1)*napt(3)*napt(4)+(k-1)*napt(4)+l
                  read(io_rdaptb) frac_apact(1,id),frac_apact(2,id)
               end do
            end do
         end do
         if ( IsMaster .and. debug ) then
           if(i==1+nint((ap_lnsig(1)-apt_st(1))/apt_add(1))) then
              write(fid_alog,*) "ap_lnsig(1), and chosen sigma:",ap_lnsig(1), &
                   apt_st(1)+real(i-1,PS_KIND)*apt_add(1)
           end if
           if(i==1+nint((ap_lnsig(2)-apt_st(1))/apt_add(1))) then
              write(fid_alog,*) "ap_lnsig(2), and chosen sigma:",ap_lnsig(2), &
                   apt_st(1)+real(i-1,PS_KIND)*apt_add(1)
           end if
         end if
      end do
      close(io_rdaptb)
      return
      end subroutine RDAPTB
```

(If it were live, the target is `frac_apact(2,8320)` plus `apt_st(4),apt_add(4),apt_max(4),napt(4)` in `common/RAPTBL`/`common/IAPTBL`, `com_amps.F90:63-70`.)

### 2.3 RDSTTB — `mod_amps_utility.F90:333-348` (complete)

```fortran
      subroutine RDSTTB()
      use com_amps
      implicit none
      real(PS) :: dum1,dum2
      integer i,j
      integer :: io_rdsttb
!     +++++++++++ read standard normal distribution dat +++++++++++++++++++++++++++++
      io_rdsttb = IO_get_available_fid()
      open(io_rdsttb, file=trim(DRSTTB)//"/stdnorm.dat",form='formatted')
      do i=1,451
         read(io_rdsttb,*) dum1,dum2,(znorm(j,i),j=1,4)
      end do
      close(io_rdsttb)

      return
      end subroutine RDSTTB
```

No header; exactly 451 lines of 6 columns; columns 1-2 discarded (col 2 is the abscissa x = 0.00…4.50 in steps of 0.01), columns 3-6 → `znorm(1:4,i)`, i.e. `znorm(4,451)` in `common/STTBL` (`com_amps.F90:87-88`). Semantics (from the commented `getznorm` at `mod_amps_utility.F90:351-367`): `znorm(which,i)` = integral of `z**(which-1) * exp(-z^2/2)` (normalized) from 0 to x.

## 3. Representative file heads (verbatim)

### `collision_data/drop_drop_Rey4.dat` (204 lines = 2 comment + 1 dims + 201 data rows of 201 values)

```
# collision efficiency between drops. Reynolds number of large drop vs ratio of radius
# nrow, ncol, xstart, dx, ystart, dy
201 201 0.004 0.005 -3.9379  0.0286
          0.0           0.0           0.0           0.0           0.0           0.0 ...
          0.0           0.0           0.0           0.0           0.0           0.0 ...
          0.0           0.0           0.0           0.0           0.0           0.0 ...
```

Header line 3 order is `nrow ncol xstart dx ystart dy` → `nr, nc, xs, dx, ys, dy`. Note the "Rey_Ec" family files have **2** comment lines then the dims line; each data row is one long physical line (space- or tab-separated).

### `collision_data/pol_frq.dat` (53 lines = 1 comment + 1 dims + 51 rows of 101 values)

```
# this is a frequency of polycrystals -70<T<-20
 51 101
0.13987755 0.15926646 0.18153173 0.20896764 0.23785380 0.26494684 0.28649564 0.30493480 ...
0.15854845 0.17851412 0.20085279 0.22547584 0.25092999 0.27511609 0.29661720 0.31543388 ...
```

Habit-frequency and map files have only **1** comment line, then `nrow ncol` (no xs/dx/ys/dy — grid is implicit).

### `statpack/stdnorm.dat` (451 lines, no header)

```
   0.0000000E+00   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000
   0.0000000E+00   1.0000000E-02   3.9893563E-03   1.9946615E-05   1.3297677E-07   9.9732247E-10
   0.0000000E+00   2.0000000E-02   7.9783137E-03   7.9780478E-05   1.0637184E-06   1.5955564E-08
   0.0000000E+00   3.0000000E-02   1.1966473E-02   1.7948364E-04   3.5895113E-06   8.0761581E-08
   0.0000000E+00   4.0000000E-02   1.5953437E-02   3.1902620E-04   8.5066847E-06   2.5518693E-07
   0.0000000E+00   5.0000000E-02   1.9938806E-02   4.9836631E-04   1.6610134E-05   6.2282811E-07
```

Also useful (other header variants):

`hex_drop_Nre_Ec.dat` (tab-separated data, 2-comment header):
```
# collision efficiency bet plate and drop. log10(Nre) versus log10(Nre)
# nrow, ncol, xstart, dx, ystart, dy
 64 71 -4.5 0.1 -4.3 0.1
0.0	0.0	0.0	0.0	...
```

`tmp_map_col.dat` (note ncol=91 < array dim 101):
```
# diagnosis of log10(phi=c/a) for columnar hexagonal ice crystal given mass (-12<log10(m)<-3) and temperature (-50<T<-1)
  50  91
-9.9066923e-02 -8.6345993e-02 ...
```

`lmt_mass_col.dat` (4 columns per row; only column 2 is kept → `lmt_mass_col(i)`):
```
# limit mass for diagnosis : columns
  50  4
-50.0 -4.44041 5.81650880E-004 1.59772
-49.0 -4.39537 6.05699080E-004 1.5867
```

## 4. Where each table lands

### com_amps common blocks (`com_amps.F90`)

```fortran
  real(PS), pointer :: bu_fd(:,:),bu_tmass(:)
  common/ECTBL/drpdrp(201,201), bbcdrp(64,71), hexdrp(64,71), &
       coldrp(62,71),&
       gp1drp(37,125), gp4drp(27,125), gp8drp(21,125),rvtm(4,3,15,9),&
       bu_fd, bu_tmass
```
plus per-table header scalars in `common/IECTBLH` (`nr_*, nc_*`) and `common/RECTBLH` (`xs_*, dx_*, ys_*, dy_*`) — `com_amps.F90:18-29`.

| File | com_amps target (declared dims) | aux (from file header) |
|---|---|---|
| drop_drop_Rey4.dat | `drpdrp(201,201)` | nr=201 nc=201 xs=0.004 dx=0.005 ys=-3.9379 dy=0.0286 |
| hex_drop_Nre_Ec.dat | `hexdrp(64,71)` | nr=64 nc=71 xs=-4.5 dx=0.1 ys=-4.3 dy=0.1 |
| bbc_drop_Nre_Ec.dat | `bbcdrp(64,71)` | from file line 3 |
| col_drop_Nre_Ec.dat | `coldrp(62,71)` | from file line 3 |
| grp01_ratNre_Ec.dat | `gp1drp(37,125)` | from file line 3 |
| grp04_ratNre_Ec.dat | `gp4drp(27,125)` | from file line 3 |
| grp08_ratNre_Ec.dat | `gp8drp(21,125)` | from file line 3 |
| pol/pla/col/ros/ppo_frq.dat | `pol_frq(51,101)` etc., `common/FRQTBL` (com_amps.F90:53-55); pol/pla/col renormalized in-place after read | nrow ncol in file (51 101) |
| tmp_map_col.dat | `mtac_map_col(:,:,1)` of `mtac_map_col(50,101,2)`, `common/RMTACTBL` (line 58-59) | nrow=50 ncol=91 |
| tmd_map_col.dat | `mtac_map_col(:,:,2)` | 〃 |
| lmt_mass_col.dat | `lmt_mass_col(50)` (col 2 of 4) | nrow=50 ncol=4 |
| tmp_map_pla.dat / tmd_map_pla.dat / lmt_mass_pla.dat | `mtac_map_pla(50,101,2)`, `lmt_mass_pla(50)` | 〃 |
| stdnorm.dat | `znorm(4,451)`, `common/STTBL` (line 87-88) | none (x = (i-1)*0.01 implicit) |
| ap_act.bin.dat | (dead) `frac_apact(2,8320)` + `apt_st/add/max(4)`, `napt(4)` | binary self-described |

After `lmt_mass_*` read, `cal_indx_lmtmass()` is called (derives `i_lmt_mass_col/pla(50)` index arrays).

### col_lut_aux type — `class_Group.F90:179-183` (verbatim)

```fortran
  type col_lut_aux
    sequence
    real(PS) :: xs,dx,ys,dy
    integer :: nr,nc
  end type col_lut_aux
```

and the IGP holder, `class_Group.F90:186-195`:

```fortran
  integer,parameter :: nok_max=23
  type vap_igp_aux
    sequence
    ! number of knots
    integer :: nok
    integer :: align_vap_igp_aux ! CHIARUI
    ! parameters for cubic approximation
    real(PS),dimension(nok_max)   :: x
    real(PS),dimension(nok_max,4) :: a,b
  end type vap_igp_aux
```

1-D LUT holders, `class_Mass_Bin.F90:281-295`:

```fortran
  integer,parameter :: nx_lut_max=100
  type data1d_lut
    sequence 
    integer :: n
    integer :: align_data1d_lut ! CHIARUI
    real(PS) :: xs,dx
    real(PS),dimension(nx_lut_max) :: y
  end type data1d_lut
  type data1d_lut_big
    sequence 
    integer :: n
    integer :: align_data1d_lut_big ! CHIARUI
    real(PS) :: xs,dx
    real(PS),dimension(501) :: y
  end type data1d_lut_big
```

### Cloud_Micro members — `class_Cloud_Micro.F90:171-224`

```fortran
     ! collisional-breakup variables
     integer :: imin_bk,imax_bk,jmin_bk,jmax_bk
! this is for 30 bins
!     real(PS) :: bu_fd(2,17835),bu_tmass(435)
! this is for 80 bins
     real(PS) :: bu_fd(2,62400),bu_tmass(780)
     ...
     ! collision LUT
     type (col_lut_aux) :: adrpdrp
     real(PS),dimension(201,201) :: drpdrp
     type (col_lut_aux) :: ahexdrp
     real(PS),dimension(64,71) :: hexdrp
     type (col_lut_aux) :: abbcdrp
     real(PS),dimension(64,71) :: bbcdrp
     type (col_lut_aux) :: acoldrp
     real(PS),dimension(62,71) :: coldrp
     type (col_lut_aux) :: agp1drp
     real(PS),dimension(37,125) :: gp1drp
     type (col_lut_aux) :: agp4drp
     real(PS),dimension(27,125) :: gp4drp
     type (col_lut_aux) :: agp8drp
     real(PS),dimension(21,125) :: gp8drp
     ...
     ! Inherent Growth parameterization
     type (vap_igp_aux) :: vigp
     ! Osmotic coefficients
     type (data1d_lut) :: osm_nhs4,osm_sdch
     ...
     ! lookup table for standard normal distribution
     type (data1d_lut_big) :: snrml
     ! lookup table for invserse standard normal distribution
     type (data1d_lut_big) :: isnrml
```

Copied from com_amps in `make_Cloud_Micro_cnfg`, `class_Cloud_Micro.F90:436-519`, e.g.:

```fortran
    i1d_pair_max=(imax_bk-1)-jmin_bk+1+(imax_bk-imin_bk)*(1+imax_bk-imin_bk)/2
    kk_max=i1d_pair_max*NRBIN
    CM%bu_fd(1:2,1:kk_max)=bu_fd(1:2,1:kk_max)
    CM%bu_tmass(1:i1d_pair_max)=bu_tmass(1:i1d_pair_max)
    ...
    CM%adrpdrp=make_col_lut(nr_drpdrp,nc_drpdrp,xs_drpdrp,ys_drpdrp &
                           ,dx_drpdrp,dy_drpdrp)
    CM%drpdrp(1:nr_drpdrp,1:nc_drpdrp)=drpdrp(1:nr_drpdrp,1:nc_drpdrp)
    ... (same pattern for ahexdrp/hexdrp, abbcdrp/bbcdrp, acoldrp/coldrp,
         agp1drp/gp1drp, agp4drp/gp4drp, agp8drp/gp8drp) ...
    CM%vigp%nok=nok_igp
    CM%vigp%x(1:nok_igp)=x_igp(1:nok_igp)
    CM%vigp%a(1:nok_igp,1:4)=a_igp(1:nok_igp,1:4)
    CM%vigp%b(1:nok_igp,1:4)=b_igp(1:nok_igp,1:4)
    CM%osm_nhs4%n=n_osm_nh42so4 ; %dx=dx_osm_nh42so4 ; %xs=xs_osm_nh42so4 ; %y(1:n)=y_osm_nh42so4(1:n)
    CM%osm_sdch%n=n_osm_sodchl  ; ... 
    CM%snrml%n=n_snrml   ; %dx=dx_snrml  ; %xs=xs_snrml  ; %y(1:n)=y_snrml(1:n)
    CM%isnrml%n=n_isnrml ; %dx=dx_isnrml ; %xs=xs_isnrml ; %y(1:n)=y_isnrml(1:n)
```

Note `make_col_lut` argument order is `(nr,nc,xs,ys,dx,dy)`. The habit-frequency (`*_frq`), `mtac_map_*`, `lmt_mass_*`, and `znorm` tables are **not** Cloud_Micro members — they stay in com_amps common blocks and are used from there. Allocation of the com_amps `bu_fd`/`bu_tmass` pointers: `mod_amps_lib.F90:443-447` — `bu_fd(2,17835)`/`bu_tmass(435)` for 30 bins, `bu_fd(2,62400)`/`bu_tmass(780)` for 80 bins.

## 5. Tables computed at init, not read from disk

All in com_amps, filled at startup (called from `mod_amps_lib.F90` init path):

1. **Osmotic coefficient LUTs** — `init_osmo_par` (`mod_amps_utility.F90:12923-12986`): tabulates `osm_ammsul` / `osm_sodchl` (hard-coded 23/24-point piecewise-linear experimental curves, `:12988-13085`) on molality grids 0:5.5:0.1 and 0:6.0:0.1 → `y_osm_nh42so4(100)`, `y_osm_sodchl(100)` with `n_*, xs_*, dx_*`.

```fortran
  subroutine init_osmo_par
  use com_amps
  implicit none

    real(PS) :: x,dx
    real(PS) :: xmin,xmax
    real(PS),dimension(100) :: y
    integer :: i,n

    ! (NH4)2 SO4
    xmin=0.0
    xmax=5.5
    dx=0.1

    n=(xmax-xmin)/dx+1
    i=1
    x=xmin
    do
      if(x>xmax) exit
      y(i)=osm_ammsul(x)
      x=x+dx
      i=i+1
    enddo
    n_osm_nh42so4=n
    xs_osm_nh42so4=xmin
    do i=1,n
      y_osm_nh42so4(i)=y(i)
    enddo
    dx_osm_nh42so4=dx
    ! ... (NaCl block identical with osm_sodchl, xmax=6.0) ...
  end subroutine init_osmo_par
```

Hard-coded curve data (osm_ammsul): `x = (/ 0.0,0.1,0.2,...,5.5 /)` (23 pts), `y = (/ 1.0,0.767,0.731,0.707,0.690,0.677,0.667,0.658,0.652,0.646,0.640,0.632,0.628,0.624,0.623,0.623,0.626,0.635,0.647,0.660,0.673,0.686,0.699 /)`; osm_sodchl: 24 pts, `y = (/ 1.0,0.932,0.925,0.922,0.920,0.921,0.923,0.926,0.929,0.932,0.936,0.943,0.951,0.962,0.972,0.983,1.013,1.045,1.080,1.116,1.153,1.192,1.231,1.271 /)`; flat extrapolation beyond ends, linear interpolation inside.

2. **Normal CDF LUT** — `init_normal_lut` (`:13087-13129`): x = 0:5:0.01 (n=501), `y_snrml(i) = 1 - Phi(x)` via `cdfnor` (DCDFLIB); stores `n_snrml, xs_snrml=0, dx_snrml=0.01, y_snrml(501)`.

3. **Inverse normal CDF LUT** — `init_inv_normal_lut` (`:13131-13187`): log10-spaced probability grid, `ymin=1e-30, ymax=0.5, n=501, dy=(log10(ymax)-log10(ymin))/(n-1)`; `x(i)=dinvnr(y, 1-y)`; stores `xs_isnrml=log10(1e-30)=-30, dx_isnrml=dy, y_isnrml(1:n)=x`. (Note loop recomputes n as `i-1`; grid marches `y=10**(log10(y)+dy)`.)

4. **Inherent Growth Parameterization splines** — `init_inherent_growth_par` (`mod_amps_utility.F90:1247-1321`), all constants hard-coded. Generating code verbatim (comments elided at `!!c` blocks only where they are pure history; the live statements complete):

```fortran
  subroutine init_inherent_growth_par
  use com_amps
  implicit none

  if ( IsMaster .and. debug ) then
    write(fid_alog,*) "Initializing inherent growth parameterization"
  end if
  ! +++ definition +++
   nok_igp=23

   x_igp = (/-60.0,-55.0, -50.0, -45.0, -40.0, -35.0, -30.0, -27.0, -25.0, -23.0, -21.0, -20.0, -17.0, &
             -15.0, -12.0, -10.0, -8.0, -6.0, -5.0, -4.0, -3.5, -2.5, -1.5 /)

    a_igp(1,:)=(/ 0.000000e+00 , 0.000000e+00 , -3.000000e-02 , 2.300000e+00 /)
    a_igp(2,:)=(/ 0.000000e+00 , 0.000000e+00 , -3.000000e-02 , 2.150000e+00 /)
    a_igp(3,:)=(/ -4.800000e-04 , 2.400000e-03 , -3.000000e-02 , 2.000000e+00 /)
    a_igp(4,:)=(/ 1.586667e-03 , -1.353333e-02 , -4.200000e-02 , 1.850000e+00 /)
    a_igp(5,:)=(/ 1.512821e-03 , -5.897436e-03 , -5.833333e-02 , 1.500000e+00 /)
    a_igp(6,:)=(/ -9.597381e-05 , 8.490998e-04 , -3.846154e-03 , 1.250000e+00 /)
    a_igp(7,:)=(/ -1.176598e-04 , 9.293226e-05 , -2.553191e-03 , 1.240000e+00 /)
    a_igp(8,:)=(/ 2.040230e-03 , -6.494253e-03 , -5.172414e-03 , 1.230000e+00 /)
    a_igp(9,:)=(/ -1.439394e-03 , 3.712121e-03 , -6.666667e-03 , 1.210000e+00 /)
    a_igp(10,:)=(/ -1.597052e-03 , -1.726044e-02 , -9.090909e-03 , 1.200000e+00 /)
    a_igp(11,:)=(/ 4.920791e-01 , -7.947818e-01 , -9.729730e-02 , 1.100000e+00 /)

    b_igp(1,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(2,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(3,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(4,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(5,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(6,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(7,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(8,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(9,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(10,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(11,:) = (/ 0.0431 , -0.1031 , 0.0 , 0.7600 /)

    a_igp(12,:) = (/-0.0012,0.0363,-0.2244,0.7000/)
    a_igp(13,:) = (/ -0.0005, 0.0141, -0.0511, 0.3200 /)
    a_igp(14,:) = (/  0.0035, 0.0109, -0.0011, 0.2700 /)
    a_igp(15,:) = (/ -0.0063, 0.0427, 0.1597, 0.4600 /)
    a_igp(16,:) = (/  0.0212, 0.0051, 0.2552, 0.9000 /)
    a_igp(17,:) = (/ -0.1109, 0.1320, 0.5294, 1.6000 /)
    a_igp(18,:) = (/  0.1061, -0.5332, -0.2729, 2.3000 /)
    a_igp(19,:) = (/  0.2758, -0.2148, -1.0209, 1.6000 /)
    a_igp(20,:) = (/ -0.0118, 0.6125, -0.6233, 0.6400 /)
    a_igp(21,:) = (/ -0.2951, 0.5948, -0.0197, 0.4800 /)
    a_igp(22,:) = (/  0.0550,   -0.0995, 0.1244, 0.7600 /)
    a_igp(23,:) = (/ -0.0000, 0.0108, 0.0906, 0.8400 /)

    b_igp(12:23,1:4)=a_igp(12:23,1:4)

  end subroutine init_inherent_growth_par
```

Targets: `x_igp(23), a_igp(23,4), b_igp(23,4)` in `common/VAPIGP` (`com_amps.F90:123-125`); copied into `CM%vigp` (vap_igp_aux). Knots are temperature (deg C); a/b rows are cubic coefficients per interval (a = gamma spline, b = second spline; identical for warm knots 12-23).

5. **Collisional-breakup fragment tables** (`bu_fd`, `bu_tmass`) — computed by `cal_breakfragment` (`mod_amps_lib.F90:1831-2017`), FULL verbatim:

```fortran
    subroutine cal_breakfragment(NRBIN,estbar,esitbar)
      use scale_prc, only: &
         PRC_abort
      use com_amps
      use class_AirGroup, only: &
         AirGroup, &
         make_AirGroup_2
      use class_Group, only: &
         Group, &
         make_group, &
         ini_group_MP
      use mod_amps_utility, only: &
         random_genvar
      use mod_amps_core, only: &
         cal_Coalescence_Efficiency, &
         cal_breakup_dis_LL, &
         diag_pq
      implicit none

  integer,intent(in) :: NRBIN

  ! saturation vapor lookup table
!tmp  real  :: estbar(150),esitbar(111)
  real(DS),intent(in)  :: estbar(150),esitbar(111)

  ! steady growth problem: phase over which saturation is calculated.
  ! phase 1: water, 2: ice
  integer    :: phase

  ! random generator vars (not used)
  type (random_genvar) :: rdsd

  ! steady growth problem: relative humidity
  real(PS)    :: RH

  ! define the vapor pressure
  real(PS)              :: e

  ! definition for reading command line
  character(len=20)      :: chbuf
  real(PS)               :: temp ! ambient temperature in Celcius


  real(PS),dimension(1) :: rv,den,pt,T,W,z


  type (Group)               :: liquid
  type (AirGroup)            :: steady


  ! coalescence efficiency
  real (PS)                 :: E_coal
  real (PS) :: D,D_L,D_S,S_T,S_C,DS_S,CKE
  ! low-diameter cut off related to the resolution of the experiments (cm)
  real(DS),parameter :: D_0=0.01

  integer                                :: i, j,n



  ! message from reality_check
!tmp  integer,pointer,dimension(:)   :: mes_rc
  integer,dimension(1)   :: mes_rc
  integer,dimension(1)   :: ID,JD,KD

  ! optional status check
  integer  :: var_Status,em

  ! coefficient: sqrt(2*pi), pi/6.0
  real(DS),parameter  :: coef1=2.506628275,coef2=0.523598776

  ! maximum number of elements for lookup table, bu_tmass and bu_fd
  integer :: i1d_pair_max,kk_max

  integer :: ierr

!tp  allocate(mes_rc(1))


  !if ( IsMaster ) then
  !  write(fid_alog,*) "Now calculating collisional breakup fragments."
  !end if

  bu_fd=0.0_PS
  bu_tmass=0.0_PS

  T=278.6795
  PT=850.0e+2_PS
  W=0.0
  phase=1
  RH=100.0_PS

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! construct liquid hydrometeor
!tmp  liquid = make_Group ( 1, 0.0, 1, NRBIN, &
!tmp       srat_r, sadd_r, minmass_r, 1,1.0,1.0,1.0,1.0,binbr)
  call make_Group (liquid,1, 0.0_PS, 1, NRBIN, &
       srat_r, sadd_r, minmass_r, 1,1.0_PS,den_aps(1),den_api(1),eps_ap(1),binbr)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! construt thermo_var object
  if ( IsMaster .and. debug ) then
     write(fid_alog,*) "ck estbar",estbar(1:10)
  end if

  steady=make_AirGroup_2(1,estbar,esitbar,RV,DEN,PT,T,W,phase,RH)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! initialize the mass and con of liquid
  call ini_group_MP(liquid)


  ! diagnose the property
  call diag_pq(liquid,steady,1,em, mes_rc,ID,JD,KD,rdsd,ihabit_gm_random &
              ,eps_ap(1),nu_aps,phi_aps,m_aps)


  ! find bin that has the minimum size for possible breakup.
  do i=1,NRBIN
     if(liquid%MS(i,1)%len>=D_0) then
        jmin_bk=i
        exit
     end if
  end do

  imin_bk=jmin_bk+1
  imax_bk=NRBIN
  jmax_bk=NRBIN-1

  i1d_pair_max=(imax_bk-1)-jmin_bk+1+(imax_bk-imin_bk)*(1+imax_bk-imin_bk)/2
  kk_max=i1d_pair_max*liquid%N_BIN

  if ( IsMaster .and. debug ) then
    write(fid_alog,*) "i1d_pair_max,kk_max",i1d_pair_max,kk_max

    write(fid_alog,*) "size of bu_tmass, bu_fd",size(bu_tmass),size(bu_fd)/2
  end if

  if(size(bu_tmass)<i1d_pair_max) then
    LOG_ERROR("cal_breakfragment",*) "Error.  Increase the bu_tmass array"
    call PRC_abort
  endif
  if(size(bu_fd)/2<kk_max) then
    LOG_ERROR("cal_breakfragment",*) "Error.  Increase the bu_fd array"
    call PRC_abort
  endif



  do i=imin_bk,imax_bk
     do j=jmin_bk,i-1
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! determine the coalescence efficiency
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        call cal_Coalescence_Efficiency(liquid,i,liquid,j,1,steady%TV(1),E_coal,&
             D_L,D_S,S_T,S_C,DS_S,CKE)

        if(CKE<=1.0e-20) cycle
!!c             call cal_breakup_dis(liquid,i,j,imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,bu_fd,&
!!c                                  1.0,D_L,D_S,S_T,S_C,dS_S,CKE)
        call cal_breakup_dis_LL(liquid,i,j,imin_bk,jmin_bk,bu_tmass,bu_fd,&
                                D_L,D_S,S_T,S_C,CKE)


     end do

     if(D_L*1.0e+2>0.4.and.D_L*1.0e+2<0.5) then
        call print_bufd(i,liquid,"check_bkup0p4.dat")
     elseif(D_L*1.0e+2>0.3.and.D_L*1.0e+2<0.4) then
        call print_bufd(i,liquid,"check_bkup0p3.dat")
     end if

  end do

!!c  write(fid_alog,*) "cal_break bu_fd"
!!c  write(fid_alog,*) bu_fd
!!c  write(fid_alog,*) "cal_break bu_tmass"
!!c  write(fid_alog,*) bu_tmass


  ! +++ deallocate memory used for the Group Object +++
!tmp  call delete_Group(liquid)
!tmp  call delete_AirGroup (steady)

!tmp  deallocate(mes_rc)
end subroutine cal_breakfragment
```

This is the Low–List (1982) breakup fragment machinery: for each colliding bin pair (i large, j small) above cutoff D_0 = 0.01 cm, `cal_Coalescence_Efficiency` gives collision kinetic energy CKE and surface energies S_T/S_C, then `cal_breakup_dis_LL` (mod_amps_core) fills `bu_tmass(pair)` (total fragment mass per pair) and `bu_fd(2, pair×bin)` (fragment number/mass distribution over destination bins). It depends on the runtime bin grid (`binbr`, `srat_r`, `sadd_r`, `minmass_r`) and on `estbar(150)`/`esitbar(111)` — which are the **saturation vapor pressure lookup tables** (water/ice) passed in from the host model, not AMPS_DATA files.

Also computed at runtime, not read: `estbar/esitbar` themselves (host-side sat-pressure tables, see `get_sat_vapor_pres_lk`, `mod_amps_lib.F90:1818-1829`), and `getznorm2` uses DCDFLIB `cdfnor` directly (`mod_amps_utility.F90:369-375`), so `stdnorm.dat`'s `znorm` is partially superseded but still read and stored.

## 6. Parser notes for the Python conversion plan

- All ASCII files are list-directed: whitespace-split works everywhere (some files space-, some tab-separated).
- Three header conventions: (a) collision-efficiency `*_Nre_Ec`/`drop_drop_Rey4`: 2 `#` comment lines + `nr nc xs dx ys dy` + nr rows × nc floats; (b) frequency/map files: 1 `#` comment line + `nrow ncol` + data (declared array may be larger than file ncol — mtac maps declared (50,**101**,2) but files have ncol=91; the tail stays uninitialized/zero, keep the file dims); (c) `stdnorm.dat`: headerless 451×6, keep cols 3-6 as (451,4), x-grid = 0.01*(row-1).
- `lmt_mass_*`: 4 columns, only col 2 retained (col 1 = temperature, cols 3-4 discarded).
- Fortran fills row-index-major here (row i = one file line), so a straight `np.loadtxt` gives the same (nr, nc) layout as `drpdrp(i,j)`; no transpose needed if you index `[row, col]`.
- Post-read transform to replicate: clip-to-zero and renormalize pol/pla/col frq by their sum (ros/ppo clipped only), RDCETB lines 226-242.
- `ap_act.bin.dat` can be skipped (RDAPTB returns immediately).
- Aux metadata per 2-D collision LUT should mirror `col_lut_aux`: `xs, dx, ys, dy, nr, nc` (x = row axis, y = column axis; e.g. drpdrp x = Reynolds bin start 0.004 step 0.005, y = log(ratio) start -3.9379 step 0.0286; hexdrp both axes log10(Nre) from -4.5/-4.3 step 0.1).
- Everything in section 5 needs no data files: osmotic and normal/inverse-normal LUTs and IGP splines are pure constants + scipy equivalents (`scipy.stats.norm.sf`, `norm.ppf`); breakup `bu_fd`/`bu_tmass` depend on the model bin grid and physics kernels and must be either ported (cal_Coalescence_Efficiency + cal_breakup_dis_LL in mod_amps_core.F90) or dumped once from Fortran.

Key file paths: `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_utility.F90` (RDCETB 96, RDAPTB 296, RDSTTB 333, init_inherent_growth_par 1247, init_osmo_par 12923, init_normal_lut 13087, init_inv_normal_lut 13131), `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_lib.F90` (readers called 626/632/638, bu allocation 443-447, cal_breakfragment 1831), `/Users/jcanton/projects/scale_amps/contrib/AMPS/com_amps.F90` (all common-block targets), `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Group.F90` (col_lut_aux 179, vap_igp_aux 187), `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Mass_Bin.F90` (data1d_lut 282), `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Cloud_Micro.F90` (members 171-224, copy 436-519).