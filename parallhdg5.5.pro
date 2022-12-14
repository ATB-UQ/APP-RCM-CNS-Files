remark   file protein-allhdg-ucl.param  version UCL  date 07-JUL-01
remark   for file protein-allhdg-ucl.top  version UCL  date 14-MAR-00
remark   for file protein-allhdg-dih-ucl.top  version UCL  date 07-JUL-01
remark   Geometric energy function parameters for distance geometry and
remark   simulated annealing.
remark   Original author: Michael Nilges, EMBL Heidelberg
remark   Modifications: Mark A. Williams, UCL London

set echo off message off end

!***********************************************************************!
! Copyright (C) 1995,1996 by Michael Nilges. All rights reserved.       !
! Copying and redistribution of this files is authorized only if etiher !
! (1) you make absolutely no changes to your copy, including name, or   !
! (2) if you do make changes, you name it something other than          !
! topallhdg.pro and topallhdg.x.xx.pro, and clearly mark the changes.   ! 
! The information in this software is subject to change without notice  !
! and should not be construed as a commitment by the EMBL, UCL, Yale    !
! University, UCL or by the authors. Neither the EMBL, Yale University, !
! UCL, nor the authors assume responsibility for the use or reliability !
! of this software. We do hope, however, to get responses from users,   !
! especially when errors have been found.                               !               !
!***********************************************************************!
! Description:                                                          !
! This parameter file was originally derived from the parameter files   !
! PARMALLH6 and PARHCSDX. It was designed specifically for the initial  !
! stages of calculating structures from NMR restraints.                 !
!***********************************************************************!
! PARHCSDX includes bond and angle parameters for non-hydrogen atoms    ! 
! derived from Cambridge Data Base model structures (R. A. Engh and R.  !
! Huber, Acta Cryst. Sect. A., 1991). Hydrogens were added with X-PLOR  !
! scripts for minimization and the PARAmeter LEARn statement. Dihedral, !
! improper and non-bonded values are from previous PARALLHDG versions,  !
! and assigned to new atom types where appropriate. Due to the          !
! minimization procedure used in the derivation, there are very small   !
! deviations from the parameter values in PARHCSDX.                     !
! Heavy atom types are exactly as in PARHCSDX, hydrogen types as in     !
! PARALLHDG.                                                            !
!***********************************************************************!
! History:
! version UCL   (07-Jul-01) : modified backbone & disulphide paprameters 
! version UCL   (05-Jul-01) : put dihedrals to maximum of 2 kcal/mol 
! version UCL   (03-Jul-01) : added new backbone dihedral parameters
! version UCL   (27-Jun-01) : restored original PROLSQ 1-4 parameters except larger aliphatic H radius
! version UCL   (02-Apr-00) : modified 1-4 PROLSQ parameters
! version UCL   (02-Apr-00) : modified aliphatic H radii
! version UCL   (14-Mar-00) : modified peptide bond and disulphide impropers/dihedrals
! version UCL   (12-Mar-00) : added or modified N and C terminus, CH3 & NH3 parameters
! version UCL   (10-Mar-00) : added CONTACT non-bond parameters
! version UCL   (09-Mar-00) : added or modified planar chi2 or chi3
! version 5.10  (24-Feb-99) : larger aliphatic hydrogens (test - not now active?)
! version 5.03  (13-Nov-98) : corrected aromatic chi2
! version 5.02  (19-Aug-98) : corrected cis peptide bond
! version 5.01  (20-Mar-98) : put side-chain dihedrals to 5 kcal/mol
! version 5.00  (16-Feb-98) : allow different nonbonded parameter options
! version 4.04  (28-Jan-98) : added missing dihedral parameters
! version 4.03  (27-Mar-97) : added missing dihedral parameters
! version 4.02  (25-Sep-96) : added missing parameters
! version 4.01  (29-Jul-96) : added missing covalent parameters
! version 4.00  (19-Jul-96) : all atom types from CSDX implemented
! version 3.00  (24-Oct-95) : mapped CSDX parmameters on parallhdg, 
!                             no changes in topallhdg
!
! previous modifications:
! proline residue modified (MN)
! added hbond acceptor and donor definitions for analysis (MN)
! all references to internal coordinates (IC's) removed (MN) 
! added stereospecific impropers for all pro-chiral centers (ATB, JK)
! modification of PARMALLH6 parameters to improve geometric consistency (JK)
! all dihedrals defining planarity converted to impropers (MN, PK)
! additional impropers at planar centers (MN)
!***********************************************************************!

   BOND NH1  H2   478.0   1.040
   BOND NH1  H3   478.0   1.040
   BOND NH1  H1   478.0   1.040
   BOND CH2E HB31 698.5   1.090
   BOND CH2E HB24 698.5   1.090
   BOND CH2E CG3  800.0   1.530
   BOND CG3  OD1  1000.8  1.250
   BOND CG3  ND2  901.2   1.340
   BOND ND2  HD211 1027.7  1.010
   BOND ND2  HD22 1027.7  1.010

   ANGLe CH1E NH1  H2   101.6   109.5
   ANGLe CH1E NH1  H3   101.6   109.5
   ANGLe CH1E NH1  H1   101.6   109.5
   ANGLe H2   NH1  H3   130.3   113.0
   ANGLe H2   NH1  H1   130.3   113.0
   ANGLe H3   NH1  H1   130.3   113.0
   ANGLe CH1E CH2E HB31 105.9   108.5
   ANGLe CH1E CH2E HB24 127.1   111.4
   ANGLe CH1E CH2E CG3  126.7   111.0
   ANGLe HB31 CH2E HB24 115.7   107.6
   ANGLe CG3  CH2E HB31 126.7   111.0
   ANGLe CG3  CH2E HB24 120.2   106.8
   ANGLe CH2E CG3  OD1  163.7   121.0
   ANGLe CH2E CG3  ND2  145.8   115.0
   ANGLe ND2  CG3  OD1  174.5   124.0
   ANGLe CG3  ND2  HD211 99.2    123.0
   ANGLe CG3  ND2  HD22 99.2    123.0
   ANGLe HD211 ND2  HD22 106.4   120.0

   DIHEdral H2   NH1  CH1E C    0.3     3  0.0
   DIHEdral CH1E CH2E CG3  ND2  0.2     6  0.0
   DIHEdral NH1  CH1E CH2E CG3  1.4     3  0.0
   DIHEdral CH2E CG3  ND2  HD211 8.0     2  180.0

   IMPRoper CG3  CH2E OD1  ND2  40.038  0  0.0
   IMPRoper ND2  CG3  HD211 HD22 40.038  0  0.0


! NBXMod=5 excludes all 1-2 and 1-3 pairs for non-bonded interactions,
! but the 1-4 nonbonded interactions are computed using the 1-4 Lennard-Jones
! parameters and the electrostatic scale factor E14Fac

   NBONds
   NBXMOD = 5
   E14Fac = 1
 END

   NONBONDED   CG3   1.025 2.81  1.025 2.53 !GROMOS LJ type: CPos
   NONBONDED   ND2   0.5969   4.07  0.5969   3.66 !GROMOS LJ type: NPri
   NONBONDED   OD1   0.4997   3.4   0.4997   3.06 !GROMOS LJ type: OEOp
   NONBONDED   H2 0  0  0  0 !GROMOS LJ type: HS14
   NONBONDED   H3 0  0  0  0 !GROMOS LJ type: HS14
   NONBONDED   H1 0  0  0  0 !GROMOS LJ type: HS14
   NONBONDED   HB31  0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HB24  0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HD211 0  0  0  0 !GROMOS LJ type: HS14
   NONBONDED   HD22  0  0  0  0 !GROMOS LJ type: HS14

!SEGMent
!   name=NTER
!   molecule number=1 name=NTER end
!END

  BOND CD3  OE2  1002.0  1.265
   BOND CD3  OE1  1002.0  1.265
   BOND CG2  CD3  358.5   1.560
   BOND CG2  HG32 699.8   1.100
   BOND CG2  HG22 698.5   1.090
   BOND CH2E CG2  800.0   1.530
   BOND CH2E HB3  698.5   1.090
   BOND CH2E HB2  698.5   1.090
   BOND C    O2   1000.8  1.250

   ANGLe OE1  CD3  OE2  184.0   126.0
   ANGLe CG2  CD3  OE2  151.8   117.0
   ANGLe CG2  CD3  OE1  151.8   117.0
   ANGLe CD3  CG2  HG32 120.2   106.8
   ANGLe CD3  CG2  HG22 107.1   109.5
   ANGLe CH2E CG2  CD3  126.7   111.0
   ANGLe HG32 CG2  HG22 120.2   106.8
   ANGLe CH2E CG2  HG32 105.9   108.5
   ANGLe CH2E CG2  HG22 126.7   111.0
   ANGLe CG2  CH2E HB3  125.2   110.3
   ANGLe CG2  CH2E HB2  121.2   107.6
   ANGLe CH1E CH2E CG2  126.7   111.0
   ANGLe HB3  CH2E HB2  120.2   106.8
   ANGLe CH1E CH2E HB3  120.2   106.8
   ANGLe CH1E CH2E HB2  107.1   109.5
   ANGLe CH1E C    O2   145.8   115.0
   ANGLe O    C    O2   184.0   126.0

   DIHEdral CH1E CH2E CG2  CD3  1.4     3  0.0
   DIHEdral CH2E CG2  CD3  OE1  0.2     6  0.0
   DIHEdral NH1  CH1E CH2E CG2  1.4     3  0.0

   IMPRoper C    CH1E O    O2   40.038  0  0.0
   IMPRoper CD3  OE2  OE1  CG2  40.038  0  0.0


! NBXMod=5 excludes all 1-2 and 1-3 pairs for non-bonded interactions,
! but the 1-4 nonbonded interactions are computed using the 1-4 Lennard-Jones
! parameters and the electrostatic scale factor E14Fac

   NBONds
   NBXMOD = 5
   E14Fac = 1
 END

   NONBONDED   CG2   0.2774   3.58  0.2774   3.22 !GROMOS LJ type: C
   NONBONDED   CD3   1.025 2.81  1.025 2.53 !GROMOS LJ type: CPos
   NONBONDED   OE1   1.725 2.63  1.725 2.36 !GROMOS LJ type: OM
   NONBONDED   OE2   1.725 2.63  1.725 2.36 !GROMOS LJ type: OM
   NONBONDED   O2 1.725 2.63  1.725 2.36 !GROMOS LJ type: OM
   NONBONDED   HB3   0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HB2   0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HG32  0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HG22  0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC

!SEGMent
!   name=CTER
!   molecule number=1 name=CTER end
!END
  BOND CH1E CB2  478.0   1.540
   BOND CH1E CB1  358.5   1.560
   BOND CB2  HB21 698.5   1.090
   BOND CB2  HB22 698.5   1.090
   BOND CB2  HB23 698.5   1.090
   BOND CB1  HB13 698.5   1.090
   BOND CB1  HB12 698.5   1.090
   BOND CB1  CG   800.0   1.530
   BOND CG   HG2  698.5   1.090
   BOND CG   HG3  698.5   1.090
   BOND CG   CD1  478.0   1.540
   BOND CD1  HD2  699.8   1.100
   BOND CD1  HD3  698.5   1.090
   BOND CD1  CE   406.3   1.510
   BOND CE   HE   698.5   1.090
   BOND CE   CZ   997.7   1.330
   BOND CE   CE   997.7   1.330
   BOND CZ   HZ   698.5   1.090
   BOND CZ   CH   900.2   1.500
   BOND CH   HH1  698.5   1.090
   BOND CH   HH2  698.5   1.090
   BOND CH   HH3  698.5   1.090

   ANGLe CB2  CH1E C    126.7   111.0
   ANGLe CB1  CH1E C    124.3   109.5
   ANGLe CB1  CH1E CB2  126.7   111.0
   ANGLe NH1  CH1E CB2  124.3   109.5
   ANGLe NH1  CH1E CB1  126.7   111.0
   ANGLe CH1E CB2  HB21 401.6   109.0
   ANGLe CH1E CB2  HB22 401.6   109.0
   ANGLe CH1E CB2  HB23 401.6   109.0
   ANGLe HB21 CB2  HB22 105.9   108.5
   ANGLe HB21 CB2  HB23 105.9   108.5
   ANGLe HB22 CB2  HB23 105.9   108.5
   ANGLe CH1E CB1  HB13 651.6   107.0
   ANGLe CH1E CB1  HB12 111.1   108.0
   ANGLe CH1E CB1  CG   126.7   111.0
   ANGLe HB13 CB1  HB12 120.2   106.8
   ANGLe CG   CB1  HB13 107.5   109.6
   ANGLe CG   CB1  HB12 107.5   109.6
   ANGLe CB1  CG   HG2  107.5   109.6
   ANGLe CB1  CG   HG3  125.2   110.3
   ANGLe CB1  CG   CD1  126.7   111.0
   ANGLe HG2  CG   HG3  414.3   106.0
   ANGLe CD1  CG   HG2  105.9   108.5
   ANGLe CD1  CG   HG3  401.6   109.0
   ANGLe CG   CD1  HD2  105.9   108.5
   ANGLe CG   CD1  HD3  68.1    109.5
   ANGLe CG   CD1  CE   126.7   111.0
   ANGLe HD2  CD1  HD3  120.2   106.8
   ANGLe CE   CD1  HD2  107.5   109.6
   ANGLe CE   CD1  HD3  68.1    110.0
   ANGLe CD1  CE   HE   120.7   120.0
   ANGLe CD1  CE   CZ   153.0   126.0
   ANGLe CD1  CE   CE   153.0   126.0
   ANGLe CZ   CE   HE   120.7   120.0
   ANGLe CE   CZ   HZ   120.7   120.0
   ANGLe CE   CZ   CH   153.0   126.0
   ANGLe CH   CZ   HZ   120.7   120.0
   ANGLe CZ   CH   HH1  151.0   111.3
   ANGLe CZ   CH   HH2  151.0   111.3
   ANGLe CZ   CH   HH3  151.0   111.3
   ANGLe HH1  CH   HH2  111.1   108.0
   ANGLe HH1  CH   HH3  111.1   108.0
   ANGLe HH2  CH   HH3  111.1   108.0

   DIHEdral CH1E CB1  CG   CD1  1.4     3  0.0
   DIHEdral CB1  CG   CD1  CE   1.4     3  0.0
   DIHEdral CG   CD1  CE   CZ   0.2     6  180.0
   DIHEdral CD1  CE   CZ   CH   10.0    2  180.0
   DIHEdral CE   CZ   CH   HH1  0.2     6  180.0
   DIHEdral NH1  CH1E CB2  HB21 1.4     3  0.0
   DIHEdral NH1  CH1E CB1  CG   1.4     3  0.0

   IMPRoper CE   CD1  HE   CZ   40.038  0  0.0
   IMPRoper CZ   CE   HZ   CH   40.038  0  0.0


! NBXMod=5 excludes all 1-2 and 1-3 pairs for non-bonded interactions,
! but the 1-4 nonbonded interactions are computed using the 1-4 Lennard-Jones
! parameters and the electrostatic scale factor E14Fac

   NBONds
   NBXMOD = 5
   E14Fac = 1
 END

   NONBONDED   CB1   0.2774   3.58  0.2774   3.22 !GROMOS LJ type: C
   NONBONDED   CG 0.2774   3.58  0.2774   3.22 !GROMOS LJ type: C
   NONBONDED   CD1   1.025 2.81  1.025 2.53 !GROMOS LJ type: CPos
   NONBONDED   CE 0.2774   3.58  0.2774   3.22 !GROMOS LJ type: C
   NONBONDED   CZ 0.2774   3.58  0.2774   3.22 !GROMOS LJ type: C
   NONBONDED   CH 0.2774   3.58  0.2774   3.22 !GROMOS LJ type: C
   NONBONDED   CB2   0.2774   3.58  0.2774   3.22 !GROMOS LJ type: C
   NONBONDED   HB13  0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HB12  0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HG2   0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HD2   0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HG3   0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HE 0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HD3   0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HZ 0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HH1   0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HH2   0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HH3   0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HB21  0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HB22  0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC
   NONBONDED   HB23  0.1184   2.37  0.1184   2.14 !GROMOS LJ type: HC

!SEGMent
!   name=PAL
!   molecule number=1 name=ST1 end
!END
 


!added by nocky
 BOND  S     C      1000.000 {sd=     0.001}      1.822
 ANGLe  CH2E S     C        500.00 {sd=     0.031}    100.00   
 ANGLe  S    C     CH1E     500.00 {sd=     0.031}    113.00
 ANGLe  S    C     O        500.00 {sd=     0.031}    120.00
 IMPRoper  HA   CH1E CH2E HA        500.00 {sd=     0.031}    0    180.0000 
 DIHEdral  O       C     CH1E    HA    5.00     3  0.0000
 DIHEdral  CH2E    S     C      CH1E   5.00     3  0.0000
 DIHEdral  CH2E    S     C      O      5.00     3  0.0000




 BOND  C    CH1E    1000.000 {sd=     0.001}      1.525
 BOND  C    CH2E    1000.000 {sd=     0.001}      1.516
 BOND  C    CH3E    1000.000 {sd=     0.001}      1.507   ! added EAB
 BOND  C    CH2G    1000.000 {sd=     0.001}      1.516
 BOND  C    N       1000.000 {sd=     0.001}      1.341
 BOND  C    NC2     1000.000 {sd=     0.001}      1.326
 BOND  C    NH1     1000.000 {sd=     0.001}      1.329
 BOND  C    NH2     1000.000 {sd=     0.001}      1.328
 BOND  C    O       1000.000 {sd=     0.001}      1.231
 BOND  C    OC      1000.000 {sd=     0.001}      1.249
 BOND  C5   CH2E    1000.000 {sd=     0.001}      1.497
 BOND  C5   CR1E    1000.000 {sd=     0.001}      1.357
 BOND  C5   CR1H    1000.000 {sd=     0.001}      1.354
 BOND  C5   NH1     1000.000 {sd=     0.001}      1.378
 BOND  C5   NR      1000.000 {sd=     0.001}      1.373
 BOND  C5W  CH2E    1000.000 {sd=     0.001}      1.498
 BOND  C5W  CR1E    1000.000 {sd=     0.001}      1.365
 BOND  C5W  CW      1000.000 {sd=     0.001}      1.433
 BOND  CF   CH2E    1000.000 {sd=     0.001}      1.502
 BOND  CF   CR1E    1000.000 {sd=     0.001}      1.385
 BOND  CH1E CH1E    1000.000 {sd=     0.001}      1.540
 BOND  CH1E CH2E    1000.000 {sd=     0.001}      1.530
 BOND  CH1E CH3E    1000.000 {sd=     0.001}      1.521
 BOND  CH1E HA      1000.000 {sd=     0.001}      1.080
 BOND  CH1E N       1000.000 {sd=     0.001}      1.466
 BOND  CH1E NH1     1000.000 {sd=     0.001}      1.458
 BOND  CH1E NH3     1000.000 {sd=     0.001}      1.486
 BOND  CH1E OH1     1000.000 {sd=     0.001}      1.433 !included 14-apr-99
 BOND  CH2E CH2E    1000.000 {sd=     0.001}      1.520
 BOND  CH2E CH2P    1000.000 {sd=     0.001}      1.490
 BOND  CH2E CH3E    1000.000 {sd=     0.001}      1.513
 BOND  CH2E CY      1000.000 {sd=     0.001}      1.512
 BOND  CH2E HA      1000.000 {sd=     0.001}      1.080
 BOND  CH2E NH1     1000.000 {sd=     0.001}      1.460
 BOND  CH2E NH3     1000.000 {sd=     0.001}      1.489
 BOND  CH2E OH1     1000.000 {sd=     0.001}      1.417
 BOND  CH2E S       1000.000 {sd=     0.001}      1.822
 BOND  CH2E SH1E    1000.000 {sd=     0.001}      1.808
 BOND  CH2E SM      1000.000 {sd=     0.001}      1.803
 BOND  CH2G HA      1000.000 {sd=     0.001}      1.080
 BOND  CH2G NH1     1000.000 {sd=     0.001}      1.451
 BOND  CH2G NH3     1000.000 {sd=     0.001}      1.489
 BOND  CH2P CH2P    1000.000 {sd=     0.001}      1.504
 BOND  CH2P HA      1000.000 {sd=     0.001}      1.080
 BOND  CH2P N       1000.000 {sd=     0.001}      1.473
 BOND  CH2P NH3     1000.000 {sd=     0.001}      1.474
 BOND  CH3E HA      1000.000 {sd=     0.001}      1.080
 BOND  CH3E SM      1000.000 {sd=     0.001}      1.791
 BOND  CR1E CR1E    1000.000 {sd=     0.001}      1.382
 BOND  CR1E CR1W    1000.000 {sd=     0.001}      1.400
 BOND  CR1E CW      1000.000 {sd=     0.001}      1.398
 BOND  CR1E CY      1000.000 {sd=     0.001}      1.389
 BOND  CR1E CY2     1000.000 {sd=     0.001}      1.379
 BOND  CR1E HA      1000.000 {sd=     0.001}      1.080
 BOND  CR1E NH1     1000.000 {sd=     0.001}      1.373
 BOND  CR1E NR      1000.000 {sd=     0.001}      1.391
 BOND  CR1H HA      1000.000 {sd=     0.001}      1.081
 BOND  CR1H NH1     1000.000 {sd=     0.001}      1.374
 BOND  CR1H NR      1000.000 {sd=     0.001}      1.374
 BOND  CR1W CR1W    1000.000 {sd=     0.001}      1.368
 BOND  CR1W CW      1000.000 {sd=     0.001}      1.393
 BOND  CR1W HA      1000.000 {sd=     0.001}      1.080
 BOND  CRH  HA      1000.000 {sd=     0.001}      1.080
 BOND  CRH  NH1     1000.000 {sd=     0.001}      1.341
 BOND  CRH  NR      1000.000 {sd=     0.001}      1.319
 BOND  CRHH HA      1000.000 {sd=     0.001}      1.080
 BOND  CRHH NH1     1000.000 {sd=     0.001}      1.321
 BOND  CRHH NR     1000.000 {sd=     0.001}      1.321
 BOND  CW   CW      1000.000 {sd=     0.001}      1.409
 BOND  CW   NH1     1000.000 {sd=     0.001}      1.370
 BOND  CY2  OH1     1000.000 {sd=     0.001}      1.376
 BOND  CY2  OHP     1000.000 {sd=     0.001}      1.376
 BOND  H    NH1     1000.000 {sd=     0.001}      0.980
 BOND  H    N       1000.000 {sd=     0.001}      0.980
 BOND  H    NH2     1000.000 {sd=     0.001}      0.980
 BOND  H    OH1     1000.000 {sd=     0.001}      0.960
 BOND  H    SH1E    1000.000 {sd=     0.001}      1.325 ! from CHARMM22
 BOND  HC   NC2     1000.000 {sd=     0.001}      1.000
 BOND  HC   NH3     1000.000 {sd=     0.001}      1.040
 BOND  S    S       1000.000 {sd=     0.001}      2.030
 BOND  CCIS  CH1E    1000.000 {sd=     0.001}      1.525 ! CCIS is for cis peptide bond
 BOND  CCIS  CH2E    1000.000 {sd=     0.001}      1.516
 BOND  CCIS  CH2G    1000.000 {sd=     0.001}      1.516
 BOND  CCIS  N       1000.000 {sd=     0.001}      1.341
 BOND  CCIS  NC2     1000.000 {sd=     0.001}      1.326
 BOND  CCIS  NH1     1000.000 {sd=     0.001}      1.329
 BOND  CCIS  NH2     1000.000 {sd=     0.001}      1.328
 BOND  CCIS  O       1000.000 {sd=     0.001}      1.231
 BOND  CCIS  OC      1000.000 {sd=     0.001}      1.249
 BOND  NR    P       1000.0  1.610
 BOND  SM    P       1000.0  2.094
 BOND  P     O3R     1000.0  1.610
 BOND  P     OHP     1000.0  1.610
 BOND  P     O1P     1000.0  1.485
 BOND  OC    CA+2    1000.0  2.500
 BOND  CH2E CH1P    1000.000 {sd=     0.001}      1.490 ! for HYP (KJR)
 BOND  CH1P CH2P    1000.000 {sd=     0.001}      1.504 ! for HYP (KJR)
 BOND  CH1P HA      1000.000 {sd=     0.001}      1.080 ! for HYP (KJR)
 BOND  CH2P OH      1000.000 {sd=     0.001}      1.417 ! for HYP (KJR)
 BOND  CH1P OH      1000.000 {sd=     0.001}      1.417 ! for HYP (KJR)
 BOND  H    OH      1000.000 {sd=     0.001}      0.960 ! for HYP (KJR)
 BOND  CH2E CQ      1000.000 {sd=     0.001}      1.507 ! for PRG when in triazole bridge (KJR)
 BOND  CQ   CL      1000.000 {sd=     0.001}      1.350 ! for PRG when in triazole birdge (KJR)
 BOND  CL   HA      1000.000 {sd=     0.001}      1.080 ! for PRG when in triazole bridge (KJR)
 BOND  CH2E NK      1000.000 {sd=     0.001}      1.400 ! for AZA/AZH when in triazole bridge (KJR)
 BOND  NK   NL      1000.000 {sd=     0.001}      1.403 ! for AZA/AZH when in triazole bridge (KJR)
 BOND  NL   NM      1000.000 {sd=     0.001}      1.288 ! for AZA/AZH when in triazole bridge (KJR)
 BOND  NM   CQ      1000.000 {sd=     0.001}      1.334 ! for triazole bridge (KJR)
 BOND  CL   NK      1000.000 {sd=     0.001}      1.334 ! for triazole bridge (KJR)
 BOND  CH1E CH2K    1000.000 {sd=     0.001}      1.526 !LJ 31/01/07
 BOND  CH2K CH2L    1000.000 {sd=     0.001}      1.526 !LJ 31/01/07
 BOND  CH2K HA      1000.000 {sd=     0.001}      1.090 !LJ 31/01/07
 BOND  CH2L CJ      1000.000 {sd=     0.001}      1.522 !LJ 31/01/07
 BOND  CH2L HA      1000.000 {sd=     0.001}      1.090 !LJ 31/01/07
 BOND  CJ   N       1000.000 {sd=     0.001}      1.335 !LJ 31/01/07
 BOND  CJ   O       1000.000 {sd=     0.001}      1.229 !LJ 31/01/07
 
 ANGLe  C    CH1E CH1E     500.00 {sd=     0.031}    109.0754
 ANGLe  C    CH1E CH2E     500.00 {sd=     0.031}    110.1094
 ANGLe  C    CH1E CH3E     500.00 {sd=     0.031}    110.4838
 ANGLe  C    CH1E HA       500.00 {sd=     0.031}    108.9914
 ANGLe  C    CH1E N        500.00 {sd=     0.031}    111.9082
 ANGLe  C    CH1E NH1      500.00 {sd=     0.031}    111.1396
 ANGLe  C    CH1E NH3      500.00 {sd=     0.031}    111.1703
 ANGLe  C    CH2E CH1E     500.00 {sd=     0.031}    112.5947
 ANGLe  C    CH2E CH2E     500.00 {sd=     0.031}    112.5943
 ANGLe  C    CH2E HA       500.00 {sd=     0.031}    108.5877
 ANGLe  C    CH3E HA       500.00 {sd=     0.031}    108.5877
 ANGLe  C    CH2G HA       500.00 {sd=     0.031}    108.8528
 ANGLe  C    CH2G NH1      500.00 {sd=     0.031}    112.4999
 ANGLe  C    CH2G NH3      500.00 {sd=     0.031}    112.4990
 ANGLe  C    N    CH1E     500.00 {sd=     0.031}    122.7632
 ANGLe  C    N    CH2P     500.00 {sd=     0.031}    125.1134
 ANGLe  C    NC2  HC       500.00 {sd=     0.031}    119.9992
 ANGLe  C    NH1  CH1E     500.00 {sd=     0.031}    121.6541
 ANGLe  C    NH1  CH2E     500.00 {sd=     0.031}    124.1226
 ANGLe  C    NH1  CH2G     500.00 {sd=     0.031}    120.5859
 ANGLe  C    NH1  H        500.00 {sd=     0.031}    119.2489
 ANGLe  C    NH2  H        500.00 {sd=     0.031}    118.1853
 ANGLe  C5   CH2E CH1E     500.00 {sd=     0.031}    113.7931
 ANGLe  C5   CH2E HA       500.00 {sd=     0.031}    108.1195
 ANGLe  C5   CR1E HA       500.00 {sd=     0.031}    126.2616
 ANGLe  C5   CR1E NH1      500.00 {sd=     0.031}    106.5015
 ANGLe  C5   CR1E NR       500.00 {sd=     0.031}    109.4272
 ANGLe  C5   CR1H HA       500.00 {sd=     0.031}    126.4031
 ANGLe  C5   CR1H NH1      500.00 {sd=     0.031}    107.1610
 ANGLe  C5   CR1H NR       500.00 {sd=     0.031}    107.1610
 ANGLe  C5   NH1  CRH      500.00 {sd=     0.031}    108.0959
 ANGLe  C5   NH1  CRHH     500.00 {sd=     0.031}    109.4352
 ANGLe  C5   NR   CRHH     500.00 {sd=     0.031}    109.4352
 ANGLe  C5   NH1  H        500.00 {sd=     0.031}    126.0497
 ANGLe  C5   NR   CRH      500.00 {sd=     0.031}    105.7163
 ANGLe  C5W  CH2E CH1E     500.00 {sd=     0.031}    113.5788
 ANGLe  C5W  CH2E HA       500.00 {sd=     0.031}    108.1831
 ANGLe  C5W  CR1E HA       500.00 {sd=     0.031}    124.5037
 ANGLe  C5W  CR1E NH1      500.00 {sd=     0.031}    110.0962
 ANGLe  C5W  CW   CR1E     500.00 {sd=     0.031}    133.9320
 ANGLe  C5W  CW   CW       500.00 {sd=     0.031}    107.2333
 ANGLe  CF   CH2E CH1E     500.00 {sd=     0.031}    113.7937
 ANGLe  CF   CH2E HA       500.00 {sd=     0.031}    108.1268
 ANGLe  CF   CR1E CR1E     500.00 {sd=     0.031}    120.7850
 ANGLe  CF   CR1E HA       500.00 {sd=     0.031}    119.4540
 ANGLe  CH1E C    N        500.00 {sd=     0.031}    116.9940
 ANGLe  CH3E C    N        500.00 {sd=     0.031}    116.9940
 ANGLe  CH1E C    NH1      500.00 {sd=     0.031}    116.1998
 ANGLe  CH3E C    NH1      500.00 {sd=     0.031}    116.1998
 ANGLe  CH1E C    O        500.00 {sd=     0.031}    120.8258
 ANGLe  CH1E C    OC       500.00 {sd=     0.031}    118.0611
 ANGLe  CH1E CH1E CH2E     500.00 {sd=     0.031}    110.3824
 ANGLe  CH1E CH1E CH3E     500.00 {sd=     0.031}    110.4882
 ANGLe  CH1E CH1E HA       500.00 {sd=     0.031}    108.2775
 ANGLe  CH1E CH1E NH1      500.00 {sd=     0.031}    111.4875
 ANGLe  CH1E CH1E NH3      500.00 {sd=     0.031}    111.4875 
 ANGLe  CH1E CH1E OH1      500.00 {sd=     0.031}    109.600  ! included 14-APR-99
 ANGLe  CH1E CH2E CH1E     500.00 {sd=     0.031}    116.0385
 ANGLe  CH1E CH2E CH2E     500.00 {sd=     0.031}    114.0589
 ANGLe  CH1E CH2E CH2P     500.00 {sd=     0.031}    104.3952
 ANGLe  CH1E CH2E CH3E     500.00 {sd=     0.031}    113.7404
 ANGLe  CH1E CH2E CY       500.00 {sd=     0.031}    113.8748
 ANGLe  CH1E CH2E HA       500.00 {sd=     0.031}    109.2833
 ANGLe  CH1E CH2E OH1      500.00 {sd=     0.031}    111.0593
 ANGLe  CH1E CH2E S        500.00 {sd=     0.031}    114.3551
 ANGLe  CH1E CH2E SH1E     500.00 {sd=     0.031}    114.3558
 ANGLe  CH1E CH2E SM       500.00 {sd=     0.031}    114.3558
 ANGLe  CH1E CH3E HA       500.00 {sd=     0.031}    109.4726
 ANGLe  CH1E N    CH2P     500.00 {sd=     0.031}    112.1234
 ANGLe  CH1E NH1  H        500.00 {sd=     0.031}    119.2367
 ANGLe  CH1E NH3  CH2P     500.00 {sd=     0.031}    110.6738
!!! ANGLe  CH1E NH3  HC       500.00 {sd=     0.031}    104.9708
 ANGLe  CH1E NH3  HC       500.00 {sd=     0.031}    109.4688 !estimated (from Gly) 14-APR-99
 ANGLe  CH1E OH1  H        500.00 {sd=     0.031}    109.5    ! included 14-APR-99
 ANGLe  CH2E C    NH2      500.00 {sd=     0.031}    116.4617
 ANGLe  CH2E C    O        500.00 {sd=     0.031}    120.9106
 ANGLe  CH3E C    O        500.00 {sd=     0.031}    120.9106 ! added EAB
 ANGLe  CH2E C    OC       500.00 {sd=     0.031}    118.4969
 ANGLe  CH2E C5   CR1E     500.00 {sd=     0.031}    129.6325
 ANGLe  CH2E C5   CR1H     500.00 {sd=     0.031}    131.2043
 ANGLe  CH2E C5   NH1      500.00 {sd=     0.031}    123.4237
 ANGLe  CH2E C5   NR       500.00 {sd=     0.031}    121.5772
 ANGLe  CH2E C5W  CR1E     500.00 {sd=     0.031}    126.9191
 ANGLe  CH2E C5W  CW       500.00 {sd=     0.031}    126.8167
 ANGLe  CH2E CF   CR1E     500.00 {sd=     0.031}    120.6527
 ANGLe  CH2E CH1E CH3E     500.00 {sd=     0.031}    110.6376
 ANGLe  CH2E CH1E HA       500.00 {sd=     0.031}    109.2487
 ANGLe  CH2E CH1E N        500.00 {sd=     0.031}    103.0552
 ANGLe  CH2E CH1E NH1      500.00 {sd=     0.031}    110.4763
 ANGLe  CH2E CH1E NH3      500.00 {sd=     0.031}    110.4763  ! N term from 108.4924 14-MAR-00
 ANGLe  CH2E CH2E CH2E     500.00 {sd=     0.031}    111.3121
 ANGLe  CH2E CH2E HA       500.00 {sd=     0.031}    108.7236
 ANGLe  CH2E CH2E NH1      500.00 {sd=     0.031}    111.9978
 ANGLe  CH2E CH2E NH3      500.00 {sd=     0.031}    111.8939
 ANGLe  CH2E CH2E SM       500.00 {sd=     0.031}    112.6822
 ANGLe  CH2E CH2P CH2P     500.00 {sd=     0.031}    106.1000
 ANGLe  CH2E CH2P HA       500.00 {sd=     0.031}    109.9548
 ANGLe  CH2E CH3E HA       500.00 {sd=     0.031}    109.4694
 ANGLe  CH2E CY   CR1E     500.00 {sd=     0.031}    120.9304
 ANGLe  CH2E NH1  H        500.00 {sd=     0.031}    118.0987
 ANGLe  CH2E NH3  HC       500.00 {sd=     0.031}    109.4693
 ANGLe  CH2E OH1  H        500.00 {sd=     0.031}    109.4969
 ANGLe  CH2E S    S        500.00 {sd=     0.031}    103.7998
 ANGLe  CH2E SH1E H        500.00 {sd=     0.031}    107.9769
 ANGLe  CH2E SM   CH3E     500.00 {sd=     0.031}    100.8987
 ANGLe  CH2E SM   P        500.00 {sd=     0.031}     96.53
 ANGLe  CH2G C    N        500.00 {sd=     0.031}    117.7918
 ANGLe  CH2G C    NH1      500.00 {sd=     0.031}    116.3225
 ANGLe  CH2G C    O        500.00 {sd=     0.031}    120.6203
 ANGLe  CH2G C    OC       500.00 {sd=     0.031}    118.4971
 ANGLe  CH2G NH1  H        500.00 {sd=     0.031}    119.7297
 ANGLe  CH2G NH3  HC       500.00 {sd=     0.031}    109.4688
 ANGLe  CH2P CH2E HA       500.00 {sd=     0.031}    111.1127
 ANGLe  CH2P CH2P HA       500.00 {sd=     0.031}    110.3818
 ANGLe  CH2P CH2P N        500.00 {sd=     0.031}    103.2695
 ANGLe  CH2P CH2P NH3      500.00 {sd=     0.031}    103.6880
 ANGLe  CH2P NH3  HC       500.00 {sd=     0.031}    123.8148
 ANGLe  CH3E CH1E CH3E     500.00 {sd=     0.031}    110.7707
 ANGLe  CH3E CH1E HA       500.00 {sd=     0.031}    108.1279
 ANGLe  CH3E CH1E NH1      500.00 {sd=     0.031}    110.3844
 ANGLe  CH3E CH1E NH3      500.00 {sd=     0.031}    110.4751
 ANGLe  CH3E CH1E OH1      500.00 {sd=     0.031}    109.300 !included 14-APR-99
 ANGLe  CH3E CH2E HA       500.00 {sd=     0.031}    108.0408
 ANGLe  CH3E CH2E OH1      500.00 {sd=     0.031}    108.0961
 ANGLe  CR1E C5   NH1      500.00 {sd=     0.031}    105.6758
 ANGLe  CR1E C5   NR       500.00 {sd=     0.031}    109.3402
 ANGLe  CR1E C5W  CW       500.00 {sd=     0.031}    106.2641
 ANGLe  CR1E CF   CR1E     500.00 {sd=     0.031}    118.6946
 ANGLe  CR1E CR1E CR1E     500.00 {sd=     0.031}    119.9118
 ANGLe  CR1E CR1E CR1W     500.00 {sd=     0.031}    121.1513
 ANGLe  CR1E CR1E CW       500.00 {sd=     0.031}    118.6734
 ANGLe  CR1E CR1E CY       500.00 {sd=     0.031}    121.1348
 ANGLe  CR1E CR1E CY2      500.00 {sd=     0.031}    119.6224
 ANGLe  CR1E CR1E HA       500.00 {sd=     0.031}    119.9433
 ANGLe  CR1E CR1W CR1W     500.00 {sd=     0.031}    121.4832
 ANGLe  CR1E CR1W HA       500.00 {sd=     0.031}    118.7598
 ANGLe  CR1E CW   CW       500.00 {sd=     0.031}    118.8347
 ANGLe  CR1E CY   CR1E     500.00 {sd=     0.031}    118.1392
 ANGLe  CR1E CY2  CR1E     500.00 {sd=     0.031}    120.3463
 ANGLe  CR1E CY2  OH1      500.00 {sd=     0.031}    119.8269
 ANGLe  CR1E CY2  OHP      500.00 {sd=     0.031}    119.8269
 ANGLe  CR1E NH1  CRH      500.00 {sd=     0.031}    106.8630
 ANGLe  CR1E NH1  CW       500.00 {sd=     0.031}    108.9983
 ANGLe  CR1E NH1  H        500.00 {sd=     0.031}    125.8235
 ANGLe  CR1E NR   CRH      500.00 {sd=     0.031}    105.7678
 ANGLe  CR1H C5   NH1      500.00 {sd=     0.031}    106.0900
 ANGLe  CR1H NH1  CRHH     500.00 {sd=     0.031}    108.9901
 ANGLe  CR1H NR   CRHH     500.00 {sd=     0.031}    108.9901
 ANGLe  CR1H NH1  H        500.00 {sd=     0.031}    125.5054
 ANGLe  CR1W CR1E HA       500.00 {sd=     0.031}    119.1706
 ANGLe  CR1W CR1W CW       500.00 {sd=     0.031}    117.4515
 ANGLe  CR1W CR1W HA       500.00 {sd=     0.031}    120.2616
 ANGLe  CR1W CW   CW       500.00 {sd=     0.031}    122.4059
 ANGLe  CR1W CW   NH1      500.00 {sd=     0.031}    130.1860
 ANGLe  CRH  NH1  H        500.00 {sd=     0.031}    126.0322
 ANGLe  CRHH NH1  H        500.00 {sd=     0.031}    125.1896
 ANGLe  CRHH NR   H        500.00 {sd=     0.031}    125.1896
 ANGLe  CW   CR1E HA       500.00 {sd=     0.031}    121.0317
 ANGLe  CW   CR1W HA       500.00 {sd=     0.031}    121.7822
 ANGLe  CW   CW   NH1      500.00 {sd=     0.031}    107.4081
 ANGLe  CW   NH1  H        500.00 {sd=     0.031}    125.9221
 ANGLe  CY   CH2E HA       500.00 {sd=     0.031}    108.0910
 ANGLe  CY   CR1E HA       500.00 {sd=     0.031}    119.1931
 ANGLe  CY2  CR1E HA       500.00 {sd=     0.031}    120.3261
 ANGLe  CY2  OH1  H        500.00 {sd=     0.031}    109.4984
 ANGLe  CY2  OHP  P        500.00 {sd=     0.031}    120.0000
 ANGLe  H    NH2  H        500.00 {sd=     0.031}    123.6294
 ANGLe  HA   CH1E N        500.00 {sd=     0.031}    111.0977
 ANGLe  HA   CH1E NH1      500.00 {sd=     0.031}    108.0508
 ANGLe  HA   CH1E NH3      500.00 {sd=     0.031}    108.5074
 ANGLe  HA   CH1E OH1      500.00 {sd=     0.031}    108.6930 !assume like CH2E
 ANGLe  HA   CH2E HA       500.00 {sd=     0.031}    109.4074
 ANGLe  HA   CH2E NH1      500.00 {sd=     0.031}    108.9030
 ANGLe  HA   CH2E NH3      500.00 {sd=     0.031}    108.9390
 ANGLe  HA   CH2E OH1      500.00 {sd=     0.031}    108.6930
 ANGLe  HA   CH2E S        500.00 {sd=     0.031}    107.9228
 ANGLe  HA   CH2E SH1E     500.00 {sd=     0.031}    107.9185
 ANGLe  HA   CH2E SM       500.00 {sd=     0.031}    108.6768
 ANGLe  HA   CH2G HA       500.00 {sd=     0.031}    108.8718
 ANGLe  HA   CH2G NH1      500.00 {sd=     0.031}    108.8510
 ANGLe  HA   CH2G NH3      500.00 {sd=     0.031}    108.8586
 ANGLe  HA   CH2P HA       500.00 {sd=     0.031}    110.4563
 ANGLe  HA   CH2P N        500.00 {sd=     0.031}    110.8278
 ANGLe  HA   CH2P NH3      500.00 {sd=     0.031}    110.7246
 ANGLe  HA   CH3E HA       500.00 {sd=     0.031}    109.4703
 ANGLe  HA   CH3E SM       500.00 {sd=     0.031}    109.4700
 ANGLe  HA   CR1E NH1      500.00 {sd=     0.031}    125.8803
 ANGLe  HA   CR1E NR       500.00 {sd=     0.031}    125.1878
 ANGLe  HA   CR1H NH1      500.00 {sd=     0.031}    126.4359
 ANGLe  HA   CR1H NR       500.00 {sd=     0.031}    126.4359
 ANGLe  HA   CRH  NH1      500.00 {sd=     0.031}    124.3534
 ANGLe  HA   CRH  NR       500.00 {sd=     0.031}    124.3404
 ANGLe  HA   CRHH NH1      500.00 {sd=     0.031}    125.8381
 ANGLe  HA   CRHH NR       500.00 {sd=     0.031}    125.8381
 ANGLe  HC   NC2  HC       500.00 {sd=     0.031}    120.0016
 ANGLe  HC   NH3  HC       500.00 {sd=     0.031}    108.1992
 ANGLe  N    C    O        500.00 {sd=     0.031}    122.1842  ! 14-MAR-00
 ANGLe  NC2  C    NC2      500.00 {sd=     0.031}    119.7933
 ANGLe  NC2  C    NH1      500.00 {sd=     0.031}    120.1034
 ANGLe  NH1  C    O        500.00 {sd=     0.031}    122.9907
 ANGLe  NH1  CRH  NR       500.00 {sd=     0.031}    111.3061
 ANGLe  NH1  CRHH NH1      500.00 {sd=     0.031}    108.3237
 ANGLe  NH1  CRHH NR       500.00 {sd=     0.031}    108.3237
 ANGLe  NH2  C    O        500.00 {sd=     0.031}    122.6277
 ANGLe  OC   C    OC       500.00 {sd=     0.031}    123.3548
 ANGLe  CCIS  CH1E CH1E     500.00 {sd=     0.031}    109.0754
 ANGLe  CCIS  CH1E CH2E     500.00 {sd=     0.031}    110.1094
 ANGLe  CCIS  CH1E CH3E     500.00 {sd=     0.031}    110.4838
 ANGLe  CCIS  CH1E HA       500.00 {sd=     0.031}    108.9914
 ANGLe  CCIS  CH1E N        500.00 {sd=     0.031}    111.9082
 ANGLe  CCIS  CH1E NH1      500.00 {sd=     0.031}    111.1396
 ANGLe  CCIS  CH1E NH3      500.00 {sd=     0.031}    111.1703
 ANGLe  CCIS  CH2E CH1E     500.00 {sd=     0.031}    112.5947
 ANGLe  CCIS  CH2E CH2E     500.00 {sd=     0.031}    112.5943
 ANGLe  CCIS  CH2E HA       500.00 {sd=     0.031}    108.5877
 ANGLe  CCIS  CH2G HA       500.00 {sd=     0.031}    108.8528
 ANGLe  CCIS  CH2G NH1      500.00 {sd=     0.031}    112.4999
 ANGLe  CCIS  CH2G NH3      500.00 {sd=     0.031}    112.4990
 ANGLe  CCIS  N    CH1E     500.00 {sd=     0.031}    122.7632
 ANGLe  CCIS  N    CH2P     500.00 {sd=     0.031}    125.1134
 ANGLe  CCIS  NC2  HC       500.00 {sd=     0.031}    119.9992
 ANGLe  CCIS  NH1  CH1E     500.00 {sd=     0.031}    121.6541
 ANGLe  CCIS  NH1  CH2E     500.00 {sd=     0.031}    124.1226
 ANGLe  CCIS  NH1  CH2G     500.00 {sd=     0.031}    120.5859
 ANGLe  CCIS  NH1  H        500.00 {sd=     0.031}    119.2489
 ANGLe  CCIS  NH2  H        500.00 {sd=     0.031}    118.1853
 ANGLe  CH1E CCIS  N        500.00 {sd=     0.031}    116.9940 ! CCIS for cis peptide
 ANGLe  CH1E CCIS  NH1      500.00 {sd=     0.031}    116.1998
 ANGLe  CH1E CCIS  O        500.00 {sd=     0.031}    120.8258
 ANGLe  CH1E CCIS  OC       500.00 {sd=     0.031}    118.0611
 ANGLe  CH2E CCIS  NH2      500.00 {sd=     0.031}    116.4617
 ANGLe  CH2E CCIS  O        500.00 {sd=     0.031}    120.9106
 ANGLe  CH2E CCIS  OC       500.00 {sd=     0.031}    118.4969
 ANGLe  CH2G CCIS  N        500.00 {sd=     0.031}    117.7918
 ANGLe  CH2G CCIS  NH1      500.00 {sd=     0.031}    116.3225
 ANGLe  CH2G CCIS  O        500.00 {sd=     0.031}    120.6203
 ANGLe  CH2G CCIS  OC       500.00 {sd=     0.031}    118.4971
 ANGLe  N    CCIS  O        500.00 {sd=     0.031}    122.0016
 ANGLe  NC2  CCIS  NC2      500.00 {sd=     0.031}    119.7933
 ANGLe  NC2  CCIS  NH1      500.00 {sd=     0.031}    120.1034
 ANGLe  NH1  CCIS  O        500.00 {sd=     0.031}    122.9907
 ANGLe  NH2  CCIS  O        500.00 {sd=     0.031}    122.6277
 ANGLe  OC   CCIS  OC       500.00 {sd=     0.031}    123.3548

 ANGLE  O3R  P     O1P      500.00  109.600
 ANGLE  O1P  P     O1P      500.00  120.000
 ANGLE  O3R  P     NR       500.00  109.600
 ANGLE  O1P  P     NR       500.00  120.000
 ANGLE  O3R  P     SM       500.00  109.600
 ANGLE  O1P  P     SM       500.00  120.000
 ANGLE  O3R  P     OHP      500.00  109.600
 ANGLE  O1P  P     OHP      500.00  120.000
 ANGLE  P    NR    CRHH     500.00  120.000
 ANGLE  P    NR    CR1H     500.00  120.000
 
 ANGLe  CH1E CH2E CH1P     500.00 {sd=     0.031}    104.3952 ! for HYP added by JR
 ANGLe  CH2E CH1P CH2P     500.00 {sd=     0.031}    106.1000 ! for HYP added by JR
 ANGLe  CH2E CH1P OH       500.00 {sd=     0.031}    110.0000 ! for HYP added by JR
 ANGLe  CH2E CH1P HA       500.00 {sd=     0.031}    109.9548 ! for HYP added by JR
 ANGLe  CH1P CH2E HA       500.00 {sd=     0.031}    111.1127 ! for HYP added by JR
 ANGLe  CH1P CH2P HA       500.00 {sd=     0.031}    110.3818 ! for HYP added by JR
 ANGLe  CH2P CH1P HA       500.00 {sd=     0.031}    110.3818 ! for HYP added by JR
 ANGLe  CH1P CH2P N        500.00 {sd=     0.031}    103.2695 ! for HYP added by JR
 ANGLe  CH1P OH   H        500.00 {sd=     0.031}    106.0000 ! for HYP added by JR
 ANGLe  HA   CH1P OH       500.00 {sd=     0.031}    104.1000 ! for HYP added by JR
 ANGLe  OH   CH1P CH2P     500.00 {sd=     0.031}    108.8000 ! for HYP added by JR
 ANGLe  CH1E CH2E CQ       500.00 {sd=     0.031}    113.7931 ! for PRG (based on His) added by KJR
 ANGLe  HA   CH2E CQ       500.00 {sd=     0.031}    108.1195 ! for PRG (based on His) added by KJR
 ANGLe  CH2E CQ   CL       500.00 {sd=     0.031}    125.4500 ! for PRG when in triazole bridge added by KJR
 ANGLe  CQ   CL   HA       500.00 {sd=     0.031}    126.6400 ! for PRG when in triazole bridge added by KJR
 ANGLe  CH1E CH2E NK       500.00 {sd=     0.031}    113.7931 ! for AZA (based on His) added by KJR
 ANGLe  CH2E CH2E NK       500.00 {sd=     0.031}    113.7931 ! for AZH (based on His) added by KJR
 ANGLe  HA   CH2E NK       500.00 {sd=     0.031}    108.1195 ! for AZA/AZH (based on His) added by KJR
 ANGLe  CH2E NK   NL       500.00 {sd=     0.031}    126.8900 ! for AZA/AZH when in triozole bridge added by KJR
 ANGLe  NK   NL   NM       500.00 {sd=     0.031}    108.0100 ! for AZA/AZH when in triazole bridge added by KJR
 ANGLe  CH1E NK   CL       500.00 {sd=     0.031}    126.8200 ! for triazole bridge added by KJR
 ANGLe  CH2E NK   CL       500.00 {sd=     0.031}    126.8200 ! for triazole bridge added by KJR
 ANGLe  NK   CL   CQ       500.00 {sd=     0.031}    106.7500 ! for triazole bridge added by KJR
 ANGLe  NL   NM   CQ       500.00 {sd=     0.031}    109.9300 ! for triazole bridge added by KJR
 ANGLe  NM   CQ   CH2E     500.00 {sd=     0.031}    125.5300 ! for triazole bridge added by KJR
 ANGLe  NK   CL   HA       500.00 {sd=     0.031}    126.6200 ! for triazole bridge added by KJR
 ANGLe  CH2K CH1E C        500.00 {sd=     0.031}    110.1000  !!!!included 30-JAN-07 LJ
 ANGLe  CH2K CH2L HA       500.00 {sd=     0.031}    109.5000     !included 30-JAN-07 LJ
 ANGLe  CH2K CH2L CJ       500.00 {sd=     0.031}    106.1000 !!!!!included 30-JAN-07 LJ
 ANGLe  CH2L CJ   O        500.00 {sd=     0.031}    120.4000  !!!!included 30-JAN-07 LJ
 ANGLe  CH2L CJ   N        500.00 {sd=     0.031}    110.0962 !!!!!included 30-JAN-07 LJ
 ANGLe  HA   CH2K HA       500.00 {sd=     0.031}    109.5000   !LJ 31/01/07
 ANGLe  HA   CH2K CH2L     500.00 {sd=     0.031}    109.5000   !LJ 31/01/07
 ANGLe  HA   CH2L HA       500.00 {sd=     0.031}    109.5000   !LJ 31/01/07
 ANGLe  HA   CH2L CJ       500.00 {sd=     0.031}    109.5000   !LJ 31/01/07
 ANGLe  N    CJ   O        500.00 {sd=     0.031}    122.9000   !!!LJ 31/01/07
 ANGLe  N    CH1E CH2K     500.00 {sd=     0.031}    109.5000  !!!!LJ 31/01/07
 ANGLe  HA   CH1E CH2K     500.00 {sd=     0.031}    109.5000   !LJ 31/01/07
 ANGLe  CH1E CH2K HA       500.00 {sd=     0.031}    109.5000    !LJ 31/01/07
 ANGLe  CH1E CH2K CH2L     500.00 {sd=     0.031}    109.5000    !!!LJ 31/01/07
 ANGLe  CH1E N    CJ       500.00 {sd=     0.031}    112.1234   !!!!Included 31-JAN-07 LJ
 ANGLe  CH1E N    H        500.00 {sd=     0.031}    124.0150   !!!!Included 31-JAN-07 LJ
 ANGLe  H    N    CJ       500.00 {sd=     0.031}    124.0150   !!!!Included 31-JAN-07 LJ
 ANGLe  CH1E C    NH2      500.00 {sd=     0.031}    116.4617 !CTN KJR

 
 IMPRoper  C    CH1E HA   HA        500.00 {sd=     0.031}    0    -70.4072
 IMPRoper  C    CH1E OC   OC        500.00 {sd=     0.031}    0      0.0210
 IMPRoper  C    CH2E HA   HA        500.00 {sd=     0.031}    0    -70.4459
 IMPRoper  HA   HA   C    HA        500.00 {sd=     0.031}    0    -70.4459 !EAB
 IMPRoper  C    CH2E O    NH2       500.00 {sd=     0.031}    0      0.0124
 IMPRoper  C    CH2E OC   OC        500.00 {sd=     0.031}    0     -0.0137
 IMPRoper  C    CH2G OC   OC        500.00 {sd=     0.031}    0      0.0223
 IMPRoper  C    NC2  H    NH1       500.00 {sd=     0.031}    0     -0.0121
 IMPRoper  C    NH1  HA   HA        500.00 {sd=     0.031}    0    -70.8745
 IMPRoper  C    NH1  NC2  NC2       500.00 {sd=     0.031}    0     -0.0088
 IMPRoper  C    NH3  HA   HA        500.00 {sd=     0.031}    0     70.6479
 IMPRoper  C5   CH1E HA   HA        500.00 {sd=     0.031}    0    -69.9815
 IMPRoper  C5   CH2E NH1  CR1E      500.00 {sd=     0.031}    0     -0.0386
 IMPRoper  C5   CH2E NH1  CR1H      500.00 {sd=     0.031}    0     -0.0237
 IMPRoper  C5   CH2E NR   CR1E      500.00 {sd=     0.031}    0      0.0178
 IMPRoper  C5   CR1E NH1  CRH       500.00 {sd=     0.031}    0      0.0072
 IMPRoper  C5   CR1E NR   CRH       500.00 {sd=     0.031}    0     -0.0138
 IMPRoper  C5   CR1H NH1  CRHH      500.00 {sd=     0.031}    0     -0.0209
 IMPRoper  C5   CR1H NR   CRHH      500.00 {sd=     0.031}    0     -0.0209
 IMPRoper  C5   NH1  CRH  NR        500.00 {sd=     0.031}    0     -0.0652
 IMPRoper  C5   NH1  CRHH NH1       500.00 {sd=     0.031}    0     -0.0660
 IMPRoper  C5   NH1  CRHH NR        500.00 {sd=     0.031}    0     -0.0660
 IMPRoper  C5   NR   CRH  NH1       500.00 {sd=     0.031}    0      0.0250
 IMPRoper  C5W  CH1E HA   HA        500.00 {sd=     0.031}    0    -70.0142
 IMPRoper  C5W  CW   CR1E CR1E      500.00 {sd=     0.031}    0   -179.9506
 IMPRoper  C5W  CW   CW   CR1W      500.00 {sd=     0.031}    0    179.9618
 IMPRoper  CF   CH1E HA   HA        500.00 {sd=     0.031}    0    -70.0169
 IMPRoper  CF   CR1E CR1E CR1E      500.00 {sd=     0.031}    0      0.0069
 IMPRoper  CF   CR1E CR1E HA        500.00 {sd=     0.031}    0    179.9729
 IMPRoper  CH1E C    CH2P N         500.00 {sd=     0.031}    0      0.0332
 IMPRoper  CH1E C    H    NH1       500.00 {sd=     0.031}    0      0.0000 !c
 IMPRoper  CH1E C    N    CH1E      500.00 {sd=     0.031}    0    179.9873
 IMPRoper  CH3E C    N    CH1E      500.00 {sd=     0.031}    0    179.9873
 IMPRoper  CH1E C    N    CH2P      500.00 {sd=     0.031}    0      0.0025
 IMPRoper  CH1E N    C    CH2G      500.00 {sd=     0.031}    0    179.9856
 IMPRoper  CH1E NH1  C    CH2G      500.00 {sd=     0.031}    0   -179.9916
 IMPRoper  CH1E C    NH1  CH1E      500.00 {sd=     0.031}    0   -180.0067
 IMPRoper  CH1E C    NH1  CH2G      500.00 {sd=     0.031}    0   -180.0018
 IMPRoper  CH1E C    NH1  H         500.00 {sd=     0.031}    0      0.0000 !c
 IMPRoper  CH1E C    NH1  HA        500.00 {sd=     0.031}    0     66.2535 
 IMPRoper  CH1E CCIS NH1  HA        500.00 {sd=     0.031}    0     66.2535 
 IMPRoper  CH1E C    NH3  HA        500.00 {sd=     0.031}    0     66.3265
 IMPRoper  CH1E CH1E HA   HA        500.00 {sd=     0.031}    0    -69.6639
 IMPRoper  CH2E C    N    HA        500.00 {sd=     0.031}    0     67.7957
 IMPRoper  CH2E CCIS N    HA        500.00 {sd=     0.031}    0     67.7957 ! cis pept
 IMPRoper  CH2E C    NH1  HA        500.00 {sd=     0.031}    0     66.1640
 IMPRoper  CB1  C    NH1  CB2       500.00 {sd=     0.031}    0     66.1640 !for PAL JR
 IMPRoper  CH2E CCIS NH1  HA        500.00 {sd=     0.031}    0     66.1640
 IMPRoper  CH2E C    NH2  H         500.00 {sd=     0.031}    0      0.0000
 IMPRoper  CH2E C    NH3  HA        500.00 {sd=     0.031}    0     66.3265
 IMPRoper  CH2E C5W  CW   CW        500.00 {sd=     0.031}    0    179.9679
 IMPRoper  CH2E CF   CR1E CR1E      500.00 {sd=     0.031}    0   -179.9993
 IMPRoper  CH2E CH1E HA   HA        500.00 {sd=     0.031}    0    -70.0781
 IMPRoper  CH2E CH2E HA   HA        500.00 {sd=     0.031}    0    -70.7825
 IMPRoper  CH2E CH3E CH1E HA        500.00 {sd=     0.031}    0    -65.2137
 IMPRoper  CH2E CY   CR1E CR1E      500.00 {sd=     0.031}    0   -179.9903
 IMPRoper  CH2G C    N    CH2P      500.00 {sd=     0.031}    0     -0.0116
 IMPRoper  CH2G C    NH1  CH2G      500.00 {sd=     0.031}    0    179.9899
 IMPRoper  CH2G C    NH1  H         500.00 {sd=     0.031}    0      0.0000 !c
 IMPRoper  CH2P CH1E HA   HA        500.00 {sd=     0.031}    0    -71.9385
 IMPRoper  CH2P CH1E HC   HC        500.00 {sd=     0.031}    0    -70.7727
 IMPRoper  CH2P CH2E HA   HA        500.00 {sd=     0.031}    0    -71.8986
 IMPRoper  CH3E C    NH1  HA        500.00 {sd=     0.031}    0     65.9907
 IMPRoper  CH3E CCIS NH1  HA        500.00 {sd=     0.031}    0     65.9907
 IMPRoper  CH3E C    NH3  HA        500.00 {sd=     0.031}    0     65.6779
 IMPRoper  CH3E CCIS NH3  HA        500.00 {sd=     0.031}    0     65.6779
 IMPRoper  CH3E CH1E HA   HA        500.00 {sd=     0.031}    0    -70.1069
 IMPRoper  CH3E CH3E CH1E HA        500.00 {sd=     0.031}    0    -65.0462
 IMPRoper  CH3E CH3E CH2E HA        500.00 {sd=     0.031}    0    -65.1424
 IMPRoper  CH3E OH1  CH1E HA        500.00 {sd=     0.031}    0     66.1521
 IMPRoper  CR1E C5   NH1  CRH       500.00 {sd=     0.031}    0      0.0557
 IMPRoper  CR1E C5   NR   CRH       500.00 {sd=     0.031}    0     -0.0198
 IMPRoper  CR1E C5W  CW   CR1E      500.00 {sd=     0.031}    0    179.9645
 IMPRoper  CR1E CF   CR1E CR1E      500.00 {sd=     0.031}    0     -0.0034
 IMPRoper  CR1E CF   CR1E HA        500.00 {sd=     0.031}    0   -179.9624
 IMPRoper  CR1E CR1E CR1E CR1E      500.00 {sd=     0.031}    0     -0.0034
 IMPRoper  CR1E CR1E CR1E HA        500.00 {sd=     0.031}    0    179.9935
 IMPRoper  CR1E CR1E CR1W CR1W      500.00 {sd=     0.031}    0     -0.0413
 IMPRoper  CR1E CR1E CR1W HA        500.00 {sd=     0.031}    0   -179.9535
 IMPRoper  CR1E CR1E CW   CW        500.00 {sd=     0.031}    0     -0.0109
 IMPRoper  CR1E CR1E CY   CR1E      500.00 {sd=     0.031}    0      0.0135
 IMPRoper  CR1E CR1E CY2  CR1E      500.00 {sd=     0.031}    0      0.0133
 IMPRoper  CR1E CR1E CY2  OH1       500.00 {sd=     0.031}    0   -179.9788
 IMPRoper  CR1E CR1E CY2  OHP       500.00 {sd=     0.031}    0   -179.9788
 IMPRoper  CR1E CR1W CR1W CW        500.00 {sd=     0.031}    0      0.0360
 IMPRoper  CR1E CR1W CR1W HA        500.00 {sd=     0.031}    0    179.9725
 IMPRoper  CR1E CW   CW   CR1W      500.00 {sd=     0.031}    0      0.0072
 IMPRoper  CR1E CW   CW   NH1       500.00 {sd=     0.031}    0   -179.9720
 IMPRoper  CR1E CY   CR1E HA        500.00 {sd=     0.031}    0   -179.9985
 IMPRoper  CR1E CY2  CR1E HA        500.00 {sd=     0.031}    0    179.9891
 IMPRoper  CR1E NH1  C5   HA        500.00 {sd=     0.031}    0     -0.0175
 IMPRoper  CR1E NH1  CRH  NR        500.00 {sd=     0.031}    0     -0.0206
 IMPRoper  CR1E NH1  CW   CR1W      500.00 {sd=     0.031}    0   -179.9685
 IMPRoper  CR1E NR   C5   HA        500.00 {sd=     0.031}    0     -0.0096
 IMPRoper  CR1E NR   CRH  NH1       500.00 {sd=     0.031}    0      0.0490
 IMPRoper  CR1H C5   NH1  CRHH      500.00 {sd=     0.031}    0      0.0496
 IMPRoper  CR1H C5   NR   CRHH      500.00 {sd=     0.031}    0      0.0496
 IMPRoper  CR1H NH1  C5   HA        500.00 {sd=     0.031}    0      0.0047
 IMPRoper  CR1H NR   C5   HA        500.00 {sd=     0.031}    0      0.0047
 IMPRoper  CR1H NH1  CRHH NH1       500.00 {sd=     0.031}    0      0.0534
 IMPRoper  CR1H NR   CRHH NH1       500.00 {sd=     0.031}    0      0.0534
 IMPRoper  CR1W CR1E CR1E CW        500.00 {sd=     0.031}    0      0.0275
 IMPRoper  CR1W CR1E CR1E HA        500.00 {sd=     0.031}    0   -179.9902
 IMPRoper  CR1W CR1W CR1E HA        500.00 {sd=     0.031}    0    179.9587
 IMPRoper  CR1W CR1W CW   CW        500.00 {sd=     0.031}    0     -0.0194
 IMPRoper  CR1W CR1W CW   NH1       500.00 {sd=     0.031}    0    179.9546
 IMPRoper  CRH  NH1  NR   HA        500.00 {sd=     0.031}    0      0.0429
 IMPRoper  CRH  NR   NH1  HA        500.00 {sd=     0.031}    0     -0.0123
 IMPRoper  CRHH NH1  NH1  HA        500.00 {sd=     0.031}    0      0.0414
 IMPRoper  CRHH NH1  NR   HA        500.00 {sd=     0.031}    0      0.0414
 IMPRoper  CW   CW   NH1  H         500.00 {sd=     0.031}    0    179.9788
 IMPRoper  CW   NH1  CR1E HA        500.00 {sd=     0.031}    0   -179.9528
 IMPRoper  CY   CH1E HA   HA        500.00 {sd=     0.031}    0    -70.0662
 IMPRoper  CY   CR1E CR1E CY2       500.00 {sd=     0.031}    0     -0.0270
 IMPRoper  CY   CR1E CR1E HA        500.00 {sd=     0.031}    0   -179.9841
 IMPRoper  CY2  CR1E CR1E HA        500.00 {sd=     0.031}    0    179.9517
 IMPRoper  H    C    CH2E NH1       500.00 {sd=     0.031}    0      0.0051
 IMPRoper  H    C5   CRH  NH1       500.00 {sd=     0.031}    0      0.0263
 IMPRoper  H    C5   CRHH NH1       500.00 {sd=     0.031}    0      0.0282
 IMPRoper  H    H    C    NH2       500.00 {sd=     0.031}    0      0.0032
 IMPRoper  HA   CH1E HA   HA        500.00 {sd=     0.031}    0    -66.5692
 IMPRoper  HA   CH2E HA   HA        500.00 {sd=     0.031}    0    -66.5934
 IMPRoper  HA   HA   CH1E OH1       500.00 {sd=     0.031}    0    -69.8494
 IMPRoper  HA   HA   CH1E S         500.00 {sd=     0.031}    0    -72.0980
 IMPRoper  HA   HA   CH1E SH1E      500.00 {sd=     0.031}    0    -72.0234
 IMPRoper  HA   HA   CH1E SM        500.00 {sd=     0.031}    0    -72.0234
 IMPRoper  HA   HA   CH2E NH1       500.00 {sd=     0.031}    0    -70.1253
 IMPRoper  HA   HA   CH2E NH3       500.00 {sd=     0.031}    0    -70.4126
 IMPRoper  HA   HA   CH2E SM        500.00 {sd=     0.031}    0    -72.4655
 IMPRoper  HA   HA   CH2P N         500.00 {sd=     0.031}    0    -72.1561
 IMPRoper  HA   HA   CH2P NH3       500.00 {sd=     0.031}    0    -71.9018
 IMPRoper  HA   HA   SM   HA        500.00 {sd=     0.031}    0    -65.1411
 IMPRoper  HC   CH1E HC   HC        500.00 {sd=     0.031}    0    -66.4313
 IMPRoper  HC   CH2E HC   HC        500.00 {sd=     0.031}    0    -66.4262
 IMPRoper  HC   CH2G HC   HC        500.00 {sd=     0.031}    0    -66.4073
 IMPRoper  HC   HC   C    NC2       500.00 {sd=     0.031}    0     -0.0094
 IMPRoper  HC   NC2  C    NH1       500.00 {sd=     0.031}    0      0.0000
 IMPRoper  NH1  C5   CR1E NR        500.00 {sd=     0.031}    0     -0.0249
 IMPRoper  NH1  C5   CR1H NH1       500.00 {sd=     0.031}    0     -0.0178
 IMPRoper  NH1  C5   CR1H NR        500.00 {sd=     0.031}    0     -0.0178
 IMPRoper  NH1  CR1E C5   NR        500.00 {sd=     0.031}    0      0.0078
 IMPRoper  NH1  CR1E C5   NR        500.00 {sd=     0.031}    0      0.0078
 IMPRoper  NR   P    CR1H CRHH      500.00 {sd=     0.031}    0      0.0078
 IMPRoper  NR   P    O3R  O1P       500.00 {sd=     0.031}    0    -35.0000
 IMPRoper  SM   P    O3R  O1P       500.00 {sd=     0.031}    0    -35.0000
 IMPRoper  OHP  P    O3R  O1P       500.00 {sd=     0.031}    0    -35.0000
 IMPRoper  CQ   CH1E HA   HA        500.00 {sd=     0.031}    0    -70.4072 !for PRG
 IMPRoper  NK   CH2E HA   HA        500.00 {sd=     0.031}    0    -70.4072 !for AZA
 IMPRoper  NL   CH1E HA   HA        500.00 {sd=     0.031}    0    -70.4072 !for AZH
 IMPRoper  NK   CH2E NL   CL        500.00 {sd=     0.031}    0       0.000 !for PRG
 IMPRoper  CQ   CH2E CL   NM        500.00 {sd=     0.031}    0       0.000 !for AZA
 IMPRoper  CL   HA   NK   CQ        500.00 {sd=     0.031}    0       0.000 !for AZH
 IMPRoper  NK   CL   CQ   NM        500.00 {sd=     0.031}    0       0.000 !for PRG
 IMPRoper  CL   CQ   NM   NL        500.00 {sd=     0.031}    0       0.000 !for AZA
 IMPRoper  CQ   NM   NL   NK        500.00 {sd=     0.031}    0       0.000 !for AZH
 IMPRoper  NM   NL   NK   CL        500.00 {sd=     0.031}    0       0.000 !for AZA
 IMPRoper  NL   NK   CL   CQ        500.00 {sd=     0.031}    0       0.000 !for AZH

! peptide group (modified 14-MAR-00)
! omega modified again (MN, 18-MAR-02)

 IMPRoper  C    CH1E NH1  O         500.00   0      0.0000 !
 IMPRoper  C    CH3E NH1  O         500.00   0      0.0000 !
 IMPRoper  C    CH1E N    O         500.00   0      0.0000 ! aaa-pro
 IMPRoper  C    CH2G N    O         500.00   0      0.0000 ! gly-pro
 IMPRoper  C    CH2G NH1  O         500.00   0      0.0000 ! gly-aaa
 IMPRoper  NH1  C    CH1E H         500.00   0      0.0000 !
 IMPRoper  N    C    CH1E CH2P      500.00   0      0.0000 ! aaa-pro
 IMPRoper  NH1  C    CH2G H         500.00   0      0.0000 ! aaa-gly

 IMPRoper  CCIS CH1E NH1  O         500.00   0      0.0000 !
 IMPRoper  CCIS CH1E N    O         500.00   0      0.0000 ! aaa-pro
 IMPRoper  CCIS CH2G N    O         500.00   0      0.0000 ! gly-pro
 IMPRoper  CCIS CH2G NH1  O         500.00   0      0.0000 ! gly-aaa
 IMPRoper  CCIS NH1  HA   HA        500.00   0    -78.8745
 IMPRoper  CCIS NH1  N    CH1E      500.00   0      0.0000 
 IMPRoper  NH1  CCIS CH1E H         500.00   0      0.0000 !
 IMPRoper  N    CCIS CH1E CH2P      500.00   0      0.0000 ! aaa-pro
 IMPRoper  NH1  CCIS CH2G H         500.00   0      0.0000 ! aaa-gly

 IMPRoper  CH1E C    NH1  CH1E      500.00   0    180.0000 
 IMPRoper  CH3E C    NH1  CH1E      500.00   0    180.0000 
 IMPRoper  CH1E C    NH1  CH2G      500.00   0    180.0000 
 IMPRoper  CH1E C    N    CH2P      500.00   0    180.0000 
 IMPRoper  CH2G C    NH1  CH1E      500.00   0    180.0000 
 IMPRoper  CH2G C    NH1  CH2G      500.00   0    180.0000 
 IMPRoper  CH2G C    N    CH2P      500.00   0    180.0000 
 IMPRoper  CH2P C    NH1  CH1E      500.00   0    180.0000 
 IMPRoper  CH2P C    NH1  CH2G      500.00   0    180.0000 
 IMPRoper  CH2P C    NH1  CH2P      500.00   0    180.0000 

 IMPRoper  CH1E CCIS NH1  CH1E      500.00   0      0.0000 
 IMPRoper  CH1E CCIS NH1  CH2G      500.00   0      0.0000 
 IMPRoper  CH1E CCIS N    CH1E      500.00   0      0.0000 
 IMPRoper  CH1E CCIS N    CH2P      500.00   0      0.0000 
 IMPRoper  CH2G CCIS NH1  CH1E      500.00   0      0.0000 
 IMPRoper  CH2G CCIS N    CH1E      500.00   0      0.0000 
 IMPRoper  CH2G CCIS NH1  CH2G      500.00   0      0.0000 
 IMPRoper  CH2G CCIS N    CH2P      500.00   0      0.0000 
 IMPRoper  CH2P CCIS NH1  CH1E      500.00   0      0.0000 
 IMPRoper  CH2P CCIS NH1  CH2G      500.00   0      0.0000 
 IMPRoper  CH2P CCIS NH1  CH2P      500.00   0      0.0000 

 IMPRoper  CH1P CH1E HA   HA        500.00 {sd=     0.031}    0    -71.9385 ! for HYP
 IMPRoper  HA   HA   CH1P N         500.00 {sd=     0.031}    0    -72.1561 ! for HYP
 IMPRoper  OH   HA   CH2E CH2P      500.00 {sd=     0.031}    0    -70.8740 ! for HYP
 IMPRoper  C    CH1E NH2  O         500.00 {sd=     0.031}    0      0.0000 !for CTN
 IMPRoper  C    NH2  O    H         500.00 {sd=     0.031}    0      0.0000 ! for CTN
 IMPRoper  NH2  H    H    C         500.00 {sd=     0.031}    0      0.0000 ! for CTN
 IMPRoper  N    CJ   H    CH1E      500.00 {sd=     0.031}    0      0.0032 !LJ
 IMPRoper  NH2  H    H    C         500.00 {sd=     0.031}    0      0.0032 !LJ
 IMPRoper  CH2L CH1E HA   HA        500.00 {sd=     0.031}    0    -71.9385 !LJ
 IMPRoper  CH2K C    N    HA        500.00 {sd=     0.031}    0     67.7957 !LJ
 IMPRoper  C    NH2  O    H         500.00 {sd=     0.031}    0      0.0000 !LJ
 IMPRoper  C    CH2K HA   HA        500.00 {sd=     0.031}    0    -70.4459 !LJ
 IMPRoper  CJ   CH2K HA   HA        500.00 {sd=     0.031}    0    -70.4459 !LJ
 IMPRoper  CJ   CH2L O    N         500.00 {sd=     0.031}    0      0.0124 !LJ
 IMPRoper  C    CH1E NH2  O         500.00 {sd=     0.031}    0      0.0124 !LJ

! CR3-CR3 sidechain dihedrals

 DIHEdral  C    CH2E CH2E CH1E   2.00     3  0.0000
 DIHEdral  CH1E CH1E CH2E CH3E   2.00     3  0.0000
 DIHEdral  CH1E CH2E CH2E CH2E   2.00     3  0.0000
 DIHEdral  CH1E CH2E CH2E SM     2.00     3  0.0000
 DIHEdral  CH2E CH1E CH1E NH1    2.00     3  0.0000
 DIHEdral  CH2E CH2E CH2E CH2E   2.00     3  0.0000
 DIHEdral  NH3  CH2E CH2E CH2E   2.00     3  0.0000 
 DIHEdral  CH2E CH2E CH2E NH1    2.00     3  0.0000
 DIHEdral  CH2E CH2E SM   CH3E   2.00     3  0.0000
 DIHEdral  CH3E CH1E CH1E NH1    2.00     3  0.0000
 DIHEdral  CH3E CH1E CH2E CH1E   2.00     3  0.0000
 DIHEdral  NH1  CH1E CH2E C      2.00     3  0.0000
 DIHEdral  NH1  CH1E CH2E C5     2.00     3  0.0000
 DIHEdral  NH1  CH1E CH2E C5W    2.00     3  0.0000
 DIHEdral  NH1  CH1E CH2E CF     2.00     3  0.0000
 DIHEdral  NH1  CH1E CH2E CH1E   2.00     3  0.0000
 DIHEdral  NH1  CH1E CH2E CH2E   2.00     3  0.0000
 DIHEdral  NH1  CH1E CH2E CY     2.00     3  0.0000
 DIHEdral  NH1  CH1E CH2E OH1    2.00     3  0.0000
 DIHEdral  NH1  CH1E CH1E OH1    2.00     3  0.0000 !added 14-APR-99
 DIHEdral  NH1  CH1E CH2E S      2.00     3  0.0000
 DIHEdral  NH1  CH1E CH2E SH1E   2.00     3  0.0000
 DIHEdral  NH1  CH1E CH2E SM     2.00     3  0.0000
 DIHEdral  CH1E CH2E CH2E NK     2.00     3  0.0000 ! for AZH
 DIHEdral  NH1  CH1E CH2E NL     2.00     3  0.0000 ! for AZA
 DIHEdral  NH1  CH1E CH2E CQ     2.00     3  0.0000 ! for AZA
 ! pyroglutamic acid
 DIHEdral  CJ   N     CH1E   CH2K   750.00 {sd=     0.031}    0      0.0000 !LJ 07-03-15 
 DIHEdral  CJ   N     CH1E   C      750.00 {sd=     0.031}    0    120.0000
 DIHEdral  CH1E N     CJ     CH2L   750.00 {sd=     0.031}    0      0.0000
 DIHEdral  CH1E N     CJ     O      750.00 {sd=     0.031}    0    180.0000
 DIHEdral  N    CH1E  CH2K   CH2L   750.00 {sd=     0.031}    0      0.0000
 DIHEdral  C    CH1E  CH2K   CH2L   750.00 {sd=     0.031}    0   -120.0000
 DIHEdral  CH2K CH1E  C      O      750.00 {sd=     0.031}    0    -90.0000
 DIHEdral  CH1E CH2K  CH2L   CJ     750.00 {sd=     0.031}    0      0.0000
 DIHEdral  CH2K CH2L  CJ     N      750.00 {sd=     0.031}    0      0.0000
 DIHEdral  CH2K CH2L  CJ     O      750.00 {sd=     0.031}    0    180.0000


! chi1 modifications at N-terminus

 DIHEdral  NH3  CH1E CH2E C      2.00     3  0.0000 
 DIHEdral  NH3  CH1E CH2E C5     2.00     3  0.0000
 DIHEdral  NH3  CH1E CH2E C5W    2.00     3  0.0000
 DIHEdral  NH3  CH1E CH2E CF     2.00     3  0.0000
 DIHEdral  NH3  CH1E CH2E CH1E   2.00     3  0.0000
 DIHEdral  NH3  CH1E CH2E CH2E   2.00     3  0.0000
 DIHEdral  NH3  CH1E CH2E CY     2.00     3  0.0000
 DIHEdral  NH3  CH1E CH2E OH1    2.00     3  0.0000
 DIHEdral  NH3  CH1E CH1E OH1    2.00     3  0.0000
 DIHEdral  NH3  CH1E CH2E S      2.00     3  0.0000
 DIHEdral  NH3  CH1E CH2E SH1E   2.00     3  0.0000
 DIHEdral  NH3  CH1E CH1E CH3E   2.00     3  0.0000
 DIHEdral  NH3  CH1E CH1E CH2E   2.00     3  0.0000

! CR3-CR2 (Planar) sidechain dihedrals (modified 09-MAR-00)

 DIHEdral  CR1E C5W  CH2E CH1E   1.00     2  0.0000 !Trp
 DIHEdral  CR1E CF   CH2E CH1E   1.00     2  0.0000 !Phe
 DIHEdral  NH1  C5   CH2E CH1E   0.60     1  0.0000 !His
 DIHEdral  CR1H C5   CH2E CH1E   0.80     2  0.0000 !His
 DIHEdral  NR   C5   CH2E CH1E   0.60     1  0.0000 !His 
 DIHEdral  CR1E C5   CH2E CH1E   0.80     2  0.0000 !His
 DIHEdral  CH1E CH2E CY   CR1E   1.00     2  0.0000 !Tyr
 DIHEDral  CR1E CY2  OH1  H      2.00     2  90.000 !Tyr
 DIHEDral  CR1E CY2  OHP  P      2.00     2  90.000 !Ptr
 DIHEdral  CH2E CH2E NH1  C      2.00     1  0.0000 !Arg 
 DIHEdral  O    C    CH2E CH1E   1.00     6  0.0000 !Asn (+ 3 60.000?) 
 DIHEdral  OC   C    CH2E CH1E   MULT 2 1.00     2  0.0000  0.5     6  0.0000 !Asp
 DIHEdral  O    C    CH2E CH2E   2.00     6  0.0000 !Gln (+ 3 60.000?) 
 DIHEdral  OC   C    CH2E CH2E   MULT 2 1.00     2  0.0000  0.5     6  0.0000 !Glu 

! sidechain CH3 and NH3 groups (added 12-MAR-00)

 DIHEdral  HC   NH3  CH2E CH2E   1.00     3  0.0000 !Lys 
 DIHEdral  HA   CH3E CH1E C      1.00     3  0.0000 !Ala 
 DIHEdral  HA   CH3E CH2E CH1E   1.00     3  0.0000 !Ile 
 DIHEdral  HA   CH3E CH1E CH1E   1.00     3  0.0000 !Ile/Thr/Val  
 DIHEdral  HA   CH3E CH1E CH2E   1.00     3  0.0000 !Leu  
 DIHEdral  HA   CH3E SM   CH2E   1.00     3  0.0000 !Met 

 ! added by Aart Nederveen July 3 2003
 DIHEdral  HA   CH3E CH1E CCIS      1.00     3  0.0000 !Ala

! N and C terminal groups (ace planarity taken care of by peptide bond added 12-MAR-00)

 DIHEdral  HC   NH3  CH1E CCIS   2.00     3  0.0000 !Non-Gly/Pro N-terminus
 DIHEdral  HC   NH3  CH1E C      2.00     3  0.0000 !Non-Gly/Pro N-terminus
 DIHEdral  HC   NH3  CH2G C      2.00     3  0.0000 !Gly N-terminus
 DIHEdral  OC   C    CH1E NH1    2.00     2  0.0000 !Non-Gly/Pro COO-terminus 
 DIHEdral  OC   C    CH1E NH3    2.00     2  0.0000 !Non-Gly/Pro COO-terminus 
 DIHEdral  OC   C    CH2G NH1    1.00     6  0.0000 !Gly COO-terminus 
 DIHEdral  OC   C    CH1E N      2.00     2  0.0000 !Pro COO-terminus
 DIHEdral  NH2  C    CH1E NH1    2.00     2  0.0000 !Non-Gly/Pro CONH2-terminus
 DIHEdral  NH2  C    CH2G NH1    1.00     6  0.0000 !Gly CONH2-terminus (+ 3 60.000?)
 DIHEdral  NH2  C    CH1E N      2.00     2  0.0000 !Pro CONH2-terminus

! backbone psi (not active with usual topology file modified 03-JUL-01)
! DIHEdral NH1    CH1E   C      NH1    MULT 3   0.175  2 0.0000   0.1  4 0.0000   0.03  6 0.0000 !aaa-aaa
! DIHEdral NH1    CH2G   C      NH1    MULT 3   0.175  2 0.0000   0.1  4 0.0000   0.03  6 0.0000 !gly-aaa
! DIHEdral NH3    CH1E   C      NH1    MULT 3   0.175  2 0.0000   0.1  4 0.0000   0.03  6 0.0000 !Nterm-aaa-aaa
! DIHEdral NH3    CH2G   C      NH1    MULT 3   0.175  2 0.0000   0.1  4 0.0000   0.03  6 0.0000 !Nterm-gly-aaa
! DIHEdral N      CH1E   C      NH1    MULT 3   0.175  2 0.0000   0.1  4 0.0000   0.03  6 0.0000 !pro-aaa
 
! DIHEdral NH1    CH1E   C      N      MULT 3   0.175  2 0.0000   0.1  4 0.0000   0.03  6 0.0000 !aaa-pro
! DIHEdral NH1    CH2G   C      N      MULT 3   0.175  2 0.0000   0.1  4 0.0000   0.03  6 0.0000 !gly-pro
! DIHEdral NH3    CH1E   C      N      MULT 3   0.175  2 0.0000   0.1  4 0.0000   0.03  6 0.0000 !Nterm-aaa-pro
! DIHEdral NH3    CH2G   C      N      MULT 3   0.175  2 0.0000   0.1  4 0.0000   0.03  6 0.0000 !Nterm-gly-pro
! DIHEdral N      CH1E   C      N      MULT 3   0.175  2 0.0000   0.1  4 0.0000   0.03  6 0.0000 !pro-pro

! backbone psi related due to O & CB (not active with usual topology file added 07-JUL-01)
 DIHEdral CH2K   CH1E   C      O      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !KJR PGL
 DIHEdral CH2E   CH1E   C      O      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !aaa-aaa
 DIHEdral CH3E   CH1E   C      O      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !ala-aaa
 DIHEdral CH1E   CH1E   C      O      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !(ile,thr,val)-aaa

! backbone psi related due to O & CB (not active with usual topology file added 07-JUL-01)
 DIHEdral CH2E   CH1E   CCIS   O      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !aaa-aaa
 DIHEdral CH3E   CH1E   CCIS   O      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !ala-aaa
 DIHEdral CH1E   CH1E   CCIS   O      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !(ile,thr,val)-aaa

! backbone psi related due to NH & CB (N.B. 2.0x the strength of  above not active with usual topology file added 03-JUL-01 modified 07-JUL-01)
 DIHEdral CH2K   CH1E   C      NH1    MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !KJR PGL
 DIHEdral CH2E   CH1E   C      NH1    MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !aaa-aaa
 DIHEdral CH3E   CH1E   C      NH1    MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !ala-aaa
 DIHEdral CH1E   CH1E   C      NH1    MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !(ile,thr,val)-aaa
 DIHEdral CH2E   CH1E   CCIS   NH1    MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !aaa-aaa
 DIHEdral CH3E   CH1E   CCIS   NH1    MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !ala-aaa
 DIHEdral CH1E   CH1E   CCIS   NH1    MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !(ile,thr,val)-aaa

 DIHEdral CH2E   CH1E   C      N      MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !aaa-pro
 DIHEdral CH3E   CH1E   C      N      MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !ala-pro
 DIHEdral CH1E   CH1E   C      N      MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !(ile,thr,val)-pro

 DIHEdral CH2E   CH1E   CCIS   N      MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !aaa-pro
 DIHEdral CH3E   CH1E   CCIS   N      MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !ala-pro
 DIHEdral CH1E   CH1E   CCIS   N      MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !(ile,thr,val)-pro

! backbone phi (not active with usual topology file modified 03-JUL-01)
 DIHEdral C      NH1    CH1E   C      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !aaa-aaa
 DIHEdral C      NH1    CH2G   C      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !aaa-gly
! backbone phi (not active with usual topology file modified 03-JUL-01)
 DIHEdral CCIS      NH1    CH1E   C      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !aaa-aaa
 DIHEdral CCIS      NH1    CH2G   C      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !aaa-gly
! backbone phi (not active with usual topology file modified 03-JUL-01)
 DIHEdral C      NH1    CH1E   CCIS      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !aaa-aaa
 DIHEdral C      NH1    CH2G   CCIS      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !aaa-gly
! backbone phi (not active with usual topology file modified 03-JUL-01)
 DIHEdral CCIS      NH1    CH1E   CCIS      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !aaa-aaa
 DIHEdral CCIS      NH1    CH2G   CCIS      MULT 6   0.375  1 0.0000 0.3  2 0.0000  0.225  3 0.0000  0.15  4 0.0000  0.075  5 0.0000  0.0375  6 0.0000 !aaa-gly

! backbone phi related due to O & CB  N.B. 1.33x the strength of phi above (not active with usual topology file added 03-JUL-01)
 DIHEdral C      NH1    CH1E   CH2E   MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !aaa-aaa
 DIHEdral C      NH1    CH1E   CH3E   MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !aaa-ala
 DIHEdral C      NH1    CH1E   CH1E   MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !aaa-(ile,thr,val)
 DIHEdral CCIS   NH1    CH1E   CH2E   MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !aaa-aaa
 DIHEdral CCIS   NH1    CH1E   CH3E   MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !aaa-ala
 DIHEdral CCIS   NH1    CH1E   CH1E   MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.3   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.05   6 0.0000 !aaa-(ile,thr,val)

! disulphide (modified 07-JUL-01)
   
 DIHEdral CH2E   S      S      CH2E   MULT 6   0.5   1 0.0000 0.4  2 0.0000  0.4   3 0.0000  0.2  4 0.0000  0.1   5 0.0000  0.075  6 0.0000 !aaa-aaa
! DIHEdral CH2E   S      S    CH2E     0.5       2 0.0000  

! Phosphorylated Histidine (AB 1-05-03)
 DIHEdral O3R    P      NR     CRHH   2.0  2  0.0

! Phosphorylated Cysteine (AB 7-05-03)
 DIHEdral O3R    P      SM     CH2E   2.0  2  0.0

! Phosphorylated Tyrosine (AB 8-05-06)
 DIHEdral O3R    P      OHP     CY2   2.0  2  0.0
 
 
if ( $exist_par_nonbonded = false ) then
  evaluate ($par_nonbonded = "OPLSX")
end if


{* nonbonding parameter section *}
{* ============================ *}

 if ($par_nonbonded = "CONTACT") then                   { added by 10-MAR-00 }

!  This uses distances consistent with those in a survey of atom-atom contacts 
!  in high resolution crystal structures (Williams et al. (1994) Protein Science 3, 1224-1235)
!  implemented with the PROLSQ form:
!
!    fVDW(R) =  16 * ( Rmin - R ) ^ 4  
!
!  The Rmin is related to the Lennard-Jones (epsilon,sigma) form by epsilon=0.1(not used) 
!  and Rmin = sigma * 2 ^ (1/6).
!
{ mandatory values:
  NBONds
    CUTNB=7.0   WMIN=1.5
    REPEl = 1.0          
    REXPonent = 4
    IREXponent = 1
    RCONst = 16.0
    TOLErance = 0.5      NBXMOD = 5
    ctonnb=5.5 ctofnb=6.0 {* for consistency only, not needed for repel *}
  END
}
 evaluate ($repel_radius = 1.0)
 evaluate ($repel_rcons = 25)
 evaluate ($repel_rexpo  = 4)
 evaluate ($repel_irexp  = 1)

!         type          sigma     
evaluate ($VR_C=        3.15)    
evaluate ($VR_N=        2.6)     {polar value to reflect H-bonding ability} 
evaluate ($VR_O=        2.4)     {polar value to reflect H-bonding ability}                       
evaluate ($VR_S=        3.1)     
evaluate ($VR_FE=       2.15)
evaluate ($VR_H=        1.85)     {large for aliphatics}                                                     
evaluate ($VR_HH=       0.65)     {small for H bonds}                      
evaluate ($VR_P=        3.4)
evaluate ($VR_I=        3.9)

! compute 1-4 sigmas 
evaluate ($VR14_C   = $VR_C    -0.27)
evaluate ($VR14_N   = $VR_N    -0.27)
evaluate ($VR14_O   = $VR_O    -0.27)
evaluate ($VR14_S   = $VR_S    -0.27)
evaluate ($VR14_FE  = $VR_FE   -0.27)
evaluate ($VR14_H   = $VR_H    -0.27)
evaluate ($VR14_HH  = $VR_H    -0.27)  {no 1-4 H bonds are possible}
evaluate ($VR14_P   = $VR_P    -0.27)
evaluate ($VR14_I   = $VR_I    -0.27)

evaluate ($VE=0.1)
!
 !                  eps     sigma              eps(1:4) sigma(1:4)
 !                  (kcal/mol) (A)             (kcal/mol) (A)
 !                   --------------------------------------------------
 NONBonded  H       $VE $VR_HH   $VE $VR14_HH 
 NONBonded  HA      $VE $VR_H    $VE $VR14_H   ! aliphatic hydrogen
 NONBonded  HC      $VE $VR_HH   $VE $VR14_HH  ! h attached to charg.
 !
 NONBonded  C       $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CCIS    $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  C5      $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  C5W     $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CF      $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CW      $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CY      $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CY2     $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CH1P    $VE $VR_C    $VE $VR14_C         ! KJR HYP
 NONBonded  CH1E    $VE $VR_C    $VE $VR14_C         ! \
 NONBonded  CH2E    $VE $VR_C    $VE $VR14_C         !  extended carbons
 NONBonded  CH2G    $VE $VR_C    $VE $VR14_C         !  extended carbons
 NONBonded  CH2P    $VE $VR_C    $VE $VR14_C         !  extended carbons
 NONBonded  CH3E    $VE $VR_C    $VE $VR14_C         ! /
 NONBonded  CR1E    $VE $VR_C    $VE $VR14_C         !  ring carbons
 NONBonded  CH2K    $VE $VR_C    $VE $VR14_C         !  !KJR PGL
 NONBonded  CH2L    $VE $VR_C    $VE $VR14_C         !  !KJR PGL
 NONBonded  CJ      $VE $VR_C    $VE $VR14_C         !  !KJR PGL
 NONBonded  CR1H    $VE $VR_C    $VE $VR14_C         !  ring carbons
 NONBonded  CR1W    $VE $VR_C    $VE $VR14_C         !  ring carbons
 NONBonded  CRHH    $VE $VR_C    $VE $VR14_C         !  ring carbons
 NONBonded  CRH     $VE $VR_C    $VE $VR14_C         !  ring carbons
		    
 !
 NONBonded  N       $VE $VR_N    $VE $VR14_N         
 NONBonded  NC2     $VE $VR_N    $VE $VR14_N
 NONBonded  NH1     $VE $VR_N    $VE $VR14_N
 NONBonded  NH2     $VE $VR_N    $VE $VR14_N
 NONBonded  NH3     $VE $VR_N    $VE $VR14_N
 NONBonded  NP      $VE $VR_N    $VE $VR14_N
 NONBonded  NR      $VE $VR_N    $VE $VR14_N
 !
 NONBonded  O       $VE $VR_O    $VE $VR14_O
 NONBonded  OC      $VE $VR_O    $VE $VR14_O
 NONBonded  OH1     $VE $VR_O    $VE $VR14_O
 NONBonded  OH      $VE $VR_O    $VE $VR14_O
 NONBonded  OHP     $VE $VR_O    $VE $VR14_O
 !
 NONBonded  S       $VE $VR_S    $VE $VR14_S
 NONBonded  SM      $VE $VR_S    $VE $VR14_S
 NONBonded  SH1E    $VE $VR_S    $VE $VR14_S


 elseif ($par_nonbonded eq "PROLSQ") then
{* nonbonding parameter section *}
{* ============================ *}
!!
!  This uses a new form of the REPEL function:
!    fVDW(R) =  RCON *( Rmin ^ IREX - R  ^ IREX ) ^ REXP
!
!  PROLSQ uses a function of the form:
!    fVDW(R) =  (1 / 0.5) ^ 4 * ( Rmin ^ 4 - R  ^ 4 ) 
!
!  The epsilon values are arbitrary since the repel function does not depend
!   on epsilon.  The sigma values come from converting the Van der Waals 
!   radii of the PROLSQ program into sigma values using the formula:
!     Rmin = sigma * 2 ^ (1/6)
!   Note:  Prolsq decrements Van der Waals radii for non-bonded contacts
!           that involve torsion angles (1:4 contacts) by .30 A, and
!           hydrogen bonds (X...Y) by .2 A (X-H...Y) by .9.  The former
!           decrement is accomplished in CNS by using the
!           1-4 nonbonded terms.  The latter decrement is accomplished by
!           decreasing the van der Waals radius of hydrogens by 0.8 
!           and that of O and N by 0.1 A.
!
{ suggested values:
  NBONds
    CUTNB=7.0   WMIN=1.5
    REPEl = 1.0          
    REXPonent = 4
    IREXponent = 1
    RCONst = 16.0
    TOLErance = 0.5      NBXMOD = 5
    ctonnb=5.5 ctofnb=6.0 {* for consistency only, not needed for repel *}
  END
}
 evaluate ($repel_radius = 1.0)
 evaluate ($repel_rcons = 25)
 evaluate ($repel_rexpo  = 4)
 evaluate ($repel_irexp  = 1)
!                  type      van der Waals radius         correction applied for hbond
evaluate ($VR_C=        3.7)   
evaluate ($VR_N=        3.0)                             {-0.1}
evaluate ($VR_O=        2.9)                             {-0.1}
evaluate ($VR_S=        3.6)
evaluate ($VR_FE=       2.4)
evaluate ($VR_H=        2.0) {from 1.6 02-APRIL-00}                             
evaluate ($VR_HH=       1.6)                             {-0.8}
evaluate ($VR_P=        3.8)
evaluate ($VR_I=        4.3)
evaluate ($VR_C_SP2=    3.4)

{ convert radii into sigmas }
!
! sigma= vdw radius / 2 ^ (1/6)
!
evaluate ($VR_C   = $VR_C    / 2^(1/6))
evaluate ($VR_N   = $VR_N    / 2^(1/6))
evaluate ($VR_O   = $VR_O    / 2^(1/6))
evaluate ($VR_S   = $VR_S    / 2^(1/6))
evaluate ($VR_FE  = $VR_FE   / 2^(1/6))
evaluate ($VR_H   = $VR_H    / 2^(1/6))
evaluate ($VR_HH  = $VR_HH   / 2^(1/6))
evaluate ($VR_P   = $VR_P    / 2^(1/6))
evaluate ($VR_I   = $VR_I    / 2^(1/6))
evaluate ($VR_C_SP2=$VR_C_SP2/ 2^(1/6))

{ compute 1-4 sigmas from reduced radii}
evaluate ($VR14_C   = $VR_C    -(0.3/ 2^(1/6)))
evaluate ($VR14_N   = $VR_N    -(0.3/ 2^(1/6)))
evaluate ($VR14_O   = $VR_O    -(0.3/ 2^(1/6)))
evaluate ($VR14_S   = $VR_S    -(0.3/ 2^(1/6)))
evaluate ($VR14_FE  = $VR_FE   -(0.3/ 2^(1/6)))
evaluate ($VR14_H   = $VR_H    -(0.3/ 2^(1/6)))
evaluate ($VR14_HH  = $VR_HH   -(0.3/ 2^(1/6)))
evaluate ($VR14_P   = $VR_P    -(0.3/ 2^(1/6)))
evaluate ($VR14_I   = $VR_I    -(0.3/ 2^(1/6)))
evaluate ($VR14_C_SP2=$VR_C_SP2-(0.3/ 2^(1/6)))

evaluate ($VE=0.1)
!
 !                  eps     sigma              eps(1:4) sigma(1:4)
 !                  (kcal/mol) (A)             (kcal/mol) (A)
 !                   --------------------------------------------------
 NONBonded  H       $VE $VR_HH   $VE $VR14_HH 
 NONBonded  HA      $VE $VR_H    $VE $VR14_H   ! aliphatic hydrogen
 NONBonded  HC      $VE $VR_HH   $VE $VR14_HH  ! h attached to charg.
 !
 NONBonded  C       $VE $VR_C_SP2    $VE $VR14_C_SP2 ! carbonyl carbon
 NONBonded  CCIS    $VE $VR_C_SP2    $VE $VR14_C_SP2 ! carbonyl carbon
 NONBonded  C5      $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  C5W     $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CF      $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CW      $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CY      $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CY2     $VE $VR_C    $VE $VR14_C         ! carbonyl carbon
 NONBonded  CH1P    $VE $VR_C    $VE $VR14_C         ! KJR HYP
 NONBonded  CH1E    $VE $VR_C    $VE $VR14_C         ! \
 NONBonded  CH2E    $VE $VR_C    $VE $VR14_C         !  extended carbons
 NONBonded  CH2G    $VE $VR_C    $VE $VR14_C         !  extended carbons
 NONBonded  CH2P    $VE $VR_C    $VE $VR14_C         !  extended carbons
 NONBonded  CH3E    $VE $VR_C    $VE $VR14_C         ! /
 NONBonded  CR1E    $VE $VR_C    $VE $VR14_C         !  ring carbons
 NONBonded  CH2K    $VE $VR_C    $VE $VR14_C         !  !KJR PGL
 NONBonded  CH2L    $VE $VR_C    $VE $VR14_C         !  !KJR PGL
 NONBonded  CJ      $VE $VR_C    $VE $VR14_C         !  !KJR PGL
 NONBonded  CR1H    $VE $VR_C    $VE $VR14_C         !  ring carbons
 NONBonded  CR1W    $VE $VR_C    $VE $VR14_C         !  ring carbons
 NONBonded  CRHH    $VE $VR_C    $VE $VR14_C         !  ring carbons
 NONBonded  CRH     $VE $VR_C    $VE $VR14_C         !  ring carbons
		    
 !
 NONBonded  N       $VE $VR_N    $VE $VR14_N         
 NONBonded  NC2     $VE $VR_N    $VE $VR14_N
 NONBonded  NH1     $VE $VR_N    $VE $VR14_N
 NONBonded  NH2     $VE $VR_N    $VE $VR14_N
 NONBonded  NH3     $VE $VR_N    $VE $VR14_N
 NONBonded  NP      $VE $VR_N    $VE $VR14_N
 NONBonded  NR      $VE $VR_N    $VE $VR14_N
 !
 NONBonded  O       $VE $VR_O    $VE $VR14_O
 NONBonded  O1P     $VE $VR_O    $VE $VR14_O
 NONBonded  O3R     $VE $VR_O    $VE $VR14_O
 NONBonded  OC      $VE $VR_O    $VE $VR14_O
 NONBonded  OH1     $VE $VR_O    $VE $VR14_O
 NONBonded  OHP     $VE $VR_O    $VE $VR14_O
 NONBonded  OH      $VE $VR_O    $VE $VR14_O
 !
 NONBonded  S       $VE $VR_S    $VE $VR14_S
 NONBonded  SM      $VE $VR_S    $VE $VR14_S
 NONBonded  SH1E    $VE $VR_S    $VE $VR14_S

 

 elseif ($par_nonbonded eq "PARMALLH6") then
{ suggested values for refinement in vacuo
  NBONds
    CUTNB=10   WMIN=1.5
    REPEl = 0.0     rdie vswitch     
    TOLErance = 0.5      NBXMOD = 5
    ctonnb=5.5 ctofnb=9.0 
  END
  HBOND 
    AEXP=4 REXP=6 HAEX=4 AAEX=2   ACCEPTORS
    CUTHB=4.5 CTOFHB=4.0 CTONHB=3.5  CUTHA=90.0  CTOFHA=70.0  CTONHA=50.0
  end

or

{ suggested values:
  NBONds
  CUTNB=7.0   WMIN=1.5
  REPEl = 0.8          
  REXPonent = 2
  IREXponent = 2
  RCONst = 5.0
  TOLErance = 0.5      NBXMOD = 5
  ctonnb=5.5 ctofnb=6.0 {* for consistency only, not needed for repel *}
END
}
}
 evaluate ($repel_radius = 0.8)
 evaluate ($repel_rcons = 5.0)
 evaluate ($repel_rexpo  = 2)
 evaluate ($repel_irexp  = 2)
 NONBonded  C       0.0903   3.2072      0.0903   3.2072
 NONBonded  CCIS    0.0903   3.2072      0.0903   3.2072
 NONBonded  CR1E    0.0903   3.2072      0.0903   3.2072
 NONBonded  CF      0.0903   3.2072      0.0903   3.2072
 NONBonded  CY      0.0903   3.2072      0.0903   3.2072
 NONBonded  CY2     0.0903   3.2072      0.0903   3.2072
 NONBonded  CR1W    0.0903   3.2072      0.0903   3.2072
 NONBonded  CW      0.0903   3.2072      0.0903   3.2072
 NONBonded  C5W     0.0903   3.2072      0.0903   3.2072
 NONBonded  CN      0.0903   3.2072      0.0903   3.2072
 NONBonded  C5      0.0903   3.2072      0.0903   3.2072
 NONBonded  CH1E    0.0903   3.2072      0.0903   3.2072
 NONBonded  CH1P    0.0903   3.2072      0.0903   3.2072 !KJR HYP
 NONBonded  CH2E    0.0903   3.2072      0.0903   3.2072
 NONBonded  CH2K    0.0903   3.2072      0.0903   3.2072 !KJR PGL
 NONBonded  CH2L    0.0903   3.2072      0.0903   3.2072 !KJR PGL
 NONBonded  CJ      0.0903   3.2072      0.0903   3.2072 !KJR PGL
 NONBonded  CH3E    0.0903   3.2072      0.0903   3.2072
 NONBonded  CH2G    0.0903   3.2072      0.0903   3.2072
 NONBonded  CH2P    0.0903   3.2072      0.0903   3.2072
 NONBonded  CRH     0.0903   3.2072      0.0903   3.2072
 NONBonded  CR1H    0.0903   3.2072      0.0903   3.2072
 NONBonded  CRHH    0.0903   3.2072      0.0903   3.2072
 NONBonded  H       0.0498   1.4254      0.0498   1.4254
 NONBonded  HA      0.0045   2.6157      0.0045   2.6157
 NONBonded  HC      0.0498   1.4254      0.0498   1.4254
 NONBonded  N       0.1592   2.7618      0.1592   2.7618
 NONBonded  NR      0.1592   2.7618      0.1592   2.7618
 NONBonded  NH1     0.1592   2.7618      0.1592   2.7618
 NONBonded  NH2     0.1592   2.7618      0.1592   2.7618
 NONBonded  NH3     0.1592   2.7618      0.1592   2.7618
 NONBonded  NC2     0.1592   2.7618      0.1592   2.7618
 NONBonded  O       0.2342   2.6406      0.2342   2.6406
 NONBonded  OC      1.0244   2.6406      1.0244   2.6406
 NONBonded  OH1     0.2342   2.6406      0.2342   2.6406
 NONBonded  OH      0.2342   2.6406      0.2342   2.6406 ! for HYP
 NONBonded  OHP     0.2342   2.6406      0.2342   2.6406
 NONBonded  S       0.0239   3.3854      0.0239   3.3854
 NONBonded  SM      0.0239   3.3854      0.0239   3.3854
 NONBonded  SH1E    0.0239   3.3854      0.0239   3.3854

!! HBOND AEXP=4 REXP=6 HAEX=4 AAEX=2   ACCEPTORS
!!     CUTHB=4.5 CTOFHB=4.0 CTONHB=3.5  CUTHA=90.0  CTOFHA=70.0  CTONHA=50.0
!
 AEXP 4
 REXP 6
 HAEX 4
 AAEX 2
!                   Emin      Rmin
!                (Kcal/mol)   (A)
 hbond N*+* N%      -3.00      3.0!  VALUES FROM VINOGRADOV AND LINELL FOR
 hbond N*+* O*      -3.50      2.9!  TYPICAL LENGTHS AND DEPTHS.
 hbond OH*  N%      -4.00      2.85
 hbond OH*  O*      -4.25      2.75
 hbond SH1E    N%      -3.00      3.0 !! added, ATB
 hbond SH1E    O*      -3.50      2.9 !! added, ATB



 elseif ($par_nonbonded eq "OPLSX") then
! these are close to the original OPLS parameters without introducing
! new atom types. The commented out lines are for atom types that 
! are unique in OPLS but would require additional atom types in 
! allhdg. Michael Nilges.
{ suggested values for refinement in H2O
  NBONds
    CUTNB=12   WMIN=1.5
    REPEl = 0.0     cdie shift     
    TOLErance = 0.5      NBXMOD = 5
    ctonnb=5.5 ctofnb=11.0 
  END
}
 evaluate ($repel_radius = 0.0)
 NONBonded  C       0.105   3.750       0.013    3.750  {OPLS C}!
 NONBonded  CCIS    0.105   3.750       0.013    3.750  {OPLS C}!
! NONBonded  C       0.110   3.750       0.014    3.750  {OPLS CAJ}!
 NONBonded  CY2     0.105   3.750       0.013    3.750  {OPLS C}
 NONBonded  CF      0.110   3.750       0.014    3.750  {OPLS CA}
 NONBonded  CY      0.110   3.750       0.014    3.750  {OPLS CA}
 NONBonded  CW      0.145   3.750       0.018    3.750  {OPLS CB}
! NONBonded  CW      0.130   3.750       0.016    3.750  {OPLS CN}
 NONBonded  CR1E    0.110   3.750       0.014    3.750  {OPLS CD}!
! NONBonded  CR1E    0.130   3.750       0.016    3.750  {OPLS CG}!
 NONBonded  CR1H    0.145   3.750       0.018    3.750  {OPLS CGJ}
 NONBonded  CR1W    0.110   3.750       0.014    3.750  {OPLS CD}!
 NONBonded  CH1E    0.080   3.800       0.010    3.800  {OPLS CH}
 NONBonded  CH1P    0.080   3.800       0.010    3.800  {OPLS CH} !KJR HYP
! NONBonded  CH1E    0.080   3.850       0.010    3.850  {OPLS CHJ}
 NONBonded  CRHH    0.145   3.750       0.018    3.750  {OPLS CP}
 NONBonded  C5W     0.145   3.750       0.018    3.750  {OPLS CS}
 NONBonded  CH2G    0.118   3.800       0.015    3.800  {OPLS C2} 
 NONBonded  CH2P    0.118   3.905       0.015    3.905  {OPLS C2J}
 NONBonded  CH2K    0.118   3.905       0.015    3.905  {OPLS C2J} !PGL
 NONBonded  CH2L    0.118   3.905       0.015    3.905  {OPLS C2J} !PGL
 NONBonded  CJ      0.105   3.750       0.013    3.905  {OPLS C2J} !PGL
 NONBonded  CQ      0.118   3.905       0.015    3.905  {OPLS C2J} !!!
 NONBonded  CL      0.118   3.905       0.015    3.905  {OPLS C2J} !!!
! NONBonded  CH2P    0.118   3.800       0.015    3.800  {OPLS C2} 
 NONBonded  CH2E    0.118   3.905       0.015    3.905  {OPLS C2J}
! NONBonded  CH3E    0.118   3.905       0.015    3.905  {OPLS C3}  
 NONBonded  CH3E    0.160   3.910       0.020    3.910  {OPLS C3J}
! NONBonded  CH3E    0.175   3.905       0.022    3.905  {OPLS C3K}
! NONBonded  CH3E    0.170   3.800       0.021    3.800  {OPLS C3L}
 NONBonded  C5      0.145   3.750       0.145    3.750  {OPLS CC}
 NONBonded  CRH     0.145   3.750       0.145    3.750  {OPLS 43}

 NONBonded  H       0.05     0.50        0.004     0.50
 NONBonded  HA      0.05     0.50        0.004     0.50
 NONBonded  HC      0.05     0.50        0.004     0.50

 NONBonded  N       0.170    3.250       0.021    3.250
 NONBonded  NR      0.170    3.250       0.021    3.250
 NONBonded  NK      0.170    3.250       0.021    3.250 !AZA/AZH (KJR)
 NONBonded  NL      0.170    3.250       0.021    3.250 !AZA/AZH (KJR)
 NONBonded  NM      0.170    3.250       0.021    3.250 !AZA/AZH (KJR)
 NONBonded  NH1     0.170    3.250       0.021    3.250
 NONBonded  NH2     0.170    3.250       0.021    3.250
 NONBonded  NH3     0.170    3.250       0.021    3.250
 NONBonded  NC2     0.170    3.250       0.021    3.250

 NONBonded  O       0.210    2.960       0.021    2.960
 NONBonded  OC      0.210    2.960       0.021    2.960
 NONBonded  O1P     0.210    2.960       0.021    2.960
 NONBonded  O3R     0.210    2.960       0.021    2.960
 NONBonded  OH1     0.170    3.070       0.021    3.070
 NONBonded  OHP     0.170    3.070       0.021    3.070
 NONBonded  OH      0.170    3.070       0.021    3.070 ! for HYP
 NONBonded  P       0.2000   3.738       0.085    3.738
 NONBonded  S       0.250    3.550       0.031    3.550
 NONBonded  SM      0.250    3.550       0.031    3.550
 NONBonded  SH1E    0.250    3.550       0.031    3.550

 NBFIX      H   H      14.013   1.67074      14.013   1.67074
 NBFIX      H   HA    277.278   4.07473     277.278   4.07473
 NBFIX      H   HC     14.013   1.67074      14.013   1.67074
 NBFIX      HA  HA   1846.41    5.76501    1846.41    5.76501  
 NBFIX      HA  HC    277.278   4.07473     277.278   4.07473
 NBFIX      HC  HC     14.013   1.67074      14.013   1.67074

 elseif ($par_nonbonded eq "PARAM19") then
! these are close to the original OPLS parameters without introducing
! new atom types. The commented out lines are for atom types that 
! are unique in OPLS but would require additional atom types in 
! allhdg. Michael Nilges.
{ suggested values for refinement in H2O
  NBONds
    CUTNB=12   WMIN=1.5
    REPEl = 0.0     cdie shift     
    TOLErance = 0.5      NBXMOD = 5
    ctonnb=5.5 ctofnb=11.0 
  END
}
 evaluate ($repel_radius = 0.89)
 evaluate ($repel_rcons = 25)
 evaluate ($repel_rexpo  = 4)
 evaluate ($repel_irexp  = 1)

 NONBonded  H       0.0498   1.4254      0.0498   1.4254
 NONBonded  HA      0.0498   1.4254      0.0450   2.6157 !- charged group.
 NONBonded  HC      0.0498   1.0691      0.0498   1.0691 !   Reduced vdw radius
 !
 NONBonded  C       0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CCIS    0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  C5      0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  C5W     0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CF      0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CW      0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CY      0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CY2     0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CH1P    0.0486   4.2140      0.1000   3.3854 ! KJR HYP
 NONBonded  CH1E    0.0486   4.2140      0.1000   3.3854 ! \
 NONBonded  CH2E    0.1142   3.9823      0.1000   3.3854 !  extended carbons
 NONBonded  CH2G    0.1142   3.9823      0.1000   3.3854 !  extended carbons
 NONBonded  CH2P    0.1142   3.9823      0.1000   3.3854 !  extended carbons
 NONBonded  CH3E    0.1811   3.8576      0.1000   3.3854 ! /
!! NONBonded  CM      0.0262   4.4367      0.1000   3.3854
 NONBonded  CR1E    0.1200   3.7418      0.1000   3.3854 !  ring carbons
 NONBonded  CR1H    0.1200   3.7418      0.1000   3.3854 !  ring carbons
 NONBonded  CR1W    0.1200   3.7418      0.1000   3.3854 !  ring carbons
 NONBonded  CRHH    0.1200   3.7418      0.1000   3.3854 !  ring carbons
 NONBonded  CRH     0.1200   3.7418      0.1000   3.3854 !  ring carbons
!! NONBonded  CT      0.0262   4.4367      0.1000   3.3854
 !
 NONBonded  N       0.2384   2.8509      0.2384   2.8509
 NONBonded  NC2     0.2384   2.8509      0.2384   2.8509
 NONBonded  NH1     0.2384   2.8509      0.2384   2.8509
 NONBonded  NH2     0.2384   2.8509      0.2384   2.8509
 NONBonded  NH3     0.2384   2.8509      0.2384   2.8509
 NONBonded  NP      0.2384   2.8509      0.2384   2.8509
 NONBonded  NR      0.2384   2.8509      0.2384   2.8509
 !
 NONBonded  O       0.1591   2.8509      0.1591   2.8509
 NONBonded  OC      0.6469   2.8509      0.6469   2.8509
 NONBonded  OH1     0.1591   2.8509      0.1591   2.8509
 NONBonded  OH      0.1591   2.8509      0.1591   2.8509 !KJR HYP
!! NONBonded  OM      0.1591   2.8509      0.1591   2.8509
 NONBonded  OS      0.1591   2.8509      0.1591   2.8509
 !
 NONBonded  S       0.0430   3.3676      0.0430   3.3676
 NONBonded  SM      0.0430   3.3676      0.0430   3.3676
 NONBonded  SH1E    0.0430   3.3676      0.0430   3.3676
 !
!! NONBONDED FE        0.0000    1.1582      0.0000 1.1582




 else {standard PARALLHDG parameters}
{ suggested values:
  NBONds
  CUTNB=7.0   WMIN=1.5
  REPEl = 0.78          
  REXPonent = 2
  IREXponent = 2
  RCONst = 5.0
  TOLErance = 0.5      NBXMOD = 5
  ctonnb=5.5 ctofnb=6.0 {* for consistency only, not needed for repel *}
END
} 
 evaluate ($repel_radius = 0.78)
 evaluate ($repel_rcons = 5.0)
 evaluate ($repel_rexpo  = 2)
 evaluate ($repel_irexp  = 2)
 NONBonded  C       0.0903   3.3409      0.0903   3.3409
 NONBonded  CCIS    0.0903   3.3409      0.0903   3.3409
 NONBonded  CR1E      0.1200   3.3409      0.1200   3.3409
 NONBonded  CF      0.1200   3.3409      0.1200   3.3409
 NONBonded  CY      0.1200   3.3409      0.1200   3.3409
 NONBonded  CY2      0.1200   3.3409      0.1200   3.3409
 NONBonded  CR1W      0.1200   3.3409      0.1200   3.3409
 NONBonded  CW      0.1450   3.3409      0.1450   3.3409
 NONBonded  C5W      0.1450   3.3409      0.1450   3.3409
 NONBonded  CN      0.1450   3.3409      0.1450   3.3409
 NONBonded  C5      0.1200   3.3409      0.1200   3.3409
 NONBonded  CH1E      0.0903   3.3409      0.0903   3.3409
 NONBonded  CH1P      0.0903   3.3409      0.0903   3.3409 !KJR HYP
 NONBonded  CH2E      0.0903   3.3409      0.0903   3.3409
 NONBonded  CH3E      0.0903   3.3409      0.0903   3.3409
 NONBonded  CH2G      0.0903   3.3409      0.0903   3.3409
 NONBonded  CH2P      0.1450   3.3409      0.1450   3.3409
 NONBonded  CH2K      0.0903   3.3409      0.0903   3.3409 !KJR PGL
 NONBonded  CH2L      0.0903   3.3409      0.0903   3.3409 !KJR PGL
 NONBonded  CJ        0.0903   3.3409      0.0903   3.3409 !KJR PGL
 NONBonded  CRH       0.1200   3.3409      0.1200   3.3409
 NONBonded  CR1H      0.1200   3.3409      0.1200   3.3409
 NONBonded  CRHH      0.1200   3.3409      0.1200   3.3409
 NONBonded  H       0.0498   2.2272      0.0498   2.2272
 NONBonded  HA      0.0045   2.2272      0.0045   2.2272
 NONBonded  HC      0.0498   2.2272      0.0498   2.2272
 NONBonded  N       0.1592   3.0068      0.1592   3.0068
 NONBonded  NR      0.1592   3.0068      0.1592   3.0068
 NONBonded  NH1      0.1592   3.0068      0.1592   3.0068
 NONBonded  NH2     0.1592   3.0068      0.1592   3.0068
 NONBonded  NH3     0.1592   3.0068      0.1592   3.0068
 NONBonded  NC2     0.1592   3.0068      0.1592   3.0068
 NONBonded  O       0.2342   2.7755      0.2342   2.7755
 NONBonded  OC      1.0244   2.7755      1.0244   2.7755
 NONBonded  OH1      0.2342   2.7755      0.2342   2.7755
 NONBonded  OH       0.2342   2.7755      0.2342   2.7755 !KJR HYP
 NONBonded  S        0.0239   3.7458      0.0239   3.7458
 NONBonded  SM       0.0239   3.7458      0.0239   3.7458
 NONBonded  SH1E       0.0239   3.7458      0.0239   3.7458

 NBFIx  H    NR         44.200        1.000          44.200        1.000
 NBFIx  H    O          44.200        1.000          44.200        1.000
 NBFIx  H    OC         44.200        1.000          44.200        1.000
 NBFIx  H    OH1         44.200        1.000          44.200        1.000
 NBFIx  HC   NR         44.200        1.000          44.200        1.000
 NBFIx  HC   O          44.200        1.000          44.200        1.000
 NBFIx  HC   OC         44.200        1.000          44.200        1.000
 NBFIx  HC   OH1         44.200        1.000          44.200        1.000

 end if

!set echo on message on end

