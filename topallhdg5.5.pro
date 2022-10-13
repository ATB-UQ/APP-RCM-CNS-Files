remark   file topallhdg.pro version 5.5 date 10-July-2012 
remark   file topallhdg.pro  version 5.4  date 17-March-06
remark      for file parallhdg.pro version 5.3 date 13-Feb-02 or later
remark   Geometric energy function parameters for distance geometry and
remark      simulated annealing.
remark   Author: Michael Nilges, EMBL Heidelberg; Institut Pasteur, Paris
remark   This file contains modifications from M. Williams, UCL London
remark   Latest modifications by Alexandre Bonvin: 
remark     adding partial charges to aromatic protons
remark     and carbons for TYR,PHE,HIS and TRP (taken from OPLSAA)

set echo off message off end

!***********************************************************************!
! Copyright (C) 1995,1996 by Michael Nilges. All rights reserved.       !
! Copying and redistribution of this files is authorized only if etiher !
! (1) you make absolutely no changes to your copy, including name, or   !
! (2) if you do make changes, you name it something other than          !
! topallhdg.pro and topallhdg.x.xx.pro, and clearly mark the changes.   ! 
! The information in this software is subject to change without notice  !
! and should not be construed as a commitment by the EMBL or by the     !
! authors. Neither the EMBL, Yale University, nor the authors assume    !
! responsibility for the use or reliability of this software.           !
! We do hope, however, to get responses from users, especially when     !
! errors have been found.                                               !
!***********************************************************************!
! Description:                                                          !
! This parameter file was originally derived from the  CHARMM parameter !
! file PARMALLH6. It was designed specifically for the initial stages   ! 
! of calculating structures from NMR restraints.                           !
!***********************************************************************!
! History:
! version 5.5  (12-July-12) :Original file topallhdg5.4.pro taken from the topppar directory of Haddock2.0
! Modified by Johan Rosengren to include: Hydroxy proline, Pyroglutamic acid, Abu, triazole-bridge, C-terminal amidation.
! version 5.3  (13-Feb-02) : merged UCL version with EMBL version
! version UCL  (07-JUL-01) : added new MULT disulphide/sidechain COO dihedral parameters
! version UCL  (03-JUL-01) : added new MULT phi/psi dihedral related parameters
! version UCL  (14-MAR-00) : modified peptide bond parameters for flexible omega
! version UCL  (12-MAR-00) : added C and N terminus and disulphide dihedrals
! version UCL  (09-MAR-00) : TYR O-H planarity & LYS NH3 stagger 
!                            & HIS dihedrals modified
! version 5.2  (18-Jun-99) : new cis peptide patch
! version 4.03 (04-APR-99) : THR CB atom type corrected
! version 4.02 (17-DEC-97) : Histidine hbond acceptor/donor corrected
! version 4.01 (29-Jul-96) : all covalent parameters
! version 4.00  (19-Jul-96) : all atom types from CSDX implemented
! version 3.00 (24-Oct-95) : mapped CSDX parmameters on parallhdg, 
!                            no changes in topallhdg
! previous modifications:
! proline residue modified, puckering enforced (MN)
! added hbond acceptor and donor definitions for analysis (MN)
! all references to internal coordinates (IC's) removed (MN) 
! added stereospecific impropers for all pro-chiral centers (ATB, JK)
! all dihedrals defining planarity converted to impropers (MN, PK)
! additional impropers at planar centers (MN)
!***********************************************************************!


set message off echo off end

autogenerate 
  angles=true
  dihedrals=false
end

mass H    1.008
mass HC   1.008
mass HA   1.008
mass C   12.011
mass CCIS  12.011
mass CH1E  12.011
mass CH2E  12.011
mass CH3E  12.011
mass CH2G  12.011
mass CH2P  12.011
mass CH1P  12.011 !HYP added by JR
mass C5W  12.011
mass CW  12.011

mass CR1E  12.011
mass C5  12.011
mass CRH  12.011
mass CR1H  12.011
mass CR1W  12.011
mass CRHH  12.011
mass CJ  12.011 !PGL added by JR
mass CH2K  12.011 !PGL added by JR
mass CH2L  12.011 !PGL added by JR
mass CQ  12.011 !PRG added by JR
mass CL  12.011 !PRG added by JR
mass CF  12.011
mass CY  12.011
mass CY2 12.011
mass N   14.007
mass NR  14.007
mass NK  14.007 !AZA/AZH added by JR
mass NL  14.007 !AZA/AZH added by JR
mass NM  14.007
mass NH1 14.007
mass NH2 14.007
mass NH3 14.007
mass NC2 14.007
mass O   15.999
mass OH  15.999
mass O1P 15.999
mass O3R 15.999
mass OC  15.999
mass OH1 15.999
mass P   30.97400
mass SH1E  32.060
mass SM  32.060
mass S   32.060

    RESIdue NTER
      ATOM N    TYPE=NH1  CHARge=-0.570   MASS=14.0067  END
      ATOM CA   TYPE=CH1E CHARge=0.200    MASS=12.0110  END
      ATOM CB   TYPE=CH2E CHARge=-0.349   MASS=12.0110  END
      ATOM CG   TYPE=CG3  CHARge=0.824    MASS=12.0110  END
      ATOM ND2  TYPE=ND2  CHARge=-1.116   MASS=14.0067  END
      ATOM OD1  TYPE=OD1  CHARge=-0.611   MASS=15.9994  END
      ATOM H2   TYPE=H2   CHARge=0.346    MASS=1.0080   END
      ATOM H3   TYPE=H3   CHARge=0.346    MASS=1.0080   END
      ATOM H1   TYPE=H1   CHARge=0.346    MASS=1.0080   END
      ATOM HA   TYPE=HA   CHARge=0.000    MASS=1.0080   END
      ATOM HB3  TYPE=HB31 CHARge=0.148    MASS=1.0080   END
      ATOM HB2  TYPE=HB24 CHARge=0.146    MASS=1.0080   END
      ATOM HD21 TYPE=HD211 CHARge=0.460    MASS=1.0080   END
      ATOM HD22 TYPE=HD22 CHARge=0.460    MASS=1.0080   END
      ATOM C    TYPE=C    CHARge=0.500    MASS=12.0110  END
      ATOM O    TYPE=O    CHARge=-0.500   MASS=15.9994  END

      BOND C    O   
      BOND CA   C   
      BOND CA   HA  
      BOND N    CA  
      BOND CA   CB  
      BOND N    H2  
      BOND N    H3  
      BOND N    H1  
      BOND CB   HB3 
      BOND CB   HB2 
      BOND CB   CG  
      BOND CG   OD1 
      BOND CG   ND2 
      BOND ND2  HD21
      BOND ND2  HD22

      ANGLe CA   C    O   
      ANGLe HA   CA   C   
      ANGLe N    CA   C   
      ANGLe CB   CA   C   
      ANGLe N    CA   HA  
      ANGLe CB   CA   HA  
      ANGLe N    CA   CB  
      ANGLe CA   N    H2  
      ANGLe CA   N    H3  
      ANGLe CA   N    H1  
      ANGLe H2   N    H3  
      ANGLe H2   N    H1  
      ANGLe H3   N    H1  
      ANGLe CA   CB   HB3 
      ANGLe CA   CB   HB2 
      ANGLe CA   CB   CG  
      ANGLe HB3  CB   HB2 
      ANGLe CG   CB   HB3 
      ANGLe CG   CB   HB2 
      ANGLe CB   CG   OD1 
      ANGLe CB   CG   ND2 
      ANGLe ND2  CG   OD1 
      ANGLe CG   ND2  HD21
      ANGLe CG   ND2  HD22
      ANGLe HD21 ND2  HD22

      DIHEdral H2   N    CA   C   
      DIHEdral CA   CB   CG   ND2 
      DIHEdral N    CA   CB   CG  
      DIHEdral CB   CG   ND2  HD21

      IMPRoper CG   CB   OD1  ND2 
      IMPRoper ND2  CG   HD21 HD22
      IMPRoper HA   N    C    CB  
   END

  RESIdue CTER
      ATOM N    TYPE=NH1  CHARge=-0.570   MASS=14.0067  END
      ATOM CA   TYPE=CH1E CHARge=0.200    MASS=12.0110  END
      ATOM C    TYPE=C    CHARge=0.500    MASS=12.0110  END
      ATOM O    TYPE=O    CHARge=-0.500   MASS=15.9994  END
      ATOM CB   TYPE=CH2E CHARge=-0.108   MASS=12.0110  END
      ATOM CG   TYPE=CG2  CHARge=-0.086   MASS=12.0110  END
      ATOM CD   TYPE=CD3  CHARge=0.431    MASS=12.0110  END
      ATOM OE1  TYPE=OE1  CHARge=-0.786   MASS=15.9994  END
      ATOM OE2  TYPE=OE2  CHARge=-0.786   MASS=15.9994  END
      ATOM HA   TYPE=HA   CHARge=0.000    MASS=1.0080   END
      ATOM O2   TYPE=O2   CHARge=-0.729   MASS=15.9994  END
      ATOM HN   TYPE=H    CHARge=0.370    MASS=1.0080   END
      ATOM HB3  TYPE=HB3  CHARge=0.028    MASS=1.0080   END
      ATOM HB2  TYPE=HB2  CHARge=0.034    MASS=1.0080   END
      ATOM HG3  TYPE=HG32 CHARge=0.003    MASS=1.0080   END
      ATOM HG2  TYPE=HG22 CHARge=-0.001   MASS=1.0080   END

      BOND CD   OE2 
      BOND CD   OE1 
      BOND CG   CD  
      BOND CG   HG3 
      BOND CG   HG2 
      BOND CB   CG  
      BOND CB   HB3 
      BOND CB   HB2 
      BOND CA   CB  
      BOND CA   HA  
      BOND CA   C   
      BOND N    CA  
      BOND C    O   
      BOND C    O2  
      BOND N    HN  

      ANGLe OE1  CD   OE2 
      ANGLe CG   CD   OE2 
      ANGLe CG   CD   OE1 
      ANGLe CD   CG   HG3 
      ANGLe CD   CG   HG2 
      ANGLe CB   CG   CD  
      ANGLe HG3  CG   HG2 
      ANGLe CB   CG   HG3 
      ANGLe CB   CG   HG2 
      ANGLe CG   CB   HB3 
      ANGLe CG   CB   HB2 
      ANGLe CA   CB   CG  
      ANGLe HB3  CB   HB2 
      ANGLe CA   CB   HB3 
      ANGLe CA   CB   HB2 
      ANGLe CB   CA   HA  
      ANGLe C    CA   CB  
      ANGLe N    CA   CB  
      ANGLe C    CA   HA  
      ANGLe N    CA   HA  
      ANGLe N    CA   C   
      ANGLe CA   C    O   
      ANGLe CA   C    O2  
      ANGLe O    C    O2  
      ANGLe CA   N    HN  

      DIHEdral CA   CB   CG   CD  
      DIHEdral CB   CG   CD   OE1 
      DIHEdral N    CA   CB   CG  

      IMPRoper CD   OE2  OE1  CG  
      IMPRoper HA   N    C    CB  
   END

 RESIdue NMA
      ATOM N    TYPE=NH1  CHARge=-0.570   MASS=14.0067  END
      ATOM HN   TYPE=H    CHARge=0.370    MASS=1.0080   END
      ATOM CA   TYPE=CH1E CHARge=0.200    MASS=12.0110  END
      ATOM CB   TYPE=CH2E CHARge=-0.172   MASS=12.0110  END
      ATOM CG   TYPE=CG   CHARge=-0.175   MASS=12.0110  END
      ATOM CD   TYPE=CD1  CHARge=-0.049   MASS=12.0110  END
      ATOM NE   TYPE=NE   CHARge=-0.278   MASS=14.0067  END
      ATOM CZ   TYPE=CZ   CHARge=0.151    MASS=12.0110  END
      ATOM NH2  TYPE=NH21 CHARge=-0.297   MASS=14.0067  END
      ATOM NH1  TYPE=NH11 CHARge=-0.362   MASS=14.0067  END
      ATOM CQ2  TYPE=CQ2  CHARge=-0.124   MASS=12.0110  END
      ATOM HA   TYPE=HA   CHARge=0.000    MASS=1.0080   END
      ATOM HB3  TYPE=HB3  CHARge=0.107    MASS=1.0080   END
      ATOM HB2  TYPE=HB2  CHARge=0.105    MASS=1.0080   END
      ATOM HG3  TYPE=HG3  CHARge=0.082    MASS=1.0080   END
      ATOM HG2  TYPE=HG2  CHARge=0.101    MASS=1.0080   END
      ATOM HD3  TYPE=HD3  CHARge=0.111    MASS=1.0080   END
      ATOM HD2  TYPE=HD2  CHARge=0.080    MASS=1.0080   END
      ATOM HE   TYPE=HE   CHARge=0.135    MASS=1.0080   END
      ATOM HH11 TYPE=HH11 CHARge=0.186    MASS=1.0080   END
      ATOM HH12 TYPE=HH12 CHARge=0.186    MASS=1.0080   END
      ATOM HQ21 TYPE=HQ21 CHARge=0.071    MASS=1.0080   END
      ATOM HQ22 TYPE=HQ22 CHARge=0.071    MASS=1.0080   END
      ATOM HQ23 TYPE=HQ23 CHARge=0.071    MASS=1.0080   END
      ATOM C    TYPE=C    CHARge=0.500    MASS=12.0110  END
      ATOM O    TYPE=O    CHARge=-0.500   MASS=15.9994  END

      BOND C    O   
      BOND CA   C   
      BOND CA   HA  
      BOND N    CA  
      BOND CA   CB  
      BOND N    HN  
      BOND CB   HB3 
      BOND CB   HB2 
      BOND CB   CG  
      BOND CG   HG3 
      BOND CG   HG2 
      BOND CG   CD  
      BOND CD   HD3 
      BOND CD   HD2 
      BOND CD   NE  
      BOND NE   HE  
      BOND NE   CZ  
      BOND CZ   NH1 
      BOND CZ   NH2 
      BOND NH1  HH11
      BOND NH1  HH12
      BOND NH2  CQ2 
      BOND CQ2  HQ21
      BOND CQ2  HQ22
      BOND CQ2  HQ23

      ANGLe CA   C    O   
      ANGLe HA   CA   C   
      ANGLe N    CA   C   
      ANGLe CB   CA   C   
      ANGLe N    CA   HA  
      ANGLe CB   CA   HA  
      ANGLe N    CA   CB  
      ANGLe HN   N    CA  
      ANGLe CA   CB   HB3 
      ANGLe CA   CB   HB2 
      ANGLe CA   CB   CG  
      ANGLe HB3  CB   HB2 
      ANGLe CG   CB   HB3 
      ANGLe CG   CB   HB2 
      ANGLe CB   CG   HG3 
      ANGLe CB   CG   HG2 
      ANGLe CB   CG   CD  
      ANGLe HG3  CG   HG2 
      ANGLe CD   CG   HG3 
      ANGLe CD   CG   HG2 
      ANGLe CG   CD   HD3 
      ANGLe CG   CD   HD2 
      ANGLe CG   CD   NE  
      ANGLe HD3  CD   HD2 
      ANGLe NE   CD   HD3 
      ANGLe NE   CD   HD2 
      ANGLe CD   NE   HE  
      ANGLe CD   NE   CZ  
      ANGLe CZ   NE   HE  
      ANGLe NE   CZ   NH1 
      ANGLe NE   CZ   NH2 
      ANGLe NH2  CZ   NH1 
      ANGLe CZ   NH1  HH11
      ANGLe CZ   NH1  HH12
      ANGLe HH11 NH1  HH12
      ANGLe CZ   NH2  CQ2 
      ANGLe NH2  CQ2  HQ21
      ANGLe NH2  CQ2  HQ22
      ANGLe NH2  CQ2  HQ23
      ANGLe HQ21 CQ2  HQ22
      ANGLe HQ21 CQ2  HQ23
      ANGLe HQ22 CQ2  HQ23

      !DIHEdral N    CA   C    O   
      DIHEdral CA   CB   CG   CD  
      DIHEdral N    CA   CB   CG  
      DIHEdral CB   CG   CD   NE  
      DIHEdral CG   CD   NE   CZ  
      DIHEdral CD   NE   CZ   NH2 
      DIHEdral NE   CZ   NH1  HH11
      DIHEdral NE   CZ   NH2  CQ2 
      DIHEdral CZ   NH2  CQ2  HQ21

      IMPRoper CZ   NE   NH1  NH2 
      IMPRoper HA   N    C    CB  
   END

  RESIdue PAL
      ATOM N    TYPE=NH1  CHARge=-0.570   MASS=14.0067  END
      ATOM HN   TYPE=H    CHARge=0.370    MASS=1.0080   END
      ATOM CA   TYPE=CH1E CHARge=0.200    MASS=12.0110  END
      ATOM CB1  TYPE=CB1  CHARge=-0.413   MASS=12.0110  END
      ATOM CG   TYPE=CG   CHARge=0.399    MASS=12.0110  END
      ATOM CD   TYPE=CD1  CHARge=-0.045   MASS=12.0110  END
      ATOM CE   TYPE=CE   CHARge=-0.329   MASS=12.0110  END
      ATOM CZ   TYPE=CZ   CHARge=-0.022   MASS=12.0110  END
      ATOM CH   TYPE=CH   CHARge=-0.270   MASS=12.0110  END
      ATOM CB2  TYPE=CB2  CHARge=-0.375   MASS=12.0110  END !edited by JR to make net charge 0
      ATOM HB13 TYPE=HB13 CHARge=0.1065   MASS=1.0080   END !averaged
      ATOM HB12 TYPE=HB12 CHARge=0.1065   MASS=1.0080   END !averaged
      ATOM HG2  TYPE=HG2  CHARge=-0.040   MASS=1.0080   END !averaged
      ATOM HD2  TYPE=HD2  CHARge=0.045    MASS=1.0080   END !averaged
      ATOM HG3  TYPE=HG3  CHARge=-0.040   MASS=1.0080   END !averaged
      ATOM HE   TYPE=HE   CHARge=0.165    MASS=1.0080   END
      ATOM HD3  TYPE=HD3  CHARge=0.045    MASS=1.0080   END !averaged
      ATOM HZ   TYPE=HZ   CHARge=0.102    MASS=1.0080   END
      ATOM HH1  TYPE=HH1  CHARge=0.088    MASS=1.0080   END
      ATOM HH2  TYPE=HH2  CHARge=0.088    MASS=1.0080   END
      ATOM HH3  TYPE=HH3  CHARge=0.088    MASS=1.0080   END
      ATOM HB21 TYPE=HB21 CHARge=0.125    MASS=1.0080   END
      ATOM HB22 TYPE=HB22 CHARge=0.125    MASS=1.0080   END
      ATOM HB23 TYPE=HB23 CHARge=0.125    MASS=1.0080   END
      ATOM C    TYPE=C    CHARge=0.500    MASS=12.0110  END
      ATOM O    TYPE=O    CHARge=-0.500   MASS=15.9994  END

      BOND C    O   
      BOND CA   C   
      BOND CA   CB2 
      BOND CA   CB1 
      BOND N    CA  
      BOND CB2  HB21
      BOND CB2  HB22
      BOND CB2  HB23
      BOND CB1  HB13
      BOND CB1  HB12
      BOND CB1  CG  
      BOND CG   HG2 
      BOND CG   HG3 
      BOND CG   CD  
      BOND CD   HD2 
      BOND CD   HD3 
      BOND CD   CE  
      BOND CE   HE  
      BOND CE   CZ  
      BOND CZ   HZ  
      BOND CZ   CH  
      BOND CH   HH1 
      BOND CH   HH2 
      BOND CH   HH3 
      BOND N    HN  

      ANGLe CA   C    O   
      ANGLe CB2  CA   C   
      ANGLe CB1  CA   C   
      ANGLe N    CA   C   
      ANGLe CB1  CA   CB2 
      ANGLe N    CA   CB2 
      ANGLe N    CA   CB1 
      ANGLe CA   CB2  HB21
      ANGLe CA   CB2  HB22
      ANGLe CA   CB2  HB23
      ANGLe HB21 CB2  HB22
      ANGLe HB21 CB2  HB23
      ANGLe HB22 CB2  HB23
      ANGLe CA   CB1  HB13
      ANGLe CA   CB1  HB12
      ANGLe CA   CB1  CG  
      ANGLe HB13 CB1  HB12
      ANGLe CG   CB1  HB13
      ANGLe CG   CB1  HB12
      ANGLe CB1  CG   HG2 
      ANGLe CB1  CG   HG3 
      ANGLe CB1  CG   CD  
      ANGLe HG2  CG   HG3 
      ANGLe CD   CG   HG2 
      ANGLe CD   CG   HG3 
      ANGLe CG   CD   HD2 
      ANGLe CG   CD   HD3 
      ANGLe CG   CD   CE  
      ANGLe HD2  CD   HD3 
      ANGLe CE   CD   HD2 
      ANGLe CE   CD   HD3 
      ANGLe CD   CE   HE  
      ANGLe CD   CE   CZ  
      ANGLe CZ   CE   HE  
      ANGLe CE   CZ   HZ  
      ANGLe CE   CZ   CH  
      ANGLe CH   CZ   HZ  
      ANGLe CZ   CH   HH1 
      ANGLe CZ   CH   HH2 
      ANGLe CZ   CH   HH3 
      ANGLe HH1  CH   HH2 
      ANGLe HH1  CH   HH3 
      ANGLe HH2  CH   HH3 
      ANGLe HN   N    CA  

      DIHEdral CA   CB1  CG   CD  
      DIHEdral CB1  CG   CD   CE  
      DIHEdral CG   CD   CE   CZ  
      DIHEdral CD   CE   CZ   CH  
      DIHEdral CE   CZ   CH   HH1 
      DIHEdral N    CA   CB2  HB21
      DIHEdral N    CA   CB1  CG  

      IMPRoper CE   CD   HE   CZ  
      IMPRoper CZ   CE   HZ   CH  
      IMPRoper CB2  N    C    CB1  
   END

residue ALA
  group
    atom N   type=NH1     charge=-0.570 end
    atom HN  type=H       charge= 0.370 end
    atom CA  type=CH1E    charge= 0.200 end
    atom HA  type=HA      charge= 0.000 end
    atom CB  type=CH3E    charge= 0.000 end
    atom HB1 type=HA      charge= 0.000 end
    atom HB2 type=HA      charge= 0.000 end
    atom HB3 type=HA      charge= 0.000 end
    atom C   type=C       charge= 0.500 end
    atom O   type=O       charge=-0.500 end

  bond N  HN       
  bond N  CA    bond CA  HA   
  bond CA  CB   bond CB  HB1   bond CB  HB2     bond CB  HB3
  bond CA  C       
  bond C   O

  improper HA  N   C   CB  ! chirality CA
  improper HB1 HB2 CA HB3  ! methyl CB

  dihedral HB1 CB  CA  C   ! methyl stagger UCL 12-MAR-00 

  DONO HN   N
  ACCE O    C
end
 

  

residue ARG
  group
    atom N    type=NH1    charge=-0.570 end
    atom HN   type=H      charge= 0.370 end
    atom CA   type=CH1E   charge= 0.200 end
    atom HA   type=HA     charge= 0.000 end
    atom CB   type=CH2E   charge= 0.000 end
    atom HB1  type=HA     charge= 0.000 end
    atom HB2  type=HA     charge= 0.000 end
    atom CG   type=CH2E   charge= 0.070 end
    atom HG1  type=HA     charge= 0.000 end
    atom HG2  type=HA     charge= 0.000 end
    atom CD   type=CH2E   charge= 0.310 end
    atom HD1  type=HA     charge= 0.000 end
    atom HD2  type=HA     charge= 0.000 end
    atom NE   type=NH1    charge=-0.700 end
    atom HE   type=H      charge= 0.440 end
    atom CZ   type=C      charge= 0.640 end
    atom NH1  type=NC2    charge=-0.800 end
    atom HH11 type=HC     charge= 0.460 end
    atom HH12 type=HC     charge= 0.460 end
    atom NH2   type=NC2   charge=-0.800 end
    atom HH21  type=HC    charge= 0.460 end
    atom HH22  type=HC    charge= 0.460 end
    atom C     type=C     charge= 0.500 end
    atom O     type=O     charge=-0.500 end
 
  bond N  HN
  bond N  CA     bond CA  HA
  bond CA CB     bond CB  HB1     bond CB  HB2
  bond CB CG     bond CG  HG1     bond CG  HG2
  bond CG CD     bond CD  HD1     bond CD  HD2
  bond CD NE     bond NE  HE
  bond NE CZ
  bond CZ NH1    bond NH1 HH11    bond NH1 HH12
  bond CZ NH2    bond NH2 HH21    bond NH2 HH22
  bond CA C
  bond C  O

  improper HA  N  C    CB  !chirality CA
  improper NE  CD CZ   HE  
  improper CZ  NE NH1  NH2
  improper NH1 CZ HH11 HH12
  improper NH2 CZ HH21 HH22
  improper NE  CZ NH1  HH11
  improper NE  CZ NH2  HH21
  improper CZ  NH2 HE NE  ! planar HE, CZ
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG
  improper HD1 HD2 CG NE  !stereo CD
  
  dihedral CG  CB  CA  N
  dihedral CD  CG  CB  CA
  dihedral NE  CD  CG  CB
  dihedral CZ  NE  CD  CG
  
  DONO HN   N
  DONO HE   NE
  DONO HH11 NH1
  DONO HH12 NH1
  DONO HH21 NH2
  DONO HH22 NH2
  ACCE O    C
end

residue ASN
  group
    atom N    type=NH1    charge=-0.570 end
    atom HN   type=H      charge= 0.370 end
    atom CA   type=CH1E   charge= 0.200 end
    atom HA   type=HA     charge= 0.000 end
    atom CB   type=CH2E   charge=-0.000 end
    atom HB1  type=HA     charge= 0.000 end
    atom HB2  type=HA     charge= 0.000 end
    atom CG   type=C      charge= 0.500 end
    atom OD1  type=O      charge=-0.500 end
    atom ND2  type=NH2    charge=-0.850 end
    atom HD21 type=H      charge= 0.425 end
    atom HD22 type=H      charge= 0.425 end
    atom C    type=C      charge= 0.500 end
    atom O    type=O      charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA  HA
  bond CA CB     bond CB  HB1     bond CB  HB2
  bond CB CG
  bond CG OD1
  bond CG ND2    bond ND2 HD21    bond ND2 HD22
  bond CA C
  bond C  O

  improper HA  N  C    CB  ! chirality CA
  improper CG  CB OD1  ND2
  improper ND2 CG HD21 HD22
  improper CB  CG ND2  HD21
  improper HB1 HB2 CA CG  ! stereo CB

  dihedral CG  CB  CA  N
  dihedral OD1 CG  CB  CA ! mult 2 ! UCL

  DONO HN   N
  DONO HD21 ND2
  DONO HD22 ND2
  ACCE OD1  CG
  ACCE O    C
end


residue ASP
  group
    atom N   type=NH1   charge=-0.570 end
    atom HN  type=H     charge= 0.370 end
    atom CA  type=CH1E  charge= 0.200 end
    atom HA  type=HA    charge= 0.000 end
    atom CB  type=CH2E  charge=-0.100 end
    atom HB1 type=HA    charge= 0.000 end
    atom HB2 type=HA    charge= 0.000 end
    atom CG  type=C     charge= 0.700 end
    atom OD1 type=OC    charge=-0.800 end
    atom OD2 type=OC    charge=-0.800 end
    atom C   type=C     charge= 0.500 end
    atom O   type=O     charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG
  bond CG OD1
  bond CG OD2
  bond CA C
  bond C  O

  improper HA N  C   CB   !chirality CA
  improper CG CB OD1 OD2
  improper HB1 HB2 CA CG  !stereo CB

  dihedral CG  CB  CA  N

  dihedral OD1 CG  CB  CA

  DONO HN   N
  ACCE OD1  CG
  ACCE OD2  CG
  ACCE O    C
end


residue ASPH
  group
    atom N   type=NH1   charge=-0.570 end
    atom HN  type=H     charge= 0.370 end
    atom CA  type=CH1E  charge= 0.200 end
    atom HA  type=HA    charge= 0.000 end
    atom CB  type=CH2E  charge= 0.000 end
    atom HB1 type=HA    charge= 0.000 end
    atom HB2 type=HA    charge= 0.000 end
    atom CG  type=C     charge= 0.500 end
    atom OD1 type=OC    charge=-0.500 end
    atom OD2 type=OH1   charge=-0.450 end
    atom HD2 type=H     charge= 0.450 end
    atom C   type=C     charge= 0.500 end
    atom O   type=O     charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG
  bond CG OD1
  bond CG OD2    bond OD2 HD2
  bond CA C
  bond C  O

  improper HA N  C   CB   !chirality CA
  improper CG CB OD1 OD2
  improper HB1 HB2 CA CG  !stereo CB

  dihedral CG  CB  CA  N

  dihedral OD1 CG  CB  CA

  DONO HN   N
  DONO HD2  OD2
  ACCE OD1  CG
  ACCE OD2  CG
  ACCE O    C
end

residue AZA !azidoalanine - charges guessed KJR.
  group
    atom N   type=NH1     charge=-0.570 end
    atom HN  type=H       charge= 0.370 end
    atom CA  type=CH1E    charge= 0.200 end
    atom HA  type=HA      charge= 0.000 end
    atom CB  type=CH2E    charge= 0.180 end
    atom HB1 type=HA      charge= 0.000 end
    atom HB2 type=HA      charge= 0.000 end
    atom NG  type=NK      charge=-0.180 end
    atom ND  type=NL      charge= 0.500 end
    atom NE  type=NM      charge=-0.500 end
    atom C   type=C       charge= 0.500 end
    atom O   type=O       charge=-0.500 end

  bond N  HN       
  bond N  CA    bond CA  HA   
  bond CA  CB   bond CB  HB1   bond CB  HB2     bond CB  NG
  bond NG  ND   bond ND  NE
  bond CA  C       
  bond C   O

  improper HA  N   C   CB  ! chirality CA
  improper HB1 HB2 CA  ND  ! methyl CB
  
  dihedral ND  CB  CA  N     

  DONO HN   N
  ACCE O    C
  ACCE NE   ND
  
end

residue AZH !azidohomoalanine - charges guessed KJR.
  group
    atom N   type=NH1     charge=-0.570 end
    atom HN  type=H       charge= 0.370 end
    atom CA  type=CH1E    charge= 0.200 end
    atom HA  type=HA      charge= 0.000 end
    atom CB  type=CH2E    charge= 0.000 end
    atom HB1 type=HA      charge= 0.000 end
    atom HB2 type=HA      charge= 0.000 end
    atom CG  type=CH2E    charge= 0.180 end
    atom HG1 type=HA      charge= 0.000 end
    atom HG2 type=HA      charge= 0.000 end
    atom ND  type=NK      charge=-0.180 end
    atom NE  type=NL      charge= 0.500 end
    atom NZ  type=NM      charge=-0.500 end
    atom C   type=C       charge= 0.500 end
    atom O   type=O       charge=-0.500 end

  bond N  HN       
  bond N  CA    bond CA  HA   
  bond CA  CB   bond CB  HB1   bond CB  HB2     bond CB  CG
  bond CG  ND   bond ND  NE    bond NE  NZ
  bond CA  C       
  bond C   O

  improper HA  N   C   CB  ! chirality CA
  improper HB1 HB2 CA  CG  ! methyl CB
  improper HG1 HG2 CB  ND  ! methyl CG
  
  dihedral CG  CB  CA  N     
  dihedral ND  CG  CB  CA

  DONO HN   N
  ACCE O    C
  ACCE NZ   NE
  
end

residue CSP
  group
    atom N   type=NH1    charge=-0.570 end
    atom HN  type=H      charge= 0.370 end
    atom CA  type=CH1E   charge= 0.200 end
    atom HA  type=HA     charge= 0.000 end
    atom CB  type=CH2E   charge= 0.080 end
    atom HB1 type=HA     charge= 0.000 end
    atom HB2 type=HA     charge= 0.000 end
    atom SG  type=SM     charge=-0.440 end
    atom P   type=P      charge= 0.520 end
    atom O1P type=O1P    charge=-0.720 end
    atom O2P type=O1P    charge=-0.720 end
    atom O3P type=O3R    charge=-0.720 end
    atom C   type=C      charge= 0.500 end
    atom O   type=O      charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB SG     bond SG P 
  bond P   O1P
  bond P   O2P
  bond P   O3P
  bond CA C
  bond C  O

  improper HA N C CB   !chirality CA
  improper HB1 HB2 CA SG  !stereo CB
  improper SG  P O3P O1P ! chirality or flatness improper   -37.36

  dihedral SG  CB  CA  N
  dihedral CB  SG  P   O3P

  DONO HN   N
  ACCE O    C
  ACCE O1P  P
  ACCE O2P  P
  ACCE O3P  P

end


residue CYS
  group
    atom N   type=NH1    charge=-0.570 end
    atom HN  type=H      charge= 0.370 end
    atom CA  type=CH1E   charge= 0.200 end
    atom HA  type=HA     charge= 0.000 end
    atom CB  type=CH2E   charge= 0.180 end
    atom HB1 type=HA     charge= 0.000 end
    atom HB2 type=HA     charge= 0.000 end
    atom SG  type=SH1E   charge=-0.450 end
    atom HG  type=H      charge= 0.270 end
    atom C   type=C      charge= 0.500 end
    atom O   type=O      charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB SG     bond SG HG 
  bond CA C
  bond C  O

  improper HA N C CB   !chirality CA
  improper HB1 HB2 CA SG  !stereo CB

  dihedral SG  CB  CA  N

  DONO HN   N
  DONO HG  SG
  !!! ACCE SG    !REMOVED, ATB
  ACCE O    C
end


residue GLN
  group
    atom N    type=NH1     charge=-0.570 end
    atom HN   type=H       charge= 0.370 end
    atom CA   type=CH1E    charge= 0.200 end
    atom HA   type=HA      charge= 0.000 end
    atom CB   type=CH2E    charge= 0.000 end
    atom HB1  type=HA      charge= 0.000 end
    atom HB2  type=HA      charge= 0.000 end
    atom CG   type=CH2E    charge= 0.000 end
    atom HG1  type=HA      charge= 0.000 end
    atom HG2  type=HA      charge= 0.000 end
    atom CD   type=C       charge= 0.500 end
    atom OE1  type=O       charge=-0.500 end
    atom NE2  type=NH2     charge=-0.850 end
    atom HE21 type=H       charge= 0.425 end
    atom HE22 type=H       charge= 0.425 end
    atom C    type=C       charge= 0.500 end
    atom O    type=O       charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA  HA
  bond CA CB     bond CB  HB1     bond CB  HB2
  bond CB CG     bond CG  HG1     bond CG  HG2
  bond CG CD
  bond CD OE1
  bond CD NE2    bond NE2 HE21    bond NE2 HE22
  bond CA C
  bond C  O

  improper HA  N  C    CB  !chirality CA
  improper CD  CG OE1  NE2
  improper NE2 CD HE21 HE22
  improper CG  CD NE2  HE21
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG

  dihedral CG  CB  CA  N
  dihedral CD  CG  CB  CA
  dihedral OE1 CD  CG  CB 

  DONO HN   N
  DONO HE21 NE2
  DONO HE22 NE2
  ACCE OE1  CD
  ACCE O    C
end


residue GLU
  group
    atom N   type=NH1    charge=-0.570 end
    atom HN  type=H      charge= 0.370 end
    atom CA  type=CH1E   charge= 0.200 end
    atom HA  type=HA     charge= 0.000 end
    atom CB  type=CH2E   charge= 0.000 end
    atom HB1 type=HA     charge= 0.000 end
    atom HB2 type=HA     charge= 0.000 end
    atom CG  type=CH2E   charge=-0.100 end
    atom HG1 type=HA     charge= 0.000 end
    atom HG2 type=HA     charge= 0.000 end
    atom CD  type=C      charge= 0.700 end
    atom OE1 type=OC     charge=-0.800 end
    atom OE2 type=OC     charge=-0.800 end
    atom C   type=C      charge= 0.500 end
    atom O   type=O      charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG1     bond CG HG2
  bond CG CD
  bond CD OE1
  bond CD OE2
  bond CA C
  bond C  O

  improper HA N  C   CB    !chirality CA
  improper CD CG OE1 OE2
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG

  dihedral CG  CB  CA  N
  dihedral CD  CG  CB  CA
  dihedral OE1 CD  CG  CB ! mult 2 ! UCL

  DONO HN   N
  ACCE OE1  CD
  ACCE OE2  CD
  ACCE O    C
end


residue GLUH
  group
    atom N   type=NH1    charge=-0.570 end
    atom HN  type=H      charge= 0.370 end
    atom CA  type=CH1E   charge= 0.200 end
    atom HA  type=HA     charge= 0.000 end
    atom CB  type=CH2E   charge= 0.000 end
    atom HB1 type=HA     charge= 0.000 end
    atom HB2 type=HA     charge= 0.000 end
    atom CG  type=CH2E   charge= 0.000 end
    atom HG1 type=HA     charge= 0.000 end
    atom HG2 type=HA     charge= 0.000 end
    atom CD  type=C      charge= 0.500 end
    atom OE1 type=OC     charge=-0.500 end
    atom OE2 type=OH1    charge=-0.450 end
    atom HE2 type=H      charge= 0.450 end
    atom C   type=C      charge= 0.500 end
    atom O   type=O      charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG1     bond CG HG2
  bond CG CD
  bond CD OE1
  bond CD OE2    bond OE2 HE2
  bond CA C
  bond C  O

  improper HA N  C   CB    !chirality CA
  improper CD CG OE1 OE2
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG

  dihedral CG  CB  CA  N
  dihedral CD  CG  CB  CA
  dihedral OE1 CD  CG  CB ! mult 2 ! UCL

  DONO HN   N
  DONO HE2  OE2
  ACCE OE1  CD
  ACCE OE2  CD
  ACCE O    C
end


residue GLY
  group
    atom N   type=NH1    charge=-0.570 end
    atom HN  type=H      charge= 0.370 end
    atom CA  type=CH2G   charge= 0.200 end
    atom HA1 type=HA     charge= 0.000 end
    atom HA2 type=HA     charge= 0.000 end
    atom C   type=C      charge= 0.500 end
    atom O   type=O      charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA1     bond CA HA2
  bond CA C
  bond C  O

  DONO HN   N
  ACCE O    C

  improper HA1 HA2 N C  !stereo CA
end

residue NEP
  group
    atom N   type=NH1     charge=-0.570 end
    atom HN  type=H       charge= 0.370 end
    atom CA  type=CH1E    charge= 0.200 end
    atom HA  type=HA      charge= 0.000 end
    atom CB  type=CH2E    charge= 0.000 end
    atom HB1 type=HA      charge= 0.000 end
    atom HB2 type=HA      charge= 0.000 end
    atom CG  type=C5      charge= 0.330 end
    atom ND1 type=NH1     charge=-0.540 end
    atom HD1 type=H       charge= 0.460 end
    atom CD2 type=CR1H    charge= 0.215 end
    atom HD2 type=HA      charge= 0.115 end
    atom CE1 type=CRHH    charge= 0.385 end
    atom HE1 type=HA      charge= 0.115 end
    atom NE2 type=NR      charge=-0.440 end
    atom P   type=P       charge= 0.520 end
    atom O1P type=O1P     charge=-0.720 end
    atom O2P type=O1P     charge=-0.720 end
    atom O3P type=O3R     charge=-0.720 end
    atom C   type=C       charge= 0.500 end
    atom O   type=O       charge=-0.500 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG
  bond CG  ND1    bond ND1 HD1
  bond ND1 CE1    bond CE1 HE1
  bond CG  CD2    bond CD2 HD2
  bond CD2 NE2    bond NE2 P  
  bond CE1 NE2
  bond P   O1P  
  bond P   O2P    
  bond P   O3P    
  bond CA  C
  bond C   O

  improper HA  N   C   CB  !chirality CA
  improper CG  CB  ND1 CD2
  improper ND1 CE1 CG  HD1
  improper CD2 NE2 CG  HD2
  improper CE1 ND1 NE2 HE1
  improper CG  ND1 CE1 NE2
  improper ND1 CE1 NE2 CD2
  improper CE1 NE2 CD2 CG
  improper NE2 CD2 CG  ND1
  improper CD2 CG  ND1 CE1
  improper HB1 HB2 CA CG  !stereo CB
  improper NE2   P     CD2   CE1  

  dihedral CG  CB  CA  N
  dihedral ND1 CG  CB  CA
  dihedral CD2 CG  CB  CA
  dihedral O3P P   NE2 CE1  

  DONO HN   N
  DONO HD1  ND1
  DONO HE2  NE2
  ACCE O    C
  ACCE O1P  P
  ACCE O2P  P
  ACCE O3P  P
end


residue HISH
  group
    atom N   type=NH1     charge=-0.570 end
    atom HN  type=H       charge= 0.370 end
    atom CA  type=CH1E    charge= 0.200 end
    atom HA  type=HA      charge= 0.000 end
    atom CB  type=CH2E    charge= 0.000 end
    atom HB1 type=HA      charge= 0.000 end
    atom HB2 type=HA      charge= 0.000 end
    atom CG  type=C5      charge= 0.330 end
    atom ND1 type=NH1     charge=-0.540 end
    atom HD1 type=H       charge= 0.460 end
    atom CD2 type=CR1H    charge= 0.215 end
    atom HD2 type=HA      charge= 0.115 end
    atom CE1 type=CRHH    charge= 0.385 end
    atom HE1 type=HA      charge= 0.115 end
    atom NE2 type=NH1     charge=-0.540 end
    atom HE2 TYPE=H       charge= 0.460 end
    atom C   type=C       charge= 0.500 end
    atom O   type=O       charge=-0.500 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG
  bond CG  ND1    bond ND1 HD1
  bond ND1 CE1    bond CE1 HE1
  bond CG  CD2    bond CD2 HD2
  bond CD2 NE2    bond NE2 HE2
  bond CE1 NE2
  bond CA  C
  bond C   O

  improper HA  N   C   CB  !chirality CA
  improper CG  CB  ND1 CD2
  improper ND1 CE1 CG  HD1
  improper CD2 NE2 CG  HD2
  improper CE1 ND1 NE2 HE1
  improper CG  ND1 CE1 NE2
  improper ND1 CE1 NE2 CD2
  improper CE1 NE2 CD2 CG
  improper NE2 CD2 CG  ND1
  improper CD2 CG  ND1 CE1
  improper HB1 HB2 CA CG  ! stereo CB

  dihedral CG  CB  CA  N
  dihedral ND1 CG  CB  CA
  dihedral CD2 CG  CB  CA ! UCL added 09-MAR-00

  DONO HN   N
  DONO HD1  ND1
  DONO HE2  NE2
  ACCE O    C
end

residue HIS
  group
    atom N   type=NH1     charge=-0.570 end
    atom HN  type=H       charge= 0.370 end
    atom CA  type=CH1E    charge= 0.200 end
    atom HA  type=HA      charge= 0.000 end
    atom CB  type=CH2E    charge= 0.000 end
    atom HB1 type=HA      charge= 0.000 end
    atom HB2 type=HA      charge= 0.000 end
    atom CG  type=C5      charge= 0.330 end
    atom ND1 type=NH1     charge=-0.540 end
    atom HD1 type=H       charge= 0.460 end
    atom CD2 type=CR1H    charge= 0.215 end
    atom HD2 type=HA      charge= 0.115 end
    atom CE1 type=CRHH    charge= 0.385 end
    atom HE1 type=HA      charge= 0.115 end
    atom NE2 type=NH1     charge=-0.540 end
    atom HE2 TYPE=H       charge= 0.460 end
    atom C   type=C       charge= 0.500 end
    atom O   type=O       charge=-0.500 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG
  bond CG  ND1    bond ND1 HD1
  bond ND1 CE1    bond CE1 HE1
  bond CG  CD2    bond CD2 HD2
  bond CD2 NE2    bond NE2 HE2
  bond CE1 NE2
  bond CA  C
  bond C   O

  improper HA  N   C   CB  !chirality CA
  improper CG  CB  ND1 CD2
  improper ND1 CE1 CG  HD1
  improper CD2 NE2 CG  HD2
  improper CE1 ND1 NE2 HE1
  improper CG  ND1 CE1 NE2
  improper ND1 CE1 NE2 CD2
  improper CE1 NE2 CD2 CG
  improper NE2 CD2 CG  ND1
  improper CD2 CG  ND1 CE1
  improper HB1 HB2 CA CG  !stereo CB

  dihedral CG  CB  CA  N
  dihedral ND1 CG  CB  CA
  dihedral CD2 CG  CB  CA ! UCL added 09-MAR-00

  DONO HN   N
  DONO HD1  ND1
  DONO HE2  NE2
  ACCE O    C
end

residue HYP   !trans-4-hydroxyproline ! written by KN modified by JR
  group
    atom N    type=N    charge=-0.570 end
    atom CA   type=CH1E charge= 0.285 end
    atom HA   type=HA   charge= 0.000 end
    atom CB   type=CH2E charge= 0.000 end
    atom HB1  type=HA   charge= 0.000 end
    atom HB2  type=HA   charge= 0.000 end
    atom CG   type=CH1P charge= 0.265 end 
    atom HG   type=HA   charge= 0.000 end
    atom CD1  type=CH2P charge= 0.285 end 
    atom HD12 type=HA   charge= 0.000 end   ! 
    atom HD11 type=HA   charge= 0.000 end   ! 
    atom OD2  type=OH   charge=-0.700 end
    atom HD2  type=H    charge= 0.435 end
    atom C    type=C    charge= 0.500 end  
    atom O    type=O    charge=-0.500 end

  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG  
  bond CG OD2    bond OD2 HD2
  bond CG CD1    bond CD1 HD11     bond CD1 HD12
  bond CD1 N
  bond CA C
  bond C  O

  improper HA N C CB       !stereo CA
  improper HB1 HB2 CA CG   !stereo CB
  improper OD2 HG CB CD1   !stereo CG
  improper HD11 HD12 CG N  !stereo CD

  DONO HD2  OD2
  ACCE OD2  " "
  ACCE O    C

end

residue ILE
  group
    atom N    type=NH1     charge=-0.570 end
    atom HN   type=H       charge= 0.370 end
    atom CA   type=CH1E    charge= 0.200 end
    atom HA   type=HA      charge= 0.000 end
    atom CB   type=CH1E    charge= 0.000 end
    atom HB   type=HA      charge= 0.000 end
    atom CG1  type=CH2E    charge= 0.000 end
    atom HG11 type=HA      charge= 0.000 excl = (HG21 HG22 HG23 HD11 HD12 HD13) end
    atom HG12 type=HA      charge= 0.000 excl = (HG21 HG22 HG23 HD11 HD12 HD13) end
    atom CG2  type=CH3E    charge= 0.000 end
    atom HG21 type=HA      charge= 0.000 excl = (HG11 HG12 HD11 HD12 HD13) end
    atom HG22 type=HA      charge= 0.000 excl = (HG11 HG12 HD11 HD12 HD13) end
    atom HG23 type=HA      charge= 0.000 excl = (HG11 HG12 HD11 HD12 HD13) end
    atom CD1  type=CH3E    charge= 0.000 end
    atom HD11 type=HA      charge= 0.000 excl = (HG21 HG22 HG23 HG11 HG12) end
    atom HD12 type=HA      charge= 0.000 excl = (HG21 HG22 HG23 HG11 HG12) end
    atom HD13 type=HA      charge= 0.000 excl = (HG21 HG22 HG23 HG11 HG12) end
    atom C    type=C       charge= 0.500 end
    atom O    type=O       charge=-0.500 end
  
  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB
  bond CB  CG1    bond CG1 HG11    bond CG1 HG12
  bond CB  CG2    bond CG2 HG21    bond CG2 HG22    bond CG2 HG23
  bond CG1 CD1    bond CD1 HD11    bond CD1 HD12    bond CD1 HD13
  bond CA  C
  bond C   O

  improper HA N  C   CB   !chirality CA
  improper HB CA CG2 CG1  !chirality CB
  
  improper HG11 HG12 CB CD1  !stereo CG1
  improper HG21 HG22 CB HG23   !methyl CG2
  improper HD11 HD12 CG1 HD13  !methyl CD1

  dihedral CG1 CB  CA  N
  dihedral CD1 CG1 CB  CA
  dihedral HD11 CD1 CG1 CB   ! UCL methyl stagger 12-MAR-00 
  dihedral HG21 CG2 CB  CA   ! UCL methyl stagger 12-MAR-00 

  DONO HN   N
  ACCE O    C
end


residue LEU
  group
    atom N    type=NH1    charge=-0.570 end
    atom HN   type=H      charge= 0.370 end
    atom CA   type=CH1E   charge= 0.200 end
    atom HA   type=HA     charge= 0.000 end
    atom CB   type=CH2E   charge= 0.000 end
    atom HB1  type=HA     charge= 0.000 end
    atom HB2  type=HA     charge= 0.000 end
    atom CG   type=CH1E   charge= 0.000 end
    atom HG   type=HA     charge= 0.000 end
    atom CD1  type=CH3E   charge= 0.000 end
    atom HD11 type=HA     charge= 0.000 excl = (HD21 HD22 HD23 HG) end
    atom HD12 type=HA     charge= 0.000 excl = (HD21 HD22 HD23 HG) end
    atom HD13 type=HA     charge= 0.000 excl = (HD21 HD22 HD23 HG) end
    atom CD2  type=CH3E   charge= 0.000 end
    atom HD21 type=HA     charge= 0.000 excl = (HD11 HD12 HD13 HG) end
    atom HD22 type=HA     charge= 0.000 excl = (HD11 HD12 HD13 HG) end
    atom HD23 type=HA     charge= 0.000 excl = (HD11 HD12 HD13 HG) end
    atom C    type=C      charge= 0.500 end
    atom O    type=O      charge=-0.500 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG     bond CG  HG
  bond CG  CD1    bond CD1 HD11    bond CD1 HD12    bond CD1 HD13
  bond CG  CD2    bond CD2 HD21    bond CD2 HD22    bond CD2 HD23
  bond CA  C
  bond C   O

  improper HA   N    C   CB    !chirality CA
  improper HG   CB   CD1 CD2   !stereo CG
  improper HB1 HB2 CA CG       !stereo CB
  improper HD11 HD12 CG HD13   !methyl CD1
  improper HD21 HD22 CG HD23   !methyl CD2
    
  dihedral CG  CB  CA  N
  dihedral CD1 CG  CB  CA
  dihedral HD11 CD1 CG  CB   ! UCL methyl stagger 12-MAR-00 
  dihedral HD21 CD2 CG  CB   ! UCL methyl stagger 12-MAR-00 

  DONO HN   N
  ACCE O    C
end


residue LYS
  group
    atom N   type=NH1    charge=-0.570 end
    atom HN  type=H      charge= 0.370 end
    atom CA  type=CH1E   charge= 0.200 end
    atom HA  type=HA     charge= 0.000 end
    atom CB  type=CH2E   charge= 0.000 end
    atom HB1 type=HA     charge= 0.000 end
    atom HB2 type=HA     charge= 0.000 end
    atom CG  type=CH2E   charge= 0.000 end
    atom HG1 type=HA     charge= 0.000 end
    atom HG2 type=HA     charge= 0.000 end
    atom CD  type=CH2E   charge= 0.000 end
    atom HD1 type=HA     charge= 0.000 end
    atom HD2 type=HA     charge= 0.000 end
    atom CE  type=CH2E   charge= 0.310 end
    atom HE1 type=HA     charge= 0.000 end
    atom HE2 type=HA     charge= 0.000 end
    atom NZ  type=NH3    charge=-0.300 end
    atom HZ1 type=HC     charge= 0.330 end
    atom HZ2 type=HC     charge= 0.330 end
    atom HZ3 type=HC     charge= 0.330 end
    atom C   type=C      charge= 0.500 end
    atom O   type=O      charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG1     bond CG HG2
  bond CG CD     bond CD HD1     bond CD HD2
  bond CD CE     bond CE HE1     bond CE HE2
  bond CE NZ     bond NZ HZ1     bond NZ HZ2     bond NZ HZ3
  bond CA C
  bond C  O

  improper HA N C CB      !chirality CA
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG
  improper HD1 HD2 CG CE  !stereo CD
  improper HE1 HE2 CD NZ  !stereo CE
  improper HZ1 HZ2 CE HZ3 !methyl NZ

  dihedral CG  CB  CA  N
  dihedral CD  CG  CB  CA
  dihedral CE  CD  CG  CB
  dihedral NZ  CE  CD  CG
  dihedral HZ1 NZ  CE  CD ! UCL stagger NH3 group 12-MAR-00

  DONO HN   N
  DONO HZ1  NZ
  DONO HZ2  NZ
  DONO HZ3  NZ
  ACCE O    C
end

residue MET
  group
    atom N   type=NH1     charge=-0.570 end
    atom HN  type=H       charge= 0.370 end
    atom CA  type=CH1E    charge= 0.200 end
    atom HA  type=HA      charge= 0.000 end
    atom CB  type=CH2E    charge=-0.000 end
    atom HB1 type=HA      charge= 0.000 end
    atom HB2 type=HA      charge= 0.000 end
    atom CG  type=CH2E    charge= 0.235 end
    atom HG1 type=HA      charge= 0.000 end
    atom HG2 type=HA      charge= 0.000 end
    atom SD  type=SM      charge=-0.470 end
    atom CE  type=CH3E    charge= 0.235 end
    atom HE1 type=HA      charge= 0.000 end
    atom HE2 type=HA      charge= 0.000 end
    atom HE3 type=HA      charge= 0.000 end
    atom C   type=C       charge= 0.500 end
    atom O   type=O       charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG1     bond CG HG2
  bond CG SD
  bond SD CE     bond CE HE1     bond CE HE2     bond CE HE3
  bond CA C
  bond C  O

  improper HA N C CB      !chirality CA
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB SD  !stereo CG
  improper HE1 HE2 SD HE3 !methyl CE

  dihedral CG  CB  CA  N
  dihedral SD  CG  CB  CA
  dihedral CE  SD  CG  CB
  dihedral HE1 CE  SD  CG ! UCL methyl stagger 12-MAR-00

  DONO HN   N
  ACCE O    C
end

residue PGL !N-terminal pyroglutamic acid
  group
    atom N   type=N     charge=-0.5700 end
    atom HN  type=H     charge= 0.3700 end
    atom CA  type=CH1E  charge= 0.3200 end
    atom HA  type=HA    charge= 0.0000 end
    atom CB  type=CH2K  charge= 0.0000 end
    atom HB1 type=HA    charge= 0.0000 end
    atom HB2 type=HA    charge= 0.0000 end
    atom CG  type=CH2L  charge=-0.1000 end
    atom HG1 type=HA    charge= 0.0000 end
    atom HG2 type=HA    charge= 0.0000 end
    atom CD  type=CJ    charge= 0.5000 end   
    atom C   type=C     charge= 0.5000 end   
    atom O   type=O     charge=-0.5000 end
    atom OE  type=O     charge=-0.5000 end

  bond N  CA     bond N HN       bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG1     bond CG HG2
  bond CG CD     bond CD OE
  bond CD N
  bond CA C
  bond C  O

  !dihedral CD N CA CB
  !dIhedral CA N CD CG
  !dihedral CA N CD O
  !dihedral CA CB CG CD 
  !dihedral CB CG CD N
  !dihedral CB CG CD OE
   
  improper HA N C CB      !chirality CA
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG
  improper N CD HN CA      !chirality
  improper CD CG OE N     !chirality 
   
  DONO HN   N
  ACCE O    C
  ACCE OE   CD
end


residue PHE
  group
    atom N   type=NH1     charge=-0.570 end
    atom HN  type=H       charge= 0.370 end
    atom CA  type=CH1E    charge= 0.200 end
    atom HA  type=HA      charge= 0.000 end
    atom CB  type=CH2E    charge= 0.000 end
    atom HB1 type=HA      charge= 0.000 end
    atom HB2 type=HA      charge= 0.000 end
    atom CG  type=CF      charge= 0.000 exclude=(CZ) end
    atom CD1 type=CR1E    charge= -0.115 exclude=(CE2) end
    atom HD1 type=HA      charge= 0.115 end
    atom CD2 type=CR1E    charge= -0.115 exclude=(CE1) end
    atom HD2 type=HA      charge= 0.115 end
    atom CE1 type=CR1E    charge= -0.115 exclude=(CD2) end
    atom HE1 type=HA      charge= 0.115 end
    atom CE2 type=CR1E    charge= -0.115 exclude=(CD1) end
    atom HE2 type=HA      charge= 0.115 end
    atom CZ  type=CR1E    charge= -0.115 exclude=(CG) end
    atom HZ  type=HA      charge= 0.115 end
    atom C   type=C       charge= 0.500 end
    atom O   type=O       charge=-0.500 end
  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG
  bond CG  CD1    bond CD1 HD1
  bond CG  CD2    bond CD2 HD2
  bond CD1 CE1    bond CE1 HE1
  bond CD2 CE2    bond CE2 HE2
  bond CE1 CZ     bond CZ  HZ
  bond CE2 CZ
  bond CA  C
  bond C   O

  improper HA  N   C   CB !chirality CA
  improper HB1 HB2 CA CG  !stereo CB

! Hs and CB around the ring
  improper HD2 CD2 CE2 CZ
  improper HE2 CE2 CZ  CE1
  improper HZ  CZ  CE1 CD1
  improper HE1 CE1 CD1 CG
  improper HD1 CD1 CG  CD2
  improper CB  CG  CD2 CE2

! around the ring
  improper CG  CD1 CE1 CZ 
  improper CD1 CE1 CZ  CE2
  improper CE1 CZ  CE2 CD2
  improper CZ  CE2 CD2 CG
  improper CE2 CD2 CG  CD1
  improper CD2 CG  CD1 CE1

  dihedral CG  CB  CA  N
  dihedral CD1 CG  CB  CA

  DONO HN   N
  ACCE O    C
end

residue PRG !propargylglycine -charges guessed KJR
  group
    atom N    type=NH1    charge=-0.570 end
    atom HN   type=H      charge= 0.370 end
    atom CA   type=CH1E   charge= 0.200 end
    atom HA   type=HA     charge= 0.000 end
    atom CB   type=CH2E   charge=-0.000 end
    atom HB1  type=HA     charge= 0.000 end
    atom HB2  type=HA     charge= 0.000 end
    atom CG   type=CQ     charge= 0.000 end
    atom CD   type=CL     charge= 0.000 end
    atom HD   type=HA     charge= 0.000 end
    atom C    type=C      charge= 0.500 end
    atom O    type=O      charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA  HA
  bond CA CB     bond CB  HB1     bond CB  HB2
  bond CB CG
  bond CG CD
  bond CD HD
  bond CA C
  bond C  O

  improper HA  N  C    CB  ! chirality CA
  improper CG  CB OD1  ND2
  improper ND2 CG HD21 HD22
  improper CB  CG ND2  HD21
  improper HB1 HB2 CA CG  ! stereo CB

  dihedral CG  CB  CA  N
  dihedral OD1 CG  CB  CA ! mult 2 ! UCL

  DONO HN   N
  ACCE O    C
end

residue PRO
  group
    atom N   type=N     charge=-0.570 end
    atom CA  type=CH1E  charge= 0.285 end
    atom HA  type=HA    charge= 0.000 end
    atom CB  type=CH2E  charge= 0.000 end
    atom HB1 type=HA    charge= 0.000 end
    atom HB2 type=HA    charge= 0.000 end
    atom CG  type=CH2P  charge= 0.000 end
    atom HG1 type=HA    charge= 0.000 end
    atom HG2 type=HA    charge= 0.000 end
    atom CD  type=CH2P  charge= 0.285 end
    atom HD1 type=HA    charge= 0.000 end   
    atom HD2 type=HA    charge= 0.000 end   
    atom C   type=C     charge= 0.500 end   
    atom O   type=O     charge=-0.500 end

  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG1     bond CG HG2
  bond CG CD     bond CD HD1     bond CD HD2
  bond CD N
  bond CA C
  bond C  O

  improper HA N C CB      !chirality CA
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG
  improper HD1 HD2 CG N   !stereo CD

  ACCE O    C
end


residue SER
  group
    atom N   type=NH1    charge=-0.570 end
    atom HN  type=H      charge= 0.370 end
    atom CA  type=CH1E   charge= 0.200 end
    atom HA  type=HA     charge= 0.000 end
    atom CB  type=CH2E   charge= 0.265 end
    atom HB1 type=HA     charge= 0.000 end
    atom HB2 type=HA     charge= 0.000 end
    atom OG  type=OH1    charge=-0.700 end
    atom HG  type=H      charge= 0.435 end
    atom C   type=C      charge= 0.500 end
    atom O   type=O      charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB OG     bond OG HG
  bond O  C
  bond C  CA

  improper HA N C CB      !chirality CA
  improper HB1 HB2 CA OG  !stereo CB

  dihedral OG  CB  CA  N

  DONO HN   N
  DONO HG   OG
  ACCE OG  " "
  ACCE O    C
end


residue THR
  group
    atom N    type=NH1 charge=-0.570 end
    atom HN   type=H   charge= 0.370 end
    atom CA   type=CH1E  charge= 0.200 end
    atom HA   type=HA  charge= 0.000 end
    atom CB   type=CH1E  charge= 0.265 end
    atom HB   type=HA  charge= 0.000 end
    atom OG1  type=OH1  charge=-0.700 end
    atom HG1  type=H   charge= 0.435 end
    atom CG2  type=CH3E  charge=-0.000 end
    atom HG21 type=HA  charge= 0.000 end
    atom HG22 type=HA  charge= 0.000 end
    atom HG23 type=HA  charge= 0.000 end
    atom C    type=C   charge= 0.500 end
    atom O    type=O   charge=-0.500 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB
  bond CB OG1    bond OG1 HG1
  bond CB CG2    bond CG2 HG21    bond CG2 HG22    bond CG2 HG23
  bond CA C
  bond C  O

  improper HA N  C   CB        !chirality CA
  improper HB CA OG1 CG2       !stereo CB
  improper HG21 HG22 CB HG23   !methyl CG2

  dihedral OG1 CB  CA  N
  dihedral HG21 CG2 CB  CA   ! UCL methyl stagger 12-MAR-00 

  DONO HN   N
  DONO HG1  OG1
  ACCE OG1 " "
  ACCE O    C
end


residue TRP
  group
    atom N   type=NH1 charge=-0.570 end
    atom HN  type=H   charge= 0.370 end
    atom CA  type=CH1E  charge= 0.200 end
    atom HA  type=HA  charge= 0.000 end
    atom CB  type=CH2E  charge= 0.000 end
    atom HB1 type=HA  charge= 0.000 end
    atom HB2 type=HA  charge= 0.000 end
    atom CG  type=C5W  charge=0.075 end
    atom CD1 type=CR1E  charge= -0.115 end
    atom HD1 type=HA  charge= 0.115 end
    atom CD2 type=CW  charge=-0.055 exclude=(CH2) end
    atom NE1 type=NH1 charge=-0.570 end
    atom HE1 type=H   charge= 0.420 end
    atom CE2 type=CW  charge= 0.130 exclude=(CZ3) end
    atom CE3 type=CR1E  charge= -0.115 exclude=(CZ2) end
    atom HE3 type=HA  charge= 0.115 end
    atom CZ2 type=CR1W charge= -0.115 exclude=(CE3) end
    atom HZ2 type=HA  charge= 0.115 end
    atom CZ3 type=CR1E  charge= -0.115 exclude=(CE2) end
    atom HZ3 type=HA  charge= 0.115 end
    atom CH2 type=CR1W charge= -0.115 exclude=(CD2) end
    atom HH2 type=HA  charge= 0.100 end
    atom C   type=C   charge= 0.500 end
    atom O   type=O   charge=-0.500 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG
  bond CG  CD1    bond CD1 HD1
  bond CG  CD2
  bond CD1 NE1    bond NE1 HE1
  bond NE1 CE2
  bond CD2 CE2
  bond CD2 CE3    bond CE3 HE3
  bond CE2 CZ2    bond CZ2 HZ2
  bond CE3 CZ3    bond CZ3 HZ3
  bond CZ2 CH2    bond CH2 HH2
  bond CZ3 CH2
  bond CA  C
  bond C   O

! chirality
  improper HA  N   C   CB !chirality CA
  improper HB1 HB2 CA CG  !stereo CB

! around the 6-ring
  improper CD2 CE2 CZ2 CH2
  improper CE2 CZ2 CH2 CZ3
  improper CZ2 CH2 CZ3 CE3
  improper CH2 CZ3 CE3 CD2
  improper CZ3 CE3 CD2 CE2
  improper CE3 CD2 CE2 CZ2

! link 5-ring to 6-ring
  improper CD1 NE1 CE2 CZ2
  improper CD1 CG  CD2 CE3
  improper NE1 CE2 CZ2 CH2
  improper NE1 CE2 CD2 CE3
  improper CG  CD2 CE3 CZ3
  improper CG  CD2 CE2 CZ2

! 6-ring hydrogens
  improper HZ2 CZ2 CH2 CZ3
  improper HH2 CH2 CZ3 CE3
  improper HZ3 CZ3 CH2 CZ2
  improper HE3 CE3 CZ3 CH2

! 5-ring hydrogens and CB
  improper HE1 NE1 CE2 CD2
  improper HD1 CD1 NE1 CE2
  improper CB  CG  CD2 CE2

  dihedral CG  CB  CA  N
  dihedral CD1 CG  CB  CA

  DONO HN   N
  DONO HE1  NE1
  ACCE O    C
end

residue TYR
  group
    atom N   type=NH1 charge=-0.570 end
    atom HN  type=H   charge= 0.370 end
    atom CA  type=CH1E  charge= 0.200 end
    atom HA  type=HA  charge= 0.000 end
    atom CB  type=CH2E  charge= 0.000 end
    atom HB1 type=HA  charge= 0.000 end
    atom HB2 type=HA  charge= 0.000 end
    atom CG  type=CY  charge= 0.000 exclude=(CZ) end
    atom CD1 type=CR1E  charge= -0.115 exclude=(CE2) end
    atom HD1 type=HA  charge= 0.115 end
    atom CD2 type=CR1E  charge= -0.115 exclude=(CE1) end
    atom HD2 type=HA  charge= 0.115 end
    atom CE1 type=CR1E  charge= -0.115 exclude=(CD2) end
    atom HE1 type=HA  charge= 0.115 end
    atom CE2 type=CR1E  charge= -0.115 exclude=(CD1) end
    atom HE2 type=HA  charge= 0.115 end
    atom CZ  type=CY2   charge= 0.150 exclude=(CG) end
    atom OH  type=OH1  charge=-0.585 end
    atom HH  type=H   charge= 0.435 end
    atom C   type=C   charge= 0.500 end
    atom O   type=O   charge=-0.500 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG
  bond CG  CD1    bond CD1 HD1
  bond CG  CD2    bond CD2 HD2
  bond CD1 CE1    bond CE1 HE1
  bond CD2 CE2    bond CE2 HE2
  bond CE1 CZ
  bond CE2 CZ
  bond CZ  OH     bond OH  HH
  bond CA  C
  bond C   O

! chirality
  improper HA  N   C   CB !chirality CA
  improper HB1 HB2 CA CG  !stereo CB

! Hs, OH, and CB around the ring
  improper HD2 CD2 CE2 CZ
  improper HE2 CE2 CZ  CE1
  improper OH  CZ  CE1 CD1
  improper HE1 CE1 CD1 CG
  improper HD1 CD1 CG  CD2
  improper CB  CG  CD2 CE2

! around the ring
  improper CG  CD1 CE1 CZ 
  improper CD1 CE1 CZ  CE2
  improper CE1 CZ  CE2 CD2
  improper CZ  CE2 CD2 CG
  improper CE2 CD2 CG  CD1
  improper CD2 CG  CD1 CE1

  dihedral CG  CB  CA  N
  dihedral CD1 CG  CB  CA
  dihedral CE2 CZ  OH  HH   ! UCL Added 12-MAR-00

  DONO HN   N
  DONO HH   OH
  ACCE OH  " "
  ACCE O    C
end


residue VAL
  group
    atom N    type=NH1 charge=-0.570 end
    atom HN   type=H   charge= 0.370 end
    atom CA   type=CH1E  charge= 0.200 end
    atom HA   type=HA  charge= 0.000 end
    atom CB   type=CH1E  charge= 0.000 end
    atom HB   type=HA  charge= 0.000 end
    atom CG1  type=CH3E  charge= 0.000 end
    atom HG11 type=HA  charge= 0.000 excl = (HG21 HG22 HG23) end
    atom HG12 type=HA  charge= 0.000 excl = (HG21 HG22 HG23) end
    atom HG13 type=HA  charge= 0.000 excl = (HG21 HG22 HG23) end
    atom CG2  type=CH3E  charge=-0.000 end
    atom HG21 type=HA  charge= 0.000 excl = (HG11 HG12 HG13) end
    atom HG22 type=HA  charge= 0.000 excl = (HG11 HG12 HG13) end
    atom HG23 type=HA  charge= 0.000 excl = (HG11 HG12 HG13) end
    atom C    type=C   charge= 0.500 end
    atom O    type=O   charge=-0.500 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB
  bond CB  CG1    bond CG1 HG11    bond CG1 HG12    bond CG1 HG13
  bond CB  CG2    bond CG2 HG21    bond CG2 HG22    bond CG2 HG23
  bond CA  C
  bond C   O

  improper HA N C CB         !chirality CA
  improper HB   CA   CG1 CG2 !stereo CB
  improper HG11 HG12 CB HG13    !methyl G1
  improper HG21 HG22 CB HG23    !methyl G2

  dihedral CG1 CB  CA  N
  dihedral HG11 CG1 CB  CA   ! UCL methyl stagger 12-MAR-00 
  dihedral HG21 CG2 CB  CA   ! UCL methyl stagger 12-MAR-00 

  DONO HN   N
  ACCE O    C
end




residue CHEX  !! ADDED BY MN
  group
    atom N   type=NH1 charge=-0.360 end
    atom HN  type=H   charge= 0.260 end
    atom CA  type=CH1E  charge= 0.000 end
    atom HA  type=HA  charge= 0.100 end
    atom CB  type=CH2E  charge=-0.200 end
    atom HB1 type=HA  charge= 0.100 end
    atom HB2 type=HA  charge= 0.100 end
    atom CG  type=CH1E  charge=-0.200 end
    atom HG  type=HA  charge= 0.100 end
    atom CD1  type=CH2E  charge=-0.200 end
    atom HD11 type=HA  charge= 0.100 end
    atom HD12 type=HA  charge= 0.100 end
    atom CD2  type=CH2E  charge=-0.200 end
    atom HD21 type=HA  charge= 0.100 end
    atom HD22 type=HA  charge= 0.100 end
    atom CE1  type=CH2E  charge=-0.200 end
    atom HE11 type=HA  charge= 0.100 end
    atom HE12 type=HA  charge= 0.100 end
    atom CE2  type=CH2E  charge=-0.200 end
    atom HE21 type=HA  charge= 0.100 end
    atom HE22 type=HA  charge= 0.100 end
    atom CZ  type=CH2E  charge=-0.200 end
    atom HZ1 type=HA  charge= 0.100 end
    atom HZ2 type=HA  charge= 0.100 end
    atom C   type=C   charge= 0.480 end
    atom O   type=O   charge=-0.480 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG CD1     bond CG CD2     bond CG HG    
  bond CD1 CE1   bond CD1 HD11    bond CD1 HD12   
  bond CD2 CE2   bond CD2 HD21    bond CD2 HD22
  bond CE1 CZ    bond CE1 HE11    bond CE1 HE12
  bond CE2 CZ    bond CE2 HE21    bond CE2 HE22
  bond CZ HZ1    bond CZ  HZ2
  bond CA C
  bond C  O

  improper HA N C CB          !chirality CA
  improper HB1 HB2 CA CG  !stereo CB
  improper HG CB CD1 CD2      !stereo CG
  improper HD11 HD12 CG CE1      !stereo CD1
  improper HD21 HD22 CG CE2      !stereo CD2
  improper HE11 HE12 CE1 CZ      !stereo CE1
  improper HE11 HE22 CE2 CZ      !stereo CE2
  improper HZ1 HZ2 CE1 CE2  !stereo CZ

  dihedral CG  CB  CA  N
  dihedral CD1  CG  CB  CA

  DONO HN   N
  ACCE O    C
end





residue ACE
  group
    atom CA  type=CH3E charge= 0.000 end
    atom HA1 type=HA charge= 0.000 end
    atom HA2 type=HA charge= 0.000 end
    atom HA3 type=HA charge= 0.000 end
    atom C   type=C  charge= 0.500 end
    atom O   type=O  charge=-0.500 end

  bond C  CA     bond CA HA1     bond CA HA2     bond CA HA3
  bond C  O

  improper HA1 HA2 C HA3  !methyl CA

  ACCE O    C
end




presidue NTER                      ! patch as "NTER - *" to any except PRO
  group
    modify    atom +N   type=NH3 charge=-0.300 end
    delete    atom +HN                        end
    add       atom +HT1 type=HC  charge= 0.330 end
    add       atom +HT2 type=HC  charge= 0.330 end
    add       atom +HT3 type=HC  charge= 0.330 end
    modify    atom +CA           charge= 0.310 end

  add bond +HT1 +N
  add bond +HT2 +N
  add bond +HT3 +N

  add angle +HT1 +N +HT2
  add angle +HT2 +N +HT3
  add angle +HT2 +N +CA
  add angle +HT1 +N +HT3
  add angle +HT1 +N +CA
  add angle +HT3 +N +CA
  add improper +HT1 +HT2 +CA +HT3  !methyl N ???

  add dihedral +C  +CA  +N  +HT1      ! UCL (12-MAR-00)

  ADD DONOr +HT1  +N
  ADD DONOr +HT2  +N
  ADD DONOr +HT3  +N
end


presidue PROP                        ! N-terminal for PRO: "PROP - PRO"
! the charges are guessed from making the residue have charge +1
! and keeping the same charge on N as in NTER (Michael Nilges)
  modify    atom +CD           charge= 0.320 end
  modify    atom +CA           charge= 0.320 end
  modify    atom +N   type=NH3 charge=-0.300 end
  add       atom +HT1 type=HC  charge= 0.330 end
  add       atom +HT2 type=HC  charge= 0.330 end

  add bond +HT1 +N
  add bond +HT2 +N

  add angle +HT1 +N +HT2
  add angle +HT2 +N +CA
  add angle +HT1 +N +CD
  add angle +HT1 +N +CA
  add angle +CD  +N +HT2
  add improper +HT1 +HT2 +CA +CD  !stereo N 
 
  ADD DONOr +HT1  +N
  ADD DONOr +HT2  +N
end
 

presidue CTER               ! C-terminal for all amino acids "* - CTER"
! charge on -C changed from 0.700 to make group -1, Michael Nilges
  group
    modify    atom -CA            charge= 0.100  end
    modify    atom -C             charge= 0.700  end
    delete    atom -O                            end
    add       atom -OT1 type=OC   charge=-0.800  end
    add       atom -OT2 type=OC   charge=-0.800  end

  add bond -C -OT1
  add bond -C -OT2

  add angle -CA  -C -OT1
  add angle -CA  -C -OT2
  add angle -OT1 -C -OT2

  add improper -C -CA -OT2 -OT1
 
  add dihedral -N -CA -C -OT1      ! UCL (12-MAR-00)

  ADD ACCEptor -OT1 -C
  ADD ACCEptor -OT2 -C
end


presidue CTN                  ! C-terminal for all, CONH2 at end "* - CTN"
! charges not consistent with rest, Michael Nilges
  group
    modify    atom -C           charge= 0.48 end
    modify    atom -O           charge=-0.48 end
    add       atom -NT type=NH2 charge=-0.52 end
    add       atom -H1 type=H   charge= 0.26 end
    add       atom -H2 type=H   charge= 0.26 end

  add bond -C  -NT
  add bond -NT -H1
  add bond -NT -H2

  add angle -CA -C  -NT
  add angle -O  -C  -NT
  add angle -CA -C  -O
  add angle -C  -NT -H1
  add angle -C  -NT -H2
  add angle -H1 -NT -H2

  add improper -C  -CA -NT -O
  add improper -C  -NT -O  -H1
  add improper -NT -H1 -H2 -C

  add dihedral -N -CA -C -NT      ! UCL (12-MAR-00)

  ADD DONOR  -H1 -NT
  ADD DONOR  -H2 -NT
end


presidue PEPT     ! PEPTide bond link, for all except the  *(-) - (+)PRO link
                  ! "*(-) - PEPT - (+)*:
  add bond -C +N

  add angle -CA -C +N
  add angle -O  -C +N
  add angle -C  +N +CA
  add angle -C  +N +HN

  add improper  -C -CA +N -O                 ! planar -C
  add improper  +N -C +CA +HN                ! planar +N
  add improper -CA -C  +N  +CA               ! angle across peptide plane

! phi/psi related topology UCL 05-JUL-01

  add dihedral  -C +N +CA +C  mult 6
  add dihedral  -C +N +CA +CB mult 6
  add dihedral  -CB -CA -C +N mult 6
  add dihedral  -CB -CA -C -O mult 6         ! corrected from +O 16-Sept-02
end

presidue PEPX     ! PEPTide bond link for - PAL + *
                  ! "*(-) - PEPT - (+)*:
  add bond -C +N

  add angle -CA -C +N
  add angle -O  -C +N
  add angle -C  +N +CA
  add angle -C  +N +HN

  add improper  -C -CA +N -O                 ! planar -C
  add improper  +N -C +CA +HN                ! planar +N
  add improper -CA -C  +N  +CA               ! angle across peptide plane

! phi/psi related topology UCL 05-JUL-01

 ! add dihedral  -C +N +CA +C  mult 6
 ! add dihedral  -C +N +CA +CB mult 6
 ! add dihedral  -CB1 -CA -C +N mult 6
 ! add dihedral  -CB1 -CA -C -O mult 6         ! corrected from +O 16-Sept-02
end

presidue PEPY     ! PEPTide bond link for - * + PAL
                  ! "*(-) - PEPT - (+)*:
  add bond -C +N

  add angle -CA -C +N
  add angle -O  -C +N
  add angle -C  +N +CA
  add angle -C  +N +HN

  add improper  -C -CA +N -O                 ! planar -C
  add improper  +N -C +CA +HN                ! planar +N
  add improper -CA -C  +N  +CA               ! angle across peptide plane

! phi/psi related topology UCL 05-JUL-01

 ! add dihedral  -C +N +CA +C  mult 6
 ! add dihedral  -C +N +CA +CB1 mult 6
 ! add dihedral  -CB -CA -C +N mult 6
 ! add dihedral  -CB -CA -C -O mult 6         ! corrected from +O 16-Sept-02
end

presidue PEPP     ! for  ...*(-) - (+)PRO  link, same as PEPT except
                  ! replacement H by CD and improper +N +CA +CD -C
  add bond -C +N

  add angle -CA -C +N
  add angle -O  -C +N
  add angle -C  +N +CA
  add angle -C  +N +CD

!  ADD DIHEdral  -C +N +CA +C  ! taken out since related to phi MN 18-MAR-02
!  ADD DIHEdral  -N -CA -C +N  ! taken out since equivalent below MN 18-MAR-02
!  ADD DIHEdral  -CA -C +N +CA ! taken out since improper below MN 18-MAR-02

  add improper  -C  -CA +N  -O               ! planar -C
  add improper  +N  -C +CA +CD               ! planar +N
  add improper -CA  -C  +N  +CA              ! angle across peptide plane

! psi related topology N.B. phi is fixed by ring bonding 05-JUL-01

  add dihedral  -CB -CA -C +N mult 6
  add dihedral  -CB -CA -C -O mult 6         ! corrected from +O 16-Sept-02
end


presidue PPGP     ! for  GLY(-) - (+)PRO  link, same as PEPT except
                  ! replacement H by CD and improper +N +CA +CD -C
  add bond -C +N

  add angle -CA -C +N
  add angle -O  -C +N
  add angle -C  +N +CA
  add angle -C  +N +CD

!  ADD DIHEdral  -C +N +CA +C  ! taken out since related to phi MN 18-MAR-02
!  ADD DIHEdral  -N -CA -C +N  ! taken out since equivalent below MN 18-MAR-02
!  ADD DIHEdral  -CA -C +N +CA ! taken out since improper below MN 18-MAR-02

  add improper  -C  -CA +N  -O               ! planar -C
  add improper  +N  -C +CA +CD               ! planar +N
  add improper -CA  -C  +N  +CA              ! angle across peptide plane

end


PRESidue PPGG { for GLY (-) - (+) GLY LINK }

  add bond -C +N

  add angle -CA -C +N
  add angle -O  -C +N
  add angle -C  +N +CA
  add angle -C  +N +HN

  add improper  -C -CA +N -O                 ! planar -C
  add improper  +N -C +CA +HN                ! planar +N
  add improper -CA -C  +N  +CA               ! angle across peptide plane

! phi/psi related topology UCL 05-JUL-01

  add dihedral  -C +N +CA +C  mult 6
end

PRESidue PPG1 { for ...*(-) - (+) GLY LINK }

  add bond -C +N

  add angle -CA -C +N
  add angle -O  -C +N
  add angle -C  +N +CA
  add angle -C  +N +HN

  add improper  -C -CA +N -O                 ! planar -C
  add improper  +N -C +CA +HN                ! planar +N
  add improper -CA -C  +N  +CA               ! angle across peptide plane

! phi/psi related topology UCL 05-JUL-01

  add dihedral  -C +N +CA +C  mult 6
  add dihedral  -CB -CA -C +N mult 6
  add dihedral  -CB -CA -C -O mult 6          ! corrected from +O 16-Sept-02

end


PRESidue PPG2 { for ... GLY(-) - (+) * LINK }
  add bond -C +N

  add angle -CA -C +N
  add angle -O  -C +N
  add angle -C  +N +CA
  add angle -C  +N +HN

  add improper  -C -CA +N -O                  ! planar -C
  add improper  +N -C +CA +HN                 ! planar +N
  add improper -CA -C  +N  +CA                ! angle across peptide plane

! phi/psi related topology UCL 05-JUL-01

  add dihedral  -C +N +CA +C  mult 6
  add dihedral  -C +N +CA +CB mult 6
end



presidue DISU                ! disulfide bridge  ...CYS - DISU - CYS...
  group
    delete    atom 1HG               end
    modify    atom 1CB              charge= 0.300 end
    MODIfy    ATOM 1SG  TYPE=S      charge=-0.300 end
  group
    delete    atom 2HG               end
    modify    atom 2CB              charge= 0.300 end
    MODIfy    ATOM 2SG  TYPE=S      charge=-0.300 end

  add bond 1SG 2SG

  add angle 1CB 1SG 2SG
  add angle 1SG 2SG 2CB

  add dihedral 1CB 1SG 2SG 2CB mult 6   ! UCL modified 12-MAR-00 & 07-JUL-01
end

presidue RCM                ! disulfide bridge  ...PAL - RCM - PAL...
  group
    delete    atom 1CH               end
    delete    atom 1CZ               end
    delete    atom 1HH1              end
    delete    atom 1HH2              end
    delete    atom 1HH3              end
    delete    atom 1HZ               end
    modify    atom 1CG              charge= 0.256 end
    modify    atom 1CD              charge= 0.040 end
    modify    atom 1CE              charge=-0.242 end
    modify    atom 1HE              charge= 0.136 end
    
   group
    delete    atom 2CH               end
    delete    atom 2CZ               end
    delete    atom 2HH1              end
    delete    atom 2HH2              end
    delete    atom 2HH3              end
    delete    atom 2HZ               end
    modify    atom 2CG              charge= 0.256 end
    modify    atom 2CD              charge= 0.040 end
    modify    atom 2CE              charge=-0.242 end
    modify    atom 2HE              charge= 0.136 end

  add bond 1CE 2CE

  add angle 1CD 1CE 2CE
  add angle 1CE 2CE 2CD

  add improper 1CB 1SG 2SG 2CB  !planar C-C double bond
end

!presidue DISU                ! triazole bridge  ...PRG - 1,4 triazole - AZH...
!  group
!    modify    atom 1CB              charge= 0.300 end
!    MODIfy    ATOM 1SG  TYPE=S      charge=-0.300 end
!  group
!    delete    atom 2HG               end
!    modify    atom 2CB              charge= 0.300 end
!    MODIfy    ATOM 2SG  TYPE=S      charge=-0.300 end

!  add bond 1ND 2CD
!  add bond 1NZ 2CG
  
!  add angle 1CG 1ND 2CD
!  add angle 1ND 2CD 2CG
!  add angle 1NE 1NZ 2CG
!  add angle 1NZ 2CG 2CB
!  add angle 1ND 2CD 2HD
  
!  add improper 1ND 1CG 1NE 2CD
!  add improper 2CG 2CB 2CD 1NZ
!  add improper 2CD 2HD 1ND 2CG
!  add improper 1ND 2CD 2CG 1NZ
!  add improper 2CD 2CG 1NZ 1NE
!  add improper 2CG 1NZ 1NE 1ND
!  add improper 1NZ 1NE 1ND 2CD
!  add improper 1NE 1ND 2CD 2CG
  
end

presidue TRIZ                ! triazole bridge  ...PRG - 1,4 triazole - AZA...
!  group
!    modify    atom 1CB              charge= 0.300 end
!    MODIfy    ATOM 1SG  TYPE=S      charge=-0.300 end
!  group
!    delete    atom 2HG               end
!    modify    atom 2CB              charge= 0.300 end
!    MODIfy    ATOM 2SG  TYPE=S      charge=-0.300 end

  add bond 1NG 2CD
  add bond 1NE 2CG
  
  add angle 1CG 1NG 2CD
  add angle 1NG 2CD 2CG
  add angle 1ND 1NE 2CG
  add angle 1NE 2CG 2CB
  add angle 1NG 2CD 2HD
  
  add improper 1CG 1NG 2CG 2CD
  add improper 1NG 2CD 2CB 2CG
  add improper 1NG 2CD 1NE 2CG
  add improper 2CD 2CG 1ND 1NE
  add improper 1NG 1ND 2CG 1NE
  
end


presidue DISN                ! disulfide bridge  ...CYS - DISU - CYS...
                             ! w/o the actual bond
  group
    delete    atom 1HG               end
    modify    atom 1CB              charge= 0.300 end
    MODIfy    ATOM 1SG  TYPE=S      charge=-0.300 end
  group
    delete    atom 2HG               end
    modify    atom 2CB              charge= 0.300 end
    MODIfy    ATOM 2SG  TYPE=S      charge=-0.300 end
end


presidue LTOD                        ! change from L to D amino acid
  delete improper HA N C CB
  add    improper HA C N CB
end

presidue CISP                        ! change from trans to cis peptide bond
  modify atom C type=CCIS end
end
presidue CIPP                        ! change from trans to cis peptide bond. not necessary to have separate patch
  modify atom C type=CCIS end
end


!--------------------------------------------------------------------------

PRESidue HISE   ! Patch to change doubly protonated HIS to singly
                ! protonated histidine (HE2)
                ! has to be patched as REFErence=NIL=<selection>
  DELETE DONOR  ND1 HD1
  MODIFY ATOM  CB   TYPE=CH2E    CHARge= 0.000  END
  MODIFY ATOM  CG   TYPE=C5      CHARge= 0.130  END !#
  MODIFY ATOM  ND1  TYPE=NR      CHARge=-0.490  END !#
  MODIFY ATOM  CE1  TYPE=CRH     CHARge= 0.295  END !#
  MODIFY ATOM  CD2  TYPE=CR1E    CHARge= -0.015  END
  MODIFY ATOM  NE2  TYPE=NH1     CHARge=-0.570  END !#
  MODIFY ATOM  HE2  TYPE=H       CHARge= 0.420  END !#
  DELETE ATOM  HD1                              END

  ADD ACCEPTOR  ND1 " "
END {HISE}

!--------------------------------------------------------------------------

PRES HISD   ! Patch to change doubly protonated HIS to singly
            ! protonated histidine (HD1)
            ! has to be patched as REFErence=NIL=<selection>
  DELETE DONOR HE2 NE2
  MODIFY ATOM  CB   TYPE=CH2E    CHARge= 0.000  END
  MODIFY ATOM  CG   TYPE=C5      CHARge= 0.130  END !#
  MODIFY ATOM  ND1  TYPE=NH1     CHARge=-0.570  END  !#
  MODIFY ATOM  HD1  TYPE=H       CHARge= 0.420  END  !#
  MODIFY ATOM  CD2  TYPE=CR1E    CHARge= -0.015  END  !#
  MODIFY ATOM  NE2  TYPE=NR      CHARge=-0.490  END  !#
  MODIFY ATOM  CE1  TYPE=CRH     CHARge= 0.295  END  !#
  DELETE ATOM  HE2                              END

  ADD ACCEPTOR NE2 " "

END {HISD}


PRESidue DCA2   ! Patch to add a single bond from an ASP COO- group to a calcium ion

  ADD BOND 1OD1 2CA+2

END


PRESidue ECA2   ! Patch to add a single bond from a GLU COO- group to a calcium ion

  ADD BOND 1OE1 2CA+2

END
  




set echo=true end
