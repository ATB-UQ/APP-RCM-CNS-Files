# APP-RCM-CNS-Files

# CNS link statement used
```
presidue RCM                ! dicarbabridge between two pentenyl alanine  ...PAL - RCM - PAL...
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
  ```
