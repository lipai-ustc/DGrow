# What is **TEAKS**?
Thermodynamic equilibirum analysis and KMC simulation (**TEAKS**) package is a collection of tools
for the thermodynamic equilibrium analysis and kinetic Monte Carlo simulations on surface
based on first principle calculation results.

# License
Copyright (C) 2020-2026 Pai Li (lipai@mail.ustc.edu.cn)

# Installation
cd source/
make

# Usage
./teaks.x-1.0.0-ifort infile

# TEAKS input file description

Mode KMC
Restart  True  # restart will read Status file to restore the model
Temperature 1300 # unit Kelvin
Vibmode 1    # 1: Vib=kT/h;   0: Vib is read from infile Vib value
Vib    10    # only read when Vibmode is 0
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize    # 2D substrate mesh size 
1000 1000

CoreSize    # 2D island size
100 100

# Total steps and print intervals
TotalStep           1
PrintSpecDenStep    1E6
SaveStatus          1E6

GasNumber 2
GasName
H2 CH4
GasPressure  # in Torr
0.001 1000.0
GasAdsorptionRate   # can be calculatd using related kinetic formula
3570 0.00874
GasAdsorbReaction
#Reactant  Product1  Product2
H2  H   H
CH4 CH3 H

MSpeciesNumber 2  # Species treated using mean-field approximation
# in which, species has no position info and only density is used in the simulation.
MSpeciesName
H  F
InitialMSpeciesDensity  
0.001 0.0001
FixMSpeciesDensity  # 0 means unfixed, 1 means fixed
0    1
MAttachBarrier  # the barrier of species attach onto an island edge-site where neighboring sites are not occupied
0.6  0.8
MAttachBarrier2 # the barrier of species attach onto an island edge-site where one neighboring site is occupied
0.7  0.9
MAttachBarrier3 # the barrier of species attach onto an island edge-site where both neighboring sites are occupied
0.8  1.0
MDetachBarrier # the barrier of species detach from an island edge-site where neighboring sites are not occupied
1.57  2.1
MDetachBarrier2 # the barrier of species detach from an island edge-site where one neighboring site is occupied
1.67  2.2
MDetachBarrier3 # the barrier of species detach from an island edge-site where both neighboring sites are occupied
1.77  2.3

SpeciesNumber  9  #species without MFA
SpeciesName
C CH CH2 CH3 C2 C2H C2H2 C3 C3H
DiffusionBarrier
#negative value indicates event that will not happen in the simulation
#Barriers for species with MFA here will not be used
#but should be listed with whatever values
0.50 0.10 0.18 -1 0.49 0.32 0.44 -1 -1

AttachBarrier  # all barriers are in eV in TEAKS 
1.27  0.44  0.19 -1 0.58  0.56  1.05 -1 -1  
AttachBarrier2
1.57  0.44  0.19 -1 0.58  0.56  1.05  -1 -1
AttachBarrier3
1.87  0.44  0.19 -1 0.58  0.56  1.05  -1 -1
DetachBarrier
1.57  1.08  2.08 -1 2.29  2.31  3.55  -1 -1 
DetachBarrier2
1.57  1.08  2.08 -1 2.29  2.31  3.55  -1 -1 
DetachBarrier3
1.57  1.08  2.08 -1 2.29  2.31  3.55  -1 -1 
#InitialSpeciesDensity
#0.02  0.01  0.03  0.03  0.00  0.00  0.00 
InitialSpeciesNumber
12 22 32 18 42 72 22 42 32 
FixSpeciesDensity
0 1 0 0 1 1 1 0 0 

MergeEventNumber  13   # Merge event

MergeReaction      
#Reactant1 Reactant2 Product  Barrier
H  H   H2   0.8     
H  C   CH   0.79    
H  CH  CH2  0.65    
H  CH2 CH3  0.68    
H  CH3 CH4  0.69    
H  C2  C2H  0.72    
H  C2H C2H2 0.85    
H  C3  C3H  0.57    
C  C   C2   0.25    
C  CH  C2H  1.27    
C  C2  C3   1.14    
CH CH  C2H2 0.14    
CH C2  C3H  0.93    

#If product is gas molecule, the merge reaction corresponds to the associative desorption
#Such as the associative desorption of H2

DecompEventNumber  11   # Decomposition events
DecompReaction
#Reactant  Product1  Product2  Barrier
CH    C     H     1.65   
CH2   CH    H     1.13   
CH3   CH2   H     1.53   
C2    C     C     2.75   
C2H   C2    H     1.38   
C2H   C     CH    2.95   
C2H2  C2H   H     1.64   
C2H2  CH    CH    2.28   
C3    C2    C     2.32   
C3H   C2    CH    1.92   
C3H   C3    H     1.66   
                       
DesorpEventNumber  1    # Desorption event
DesorpReaction
C2H2    1.83
