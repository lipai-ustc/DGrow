![image](https://github.com/user-attachments/assets/250d9d17-0023-47ed-90d5-4eef6ac6a42e)  
lipai@mail.ustc.edu.cn

### 1. Introduction to DGrow

#### 1.1 What is DGrow

DGrow, fully named as Two-Dimensional Material Growth Simulation Software, is an open-source tool that employs kinetic Monte Carlo (KMC) simulation based on elementary reaction barrier information to study the growth of two-dimensional materials on substrates. It simulates processes such as gas molecule adsorption and desorption, diffusion, decomposition, and recombination reactions of multiple species on the substrate surface, attachment and detachment at the edges of two-dimensional materials, and applies mean-field approximation for species with very low diffusion barriers to optimize simulation efficiency. The software is written in Fortran and does not require additional libraries, making it compatible with various computer operating systems.

#### 1.2 Development Status of Kinetic Monte Carlo Simulations

Atomic-scale simulations of material growth are constrained by spatial and temporal limitations. While ab initio molecular dynamics (AIMD) offers highly accurate energy and force calculations, its computational efficiency is limited, restricting the simulated system size to around a hundred atoms. Empirical potential molecular dynamics (EPMD) can simulate up to millions of atoms but lacks precision, especially in describing bond formation and breaking. Neural network (NN) potentials provide a promising simulation method with near-density functional theory (DFT) accuracy and computation speeds slightly slower than empirical potentials. Even with EPMD, rare events in kinetic simulations are challenging to observe without enhanced sampling methods like umbrella sampling. For atomic-scale simulations extending into milliseconds, KMC might be the only viable option. Due to the customization required for different systems, developing a universal KMC code that satisfies all simulation needs remains challenging. Several academic publications have introduced KMC-related programs, each tailored to specific applications.

#### 1.3 Theoretical Foundations and Advantages of DGrow

In environments where two-dimensional materials grow, the substrate surface hosts numerous species with varying properties. Species with very low diffusion barriers are treated using mean-field approximation, described by density scalars rather than individual particle positions, which significantly reduces computational time. Species with higher diffusion barriers are handled normally within the KMC framework, tracking their positions and movements. Unlike molecular catalytic reactions, species on the substrate during growth can be sparse, leading to challenges in reaching equilibrium or steady states. Efficient algorithms must be carefully chosen for KMC models to address these issues. DGrow's KMC framework treats particles as statistical units for traversal, optimizing efficiency for sparse surfaces.

#### 1.4 Algorithm Framework of DGrow

The core algorithm of DGrow is based on kinetic Monte Carlo methods. Its workflow diagram is provided separately.
<img src="https://github.com/user-attachments/assets/40de03b5-0418-475f-928a-b3f1c6c9df42" width="500">


#### 1.5 Main Modules of DGrow

- **Main Program:** Controls overall execution, including file reading, KMC simulation initiation, and basic output.
- **Input Parsing Module:** Parses input files and saves data to the InputData structure.
- **Mesh Class Definition:** Initializes grid structures and related functions, such as obtaining nearest neighbor information.
- **Species Class Definition:** Stores information about surface species and those approximated by mean field, along with relevant reaction functions.
- **Gas Species Class Definition:** Handles initialization and adsorption functions for gaseous species.

#### 1.6 Primary Classes and Attributes of DGrow

Refer to the detailed documentation for class definitions and attributes.  
<img src="https://github.com/user-attachments/assets/46d43d6b-7977-44a0-9007-6207f23cf1b8" width="500">


#### 1.7 Physical Model of DGrow

DGrow focuses on the growth of two-dimensional materials on substrates, represented by a periodic structure where particles exiting one edge re-enter from the opposite side. Particle diffusion events are modeled by jumps between grid points according to the shortest distance.  
<img src="https://github.com/user-attachments/assets/03462e09-71bc-43a2-8ed6-523744dcd975" width="200">


#### 1.8 Events Included in DGrow Simulations

Common species can undergo events such as gas adsorption, surface diffusion, decomposition, recombination, desorption, edge attachment, and detachment. Event probabilities follow the formula \(P=v\times\exp(-E/(k_B T))\), where \(v\) is the pre-exponential factor, \(E\) is the elementary reaction barrier, \(k_B\) is the Boltzmann constant, and \(T\) is temperature.

#### 1.9 Mean Field Approximation for Particles

Particles with extremely low diffusion barriers are approximated using mean field, preserving only total particle numbers and concentrations. They participate in decomposition, recombination, desorption, edge attachment, and detachment events but not diffusion. Reaction frequencies are calculated based on particle densities.

#### 1.10 References

For theoretical details on KMC simulations, refer to:
[1] The Journal of Physical Chemistry C 121 (46), 25949-25955
[2] The Journal of Physical Chemistry C 124 (30), 16233-16247
[3] Chinese Science Bulletin 63 (33), 3419-3426

### 2. Compilation and Execution

#### 2.1 Compilation

To compile the source code, extract the tarball, enter the source directory, adjust the makefile as needed for your compilation environment, and run `make` to produce the executable `dgrow.x`. For quick installation:

```bash
cd source
make
```

#### 2.2 Execution

Place the compiled `dgrow.x` in the `~/bin` folder or add its path to the PATH environment variable in `~/.bashrc`:

```bash
export PATH=$PATH:/your_dgrow_path/
```

Replace `/your_dgrow_path/` with the actual path to `dgrow.x`. To run DGrow, ensure an input file exists in the current directory and use:

```bash
dgrow.x infile
```

### 3. Detailed Input Parameters

Input parameters for DGrow are managed through the `infile`, containing single-value settings and multi-value configurations. Comments start with `#`.

#### 3.1 Basic Settings for Kinetic Monte Carlo (KMC)

- Mode: Currently supports KMC simulation (`Mode KMC`)
- Restart Option: Enables continuation after interruption (`Restart True`)
- Temperature Setting: Units in Kelvin (`Temperature 1300`)
- Vibration Frequency Setting: Determines reaction rate pre-factor (`Vibmode 1`)
- Total Steps for KMC Simulation (`TotalStep 1E8`)
- Step Interval for Screen Output (`PrintSpecDenStep 1E6`)
- Step Interval for Saving State (`SaveStatus 1E6`)

#### 3.2 Mesh Configuration

- Substrate Grid Type (`MeshType 1`)
- Substrate Grid Size (`MeshSize 1000 1000`)
- Two-Dimensional Island Size (`CoreSize 100 100`)

#### 3.3 Gaseous Species Configuration

- Number of Gas Types (`GasNumber 2`)
- Gas Names (`GasName H2 CH4`)
- Gas Pressures (`GasPressure 0.001 1000.0`)
- Adsorption Rates (`GasAdsorptionRate 3570 0.00874`)
- Dissociative Adsorption Reactions (`GasAdsorbReaction`)

#### 3.4 Mean Field Approximation (MFA) Species Configuration

- Number of MFA Species (`MSpeciesNumber 2`)
- Species Names (`MSpeciesName H F`)
- Initial Concentrations (`InitialMSpeciesDensity 0.001 0.0001`)
- Fixed Concentration Flag (`FixMSpeciesDensity 0 1`)
- Attachment and Detachment Barriers for Edge Sites (`MAttachBarrier`, `MDetachBarrier`)

#### 3.5 Regular Species Configuration Without MFA

- Number of Species (`SpeciesNumber 9`)
- Species Names (`SpeciesName C CH CH2 CH3 C2 C2H C2H2 C3 C3H`)
- Diffusion Barriers (`DiffusionBarrier`)
- Attachment and Detachment Barriers (`AttachBarrier`, `DetachBarrier`)
- Initial Concentrations (`InitialSpeciesDensity`)
- Initial Particle Numbers (`InitialSpeciesNumber`)
- Fixed Concentration Flags (`FixSpeciesDensity`)

### 3.6 Reaction Event Settings

#### Merge Reaction Events

- Number of Merge Events (`MergeEventNumber 5`)
- Information on Reactants, Products, and Barriers for Merge Reactions (`MergeReaction`)

```plaintext
#Reactant1 Reactant2 Product  Barrier
H  H   H2   0.8     
H  C   CH   0.79    
H  CH  CH2  0.65    
H  CH2 CH3  0.68    
H  CH3 CH4  0.69    
```

If the product is a gas molecule, the merge reaction corresponds to associative desorption.

#### Decomposition Reaction Events

- Number of Decomposition Events (`DecompEventNumber 3`)
- Information on Reactants, Products, and Barriers for Decomposition Reactions (`DecompReaction`)

```plaintext
#Reactant  Product1  Product2  Barrier
CH    C     H     1.65   
CH2   CH    H     1.13   
CH3   CH2   H     1.53   
```

#### Desorption Events

- Number of Desorption Events (`DesorpEventNumber 1`)
- Species and Desorption Barriers for Desorption Reactions (`DesorpReaction`)

```plaintext
C2H2    1.83
```

### 3.7 Example Input Files

#### Example 1: Equilibrium of Carbon Monomers and Dimers

- **Files:** See `examples/C-C2`
- **Included Events:**
  1. Diffusion of C and C2
  2. C + C -> C2 (merge reaction)
  3. C2 -> C + C (decomposition reaction)

**Input File (`infile`) Content:**

```plaintext
Mode KMC
Temperature 1300 # unit Kelvin
Vibmode 1
Vib    10
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize
2000 2000

TotalStep           1E8

SpeciesNumber  2  #species without MFA
SpeciesName
C  C2
DiffusionBarrier
0.50 0.49
InitialSpeciesNumber
100 100

MergeEventNumber  1
MergeReaction
C  C   C2   0.90

DecompEventNumber  1
DecompReaction
C2    C     C     2.80
```

#### Example 2: Concentration of H Species on Copper Surface Under Different Hydrogen Pressures

- **Files:** See `examples/H-no_MFA` and `examples/H-saturation`
- **Included Events:**
  1. Dissociative adsorption of H2 -> H + H
  2. Diffusion of H
  3. H + H -> H2 (associative desorption)

**Input File (`infile`) Content:**

```plaintext
Mode KMC
Temperature 1300 # unit Kelvin
Vibmode 1
Vib    10
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize
500 500

TotalStep           1E10
PrintSpecDenStep    1E8
SaveStatus          1E8

GasNumber 1
GasName
H2
GasPressure  # in Torr
10
GasAdsorptionRate
3570
GasAdsorbReaction
#Reactant  Product1  Product2
H2  H   H

SpeciesNumber 1
SpeciesName
H
DiffusionBarrier
0.15

MergeEventNumber  1
MergeReaction
H  H   H2   0.8
```

#### Example 3: Same as Example 2 but with Mean Field Approximation for H Species

**Input File (`infile`) Content:**

```plaintext
Mode KMC
Temperature 1300 # unit Kelvin
Vibmode 1
Vib    10
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize
1000 1000

TotalStep           2E5

GasNumber 1
GasName
H2
GasPressure  # in Torr
0.1
GasAdsorptionRate
3570
GasAdsorbReaction
#Reactant  Product1  Product2
H2  H   H

MSpeciesNumber 1
MSpeciesName
H
InitialMSpeciesDensity
0.0001
FixMSpeciesDensity
0

MergeEventNumber  1
MergeReaction
H  H   H2   0.8
```

#### Example 4: Hydrogen Saturation at Graphene Edges Under Different Hydrogen Pressures

- **Files:** See `examples/H-no_MFA` and `examples/H-saturation`
- **Included Events:**
  1. Dissociative adsorption of H2 -> H + H
  2. Diffusion of H
  3. H + H -> H2 (associative desorption)
  4. Attachment of H at three different sites on graphene edges
  5. Detachment of H from three different sites on graphene edges

**Input File (`infile`) Content:**

```plaintext
Mode KMC
Temperature 1300 # unit Kelvin
Vibmode 1
Vib    10
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize
500 500
CoreSize
100 100

TotalStep           2E9
PrintSpecDenStep    1E7
SaveStatus          1E7

GasNumber 1
GasName
H2
GasPressure  # in Torr
0.1
GasAdsorptionRate
3570
GasAdsorbReaction
#Reactant  Product1  Product2
H2  H   H

SpeciesNumber 1
SpeciesName
H
DiffusionBarrier
0.15

AttachBarrier
0.75
AttachBarrier2
0.53
AttachBarrier3
0.84
DetachBarrier
1.74
DetachBarrier2
1.61
DetachBarrier3
1.99

MergeEventNumber  1
MergeReaction
H  H   H2   0.8
```

#### Example 5: Reaction Network of Multiple Carbon-Hydrogen Species with Mean Field Approximation for Species H

- **Files:** See `examples/concentration`
- **Included Events:**
  1. Dissociative adsorption of hydrogen gas: H2 -> H + H
  2. Dissociative adsorption of methane gas: CH4 -> CH3 + H
  3. Diffusion of species other than H
  4. Merge reactions including:
     - H+H->H2 (hydrogen recombination into gas)
     - Various combinations leading to new species or gas release
  5. Decomposition reactions for various species
  6. Desorption of C2H2

**Input File (`infile`) Content:**

```plaintext
Mode KMC
Temperature 1300 # unit Kelvin
Vibmode 1
Vib    10
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize
1000 1000

TotalStep           1E10

GasNumber 2
GasName
H2 CH4
GasPressure  # in Torr
0.001 1000.0
GasAdsorptionRate
3570 0.00874
GasAdsorbReaction
#Reactant  Product1  Product2
H2  H   H
CH4 CH3 H

MSpeciesNumber 1
MSpeciesName
H

SpeciesNumber  9  #species without MFA
SpeciesName
C CH CH2 CH3 C2 C2H C2H2 C3 C3H
DiffusionBarrier
0.50 0.10 0.18 -1 0.49 0.32 0.44 -0.1 -0.1

MergeEventNumber  13
MergeReaction
#Reactant1 Reactant2 Product  Barrier
H  H   H2   0.8     1
H  C   CH   0.79    2
H  CH  CH2  0.65    3
H  CH2 CH3  0.68    4
H  CH3 CH4  0.69    5
H  C2  C2H  0.72    6
H  C2H C2H2 0.85    7
H  C3  C3H  0.57    8
C  C   C2   0.25    9
C  CH  C2H  1.27    10
C  C2  C3   1.14    11
CH CH  C2H2 0.14    12
CH C2  C3H  0.93    13

DecompEventNumber  11
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

DesorpEventNumber  1
DesorpReaction
C2H2    1.83
```

### 4. Detailed Output Explanation

DGrow outputs information both to the screen and to files. Screen output includes summaries of calculation parameters and basic KMC run information. File outputs contain the main results of the KMC simulation. Below are descriptions of each output file:

#### 4.1 SpecNum

The `SpecNum` file records the number of particles for all species at each time point.

#### 4.2 Density

The `Density` file records the particle density for all species at each time point.

#### 4.3 Statistics

The `Statistics` file records the occurrence frequency of all KMC events.### 3.6 Reaction Event Settings

#### Merge Reaction Events

- Number of Merge Events (`MergeEventNumber 5`)
- Information on Reactants, Products, and Barriers for Merge Reactions (`MergeReaction`)

```plaintext
#Reactant1 Reactant2 Product  Barrier
H  H   H2   0.8     
H  C   CH   0.79    
H  CH  CH2  0.65    
H  CH2 CH3  0.68    
H  CH3 CH4  0.69    
```

If the product is a gas molecule, the merge reaction corresponds to associative desorption.

#### Decomposition Reaction Events

- Number of Decomposition Events (`DecompEventNumber 3`)
- Information on Reactants, Products, and Barriers for Decomposition Reactions (`DecompReaction`)

```plaintext
#Reactant  Product1  Product2  Barrier
CH    C     H     1.65   
CH2   CH    H     1.13   
CH3   CH2   H     1.53   
```

#### Desorption Events

- Number of Desorption Events (`DesorpEventNumber 1`)
- Species and Desorption Barriers for Desorption Reactions (`DesorpReaction`)

```plaintext
C2H2    1.83
```

### 3.7 Example Input Files

#### Example 1: Equilibrium of Carbon Monomers and Dimers

- **Files:** See `examples/C-C2`
- **Included Events:**
  1. Diffusion of C and C2
  2. C + C -> C2 (merge reaction)
  3. C2 -> C + C (decomposition reaction)

**Input File (`infile`) Content:**

```plaintext
Mode KMC
Temperature 1300 # unit Kelvin
Vibmode 1
Vib    10
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize
2000 2000

TotalStep           1E8

SpeciesNumber  2  #species without MFA
SpeciesName
C  C2
DiffusionBarrier
0.50 0.49
InitialSpeciesNumber
100 100

MergeEventNumber  1
MergeReaction
C  C   C2   0.90

DecompEventNumber  1
DecompReaction
C2    C     C     2.80
```

#### Example 2: Concentration of H Species on Copper Surface Under Different Hydrogen Pressures

- **Files:** See `examples/H-no_MFA` and `examples/H-saturation`
- **Included Events:**
  1. Dissociative adsorption of H2 -> H + H
  2. Diffusion of H
  3. H + H -> H2 (associative desorption)

**Input File (`infile`) Content:**

```plaintext
Mode KMC
Temperature 1300 # unit Kelvin
Vibmode 1
Vib    10
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize
500 500

TotalStep           1E10
PrintSpecDenStep    1E8
SaveStatus          1E8

GasNumber 1
GasName
H2
GasPressure  # in Torr
10
GasAdsorptionRate
3570
GasAdsorbReaction
#Reactant  Product1  Product2
H2  H   H

SpeciesNumber 1
SpeciesName
H
DiffusionBarrier
0.15

MergeEventNumber  1
MergeReaction
H  H   H2   0.8
```

#### Example 3: Same as Example 2 but with Mean Field Approximation for H Species

**Input File (`infile`) Content:**

```plaintext
Mode KMC
Temperature 1300 # unit Kelvin
Vibmode 1
Vib    10
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize
1000 1000

TotalStep           2E5

GasNumber 1
GasName
H2
GasPressure  # in Torr
0.1
GasAdsorptionRate
3570
GasAdsorbReaction
#Reactant  Product1  Product2
H2  H   H

MSpeciesNumber 1
MSpeciesName
H
InitialMSpeciesDensity
0.0001
FixMSpeciesDensity
0

MergeEventNumber  1
MergeReaction
H  H   H2   0.8
```

#### Example 4: Hydrogen Saturation at Graphene Edges Under Different Hydrogen Pressures

- **Files:** See `examples/H-no_MFA` and `examples/H-saturation`
- **Included Events:**
  1. Dissociative adsorption of H2 -> H + H
  2. Diffusion of H
  3. H + H -> H2 (associative desorption)
  4. Attachment of H at three different sites on graphene edges
  5. Detachment of H from three different sites on graphene edges

**Input File (`infile`) Content:**

```plaintext
Mode KMC
Temperature 1300 # unit Kelvin
Vibmode 1
Vib    10
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize
500 500
CoreSize
100 100

TotalStep           2E9
PrintSpecDenStep    1E7
SaveStatus          1E7

GasNumber 1
GasName
H2
GasPressure  # in Torr
0.1
GasAdsorptionRate
3570
GasAdsorbReaction
#Reactant  Product1  Product2
H2  H   H

SpeciesNumber 1
SpeciesName
H
DiffusionBarrier
0.15

AttachBarrier
0.75
AttachBarrier2
0.53
AttachBarrier3
0.84
DetachBarrier
1.74
DetachBarrier2
1.61
DetachBarrier3
1.99

MergeEventNumber  1
MergeReaction
H  H   H2   0.8
```

#### Example 5: Reaction Network of Multiple Carbon-Hydrogen Species with Mean Field Approximation for Species H

- **Files:** See `examples/concentration`
- **Included Events:**
  1. Dissociative adsorption of hydrogen gas: H2 -> H + H
  2. Dissociative adsorption of methane gas: CH4 -> CH3 + H
  3. Diffusion of species other than H
  4. Merge reactions including:
     - H+H->H2 (hydrogen recombination into gas)
     - Various combinations leading to new species or gas release
  5. Decomposition reactions for various species
  6. Desorption of C2H2

**Input File (`infile`) Content:**

```plaintext
Mode KMC
Temperature 1300 # unit Kelvin
Vibmode 1
Vib    10
MeshType 1  # 1 for hexagonal grid, 2 for square grid

MeshSize
1000 1000

TotalStep           1E10

GasNumber 2
GasName
H2 CH4
GasPressure  # in Torr
0.001 1000.0
GasAdsorptionRate
3570 0.00874
GasAdsorbReaction
#Reactant  Product1  Product2
H2  H   H
CH4 CH3 H

MSpeciesNumber 1
MSpeciesName
H

SpeciesNumber  9  #species without MFA
SpeciesName
C CH CH2 CH3 C2 C2H C2H2 C3 C3H
DiffusionBarrier
0.50 0.10 0.18 -1 0.49 0.32 0.44 -0.1 -0.1

MergeEventNumber  13
MergeReaction
#Reactant1 Reactant2 Product  Barrier
H  H   H2   0.8     1
H  C   CH   0.79    2
H  CH  CH2  0.65    3
H  CH2 CH3  0.68    4
H  CH3 CH4  0.69    5
H  C2  C2H  0.72    6
H  C2H C2H2 0.85    7
H  C3  C3H  0.57    8
C  C   C2   0.25    9
C  CH  C2H  1.27    10
C  C2  C3   1.14    11
CH CH  C2H2 0.14    12
CH C2  C3H  0.93    13

DecompEventNumber  11
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

DesorpEventNumber  1
DesorpReaction
C2H2    1.83
```

### 4. Detailed Output Explanation

DGrow outputs information both to the screen and to files. Screen output includes summaries of calculation parameters and basic KMC run information. File outputs contain the main results of the KMC simulation. Below are descriptions of each output file:

#### 4.1 SpecNum

The `SpecNum` file records the number of particles for all species at each time point.

#### 4.2 Density

The `Density` file records the particle density for all species at each time point.

#### 4.3 Statistics

The `Statistics` file records the occurrence frequency of all KMC events.

#### 4.4 EdgeInfo

The `EdgeInfo` file documents the adsorption status of particles at the edges of two-dimensional materials. This information is crucial for understanding the interaction dynamics specifically occurring at the boundaries of the growing material.

#### 4.5 Status

The `Status` file records the occupation status of particles on each grid point of both the substrate surface and the edges of the two-dimensional material. It saves all current state information, which is essential for resuming simulations from a previous state. This file allows users to restore the exact configuration of the simulation environment, including particle positions and states, ensuring continuity in computational studies.

#### 4.6 EventInfo

The `EventInfo` file contains detailed information about all kinetic Monte Carlo (KMC) events that have occurred during the simulation. This data serves as a comprehensive log for users to review and analyze the sequence of events, aiding in validation and troubleshooting. Additionally, this file is invaluable when resuming simulations, as it provides insights into past events and helps set up accurate continuation parameters. 

These output files collectively provide a thorough record of the simulation process and outcomes, supporting both immediate analysis and future reference or resumption of the simulation.
