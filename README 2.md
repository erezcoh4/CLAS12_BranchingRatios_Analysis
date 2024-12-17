# CLAS12_BranchingRatios_Analysis


Last update Dec-16, 2024
    


# **(1) General**
--------------------------------------------------------
Data analysis of CLAS12/RGB data for brnaching ratio study
  
## (1.1) Ratios of interest to be studied 
--------------------------------------------------------
We look for the following ratio 

    N(𝛾*D → π⁰pX) / N(𝛾*D → ηpX),

as a function of the momentum of the struck nucleon.

The analysis is limited to cases where the meson decays to 2 gamma-rays, i.e
    
    π⁰ → 2𝛾 and η → 2𝛾
    
[π⁰ ~ (uu¯ - dd¯)√2, η ~ (uu¯ + dd¯ - 2ss¯)√6] 



### (1.1.1) Desired event class 
--------------------------------------------------------
We are looking for a clas of events with exactly one proton and two photons in the final state:  

    D(e,e'p2𝛾X)

where the invariant mass of the two photons can be reconstructed as a meson production:
 
    For π⁰ production, D(e,e'pπ⁰X)
    For η production,  D(e,e'pηX)     
  
  


# **(2) Software**
--------------------------------------------------------

# Prequisits
---------------------------------------
### Load modules on ifarm:

module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module purge
module load sqlite/dev
module load clas12


# **(2) Comments**
--------------------------------------------------------

## (2.1) $\gamma$ classification
1. g1 - photon with higher energy
2. g2 - photon with smaller energy

## (2.2) observables and output variables

**Name** calculation                definition

**M2g**  (g1_p4 + g2_p4).M();       invariant mass of the two photons M2g = |g1_p4 + g2_p4|



# **(3) ToDo**
--------------------------------------------------------

1. We start only with particle detection in the forward detector - ?
2. Learn how to identify protons in the forward detector
3. Learn how to identify and isolate photons (g1 and g2) in the forward detector
4. Add cuts in a dedicated file cuts/RGBcutValues.csv



# **(4) ToDo**
--------------------------------------------------------
    git pull && python ./macros/execute_skim --fdebug=4 --no-email --NeventsMax=44
