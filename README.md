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

1. module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
2. module purge
3. module load sqlite/dev
4. module load clas12


# **(2) Comments**
--------------------------------------------------------

## (2.1) $\gamma$ classification
1. g1 - photon with higher energy
2. g2 - photon with smaller energy

## (2.2) observables and output variables

**Name** calculation                definition


**status**
**runnum**
**evnum**
**e_P**
**e_Theta**
**e_Phi**
**e_Vz**
**e_DC_sector**
**e_DC_Chi2N**
**p_P**         p_p4.P()                proton momentum         [GeV/c]
**p_Theta**     p_p4.Theta()            proton polar angle      [rad.]
**p_Phi**       p_p4.Phi()              proton azimuthal angle  [rad.]
**p_Vz**        Vp.Z()                  proton vertex z position [cm]                 
**p_DC_sector**                         proton DC sector
**p_DC_Chi2N**                          proton Chi2/NDF  
**g1_E**        g1_p4.E()               energy of the first gamma-ray
**g1_Theta**
**g1_Phi**
**g1_Vz**
**g1_DC_sector**
**g1_DC_Chi2N**
**g2_E**        g2_p4.E()               energy of the second gamma-ray
**g2_Theta**
**g2_Phi**
**g2_Vz**
**g2_DC_sector**
**g2_DC_Chi2N**
**Q2**          -q_p4.Mag2()                                            momentum transfer
**xB**          Q2/(2*Mp*q_p4.E())                                      Bjorken x
**omega**       q_p4.E()                                                energy transfer in the (e,e') reaction
**q**           q_p4.P()                                                momentum transfer in the (e,e') reaction
**W**           sqrt((p_rest_p4 + q_p4).Mag2());                        invariant mass  of the (e,e'p) reaction
**M2g**         (g1_p4 + g2_p4).M();                                    invariant mass of the two photons M2g = |g1_p4 + g2_p4|, in [GeV/c]
**M_x_peep**    ((q_p4 + p_rest_p4) - (p_p4)).Mag();                    missing mass of the p(e,e'p) reaction
**M_x_deep**    ((q_p4 + d_rest_p4) - (p_p4)).Mag();                    missing mass of the d(e,e'p) reaction
**M_x_deep2g**  ((q_p4 + d_rest_p4) - (p_p4 + g1_p4 + g2_p4)).Mag();    missing mass of the d(e,e'p𝛾𝛾) reaction

# **(3) ToDo**
--------------------------------------------------------

1. Restrict particle detection to only the forward detector
2. Add kinematic cuts - can we learn from what has been done in GlueX?
3. Detector cuts on proton for PID refinement?
4. Add fiducial cuts for proton
5. Study what is the large background of M(gamma-gamma) = 0



# **(4) ToDo**
--------------------------------------------------------
    git pull && python ./macros/execute_skim --fdebug=4 --no-email --NeventsMax=44
