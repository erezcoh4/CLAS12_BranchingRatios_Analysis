# CLAS12_BranchingRatios_Analysis


Last update Dec-26, 2024
    


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

## (2.1) Prequisits
---------------------------------------
### Load modules on ifarm:

1. module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
2. module purge
3. module load sqlite/dev
4. module load clas12


## (2.2) Execute
---------------------------------------
    git pull && python ./macros/execute_skim --fdebug=4 --no-email --NeventsMax=44


# **(3) Comments**
--------------------------------------------------------

## (3.1) $\gamma$ classification
1. g1 - photon with higher energy
2. g2 - photon with smaller energy

## (3.2) observables and output variables

**Name**        Computation formula                     definition

**status**

**runnum**

**evnum**

**e_P**

**e_Theta**

**e_Phi**

**e_Vz**

**e_DC_sector**

**e_DC_Chi2N**

**p_P**         p_p4.P()                                            proton momentum         [GeV/c]

**p_Theta**     p_p4.Theta()                                        proton polar angle      [rad.]

**p_Phi**       p_p4.Phi()                                          proton azimuthal angle  [rad.]

**p_Vz**        Vp.Z()                                              proton vertex z position [cm]                 

**p_DC_sector**                                                     proton DC sector

**p_DC_Chi2N**                                                      proton Chi2/NDF  

**g1_E**        g1_p4.E()                                           energy of the first gamma-ray

**g1_Theta**

**g1_Phi**

**g1_Vz**

**g1_DC_sector**

**g1_DC_Chi2N**

**g1_beta**     β(𝛾1)

**g2_E**        g2_p4.E()                                           energy of the second gamma-ray

**g2_Theta**

**g2_Phi**

**g2_Vz**

**g2_DC_sector**

**g2_DC_Chi2N**

**g2_E_PCAL**                                                           energy deposited in CTOF

**g2_E_ECIN**                                                           energy deposited in EC-in

**g2_E_ECOUT**                                                          energy deposited in EC-out

**g2_E_EC**     g2_E_PCAL + g2_E_ECIN + g2_E_ECOUT                      energy deposited in PCAL and ECAL

**g2_E_CTOF**                                                           energy deposited in CTOF

**g2_E_CND1**                                                           energy deposited in CND1

**g2_E_CND2**                                                           energy deposited in CND2

**g2_E_CND3**                                                           energy deposited in CND3

**g2_E_CN**     g2_E_CTOF + g2_E_CND1 + g2_E_CND2 + g2_E_CND3           energy deposited in CTOF and CND

**g2_beta**     β(𝛾1)

**Q2**          -q_p4.Mag2()                                            momentum transfer

**xB**          Q2/(2*Mp*q_p4.E())                                      Bjorken x

**omega**       q_p4.E()                                                energy transfer in the (e,e') reaction

**q**           q_p4.P()                                                momentum transfer in the (e,e') reaction

**W**           sqrt((p_rest_p4 + q_p4).Mag2());                        invariant mass  of the (e,e'p) reaction

**M2g**         (g1_p4 + g2_p4).M();                                    invariant mass of the two photons M2g = |g1_p4 + g2_p4|, in [GeV/c]

**M_x_peep**    ((q_p4 + p_rest_p4) - (p_p4)).Mag();                    missing mass of the p(e,e'p) reaction

**M_x_deep**    ((q_p4 + d_rest_p4) - (p_p4)).Mag();                    missing mass of the d(e,e'p) reaction

**M_x_deep2g**  ((q_p4 + d_rest_p4) - (p_p4 + g1_p4 + g2_p4)).Mag();    missing mass of the d(e,e'p𝛾𝛾) reaction


# **(4) ToDo**
--------------------------------------------------------
[1] Find why PCAL coordinates for the gammas (lv and lw) are so small (we normally want > 14 or so) - is the problem with them? 
[2] Is there a study in CLAS12 that looked for eta mesons? If yes - how did they detect it in the M(𝛾𝛾)?
 => Found in [2] 
[5] Find what is the ML cut and how to apply it
[4] Study the impact of a minimal Edep cut for gamma rays in the EC on the Mgg distribution
[6] Study what is the large background of M(gamma-gamma) = 0
[7] Add kinematic cuts - can we learn from previous studies in CLAS12 of diphoton and / or GlueX?
[8] Detector cuts on proton for PID refinement?
[9] Add fiducial cuts for proton as a requirement
[v] Debug the vanishing values g1_E_EC, g1_E_CN, and g2_E_CN. Done Dec-2024  
[v] Are there events in which the energy of gamma rays is smaller than the energy they deposited in the EC? - No. Verified on Dec-2024.     
[v] Correct gamma momentum to energy 
[v] Add gamma velocity, β(𝛾), which should be chosen in the range 0.9 - 1.1 following [1]
[v] Restrict particle detection to only the forward detector



# **(4) References**
--------------------------------------------------------

[1] G. Matousek "Measurements of Beam Spin Asymmetries of π±π0 dihadrons at CLAS12 using Gradient Boosted Trees", https://www.jlab.org/Hall-B/shifts/admin/paper_reviews/2024/PiPi0_AnalysisNote_(3)_3.pdf-7600066-2024-03-09-v1.pdf

[2] Celentano "Analysis of the reaction γp → pπ0η with the g12 dataset and extraction of the γp → pa2(1320) differential cross-section", https://www.jlab.org/Hall-B/shifts/admin/paper_reviews/2019/Analysis_of_the_reaction_gamma_p____p_pi0_eta_with_the_g12_dataset_3-7164322-2019-10-25-v8.pdf
