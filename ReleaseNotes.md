

## **Jan-1, 2025**
## git commit 
-------------------------------------------------------- 
1. Corrected an error in M_x_deep2g (it accidently counted g1 twice)


## **Dec-27, 2024**
## git commit 1fb8ef9  
-------------------------------------------------------- 
1. retained PCAL energy deposit: g1_E_PCAL and g2_E_PCAL (following Matousek in [1])
2. Removed g1_E_EC and g2_E_EC
3. Corrected gamma momentum to energy
4. Added PCAL coordinates g_PCAL_W and g_PCAL_V to add a cut on them (following Matousek in [1])
5. Added gamma velocity, β(𝛾) - g1_beta and g2_beta, which should be chosen in the range 0.9 - 1.1 following [1] 
6. Corrected a bug with gamma PCPAL coordinates V,W


# **Dec-20, 2024**
## git commit b85822a
--------------------------------------------------------
1. Added Electromagnetic Calorimeters (ECAL) and CN (CTOF, CND) features of the two gamma-rays: g_E_PCAL, g_E_ECIN, g_E_ECOUT, g_E_CTOF, g_E_CND...
2. Added kinematic variables:
    2.a. M_x_peep - invariant mass of the p(e,e'p) reaction       |(q_p4 + p_rest_p4) - (p_p4)|
    2.b. M_x_deep - invariant mass of the d(e,e'p) reaction       |(q_p4 + d_rest_p4) - (p_p4)|
    2.c. M_x_deep2g - invariant mass of the d(e,e'p𝛾𝛾) reaction   |(q_p4 + d_rest_p4) - (p_p4 + g1_p4 + g2_p4)|


# **Dec-17, 2024**
## git commit e1d78c 
--------------------------------------------------------
1. Added a calculation of Mx for the (e,e'p2g)X event
2. Added extraction of the data for the two gamma-rays (energy and momomentum, and vertex)
3. Added a cut on the vertex z-difference between the gamma rays and the electron
4. Solved problem in ExtractElectronInformation(), it was due to a missing definition of int DC_layers[3] = {6,18,36};


# **Nov-24, 2024**
## git commit fb5ffa
--------------------------------------------------------
Created repo, and added some content 
  


# **References**
--------------------------------------------------------
[1] G. Matousek Measurements of Beam Spin Asymmetries of π±π0 dihadrons at CLAS12 using Gradient Boosted Trees https://www.jlab.org/Hall-B/shifts/admin/paper_reviews/2024/PiPi0_AnalysisNote_(3)_3.pdf-7600066-2024-03-09-v1.pdf

