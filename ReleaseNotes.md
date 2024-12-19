

# **Dec-19, 2024**
## git commit 260d307 
--------------------------------------------------------
1. Added Electromagnetic Calorimeters (ECAL) features of the two gamma-rays: g_E_PCAL, g_E_ECIN, g_E_ECOUT.
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
  
