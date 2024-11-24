// last edit Nov-24, 2023

#include "Auxiliary/BranchingRatio_CLAS12_Auxiliary.h"

#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>+;
#endif

using namespace clas12;

// Results in CSV file (e,e'p2𝛾X)
TString csvheader = ( (TString)"status,runnum,evnum,beam_helicity,"
                     // e'
                     +(TString)"e_P,e_Theta,e_Phi,e_Vz,e_DC_sector,"
                     // p
                     +(TString)"p_P,p_Theta,p_Phi,p_Vz,p_DC_sector,"
                     // two photons
                     +(TString)"g1_E,g1_Theta,g1_Phi,g1_Vz,g1_DC_sector,"
                     +(TString)"g2_E,g2_Theta,g2_Phi,g2_Vz,g2_DC_sector,"
                     // kinematics
                     +(TString)"Q2,xB,omega,y,W,M_x,"
                     );

std::vector<int> csvprecisions = {
    0,0,0,0,
    9,9,9,9,0,
    9,9,9,9,0,
    9,9,9,9,0,
    9,9,9,9,0,
    9,9,9,9,9,9    
};

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// globals
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
int fdeubg = 2;



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void Debug(int _fdebug_, const char* fmt, ...) {
    va_list arg;
    va_start(arg, fmt);

    if (fdebug > _fdebug_) {
        vprintf(fmt, arg);
        std::cout << std::endl;
    }
    //vprintf can be replaced with vsprintf (for sprintf behavior)
    //or any other printf function preceded by a v
    va_end(arg);
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// main
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void c12rSkimmer_BranchingRatios(int RunNumber    = 6420)
{
    Debug(1, "Begin main");
    Debug(1, "Done.");
}
