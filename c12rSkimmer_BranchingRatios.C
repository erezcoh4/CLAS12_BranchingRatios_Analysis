// last edit Dec-17, 2024

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include "clas12reader.h"
#include "Auxiliary/csv_reader.h"
#include "Auxiliary/BranchingRatio_CLAS12_Auxiliary.cpp"
#include "Auxiliary/DCfid_SIDIS.cpp"

#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>+;
#endif

using namespace clas12;
BranchingRatio_CLAS12_Auxiliary aux;
auto db = TDatabasePDG::Instance();
DCfid_SIDIS dcfid;

//std::pair<TString,int> csv_var_precision;

// Results in CSV file d(e,e'p2𝛾X)
TString csvheader = ( (TString)"status,runnum,evnum,"
                     +(TString)"e_P,e_Theta,e_Phi,e_Vz,e_DC_sector,e_DC_Chi2N,"         // e
                     +(TString)"p_P,p_Theta,p_Phi,p_Vz,p_DC_sector,p_DC_Chi2N,"         // p
                     +(TString)"g1_E,g1_Theta,g1_Phi,g1_Vz,g1_DC_sector,g1_DC_Chi2N,"   // photon-1
                     +(TString)"g2_E,g2_Theta,g2_Phi,g2_Vz,g2_DC_sector,g2_DC_Chi2N,"   // photon-2
                     +(TString)"Q2,xB,omega,W,M_x,q,M2g,"                               // kinematics
                     );

std::vector<int> csvprecisions = {
    0,0,0,
    4,4,4,4,0,4,
    4,4,4,4,0,4,
    4,4,4,4,0,4,
    4,4,4,4,0,4,
    4,4,4,4,4,4,4
};

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// globals
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
int verbosity = 2;

// 4-vectors for the reaction d(e,e'p2𝛾X)
TLorentzVector              p_rest_p4;
TLorentzVector     Beam_p4, target_p4;
TLorentzVector       e_p4, q_p4, p_p4;
TLorentzVector           g1_p4, g2_p4; // gamma 1 and gamma 2
TLorentzVector*    reco_pi0_p4 = NULL; // best fit pi0
TLorentzVector*    reco_eta_p4 = NULL; // best fit eta
std::vector<region_part_ptr>  electrons, protons, gammas;
int                Ne, Np, Ngammas;
int          Nevents_processed = 0;
int                   evnum, runnum;
int               torusBending = -1; // -1 for In-bending, +1 for Out-bending
int bending  = 1 ? (torusBending==-1) : 0; // bending: 0(out)/1(in)

float                         Ebeam;
TString               Skimming = "";
TString                 prefix = "";
TString               DataPath = "";
TString infilename, outfilepath, outfilename;
TString    full_outcsvfilename = "";
ofstream          outcsvfile_eep2gX;

// detector features
int               DC_layer, status;
int               DC_layers[3] = {6,18,36};// Region 1 is denoted at DC detector 6, Region 2 is denoted 18, Region 3 - as 36
int      Nevents_passed_e_cuts = 0;
int      Nevents_passed_p_cuts = 0;
int     Nevents_passed_g1_cuts = 0;
int     Nevents_passed_g2_cuts = 0;
int    Nevents_passed_eep_cuts = 0;
int  Nevents_passed_eep2g_cuts = 0;

// leading electron
// electron energy deposit in PCAL [GeV], in ECAL_in [GeV], in ECAL_out [GeV]...
double e_E_PCAL, e_E_ECIN, e_E_ECOUT, e_PCAL_W, e_PCAL_V, e_PCAL_x, e_PCAL_y, e_PCAL_z;
double e_PCAL_sector, e_DC_sector;
double e_DC_Chi2N, e_DC_x[3], e_DC_y[3], e_DC_z[3];
TVector3 Ve;
bool     ePastCutsInEvent;

// proton
double p_E_PCAL, p_E_ECIN, p_E_ECOUT, p_PCAL_W, p_PCAL_V, p_PCAL_x, p_PCAL_y, p_PCAL_z;
double p_PCAL_sector, p_DC_sector, p_DC_Chi2N, p_DC_x[3], p_DC_y[3], p_DC_z[3];
TVector3 Vp;
bool     pPastCutsInEvent;

// two photons
double g1_E_PCAL, g1_E_ECIN, g1_E_ECOUT, g1_PCAL_W, g1_PCAL_V, g1_PCAL_x, g1_PCAL_y, g1_PCAL_z;
double g1_PCAL_sector, g1_DC_sector, g1_DC_Chi2N, g1_DC_x[3], g1_DC_y[3], g1_DC_z[3];
TVector3 Vg1;
bool     g1PastCutsInEvent;

double g2_E_PCAL, g2_E_ECIN, g2_E_ECOUT, g2_PCAL_W, g2_PCAL_V, g2_PCAL_x, g2_PCAL_y, g2_PCAL_z;
double g2_PCAL_sector, g2_DC_sector, g2_DC_Chi2N, g2_DC_x[3], g2_DC_y[3], g2_DC_z[3];
TVector3 Vg2;
bool     g2PastCutsInEvent;


// kinematics
double xB, Q2, omega, W, M_x;
double M2g; // invariant mass of the two photons M2g = |g1_p4 + g2_p4|
double Pe_phi, Pe_theta;
double q_phi,   q_theta;
double Pp_phi, Pp_theta;
bool    eepPastCutsInEvent;


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// Routines
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void DEBUG(int v, const char* fmt, ...) {
    va_list arg;
    va_start(arg, fmt);
    //    std::cout << "Debug(), verbosity: " << verbosity << ", v: " << v  << ", " << fmt << std::endl;
    if (verbosity > v) {
        vprintf(fmt, arg);
        std::cout << std::endl;
    }
    //vprintf can be replaced with vsprintf (for sprintf behavior)
    //or any other printf function preceded by a v
    va_end(arg);
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void ConfrontValueWithCut(TString varlabel, double var, double cutValue){
    DEBUG(3,"%s: %.2f (cut value %.2f)",varlabel.Data(),var, cutValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TString GetRunNumberSTR( int RunNumber, TString fSkimming ){
    char RunNumberStr[20];
    // sprintf( RunNumberStr, "00%d", RunNumber );
    
    if(fSkimming == "p_uniform_distribution"){
        // "white" GEMC simulation runs
        sprintf( RunNumberStr, "%d", RunNumber );
    } else {
        sprintf( RunNumberStr, "%06d", RunNumber );
    }
    if (verbosity>1) std::cout << "(SIDIS) skimming run " << RunNumberStr << std::endl;
    return (TString)RunNumberStr;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
TVector3 GetParticleVertex(clas12::region_part_ptr rp){
    TVector3 V(rp->par()->getVx(),
               rp->par()->getVy(),
               rp->par()->getVz());
    return V;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetLorentzVector (TLorentzVector &p4,clas12::region_part_ptr rp){
    p4.SetXYZM(rp->par()->getPx(),
               rp->par()->getPy(),
               rp->par()->getPz(),
               p4.M());
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetDataPath (TString fDataPath, Double_t fEbeam) {
    DEBUG(2,"SetDataPath('%s',%f)",fDataPath.Data(),fEbeam);
    prefix   = "sidisdvcs_"; // default
    
    if (fDataPath=="" || fDataPath=="sidisdvcs" || fDataPath=="sidis dvcs"){
        
        // sidis-dvcs train files, used since July 2022
        // (the 'usual' train files)
        if        (fabs(fEbeam-10.2)<0.01){
            DataPath = "/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/";
        } else if (fabs(fEbeam-10.4)<0.01){
            DataPath = "/cache/clas12/rg-b/production/recon/spring2020/torus-1/pass1/v1/dst/train/sidisdvcs/";
        } else if (fabs(fEbeam-10.6)<0.01){
            DataPath = "/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/";
        }
        prefix   = "sidisdvcs_";
    }
    else if (fDataPath=="inclusive" || fDataPath=="inc"){
        // inclusive train files, used until July 2022
        // (inclusive train files were only generated in the beginning of RGB without any backup)
        DataPath = "/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train_20200610/inc/";
        prefix   = "inc_";
    }
    else if (fDataPath=="nSidis" || fDataPath=="nsidis"){
        // free-p data from RGA data
        // For RGA we use nSidis, they key difference is sidisdvcs has e_p > 1 GeV and nSidis has e_p > 2 GeV.
        DataPath = "/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass1/v1/dst/train/nSidis/";
        prefix   = "nSidis_";
    }
    else if (fDataPath=="AcceptanceCorrection"){
        // GEMC simulations of "white" spectra
        // i.e. (e,e'π) events with no physics generator
        DataPath = "/volatile/clas12/users/ecohen/GEMC/hipo/10.2/AcceptanceCorrection/";
        prefix = "p_uniform_distribution";
    }
    
    DEBUG(2,DataPath);
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetVerbosity( int v ){
    verbosity = v;
    aux.SetVerbosity(v);
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetEbeam (double fEbeam=10.2) { // [GeV]
    // RGA the enrgy was 10.6
    // RGB Spring-2019 the enrgy was 10.2
    // RGB Fall-2019 the enrgy was 10.4096
    Ebeam = fEbeam;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetTarget (){
    p_rest_p4.SetXYZM (0, 0, 0, aux.Mp);
    target_p4.SetXYZM (0, 0, 0, aux.Md);
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetGlobals(float fEbeam=10.2, TString fDataPath = "sidisdvcs") {
    SetEbeam            (fEbeam);
    SetTarget           ();
    SetDataPath         (fDataPath, fEbeam);
    //    SetSkimming         ( fSkimming  );
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetFileNames(int RunNumber=6164) {
    TString RunNumberStr = GetRunNumberSTR(RunNumber, Skimming);
    // define input filename
    
    infilename  = DataPath + prefix + RunNumberStr + ".hipo";
    outfilepath = "/volatile/clas12/users/ecohen/RGB/" + Skimming + "/";
    outfilename = "skimmed_BranchingRatios_" + prefix + RunNumberStr;
    
    if (verbosity>1){
        std::cout
        << "Input file name: " << std::endl
        << infilename          << std::endl
        << "Output file name: "<< std::endl
        << outfilepath + outfilename
        << std::endl;
    }
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void LoadCutValues(TString cutFilename = "cuts/cutValues_16Dec2024.csv") {
    // read cut values
    std::string cutFilename_string(cutFilename.Data());
    aux.loadCutValues(cutFilename_string,torusBending);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GetParticlesByType (){
    // get particles by type
    Ne      = electrons .size();
    Np      = protons   .size();
    Ngammas = gammas    .size();
    DEBUG(1,"N(e):%d, N(p):%d, N(g):%d ",Ne,Np,Ngammas);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpenResultFiles(){
    
    full_outcsvfilename = outfilepath + outfilename + "_eep2gX.csv";
    outcsvfile_eep2gX.open( full_outcsvfilename );
    outcsvfile_eep2gX << csvheader << std::endl;
    
    if (verbosity>1) std::cout << "Done OpenOutputFiles( " << outfilename << ")" << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeVariables(){
    DEBUG(5, "InitializeVariables()");
    Beam_p4    .SetPxPyPzE (0, 0, Ebeam, Ebeam );

    
    electrons   .clear();
    protons     .clear();
    gammas      .clear();
    
    DC_layer    = -9999;
    status      = 1; // 0 is good...
    
    DEBUG(4, "Initialize electron...");
    e_p4 = TLorentzVector(0,0,0,aux.Me);
    xB  = Q2  = omega     = -9999;
    W   = M_x             = -9999;
    e_E_ECIN    = e_E_ECOUT = e_E_PCAL  = -9999;
    e_PCAL_W    = e_PCAL_V              = -9999;
    e_PCAL_x    = e_PCAL_y  = e_PCAL_z  = -9999;
    e_PCAL_sector                       = -9999;
    e_DC_sector = e_DC_Chi2N            = -9999;
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        e_DC_x[regionIdx]               = -9999;
        e_DC_y[regionIdx]               = -9999;
        e_DC_z[regionIdx]               = -9999;
    }
    Pe_phi = q_phi = q_theta            = 0;
    Ve                                  = TVector3();
    ePastCutsInEvent                    = false;
    
    DEBUG(4, "Initialize proton...");
    p_p4 = TLorentzVector(0,0,0,aux.Mp);
    p_E_ECIN    = p_E_ECOUT = p_E_PCAL  = -9999;
    p_PCAL_W    = p_PCAL_V              = -9999;
    p_PCAL_x    = p_PCAL_y  = p_PCAL_z  = -9999;
    p_PCAL_sector                       = -9999;
    p_DC_sector = p_DC_Chi2N            = -9999;
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        p_DC_x[regionIdx]               = -9999;
        p_DC_y[regionIdx]               = -9999;
        p_DC_z[regionIdx]               = -9999;
    }
    Pp_phi                              = 0;
    Vp                                  = TVector3();
    pPastCutsInEvent                    = false;
    
    DEBUG(5, "Initialize (e,e'p)");
    eepPastCutsInEvent                  = false;
    
    DEBUG(5, "Done InitializeVariables()");
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfElectronPassedSelectionCuts(){
    DEBUG(3,"CheckIfElectronPassedSelectionCuts()");
    // decide if electron in event passes event selection cuts
    
    // DC - fiducial cuts on DC
    // from bandsoft_tools/skimmers/electrons.cpp,
    // where eHit.getDC_x1() - x position in first region of the drift chamber
    // same for y1,x2,y2,...
    // eHit.getDC_sector() - sector
    // checking DC Fiducials
    // torusBending         torus magnet bending:   ( 1 = inbeding, -1 = outbending    )
    
    // sometimes the readout-sector is 0. This is funny
    // Justin B. Estee (June-21): I also had this issue. I am throwing away sector 0. The way you check is plot the (x,y) coordinates of the sector and you will not see any thing. Double check me but I think it is 0.
    DEBUG(3,"electron DC sector: %.0f, bending: %d",e_DC_sector, bending);
    if (e_DC_sector == 0) return false;
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        // DC_e_fid: sector 1-6, layer 1-3
        DEBUG(3,"\t(x=%.1f,y=%.1f), sector: %.0f",e_DC_x[regionIdx], e_DC_y[regionIdx],e_DC_sector);
        bool DC_fid  = dcfid.DC_fid_xy_sidis(11,                 // particle PID,
                                             e_DC_x[regionIdx],  // x
                                             e_DC_y[regionIdx],  // y
                                             e_DC_sector,        // sector
                                             regionIdx+1,        // layer
                                             bending);           // torus bending
        DEBUG(3,"\tDC fid (region %d): %d",regionIdx, DC_fid);
        if (DC_fid == false) {
            return false;
        }
    }
    
    ConfrontValueWithCut("e PCAL(W)",e_PCAL_W,aux.cutValue_e_PCAL_W);
    ConfrontValueWithCut("e PCAL(V)",e_PCAL_V,aux.cutValue_e_PCAL_V);
    ConfrontValueWithCut("e E-PCAL",e_E_PCAL,aux.cutValue_e_E_PCAL);
    ConfrontValueWithCut("Sampling Fraction minimum",(e_E_PCAL + e_E_ECIN + e_E_ECOUT)/e_p4.P(), aux.cutValue_SamplingFraction_min);
    ConfrontValueWithCut("e PCAL ECIN SF min",e_E_ECIN/e_p4.P(),aux.cutValue_PCAL_ECIN_SF_min - e_E_PCAL/e_p4.P());
    ConfrontValueWithCut("V(e)-z", Ve.Z(), aux.cutValue_Vz_min );
    
    if(!(true
         // fiducial cuts on PCAL
         //fabs(e_PCAL_x)>0
         //&&  fabs(e_PCAL_y)>0
         &&  e_PCAL_W > aux.cutValue_e_PCAL_W
         &&  e_PCAL_V > aux.cutValue_e_PCAL_V
         
         // Electron Identification Refinement  - PCAL Minimum Energy Deposition Cut
         &&  e_E_PCAL > aux.cutValue_e_E_PCAL
         
         // Sampling fraction cut
         && ((e_E_PCAL + e_E_ECIN + e_E_ECOUT)/e_p4.P()) > aux.cutValue_SamplingFraction_min
         && (e_E_ECIN/e_p4.P() > aux.cutValue_PCAL_ECIN_SF_min - e_E_PCAL/e_p4.P()) // RGA AN puts "<" here mistakenly
         
         // Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
         // Spring 19 and Spring 2020 in-bending.
         // Fall 2019 (without low-energy-run) was out-bending.
         &&  ((aux.cutValue_Vz_min < Ve.Z()) && (Ve.Z() < aux.cutValue_Vz_max))
         )) return false;
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractElectronInformation(){
    DEBUG(3,"ExtractElectronInformation()");
    // ------------------------------------------------------------------------------------------------
    // extract electron information
    // ------------------------------------------------------------------------------------------------
    // find leading electron as the one with highest energy
    if (Ne == 0) return;
    double  leading_e_E;
    int     leading_e_index = 0;
    SetLorentzVector(e_p4,electrons[0]);
    TLorentzVector e_tmp(0,0,0,aux.Me);
    for (int eIdx=0; eIdx < Ne; eIdx++) {
        SetLorentzVector(e_tmp  ,electrons[eIdx]);
        double Ee = e_tmp.E();
        if (Ee > leading_e_E) {
            leading_e_index = eIdx;
            leading_e_E     = Ee;
        }
    }
    DEBUG(3,"lead index: %d",leading_e_index);
    
    // set leading electron 4-momentum
    SetLorentzVector(e_p4 , electrons[leading_e_index]);
    // set leading electron vertex
    Ve              = GetParticleVertex( electrons[leading_e_index] );
    
    // detector information on electron
    auto e_PCAL_info= electrons[leading_e_index]->cal(PCAL);
    e_E_PCAL        = e_PCAL_info->getEnergy();
    e_PCAL_sector   = e_PCAL_info->getSector();
    e_PCAL_V        = e_PCAL_info->getLv();
    e_PCAL_W        = e_PCAL_info->getLw();
    e_E_ECIN        = electrons[leading_e_index]->cal(ECIN)->getEnergy();
    e_E_ECOUT       = electrons[leading_e_index]->cal(ECOUT)->getEnergy();
    DEBUG(3,"PCAL sector %.0f, (V=%.1f,W=%.1f)",e_PCAL_sector, e_PCAL_V,e_PCAL_W);
    
    // hit position in PCAL
    e_PCAL_x        = e_PCAL_info->getX();
    e_PCAL_y        = e_PCAL_info->getY();
    e_PCAL_z        = e_PCAL_info->getZ();
    
    // Drift Chamber tracking system
    auto e_DC_info  = electrons[leading_e_index]->trk(DC);
    e_DC_sector     = e_DC_info->getSector(); // tracking sector
    e_DC_Chi2N      = e_DC_info->getChi2N();  // tracking chi^2/NDF
    DEBUG(3,"DC: sector %.0f, Chi2/N %.1f",e_DC_sector,e_DC_Chi2N);
    
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        int DC_layer = DC_layers[regionIdx];
        e_DC_x[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getX();
        e_DC_y[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getY();
        e_DC_z[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getZ();
        DEBUG(3,"\tRegion %d, DC_layer %d: (x=%.3f,y=%.3f)",regionIdx, DC_layer, e_DC_x[regionIdx], e_DC_y[regionIdx]);
    }
    
    
    DEBUG(2,"Extracted electron information");
    // ------------------------------------------------------------------------------------------------
    // now, check if electron passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    ePastCutsInEvent = CheckIfElectronPassedSelectionCuts();
    if (ePastCutsInEvent) {DEBUG(2, "** electron succesfully past cuts **"); Nevents_passed_e_cuts++ ;}
    else                  {DEBUG(2, "** electron did not pass cuts succesfully **");}
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfProtonPassedSelectionCuts(){
    DEBUG(3,"CheckIfProtonPassedSelectionCuts()");
    
    DEBUG(3,"proton DC sector: %.0f, bending: %d",p_DC_sector, bending);
    if (p_DC_sector == 0) return false;
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        DEBUG(3,"\t(x=%.1f,y=%.1f), sector: %.0f",p_DC_x[regionIdx], p_DC_y[regionIdx], p_DC_sector);
        bool DC_fid  = dcfid.DC_fid_xy_sidis(2212,               // particle PID,
                                             p_DC_x[regionIdx],  // x
                                             p_DC_y[regionIdx],  // y
                                             p_DC_sector,        // sector
                                             regionIdx+1,        // layer
                                             bending);           // torus bending
        DEBUG(3,"\tDC fid (region %d): %d",regionIdx, DC_fid);
        if (DC_fid == false) {
            return false;
        }
    }
    
    //    ConfrontValueWithCut("p PCAL(W)",p_PCAL_W,aux.cutValue_p_PCAL_W);
    //    ConfrontValueWithCut("p PCAL(V)",p_PCAL_V,aux.cutValue_p_PCAL_V);
    //    ConfrontValueWithCut("p E-PCAL",p_E_PCAL,aux.cutValue_p_E_PCAL);
    //    ConfrontValueWithCut("Sampling Fraction minimum",(p_E_PCAL + p_E_ECIN + p_E_ECOUT)/p_p4.P(), aux.cutValue_SamplingFraction_min);
    //    ConfrontValueWithCut("p PCAL ECIN SF min",p_E_ECIN/p_p4.P(),aux.cutValue_PCAL_ECIN_SF_min - p_E_PCAL/p_p4.P());
    //    ConfrontValueWithCut("V(p)-z", Vp.Z(), aux.cutValue_Vz_min );
    ConfrontValueWithCut("|Ve(z) - Vp(z)|", fabs((Ve-Vp).Z()), aux.cutValue_Ve_Vp_dz_max );
    
    if(!(true
         // Cut on the z-Vertex Difference Between Electrons and Hadrons
         &&  ( fabs((Ve-Vp).Z()) < aux.cutValue_Ve_Vp_dz_max )
         )) return false;
    
    return true;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfGammaPassedSelectionCuts(TVector3 Vg){
    DEBUG(3,"CheckIfGammaPassedSelectionCuts()");
        
    ConfrontValueWithCut("|Ve(z) - Vg(z)|", fabs((Ve-Vg).Z()), aux.cutValue_Ve_Vg_dz_max );
    
    if(!(true
         // Cut on the z-Vertex Difference Between Electrons and Hadrons
         &&  ( fabs((Ve-Vg).Z()) < aux.cutValue_Ve_Vg_dz_max )
         )) return false;
    
    return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractProtonInformation(){
    // ------------------------------------------------------------------------------------------------
    // extract proton info
    // ------------------------------------------------------------------------------------------------
    // find leading proton as the one with highest energy
    if (Np == 0) return;
    double  leading_p_E;
    int     leading_p_index = 0;
    SetLorentzVector(p_p4,protons[0]);
    TLorentzVector p_tmp(0,0,0,aux.Mp);
    for (int pIdx=0; pIdx < Np; pIdx++) {
        SetLorentzVector(p_tmp  ,protons[pIdx]);
        double Ep = p_tmp.E();
        if (Ep > leading_p_E) {
            leading_p_index = pIdx;
            leading_p_E     = Ep;
        }
    }
    // set leading proton 4-momentum
    SetLorentzVector(p_p4 , protons[leading_p_index]);
    // set leading proton vertex
    Vp              = GetParticleVertex( protons[leading_p_index] );
    
    // detector information on proton
    auto p_PCAL_info= protons[leading_p_index]->cal(PCAL);
    p_E_PCAL        = p_PCAL_info->getEnergy();
    p_PCAL_sector   = p_PCAL_info->getSector();
    p_PCAL_V        = p_PCAL_info->getLv();
    p_PCAL_W        = p_PCAL_info->getLw();
    p_E_ECIN        = protons[leading_p_index]->cal(ECIN)->getEnergy();
    p_E_ECOUT       = protons[leading_p_index]->cal(ECOUT)->getEnergy();
    
    // hit position in PCAL
    p_PCAL_x        = p_PCAL_info->getX();
    p_PCAL_y        = p_PCAL_info->getY();
    p_PCAL_z        = p_PCAL_info->getZ();
    
    // Drift Chamber tracking system
    auto p_DC_info  = protons[leading_p_index]->trk(DC);
    p_DC_sector     = p_DC_info->getSector(); // tracking sector
    p_DC_Chi2N      = p_DC_info->getChi2N();  // tracking chi^2/NDF
    
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        int DC_layer = DC_layers[regionIdx];
        p_DC_x[regionIdx] = protons[leading_p_index]->traj(DC,DC_layer)->getX();
        p_DC_y[regionIdx] = protons[leading_p_index]->traj(DC,DC_layer)->getY();
        p_DC_z[regionIdx] = protons[leading_p_index]->traj(DC,DC_layer)->getZ();
    }
    DEBUG(2,"Extracted proton information");
    
    // ------------------------------------------------------------------------------------------------
    // now, check if proton passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    pPastCutsInEvent = CheckIfProtonPassedSelectionCuts();
    if (pPastCutsInEvent) {DEBUG(2, "** proton succesfully past cuts **"); Nevents_passed_p_cuts++ ;}
    else                  {DEBUG(2, "** proton did not pass cuts succesfully **");}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractGammasInformation(){
    // ------------------------------------------------------------------------------------------------
    // extract photon info
    // ------------------------------------------------------------------------------------------------
    // Gamma classification:
    // we are restricting our analysis to events with one proton and two detected photons,
    // where we note
    // g1 - photon with higher energy
    // g2 - photon with smaller energy
    DEBUG(3, "ExtractGammasInformation()");
    if (Ngammas != 2) return;
    SetLorentzVector(g1_p4,  gammas[0]);
    Vg1 = GetParticleVertex( gammas[0] );
    SetLorentzVector(g2_p4,  gammas[1]);
    Vg2 = GetParticleVertex( gammas[1] );
    
    DEBUG(5, "g1_p4.E(): %.3f GeV, g2_p4.E(): %.3f GeV, Vg1.Z(): %.3f cm, Vg2.Z(): %.3f cm",g1_p4.E(),g2_p4.E(),Vg1.Z(),Vg2.Z());
    
    auto g1_DC_info  = gammas[0]->trk(DC);
    g1_DC_Chi2N      = g1_DC_info->getChi2N();  // tracking chi^2/NDF
    auto g2_DC_info  = gammas[1]->trk(DC);
    g2_DC_Chi2N      = g2_DC_info->getChi2N();  // tracking chi^2/NDF
    
    
    TLorentzVector g_tmp(0,0,0,0);
    TVector3      Vg_tmp  (0,0,0);
    double         g_tmp_DC_Chi2N;
    
    if (g2_p4.E() > g1_p4.E()) {
        g_tmp  = g1_p4;
        Vg_tmp = Vg1;
        g_tmp_DC_Chi2N = g1_DC_Chi2N;
        
        g1_p4 = g2_p4;
        Vg1   = Vg2;
        g1_DC_Chi2N = g2_DC_Chi2N;
        
        g2_p4 = g_tmp;
        Vg2   = Vg1;
        g2_DC_Chi2N = g_tmp_DC_Chi2N;
        DEBUG(5, "after swap g1_p4.E(): %.3f GeV, g2_p4.E(): %.3f GeV, Vg1.Z(): %.3f cm, Vg2.Z(): %.3f cm",g1_p4.E(),g2_p4.E(),Vg1.Z(),Vg2.Z());
    }
    DEBUG(2,"Extracted gamma information");
    
    g1PastCutsInEvent = CheckIfGammaPassedSelectionCuts(Vg1);
    if (g1PastCutsInEvent){DEBUG(2, "** gamma-1 succesfully past cuts **"); Nevents_passed_g1_cuts++ ;}
    else                  {DEBUG(2, "** gamma-1 did not pass cuts succesfully **");}
    
    g2PastCutsInEvent = CheckIfGammaPassedSelectionCuts(Vg2);
    if (g2PastCutsInEvent){DEBUG(2, "** gamma-2 succesfully past cuts **"); Nevents_passed_g2_cuts++ ;}
    else                  {DEBUG(2, "** gamma-2 did not pass cuts succesfully **");}
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ComputeKinematics(){
    // compute event kinematics
    DEBUG(3, "ComputeKinematics()");
    q_p4    = Beam_p4 - e_p4;
    Q2      = -q_p4.Mag2();
    omega   = q_p4.E();
    xB      = Q2/(2. * aux.Mp * q_p4.E());
    W       = sqrt((p_rest_p4 + q_p4).Mag2());
    M_x     = ( (q_p4 + p_rest_p4) - (p_p4 + g1_p4 + g2_p4) ).Mag(); // Mx_eep2gX
    M2g     = (g1_p4 + g2_p4).M();
    DEBUG(5, "q: %.1f, omega: %.1f, Q2: %.1f",q_p4.P(), omega, Q2);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrintVariables(){
    std::cout <<
    "run "          << runnum                   << ", "
    "event "        << evnum                    << ", "
    << std::endl << "electron: " << std::endl   <<
    "p: "           << e_p4.P()                 << " GeV/c, "
    "𝜃: "           << e_p4.Theta()*180./3.14   << " deg., "
    "ϕ: "           << e_p4.Phi()*180./3.14     << " deg., "
    "V(z) "         << Ve.Z()                   << " cm, "
    "DC-sector: "   << e_DC_sector              << ", "
    "χ2/NDF "       << e_DC_Chi2N               << ", "
    << std::endl << "proton: " << std::endl     <<
    "p: "           << p_p4.P()                 << " GeV/c ,"
    "𝜃: "           << p_p4.Theta()*180./3.14   << " deg., "
    "ϕ: "           << p_p4.Phi()*180./3.14     << " deg., "
    "V(z) "         << Vp.Z()                   << " cm, "
    "DC-sector: "   << p_DC_sector              << ", "
    "χ2/NDF "       << p_DC_Chi2N               << ", "
    << std::endl << "g1: " << std::endl         <<
    "p: "           << g1_p4.P()                 << " GeV/c, "
    "𝜃: "           << g1_p4.Theta()*180./3.14   << " deg., "
    "ϕ: "           << g1_p4.Phi()*180./3.14     << " deg., "
    "V(z) "         << Vg1.Z()                   << " cm, "
    "DC-sector: "   << g1_DC_sector              << ", "
    "χ2/NDF "       << g1_DC_Chi2N               << ", "
    << std::endl << "g2: " << std::endl         <<
    "p: "           << g2_p4.P()                 << " GeV/c,"
    "𝜃: "           << g2_p4.Theta()*180./3.14   << " deg.,"
    "ϕ: "           << g2_p4.Phi()*180./3.14     << " deg.,"
    "V(z) "         << Vg2.Z()                   << " cm,"
    "DC-sector: "   << g2_DC_sector              << ","
    "χ2/NDF "       << g2_DC_Chi2N               << ","
    << std::endl    <<
    "Q2: "          << Q2                       << " (GeV/c)2, "
    "xB: "          << xB                       << " , "
    "ω: "           << omega                    << " GeV, "
    "W: "           << W                        << " GeV/c2, "
    "Mx: "          << M_x                      << " GeV/c2, "
    "q: "           << q_p4.P()                 << " GeV/c, "
    << std::endl    <<
    "M2g: "         << M2g                      << " GeV/c2, "
    << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WriteEventToOutput(){
    // (Maybe) write this event to "selected events csv-file"
    bool IsSelected = false;
    
    if (ePastCutsInEvent && pPastCutsInEvent && g1PastCutsInEvent && g2PastCutsInEvent) {
        IsSelected = true;
        Nevents_passed_eep_cuts ++ ;
        
        std::vector<double> variables =
        {   (double)status, (double)runnum,     (double)evnum,
            e_p4.P(),       e_p4.Theta(),       e_p4.Phi(),         Ve.Z(),
            (double)e_DC_sector, e_DC_Chi2N,
            p_p4.P(),       p_p4.Theta(),       p_p4.Phi(),         Vp.Z(),
            (double)p_DC_sector, p_DC_Chi2N,
            g1_p4.P(),      g1_p4.Theta(),      g1_p4.Phi(),        Vg1.Z(),
            (double)g1_DC_sector, g1_DC_Chi2N,
            g2_p4.P(),      g2_p4.Theta(),      g2_p4.Phi(),        Vg2.Z(),
            (double)g2_DC_sector, g2_DC_Chi2N,
            Q2, xB, omega,  W, M_x, q_p4.P(),
            M2g,
        };
        DEBUG(3,"--- -- - electron, proton, 𝛾1 and 𝛾2 passed cuts, writing (e,e'p2𝛾)X event - -- ---");
        if (verbosity > 4) PrintVariables();
        aux.StreamToCSVfile(outcsvfile_eep2gX,
                            variables,
                            csvprecisions );
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FinishProgram(){
    DEBUG(3, "FinishProgram()");
    outcsvfile_eep2gX.close();
    DEBUG(0,"See results at (e,e'p2𝛾)X csv file: %s",full_outcsvfilename.Data());
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// main
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void c12rSkimmer_BranchingRatios(int            RunNumber = 6164,
                                 int           FirstEvent = 0,
                                 int  NeventsMaxToProcess = -1,
                                 int        PrintProgress = 100,
                                 TString        fDataPath = "sidisdvcs",
                                 float             fEbeam = 10.2,
                                 int               fdebug = 0
                                 ){
    SetVerbosity   (fdebug);
    DEBUG(1, "Begin main");
    
    SetGlobals     (fEbeam, fDataPath );
    LoadCutValues  ();
    SetFileNames   ();
    
    
    // open input file and get the hipo data
    TChain fake("hipo");
    fake.Add(infilename.Data());
    auto files = fake.GetListOfFiles();
    
    // open output files
    OpenResultFiles();
    
    // start analysis
    // step over events and extract information....
    for(Int_t i=0;i<files->GetEntries();i++){
        
        //create the event reader
        DEBUG(2, "file %d", i);
        clas12reader c12(files->At(i)->GetTitle(),{0});
        DEBUG(2, "Read title of file %d", i);
        int event = 0;
        
        // process the events...
        while((c12.next()==true) && (event < (FirstEvent + NeventsMaxToProcess))){
            InitializeVariables();
            event++;
            DEBUG(3, "hipo entry %d", event);
            
            if ((event%PrintProgress==0) && (event > FirstEvent))
                DEBUG(3,"Start processing %d/%d (run %d, event %d)",
                      (event-FirstEvent),NeventsMaxToProcess,runnum,evnum);
            
            if (event > FirstEvent) {
                
                runnum = c12.runconfig()->getRun();
                evnum  = c12.runconfig()->getEvent();
                
                // Get Particles By Type
                electrons   = c12.getByID( 11   );
                protons     = c12.getByID( 2212 );
                gammas      = c12.getByID( 22   );
                GetParticlesByType ();
                
                // filter events, extract information, and compute event kinematics
                //
                // we are restricting our analysis to events
                // with exactly one proton and exactly two detected photons
                //
                if(( Ne > 0) && ( Np == 1 ) && ( Ngammas == 2 ))   {
                    
                    DEBUG(2,"Extracting information...");
                    ExtractElectronInformation  ();
                    ExtractProtonInformation    ();
                    ExtractGammasInformation    ();
                    ComputeKinematics           ();
                    WriteEventToOutput          ();
                    DEBUG(2,"Done extracting information...");
                    
                } else {
                    DEBUG(2,"Skipped computation, since N(e)=%d, N(p)=%d, N(gamma)=%d",Ne,Np,Ngammas);
                }
                //                Nevents_processed++;
            } // end if (event > FirstEvent)
            if ((event%PrintProgress==0) && (event > FirstEvent)){
                DEBUG(0,"Done %d/%d",(event-FirstEvent),NeventsMaxToProcess);
                DEBUG(3,"----------------------------------------------------------");
            } // end if (event%PrintProgress==0 && (event > FirstEvent))
        }// end event loop
    } // end file loop
    
    FinishProgram();
    DEBUG(1, "\nDone main.\n");
} // end main

