// last edit Nov-24, 2023

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
int verbosity = 2;

// 4-vectors for the reaction d(e,e'p2𝛾X)
TLorentzVector   Beam_p4, target_p4, e_p4, q_p4, p_p4;
TLorentzVector   p_rest;
TLorentzVector   g1_p4, g2_p4; // gamma 1 and gamma 2
TLorentzVector*    reco_pi0_p4 = NULL; // best fit pi0
TLorentzVector*    reco_eta_p4 = NULL; // best fit eta
std::vector<region_part_ptr>  electrons, protons, gammas;
int                Ne, Np, Ngammas;
int          Nevents_processed = 0;
int                   evnum, runnum;
int               torusBending = -1; // -1 for In-bending, +1 for Out-bending

float                         Ebeam;
TString               Skimming = "";
TString                 prefix = "";
TString               DataPath = "";
TString infilename, outfilepath, outfilename;
ofstream          outcsvfile_eep2gX;

// detector features
int    DC_layer, status, DC_layers[3];
int      Nevents_passed_e_cuts = 0;
int      Nevents_passed_p_cuts = 0;
int    Nevents_passed_eep_cuts = 0;

// leading electron
// electron energy deposit in PCAL [GeV], in ECAL_in [GeV], in ECAL_out [GeV]...
double e_E_PCAL, e_E_ECIN, e_E_ECOUT, e_PCAL_W, e_PCAL_V, e_PCAL_x, e_PCAL_y, e_PCAL_z;
double e_PCAL_sector, e_DC_sector, e_DC_Chi2N, e_DC_x[3], e_DC_y[3], e_DC_z[3];
TVector3 Ve;
bool     ePastCutsInEvent;

// proton
double p_E_PCAL, p_E_ECIN, p_E_ECOUT, p_PCAL_W, p_PCAL_V, p_PCAL_x, p_PCAL_y, p_PCAL_z;
double p_PCAL_sector, p_DC_sector, p_DC_Chi2N, p_DC_x[3], p_DC_y[3], p_DC_z[3];
TVector3 Vp;
bool     pPastCutsInEvent;

// kinematics
double xB, Q2, omega, W, M_x;
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
    p_rest.SetXYZM (0, 0, 0, aux.Mp);
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
void LoadCutValues() {
    // read cut values
    aux.loadCutValues("cuts/cutValues_5Dec2024.csv",torusBending);
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
    
    
    outcsvfile_eep2gX.open( outfilename + "_eep2gX.csv" );
    outcsvfile_eep2gX << csvheader << std::endl;
    
    if (verbosity>1) std::cout << "Done OpenOutputFiles( " << outfilename << ")" << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InitializeVariables(){
    DEBUG(5, "InitializeVariables()");
    electrons   .clear();
    protons     .clear();
    gammas      .clear();
    
    DC_layer    = -9999;
    status      = 1; // 0 is good...
    
    DEBUG(3, "electron...");
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
    
    DEBUG(3, "proton...");
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
    
    DEBUG(5, "eep");
    eepPastCutsInEvent                  = false;
    
    DEBUG(5, "Done InitializeVariables()");
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfProtonPassedSelectionCuts(Double_t p_PCAL_x, Double_t p_PCAL_y,
                                        Double_t p_PCAL_W, Double_t p_PCAL_V,
                                        Double_t p_E_PCAL,
                                        Double_t p_E_ECIN, Double_t p_E_ECOUT,
                                        TLorentzVector p,
                                        TVector3 Vp,
                                        Double_t p_DC_sector,
                                        Double_t p_DC_x[3],
                                        Double_t p_DC_y[3],
                                        Double_t p_DC_z[3],
                                        int torusBending){
    
    if (p_DC_sector == 0) return false;
    
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        int bending  = 1 ? (torusBending==-1) : 0;
        bool DC_fid  = dcfid.DC_fid_xy_sidis(11,                 // particle PID,
                                             p_DC_x[regionIdx],  // x
                                             p_DC_y[regionIdx],  // y
                                             p_DC_sector,        // sector
                                             regionIdx+1,        // layer
                                             bending);           // torus bending
        if (DC_fid == false) {
            return false;
        }
    }
    if(!(true
         &&  p_PCAL_W > aux.cutValue_p_PCAL_W
         &&  p_PCAL_V > aux.cutValue_p_PCAL_V
         &&  p_E_PCAL > aux.cutValue_p_E_PCAL
         && ((p_E_PCAL + p_E_ECIN + p_E_ECOUT)/p_p4.P()) > aux.cutValue_SamplingFraction_min
         && (p_E_ECIN/p_p4.P() > aux.cutValue_PCAL_ECIN_SF_min - p_E_PCAL/p_p4.P())
         &&  ((aux.cutValue_Vz_min < Vp.Z()) && (Vp.Z() < aux.cutValue_Vz_max))
         )) return false;
    
    return true;
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfElectronPassedSelectionCuts(Double_t e_PCAL_x, Double_t e_PCAL_y,
                                        Double_t e_PCAL_W, Double_t e_PCAL_V,
                                        Double_t e_E_PCAL,
                                        Double_t e_E_ECIN, Double_t e_E_ECOUT,
                                        TLorentzVector e,
                                        TVector3 Ve,
                                        Double_t e_DC_sector,
                                        Double_t e_DC_x[3],
                                        Double_t e_DC_y[3],
                                        Double_t e_DC_z[3],
                                        int torusBending){
    
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
    if (e_DC_sector == 0) return false;
    
    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        // DC_e_fid:
        // sector:  1-6
        // layer:   1-3
        // bending: 0(out)/1(in)
        // std::cout << "e_DC_sector: " << e_DC_sector << ", regionIdx: " << regionIdx << std::endl;
        int bending  = 1 ? (torusBending==-1) : 0;
        bool DC_fid  = dcfid.DC_fid_xy_sidis(11,                 // particle PID,
                                             e_DC_x[regionIdx],  // x
                                             e_DC_y[regionIdx],  // y
                                             e_DC_sector,        // sector
                                             regionIdx+1,        // layer
                                             bending);           // torus bending
        if (DC_fid == false) {
            return false;
        }
    }
    
    
    if(!(true
         // fiducial cuts on PCAL
         //fabs(e_PCAL_x)>0
         //&&  fabs(e_PCAL_y)>0
         &&  e_PCAL_W > aux.cutValue_e_PCAL_W
         &&  e_PCAL_V > aux.cutValue_e_PCAL_V
         
         // Electron Identification Refinement  - PCAL Minimum Energy Deposition Cut
         &&  e_E_PCAL > aux.cutValue_e_E_PCAL
         
         // Sampling fraction cut
         && ((e_E_PCAL + e_E_ECIN + e_E_ECOUT)/e.P()) > aux.cutValue_SamplingFraction_min
         && (e_E_ECIN/e.P() > aux.cutValue_PCAL_ECIN_SF_min - e_E_PCAL/e.P()) // RGA AN puts "<" here mistakenly
         
         // Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
         // Spring 19 and Spring 2020 in-bending.
         // Fall 2019 (without low-energy-run) was out-bending.
         &&  ((aux.cutValue_Vz_min < Ve.Z()) && (Ve.Z() < aux.cutValue_Vz_max))
         )) return false;
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractElectronInformation(){
    // ------------------------------------------------------------------------------------------------
    // extract electron information
    // ------------------------------------------------------------------------------------------------
    // find leading electron as the one with highest energy
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

    // hit position in PCAL
    e_PCAL_x        = e_PCAL_info->getX();
    e_PCAL_y        = e_PCAL_info->getY();
    e_PCAL_z        = e_PCAL_info->getZ();

    // Drift Chamber tracking system
    auto e_DC_info  = electrons[leading_e_index]->trk(DC);
    e_DC_sector     = e_DC_info->getSector(); // tracking sector
    e_DC_Chi2N      = e_DC_info->getChi2N();  // tracking chi^2/NDF

    for (int regionIdx=0; regionIdx<3; regionIdx++) {
        int DC_layer = DC_layers[regionIdx];
        e_DC_x[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getX();
        e_DC_y[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getY();
        e_DC_z[regionIdx] = electrons[leading_e_index]->traj(DC,DC_layer)->getZ();
    }
    DEBUG(2,"extracted electron information");

    // ------------------------------------------------------------------------------------------------
    // now, check if electron passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    ePastCutsInEvent = CheckIfElectronPassedSelectionCuts(e_PCAL_x, e_PCAL_y,
                                                          e_PCAL_W, e_PCAL_V,
                                                          e_E_PCAL, e_E_ECIN,
                                                          e_E_ECOUT,
                                                          e_p4, Ve,
                                                          e_PCAL_sector,
                                                          e_DC_x, e_DC_y, e_DC_z,
                                                          torusBending );
    if (ePastCutsInEvent)  Nevents_passed_e_cuts++ ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExtractProtonInformation(){
    // ------------------------------------------------------------------------------------------------
    // extract electron information
    // ------------------------------------------------------------------------------------------------
    // find leading electron as the one with highest energy
    double  leading_p_E;
    int     leading_p_index = 0;
    SetLorentzVector(p_p4,protons[0]);
    TLorentzVector p_tmp(0,0,0,db->GetParticle(11)->Mass());
    for (int pIdx=0; pIdx < Np; pIdx++) {
        SetLorentzVector(p_tmp  ,protons[pIdx]);
        double Ep = p_tmp.E();
        if (Ep > leading_p_E) {
            leading_p_index = pIdx;
            leading_p_E     = Ep;
        }
    }
    // set leading electron 4-momentum
    SetLorentzVector(p_p4 , protons[leading_p_index]);
    // set leading proton vertex
    Vp              = GetParticleVertex( protons[leading_p_index] );
    
    // detector information on electron
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
    DEBUG(2,"extracted proton information");
    
    // ------------------------------------------------------------------------------------------------
    // now, check if proton passed event selection requirements
    // ------------------------------------------------------------------------------------------------
    pPastCutsInEvent = CheckIfProtonPassedSelectionCuts(p_PCAL_x, p_PCAL_y,
                                                        p_PCAL_W, p_PCAL_V,
                                                        p_E_PCAL, p_E_ECIN,
                                                        p_E_ECOUT,
                                                        p_p4, Vp,
                                                        p_PCAL_sector,
                                                        p_DC_x, p_DC_y, p_DC_z,
                                                        torusBending );
    if (pPastCutsInEvent)  Nevents_passed_p_cuts++ ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ComputeElectronKinematics(){
    // compute event kinematics (from e-only information)
    q_p4    = Beam_p4 - e_p4;
    Q2      = -q_p4.Mag2();
    omega   = q_p4.E();
    xB      = Q2/(2. * aux.Mp * q_p4.E());
    W       = sqrt((p_rest + q_p4).Mag2());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WriteEventToOutput(){
    // (Maybe) write this event to "selected events csv-file"
    bool IsSelected = false;
    
    if (ePastCutsInEvent && pPastCutsInEvent) {
        IsSelected = true;
        Nevents_passed_eep_cuts ++ ;
        
        std::vector<double> variables =
        {   (double)status, (double)runnum,     (double)evnum,
            e_p4.P(),       e_p4.Theta(),       e_p4.Phi(),         Ve.Z(),
            Q2,             xB,                 omega,
            (double)e_DC_sector,                (double)p_DC_sector,
            q_p4.P(),
        };
        
        aux.StreamToCSVfile(outcsvfile_eep2gX,
                            variables,
                            csvprecisions );
        
        DEBUG(3,"writing (e,e'p) event");
        
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FinishProgram(){
    DEBUG(2, "FinishProgram()");
    outcsvfile_eep2gX.close();
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

            if (event%PrintProgress==0 && (event > FirstEvent))
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
                if(( 0 < Ne ) && ( Np == 1 ) && ( Ngammas == 2 ))   {

                    DEBUG(2,"Extracting information...");
                    ExtractElectronInformation  ();
//                    ComputeElectronKinematics   ();
//                    ExtractProtonInformation    ();
//                    WriteEventToOutput          ();
                    DEBUG(2,"Done extracting information...");

                } else {
                    DEBUG(2,"Skipped computation, since N(e)=%d, N(p)=%d, N(gamma)=%d",Ne,Np,Ngammas);
                }
                //                Nevents_processed++;
            } // end if (event > FirstEvent)
            if (event%PrintProgress==0 && (event > FirstEvent)){
                DEBUG(1,"Done %d/%d",(event-FirstEvent),NeventsMaxToProcess);
                DEBUG(3,"----------------------------------------------------------");
            } // end if (event%PrintProgress==0 && (event > FirstEvent))
        }// end event loop
    } // end file loop
                
    FinishProgram();
    DEBUG(1, "\nDone main.\n");
} // end main

