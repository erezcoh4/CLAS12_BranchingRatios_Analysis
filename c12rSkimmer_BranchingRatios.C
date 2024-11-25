// last edit Nov-24, 2023

#include "Auxiliary/BranchingRatio_CLAS12_Auxiliary.h"

#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>+;
#endif

using namespace clas12;
aux;

// Results in CSV file d(e,e'p2𝛾X)
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
int verbosity = 2;

// 4-vectors for the reaction d(e,e'p2𝛾X)
TLorentzVector*        beam_p4 = NULL; // e-beam
TLorentzVector*           e_p4 = NULL; // e'
TLorentzVector*           q_p4 = NULL; // q 4 vector
TLorentzVector*           p_p4 = NULL; // q 4 vector
TLorentzVector*          g1_p4 = NULL; // gamma 1
TLorentzVector*          g2_p4 = NULL; // gamma 2
TLorentzVector*    reco_pi0_p4 = NULL; // best fit pi0
TLorentzVector*    reco_eta_p4 = NULL; // best fit eta
std::vector<region_part_ptr>  electrons, protons, gammas;


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void Debug(int v, const char* fmt, ...) {
    va_list arg;
    va_start(arg, fmt);
    
    if (verbosity > v) {
        vprintf(fmt, arg);
        std::cout << std::endl;
    }
    //vprintf can be replaced with vsprintf (for sprintf behavior)
    //or any other printf function preceded by a v
    va_end(arg);
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetVerbosity( int v ){
    verbosity = v;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetEbeam (double fEbeam) { // [GeV]
    // RGA the enrgy was 10.6
    // RGB Spring-2019 the enrgy was 10.2
    // RGB Fall-2019 the enrgy was 10.4096
    Ebeam = fEbeam;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetGlobals(int v, float fEbeam) {
    SetVerbosity        ( v          );
//    SetDataPath         ( fDataPath, fEbeam );
//    SetSkimming         ( fSkimming  );
    SetEbeam            ( fEbeam     );
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void SetFileNames() {
    TString RunNumberStr = aux.GetRunNumberSTR(RunNumber, Skimming);
    // define input filename
    TString infilename, outfilepath, outfilename;
    
    infilename  = DataPath + prefix + RunNumberStr + ".hipo";
    outfilepath = "/volatile/clas12/users/ecohen/RGB/" + Skimming + "/";
    outfilename = "skimmed_BranchingRatios_" + prefix + RunNumberStr;
    
    if (fdebug>1){
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
    aux.loadCutValues("cuts/RGBcutValues.csv",torusBending);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GetParticlesByType (){
    // get particles by type
    Ne      = electrons .size();
    Np      = protons   .size();
    Ngammas = gammas    .size();
}


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
// main
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void c12rSkimmer_BranchingRatios(int RunNumber    = 6420)
{
    Debug(1, "Begin main");
    
//    SetGlobals     ();
//    LoadCutValues  ();
//    SetFileNames   ();
    
    
    // open input file and get the hipo data
    TChain fake("hipo");
    fake.Add(infilename.Data());
    auto files = fake.GetListOfFiles();
    
    // open output files
    OpenResultFiles( outfilepath, outfilename );
    
    // start analysis
    // step over events and extract information....
    for(Int_t i=0;i<files->GetEntries();i++){
        
        //create the event reader
        clas12reader c12(files->At(i)->GetTitle(),{0});
        //        InitializeFileReading( NeventsMax, c12.getReader().getEntries(), fdebug );
        int event = 0;
        
        // process the events...
        while((c12.next()==true) && (event < NeventsMaxToProcess)){
            event++;
            
            if (event > FirstEvent) {
                
                runnum = c12.runconfig()->getRun();
                evnum  = c12.runconfig()->getEvent();
                
                InitializeVariables();
                // Get Particles By Type
                electrons   = c12.getByID( 11   );
                protons     = c12.getByID( 2212 );
                gammas      = c12.getByID( 22   );
                GetParticlesByType ( evnum, fdebug );
                
                Debug(1,"N(e):%d, N(p):%d, N(g):%d ",Ne,Np,Ngammas)
                
                // filter events, extract information, and compute event kinematics:
                // ....
                
                Nevents_processed++;
            }
            if (event%PrintProgress==0 && (event > FirstEvent)) Debug(1,"%d/%d",event,NeventsMaxToProcess);
            
        } // end event loop
    } // end file loop
    Debug(1, "\nDone.\n");
}
