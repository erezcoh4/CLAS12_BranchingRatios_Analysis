
// last edit Dec-5, 2024 (EOC)

#include "BranchingRatio_CLAS12_Auxiliary.h"
#define r2d 180./3.1415 // radians to degrees


BranchingRatio_CLAS12_Auxiliary::BranchingRatio_CLAS12_Auxiliary(int _fdebug_, int _torusBending_){
//    SetVerbosity    (_fdebug_);
//    SetTorusBending (_torusBending_);
}

BranchingRatio_CLAS12_Auxiliary::~BranchingRatio_CLAS12_Auxiliary(){}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BranchingRatio_CLAS12_Auxiliary::loadCutValues(std::string cutValuesFilename,
                                                    int torusBending){
    if (fdebug>2) {
        std::cout << "BranchingRatio_CLAS12_Auxiliary::loadCutValues()" << std::endl;
    }
    // read cut values csv file
    csv_reader csvr;
    cutValues = csvr.read_csv(cutValuesFilename);

    // assign specific cut values - to speed things up
    // by avoiding recalling FindCutValue() on every event

    // Cut on z-vertex position:
    // based on RGA_Analysis_Overview_and_Procedures_Nov_4_2020-6245173-2020-12-09-v3.pdf
    // p. 71
    if (torusBending==-1){ // in-bending torus field
        // Spring 19 and Spring 2020 in-bending.
        cutValue_Vz_min = FindCutValue("Vz_e_min_inbending");
        cutValue_Vz_max = FindCutValue("Vz_e_max_inbending");
    } else if (torusBending==1){ // Out-bending torus field
        // Fall 2019 (without low-energy-run) was out-bending.
        cutValue_Vz_min = FindCutValue("Vz_e_min_outbending");
        cutValue_Vz_max = FindCutValue("Vz_e_max_outbending");

    } else {
        std::cout
        << "Un-identified torus bending "
        << torusBending
        << ", return" << std::endl;
        return;
    }

    cutValue_e_PCAL_W               = FindCutValue("e_PCAL_W_min");
    cutValue_e_PCAL_V               = FindCutValue("e_PCAL_V_min");
    cutValue_e_E_PCAL               = FindCutValue("e_E_PCAL_min");
    cutValue_SamplingFraction_min   = FindCutValue("SamplingFraction_min");
    cutValue_PCAL_ECIN_SF_min       = FindCutValue("PCAL_ECIN_SF_min");
    cutValue_Ve_Vpi_dz_max          = FindCutValue("(Ve-Vpi)_z_max");
    cutValue_Q2_min                 = FindCutValue("Q2_min");
    cutValue_Q2_max                 = FindCutValue("Q2_max");
    cutValue_W_min                  = FindCutValue("W_min");
    cutValue_y_max                  = FindCutValue("y_max");
    cutValue_e_theta_min            = FindCutValue("e_theta_min");
    cutValue_e_theta_max            = FindCutValue("e_theta_max");
    cutValue_pi_theta_min           = FindCutValue("pi_theta_min");
    cutValue_pi_theta_max           = FindCutValue("pi_theta_max");
    cutValue_Ppi_min                = FindCutValue("Ppi_min");
    cutValue_Ppi_max                = FindCutValue("Ppi_max");
    cutValue_Pe_min                 = FindCutValue("Pe_min");
    cutValue_Pe_max                 = FindCutValue("Pe_max");
    cutValue_Zpi_min                = FindCutValue("Zpi_min");
    cutValue_Zpi_max                = FindCutValue("Zpi_max");

    // Kaons
    cutValue_Ve_VK_dz_max          = FindCutValue("(Ve-VK)_z_max");
    cutValue_K_theta_min           = FindCutValue("K_theta_min");
    cutValue_K_theta_max           = FindCutValue("K_theta_max");
    cutValue_PK_min                = FindCutValue("PK_min");
    cutValue_PK_max                = FindCutValue("PK_max");
    cutValue_ZK_min                = FindCutValue("ZK_min");
    cutValue_ZK_max                = FindCutValue("ZK_max");


    if (fdebug>2) { printCutValues(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BranchingRatio_CLAS12_Auxiliary::printCutValues(){
    std::cout << "Using the following cut values:" << std::endl;
    std::cout <<
    "Vz_min: "                  << cutValue_Vz_min                  << ", " << std::endl <<
    "Vz_max: "                  << cutValue_Vz_max                  << ", " << std::endl <<
    "e_PCAL_W: "                << cutValue_e_PCAL_W                << ", " << std::endl <<
    "e_PCAL_V: "                << cutValue_e_PCAL_V                << ", " << std::endl <<
    "e_E_PCAL: "                << cutValue_e_E_PCAL                << ", " << std::endl <<
    "SamplingFraction_min: "    << cutValue_SamplingFraction_min    << ", " << std::endl <<
    "PCAL_ECIN_SF_min: "        << cutValue_PCAL_ECIN_SF_min        << ", " << std::endl <<
    "Ve_Vpi_dz_max: "           << cutValue_Ve_Vpi_dz_max           << ", " << std::endl <<
    "Q2_min: "                  << cutValue_Q2_min                  << ", " << std::endl <<
    "Q2_max: "                  << cutValue_Q2_max                  << ", " << std::endl <<
    "W_min: "                   << cutValue_W_min                   << ", " << std::endl <<
    "y_max: "                   << cutValue_y_max                   << ", " << std::endl <<
    "e_theta_min: "             << cutValue_e_theta_min             << ", " << std::endl <<
    "e_theta_max: "             << cutValue_e_theta_max             << ", " << std::endl <<
    "pi_theta_min: "            << cutValue_pi_theta_min            << ", " << std::endl <<
    "pi_theta_max: "            << cutValue_pi_theta_max            << ", " << std::endl <<
    "Ppi_min: "                 << cutValue_Ppi_min                 << ", " << std::endl <<
    "Ppi_max: "                 << cutValue_Ppi_max                 << ", " << std::endl <<
    "Pe_min: "                  << cutValue_Pe_min                  << ", " << std::endl <<
    "Pe_max: "                  << cutValue_Pe_max                  << ", " << std::endl <<
    "Zpi_min: "                 << cutValue_Zpi_min                 << ", " << std::endl <<
    "Zpi_max: "                 << cutValue_Zpi_max                 << ", " << std::endl <<
    "Ve_VK_dz_max: "            << cutValue_Ve_VK_dz_max            << ", " << std::endl <<
    "K_theta_min: "             << cutValue_K_theta_min             << ", " << std::endl <<
    "K_theta_max: "             << cutValue_K_theta_max             << ", " << std::endl <<
    "PK_min: "                  << cutValue_PK_min                  << ", " << std::endl <<
    "PK_max: "                  << cutValue_PK_max                  << ", " << std::endl <<
    "ZK_min: "                  << cutValue_ZK_min                  << ", " << std::endl <<
    "ZK_max: "                  << cutValue_ZK_max                  << ", " << std::endl <<
    std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double BranchingRatio_CLAS12_Auxiliary::FindCutValue( std::string cutName ){
    for (auto cut: cutValues) {
        if (strcmp(cut.first.c_str(),cutName.c_str())==0){
            return cut.second;
        }
    }
    return -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TString BranchingRatio_CLAS12_Auxiliary::GetRunNumberSTR( int RunNumber, TString fSkimming ){
    char RunNumberStr[20];
    // sprintf( RunNumberStr, "00%d", RunNumber );

    if(fSkimming == "p_uniform_distribution"){
        // "white" GEMC simulation runs
        sprintf( RunNumberStr, "%d", RunNumber );
    } else {
        sprintf( RunNumberStr, "%06d", RunNumber );
    }
    if (fdebug>1) std::cout << "(SIDIS) skimming run " << RunNumberStr << std::endl;
    return (TString)RunNumberStr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Double_t BranchingRatio_CLAS12_Auxiliary::GetEbeamFromRunNumber ( Int_t RunNumber ){
    if (6420 <= RunNumber && RunNumber <= 6598){
        return 10.2; // GeV
    }
    else if (11362 <= RunNumber && RunNumber <= 11571){
        return 10.4; // GeV
    }
    else if (6164 <= RunNumber && RunNumber <= 6399){
        return 10.6; // GeV
    }
    else{
        return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BranchingRatio_CLAS12_Auxiliary::SetTorusBendingFromRunNumber ( Int_t RunNumber ){
    // -1 for In-bending, +1 for Out-bending
    // For BAND data
    // Spring 19 and Spring 2020 was in-bending
    // Fall 2019 (without low-energy-run) was out-bending

    if (6420 <= RunNumber && RunNumber <= 6598){
        this->torusBending = -1;
    }
    else if (11362 <= RunNumber && RunNumber <= 11571){
        this->torusBending = -1;
    }
    else if (6164 <= RunNumber && RunNumber <= 6399){
        this->torusBending = +1;
    }
    else{
        this->torusBending = 0;
    }
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void BranchingRatio_CLAS12_Auxiliary::StreamToCSVfile (std::ofstream&         csvfile,
                                                       std::vector<Double_t>  observables,
                                                       std::vector<int>       precisions){

    //    for (auto v:observables) csvfile << v << ",";
    //    csvfile << std::endl;

    for (int j=0; j < observables.size(); j++){
        int precision = 9;
        if (j < precisions.size()) precision = precisions.at(j);
        auto v = observables.at(j);
        csvfile << std::setprecision(precision) << std::fixed << v << ",";
    }
    csvfile << std::endl;


    if (fdebug>3) {
        std::cout << "StreamToEventCSVfile()" << std::endl;
        for (auto v:observables) std::cout << v << ",";
        std::cout << std::endl;
    }
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void BranchingRatio_CLAS12_Auxiliary::OpenCSVfile (std::ofstream& csvfile,
                                                   TString filename,
                                                   std::string header){
    csvfile.open( filename );
    csvfile << header << "," << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BranchingRatio_CLAS12_Auxiliary::Print4Vector( TLorentzVector v, std::string label ){
    std::cout << label << " 4-vector:"<<std::endl;
    std::cout
    << std::setprecision(2)
    << "(Px = " << v.Px()
    << ", Py = " << v.Py()
    << ", Pz = " << v.Pz()
    << ", E = " << v.E()
    << "), M = " << v.Mag()
    << " GeV"
    << std::endl;
}




// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
double BranchingRatio_CLAS12_Auxiliary::ComputeLightConeFraction( TLorentzVector p ){
    // compute light-cone momentum fraction
    double m = p.Mag();
    double alpha = (p.E() - p.Z())/m;
    return alpha;
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
void BranchingRatio_CLAS12_Auxiliary::PrintMonitorHello(){
    std::cout << "Hello..." << std::endl;
    std::cout << "Is it me you're looking for?..." << std::endl;
    std::cout << "I can see it in your eyes..." << std::endl;
    std::cout << "I can see it in your smile..." << std::endl;
    std::cout << "You're all I've ever wanted, and your arms are open wide..." << std::endl;
    std::cout << "Cause you know just what to say, and you know just what to do..." << std::endl;
    std::cout << "And I want to tell you so much, I love you" << std::endl;
    std::cout << "I long to see the sunlight in your hair" << std::endl;
}




