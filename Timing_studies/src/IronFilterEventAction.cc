
/// \file IronFilterEventAction.cc
/// \brief Implementation of the IronFilterEventAction class

#include "IronFilterEventAction.hh"
#include "IronFilterRunAction.hh"
//#include "IronFilterAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include "StepInfo.hh"
#include "G4ThreeVector.hh"

#include "TTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterEventAction::IronFilterEventAction( IronFilterRunAction* input_run_action )
 : G4UserEventAction(),
   stepCollection(),
   run_action(input_run_action),
   eventID(0),
   trackID(0),
   stepID(0),
   parentID(0),
   //particle_name(""),
   //volume_name(""),
   volume_copy_number(0),
   Eki(0),
   Ekf(0),
   edep(0),
   position(0),
   momentum(0),
   x(0),
   y(0),
   z(0),
   px(0),
   py(0),
   pz(0),
   global_time(0),
   //process_name("")
   tmp_particle_name(""),
   tmp_volume_name(""),
   tmp_process_name("")
{
    max_char_len = 15;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


IronFilterEventAction::~IronFilterEventAction(){}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void IronFilterEventAction::PrintEventStatistics() const {}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void IronFilterEventAction::BeginOfEventAction(const G4Event*){

    // If pointer to ROOT tree is empty, then ask RunAction to create the ROOT tree
    // and assign address of variables for output.
    if( data_tree==0 ){

        data_tree = run_action->GetDataTree();

        // Proceed only if data output is enabled.
        if( data_tree!=0 ){
            // information about its order in the event/run sequence
            data_tree->Branch("eventID", &eventID, "eventID/I");
            data_tree->Branch("trackID", &trackID, "trackID/I");
            data_tree->Branch("stepID", &stepID, "stepID/I");

            // information about its idenity
            data_tree->Branch("particle", particle_name, "particle[16]/C");
            data_tree->Branch("parentID", &parentID, "parentID/I");

            // geometric information
            data_tree->Branch("volume", volume_name, "volume[16]/C");
            data_tree->Branch("copy_n", &volume_copy_number, "copy_n/I");
            data_tree->Branch("x", &x, "x/D");
            data_tree->Branch("y", &y, "y/D");
            data_tree->Branch("z", &z, "z/D");
            data_tree->Branch("px", &px, "px/D");
            data_tree->Branch("py", &py, "py/D");
            data_tree->Branch("pz", &pz, "pz/D");

            // dynamic information
            data_tree->Branch("t", &global_time, "t/D");
            data_tree->Branch("Eki", &Eki, "Eki/D"); // initial kinetic energy before the step
            data_tree->Branch("Ekf", &Ekf, "Ekf/D"); // final kinetic energy after the step
            data_tree->Branch("Edep", &edep, "Edep/D"); // energy deposit calculated by Geant4
            data_tree->Branch("process", process_name, "process[16]/C");
        }
    }
    //At the beginning of the event, insert a special flag.
    StepInfo stepinfo;
    stepinfo.SetProcessName( "newEvent" );
    GetStepCollection().push_back(stepinfo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IronFilterEventAction::EndOfEventAction(const G4Event* event){

    if( data_tree!=0 ){

        // Print per event (modulo n)
        G4int evtID = event->GetEventID();
        G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
        if ( ( printModulo > 0 ) && ( evtID % printModulo == 0 ) ) {
            G4cout << "---> End of event: " << evtID << G4endl;
        }

        //if_helium = 0;
        if_first = 0;
        if_second =0;
        if_third =0;
        if_fourth = 0;
        if_fifth = 0;
        if_sixth = 0;

        // Select the tracks of interest
        for( size_t i=0; i < stepCollection.size(); ++i ){

            trackID = stepCollection[i].GetTrackID();
            tmp_volume_name = stepCollection[i].GetVolumeName();
            edep = stepCollection[i].GetDepositedEnergy();
            tmp_particle_name = stepCollection[i].GetParticleName();

            //cout<<tmp_volume_name<<"  "<<tmp_particle_name<<endl;


            //if( (tmp_volume_name=="helium")&& edep!=0){
            //    if_helium = 1;
            //}
            //if( (tmp_volume_name=="CPD")&& edep!=0){
            //    if_helium = 1;
            //}
            if( (tmp_volume_name=="1st_Leadlayer_A"||tmp_volume_name=="1st_FeCap_A")&& tmp_particle_name=="alpha"){
                if_first = 1;
            }
            if( (tmp_volume_name=="2nd_Leadlayer_A"||tmp_volume_name=="2nd_FeCap_A")&& tmp_particle_name=="alpha"){
                if_second = 1;
            }
            if( (tmp_volume_name=="3rd_Leadlayer_A"||tmp_volume_name=="3rd_FeCap_A")&& tmp_particle_name=="alpha"){
                if_third = 1;
            }
            if( (tmp_volume_name=="4th_Leadlayer_B"||tmp_volume_name=="4th_FeCap_B")&& tmp_particle_name=="alpha"){
                if_fourth = 1;
            }
            if( (tmp_volume_name=="5th_Leadlayer_B"||tmp_volume_name=="5th_FeCap_B")&& tmp_particle_name=="alpha"){
                if_fifth = 1;
            }
            if( (tmp_volume_name=="6th_Leadlayer_B"||tmp_volume_name=="6th_FeCap_B")&& tmp_particle_name=="alpha"){
                if_sixth = 1;
            }
      }

        // There is coincidence. Fill the wanted tracks
        //if(if_helium==1){
        //if(if_helium==1 && (if_fourth == 1||if_second == 1)) {
        if(if_first == 1||if_second == 1||if_third == 1||if_fourth == 1||if_fifth == 1||if_sixth == 1) {

            for( size_t i=0; i < stepCollection.size(); ++i ){
              //tmp_volume_name = stepCollection[i].GetVolumeName();
              //if (  (  ( (stepCollection[i].GetParticleName()== "gamma")  || (stepCollection[i].GetParticleName()== "neutron") ) && (stepCollection[i].GetVolumeName()=="helium")) || (stepCollection[i].GetProcessName()=="newEvent") ) {
              //if (  (  ( (stepCollection[i].GetParticleName()== "gamma")  || (stepCollection[i].GetParticleName()== "neutron") ) && (stepCollection[i].GetVolumeName()=="CPD")) || (stepCollection[i].GetProcessName()=="newEvent") ) {
              //if (     (stepCollection[i].GetParticleName()== "gamma")  || (stepCollection[i].GetProcessName()=="newEvent") ) {
              //if ((stepCollection[i].GetParticleName()== "alpha"&& (stepCollection[i].GetVolumeName()!="helium") && (stepCollection[i].GetDepositedEnergy()!=0) )
              //||(stepCollection[i].GetParticleName()== "neutron" && (stepCollection[i].GetVolumeName()=="helium") )) {
              if ( ( (stepCollection[i].GetParticleName()== "alpha") && ( isdigit(stepCollection[i].GetVolumeName()[0]) ) && (stepCollection[i].GetDepositedEnergy() !=0)  ) || (stepCollection[i].GetParticleName()== "neutron")) {
                eventID = stepCollection[i].GetEventID();
                trackID = stepCollection[i].GetTrackID();
                stepID = stepCollection[i].GetStepID();
                edep = stepCollection[i].GetDepositedEnergy();

                parentID = stepCollection[i].GetParentID();

                tmp_particle_name = stepCollection[i].GetParticleName();
                tmp_volume_name = stepCollection[i].GetVolumeName();
                tmp_process_name = stepCollection[i].GetProcessName();


                strncpy( particle_name, tmp_particle_name.c_str(), max_char_len);
                strncpy( process_name, tmp_process_name.c_str(), max_char_len);
                strncpy( volume_name, tmp_volume_name.c_str(), max_char_len);

                volume_copy_number = stepCollection[i].GetVolumeCopyNumber();
                Eki = stepCollection[i].GetEki();
                Ekf = stepCollection[i].GetEkf();

                position = stepCollection[i].GetPosition();
                momentum = stepCollection[i].GetMomentumDirection();

                x = position.x();
                y = position.y();
                z = position.z();

                px = momentum.x();
                py = momentum.y();
                pz = momentum.z();

                global_time = stepCollection[i].GetGlobalTime();

                //G4cout<<eventID<<"   "<<trackID
                //                          <<"  "<<stepID
                //                          <<"  "<<tmp_particle_name
                //                          <<"  "<<tmp_volume_name
                //                          <<"  "<<Eki
                //                          <<"   "<<global_time/1000
                //                          <<"   "<<edep<<G4endl; //time in µs

                data_tree->Fill();
              }

            }
        }

        if_helium = 0;
        if_first = 0;
        if_second =0;
        if_third =0;
        if_fourth = 0;
        if_fifth = 0;
        if_sixth = 0;

    }

    stepCollection.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

vector<StepInfo>& IronFilterEventAction::GetStepCollection(){
    return stepCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
