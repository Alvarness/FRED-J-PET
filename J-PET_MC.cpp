using namespace std;

double ScatteredGammaEnergy(double InitialEnergy, double ComptonScatteringAngle){

    double denominator = 1 + (InitialEnergy/511)*(1 - cos(ComptonScatteringAngle));

    return InitialEnergy/denominator;


}

double KleinNishinaCrossSection(double InitialEnergy /* In keV*/, double ComptonScatteringAngle, double EtaAngle){

    double result = pow((ScatteredGammaEnergy( InitialEnergy, ComptonScatteringAngle)/InitialEnergy), 2) * (
            InitialEnergy/ScatteredGammaEnergy(InitialEnergy, ComptonScatteringAngle) + ScatteredGammaEnergy( InitialEnergy, ComptonScatteringAngle)/InitialEnergy -
            2*pow(sin(ComptonScatteringAngle),2)*pow(cos(EtaAngle), 2));

    return result;

}

double Ratio(double InitialEnergy, double ComptonScatteringAngle){

   return  KleinNishinaCrossSection(InitialEnergy, TMath::PiOver2(), ComptonScatteringAngle)/KleinNishinaCrossSection(InitialEnergy, 0, ComptonScatteringAngle);

}


double gamma(double momentum, double theta) {

    return ScatteredGammaEnergy(momentum, theta)/momentum + momentum/ScatteredGammaEnergy(momentum, theta);
}


double Visibility(double momentum, double theta){

    return sin(theta)*sin(theta) / (gamma(momentum, theta) - sin(theta)*sin(theta));


}

double NKleinNishinaCrossSection(double InitialEnergy /* In keV*/, double ComptonScatteringAngle){

    double result = 0.00781015*pow((ScatteredGammaEnergy( InitialEnergy, ComptonScatteringAngle)/InitialEnergy), 2) * (
            InitialEnergy/ScatteredGammaEnergy(InitialEnergy, ComptonScatteringAngle) + ScatteredGammaEnergy( InitialEnergy, ComptonScatteringAngle)/InitialEnergy -
            pow(sin(ComptonScatteringAngle),2));

    return result;

}

double IsInsideCircle(double x, double y, double r){

    if(pow(x*180/TMath::Pi() - 81.66, 2) + pow(y*180/TMath::Pi() - 81.66, 2) < r*r) return true;
    else return false;


}



void EntanglementMCanalysis(){

    TRandom3 generator(0);
// for(int file = 1; file < 2; file++){

    TFile input_file("/home/niki/work/MC/J-PET-geant4/build/mcGeant.root");
    // TFile input_file("/home/niki/work/MC/J-PET-geant4/build/bin/mcGeant.root", "READ");
    // TFile input_file("/home/niki/work/run5_mc/out.mcGeant.root", "READ");

    // cout << Form("/home/niki/work/run5_mc/%i.mcGeant.root", file) << endl;

    // TFile output(Form(   "/net/scratch/people/plgnikodem/%i/%i_out.mcGeant.root", file, file),"recreate");

    TCanvas *canvas = new TCanvas("Canvas_1", "Canvas_1", 533, 76, 1383, 852);
    TH1F *compton = new TH1F("compton", "compton", 180, 0, 180);


    TTree* tree = (TTree*)input_file.Get("T");

    double speed_of_light = 0.0299792458; // cm ps^-1
    double electron_mass = 510.998928; // electron mass in keV

    JPetGeantEventPack* eventPack = new JPetGeantEventPack();
    JPetGeantEventInformation* evtInfo = new JPetGeantEventInformation();
    JPetGeantScinHits *MCHit = new JPetGeantScinHits();



    tree->SetBranchAddress("eventPack", &eventPack);
    int nevent = tree->GetEntries();                 ////////////////////
    cout << "events to process " << nevent << endl;


    int toprint = 0;
    double minEdep = 0.0;
    int counter = 0;
    int counter2 = 0;

    for(int i = 0; i < nevent; i++){

        // if (100.0*i/nevent >= (toprint + 1.)) {
        //     toprint = toprint + 1.;
        //     cout << "Running......" << 100.0 * i/nevent << "%........." << endl;
        // }

        tree->GetEvent(i);

        evtInfo = eventPack->GetEventInformation();
        int NumberOfHits = eventPack->GetNumberOfHits();

        for(int j = 0; j < NumberOfHits; j++){
            if(eventPack->GetHit(j)->GetGenGammaMultiplicity() == 2){
                    counter2++;
                    // cout << "dupa" << endl;
                    TVector3 mom_in = eventPack->GetHit(j)->GetMomentumIn();
                    TVector3 mom_out = eventPack->GetHit(j)->GetMomentumOut();
                    TVector3 pol_in = eventPack->GetHit(j)->GetPolarizationIn();
                    TVector3 pol_out = eventPack->GetHit(j)->GetPolarizationOut();
                    compton->Fill(pol_in.Angle(pol_out)*TMath::RadToDeg() );
                    if(abs(mom_in.Angle(pol_in)*TMath::RadToDeg() - 90) > pow(10, -7)){
                        counter++;
                       cout << pol_in.Angle(mom_in)*TMath::RadToDeg() << endl;
                    }


            }
        }




    }
    cout << counter << " " << counter2<< endl;
    compton->Draw("H");
    canvas->SaveAs("S.png");
}
// }


