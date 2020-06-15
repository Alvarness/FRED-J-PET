double ScatteredGammaEnergy(double InitialEnergy, double ComptonScatteringAngle){

    double denominator = 1 + (InitialEnergy/511)*(1 - cos(ComptonScatteringAngle));

    return InitialEnergy/denominator;


}

double KleinNishinaCrossSection(double InitialEnergy /* In keV*/, double ComptonScatteringAngle, double EtaAngle){

    double result = pow((ScatteredGammaEnergy( InitialEnergy, ComptonScatteringAngle)/InitialEnergy), 2) * \
    (InitialEnergy/ScatteredGammaEnergy(InitialEnergy, ComptonScatteringAngle) + ScatteredGammaEnergy( InitialEnergy, ComptonScatteringAngle)/InitialEnergy -
            2 + 4*pow(cos(EtaAngle),2));

    return result;

}


int main(){


    TF1 *f = new TF1("f", "KleinNishinaCrossSection(511, 81*TMath::DegToRad(), TMath::DegToRad()*x)", 0, 180);

    f->Draw();


    return 0;
}

