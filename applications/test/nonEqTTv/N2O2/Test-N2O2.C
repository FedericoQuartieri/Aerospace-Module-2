// heat bath 0D di N2 + O2 (fig. 6 del paper): raffreddamento vibrazionale con
// due temperature vibrazionali, scambio V-T e, a scelta, scambio V-V
//
// Mutation++ ha una sola temperatura vibrazionale e non ha lo scambio V-V,
// quindi qui le sorgenti sono scritte a mano seguendo il paper (eq. 8 e 18);
// da Mutation++ si prendono le energie delle specie e i tempi di
// rilassamento di Millikan-White con la correzione di Park (eq. 9-17)
//
// uso: Test-N2O2 <VV: on|off> <t_fine> <file_csv>
//   fig 6: Test-N2O2 off 3e-6 output/fig6-noVV.csv
//          Test-N2O2 on  3e-6 output/fig6-VV.csv

#include "mutation++.h"
#include "HarmonicOscillator.h"
#include "MillikanWhite.h"
#include "heatBath.H"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// energia vibrazionale e_v ed elettronica e_el (J/kg) della specie s alla
// temperatura Tv, dalle entalpie adimensionali di Mutation++ (eq. 5 e 7)
static void speciesVibElEnergies
(
    Mutation::Mixture& mix,
    const int s,
    const double Tv,
    double& e_v,
    double& e_el
)
{
    const int ns = mix.nSpecies();
    std::vector<double> h_v(ns), h_el(ns);
    // tutte le temperature uguali a Tv: h_v/RT e h_el/RT dipendono solo da Tv
    mix.speciesHOverRT(Tv, Tv, Tv, Tv, Tv, nullptr, nullptr, nullptr,
                       h_v.data(), h_el.data(), nullptr);
    const double RT_over_M = Mutation::RU * Tv / mix.speciesMw(s);
    e_v = h_v[s] * RT_over_M;
    e_el = h_el[s] * RT_over_M;
}

// temperatura vibro-elettronica della specie s dalla sua energia e_ve (J/kg),
// per bisezione (e_ve cresce con Tv)
static double TvFromEve(Mutation::Mixture& mix, const int s, const double e_ve)
{
    double Tlow = 100.0;
    double Thigh = 100000.0;
    for (int iter = 0; iter < 60; iter++)
    {
        const double Tv = 0.5 * (Tlow + Thigh);
        double e_v, e_el;
        speciesVibElEnergies(mix, s, Tv, e_v, e_el);
        if (e_v + e_el > e_ve)
        {
            Thigh = Tv;
        }
        else
        {
            Tlow = Tv;
        }
    }
    return 0.5 * (Tlow + Thigh);
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cerr << "uso: Test-N2O2 <VV: on|off> <t_fine> <file_csv>" << std::endl;
        return 1;
    }
    const bool withVV = (std::string(argv[1]) == "on");
    const double t_end = std::atof(argv[2]);
    const char* csvName = argv[3];

    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism("none");
    Mutation::Mixture mix(opts);

    const int ns = mix.nSpecies();
    const int iN2 = mix.speciesIndex("N2");
    const int iO2 = mix.speciesIndex("O2");
    const double M_N2 = mix.speciesMw(iN2);
    const double M_O2 = mix.speciesMw(iO2);
    const double R_N2 = Mutation::RU / M_N2;
    const double R_O2 = Mutation::RU / M_O2;

    // ---- stato iniziale del paper: T_tr = 5000 K, T_v = 30000 K, 1 atm,
    // N2 e O2 in parti uguali (stesso numero di moli)
    double T = 5000.0;
    double Tv_N2 = 30000.0;
    double Tv_O2 = 30000.0;

    std::vector<double> Y(ns, 0.0);
    Y[iN2] = M_N2 / (M_N2 + M_O2);
    Y[iO2] = M_O2 / (M_N2 + M_O2);
    const double P_T_Tv[3] = {Mutation::ONEATM, T, Tv_N2};
    mix.setState(Y.data(), P_T_Tv, 2);

    std::vector<double> rho_s(ns);
    mix.densities(rho_s.data());
    const double rho_N2 = rho_s[iN2];
    const double rho_O2 = rho_s[iO2];

    // energie per unita' di volume (eq. 23): quella totale E si conserva,
    // quelle vibro-elettroniche dei due serbatoi cambiano nel tempo
    const double E = totalEnergy(mix, rho_s);
    double e_v, e_el;
    speciesVibElEnergies(mix, iN2, Tv_N2, e_v, e_el);
    double Eve_N2 = rho_N2 * (e_v + e_el);
    speciesVibElEnergies(mix, iO2, Tv_O2, e_v, e_el);
    double Eve_O2 = rho_O2 * (e_v + e_el);

    // energia traslazionale-rotazionale iniziale: cresce con T come
    // 2.5 R T per molecola (eq. 3 e 4), e' il modo per ricavare T dopo
    const double T0 = T;
    const double Etr0 = E - Eve_N2 - Eve_O2;
    const double rhoCvTr = 2.5 * (rho_N2 * R_N2 + rho_O2 * R_O2);
    // ---- fine stato iniziale

    // T finale dalla conservazione dell'energia (il paper la legge in figura)
    const double T_eq = equilibriumTemperature(mix, rho_s, E);

    // ---- tempi di rilassamento V-T di Mutation++ (Millikan-White + Park,
    // eq. 9-17, con le costanti del suo file VT.xml)
    Mutation::Thermodynamics::HarmonicOscillatorDB hoDB;
    Mutation::Transfer::MillikanWhiteModelDB mwDB(mix);
    Mutation::Transfer::MillikanWhiteModel tauModel_N2 =
        mwDB.create("N2", hoDB.create("N2").characteristicTemperatures()[0]);
    Mutation::Transfer::MillikanWhiteModel tauModel_O2 =
        mwDB.create("O2", hoDB.create("O2").characteristicTemperatures()[0]);
    // ---- fine tempi di rilassamento

    // ---- costanti dello scambio V-V (eq. 18)
    // probabilita' di scambio consigliata dal paper
    const double P_N2O2 = 0.01;
    // sezione d'urto della coppia N2-O2: il paper non da' il numero, questo e'
    // quello usato da hy2Foam, il codice del paper (file thermo2TModel,
    // KnabCoefficients N2_O2, github.com/vincentcasseau/hyStrath, commit 984e3000a5f8)
    const double sigma_N2O2 = 2.667e-19;
    // massa molare ridotta della coppia (eq. 14), kg/mol
    const double M_N2O2 = M_N2 * M_O2 / (M_N2 + M_O2);
    // ---- fine costanti V-V

    std::ofstream csv(csvName);
    csv << "# T_eq = " << T_eq << " K (conservazione dell'energia)\n";
    csv << "t,Ttr,Tv_N2,Tv_O2\n";
    csv << "0," << T << "," << Tv_N2 << "," << Tv_O2 << "\n";

    // ---- ciclo nel tempo: passo di 1 ns come il paper
    const double dt = 1.0e-9;
    double t = 0.0;
    int step = 0;

    while (t < t_end)
    {
        // stato di Mutation++ a T_tr: serve per i tempi di rilassamento
        const double temps[2] = {T, Tv_N2};
        mix.setState(rho_s.data(), temps, 1);
        const double tau_N2 = tauModel_N2.relaxationTime(mix);
        const double tau_O2 = tauModel_O2.relaxationTime(mix);

        // energie vibro-elettroniche a T_tr (obiettivo) e alle T_v attuali
        double ev_N2_T, eel_N2_T, ev_O2_T, eel_O2_T;
        speciesVibElEnergies(mix, iN2, T, ev_N2_T, eel_N2_T);
        speciesVibElEnergies(mix, iO2, T, ev_O2_T, eel_O2_T);
        double ev_N2, eel_N2, ev_O2, eel_O2;
        speciesVibElEnergies(mix, iN2, Tv_N2, ev_N2, eel_N2);
        speciesVibElEnergies(mix, iO2, Tv_O2, ev_O2, eel_O2);

        // scambio V-T di ogni molecola, Landau-Teller (eq. 8)
        const double Q_N2_VT = rho_N2 * (ev_N2_T + eel_N2_T - ev_N2 - eel_N2) / tau_N2;
        const double Q_O2_VT = rho_O2 * (ev_O2_T + eel_O2_T - ev_O2 - eel_O2) / tau_O2;

        // scambio V-V fra N2 e O2 (eq. 18, con le sole energie vibrazionali);
        // quello che entra in N2 esce da O2
        double Q_N2_VV = 0.0;
        if (withVV)
        {
            Q_N2_VV = Mutation::NA * sigma_N2O2 * P_N2O2
                    * std::sqrt(8.0 * Mutation::RU * T / (Mutation::PI * M_N2O2))
                    * (rho_O2 / M_O2) * rho_N2
                    * (ev_N2_T * ev_O2 / ev_O2_T - ev_N2);
        }
        const double Q_O2_VV = -Q_N2_VV;

        // eq. 22: cambiano i due serbatoi vibro-elettronici, E si conserva
        Eve_N2 += (Q_N2_VT + Q_N2_VV) * dt;
        Eve_O2 += (Q_O2_VT + Q_O2_VV) * dt;
        t += dt;
        step++;

        // nuove temperature: T_v dai serbatoi, T_tr dall'energia che resta
        Tv_N2 = TvFromEve(mix, iN2, Eve_N2 / rho_N2);
        Tv_O2 = TvFromEve(mix, iO2, Eve_O2 / rho_O2);
        T = T0 + (E - Eve_N2 - Eve_O2 - Etr0) / rhoCvTr;

        if (step % 10 == 0)
        {
            csv << t << "," << T << "," << Tv_N2 << "," << Tv_O2 << "\n";
        }
    }
    // ---- fine ciclo nel tempo

    std::cout << csvName << ": T_tr = " << T << " K, T_v,N2 = " << Tv_N2
              << " K, T_v,O2 = " << Tv_O2 << " K a t = " << t
              << " s; T_eq attesa = " << T_eq << " K" << std::endl;

    return 0;
}
