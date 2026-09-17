// heat bath 0D di N2 + O2 (fig. 6 del paper): raffreddamento vibrazionale con
// due temperature vibrazionali, scambio V-T e, a scelta, scambio V-V
//
// Mutation++ ha una sola temperatura vibrazionale e non ha lo scambio V-V,
// quindi qui il ciclo tiene un serbatoio per molecola e usa le funzioni di
// mutationSources.H specie per specie (energie, tempi di rilassamento);
// lo scambio V-V e' scritto seguendo l'eq. 18 del paper
//
// uso: Test-N2O2 <VV: on|off> <t_fine> <file_csv>
//   fig 6: Test-N2O2 off 3e-6 output/fig6-noVV.csv
//          Test-N2O2 on  3e-6 output/fig6-VV.csv

#include "mutation++.h"
#include "mutationSources.H"
#include "heatBath.H"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// temperatura vibro-elettronica della specie s dalla sua energia e_ve (J/kg),
// per bisezione (e_ve cresce con Tv)
static double TvFromEve(Mutation::Mixture& mix, const int s, const double e_ve)
{
    double Tlow = 100.0;
    double Thigh = 100000.0;
    for (int iter = 0; iter < 60; iter++)
    {
        const double Tv = 0.5 * (Tlow + Thigh);
        if (speciesEve(mix, s, Tv) > e_ve)
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
    double Eve_N2 = rho_N2 * speciesEve(mix, iN2, Tv_N2);
    double Eve_O2 = rho_O2 * speciesEve(mix, iO2, Tv_O2);

    // energia traslazionale-rotazionale iniziale: cresce con T come
    // 2.5 R T per molecola (eq. 3 e 4), e' il modo per ricavare T dopo
    const double T0 = T;
    const double Etr0 = E - Eve_N2 - Eve_O2;
    const double rhoCvTr = 2.5 * (rho_N2 * R_N2 + rho_O2 * R_O2);
    // ---- fine stato iniziale

    // T finale dalla conservazione dell'energia (il paper la legge in figura)
    const double T_eq = equilibriumTemperature(mix, rho_s, E);

    // ---- tempi di rilassamento V-T di Mutation++ (eq. 9-17), uno per molecola
    const std::vector<Vibrator> vibrators = makeVibrators(mix);
    const Vibrator* vibN2 = nullptr;
    const Vibrator* vibO2 = nullptr;
    for (const Vibrator& v : vibrators)
    {
        if (v.species == iN2) vibN2 = &v;
        if (v.species == iO2) vibO2 = &v;
    }
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
    const int nSteps = int(t_end / dt + 0.5);
    const int writeEvery = std::max(10, nSteps / 10000);
    double t = 0.0;

    for (int step = 1; step <= nSteps; step++)
    {
        // stato di Mutation++ a T_tr: serve per i tempi di rilassamento
        const double temps[2] = {T, Tv_N2};
        mix.setState(rho_s.data(), temps, 1);

        // scambio V-T di ogni molecola, Landau-Teller (eq. 8), con la forza
        // motrice e_ve(T_tr) - e_ve(T_v) della molecola stessa
        const double Q_N2_VT = rho_N2 * (speciesEve(mix, iN2, T) - speciesEve(mix, iN2, Tv_N2))
                             / vibN2->tau.relaxationTime(mix);
        const double Q_O2_VT = rho_O2 * (speciesEve(mix, iO2, T) - speciesEve(mix, iO2, Tv_O2))
                             / vibO2->tau.relaxationTime(mix);

        // scambio V-V fra N2 e O2 (eq. 18, con le sole energie vibrazionali);
        // quello che entra in N2 esce da O2
        double Q_N2_VV = 0.0;
        if (withVV)
        {
            const double ev_N2_T = speciesEv(mix, iN2, T);
            const double ev_O2_T = speciesEv(mix, iO2, T);
            const double ev_N2 = speciesEv(mix, iN2, Tv_N2);
            const double ev_O2 = speciesEv(mix, iO2, Tv_O2);
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

        // nuove temperature: T_v dai serbatoi, T_tr dall'energia che resta
        Tv_N2 = TvFromEve(mix, iN2, Eve_N2 / rho_N2);
        Tv_O2 = TvFromEve(mix, iO2, Eve_O2 / rho_O2);
        T = T0 + (E - Eve_N2 - Eve_O2 - Etr0) / rhoCvTr;

        if (step % writeEvery == 0)
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
