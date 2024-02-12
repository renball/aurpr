#ifndef NEYRON_H
#define NEYRON_H

#include <QObject>
#include <cmath>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#define Nn 4
#define V f[0]
#define m f[1]
#define N f[2]
#define h f[3]

using namespace std;

class Neyron : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int value READ getValue WRITE setValue NOTIFY IsynChanged)
private:

    const double C = 1; // muF/cm^2

    const double g_K = 35; // mS/cm^2
    const double g_Na = 40; // mS/cm^2
    const double g_L = 0.3; // mS/cm^2

    const double E_K = -77; // mV
    const double E_Na = 55; // mV
    const double E_L = -65; // mV

    double Isyn;
    double E_syn = 0;
    double Theta_syn = 0;
    double k_syn = 0.2;
    double g_syn = 0.1;

public:
    explicit Neyron(QObject *parent = nullptr);



    double f[4];


    double getValue (){
        return V;
    }

    void setValue(double I){
        Isyn = I;
    }

    double alpha_n(double f[4])
    {
        return 0.02 * (V - 25) / (1 - exp(-(V - 25) / 9));
    }

    double beta_n(double f[4])
    {
        return -0.002 * (V - 25) / (1 - exp((V - 25) / 9));
    }

    double alpha_m(double f[4])
    {
        return 0.182 * (V + 35) / (1 - exp(-(V + 35) / 9));
    }

    double beta_m(double f[4])
    {
        return -0.124 * (V + 35) / (1 - exp((V + 35) / 9));
    }

    double alpha_h(double f[4])
    {
        return 0.25 * exp(-(V + 90) / 12);
    }

    double beta_h(double f[4])
    {
        return 0.25 * exp((V + 62) / 6) / exp((V + 90) / 12);
    }

    double HodgkinHuxley(int i, double f[4], double I, int n_ind, double V_1)
    {
        switch (i)
        {
        case 0:
            if ((n_ind % 2) == 0) {
                return 1000 * ((g_Na * pow(m, 3) * h * (E_Na - V) + g_K * N * (E_K - V) + g_L * (E_L - V) + I) / C);
            }
            else {
                Isyn = g_syn * (E_syn - V) / (1 + exp(-((V_1)-Theta_syn) / k_syn));
                return 1000 * ((g_Na * pow(m, 3) * h * (E_Na - V) + g_K * N * (E_K - V) + g_L * (E_L - V) + I + Isyn) / C);
            }
        case 1:
            return 1000 * (alpha_m(f) * (1 - m) - beta_m(f) * m);

        case 2:
            return 1000 * (alpha_n(f) * (1 - N) - beta_n(f) * N);

        case 3:
            return 1000 * (alpha_h(f) * (1 - h) - beta_h(f) * h);
        }
        return 0;
    }

    void RungeKutta(double dt, double f[4], double f_next[4], double I, int n_ind, double V_1)
    {
        double k[4][4];

        // k1
        for (int i = 0; i < 4; i++)
            k[i][0] = HodgkinHuxley(i, f, I, n_ind, V_1) * dt;

        double phi_k1[4];
        for (int i = 0; i < 4; i++)
            phi_k1[i] = f[i] + k[i][0] / 2;

        // k2
        for (int i = 0; i < 4; i++)
            k[i][1] = HodgkinHuxley(i, phi_k1, I, n_ind, V_1) * dt;

        double phi_k2[4];
        for (int i = 0; i < 4; i++)
            phi_k2[i] = f[i] + k[i][1] / 2;

        // k3
        for (int i = 0; i < 4; i++)
            k[i][2] = HodgkinHuxley(i, phi_k2, I, n_ind, V_1) * dt;

        double phi_k3[4];
        for (int i = 0; i < 4; i++)
            phi_k3[i] = f[i] + k[i][2] / 2;

        // k4
        for (int i = 0; i < 4; i++)
            k[i][3] = HodgkinHuxley(i, phi_k3, I, n_ind, V_1) * dt;

        for (int i = 0; i < 4; i++)
            f_next[i] = f[i] + (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
    }

    double HodgkinHuxley_default(int i, double f[4], double I)
    {
        switch (i)
        {
        case 0:
            return 1000 * ((g_Na * pow(m, 3) * h * (E_Na - V) + g_K * N * (E_K - V) + g_L * (E_L - V) + I) / C);

        case 1:
            return 1000 * (alpha_m(f) * (1 - m) - beta_m(f) * m);

        case 2:
            return 1000 * (alpha_n(f) * (1 - N) - beta_n(f) * N);

        case 3:
            return 1000 * (alpha_h(f) * (1 - h) - beta_h(f) * h);
        }
        return 0;
    }

    void RungeKutta_default(double dt, double f[4], double f_next[4], double I)
    {
        double k[4][4];

        // k1
        for (int i = 0; i < 4; i++)
            k[i][0] = HodgkinHuxley_default(i, f, I) * dt;

        double phi_k1[4];
        for (int i = 0; i < 4; i++)
            phi_k1[i] = f[i] + k[i][0] / 2;

        // k2
        for (int i = 0; i < 4; i++)
            k[i][1] = HodgkinHuxley_default(i, phi_k1, I) * dt;

        double phi_k2[4];
        for (int i = 0; i < 4; i++)
            phi_k2[i] = f[i] + k[i][1] / 2;

        // k3
        for (int i = 0; i < 4; i++)
            k[i][2] = HodgkinHuxley_default(i, phi_k2, I) * dt;

        double phi_k3[4];
        for (int i = 0; i < 4; i++)
            phi_k3[i] = f[i] + k[i][2] / 2;

        // k4
        for (int i = 0; i < 4; i++)
            k[i][3] = HodgkinHuxley_default(i, phi_k3, I) * dt;

        for
                (int i = 0; i < 4; i++)
            f_next[i] = f[i] + (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
    }

    void CopyArray(double source[4], double target[4])
    {
        for (int i = 0; i < 4; i++)
            target[i] = source[i];
    }

    Q_INVOKABLE void neyro_dif() {
        const int count=2;
        Neyron neyrons[count];

        for (int i = 0; i < count; i++) {
            neyrons[i].V = -58.7085; neyrons[i].m = 0.0953; neyrons[i].N = 0.000913; neyrons[i].h = 0.3662;		//-58.7085 0.0953 0.000913 0.3662
        }

        double tmax = 1.0; double dt = 0.00005;
        double t = 0;
        int counter = 0;
        double I = 0.0;
        while (I < 1.1) {
            t = 0.0;
            while (t <= tmax) {
                double fn_next[Nn];
                RungeKutta_default(dt, f, fn_next, I);
                CopyArray(fn_next, f);
                cout << I << " " << t << " " << V << " " << m << " " << N << " " << h << endl;
                if(V>10){
                    counter++;
                    cout << "**************************************************************" << endl;
                }
                t += dt;
            }

            I += 0.01;
        }
        cout << counter;
}

    Q_INVOKABLE void neyrorest(){
        f[0]=0.0;
        f[1]=0.0;
        f[2]=0.0;
        f[3]=0.0;
    }


    Q_INVOKABLE double neyro_dif_2(double t){
        double I=1.1;
        double fn_next[Nn];
        RungeKutta_default(0.00005, f, fn_next, I);
        CopyArray(fn_next, f);
        return V;
    }

    Q_INVOKABLE double getIsyn(){
        return g_K;
    }

    Q_INVOKABLE double getI_V(){
        double I = 0.0;
        double t=0.0;
        const int count=2;
        Neyron neyrons[count];

        for (int i = 0; i < count; i++) {
            neyrons[i].V = -58.7085; neyrons[i].m = 0.0953; neyrons[i].N = 0.000913; neyrons[i].h = 0.3662;		//-58.7085 0.0953 0.000913 0.3662
        }

        for (int i = 0; i < count; i++) {
            t++;
            cout << t << " " << neyrons[i].V << " " << endl;
        }

        return neyrons[0].V;

    }

signals:
    void IsynChanged();
};

#endif // NEYRON_H
