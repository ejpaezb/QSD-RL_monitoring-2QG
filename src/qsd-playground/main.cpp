/*
 * Copyright (c) 2019, Seyed Shakib Vedaie & Eduardo J. Paez
 */

#include <iostream>
#include <vector>
#include <chrono>

#include <qsd/ACG.h>
#include <qsd/Traject.h>
#include <qsd/State.h>
#include <qsd/Operator.h>
#include <qsd/AtomOp.h>
#include <qsd/FieldOp.h>
#include <qsd/SpinOp.h>
#include <qsd/Complex.h>

int main(int argc, char *argv[])
{
    Complex I{0.0, 1.0};

    qsd::SigmaX sx_1{0};
    qsd::SigmaX sx_2{1};
    qsd::SigmaX sx_3{2};

    qsd::AnnihilationOperator a_1{1};

    /*
    * The initial state
    * */

    std::vector<qsd::State> psilist;

    psilist.emplace_back(2, qsd::SPIN);
    psilist.emplace_back(5, qsd::FIELD);
    // // psilist.emplace_back(2, qsd::SPIN);
    // // psilist.emplace_back(2, qsd::SPIN);

    qsd::State psi0 = qsd::State(psilist.size(), psilist.data());

    psi0 *= sx_1;
    // psi0 *= sx_2;
    // psi0 *= sx_3;

    psi0 *= a_1.hc();
    // psi0 *= a_1.hc();
    // psi0 *= a_1.hc();
    // psi0 *= a_1.hc();

    std::cout << psi0;
}

int main_1(int argc, char *argv[])
{
    Complex I{0.0, 1.0};

    qsd::SigmaZ sz_1{0};
    qsd::SigmaX sx_1{0};

    qsd::SigmaPlus sp_1{0};
    qsd::SigmaPlus sp_2{1};

    qsd::IdentityOperator id_s1{0};

    /*
    * The initial state
    * */
    std::vector<int> phonon_cutoffs{80,80,80,80};

    std::vector<qsd::State> psilist;
    for (int i = 0; i < 4 + 2; i++) {
        if (i < 2) psilist.emplace_back(2, qsd::SPIN);
        else {
            psilist.emplace_back(phonon_cutoffs[i - 2], 0, qsd::FIELD);
        }
    }

    qsd::State psi0 = qsd::State(psilist.size(), psilist.data());

    /*
    psi0 *= sp_1;
    psi0 *= sp_2;

    State psi_0_exact(2, SPIN);

    State psi_1(2, SPIN);
    State psi_2(2, SPIN);

    psi_2 *= sp_1;
    */

    qsd::State psi1 = psi0;
    qsd::State psi2 = psi0;

    psi2 *= sp_1;

    auto begin = std::chrono::high_resolution_clock::now();

    psi0 *= sx_1;

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Exact method:\n";
    std::cout << psi0 * psi1 << std::endl;
    std::cout << psi0 * psi2 << std::endl;

    std::cout << "Total time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms\n";

}

int main_2(int argc, char *argv[])
{
    std::complex<double> I{0.0, 1.0};

    qsd::SigmaZ sz_1{0};
    qsd::SigmaX sx_1{0};
    qsd::SigmaPlus sp_1{0};
    qsd::IdentityOperator id_s1{0};

    // exp^(-i pi/2 X) = - i X

    // I - iX * dt
    // I + (- iX * dt) - 1/2 * I * dt^2

    /*
    * The initial state
    * */

    qsd::State psi_0_exact(2, qsd::SPIN);
    qsd::State psi_0_appx(2, qsd::SPIN);

    qsd::State psi_1(2, qsd::SPIN);
    qsd::State psi_2(2, qsd::SPIN);

    psi_2 *= sp_1;

    auto begin = std::chrono::high_resolution_clock::now();

    psi_0_exact *= -sx_1;
    psi_0_exact *= I;

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Exact method:\n";
    std::cout << psi_0_exact * psi_1 << std::endl;
    std::cout << psi_0_exact * psi_2 << std::endl;

    std::cout << "Total time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " us\n";

    double t = 1.570796326794897;
    int steps = 1000;
    double dt = t / steps;

    begin = std::chrono::high_resolution_clock::now();

    qsd::Operator U_appx = id_s1 - I * sx_1 * dt - 0.5 * id_s1 * dt * dt;
    for (int i = 0; i < steps; i++) {
        psi_0_appx *= U_appx;
    }

    end = std::chrono::high_resolution_clock::now();

    std::cout << "Approximate method:\n";
    std::cout << psi_0_appx * psi_1 << std::endl;
    std::cout << psi_0_appx * psi_2 << std::endl;

    std::cout << "Total time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " us\n";

    return 0;
}