//   Traject.h -*- C++ -*- Stochastic simulation of QSD trajectories.
//     
//   Copyright (C) 1995  Todd Brun and Ruediger Schack
//   
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//   
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//   
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
//   ----------------------------------------------------------------------
//   If you improve the code or make additions to it, or if you have
//   comments or suggestions, please contact us:
//
//   Dr. Todd Brun			        Tel    +44 (0)171 775 3292
//   Department of Physics                      FAX    +44 (0)181 981 9465
//   Queen Mary and Westfield College           email  t.brun@qmw.ac.uk
//   Mile End Road, London E1 4NS, UK
//
//   Dr. Ruediger Schack                        Tel    +44 (0)1784 443097
//   Department of Mathematics                  FAX    +44 (0)1784 430766
//   Royal Holloway, University of London       email  r.schack@rhbnc.ac.uk
//   Egham, Surrey TW20 0EX, UK
/////////////////////////////////////////////////////////////////////////////

#ifndef _Traject_hhh
#define _Traject_hhh 1

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <map>
#include <chrono>
#include <random>
#include <vector>
#include <fstream>
#include <complex>
#include <iostream>
#include <functional>

#include "State.h"
#include "Operator.h"
#include "CmplxRan.h"

namespace qsd {

	struct Expectation {

		std::complex<double> mean;
	  	std::complex<double> var;

	  	Expectation(double mean_real, double mean_imag, double var_real, double var_imag) {

			  mean.real(mean_real);
			  mean.imag(mean_imag);

			  var.real(var_real);
			  var.imag(var_imag);
	  	}

	  	/*
		Expectation(int nOfSteps, int nX) {

			ExpR.resize(nOfSteps + 1);
			ExpI.resize(nOfSteps + 1);

			for (int j = 0; j <= nOfSteps; j++) ExpR[j].resize(2 * nX);
			for (int j = 0; j <= nOfSteps; j++) ExpI[j].resize(2 * nX);
		}

		std::vector<std::vector<double>> ExpR; // table of real expectations values
		std::vector<std::vector<double>> ExpI; // table of imag expectations values
		*/
	};

    struct TrajectoryResult {

	  bool error;
	  bool timeout;

	  std::vector<double> t; // time that is shared between all of the trajectories
	  State state; // the final state per trajectory
      std::vector<State> state_v;

	  std::map<std::string, std::vector<Expectation>> observables; // operators expectation values per trajectory
	  // std::vector<std::vector<std::vector<std::vector<double>>>> data;
    };

    struct TrajectoryConfig {
	  	std::string log_dir;

	  	bool log_progress = false;
        bool save_on_disk = false;
        bool save_per_traj = false;

        bool dynamic_cutoff = false;
        std::vector<int> dynamic_degrees;
        double cutoff_epsilon = 1e-4;
        int cutoff_pad_size = 2;

        bool phonon_moving_basis = false;
        double shift_accuracy = 1e-4;

        int move = 0;
        double moveEps = 1e-4;

        int nTrajSave = 0;
        int ReadFile = 0;
    };

    class IntegrationStep {
    public:
        virtual void operator()(State &psi, double t, double dt,
                                double &dtlast, ComplexRandom *rndm, std::mt19937 *rnd_gen) = 0;    // perform one step
        virtual ~IntegrationStep();            // destructor
        inline int nLindblads() { return nL; }    // return # of Lindblad Ops
        int getNumdts();            // get # of dts used since last call
    protected:
        Operator H;        // Hamiltonian
        int nL;        // # of Lindblad Operators
        Operator *L;        // List of Lindblad ops
        Operator *Ldag;    // Conjugates of Lindblad ops
        int stochasticFlag;    // Flag used in calculating stochastic component
        std::complex<double> *dxi;        // Array for noise
        int numdtsUsed;    // # of deterministic dts used

        State temp0;        // Temporary states
        State temp1;
        State newsum;

        virtual void derivs(double, State &, State &);

        void error(const char *message);
    };

    class Order2Step : public IntegrationStep {
    public:
        Order2Step(const State &psi, const Operator &theH, int theNL,
                   const Operator *theL);

        virtual void operator()(State &psi, double t, double dt,
                                double &dtlast, ComplexRandom *rndm);

    private:
        State psi2;
        State dpsi;
    };

    class Order4Step : public IntegrationStep {
    public:
        Order4Step(const State &psi, const Operator &theH, int theNL,
                   const Operator *theL);

        virtual void operator()(State &psi, double t, double dt,
                                double &dtlast, ComplexRandom *rndm);

    private:
        State psi0;
        State psi2;
        State dpsi;
    };

    class AdaptiveStep : public IntegrationStep {
    public:
        AdaptiveStep(const State &psi, const Operator &theH,
                     std::vector<Operator> theL, double theEpsilon = 0.000001);

        ~AdaptiveStep();

        virtual void operator()(State &psi, double t, double dt,
                                double &dtlast, ComplexRandom *rndm, std::mt19937 *rnd_gen);

        void dtListSet(int theDim);

        void dtListRead();

        int dtListElem(int theElem);

        void dtListClear();

        void dtListReset();

    protected:
        double epsilon;    // Integration accuracy

        int *dtNumList;    // Array of # of timesteps per dt
        int listDim;        // Dimension of array dtNumList
        int listPtr;        // Current position in array dtNumList

        State dydt;        // Temporary states
        State ak2;
        State ak3;
        State ak4;
        State ak5;
        State ytemp;
        State yout;
        State yerr;
        State y;

        void rkck(double t, double h);

        void rkqs(double &t, double htry, double eps,
                  double &hdid, double &hnext);

        virtual void odeint(State &ystart, double t1, double t2, double eps,
                            double &h1, double hmin, int &nok, int &nbad, std::complex<double> *dxi);
    };

    class AdaptiveStochStep : public AdaptiveStep {
    public:
        AdaptiveStochStep(const State &psi, const Operator &theH,
                           std::vector<Operator> theL, double theEpsilon = 0.000001)
                : AdaptiveStep(psi, theH, theL, theEpsilon) {}

        virtual void operator()(State &psi, double t, double dt,
                                double &dtlast, ComplexRandom *rndm);

    private:
        virtual void odeint(State &ystart, double t1, double t2, double eps,
                            double &h1, double hmin, int &nok, int &nbad, ComplexRandom *rndm);
    };

    class AdaptiveJump : public AdaptiveStep {
    public:
        AdaptiveJump(const State &psi, const Operator &theH,
                     std::vector<Operator> theL, double theEpsilon = 0.000001,
                     char *theFilename = 0)
                : AdaptiveStep(psi, theH, theL, theEpsilon) {
            lindVal = new double[theL.size()];
            outFile = 0;
            if (theFilename != 0)
                outFile = fopen(theFilename, "w");
        }

        ~AdaptiveJump();

        virtual void operator()(State &psi, double t, double dt,
                                double &dtlast, ComplexRandom *rndm, std::mt19937 *rnd_gen);

    private:
        double *lindVal;    // array of expectation values <Ldag L>
        FILE *outFile;    // output file of jump times

        virtual void derivs(double, State &, State &);

        virtual void odeint(State &ystart, double t1, double t2, double eps,
                            double &h1, double hmin, int &nok, int &nbad, ComplexRandom *rndm, std::mt19937 *rnd_gen);
    };

    class AdaptiveOrthoJump : public AdaptiveStep {
    public:
        AdaptiveOrthoJump(const State &psi, const Operator &theH,
                          std::vector<Operator> theL, double theEpsilon = 0.000001,
                          char *theFilename = 0)
                : AdaptiveStep(psi, theH, theL, theEpsilon) {
            lindVal = new double[theL.size()];
            outFile = 0;
            if (theFilename != 0)
                outFile = fopen(theFilename, "w");
        }

        ~AdaptiveOrthoJump();

        virtual void operator()(State &psi, double t, double dt,
                                double &dtlast, ComplexRandom *rndm);

    private:
        double *lindVal;    // array of expectation values <Ldag L>-<Ldag><L>
        FILE *outFile;    // output file of jump times

        virtual void derivs(double, State &, State &);

        virtual void odeint(State &ystart, double t1, double t2, double eps,
                            double &h1, double hmin, int &nok, int &nbad, ComplexRandom *rndm);
    };

    class Trajectory {
    public:

        Trajectory();

        Trajectory(const State &thePsiIni, double thedt, IntegrationStep &theStepper,
                   ComplexRandom *theRand = 0, std::mt19937 *rndGen = nullptr, double theT0 = 0.0);
        // Define a trajectory without computing it.
        // thedt specifies the integration stepsize
        // theRand is a ComplexRandom process; it can be of any of several types.

        void plotExp_obsolete(int nX, const Operator *X, FILE **fp, int *pipe,
                              int dtsperStep, int nOfSteps, int move = 0,
                              double delta = 1e-4, int width = 2, double moveEps = 1e-4);
        // Same as plotExp but:
        // Output expectation and variance of the operators `X[i]'
        // to the files pointed at by `files[i]', i<nX.

        void plotExp(int nX, const Operator *X, char **fname, int *pipe,
                     int dtsperStep, int nOfSteps, int move = 0,
                     double delta = 1e-4, int width = 2, double moveEps = 1e-4, char *savedState = 0);
        // Numerical computation of trajectory (Second order for the deterministic
        // part, Euler for the stochastic part).
        // Time-independent Lindblads only; H can be time dependent. Time dependent Lindblads are implemented via the ion trap interface.
        // Output expectation and variance of the operators `X[i]'
        // to the files named `fname[i]', i<nX.
        // 'pipe' is an array of 4 positive integers that determine which
        // expectation values go to standard output.
        // The first 'move' freedoms are recentered, and their cutoffs adjusted.
        // The cutoff is calculated by removing the top (N - width) states,
        // where the top N states contain total probability less than the
        // threshold delta.
        // Unless the string variable `savedState' is equal to the Null pointer,
        // the final state is saved to the file `savedState'.

		std::vector<TrajectoryResult> sumExp_vector(std::map<std::string, Operator> X, std::vector<std::string> fname, int dtsperStep, int nOfSteps,
				double timeout, int nTrajectory, int worker_id, TrajectoryConfig cfg,
				std::function<void(qsd::TrajectoryResult&, int, double)> callback);
        //
        // Same as sumExp but stores information in a STD vector.
        //
        // Numerical computation of trajectories
        // Same as plotExp, but computing the mean over `nTrajectory' trajectories.
        // No expectation values are written to standard output.
        // Mean expectation and variance of the operators `X[i]' are written
        //   to the files named `fname[i]', i<nX.
        // The means are saved to the files every `nTrajSave' trajectory;
        //   if nTrajSave=0, the results are saved only once at the end.
        // If `ReadFile' is non-zero and if the files `fname[i]' exist,
        //   previously computed means are read from the files `fname[i]' and
        //   taken into account.
        //   IMPORTANT:  Change the seed of the random number generator if
        //   combining mean values in this way!

        State getState();
        // returns the stored value of psi

    private:

        State psi;               // Initial state.
        double dt;               // Basic integration stepsize = step/dtsPerStep.
        ComplexRandom *rndm;     // Pointer to random number generator.
		std::mt19937 *rnd_gen;
        double t0;               // Initial time.
        IntegrationStep *stepper;    // Basic integrator step

        void error(const char *message);    // print error message and exit
    };
}
#endif
