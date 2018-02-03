//------------------------------------------------------------------------------
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------
#include <iostream>
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <string.h>
#include <sstream>

#include "qureg/qureg.hpp"
#include "interface_api_qubitid.h"
#include "interface_api_version.h"
#include "interface_api_memory.h"

#if (defined(__ICC) || defined(__INTEL_COMPILER))
#include <mkl.h>
#if defined(OPENQU_HAVE_MPI)
#include <mkl_cdft.h>
#endif
#endif

template<typename Type>
void qft(NoisyQureg<Type> *psi)
{
  int n = openqu::ilog2(psi->size());

  // main computation
  for (int i = n - 1; i >= 0; i--) {
  	for (int j = n - 1; j > i; j--) {
  		int k = j - i;
		openqu::TinyMatrix<Type, 2, 2, 32> phaseshift;
			phaseshift(0, 0) = {1, 0};
			phaseshift(0, 1) = {0, 0};
			phaseshift(1, 0) = {0, 0};
			phaseshift(1, 1) = Type(cos(M_PI / D(UL(1) << UL(k))), sin(M_PI / D(UL(1) << UL(k))));
			psi->applyControlled1QubitGate(j, i, phaseshift);
	}
	psi->applyHadamard(i);
   }
  // perform swapping
  for (int i = 0; i < (n / 2); i++) {
  	psi->applySwap(i, n - 1 - i);
  }
}

using namespace std;


using Type = ComplexDP;
extern NoisyQureg<Type> *psi1;


// Constant defining the rotational angle of a T-dagger gate. Basically, -(pi/4).
#define TDAG_THETA -0.785398163397448


unsigned long unk(string args) {
    return 1;
}


unsigned long S_handler(string args) {
    cout << "S"<< " [" << args << "]" <<endl;
    psi1->applyPauliSqrtZ(query_qubit_id(args));
    return 0;
}


unsigned long X_handler(string args) {
    cout << "X"<< " [" << args << "]" <<endl;
    psi1->applyPauliX(query_qubit_id(args));
    return 0;
}


unsigned long T_handler(string args) {
    cout << "T"<< " [" << args << "]" <<endl;
    psi1->applyT(query_qubit_id(args));
    return 0;
}

unsigned long R_handler(string args){
   stringstream ss(args);
   string angle, trueargs, token;
   int flg = 0;
   while (getline(ss, token, ',')) {
      if(flg == 0){
         angle += token;
         flg = 1;
      }else{
         trueargs += token + ",";
      }
   }
   //cout << "atruargs: " << trueargs << " | angle: " << angle <<endl;
   cout<< "R"<< " [" << args << "] " <<endl;
   psi1->applyRotationZ(query_qubit_id(trueargs), stod(angle));
   return 0;
}

unsigned long QFT_handler(string args){
   cout << "QFT " << " [" << args << "] " << endl;
   using Type = ComplexDP;
   qft<Type>(psi1);
   return 0;
}

unsigned long Tdag_handler(string args) {
    cout << "Tdag"<< " [" << args << "]" <<endl;
    psi1->applyRotationZ(query_qubit_id(args),TDAG_THETA);
    return 0;
}


unsigned long CNOT_handler(string args) {
    int qubit1,
        qubit2;
    int token_end = args.find_first_of(',');

    qubit1 = query_qubit_id(args.substr(0,token_end));
    qubit2 = query_qubit_id(args.substr(token_end+1,args.length()));

    cout << "CNOT"<< " [" << args << "]" <<endl;
    psi1->applyCPauliX(qubit1,qubit2);
    return 0;
}


unsigned long H_handler(string args) {
    cout << "H"<< " [" << args << "]" <<endl;
    psi1->applyHadamard(query_qubit_id(args));
    return 0;
}

unsigned long Noise_handler(string args){
    cout << "Noise"<< " ["<< args << "] " << endl;
    psi1->apply_noise_gates_on_all_qubits();
    return 0;
}

unsigned long MeasZ_handler(string args) {
    using Type = ComplexDP;
    Type measurement = 0.0;
    
    cout << "MeasZ"<< " [" << args << "]" <<endl;
    measurement = psi1->getProbability(query_qubit_id(args));
    cout << measurement << endl;
    return 0;
}


unsigned long PrepZ_handler(string args) {
    cout << "PrepZ"<< " [" << args << "]" <<endl;
    return 0;
}


// Hash table containing the QASM operation string and the function to call to
// handle the operation with the qHiPSTER simulation.
//
unordered_map<string, function<long(string)>> qufun_table = {\
                                                {".malloc", qumalloc},
                                                {".free", qufree},
                                                {".iversion",quiversion},
                                                {".version",quversion},
                                                {"H", H_handler},
                                                {"CNOT", CNOT_handler},
                                                {"PrepZ",PrepZ_handler},
                                                {"T", T_handler},
                                                {"X", X_handler},
                                                {"Tdag", Tdag_handler},
                                                {"S", S_handler},
                                                {"MeasZ", MeasZ_handler},
						{"Noise", Noise_handler},
						{"QFT", QFT_handler},
						{"R", R_handler},
                                                {"*", unk},
};



unsigned long ExecuteHandler(string op, string args) {

    unsigned long result = 1;

    function<long(string)> func = qufun_table[op];

    if(func) {
        result = func(args);
    }

    return result;
}


