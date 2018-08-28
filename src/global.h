// -*- c++ -*-
#ifndef GLOBAL_H
#define GLOBAL_H

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <assert.h>
#include <cassert>
#include <random>
#include <algorithm>
#include <functional>
#include <iterator>
#include <sys/stat.h>
#include <errno.h>    // errno, ENOENT, EEXIST
#if defined(_WIN32)
#include <direct.h>   // _mkdir
#endif

#include "rng.h"

/*SodaPop
Copyright (C) 2017 Louis Gauthier

    SodaPop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    SodaPop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with SodaPop.  If not, see <http://www.gnu.org/licenses/>.
 */

#define POPSIZEMAX 	1000000
#define GENECOUNTMAX 	100

/******** CONTAINER TYPEDEFS ********/
typedef std::vector<int> VectInt;
typedef std::vector<double> VectDouble;
typedef std::vector<std::string> VectStr;


/******* CONSTANTS *********/

/*****
N.B. The physically allowed value for mutational DDG is DGG_min to DGG_max.
If the estimated energy is out of this range, the mutation is ignored.
*****/
const double ddG_min = -7.5;
const double ddG_max = 99;
const double CONC_MAX = 1e15;
const double kT = 0.5922; //defines the energy units
const double COST = 1e-4; // see Geiler-Samerotte et al. 2011
const double A_FACTOR = 0.025; // see Serohijos & Shakhnovich 2013
                                // A_factor = sum_i{1/abundance_i}
const double fNs = 0.775956284; //fraction of non-synonymous substitutions in a typical protein
extern double X_FACTOR;
extern double DRUG_CONCENTRATION;
extern double DRUG_INCREASE_FACTOR;
extern double DRUG_DECREASE_FACTOR;


// exponent values are precalculated to be used readily
const double DDG_min = exp(-1*(ddG_min)/kT);
const double DDG_max = exp(-1*(ddG_max)/kT);
const int Bigbuffer_max = 80;
constexpr double PI  = 3.141592653589793238463;
constexpr double AVOGADRO = 6.02214e23;
constexpr double CELL_VOLUME = 7e-16; // in liters (L)


// If the mutation is to a stop codon
// DG_mutant is set to 99 kcal/mol 
// -> all copies are effectively aggregated
const double DG_STOP = exp(-1*(99)/kT);

// Create a 3D matrix for fitness landscape
const int max_gene = 1200;
const int max_resi = 640;

/******** Global variables (non-const) ********/

extern VectStr PrimordialAASeq;
extern double matrix[max_gene][max_resi][20];

/******* FUNCTION DECLARATIONS *******/
int GetIndexFromCodon(std::string);
std::string GetProtFromNuc(std::string);
int GetIndexFromAA(std::string);
int GetIndexFromAA(char);
int CheckBP(std::string);
std::string AdjacentBP(std::string, int);
std::string n3_to_n3(std::string, std::string, int);
std::string getBarcode();
void InitMatrix();
void ExtractPDDGMatrix(std::string);
void ExtractDMSMatrix(std::string);
double Ran_Gaussian(const double mean, const double sigma);
void LoadPrimordialGenes(const std::string&,const std::string&);
void qread_Cell(std::fstream& IN, std::fstream& OUT);
void read_Cell(std::fstream& IN, std::fstream& OUT);
int StringDiff(const std::string&, const std::string&);
std::string trim(const std::string&);
bool isDirExist(const std::string&);
bool makePath(const std::string&);

/******* GENETIC CODE MAPPINGS *******/
// these const mappings are hard-coded

struct codon_to_num{
    static  std::map<std::string,int> create_map();
    static const std::map<std::string,int> cnum;
};

struct codon_to_prot{
    static std::map<std::string,std::string> create_map();
    static const std::map<std::string,std::string> cprot;
};

struct prot_to_num{
    static std::map<std::string,int> create_map();
    static const std::map<std::string,int> pnum;
};

#endif  // GLOBAL_H
