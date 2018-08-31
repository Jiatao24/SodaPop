// -*- c++ -*-
#ifndef POLYCELL_H
#define POLYCELL_H

#include "Cell.h"

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

class PolyCell : public Cell
{
    typedef double (PolyCell::*funcPtr)(void);

public:
    static int ff_;
    static bool useDist_;
    static bool fromS_;

    PolyCell();
    PolyCell(std::fstream&);			    
    PolyCell(std::fstream& binaryInfile, const std::string& genePath);

    // Fitness functions
    void selectFitness();
    double flux();
    double toxicity();
    double metabolicOutput();
    double multiplicative();
    double neutral();
    double stochasticExpression();
    double enzymaticOutput();
    double stochasticEnzymaticOutput();
    double fitness();           //!< Getter for fitness_
    void UpdateRates();		//!< This function updates fitness_.

    void ranmut_Gene();
    void ranmut_Gene(std::ofstream&, int);
    void mutGene(int g_num, int resid, std::string resname);
    void dump(std::fstream&, int);
    void dumpShort(std::fstream&);
    void PrintCell(int);
    void UpdateNsNa();
    int Na(){return Total_Na_;}
    int Ns(){return Total_Ns_;}

protected:
    double fitness_;
    int Total_Ns_;
    int Total_Na_;
    funcPtr fit;

};

#endif  // POLYCELL_H
