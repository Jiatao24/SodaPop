#include <unistd.h>
#include <cmath>

#include <tclap/CmdLine.h>

#include "nlohmann/json.hpp"

#include "rng.h"
#include "global.h"
#include "PolyCell.h"

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

git branch moran_model
 */

using json = nlohmann::json;

void countGenotypes(std::vector<PolyCell> cells,
                    std::map<std::string, unsigned int>& genotypeCounts)
{
    for (auto& cell : cells)
    {
        std::string genotype = cell.get_genotype();
        auto search = genotypeCounts.find(genotype);
        if (search == genotypeCounts.end())
        {
            genotypeCounts.insert({genotype, 1});
        }
        else
        {
            search->second++;
        }
    }
    return;
}


//! Save a more succint json-based snapshot of genotype populations.
void saveSnapshot(char* buffer, std::vector<PolyCell>& cells,
                  int generationNumber)
{
    // Do genotype counting
    std::map<std::string, unsigned int> genotypeCounts;
    countGenotypes(cells, genotypeCounts);

    // We also want fitnesses
    // We do repeat ourselves here
    // (fitnesses are also calculated in main)
    std::map<std::string, double> genotypeFitnesses;

    // Fill genotypeFitnesses and count genotype populations.
    for (auto& cell : cells)
    {
        std::string genotype = cell.get_genotype();
        double fitness = cell.fitness();
        auto search = genotypeFitnesses.find(genotype);
        if (search == genotypeFitnesses.end())
            genotypeFitnesses.insert({genotype, fitness});
        else
            search->second += fitness;
    }

    // Get average for each genotype by dividing by count
    for (auto& it : genotypeFitnesses)
    {
        it.second /= (double)genotypeCounts[it.first];
    }

    // Calculate the overall average fitness.
    double averageFitness(0);
    for (auto& it : genotypeFitnesses)
    {
        averageFitness += it.second * genotypeCounts[it.first];
    }
    averageFitness /= (double)cells.size();

    // Let's get abundances
    int gene_count = cells.front().gene_count();
    std::vector<double> abundances(gene_count, 0);
    std::map<std::string, std::vector<double>> genotypeAbundances;

    // Fill abundance data only if fitness is stochastic concentration.
    if (PolyCell::ff_ == 8)
    {
        for (auto& cell : cells)
        {
            std::vector<double> cell_abundances = cell.get_abundances();
            std::string genotype = cell.get_genotype();

            auto search = genotypeAbundances.find(genotype);
            if (search == genotypeAbundances.end())
            {
                genotypeAbundances.insert(
                    {genotype, std::vector<double>(gene_count, 0)});
            }

            // Instead of getting mean by dividing by cell count
            //   later, just divide each value by count right now.
            //   who knows if it's worse performance anyway?
            for (int i=0; i<gene_count; i++)
            {
                abundances[i] += cell_abundances[i] / (double)cells.size();
                genotypeAbundances[genotype][i] +=
                    cell_abundances[i] / (double)genotypeCounts[genotype];
            }
        }
    }

    // Open snapshot file
    std::fstream OUT2(buffer, std::ios::out);
    if (!OUT2.is_open())
    {
         std::cerr << "Output file could not be opened";
         exit(1);
    }

    json output;
    output["generation"] = generationNumber;
    output["population"] = cells.size();
    output["drug level (nM)"] = DRUG_CONCENTRATION;
    output["average fitness"] = averageFitness;
    output["genotypes"] = genotypeCounts;
    output["genotype fitnesses"] = genotypeFitnesses;
    if (PolyCell::ff_ == 8)
    {
        output["average abundances"] = abundances;
        output["genotype abundances"] = genotypeAbundances;
    }
    OUT2 << std::setw(2) << output << std::endl;

    OUT2.close();   
    
    return;
}


int main(int argc, char *argv[])
{
    // these variables will hold the parameters input (or not) by the user
    unsigned int generationNumber = 0;
    unsigned int generationMax = generationNumber + 1;
    unsigned int mutationCount = 0;
    unsigned int equilibrationGens = 0;
    unsigned int populationSize = 1;
    unsigned int populationMax = 1;
    unsigned int outputFreq = 1;
    unsigned int snapoutputFreq = 0;
    double TIME = 0;
    char buffer[200];
    bool trackMutations = false;
    bool rampingDrug = false;
    bool bottleneck = false;
    unsigned int bottleneckSize = 1;
    unsigned int bottleneckInterval = 1;

    std::string geneListFile, genesPath;
    std::string outDir, startSnapFile, matrixFile;

    // Wrap everything in a try block
    // errors in input are caught and explained to user
    try
    { 
        // Define the command line object
        TCLAP::CmdLine cmd("SodaPop: a multi-scale model of molecular evolution", ' ', "SodaPop/0.1.3b-biophysical");

        // Define value arguments
        TCLAP::ValueArg<unsigned int> maxArg("m","maxgen","Maximum number of generations",false,10000,"int");
        TCLAP::ValueArg<unsigned int> popArg("n","size","Initial population size",false,1,"positive int");
        TCLAP::ValueArg<unsigned int> popMaxArg("", "max-size","Maximum population size",false,1,"int");
        TCLAP::ValueArg<unsigned int> dtArg("t","dt","Time interval for json snapshots",false,1,"int");
        TCLAP::ValueArg<unsigned int> fulldtArg("","fullsnap-dt","Time interval for full snapshots",false,0,"int");

        //files
        TCLAP::ValueArg<std::string> prefixArg("o","prefix","Prefix to be used for snapshot files",false,"sim","filename");
        TCLAP::ValueArg<std::string> geneArg("g","gene-list","Gene list file",true,"null","filename");
        TCLAP::ValueArg<std::string> startArg("p","pop-desc","Population description file",true,"null","filename");
        TCLAP::ValueArg<std::string> libArg("l","gene-lib","Gene library directory",false,"files/genes/","filename");

        // fitness function
        TCLAP::ValueArg<int> fitArg("f","fitness","Fitness function",false,1,"integer ID");
    
        // boolean switch to track mutations
        TCLAP::SwitchArg eventsArg("e","track-events","Track mutation events", cmd, false);
    
        // RNG seed
        TCLAP::ValueArg<unsigned long> seedArg("", "seed", "Seed value for RNG.", false, 0, "unsigned int (64-bit)");

        // RNG stream
        TCLAP::ValueArg<unsigned long> streamArg(
            "", "stream", "Stream value for RNG (set different streams for different runs to assure different "
            "random numbers for each run. But no two runs with the same stream will be identical unless the same seed "
            "is used as well.)",
            false, 0, "unsigned int (64-bit)");

        // X_FACTOR
        TCLAP::ValueArg<double> xfactorArg("x", "x-factor", "Enzymatic output factor.", false, 0, "positive double");

        TCLAP::ValueArg<double> concentrationArg("", "concentration", "Amount of drug present (nM)", false, 0, "non-negative double");

        TCLAP::SwitchArg rampingArg("", "ramping", "Drug concentration adjusts to population fitness.", cmd, false);

        TCLAP::ValueArg<unsigned int> equilArg("", "equil", "Time before mutation", false, 0, "nonnegative int");

        TCLAP::ValueArg<unsigned int> bottleneckSizeArg("", "bottleneck-size", "how much to reduce population", false, 0, "positive int");

        TCLAP::ValueArg<unsigned int> bottleneckIntervalArg("", "bottleneck-interval", "how often to reduce population", false, 0, "positive int");

        // Add the arguments to the CmdLine object.
        cmd.add(seedArg);
        cmd.add(streamArg);
        cmd.add(maxArg);
        cmd.add(popArg);
        cmd.add(popMaxArg);
        cmd.add(dtArg);
        cmd.add(fulldtArg);
        cmd.add(prefixArg);
        cmd.add(geneArg);
        cmd.add(startArg);
        cmd.add(libArg);
        cmd.add(fitArg);
        cmd.add(xfactorArg);
        cmd.add(concentrationArg);
        cmd.add(equilArg);
        cmd.add(bottleneckSizeArg);
        cmd.add(bottleneckIntervalArg);

        // Parse the argv array.
        cmd.parse(argc, argv);

        // Get values from args. 
        generationMax = maxArg.getValue();
        populationSize = popArg.getValue();
        populationMax = popMaxArg.getValue();
        outputFreq = dtArg.getValue();
        snapoutputFreq = fulldtArg.getValue();

        geneListFile = geneArg.getValue();
        outDir = prefixArg.getValue();
        startSnapFile = startArg.getValue();
        genesPath = libArg.getValue();

        equilibrationGens = equilArg.getValue();

        if (seedArg.isSet())
            setRngSeed(seedArg.getValue());

        if (streamArg.isSet())
            setRngStream(streamArg.getValue());

        std::cout << "Begin ... " << std::endl;

        PolyCell::ff_ = fitArg.getValue();

        trackMutations = eventsArg.getValue();
        rampingDrug = rampingArg.getValue();

        if (xfactorArg.isSet())
        {
            X_FACTOR = xfactorArg.getValue();
            if (X_FACTOR <= 0)
            {
                std::cerr << "error: --x-factor argument not positive ("
                          << X_FACTOR << ")\n";
                exit(1);
            }
        }

        if (concentrationArg.isSet())
        {
            DRUG_CONCENTRATION = concentrationArg.getValue();
            if (DRUG_CONCENTRATION < 0)
            {
                std::cerr << "error: --concentration argument < 0 ("
                          << DRUG_CONCENTRATION << ")\n";
                exit(1);
            }
        }

        if (bottleneckSizeArg.isSet() || bottleneckIntervalArg.isSet())
        {
            if (!(bottleneckSizeArg.isSet() && bottleneckIntervalArg.isSet()))
            {
                std::cerr << "error: if one of --bottleneck-size, "
                          << "--bottleneck-interval is set, then both "
                          << "must be set\n.";
                exit(1);
            }

            bottleneck = true;
            bottleneckSize = bottleneckSizeArg.getValue();
            bottleneckInterval = bottleneckIntervalArg.getValue();

            if (bottleneckSize < 1)
            {
                std::cerr << "error: --bottleneck-size < 1 ("
                          << bottleneckSize << ")\n";
                exit(1);
            }
            if (bottleneckInterval < 1)
            {
                std::cerr << "error: --bottleneck-interval < 1 ("
                          << bottleneckInterval << ")\n";
                exit(1);
            }
        }

    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId()
                  << std::endl;
    }

    std::cout << "Opening starting population snapshot ..." << std::endl;
    std::fstream startsnap(startSnapFile.c_str(), std::ios::in|std::ios::binary);
    if (!startsnap.is_open())
    {
        std::cerr << "File could not be open: "<< startSnapFile << std::endl;
        exit(2);
    }
    
    // header
    int Total_Cell_Count;
    double frame_time;
    startsnap.read((char*)(&frame_time), sizeof(double));
    startsnap.read((char*)(&TIME), sizeof(double));
    startsnap.read((char*)(&Total_Cell_Count), sizeof(int));

    sprintf(buffer, "out/%s/snapshots", outDir.c_str());
    std::string outPath = buffer;
    std::cout << "Creating directory " << outPath << " ... "
              << (makePath(outPath) ? "OK" : "failed") << std::endl;

    // CREATE VECTOR WITH populationSize CELLS
    std::cout << "Creating a population of " << populationSize << " cells ..."
              << std::endl;
    PolyCell A(startsnap, genesPath);
    std::vector<PolyCell> Cell_arr(populationSize, A);
    for (auto it = Cell_arr.begin(); it != Cell_arr.end(); ++it)
    {
        it->ch_barcode(getBarcode());
        it->init_gene_stochastic_concentrations();
        it->UpdateRates(0);
    } 
    startsnap.close();

    std::cout << "Saving initial genotype counts ... " << std::endl;
    // save initial population snapshot
    sprintf(buffer, "%s/%s.gen%010d.json", outPath.c_str(),
            outDir.c_str(), generationNumber); 
    saveSnapshot(buffer, Cell_arr, (int)frame_time);

    std::ofstream MUTATIONLOG;
    if (trackMutations)
    {
        // Open MUTATION LOG
        sprintf(buffer, "out/%s/MUTATION_LOG", outDir.c_str());
        MUTATIONLOG.open(buffer);
        if (!MUTATIONLOG.is_open())
        {
            std::cerr << "Mutation log file could not be opened";
            exit(1);
        }
    }
    
    std::cout << "Starting evolution ..." << std::endl;

    // MORAN PROCESS
    while (generationNumber < generationMax)
    {
        std::vector<double> cumulativeProbability(Cell_arr.size() + 1);
        // Probability for death event
        cumulativeProbability[0] = 0.5 * Cell_arr.size();

        // Iterate through cells and gather cumulative fitnesses into vector.
        for (unsigned int i = 1; i < cumulativeProbability.size(); ++i)
        {
            cumulativeProbability[i] = cumulativeProbability[i - 1]
                + Cell_arr[i - 1].fitness();
        }

        double &totalProbability = cumulativeProbability.back();
        double dt = (-1 / totalProbability) * log(randomNumber());
        TIME += dt;

        double r2 = totalProbability * randomNumber();
        unsigned int cellIndex;
        for (cellIndex = 0; cellIndex < cumulativeProbability.size(); ++cellIndex)
        {
            if (r2 < cumulativeProbability[cellIndex])
            {
                break;
            }
        }

        if (cellIndex > 0)
        {
            // reproduction event
            --cellIndex;
            Cell_arr.push_back(Cell_arr[cellIndex]);
            // Possibly mutate both daughter and "parent"
            if (generationNumber >= equilibrationGens)
            {
                // Potentially mutate parent
                PolyCell &cell1 = Cell_arr[cellIndex];
                if (cell1.mrate() * cell1.genome_size() > randomNumber())
                {
                    mutationCount++;
                    if (trackMutations)
                    {
                        cell1.ranmut_Gene(MUTATIONLOG, generationNumber);
                    }
                    else
                    {
                        cell1.ranmut_Gene();
                    }
                }
                // Potentially mutate daughter
                PolyCell &cell2 = Cell_arr.back();
                if (cell2.mrate() * cell2.genome_size() > randomNumber())
                {
                    mutationCount++;
                    if (trackMutations)
                    {
                        cell2.ranmut_Gene(MUTATIONLOG, generationNumber);
                    }
                    else
                    {
                        cell2.ranmut_Gene();
                    }
                }
            }
        }
        else
        {
            // death event
            cellIndex = (int)(randomNumber() * Cell_arr.size());
            Cell_arr.erase(Cell_arr.begin() + cellIndex);
        }

        // Keep population under max limit
        while (Cell_arr.size() > populationMax)
        {
            cellIndex = (int)(randomNumber() * Cell_arr.size());
            Cell_arr.erase(Cell_arr.begin() + cellIndex);
        }

        // Still need to update rates
        for (auto &cell : Cell_arr)
        {
            cell.UpdateRates(dt);
        }

        // update generation counter
        if (TIME - 1 > generationNumber)
        {
            ++generationNumber;

            // save population snapshot every outputFreq generations
            if ((generationNumber % outputFreq) == 0)
            {
                sprintf(buffer, "%s/%s.gen%010d.json", outPath.c_str(),
                        outDir.c_str(), generationNumber); 

                saveSnapshot(buffer, Cell_arr, generationNumber);
                std::cout << "Generation: " << generationNumber << std::endl;
            }

            // sometimes, we save a full dump of population
            if ((snapoutputFreq > 0) && ((generationNumber % snapoutputFreq) == 0))
            {
                sprintf(buffer, "%s/%s.gen%010d.snap", outPath.c_str(),
                        outDir.c_str(), generationNumber); 

                std::fstream OUT2(buffer, std::ios::out | std::ios::binary);
                if (!OUT2.is_open())
                {
                    std::cerr << "Snapshot file could not be opened";
                    exit(1);
                }
      
                double frame_time = generationNumber;
                OUT2.write((char*)(&frame_time),sizeof(double));
                OUT2.write((char*)(&TIME),sizeof(double));
                OUT2.write((char*)(&Total_Cell_Count),sizeof(int));

                int l = 1;
                for (auto k = Cell_arr.begin(); k != Cell_arr.end(); ++k)
                {
                    k->dump(OUT2, l);
                    l++;
                }
                OUT2.close();
            }
        }
     

        if (Cell_arr.size() == 0)
        {
            // extinction!
            std::cout << "Population down to 0. "<< std::endl;
            break;
        }

        if (rampingDrug)
        {
            // Calculate a fitness average
            double fitnessTotal(0);
            for (auto& cell : Cell_arr)
            {
                fitnessTotal += cell.fitness();
            }
            double fitnessAverage = fitnessTotal / (double)Cell_arr.size();

            if (fitnessAverage > 0.5)
            {
                DRUG_CONCENTRATION *= DRUG_INCREASE_FACTOR;
            }
            else
            {
                DRUG_CONCENTRATION /= DRUG_INCREASE_FACTOR;
            }
            // I suppose there could be some checks to make sure
            // DRUG_CONCENTRATION stays in some bounds.
        }
    }

    MUTATIONLOG.close();
    std::cout << "Done." << std::endl;
    std::cout << "Total number of mutation events: " << mutationCount << std::endl;

    return 0;
}
