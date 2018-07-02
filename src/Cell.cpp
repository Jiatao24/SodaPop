#include "Cell.h"

Cell::Cell() : barcode_(getBarcode()),
	       ID_(0),
	       o_mrate_(0),
	       c_mrate_(0)
{
    Gene_L_.reserve(GENECOUNTMAX);
    Gene_arr_.reserve(GENECOUNTMAX);
}


// Construct from cell file (text description).
Cell::Cell(std::fstream& cell_in)
{
    char buffer[140];
    Gene_L_.reserve(GENECOUNTMAX);
    Gene_arr_.reserve(GENECOUNTMAX);
    ch_barcode(getBarcode());
    std::string line;
    std::string genesPath = "files/genes/";
    while (!cell_in.eof())
    {
        getline(cell_in,line);
        std::string word;
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;
        if (word == "genes_path") {
	    iss >> word; 
	    genesPath = word.c_str();
        }
        if (word == "org_id") {
	    iss >> word; 
	    ID_ = atoi(word.c_str());
        }
        else if (word == "mrate" ) {
	    iss >> word; 
	    o_mrate_ = atof(word.c_str());
	    c_mrate_ = atof(word.c_str());
        }       
        else if (word == "G") {
        // Reading gene files; 
        // Concentration and stability from gene file take precedence.
              iss >> word;
              // Open gene file.
              sprintf(buffer, "%s%s.gene", genesPath.c_str(), word.c_str());
              std::fstream temp (buffer);
              if (!temp.is_open()) {
		  std::cerr << "File could not be open: " << buffer << std::endl;
		  exit(2);
              }
              Gene_arr_.emplace_back(temp);
              std::cout << "Inserted: "<< word << std::endl;
              
              // Check if gene is correctly inserted.
              std::vector<Gene>::iterator i = Gene_arr_.end();
              i--;
              std::cout << (*i).nseq() << std::endl;        
              temp.close();
        }
    }
}


// Constructs from a unit cell stored in binary.
Cell::Cell(std::fstream& IN, const std::string& genesPath)
{
    Gene_L_.reserve(GENECOUNTMAX);
    Gene_arr_.reserve(GENECOUNTMAX);

    char buffer[140];
    int cell_id, cell_index, gene_size;
    double m,m0;
   
    IN.read((char*)(&cell_index),sizeof(int));  
    IN.read((char*)(&cell_id),sizeof(int));
    ID_ = cell_id;

    //read barcode
    int l;
    IN.read((char*)&l, sizeof(int));
    //construct vector container with nl elements
    std::vector<char> buf(l);
    IN.read(&buf[0], l);
    ch_barcode(std::string().assign(buf.begin(), buf.end()));

    IN.read((char*)(&m0),sizeof(double));
    o_mrate_ = m0;
    IN.read((char*)(&m),sizeof(double));
    c_mrate_ = m;
    IN.read((char*)(&gene_size),sizeof(int));
    
    //read gene info
    for (int j = 0; j < gene_size; j++)
    {
        double e, c, dg, f, stochastic_concentration;
        int gene_nid, Ns, Na;
        std::string DNAsequence;
     
        IN.read((char*)(&gene_nid),sizeof(int));   
        IN.read((char*)(&e),sizeof(double));
        IN.read((char*)(&c),sizeof(double));
        IN.read((char*)(&stochastic_concentration),sizeof(double));
        IN.read((char*)(&dg),sizeof(double));
        IN.read((char*)(&f),sizeof(double));
        IN.read((char*)(&Ns),sizeof(int));
        IN.read((char*)(&Na),sizeof(int));

        //read DNA sequence
        int nl;
        IN.read((char*)&nl, sizeof(int));
        //construct vector container with nl elements
        std::vector<char> buf(nl);
        IN.read(&buf[0], nl);
        DNAsequence.assign(buf.begin(), buf.end());
     
        sprintf(buffer, "%s%d.gene", genesPath.c_str(), gene_nid);
        std::fstream gene_data(buffer);
        if (!gene_data.is_open())
        {
            std::cerr << "ERROR: Cannot open gene file " << buffer << std::endl;
            exit(2);
        }
        Gene G(gene_data);   
        //update gene information
        dg = exp(-dg/kT);
        // Ideally we should have a Gene constructor instead of
        //   setters for these private variables.
        G.ch_dg(dg);
        G.ch_conc(c);
        G.ch_f(f);
        G.Update_Sequences(DNAsequence);
        G.ch_Na(Na);
        G.ch_Ns(Ns); 
        Gene_arr_.push_back(G);
    }
}


// Initialize the cumulative gene length array.
void Cell::FillGene_L()
{
    int sum = 0;
    for (auto it = Gene_arr_.begin(); it != Gene_arr_.end(); ++it)
    {
        sum += it->length();
        Gene_L_.push_back(sum);
    }
}

// Return total mutation count
// spec:
//  - 0, Ns+Na
//  - 1, Ns
//  - 2. Na
int Cell::total_mutations(const int& spec)
{
    assert( (spec < 3) &&  (spec >= 0) );

    int sa = 0;
    int s = 0;
    int a = 0;

    for (auto i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i)
    {
        int Ns = i->Ns();
        int Na = i->Na();
        s += Ns;
        a += Na;
        sa += (Ns+Na);
    }

    if(spec == 0) return sa;
    else if (spec == 1) return s;
    else return a;
}


void Cell::init_gene_stochastic_concentrations()
{
    for (auto& gene : Gene_arr_)
    {
        gene.init_stochastic_conc();
    }
}
