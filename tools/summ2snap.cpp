#include "../src/PolyCell.h"

/*
DESCRIPTION: Converts population summary to a snap file.
*/

int main(int argc, char *argv[]){
    if(argc != 3){
        std::cerr <<"sodasumm <population summary> [ 0-full | 1-single cell ]\n";
        exit(1);
    }

    int flag = atoi(argv[2]);
    assert((flag==0) | (flag==1));

    char buffer[200];
    std::vector <PolyCell> Cell_arr;
    

    // open stream to read population summary
    std::fstream popf(argv[1]);
    if(!popf.is_open()){
        std::cerr << "File could not be opened";
        exit(1);
    }

    std::string line;
    Cell_arr.reserve(POPSIZEMAX);
    while (!popf.eof())
    {
        getline(popf,line);
        std::string word;
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;  
        if ( word == "C" ){
            iss >> word;//cell count
            int count = atoi(word.c_str());   
            iss >> word;//cell files
            std::fstream temp (word.c_str());//convert std::string to char
            if(!temp.is_open()){
                std::cerr << "File could not be open: "<< word <<std::endl;
                exit(1);
            }
            PolyCell A(temp);
            if(flag)
            {
                A.ch_barcode(getBarcode());
                Cell_arr.push_back(A);
                temp.close();
                break;
            }
            for(int i=0; i<count; i++){
                Cell_arr.push_back(A);
                Cell_arr.back().ch_barcode(getBarcode());
            }
            temp.close();
        }
    }
    popf.close();

    //population snapshot
    int Total_Cell_Count = (int)(Cell_arr.size());

    sprintf(buffer,"population.snap");

    // Open stream to write snapshot
    std::fstream OUT(buffer, std::ios::out);
    if ( !OUT.is_open() ) {
      std::cerr << "Snapshot file could not be opened";
      exit(1);
    }
    
    double T0 =0;
    double frame_time =0;
    OUT.write((char*)(&frame_time),sizeof(double));
    OUT.write((char*)(&T0),sizeof(double));
    OUT.write((char*)(&Total_Cell_Count),sizeof(int));

    int l=0;
    for (auto k = Cell_arr.begin(); k != Cell_arr.end(); ++k){
        ++l;
        (*k).dump(OUT,l);
    }
    OUT.close();
}
