#include "../src/global.h"
#include "../src/rng.hpp"

// Output into ASCII:
// Header: frame_time, t0, total_cell_count
//   - frame_time :: generation count
//   - t0 :: starting time. I think it should be 0

int main(int argc, char *argv[]){
    if(argc != 4){
        std::cerr <<"sodasnap <snap-binary> <out-ascii> [ 0-full | 1-minimal ]\n";
        exit(1);
    }

    int flag = atoi(argv[3]);
    assert((flag==0) | (flag==1));

    initialize_rng();

    int Total_Cell_Count;
    double frame_time, T0;

    //open binary file 
    std::fstream IN(argv[1], std::ios::in|std::ios::binary);
    if (!IN.is_open()){
        std::cerr << "Binary file could not be opened.\n";
        exit(1);
    } 

    std::fstream OUT(argv[2], std::ios::out);
    if (!OUT.is_open()){
        std::cerr << "Ascii file could not be opened.\n";
        exit(1);
    }

    //Read header
    IN.read((char*)(&frame_time),sizeof(double));
    IN.read((char*)(&T0),sizeof(double));
    IN.read((char*)(&Total_Cell_Count),sizeof(int));
    OUT << frame_time << " " << T0 << " " << Total_Cell_Count <<  std::endl;

    //Read Cell array
    for(int i=0; i < Total_Cell_Count; i++){
      OUT << std::endl;
      if(flag)
      {
          qread_Cell(IN,OUT);//numbers are DUMMY time parameters
      }
      else
      {
          read_Cell(IN,OUT);
      }
    }
    IN.close();
    OUT.close();  
}

