#include <iostream>
#include <cstdlib>
#include <typeinfo>

#include "settings.h"
#include "computation/1_computation_parallel.h"

int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }
  std::cout << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Parsing paramter file..." << std::endl;
  std::cout << std::endl;

  // start MPI
  MPI_Init(&argc, &argv);
  
  // construct computation obj: parses parameter file and prints settings
  ComputationParallel computation = ComputationParallel();
  std::cout << typeid(computation).name() << std::endl;

  computation.initialize(argc, argv);
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Starting parallel simulation ..."  << std::endl;
  std::cout << std::endl;
  
  // // start simulating
  computation.runSimulation();

  std::cout << "Back to main" << std::endl;
  MPI_Finalize();
  std::cout << "-----------------------------" << std::endl;
  std::cout << "Finalized MPI" << std::endl;

  // test some stuff
  //computation.runTest();

  std::cout << "-------------------------------------------------" << std::endl;
  
  return EXIT_SUCCESS;
}
