#include <iostream>
#include <cstdlib>

#include "settings.h"
#include "computation.h"

int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }
  
  // get path of parameter file
  std::string parameterFile = argv[1];
  
  
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Parsing paramter file..." << std::endl;
  std::cout << std::endl;
  
  
  // construct computation obj: parses parameter file and prints settings
  Computation computation = Computation();
  computation.initialize(argc, argv);
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Starting driven cavitiy simulation ..." << std::endl;
  std::cout << std::endl;
  
  // start simulating
  computation.runSimulation();
  
  return EXIT_SUCCESS;
}
