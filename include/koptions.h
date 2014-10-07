// koptions.h
//
// author:  Michael Mulhearn
//
// The k in koptions stands for "Keep it simple, stupid!"
// 
// FEATURES:
// 
// A simple command line option parsing utility..
// 
// INSTALLATION:
// 
// Use the kutil library.  Or include header "koptions.h" and link
// "koptions.cpp".
// 
// USAGE:
// 
// int main(int argc, char * argv[]){
//    char * opt;
//    char * file;
//    int max;
//    
//    koptions options(argc, argv);
//    
//    //check for the --help option:
//    if ( options.find("--help") ) { usage(); }
//    
//    //check for the --max=<MAX> option, setting a variable if found:    
//    options.set("--max", max);
//    
//    //check for the --verbose option:
//    int verbose = options.find("--verbose");
// 
//    //check for unrecognized options (beginning with -- or -)
//    //  and handle this however you like.
//    while(opt = options.unrecognized()){
//       printf("WARNING: unrecognized option %s\n", opt)
//    } 
//     
//    //assume remaining options are files:
//    while(file = options.shift()){
//       //...
//    }
//    //...
// }
// 

#ifndef __KOPTIONS_H__
#define __KOPTIONS_H__

#include <string>
#include <vector>

namespace kutil {
  
  class koptions{   
  public:     
    //return convention:
    enum {PASS=1, FAIL=0, DONE=0};
    
    koptions(int argc, char * argv[]); 
    
    //find:  was the option asserted?
    int find (const std::string & optstr);
    
    //set:  if option asserted, set value to trailing characters.
    int set  (const std::string & optstr, std::string &  value);
    int set  (const std::string & optstr, int  &         value);
    int set  (const std::string & optstr, double &       value);
    
    //set_next:  if option asserted, set value to the next option(s):
    int set_next (const std::string & optstr, std::string & follows);
    int set_next (const std::string & optstr, int num, std::vector<std::string> & follows);
    
    //get remaining arguments that beginning with "-" or "--", 
    //  return 0 when finished.
    int unrecognized_option (std::string & value);   
    
    //size:  how many arguments remain in the list?
    int size() { return _arglist.size(); }  
    
    //get remaining arguments (even if they are options!)
    //  return 0 when finished.
    int shift_argument (std::string & value);
    
  private:
    std::vector<std::string> _arglist;    
  };

}

#endif
