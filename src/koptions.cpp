// koptions.cpp

#include "koptions.h"
#include <stdlib.h>
#include <vector>

using namespace std;

namespace kutil {
  
  koptions::koptions(int argc, char * argv[]){
    _arglist.clear();
  
    for (int i=1; i<argc; i++){
      _arglist.push_back(string(argv[i]));	
    }
  }

  int koptions::find(const string & optstr){
    int n = optstr.size();
    vector<string>::iterator it;
    
    for (it=_arglist.begin(); it != _arglist.end(); it++){
      if ( it->substr(0,n) == optstr){	
	_arglist.erase(it);
	return PASS;
      }
    }
    return FAIL;
  }

  int koptions::set(const string & optstr, string & value){    
    int n = optstr.size();
    vector<string>::iterator it;
    
    for (it=_arglist.begin(); it != _arglist.end(); it++){
      if ( it->substr(0,n) == optstr){	
	value = it->erase(0,n);
	_arglist.erase(it);
	return PASS;
      }
    }
    return FAIL;
  }

  int koptions::set(const string & optstr, int & value){    
    string strval;
    
    if (set(optstr, strval)){
      value = atoi(strval.c_str());
      return PASS;
    }
    return FAIL;
  }

  int koptions::set(const string & optstr, double & value){    
    string strval;
    
    if (set(optstr, strval)){
      value = atof(strval.c_str());
      return PASS;
    }
    return FAIL;
  }

  int koptions::set_next(const string & optstr, int num, vector<string> & follows){    
    int n = optstr.size();
    vector<string>::iterator it;
    
    for (it=_arglist.begin(); it != _arglist.end(); it++){
      if ( it->substr(0,n) == optstr){
	follows.clear();
	_arglist.erase(it);
	for (int i=0; i<num; i++){
	  if (it == _arglist.end()) return FAIL;	  
	  follows.push_back(*it);
	  _arglist.erase(it);
	}
	return PASS;
      }
    }
    return FAIL;
  }
  
  int koptions::set_next(const string & optstr, string & follows){    
    vector<string> temp;
    if (set_next(optstr,1,temp)){
      follows = temp[0];
      return PASS;
    }
    return FAIL;
  }
  
  int koptions::unrecognized_option(string & value){
    vector<string>::iterator it;
    
    for (it=_arglist.begin(); it != _arglist.end(); it++){	
      if ( (*it)[0] == '-' ){	  
	value = *it;
	_arglist.erase(it);
	return PASS;
      }
    }    
    return FAIL;
  }
  
  int koptions::shift_argument(string & value){
    vector<string>::iterator it;
    
    for (it=_arglist.begin(); it != _arglist.end(); it++){	
      value = *it;
      _arglist.erase(it);
      return PASS;
    }        
    return FAIL;
  }
  
}
