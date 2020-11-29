#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <string>
#include <map>
#include <list>

struct Arg{
    std::string key;                             //! argument key
    std::string description;                     //! description of the argument content
    std::string value;                           //! argument value
    Arg();                                  //! default constructor
    Arg(std::string k,std::string d, std::string v="");
};

struct CmdLine{
    std::map<std::string, Arg> p;      //! map of arguments
    std::string usage_;                //! description of the usage
    CmdLine();                         //! 
    CmdLine(std::list<Arg> args, std::string usage="");           
    void add(Arg arg);                   
    void parse(int argc,char **argv);   
    void usage();                       
    void recap(bool values=true);       
    std::string getString(std::string key);       
    int getInt(std::string key);             
    double getReal(std::string key);
};

std::string splitKV(std::string s, std::string &k, std::string &v);

#endif // ARGUMENTS_H
