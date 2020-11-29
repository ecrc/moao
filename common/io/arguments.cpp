#ifndef ARGUMENTS_IMPL_
#define ARGUMENTS_IMPL_

#include "arguments.hpp"

#include <iostream>
#include <stdexcept>

Arg::Arg():Arg("","",""){}
Arg::Arg(std::string k,std::string d, std::string v):key(k),description(d),value(v){
    //TODO remove white spaces?
}



CmdLine::CmdLine(){
    CmdLine(std::list<Arg>());
}

/*!\brief constructor
 *
 * initialize the map with a list of arguments
 * @param[in] args : list<Arg> : list of arguments
 */
CmdLine::CmdLine(std::list<Arg> args, std::string usage):usage_(usage){
    for(std::list<Arg>::iterator it=args.begin();it!=args.end();it++){
        add(*it);
    }
}

/*!\brief Add an option to the option map
 *
 * @param[in] arg : Arg : argument to add
 */
void CmdLine::add(Arg arg){
    std::string k=std::string("--"+arg.key+"=");
    p[k]=arg;
}

/*!\brief parse the command line
 *
 * @param[in] argc : int    : number of argument in the string
 * @param[in] argv : char** : string to parse
 */
void CmdLine::parse(int argc,char **argv){
    std::string err="";
    for(int i=1;i<argc;i++){
        std::string s=argv[i];
        std::string k;
        std::string v;
        err+=splitKV(s,k,v);
        if(p.find(k)== p.end() && !k.empty())
            {err+="unknown argument: "+k+"\n";}
        else
            {p[k].value=v;}
    }
    if(!err.empty()){
        usage();
        throw std::runtime_error("\n"+err);
    }
}
/*!\brief display the usage and list the options
 */
void CmdLine::usage(){
    recap(false);
    std::cout<<"Usage:\n"<<usage_;
    std::cout<<"--------------------\n"<<std::endl;
}

/*!\brief summarize the options
 *
 * display the list of options available and additional information:
 * either the values or the description of the arguments
 *@param values : bool : display values if set to true, the description otherwise
 */
void CmdLine::recap(bool values){
    if(values)
        {std::cout<<"Using:";}
    else
        {std::cout<<"Available arguments:";}
    std::cout<<"\n";
    for(std::map<std::string,Arg>::iterator it=p.begin();it!=p.end();it++){
        std::string  second;
        if(values)
            {second = it->second.value;}
        else
            {second = it->second.description;}
        std::cout<<it->second.key<<" : "<< second<<"\n";
    }
    std::cout<<"--------------------\n"<<std::endl;

}

/*!\brief retieve an the argument value
 *
 * argument are refer to by their key 
 * @param key : string : argument key
 * @return    : string : argument value
 */
std::string CmdLine::getString(std::string key){
    key="--"+key+"=";
    //does key exists?
    if(p.find(key)== p.end()){
        std::cerr<<"reading unknown option: "<<key<<std::endl;;
        return "";
    }
    //return value
    return p[key].value;
}

/*!\brief retieve an the argument value
 *
 * argument are refer to by their key the value is cast from string to int
 * warning: no exception is handle
 * @param key : string : argument key
 * @return    : int    : argument value
 */
int CmdLine::getInt(std::string key){
    key="--"+key+"=";
    //does key exists?
    if(p.find(key)== p.end()){
        std::cerr<<"reading unknown option: "<<key<<std::endl;;
        return 0;
    }
    //return value
    return atoi(p[key].value.c_str());
}

/*!\brief retieve an the argument value
 *
 * argument are refer to by their key the value is cast from string to double
 * warning: no exception is handle
 * @param key : string : argument key
 * @return    : double : argument value
 */
double CmdLine::getReal(std::string key){
    key="--"+key+"=";
    //does key exists?
    if(p.find(key)== p.end()){
        std::cerr<<"reading unknown option: "<<key<<std::endl;;
        return 0.;
    }
    //return value
    return atof(p[key].value.c_str());
}

/*!\brief split a string in two strings: key and value
 *
 * @param[in]  s    : string : input string
 * @param[out] k    : string : key starts with "--" and end with "="
 * @param[out] v    : string : value substring from the key
 */
std::string splitKV(std::string s, std::string &k, std::string &v){
    //locate key (start and end
    int key_s = s.find("--");
    int key_e = s.find("=")+1;
    //delimiters found?
    std::string err="";
    if(key_s<0){
        err+="\nmissing '--' delimiter for argument :"+s+".";
    }
    if(key_e<1){
        err+="\nmissing '='  delimiter for argument :"+s+".";
    }
    if(err!="")
        {return err+"\n";}
    //extract key and value
    k=s.substr(key_s,key_e);
    v=s.substr(key_e,s.size());
    return err;
}
#endif // ARGUMENTS_IMPL_
