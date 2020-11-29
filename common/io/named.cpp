#ifndef NAMED_IMPL_
#define NAMED_IMPL_

#include "named.hpp"
#include <iostream>

template<typename T>
Named<T>::Named(T & t, std::string & err, std::string name):
            t_(t), err_ (err), name_(std::move(name))
            {}

template<typename T>
Named<T> named(T & t, std::string & err, std::string name)
{
    return Named<T>(t, err, std::move(name));
}

template<typename T>
std::istream& operator>>(std::istream& in, Named<T> & n)
{
    std::string s;
    std::getline(in, s);
    if(s.compare(0,n.name_.size(), n.name_ ) ){
        n.err_.append("Line does not start with: "+n.name_ +"\n(got:"+s+")\n\n");
    }
    in >> n.t_;
    in.get();

    return in;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Named<T> & n)
{    
    os<<n.name_<<'\n'<<n.t_<<'\n';
    return os;
}


template<typename T>
IO namedIO(T & t, std::string & err, std::string name)
{
    return IO(Named<T>(t, err, std::move(name)));
}

#endif // NAMED_IMPL_
