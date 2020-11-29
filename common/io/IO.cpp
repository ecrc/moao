#ifndef IO_IMPL_
#define IO_IMPL_

#include "IO.hpp"


void IO::read(std::istream & in)
{
    ptr_->read(in);
}

void IO::write(std::ostream & os) const
{
    ptr_->write(os);
}


#endif //IO_IMPL_
