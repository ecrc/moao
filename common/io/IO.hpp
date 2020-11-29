#ifndef IO_H
#define IO_H

#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <memory>

 /*! \brief Open a file 
 *
 * throw exception if the file could not be opened
 * param fn file name
 * param om open mode
 */
template<typename T>
T openFile(std::string fn,std::ios_base::openmode om);

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v);
template<typename T>
std::istream& operator>>(std::istream& in, std::vector<T>& v);

/*!\brief basic IO interface 
 */
struct IOInterface
{
    virtual void read(std::istream & in) = 0;
    virtual void write(std::ostream & os) const = 0;
    virtual ~IOInterface() = default;
};

/*!\brief IO interface implementation
 */
template<typename T>
struct IOImpl : IOInterface
{
    T t_;
    IOImpl(T t);
    void read(std::istream & in) override;
    void write(std::ostream & os) const override;
};

/*!
 */
struct IO
{
    std::unique_ptr<IOInterface> ptr_;
    template<typename T>
    IO(T t);
    void read(std::istream & in);
    void write(std::ostream & os) const;
};

template<typename T>
void readIO(std::string filename, T & t);

template<typename T>
void writeIO(std::string filename, T & t);

#include "IO.cpp.in"

#endif //IO_H
