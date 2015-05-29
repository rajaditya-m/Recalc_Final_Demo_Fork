//**************************************************************************************
// Copyright 2004, Huamin Wang.
//**************************************************************************************
// IO_FUNCTIONS
// Only consider little endian here. Not for big endian architecture, like
// SGI, Mac, Sun...
//**************************************************************************************
#ifndef __IO_FUNCTIONS_H__
#define __IO_FUNCTIONS_H__

#include <iostream>


//**************************************************************************************
// Reading Functions.
//**************************************************************************************
template<class T>
inline void Read_Binaries(std::istream &input, T *v, const int size)
{
    input.read((char*)v, sizeof(T)*size);
}

template<class T>
inline void Read_Binary(std::istream &input, T &v)
{
    input.read((char*)&v, sizeof(T));
}

template<class T1,class T2>
inline void Read_Binary(std::istream &input, T1 &v1,T2 &v2)
{
    Read_Binary(input,v1);
    Read_Binary(input,v2);
}

template<class T1,class T2,class T3>
inline void Read_Binary(std::istream &input, T1 &v1,T2 &v2,T3 &v3)
{
    Read_Binary(input,v1);
    Read_Binary(input,v2);
    Read_Binary(input,v3);
}

template<class T1,class T2,class T3,class T4>
inline void Read_Binary(std::istream &input, T1 &v1,T2 &v2,T3 &v3, T4 &v4)
{
    Read_Binary(input,v1);
    Read_Binary(input,v2);
    Read_Binary(input,v3);
    Read_Binary(input,v4);
}


//**************************************************************************************
// Writing Functions.
//**************************************************************************************
template<class T>
inline void Write_Binaries(std::ostream &output,const T *v, const int size)
{
    output.write((const char*)v, sizeof(T)*size);
}

template<class T>
inline void Write_Binary(std::ostream &output,const T& v)
{
    output.write((const char*)&v, sizeof(T));
}

template<class T1,class T2>
inline void Write_Binary(std::ostream &output,const T1& v1,const T2& v2)
{
    Write_Binary(output,v1);
    Write_Binary(output,v2);
}

template<class T1,class T2,class T3>
inline void Write_Binary(std::ostream &output,const T1& v1,const T2& v2,const T3& v3)
{
    Write_Binary(output,v1);
    Write_Binary(output,v2);
    Write_Binary(output,v3);
}

template<class T1,class T2,class T3,class T4>
inline void Write_Binary(std::ostream &output,const T1& v1,const T2& v2,const T3& v3,const T4& v4)
{
    Write_Binary(output,v1);
    Write_Binary(output,v2);
    Write_Binary(output,v3);
    Write_Binary(output,v4);
}

//**************************************************************************************
// Enforced double/float Type Reading Functions.
//**************************************************************************************
template<class T>
inline void Read_Binary_Double(std::istream &input,T &v)
{
    double temp;
    Read_Binary(input,temp);
    v=(T)temp;
}

template<class T>
inline void Read_Binary_Double(std::istream &input,T &v1,T &v2)
{
    Read_Binary_Double(input,v1);
    Read_Binary_Double(input,v2);
}

template<class T>
inline void Read_Binary_Double(std::istream &input,T &v1,T &v2,T &v3)
{
    Read_Binary_Double(input,v1);
    Read_Binary_Double(input,v2);
    Read_Binary_Double(input,v3);
}

template<class T>
inline void Read_Binary_Double(std::istream &input,T &v1,T &v2,T &v3, T &v4)
{
    Read_Binary_Double(input,v1);
    Read_Binary_Double(input,v2);
    Read_Binary_Double(input,v3);
    Read_Binary_Double(input,v4);
}


template<class T>
inline void Read_Binary_Float(std::istream &input,T &v)
{
    float temp;
    Read_Binary(input, temp);
    v=(T)temp;
}

template<class T>
inline void Read_Binary_Float(std::istream &input,T &v1,T& v2)
{
    Read_Binary_Float(input,v1);
    Read_Binary_Float(input,v2);
}

template<class T>
inline void Read_Binary_Float(std::istream &input,T &v1,T &v2,T &v3)
{
    Read_Binary_Float(input,v1);
    Read_Binary_Float(input,v2);
    Read_Binary_Float(input,v3);
}



template<class T>
inline void Write_Binary_Double(std::ostream &output,const T &v)
{
    double temp=(double)v;
    Write_Binary(output,temp);
}

template<class T>
inline void Write_Binary_Double(std::ostream &output,const T &v1,const T &v2)
{
    Write_Binary_Double(output,v1);
    Write_Binary_Double(output,v2);
}

template<class T>
inline void Write_Binary_Double(std::ostream &output,const T &v1,const T &v2,const T &v3)
{
    Write_Binary_Double(output,v1);
    Write_Binary_Double(output,v2);
    Write_Binary_Double(output,v3);
}


template<class T>
inline void Write_Binary_Float(std::ostream &output,const T &v)
{
    float temp=(float)v;
    Write_Binary(output, temp);
}

template<class T>
inline void Write_Binary_Float(std::ostream &output,const T &v1,const T &v2)
{
    Write_Binary_Float(output,v1);
    Write_Binary_Float(output,v2);
}

template<class T>
inline void Write_Binary_Float(std::ostream &output,const T &v1,const T &v2,const T &v3)
{
    Write_Binary_Float(output,v1);
    Write_Binary_Float(output,v2);
    Write_Binary_Float(output,v3);
}

#endif
