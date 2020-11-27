/*
 * Copyright (c) 2020 <copyright holder> <email>
 * 
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "CartesianCoordinate.h"

CCartesianCoordinate::CCartesianCoordinate()
{
    // do nothing: default constructor delivers no coordinates   
        
}

CCartesianCoordinate::CCartesianCoordinate(const CCartesianCoordinate& other)
{
    (*this)=other;

}

CCartesianCoordinate::~CCartesianCoordinate()
{
    clear();
}

///////////// arithmetic operators ///////////

CCartesianCoordinate& CCartesianCoordinate::operator*=(const CDigFloat& other)
{
    for(auto icoord = begin(); icoord != end(); icoord++)
        (*icoord)*=other;
    
    return *this;
}

CCartesianCoordinate CCartesianCoordinate::operator*(const CDigFloat& other)
{
    CCartesianCoordinate ccResult = *this;
    ccResult *= other;
    return ccResult;
}

CCartesianCoordinate & CCartesianCoordinate::operator/=(const CDigFloat& other)
{
    for(auto icoord = begin(); icoord != end(); icoord++)
        (*icoord)/=other;
    
    return *this;
}

CCartesianCoordinate CCartesianCoordinate::operator/(const CDigFloat& other)
{
    CCartesianCoordinate ccResult = *this;
    ccResult /= other;
    return ccResult;
}

bool CCartesianCoordinate::operator==(const CCartesianCoordinate& other) const
{
    // sizes (no. of dimensions) must be the same
    if(other.size() != size())
        return false;
 
    // iterate: return false immedeately in case of inequality
    for(int icoord = 0; icoord < size(); icoord++)
        if((*this)[icoord] != other[icoord])
            return false;
        
    return true;
}

bool CCartesianCoordinate::operator!=(const CCartesianCoordinate& other) const
{
    return !operator==(other);
}

bool CCartesianCoordinate::operator<(const CCartesianCoordinate& other) const
{
    // iterate backward: highest coordinate component rules
    for(int icoord = size(); icoord >=0 ; icoord--)
        if((*this)[icoord] < other[icoord])
            return true;
        
    return false;
}

bool CCartesianCoordinate::operator<=(const CCartesianCoordinate& other) const
{
   
    return (*this)<other || (*this) == other;
}

bool CCartesianCoordinate::operator>(const CCartesianCoordinate& other) const
{
    // iterate backward: highest coordinate component rules
    for(int icoord = size(); icoord >=0 ; icoord--)
        if((*this)[icoord] > other[icoord])
            return true;
        
    return false;
}   

bool CCartesianCoordinate::operator>=(const CCartesianCoordinate& other) const
{
   
    return (*this) > other || (*this) == other;
}   
/////////////////////////////////////
// functions
/////////////////////////////////////
CDigFloat CCartesianCoordinate::Distance()
{
    CDigFloat dfDistance = 0;
    for(auto icoord = begin(); icoord != end(); icoord++)
        dfDistance += (pow(*icoord,2));
    
    dfDistance = sqrt(dfDistance);
    
    return dfDistance;
}

std::__cxx11::string CCartesianCoordinate::Print(bool bWithError /*= false*/)
{
    ostringstream oss;
    
    // ostringstream cant handle negative precision
    auto icoord = begin();
    oss << "(" << icoord->Print(bWithError);
    for(icoord = begin()+1; icoord != end(); icoord++)
        oss << "," << icoord->Print(bWithError);
    oss << ")";
    return oss.str();
}


/////////////////////////////////////
// external functions
/////////////////////////////////////
CDigFloat ScalarProduct(CCartesianCoordinate& one, CCartesianCoordinate& other)
{
    CDigFloat dfScalarProduct;
    
    // check sizes: return invalid number in case the sizes do not fit
    if(other.size() != one.size())
        return dfScalarProduct;
    
    dfScalarProduct = 0;
    for(int icoord = 0; icoord < one.size(); icoord++)
        dfScalarProduct += one[icoord] *  other[icoord];
    
    return sqrt(dfScalarProduct);
}

CCartesianCoordinate VectorProduct(CCartesianCoordinate& one, CCartesianCoordinate& other)
{
       CCartesianCoordinate ccVectorProduct;
    
    // check sizes: return invalid number in case the sizes do not fit
    if(other.size() != one.size())
        return ccVectorProduct;
    
    int nSize = one.size();
    for(int icoord = 0; icoord < nSize; icoord++)
    {
        int nFirst = (icoord+1)%nSize;
        int nSecond =(icoord+2)%nSize;
        ccVectorProduct.push_back(one[nFirst] *  other[nSecond] - one[nSecond] * other[nFirst]);
    }
    
    return ccVectorProduct;
}
