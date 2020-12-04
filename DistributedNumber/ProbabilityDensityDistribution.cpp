/*
 * 
 * MIT License
 *  
 * Copyright (c) 2020 Michael von Mengershausen
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * 
 */

#include "ProbabilityDensityDistribution.h"
////////////////////////////////////////////
// construction / destruction
////////////////////////////////////////////
CProbabilityDensityDistribution::CProbabilityDensityDistribution()
{
    _Init();
}

CProbabilityDensityDistribution::CProbabilityDensityDistribution(const CProbabilityDensityDistribution& other)
{
    operator=(other);
}

CProbabilityDensityDistribution::CProbabilityDensityDistribution(const CDistribution& other)
{
    operator=(other);
}

CProbabilityDensityDistribution::~CProbabilityDensityDistribution()
{
    _Init();
}
////////////////////////////////////////////
// assignment operators
////////////////////////////////////////////
CProbabilityDensityDistribution& CProbabilityDensityDistribution::operator=(const CProbabilityDensityDistribution& other)
{
    CDistribution::operator=(other);
    bNormalized = other.bNormalized;
    dfAbsIntegral = other.dfAbsIntegral;
    return *this;
}
CProbabilityDensityDistribution& CProbabilityDensityDistribution::operator=(const CDistribution& other)
{
    _Init();
    CDistribution::operator=(other);
    return *this;
}
////////////////////////////////////////////
// binary operator
////////////////////////////////////////////
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator+=(const CDigFloat& Value)
{
    // _DeNormalize
    _DeNormalize();
    
    CDistribution::operator+=(Value);
    
    // normalize again
    _Normalize();
    
    return *this;
}
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator-=(const CDigFloat& Value)
{
    // _DeNormalize
    _DeNormalize();
    
    CDistribution::operator-=(Value);
    
    // normalize again
    _Normalize();
    
    return *this;
}
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator*=(const CDigFloat& Value)
{
    // simply do nothing: due to normalization the factor will be cancelled    
    return *this;
}
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator/=(const CDigFloat& Value)
{
    // simply do nothing: due to normalization the divisor will be cancelled
    return *this;
}

////////////////////////////////////////////
// comparison operators
////////////////////////////////////////////
bool CProbabilityDensityDistribution::operator==(const CProbabilityDensityDistribution& other) const
{
    // at the moment the distribution is the relevant part:
    // dfAbsIntegral can be calculated any time and is the same if the distributions are the same
    // bNormalized is hidden from user and will be performed when necessary
    return CDistribution::operator==(other);
}

bool CProbabilityDensityDistribution::operator!=(const CProbabilityDensityDistribution& other) const
{
    return !operator==(other);
}
///////////////////////////////////
// function on all distribution variables
///////////////////////////////////    
void CProbabilityDensityDistribution::Shift(const CDigFloat& shift)
{
    // this can simply applied as the normalization is not affected 
    CDistribution::Shift(shift);
    
}
void CProbabilityDensityDistribution::Scale(const CDigFloat& scale)
{
    // scaling needs denormalization and normalization 
    _DeNormalize();
    
    CDistribution::Scale(scale);
    
    _Normalize();
}


////////////////////////////////////////////
// functions
////////////////////////////////////////////
CDigFloat CProbabilityDensityDistribution::DistriValue(const CDigFloat& variable)
{
    if(!bNormalized)
        _Normalize();
    return CDistribution::DistriValue(variable);
}
CDigFloat CProbabilityDensityDistribution::AbsIntegral(const CDigFloat& variableLeft, const CDigFloat& variableRight, int nthOrder)
{
    if(!bNormalized)
        _Normalize();
    return CDistribution::AbsIntegral(variableLeft,variableRight,nthOrder);
}
CDigFloat CProbabilityDensityDistribution::AbsIntegral(const int& nthOrder)
{
    if(!bNormalized)
        _Normalize();
    return CDistribution::AbsIntegral(nthOrder);
}
void CProbabilityDensityDistribution::GetIntervall(const CDigFloat& variable, M_DFDF::const_iterator& VariableLeft, M_DFDF::const_iterator& VariableRight)
{
    if(!bNormalized)
        _Normalize();
    return CDistribution::GetIntervall(variable,VariableLeft,VariableRight);
}
CDigFloat CProbabilityDensityDistribution::Mean(int nthOrder)
{
    if(!bNormalized)
        _Normalize();
    return CDistribution::Mean(nthOrder);
}
void CProbabilityDensityDistribution::Reset()
{
    _Init();
}
bool CProbabilityDensityDistribution::Add(const CDigFloat& variable, const CDigFloat value, bool bLastElement /*== false*/)
{
    // first check if the distribution is already normalized
    if(bNormalized)
        _DeNormalize();
    
    // add the new value
    bool bSuccess = CDistribution::Add(variable,value);
    
    // check if this is the last element: if so --> normalize
    if(bLastElement)
        _Normalize();
    
    return bSuccess;
}

void CProbabilityDensityDistribution::_DeNormalize()
{
    // do not normalize again ... only once
    if(!bNormalized)
        return;
 
    // revert normalization
    CDistribution::operator*=(dfAbsIntegral);
    bNormalized = false;
    
}

void CProbabilityDensityDistribution::_Normalize()
{
    // do not normalize again ... only once
    if(bNormalized)
        return;
    
    dfAbsIntegral = CDistribution::AbsIntegral();
    CDistribution::operator/=(dfAbsIntegral);
    
    // remember we have been normalized
    bNormalized = true;
    
}
void CProbabilityDensityDistribution::_Init()
{
    dfAbsIntegral = 0;
    bNormalized = false;
    CDistribution::_Init();
}

