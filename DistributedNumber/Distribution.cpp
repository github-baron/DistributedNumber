
/*
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
 */



#include "Distribution.h"

// DEBUG
#include <iostream>

CDistribution::CDistribution()
{
    _Init();
}

CDistribution::CDistribution(const CDistribution& other)
{
    mDistribution = other.mDistribution;
}

CDistribution::~CDistribution()
{
    mDistribution.clear();
}

CDistribution& CDistribution::operator=(const CDistribution& other)
{
    mDistribution = other.mDistribution;
    return *this;
}

bool CDistribution::operator==(const CDistribution& other) const
{
    return mDistribution == other.mDistribution;
}

bool CDistribution::operator!=(const CDistribution& other) const
{
    return mDistribution == other.mDistribution;
}

    ///////////////////////////////////
    // binary operators
    /////////////////////////////////// 
CDistribution & CDistribution::operator*=(const CDigFloat& Value)
{
    for(auto iel = mDistribution.begin(); iel != mDistribution.end(); iel++) 
        iel->second*=Value; 
    return *this;
}
CDistribution & CDistribution::operator/=(const CDigFloat& Value)
{
    for(auto iel = mDistribution.begin(); iel != mDistribution.end(); iel++) 
        iel->second/=Value; 
    return *this;
}

CDistribution & CDistribution::operator+=(const CDigFloat& Value)
{
    for(auto iel = mDistribution.begin(); iel != mDistribution.end(); iel++) 
        iel->second+=Value; 
    return *this;
}
CDistribution & CDistribution::operator-=(const CDigFloat& Value)
{
    for(auto iel = mDistribution.begin(); iel != mDistribution.end(); iel++) 
        iel->second-=Value; 
    return *this;
}

///////////////////////////////////
// function on all distribution variables
///////////////////////////////////    
void CDistribution::Shift(const CDigFloat& shift)
{
    // generate map to copy from
    M_DFDF ShiftedDistribution;
    for(auto iel: mDistribution)
        ShiftedDistribution[CDigFloat(shift)+iel.first] = iel.second;
    
    // copy
    mDistribution.clear();
    for(auto iel: ShiftedDistribution)
        mDistribution[iel.first] = iel.second;
    
    // clear temporary map
    ShiftedDistribution.clear();    
}

void CDistribution::Scale(const CDigFloat& scale)
{
    // generate map to copy from
    M_DFDF ShiftedDistribution;
    for(auto iel: mDistribution)
        ShiftedDistribution[CDigFloat(scale)*iel.first] = iel.second;
    
    // copy
    mDistribution.clear();
    for(auto iel: ShiftedDistribution)
        mDistribution[iel.first] = iel.second;
    
    // clear temporary map
    ShiftedDistribution.clear();   
}
////////////////////////////////////////
// functions
////////////////////////////////////////
bool CDistribution::Add(const CDigFloat& variable, const CDigFloat value)
{   
    // negative values are not allowed
    if( value < 0 )
        return false;
    
    // now add the value
    mDistribution[variable] = value;
    return true;
}

void CDistribution::GetIntervall(const CDigFloat& variable, M_DFDF::const_iterator& VariableLeft, M_DFDF::const_iterator& VariableRight)
{
    
    int nIntervall = mDistribution.size()-1;    
    VariableLeft = mDistribution.begin();
    
    
    // set imax to the end
    VariableRight = VariableLeft;
    
    // if the map is of size zero: leave now
    if(nIntervall < 0)
        return;
    
    advance( VariableRight ,nIntervall);
 
    // if variable is left outside ...
    if(variable < VariableLeft->first)
    {
        // ... set right variable and left variable to begin() and leave
        VariableRight = VariableLeft;
        return;
    }        
        
    // if variable is right outside ...
    if(variable > VariableRight->first)
    {
        // ... set right variable and left variable to end() and leave
        VariableRight = mDistribution.end();
        VariableLeft = VariableRight;
        return;
    }
    
    // now start binary search:
    M_DFDF::const_iterator iact;
    while( nIntervall > 1 )
    {
        // halfen the intervall
        nIntervall = (nIntervall+1)/2;
        
        // set actual value to new middle
        iact = VariableLeft;        
        advance(iact, nIntervall);
        
        // set new intervall limits
        if(variable <= iact->first)
            VariableRight = iact;
        else        
            VariableLeft = iact;
            
    }   
}


CDigFloat CDistribution::DistriValue(const CDigFloat& variable)
{
    // init result with zero
    CDigFloat dfDistriValue(0);
    
    M_DFDF::const_iterator imin, imax;
    GetIntervall(variable, imin, imax);
    
    // in case variable is left out of bounds: GetIntervall did not find a matching intervall --> imin = mDistribution.begin()
    if(variable < imin->first)
        return dfDistriValue;
    
    // in case the variable is out of bound on the right: imin == imax == mDistribution.end()
    if(imin ==imax && imax == mDistribution.end())
        return dfDistriValue;
    
    // now we got the min max --> get the linearily interpolated value
    dfDistriValue = CDigFloat(imin->second) + (CDigFloat(imax->second) - imin->second) / (CDigFloat(imax->first) - imin->first) * (CDigFloat(variable) - imin->first);
 
    return dfDistriValue;
    
}

CDigFloat CDistribution::AbsIntegral(const CDigFloat& variableLeft, const CDigFloat& variableRight, int nthOrder /*= 0*/)
{
    M_DFDF::const_iterator iminLeft, imaxLeft, iminRight, imaxRight;
    GetIntervall(variableLeft, iminLeft, imaxLeft);
    GetIntervall(variableRight, iminRight, imaxRight);
    
    // check for outside value left
    CDigFloat variableLeftCorr = variableLeft;
    if(iminLeft == imaxLeft)
        variableLeftCorr = iminLeft->first;
    
    // check for outside value left: < begin()
    CDigFloat variableRightCorr = variableRight;
    if(iminRight == imaxRight && imaxRight == mDistribution.end())
        variableRightCorr = prev(mDistribution.end())->first;
    
    // DEBUG
//     cout << endl << "AbsIntegral(" << variableLeft.RawPrint(10) << "," << variableRight.RawPrint(10) << "): " << endl;
//     cout << "left intervall: " << iminLeft->first.RawPrint(10) << "," << imaxLeft->first.RawPrint(10) << endl;
//     cout << "right intervall: " << iminRight->first.RawPrint(10) << "," << imaxRight->first.RawPrint(10) << endl;
//     cout << "left corr. variable: " << variableLeftCorr.RawPrint(10)  << endl;
//     cout << "right corr. variable: " << variableRightCorr.RawPrint(10) << endl;
    
    // check for same intervall to avoid double calculations
    // in this case no calculation is done for 
    // complete elements
    // right part
    bool bSameIntervall = (imaxLeft->first == imaxRight->first);
    
    // set the limit for the left integral:
    CDigFloat dfLeftIntegralMax = min(imaxLeft->first.RawValue(), variableRightCorr.RawValue());
    
    // init variable for result
    CDigFloat dfIntegral(0);    
    
    // calculate the integral of complete elements: only if the
    // variable arguments of this function are NOT in the same intervall
    if( !bSameIntervall )
        for(auto iel = imaxLeft; iel!=iminRight; iel++)
            dfIntegral += _IntegralConsecutiveElements(iel,nthOrder);
    
    // DEBUG    
//     cout << "Integral (1): " << dfIntegral.RawPrint(10) << endl;
    // add the part on the left:
    dfIntegral += _IntegralOfTwoPoints( variableLeftCorr, DistriValue(variableLeftCorr), dfLeftIntegralMax, DistriValue(dfLeftIntegralMax),nthOrder);
    
    // DEBUG    
//  cout << "Integral (2): " << dfIntegral.RawPrint(10) << endl;
    // add the part on the right:
    if(iminRight != mDistribution.end() && !bSameIntervall)
        dfIntegral += _IntegralOfTwoPoints(iminRight->first, iminRight->second, variableRightCorr, DistriValue(variableRightCorr),nthOrder);
    
    // DEBUG    
//  cout << "Integral (3): " << dfIntegral.RawPrint(10) << endl;
    return dfIntegral;
}

 CDigFloat CDistribution::Median(int nthOrder)
{
    CDigFloat dfMedian;
    if(!CoverageFromTo(50,mDistribution.begin()->first,dfMedian,false,nthOrder))
        dfMedian = myNAN;    
    
    return dfMedian;
}

bool CDistribution::CoverageFromTo(const CDigFloat& dfCoveragePercent, const CDigFloat& dfFrom, CDigFloat& dfTo, bool bReverse, int nthOrder)
{
    CDigFloat dfStart = dfFrom;
    dfTo = bReverse ? mDistribution.begin()->first : prev(mDistribution.end())->first;
 /*   cout << "CoverageFromTo(" << dfCoveragePercent.RawPrint(10) << ", " << dfFrom.RawPrint(10) << ", " << dfTo.RawPrint(10) << ", " << bReverse << ", " << nthOrder << "):" << endl;
 */   
    // the sign of the intervall controls the direction
    CDigFloat dfIntervall = CDigFloat(dfTo) - dfFrom;
    
    CDigFloat dfTargetCoverage = AbsIntegral(nthOrder)* dfCoveragePercent/100.;
    CDigFloat dfActualCoverage = AbsIntegral(bReverse ? dfTo : dfFrom ,bReverse ? dfFrom : dfTo, nthOrder);
    
    // reset errors for later break condition in while loop
    dfTargetCoverage.ResetError();
    dfActualCoverage.ResetError();
    
    // DEBUG:
//     cout << "dfTargetCoverage: " << dfTargetCoverage.RawPrint(10) << endl;
//     cout << "dfActualCoverage: " << dfActualCoverage.RawPrint(10) << endl;
//     cout << "dfIntervall: " << dfIntervall.RawPrint(10) << endl;
    
    // give up if the dfActualCoverage (this is max right now) is less than the dfTargetCoverage
    if(dfActualCoverage < dfTargetCoverage)
        return false;
    
    // do a binary search shrinking the actual to the target coverage
    // break condition: do until the intervall for the binary search is
    // less than the error of the demanded limit
    CDigFloat dfToOld = dfTo+1.;
//     while( dfActualCoverage != dfTargetCoverage  )
    while( dfToOld != dfTo )
//     while( fabs(dfIntervall.RawValue()) > dfTo.RawError()/100.  )
    {
        // remember old value
        dfToOld = dfTo;
        
        // halfen intervall
        dfIntervall /= 2.;
                
        if(dfActualCoverage < dfTargetCoverage)
        {
            dfTo += dfIntervall;
            
        }   // endif(dfActualCoverage < dfTargetCoverage)
        else
        {
            dfTo -= dfIntervall;
            
        }   // endelse(dfActualCoverage < dfTargetCoverage)
        
        // calculate the actual coverage: take raw value is the same as
        // handing over CDigFloat and then to call ResetError()
        dfActualCoverage = AbsIntegral(bReverse ? dfTo.RawValue() : dfFrom.RawValue() ,bReverse ? dfFrom.RawValue() : dfTo.RawValue(), nthOrder).RawValue();
      
    
    }   // endwhile( dfActualCoverage != dfTargetCoverage)
        //DEBUG
//         cout << "dfActualCoverage: " << dfActualCoverage.RawPrint(50) << endl;
//         cout << "dfTargetCoverage: " << dfTargetCoverage.RawPrint(50) << endl;
//         cout << "dfIntervall: " << dfIntervall.RawPrint(50) << endl;
//         cout << "dfToOld: " << dfToOld.RawPrint(50) << endl;
//         cout << "dfTo: " << dfTo.RawPrint(50) << endl;
//      
//     cout << "returning dfTo" << dfTo.RawPrint(50) << endl;
    
    
    return true;
}

bool CDistribution::CoverageIntervall(const CDigFloat& dfCoveragePercent, CDigFloat& dfMin, CDigFloat& dfMax, int nthOrder)
{
    
    // at least the mean of first order ... zeroth order will always yield 1
    CDigFloat dfMean = Median(nthOrder);
    assert(dfMean.RawValue() >= mDistribution.begin()->first.RawValue());
    assert(dfMean.RawValue() <= prev(mDistribution.end())->first.RawValue());
    
    // simply calculate the CoverageFromTo intervalls for min and max from median
    if( ! CoverageFromTo(CDigFloat( dfCoveragePercent)/2.,dfMean, dfMin,true, nthOrder) )
        return false;
    if( ! CoverageFromTo(CDigFloat( dfCoveragePercent)/2.,dfMean, dfMax,false, nthOrder) )
        return false;
    
    return true;
    
    // init the output arguments with the mean
    dfMin = dfMean;
    dfMax = dfMean;
    
    // check coverage argument 
    if( dfCoveragePercent < 0  || dfCoveragePercent > 100)
        return false;
    
    // get the target value of the integral for this order and given coverage 
    CDigFloat dfTargetCoverage = AbsIntegral(nthOrder) * dfCoveragePercent / 100.;
    
    // do a binary search for the given covery
    CDigFloat dfActualCoverage = 0;
    CDigFloat dfIntervall = max(dfMean.RawValue() - mDistribution.begin()->first.RawValue(), prev(mDistribution.end())->first.RawValue() - dfMean.RawValue());
    assert(dfIntervall.RawValue() > 0);
    dfMin = dfMean - dfIntervall;
    dfMax = dfMean + dfIntervall;

    dfActualCoverage = AbsIntegral(dfMin, dfMax,nthOrder);

    while( dfActualCoverage != dfTargetCoverage && dfIntervall > 0 )
    {
        
        dfIntervall /= 2.;
        if(dfActualCoverage > dfTargetCoverage )
        {            
            dfMin += dfIntervall;
            dfMax -= dfIntervall;
            
        }   // endif(dfActualCoverage > dfTargetCoverage )
        else
        {
            dfMin -= dfIntervall;
            dfMax += dfIntervall;
        }
        dfActualCoverage = AbsIntegral(dfMin, dfMax,nthOrder);
        
        // DEBUG
//         cout << "dfActualCoverage: " << dfActualCoverage.RawPrint(10) << endl;
//         cout << "dfIntervall: " << dfIntervall.RawPrint(10) << endl;
//         cout << "dfMin: " << dfMin.RawPrint(10) << endl;
//         cout << "dfMax: " << dfMax.RawPrint(10) << endl;
        
    }   // endwhile( dfActualCoverage != dfTargetCoverage && dfIntervall > 0)    
    
    return true;
    
}

void CDistribution::Reset()
{
    _Init();
}

std::__cxx11::string CDistribution::Print(int nPrecision)
{
    ostringstream oss;
    for(auto iel : mDistribution)
        oss << iel.first.Print(nPrecision) << " , " << iel.second.Print(nPrecision) << endl;
    
    return oss.str();
}

CDigFloat CDistribution::AbsIntegral(const int& nthOrder /*=0*/)
{
    CDigFloat dfIntegral = 0;
    for(auto iel = mDistribution.begin(); iel != mDistribution.end(); iel++)
        // now get the next integral part between two consecutive elements of this distribution
        dfIntegral += _IntegralConsecutiveElements(iel,nthOrder);
    
    return dfIntegral;
}

CDigFloat CDistribution::_IntegralConsecutiveElements(const M_DFDF::const_iterator& iel, const int& nthOrder /*=0*/)
{
    // init integral
    CDigFloat dfIntegral = 0;
    
    // set the next element
    M_DFDF::const_iterator ielNext = iel;
    ielNext++;
    
    // if we are at the end --> return zero
    if( ielNext == mDistribution.end())
        return dfIntegral;
    
    dfIntegral = _IntegralOfTwoPoints(iel->first, iel->second, ielNext->first, ielNext->second, nthOrder);
    
    return dfIntegral;
    
}

CDigFloat CDistribution::_IntegralOfTwoPoints(const CDigFloat& x1, const CDigFloat& y1, const CDigFloat& x2, const CDigFloat& y2, const int &nthOrder /*=0*/)
{
    // use in order:
    assert(x2>=x1);
    
    CDigFloat dfIntegral = 0; 
  
    // calculation of the weighted integral of a linear function y = a*x + b
    // calculate : int( x^n * (a*x + b) )dx within the limits [x1, x2]
    // solution: [ a/(n+2)*x^(n+2) + b/2*x² ] (x=x1) - [](x=x2)
    //           where a = (y2-y1)/(x2-x1)
    //           where b = y1 + a*(-x1)
    if(x1!=x2)
    {
        CDigFloat a = (CDigFloat(y1)-y2)/(CDigFloat(x1)-x2);
        CDigFloat b = CDigFloat(y1) - a * x1;
        
        dfIntegral = _nthOrderWeightedPrimitiveIntegral(x2,a,b,nthOrder)  - _nthOrderWeightedPrimitiveIntegral(x1,a,b,nthOrder);
        
//         DEBUG
//         cout << "_IntegralOfTwoPoints( " << x1.RawPrint(3) << "," << y1.RawPrint(3) <<"," << x2.RawPrint(3) <<"," << y2.RawPrint(3) << "," << to_string(nthOrder) << "):" << endl;
//         cout << "a=" << a.RawPrint(10) << endl;
//         cout << "b=" << b.RawPrint(10) << endl;
//         cout << "upper limit integral: " << _nthOrderWeightedPrimitiveIntegral(x2,a,b,nthOrder).RawPrint(10) << endl;
//         cout << "lower limit integral: " << _nthOrderWeightedPrimitiveIntegral(x1,a,b,nthOrder).RawPrint(10) << endl;
    }
    
    return dfIntegral;    
    
}

CDigFloat CDistribution::_nthOrderWeightedPrimitiveIntegral(const CDigFloat& x, const CDigFloat& slope, const CDigFloat& offset, const int& nthOrder)
{
        // calculation of the weighted integral of a linear function y = a*x + b
        // calculates : int( x^n * (a*x + b) )dx within the limits [x1, x2]
        // solution: [ a/(n+2)*x^(n+2) + b/2*x² ] 
    return CDigFloat(slope)/(nthOrder+2.) * pow(x,nthOrder+2) + CDigFloat(offset)/(nthOrder+1.) * pow(x,nthOrder+1);
}


void CDistribution::_Init()
{
    mDistribution.clear();
}


