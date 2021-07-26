
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
    LOGTRACE(LS_Dist+"Shift", "called with arg:");
    LOGTRACE(LS_Dist+"Shift", string("CDigFloat shift = ") + shift.RawPrint(10));
    CDigFloat dfShiftReset(shift.RawValue());
    LOGTRACE(LS_Dist+"Shift", string("CDigFloat shift reset = ") + dfShiftReset.RawPrint(10));
    LOGTRACE(LS_Dist+"Shift", string("distri before \n") + Print(10));
    
     // generate map to copy from
    MapDFDFType ShiftedDistribution;
    ShiftedDistribution.clear();
    
    for(auto iel: mDistribution)
    {
        CDigFloat dfIndex = iel.first + dfShiftReset;
        ShiftedDistribution[dfIndex] = iel.second;
        LOGTRACE(LS_Dist+"Shift", string("dfIndex = ") + dfIndex.RawPrint(10));
        LOGTRACE(LS_Dist+"Shift", string("ShiftedDistribution size = ") + to_string(ShiftedDistribution.size() ));
    
    }
    
    // copy
    mDistribution.clear();
 
 
    for(auto iel: ShiftedDistribution)
    {
        mDistribution[iel.first] = iel.second;
        
        LOGTRACE(LS_Dist+"Shift", string("copy ShiftedDistribution back [") + iel.first.RawPrint(10) + "] = " + mDistribution[iel.first].RawPrint(10) );
    }
 
    
    LOGTRACE(LS_Dist+"Shift", string("distri after \n") + Print(10));
    
    // clear temporary map
    ShiftedDistribution.clear();    
}

void CDistribution::Scale(const CDigFloat& scale)
{
    // generate map to copy from
    MapDFDFType ScaledDistribution;
    for(auto iel: mDistribution)
        ScaledDistribution[CDigFloat(scale)*iel.first] = iel.second;
    
    // copy
    mDistribution.clear();
    for(auto iel: ScaledDistribution)
        mDistribution[iel.first] = iel.second;
    
    // clear temporary map
    ScaledDistribution.clear();   
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

void CDistribution::GetInterval(const CDigFloat& variable, MapDFDFType::const_iterator& VariableLeft, MapDFDFType::const_iterator& VariableRight)
{
    LOGTRACE(LS_Dist + "GetInterval", string("called with argument :"));
    LOGTRACE(LS_Dist + "GetInterval", string("CDigFloat:") + variable.RawPrint(20));
    LOGTRACE(LS_Dist + "GetInterval", string("this Distri:") + PrintMetaInfo());
    
    int nInterval = mDistribution.size()-1;    
    VariableLeft = mDistribution.begin();
    
    
    // set imax to the end
    VariableRight = VariableLeft;
    
    // if the map is of size zero: leave now
    if(nInterval < 0)
        return;
    
    advance( VariableRight ,nInterval);
 
    // if variable is left outside ...
    if(variable < VariableLeft->first)
    {
        // trace: can be controlled by project::class::function
        LOGTRACE(LS_Dist + "GetInterval", string("variable less than left limit: setting both limits to left limit of distri"));
    
        
        // ... set right variable and left variable to begin() and leave
        VariableRight = VariableLeft;
        
        
        // trace: can be controlled by project::class::function
        LOGTRACE(LS_Dist + "GetInterval", string("VariableLeft = ") + ::Print(VariableLeft));
        LOGTRACE(LS_Dist + "GetInterval", string("VariableLeft = ") + ::Print(VariableRight));
        return;
    }        
        
    // if variable is right outside ...
    if(variable > VariableRight->first)
    {
        // trace: can be controlled by project::class::function
        LOGTRACE(LS_Dist + "GetInterval", string("variable greater than right limit: setting both limits to right limit of distri"));
        
        
        // ... set right variable and left variable to end() and leave
        VariableRight = mDistribution.end();
        VariableLeft = VariableRight;
        
        // trace: can be controlled by project::class::function
        LOGTRACE(LS_Dist + "GetInterval", string("VariableLeft = ") + ::Print(VariableLeft));
        LOGTRACE(LS_Dist + "GetInterval", string("VariableRight = ") + ::Print(VariableRight));
        
        // leave
        return;
    }
    
    // now start binary search:
    // trace: can be controlled by project::class::function
    LOGTRACE(LS_Dist + "GetInterval", string("starting binary search ..."));
    MapDFDFType::const_iterator iact = VariableLeft;
    while( nInterval > 1 )
    {
        
        // trace: can be controlled by project::class::function
        LOGTRACE(LS_Dist + "GetInterval", string("... interval(begin)      = ") + to_string(nInterval));
        LOGTRACE(LS_Dist + "GetInterval", string("... VariableLeft(begin)  = ") + ::Print(VariableLeft));
        LOGTRACE(LS_Dist + "GetInterval", string("... VariableRight(begin) = ") + ::Print(VariableRight));
        LOGTRACE(LS_Dist + "GetInterval", string("... iact(begin)          = ") + ::Print(iact));
        
        // halfen the interval: do this conservatively (nInterval+1)
        // in case you only have 3 intervals (i.e. 4 points) you'll leave
        // out one intervall
        nInterval = (nInterval+1)/2;
        
        // set actual value to new middle
        iact = VariableLeft;
        
        // secure way to increment the actual iterator
        for(int istep = 0; istep<nInterval; istep++)
        {
            iact++;
            if(iact == Distribution().end())
            {
                LOGWARN("WarningLogger", LS_Dist + "GetInterval:\niterator exceeds distribution!!");
                
                // leave loop if we were at the end already 
                if(istep == 0)
                    break;
                
                // halfen the interval
                nInterval = (istep+1) / 2.;
                
                //decrement 
                for(int backStep = 0; backStep < nInterval; backStep++)
                    iact--;
                
                
                LOGTRACE(LS_Dist + "GetInterval", string("... end of distri reached: backStep = ") + to_string(nInterval));
                
            }   // endif(iact == Distribution().end())
        }
        
        // set new interval limits
        if(variable <= iact->first)
            VariableRight = iact;
        else        
            VariableLeft = iact;
        
        
        // trace: can be controlled by project::class::function
        LOGTRACE(LS_Dist + "GetInterval", string("... interval(after)      = ") + to_string(nInterval));
        LOGTRACE(LS_Dist + "GetInterval", string("... VariableLeft(after)  = ") + ::Print(VariableLeft));
        LOGTRACE(LS_Dist + "GetInterval", string("... VariableRight(after) = ") + ::Print(VariableRight));
        LOGTRACE(LS_Dist + "GetInterval", string("... iact(after)          = ") + ::Print(iact));
            
    }   
    
    // do a final iteration: check if left could be incremented by one and is still larger equal than the variable
    iact = VariableLeft; iact++;
    
    // only do if not already the tightest interval
    if(variable >= iact->first && VariableRight->first > iact->first)
        VariableLeft = iact;
 
    // do a final iteration: check if right could be decremented by one and is still less equal than the variable
    iact = VariableRight; iact--;
    if(variable <= iact->first && VariableLeft->first < iact->first)
        VariableRight = iact;
    
    // trace: can be controlled by project::class::function
    LOGTRACE(LS_Dist + "GetInterval", string("... VariableLeft(final)  = ") + ::Print(VariableLeft));
    LOGTRACE(LS_Dist + "GetInterval", string("... VariableRight(final) = ") + ::Print(VariableRight));
}


CDigFloat CDistribution::DistriValue(const CDigFloat& variable)
{   
    
    LOGTRACE(LS_Dist + "DistriValue", string("called with argument :"));
    LOGTRACE(LS_Dist + "DistriValue", string("CDigFloat =") + variable.RawPrint(20));

    // init result with zero
    CDigFloat dfDistriValue(0);
    
    MapDFDFType::const_iterator imin, imax;
    GetInterval(variable, imin, imax);
    
    LOGTRACE(LS_Dist + "DistriValue", string("imin =") + ::Print(imin));
    LOGTRACE(LS_Dist + "DistriValue", string("imax =") + ::Print(imax));

    
    // in case variable is left out of bounds: GetInterval did not find a matching interval --> imin = mDistribution.begin()
    if( (variable < imin->first ) || 
        (imin ==imax && imax == mDistribution.end())
      )
    {
        LOGWARN("WarningLogger", LS_Dist + "DistriValue: value is out of bounds");
    }    
    else        
    {
        if(imin->first == imax->first)
            dfDistriValue = imin->second;
        else
            // now we got the min max --> get the linearily interpolated value
            dfDistriValue = CDigFloat(imin->second) + (CDigFloat(imax->second) - imin->second) / (CDigFloat(imax->first) - imin->first) * (CDigFloat(variable) - imin->first);
    }
    
    LOGTRACE(LS_Dist + "DistriValue", string("interpolated value =") + dfDistriValue.RawPrint(20));
 
    return dfDistriValue;
    
}

CDigFloat CDistribution::LinearOffset(const CDigFloat& Variable)
{
    // trace: can be controlled by project::class::function
    LOGTRACE(LS_Dist + "LinearOffset", string("called with argument :"));
    LOGTRACE(LS_Dist + "LinearOffset", string("CDigFloat:") + Variable.RawPrint(20)); 
    
    MapDFDFType::const_iterator Left, Right;
    GetInterval(Variable, Left, Right);
    
    // trace: can be controlled by project::class::function
    LOGTRACE(LS_Dist + "LinearOffset", string("got the interval :"));
    LOGTRACE(LS_Dist + "LinearOffset", string("left iterator :") + ::Print(Left));
    LOGTRACE(LS_Dist + "LinearOffset", string("right iterator :") + ::Print(Right));
    
    return _LinearOffset(Left);
}

CDigFloat CDistribution::LinearSlope(const CDigFloat& Variable)
{
    // trace: can be controlled by project::class::function
    LOGTRACE(LS_Dist + "LinearSlope", string("called with argument :"));
    LOGTRACE(LS_Dist + "LinearSlope", string("CDigFloat:") + Variable.RawPrint(20)); 
    
    MapDFDFType::const_iterator Left, Right;
    GetInterval(Variable, Left, Right);

    return _LinearSlope(Left);
}


CDigFloat CDistribution::AbsIntegral(const CDigFloat& variableLeft, const CDigFloat& variableRight, int nthOrder /*= 0*/)
{
    LOGTRACE(LS_Dist+"AbsIntegral", "called with arguments:");
    LOGTRACE(LS_Dist+"AbsIntegral", string("CDigFloat Left = ")+ variableLeft.RawPrint(20));
    LOGTRACE(LS_Dist+"AbsIntegral", string("CDigFloat Right= ")+ variableRight.RawPrint(20));
    LOGTRACE(LS_Dist+"AbsIntegral", string("int nthOrder   = ")+ to_string(nthOrder));
    LOGTRACE(LS_Dist+"AbsIntegral", string("this distri: ") + PrintMetaInfo());
    
    // must be in order: left <= right
    assert(variableLeft <= variableRight);    
    
    // init variable for result
    CDigFloat dfIntegral(0);    
    
    // check for being within limits otherwise return 0
    // in case of equality --> return zero, too 
    if( (variableLeft > prev(Distribution().end())->first || variableRight < Distribution().begin()->first) ||
        ( variableLeft == variableRight)
    )
    {
        LOGTRACE(LS_Dist+"AbsIntegral", string("returning zero because \"out of limits\" or \"equality of limits\" "));
        return dfIntegral;
    }
    
    // correct exceeding limits --> set to edges
    CDigFloat variableLeftCorr = variableLeft;
    CDigFloat variableRightCorr = variableRight;
    if( variableLeft <= Distribution().begin()->first)
        variableLeftCorr = Distribution().begin()->first;
    if( variableRight >= prev(Distribution().end())->first)
        variableRightCorr = prev(Distribution().end())->first;
    
    LOGTRACE(LS_Dist+"AbsIntegral", string("CDigFloat Left corrected  = ")+ variableLeftCorr.RawPrint(20));
    LOGTRACE(LS_Dist+"AbsIntegral", string("CDigFloat Right corrected = ")+ variableRightCorr.RawPrint(20));
    
    bool bSameInterval = false;
    MapDFDFType::const_iterator iminLeft, imaxLeft, iminRight, imaxRight;
    GetInterval(variableLeftCorr, iminLeft, imaxLeft);
    
    // set limits for right value 
    if( variableRightCorr < iminLeft->first || variableRightCorr > imaxLeft->first)
        GetInterval(variableRightCorr, iminRight, imaxRight);
    else
    {
        iminRight = iminLeft;
        imaxRight = imaxLeft;
        bSameInterval = true;
    }
    
    
    LOGTRACE(LS_Dist+"AbsIntegral", string("interval left: [ ") + iminLeft->first.RawPrint(20) + ", " + imaxLeft->first.RawPrint(20) );
    LOGTRACE(LS_Dist+"AbsIntegral", string("interval right: [ ") + iminRight->first.RawPrint(20) + ", " + imaxRight->first.RawPrint(20) );
    LOGTRACE(LS_Dist+"AbsIntegral", string("bSameInterval: ") + Bool2String(bSameInterval));
    
    // check for same interval to avoid double calculations
    // in this case no calculation is done for 
    // complete elements
    // right part
    
    // set the limit for the left integral:
    CDigFloat dfLeftIntegralMax = min(imaxLeft->first, variableRightCorr);
    
    // calculate the integral of complete elements: only if the
    // variable arguments of this function are NOT in the same interval
    if( !bSameInterval )
        for(auto iel = imaxLeft; iel!=iminRight; iel++)
            dfIntegral += _IntegralConsecutiveElements(iel,nthOrder);
        
        
    LOGTRACE(LS_Dist+"AbsIntegral", string("Integral (consecutive elements) = ")+ dfIntegral.RawPrint(20));
    
    // add the part on the left:
    dfIntegral += _IntegralOfTwoPoints( variableLeftCorr, DistriValue(variableLeftCorr), dfLeftIntegralMax, DistriValue(dfLeftIntegralMax),nthOrder);
    
    LOGTRACE(LS_Dist+"AbsIntegral", string("Integral (left part added) = ")+ dfIntegral.RawPrint(20));

    // add the part on the right:
    if(iminRight != mDistribution.end() && !bSameInterval)
        dfIntegral += _IntegralOfTwoPoints(iminRight->first, iminRight->second, variableRightCorr, DistriValue(variableRightCorr),nthOrder);

    LOGTRACE(LS_Dist+"AbsIntegral", string("Integral (right part added: final) = ")+ dfIntegral.RawPrint(20));    
    
    return CDigFloat( dfIntegral );
}
CDigFloat CDistribution::Min()
{
    CDigFloat dfMin = front()->second;
    for(auto iel: Distribution())
        if(iel.second < dfMin)
            dfMin = iel.second;
    return dfMin;
}

CDigFloat CDistribution::Max()
{
    CDigFloat dfMax = front()->second;
    for(auto iel: Distribution())
        if(iel.second > dfMax)
            dfMax = iel.second;
    return dfMax;
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
    LOGTRACE(LS_Dist+"CoverageFromTo", "called with args:");
    LOGTRACE(LS_Dist+"CoverageFromTo", string("CDigFloat coverage[%] = ")+ dfCoveragePercent.RawPrint(20));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("CDigFloat From        = ")+ dfFrom.RawPrint(20));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("CDigFloat To[out]     = ")+ Address2String(&dfTo));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("bool reverse          = ")+ Bool2String(bReverse));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("int order             = ")+ to_string(nthOrder));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("this distri: ")+ PrintMetaInfo());

    // init the output variable with the corresponding starting point depending on the direction
    CDigFloat dfStart = dfFrom;
    dfTo = bReverse ? mDistribution.begin()->first : prev(mDistribution.end())->first;
 
    // the sign of the interval controls the direction
    CDigFloat dfInterval = CDigFloat(dfTo) - dfFrom;
 
    
    CDigFloat dfTargetCoverage = AbsIntegral(nthOrder)* dfCoveragePercent/100.;
    CDigFloat dfActualCoverage = AbsIntegral(bReverse ? dfTo : dfFrom ,bReverse ? dfFrom : dfTo, nthOrder);
    
    
    // reset errors for later break condition in while loop
    dfTargetCoverage.ResetError();
    dfActualCoverage.ResetError();
    
    // give up if the dfActualCoverage (this is max right now) is less than the dfTargetCoverage
    if(dfActualCoverage < dfTargetCoverage)
    {
        LOGWARN("WarningLogger", LS_Dist+"CoverageFromTo: coverage cannot be achieved with given parameters:");
        LOGWARN("WarningLogger", LS_Dist+"CoverageFromTo: calculated target coverage = " + dfTargetCoverage.RawPrint(20));
        LOGWARN("WarningLogger", LS_Dist+"CoverageFromTo: calculated actual coverage = " + dfActualCoverage.RawPrint(20));
        LOGWARN("WarningLogger", LS_Dist+"CoverageFromTo called with args:");
        LOGWARN("WarningLogger", LS_Dist+"CoverageFromTo: CDigFloat coverage[%] = " + dfCoveragePercent.RawPrint(20));
        LOGWARN("WarningLogger", LS_Dist+"CoverageFromTo: CDigFloat From        = " + dfFrom.RawPrint(20));
        LOGWARN("WarningLogger", LS_Dist+"CoverageFromTo: CDigFloat To[out]     = " + Address2String(&dfTo));
        LOGWARN("WarningLogger", LS_Dist+"CoverageFromTo: bool reverse          = " + Bool2String(bReverse));
        LOGWARN("WarningLogger", LS_Dist+"CoverageFromTo: int order             = " + to_string(nthOrder));
        LOGWARN("WarningLogger", LS_Dist+"CoverageFromTo: this distri: " + PrintMetaInfo());
        return false;
    }
    // do a binary search shrinking the actual to the target coverage
    // break condition: do until the interval for the binary search is
    // less than the error of the demanded limit
    
    CDigFloat dfToOld = dfTo+1.;
    int count = 0;
    
    LOGTRACE(LS_Dist+"CoverageFromTo", "preparing binary search with args:");
    LOGTRACE(LS_Dist+"CoverageFromTo", string("target coverage   = ")+ dfTargetCoverage.RawPrint(20));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("actual coverage   = ")+ dfActualCoverage.RawPrint(20));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("To                = ")+ dfTo.RawPrint(20));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("ToOld             = ")+ dfToOld.RawPrint(20));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("Interval(reverse) = ")+ dfInterval.RawPrint(20));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("count             = ")+ to_string(count));
    LOGTRACE(LS_Dist+"CoverageFromTo", string("max. iterations   = ")+ to_string(MaxBinarySearchIterations()));
    LOGTRACE(LS_Dist+"CoverageFromTo", "starting binary search ... ");
    while( dfActualCoverage != dfTargetCoverage && count < MaxBinarySearchIterations() )
    {
        count++;
        // remember old value
        dfToOld = dfTo.RawValue();
        
        // halfen interval
        dfInterval /= 2.;
                
        if(dfActualCoverage < dfTargetCoverage)
        {
            dfTo += dfInterval;
            
        }   // endif(dfActualCoverage < dfTargetCoverage)
        else
        {
            dfTo -= dfInterval;
            
        }   // endelse(dfActualCoverage < dfTargetCoverage)
        
        // the error might accumulate and the break condition might be true earlier due to the accumulated errors
        // dfTo.ResetError();
        
        // calculate the actual coverage: take raw value is the same as
        // handing over CDigFloat and then to call ResetError()
        dfActualCoverage = AbsIntegral(bReverse ? dfTo.RawValue() : dfFrom.RawValue() ,bReverse ? dfFrom.RawValue() : dfTo.RawValue(), nthOrder).RawValue();
        
        LOGTRACE(LS_Dist+"CoverageFromTo", string("... target coverage   = ")+ dfTargetCoverage.RawPrint(20));
        LOGTRACE(LS_Dist+"CoverageFromTo", string("... actual coverage   = ")+ dfActualCoverage.RawPrint(20));
        LOGTRACE(LS_Dist+"CoverageFromTo", string("... To                = ")+ dfTo.RawPrint(20));
        LOGTRACE(LS_Dist+"CoverageFromTo", string("... ToOld             = ")+ dfToOld.RawPrint(20));
        LOGTRACE(LS_Dist+"CoverageFromTo", string("... Interval(reverse) = ")+ dfInterval.RawPrint(20));
        LOGTRACE(LS_Dist+"CoverageFromTo", string("... count             = ")+ to_string(count));
      
    
    }   // endwhile( dfActualCoverage != dfTargetCoverage)
    
    return true;
}

bool CDistribution::CoverageInterval(const CDigFloat& dfCoveragePercent, CDigFloat& dfMin, CDigFloat& dfMax, int nthOrder)
{
    LOGTRACE(LS_Dist+"CoverageInterval", "called with args:");
    LOGTRACE(LS_Dist+"CoverageInterval", string("CDigFloat coverage[%] = ")+ dfCoveragePercent.RawPrint(20));
    LOGTRACE(LS_Dist+"CoverageInterval", string("CDigFloat Min[out]    = ")+ Address2String(&dfMin));
    LOGTRACE(LS_Dist+"CoverageInterval", string("CDigFloat Max[out]    = ")+ Address2String(&dfMax));
    LOGTRACE(LS_Dist+"CoverageInterval", string("int order             = ")+ to_string(nthOrder));
    
    // at least the mean of first order ... zeroth order will always yield 1
    CDigFloat dfMean = Median(nthOrder);
    assert(dfMean.RawValue() >= mDistribution.begin()->first.RawValue());
    assert(dfMean.RawValue() <= prev(mDistribution.end())->first.RawValue());
    
    
    LOGTRACE(LS_Dist+"CoverageInterval", string("caclulated median = ")+  dfMean.RawPrint(20));
    
    // simply calculate the CoverageFromTo intervals for min and max from median
    if( ! CoverageFromTo(CDigFloat( dfCoveragePercent)/2.,dfMean, dfMin,true, nthOrder) )
        return false;
    if( ! CoverageFromTo(CDigFloat( dfCoveragePercent)/2.,dfMean, dfMax,false, nthOrder) )
        return false;
    
    
    LOGTRACE(LS_Dist+"CoverageInterval", string("success: "));
    LOGTRACE(LS_Dist+"CoverageInterval", string("Min = ") + dfMin.RawPrint(20));
    LOGTRACE(LS_Dist+"CoverageInterval", string("Max = ") + dfMax.RawPrint(20));
    
    return true;
    
}

void CDistribution::Reset()
{
    _Init();
}
std::__cxx11::string CDistribution::PrintMetaInfo()
{
    return string("Distribution: nPoints = ") + to_string(Distribution().size()) + ", x-interval:[ " + firstVariable().RawPrint(10) + ", " + lastVariable().RawPrint(10) + " ]";
    
}

std::__cxx11::string CDistribution::Print(int nPrecision, bool bWithError)
{
    ostringstream oss;
    for(auto iel : mDistribution)
        oss << iel.first.RawPrint(nPrecision,bWithError) << " , " << iel.second.RawPrint(nPrecision,bWithError) << endl;
    
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

CDigFloat CDistribution::_IntegralConsecutiveElements(const MapDFDFType::const_iterator& iel, const int& nthOrder /*=0*/)
{
    
    LOGTRACE(LS_Dist+"_IntegralConsecutiveElements", "called with args:");
    LOGTRACE(LS_Dist+"_IntegralConsecutiveElements", string("MapDFDFType::const_iterator:") + ::Print(iel));
    LOGTRACE(LS_Dist+"_IntegralConsecutiveElements", string("int order:") + to_string(nthOrder));
    
    // init integral
    CDigFloat dfIntegral = 0;
    
    // set the next element
    MapDFDFType::const_iterator ielNext = iel;
    ielNext++;
    
    // if we are at the end --> return zero
    if( ielNext == mDistribution.end())
        return dfIntegral;
    
    dfIntegral = _IntegralOfTwoPoints(iel->first, iel->second, ielNext->first, ielNext->second, nthOrder);

    return dfIntegral;
    
}

CDigFloat CDistribution::_IntegralOfTwoPoints(const CDigFloat& x1, const CDigFloat& y1, const CDigFloat& x2, const CDigFloat& y2, const int &nthOrder /*=0*/)
{
    
    LOGTRACE(LS_Dist + "_IntegralOfTwoPoints", "called with args:");
    LOGTRACE(LS_Dist + "_IntegralOfTwoPoints", string("CDigFloat x1 =") + x1.RawPrint(20));
    LOGTRACE(LS_Dist + "_IntegralOfTwoPoints", string("CDigFloat y1 =") + y1.RawPrint(20));
    LOGTRACE(LS_Dist + "_IntegralOfTwoPoints", string("CDigFloat x2 =") + x2.RawPrint(20));
    LOGTRACE(LS_Dist + "_IntegralOfTwoPoints", string("CDigFloat y2 =") + y2.RawPrint(20));
    LOGTRACE(LS_Dist + "_IntegralOfTwoPoints", string("int order    =") + to_string(nthOrder));
    
    CDigFloat dfIntegral = 0; 
  
    // calculation of the weighted integral of a linear function y = a*x + b
    // calculate : int( x^n * (a*x + b) )dx within the limits [x1, x2]
    // solution: [ a/(n+2)*x^(n+2) + b/2*x² ] (x=x1) - [](x=x2)
    //           where a = (y2-y1)/(x2-x1)
    //           where b = (y2*x1 - y1*x2)/(x1-x2)
    if(x1!=x2)
    {
        // use in order:
        assert(x2>=x1);
        CDigFloat a = (y1-y2)/(x1-x2);
        CDigFloat b = (y2*x1 - y1*x2)/(x1-x2);
        
        LOGTRACE(LS_Dist + "_IntegralOfTwoPoints", string("calculated slope =") + a.RawPrint(20));
        LOGTRACE(LS_Dist + "_IntegralOfTwoPoints", string("calculated offset=") + b.RawPrint(20));
        
        dfIntegral = _nthOrderWeightedIntegral(x1,x2,a,b,nthOrder);
        
        LOGTRACE(LS_Dist + "_IntegralOfTwoPoints", string("calculated specific integral=") + dfIntegral.RawPrint(20));
        
    }
    
    return dfIntegral;    
    
}
CDigFloat CDistribution::_nthOrderWeightedIntegral(const CDigFloat& dfX1, const CDigFloat& dfX2, const CDigFloat& slope, const CDigFloat& offset, const int& nthOrder)
{
    LOGTRACE(LS_Dist + "_nthOrderWeightedPrimitiveIntegral", "called with args:");
    LOGTRACE(LS_Dist + "_nthOrderWeightedPrimitiveIntegral", string("CDigFloat dfX1 =") + dfX1.RawPrint(20));
    LOGTRACE(LS_Dist + "_nthOrderWeightedPrimitiveIntegral", string("CDigFloat dfX2 =") + dfX2.RawPrint(20));
    LOGTRACE(LS_Dist + "_nthOrderWeightedPrimitiveIntegral", string("CDigFloat slope  =") + slope.RawPrint(20));
    LOGTRACE(LS_Dist + "_nthOrderWeightedPrimitiveIntegral", string("CDigFloat offset =") + offset.RawPrint(20));
    LOGTRACE(LS_Dist + "_nthOrderWeightedPrimitiveIntegral", string("int order        =") + to_string(nthOrder));
    LOGTRACE(LS_Dist + "_nthOrderWeightedPrimitiveIntegral", string("result =") + (slope/(nthOrder+2.) * DifferenceNthOrder(dfX1, dfX2,nthOrder+2) + offset/(nthOrder+1.) * DifferenceNthOrder(dfX1,dfX2 ,nthOrder+1)).RawPrint(20) );
    
        // calculation of the weighted integral of a linear function y = a*x + b
        // calculates : int( x^n * (a*x + b) )dx within the limits [x1, x2]
        // solution: [ a/(n+2)*x^(n+2) + b/(n+1)*x^(n+1) ] 
    return slope/(nthOrder+2.) * DifferenceNthOrder(dfX1, dfX2,nthOrder+2) + offset/(nthOrder+1.) * DifferenceNthOrder(dfX1,dfX2 ,nthOrder+1);
}

CDigFloat CDistribution::_LinearOffset(MapDFDFType::const_iterator& Left)
{   
     // trace: can be controlled by project::class::function
    LOGTRACE(LS_Dist + "_LinearOffset", string("called with argument :"));
    LOGTRACE(LS_Dist + "_LinearOffset", string("left iterator  (x1,y1) :") + ::Print(Left));
    
    // check if is outside the distri 
    if(Left == Distribution().end() )
    {
            LOGTRACE(LS_Dist + "_LinearOffset", string("iterator exceeds distri limits: returning 0"));
            return 0;

    }   // endif(Left == Distribution().end())
    
    MapDFDFType::const_iterator Right = Left;
    Right++;
 
    // check if is outside the distri 
    if(Right == Distribution().end() )
    {
            LOGTRACE(LS_Dist + "_LinearOffset", string("iterator exceeds distri limits: returning 0"));
            return 0;

    }   // endif(Left == Distribution().end())
  
    // trace: can be controlled by project::class::function
    LOGTRACE(LS_Dist + "_LinearOffset", "calculation offset: (x2*y1 - x1*y2) / (x2 - x1)");
    LOGTRACE(LS_Dist + "_LinearOffset", string("(x2*y1 - x1*y2) = ") + (Right->first*Left->second - Left->first*Right->second).RawPrint(20)); 
    LOGTRACE(LS_Dist + "_LinearOffset", string("(x2 - x1) = ") + (Right->first - Left->first).RawPrint(20)); 
    LOGTRACE(LS_Dist + "_LinearOffset", string(" result = ") + ((Right->first*Left->second - Left->first*Right->second) / (Right->first - Left->first)        
    ).RawPrint(20)); 
    
    
    return (Right->first*Left->second - Left->first*Right->second) / (Right->first - Left->first);
    
}


CDigFloat CDistribution::_LinearSlope(MapDFDFType::const_iterator& Left)
{
    // trace: can be controlled by project::class::function
    LOGTRACE(LS_Dist + "_LinearSlope", string("called with argument :"));
    LOGTRACE(LS_Dist + "_LinearSlope", string("left iterator  (x1,y1) :") + ::Print(Left));
    
    // check if is outside the distri 
    if(Left == Distribution().end() )
    {
            LOGTRACE(LS_Dist + "_LinearOffset", string("iterator exceeds distri limits: returning 0"));
            return 0;

    }   // endif(Left == Distribution().end())
    
    MapDFDFType::const_iterator Right = Left;
    Right++;
 
    // check if is outside the distri 
    if(Right == Distribution().end() )
    {
            LOGTRACE(LS_Dist + "_LinearOffset", string("iterator exceeds distri limits: returning 0"));
            return 0;

    }   // endif(Left == Distribution().end())
    
    // trace: can be controlled by project::class::function
    LOGTRACE(LS_Dist + "_LinearSlope", "calculation offset: (y2-y1) / (x2 - x1)");
    LOGTRACE(LS_Dist + "_LinearSlope", string("(y2 - y1) = ") + (Right->second - Left->second).RawPrint(20)); 
    LOGTRACE(LS_Dist + "_LinearSlope", string("(x2 - x1) = ") + (Right->first - Left->first).RawPrint(20)); 
    LOGTRACE(LS_Dist + "_LinearSlope", string(" result = ") + ((Right->second - Left->second) / (Right->first - Left->first)).RawPrint(20)); 
    
    return (Right->second - Left->second) / (Right->first - Left->first);

}

void CDistribution::_Init()
{
    mDistribution.clear();
    MaxBinarySearchIterations(MAX_BINARY_SEARCH_ITERATIONS);
}


/////////////////////////////////////////
// external functions
//////////////////////////////////////
string Print(const MapDFDFType::const_iterator& it)
{
    return string("( ") + it->first.RawPrint(20) + ", " + it->second.RawPrint(20) + " )";
}
