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
#include "iostream"
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
// arithmetic operator
////////////////////////////////////////////
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator+=(CProbabilityDensityDistribution& other)
{
    
//     return _GeneralOperatorsFunction(other, ProbDistOp::pdoPlus );
    return _Convolution4GeneralOperators(other, ProbDistOp::pdoPlus );
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator+(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult += other;
    return ProbDistResult;
}
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator*=(CProbabilityDensityDistribution& other)
{
//     return _GeneralOperatorsFunction(other, ProbDistOp::pdoMult );
    return _Convolution4GeneralOperators(other, ProbDistOp::pdoMult );
}
CProbabilityDensityDistribution CProbabilityDensityDistribution::operator*(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult *= other;
    return ProbDistResult;
}

CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator/=(CProbabilityDensityDistribution& other)
{
    return _Convolution4GeneralOperators(other, ProbDistOp::pdoDiv);
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator/(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult /= other;
    return ProbDistResult;
}



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
    _Normalize();
    return CDistribution::DistriValue(variable);
}
CDigFloat CProbabilityDensityDistribution::AbsIntegral(const CDigFloat& variableLeft, const CDigFloat& variableRight, int nthOrder)
{
    _Normalize();
    return CDistribution::AbsIntegral(variableLeft,variableRight,nthOrder);
}
CDigFloat CProbabilityDensityDistribution::AbsIntegral(const int& nthOrder)
{
    _Normalize();
    return CDistribution::AbsIntegral(nthOrder);
}
void CProbabilityDensityDistribution::GetIntervall(const CDigFloat& variable, M_DFDF::const_iterator& VariableLeft, M_DFDF::const_iterator& VariableRight)
{
    _Normalize();
    return CDistribution::GetIntervall(variable,VariableLeft,VariableRight);
}
CDigFloat CProbabilityDensityDistribution::Mean(int nthOrder)
{
    _Normalize();
    return CDistribution::Mean(nthOrder);
}
CDigFloat CProbabilityDensityDistribution::OrigIntegral()
{
    _Normalize();
    return dfAbsIntegral;
}

void CProbabilityDensityDistribution::Reset()
{
    _Init();
}
bool CProbabilityDensityDistribution::Add(const CDigFloat& variable, const CDigFloat value, bool bLastElement /*== false*/)
{
    // first check if the distribution is already normalized
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

void CProbabilityDensityDistribution::_GetIntegralLimits4DistriIterator(const M_DFDF::const_iterator& iel, const M_DFDF& DistriMap, CDigFloat& IntegralLowerLimit, CDigFloat& IntegralUpperLimit)
{    
    // calculate the lower limit: 
    IntegralLowerLimit = iel->first;    
    // if iel is not start:  --> half way back to previous distri point
    if(iel != DistriMap.begin()) 
        IntegralLowerLimit -= (CDigFloat(iel->first) - prev(iel)->first)/2.;
    
    // calculate the upper limit: 
    IntegralUpperLimit = iel->first;
    // if iel is not the last element:  --> half way forward to next distri point
    if(next(iel) != DistriMap.end())
        IntegralUpperLimit+=(CDigFloat(next(iel)->first) - iel->first)/2.;
    
    
}

CDigFloat CProbabilityDensityDistribution::_GetMinVariableDifference(const M_DFDF& Distri)
{
       
        // get the min distance for convolution stepping
    CDigFloat MinDifference;
    for(auto iel = Distri.begin(); iel != Distri.end(); iel++)
        if( next(iel) != Distri.end())
        {
            CDigFloat ActualStep = CDigFloat(next(iel)->first) - iel->first;
            if(iel == Distri.begin())
                MinDifference = ActualStep;
            else
                if(MinDifference > ActualStep)
                    MinDifference = ActualStep;
        }
        
    return MinDifference;
}

CProbabilityDensityDistribution& CProbabilityDensityDistribution::_GeneralOperatorsFunction(CProbabilityDensityDistribution& other, const ProbDistOp Operation)
{
       // not nice but necessary: normalize the other, too: this is the state all the distribution should be in
    _Normalize();
    other._Normalize();
    
    // get the min distance for convolution stepping
    CDigFloat MinVarStepThis = _GetMinVariableDifference(Distribution());
    CDigFloat MinVarStepOther = _GetMinVariableDifference(other.Distribution());
    
//     cout << endl << "_GeneralOperatorsFunction" << endl;
//     cout << "MinVarStepThis " << MinVarStepThis.RawPrint(10) << endl;
//     cout << "MinVarStepOther " << MinVarStepOther.RawPrint(10) << endl;
    
    // calculate the number of integration steps within the distribution
    CDigFloat nIntStepsThis = VariableRange() / MinVarStepThis;
    CDigFloat nIntStepsOther = other.VariableRange() / MinVarStepOther ;
//     cout << "nIntStepsThis " << nIntStepsThis.RawPrint(10) << endl;
//     cout << "nIntStepsOther " << nIntStepsOther.RawPrint(10) << endl;
    
    // declare the half integral intervalls:
    CDigFloat IntegrationIntervallThis = VariableRange()/nIntStepsThis;
    CDigFloat IntegrationIntervallOther = other.VariableRange()/nIntStepsOther;
    
    // declare the result distri
    M_DFDF ResultDistri;
    
    // declare the variables for both distris to step through
    CDigFloat VarThis = CDigFloat( Distribution().begin()->first );
    CDigFloat VarOther;
    
    // the final variable values: end outside the distribution 
    CDigFloat VarEndThis = CDigFloat( prev(Distribution().end())->first );
    CDigFloat VarEndOther = CDigFloat( prev(other.Distribution().end())->first );
    
    // while looping until the end
    while( VarThis <= VarEndThis   )
    {        
        // reset VarThis
        VarOther = CDigFloat( other.Distribution().begin()->first ) ;
        
        while( VarOther <= VarEndOther )
        {
            CDigFloat variable;
            switch(Operation)
            {
                case ProbDistOp::pdoPlus:
                    variable = VarThis + VarOther;
                    break;

                case ProbDistOp::pdoMinus:
                    variable = VarThis - VarOther;
                    break;

                case ProbDistOp::pdoMult:
                    variable = VarThis * VarOther;
                    break;
                    
                case ProbDistOp::pdoDiv:
                    variable = VarThis / VarOther;
                    break;
                    
                case ProbDistOp::pdoUnknown:
                case ProbDistOp::pdoLast:                    
                    break;
            }
            
            
//             cout << "VarThis " << VarThis.RawPrint(10) << endl;
//             cout << "VarOther " << VarOther.RawPrint(10) << endl;
//             cout << "VarEndThis " << VarEndThis.RawPrint(10) << endl;
//             cout << "VarEndOther " << VarEndOther.RawPrint(10) << endl;
//             cout << "IntegrationIntervallThis " << IntegrationIntervallThis.RawPrint(10) << endl;
//             cout << "IntegrationIntervallOther " << IntegrationIntervallOther.RawPrint(10) << endl;
            
            
            // calculate the probability
            CDigFloat probThis, probOther;
            probThis = AbsIntegral(VarThis-IntegrationIntervallThis/2., VarThis+IntegrationIntervallThis/2.);
            probOther = other.AbsIntegral(VarOther- IntegrationIntervallOther/2., VarOther+IntegrationIntervallOther/2.);
//             cout << "probThis " << probThis.RawPrint(10) << endl;
//             cout << "probOther " << probOther.RawPrint(10) << endl;
//             cout << "[" << variable.RawPrint(10) << "]=" << (probOther*probThis).RawPrint(10) << endl;
            
            // depending on the fact if the variable already exists: 
            // a) if exists:  sum new value to the existing
            // b) if not exists: add as new element
            if( ResultDistri.find(variable) == ResultDistri.end() )
                ResultDistri[variable] = CDigFloat(probThis)*probOther;
            else 
                ResultDistri[variable] += CDigFloat(probThis)*probOther;
            
            // increment variable by integration intervall for other distri
            VarOther += IntegrationIntervallOther;
            
        }   // endwhile( VarOther < VarEndOther )
        
        // increment variable by integration intervall for this distri
        VarThis += IntegrationIntervallThis;
        
    }   // endwhile( VarThis < VarEndThis )
      
//     cout << "ResultDistr: " << endl;
//     for(auto iel: ResultDistri)
//         cout << "[" << iel.first.RawPrint(10) << "] =" << iel.second.RawPrint(10) << endl;
    
    
    // normalize again:
    // --> dfAbsIntegral == 1
    mDistribution.clear();
    for(auto iel: ResultDistri)
        mDistribution[iel.first] = iel.second;
    bNormalized=false;
    _Normalize();
    
    assert(AbsIntegral() == 1);
    
//     cout << Print(10);
    
    return *this;

}
CProbabilityDensityDistribution& CProbabilityDensityDistribution::_Convolution4GeneralOperators(CProbabilityDensityDistribution& other, const ProbDistOp Operation)
{
    // not nice but necessary: normalize the other, too: this is the state all the distribution should be in
    _Normalize();
    other._Normalize();
    
    // get the min distance for convolution stepping
    CDigFloat MinVarStepThis = _GetMinVariableDifference(Distribution());
    CDigFloat MinVarStepOther = _GetMinVariableDifference(other.Distribution());
    
    // declare the result distri
    M_DFDF TargetDistri;
    
    // declare the variables for both distris to step through
    CDigFloat VarThis = CDigFloat( Distribution().begin()->first );
    CDigFloat VarOther= CDigFloat( other.Distribution().begin()->first );
    
    // the final variable values: end outside the distribution 
    CDigFloat VarEndThis = CDigFloat( prev(Distribution().end())->first );
    CDigFloat VarEndOther = CDigFloat( prev(other.Distribution().end())->first );
    
    // calculate the target range for the resulting variable
    CDigFloat TargetRangeStart, TargetRangeEnd;
    _GetRangeFromDistributionOperation(other.Distribution(), Operation, TargetRangeStart, TargetRangeEnd);
    
    // calculate the target step size: 
    int nTargetSteps = (int) (TargetRangeEnd.RawValue() - TargetRangeStart.RawValue()) / min(MinVarStepOther.RawValue(), MinVarStepThis.RawValue()) +1;
    CDigFloat TargetStepSize = (TargetRangeEnd - TargetRangeStart) / nTargetSteps;  
    CDigFloat IntegrationStepSize = TargetStepSize / 10.;
    CDigFloat TargetVariable = TargetRangeStart;
    
//     cout << "_Convolution4GeneralOperators( other distri:" << endl << other.Print(10) << endl <<  ", this distri:" << endl << Print(10) << endl << ", operator: " << (int)Operation << " ) " << endl;
//     cout << "TargetRangeStart " << TargetRangeStart.RawPrint(10) << endl;
//     cout << "TargetRangeEnd " << TargetRangeEnd.RawPrint(10) << endl;
//     cout << "nTargetSteps " << nTargetSteps << endl;
//     cout << "TargetStepSize " << TargetStepSize.RawPrint(10) << endl;
//     cout << "IntegrationStepSize " << IntegrationStepSize.RawPrint(10) << endl;
    
    // do only continue if there is no division by zero
    if(!(TargetRangeStart == TargetRangeEnd && TargetRangeStart == 0))
    {
        while(TargetVariable <= TargetRangeEnd)
        {
            // init the probability density
            TargetDistri[TargetVariable]=0;
            
            // initialize integration variable
            VarOther = other.Distribution().begin()->first;
            
            // integrate over the multiplied distributions with given stepsize
            while(VarOther <= VarEndOther)
            {
                // set the integration factor due to the parametrization of the integral:
                // e.g. for multiplication: x1*x2 = const (=integration varibale of target distribution)
                // --> t = const --> x1 = t / x2
                // dx1 / dt = 1/x2 --> dx1 = dt / x2 
                // --> the integration factor in this case is 1/x2
                CDigFloat IntegrationFactor = 1;
                
                // calculate the integration deltas for correct integration limits
                // simple case: x1 +- IntegrationStepSize /2.
                CDigFloat NegIntegrationDeltaOther = IntegrationStepSize/2.;
                CDigFloat PosIntegrationDeltaOther = IntegrationStepSize/2.;
                CDigFloat NegIntegrationDeltaThis = IntegrationStepSize/2.;
                CDigFloat PosIntegrationDeltaThis = IntegrationStepSize/2.;
                // set the integration variable for this distribution
                switch(Operation)
                {
                    case ProbDistOp::pdoPlus:
                        VarThis = TargetVariable - VarOther;
                        break;

                    case ProbDistOp::pdoMinus:
                        VarThis = TargetVariable + VarOther;
                        break;

                    case ProbDistOp::pdoMult:
                        // avoid division by zero
                        if( abs(VarOther) > IntegrationStepSize/2. )
                        {
                            // this is the "normale case" no division by zero
                            VarThis = TargetVariable / VarOther;
                            IntegrationFactor = abs(CDigFloat(1.)/VarOther);
                            
                            // the integration limits are non-linear due to the dependency 1/VarOther of the integration variable for this distribution
                            NegIntegrationDeltaThis = abs(TargetVariable / ( VarOther - IntegrationStepSize /2.) - TargetVariable / VarOther );
                            PosIntegrationDeltaThis = abs(TargetVariable / VarOther - TargetVariable / ( VarOther + IntegrationStepSize /2.));
                        }
                        else
                        {
                            // in case of VarOther being zero: 
                            // the TargetVariable is zero (t  = VarThis * VarOther ) : if not --> continue
                            // the probability added to the target distribution[0] += :
                            // is the integral around 0 of other distribution multiplied by the total distribution integral of other (==1):
                            // this is due to the fact that all variables of this distribution multiplied with VarOther == 0 --> 0 which is 
                            // the value of the target variable
                            if(abs(TargetVariable) > IntegrationStepSize/2.)
                            {
                                // incerment: goto next value
                                VarOther += IntegrationStepSize;
                                continue;
                            }
                            
                            // coming here means: we are around zero with VarOther and TargetVariable:
                            // so the effective probability is This.Distribution(0) * 1
                            VarThis = 0;
                            NegIntegrationDeltaOther = VarOther - other.Distribution().begin()->first;
                            PosIntegrationDeltaOther = VarEndOther - VarOther;
                        }                     
                        break;
                        
                    case ProbDistOp::pdoDiv:
                        VarThis = TargetVariable * VarOther;
                        IntegrationFactor = abs(VarOther);
                        break;
                        
                    case ProbDistOp::pdoUnknown:
                    case ProbDistOp::pdoLast:                    
                        break;
                }
                
                
                // check if this is out of bounds ...
                if(OutOfBounds(VarThis-NegIntegrationDeltaThis) && OutOfBounds(VarThis + PosIntegrationDeltaThis))
                {
                    // debug:
//                     cout << "VarThis out of bounds ... leaving " << VarThis.RawPrint(10) <<  endl;
                    // ... leave if so : 
                    VarOther += IntegrationStepSize;
                    continue;
                }
                
                // DEBUG
//                 cout << "TargetVariable " << TargetVariable.RawPrint(10) << endl;
//                 cout << "VarThis " << VarThis.RawPrint(10) << endl;
//                 cout << "VarOther " << VarOther.RawPrint(10) << endl;
//                 cout << "NegIntegrationDeltaThis " << NegIntegrationDeltaThis.RawPrint(10) << endl;
//                 cout << "PosIntegrationDeltaThis " << PosIntegrationDeltaThis.RawPrint(10) << endl;
//                 cout << "NegIntegrationDeltaOther " << NegIntegrationDeltaOther.RawPrint(10) << endl;
//                 cout << "PosIntegrationDeltaOther " << PosIntegrationDeltaOther.RawPrint(10) << endl;
                
                // calculate the probability
                CDigFloat probThis, probOther;
                probThis = AbsIntegral(VarThis-NegIntegrationDeltaThis, VarThis+PosIntegrationDeltaThis);
                probOther = other.AbsIntegral(VarOther- NegIntegrationDeltaOther, VarOther+PosIntegrationDeltaOther);
                
                // add to result
                TargetDistri[TargetVariable] += IntegrationFactor * probOther*probThis;
                
                // DEBUG
//                 cout << "probThis " << probThis.RawPrint(10) << endl;
//                 cout << "probOther " << probOther.RawPrint(10) << endl;
//                 cout << "IntegrationFactor " << IntegrationFactor.RawPrint(10) << endl;
//                 cout << "TargetDistri[TargetVariable] " << TargetDistri[TargetVariable].RawPrint(10) << endl;
                
                // increment integration variable by step size
                VarOther += IntegrationStepSize;
                
            }   // endwhile(VarOther <= VarEndOther)
        // increment the target variable
        TargetVariable += TargetStepSize;
        
        }   //endwhile(TargetVariable <= TargetRangeEnd)
        
    }   // endif(!(TargetRangeStart == TargetRangeEnd && TargetRangeStart == 0))
    
    
    // normalize again:
    // --> dfAbsIntegral == 1
    _Init();
    
    //
    if(!(TargetRangeStart == TargetRangeEnd && TargetRangeStart == 0))
    {
        for(auto iel: TargetDistri)
            mDistribution[iel.first] = iel.second;
        
        
//         cout << "before Normalization: " << endl << Print(10);
        
        
        bNormalized=false;
        _Normalize();
        
//         cout << "after Normalization: " << endl << Print(10);
        assert(AbsIntegral() == 1);
        
    }   // endif(!(TargetRangeStart == TargetRangeEnd && TargetRangeStart == 0))
    
    return *this;
}

void CProbabilityDensityDistribution::_GetRangeFromDistributionOperation(const M_DFDF& other, const ProbDistOp Operation, CDigFloat& TargetRangeStart, CDigFloat& TargetRangeEnd)
{
    CDigFloat ThisStart = Distribution().begin()->first;
    CDigFloat ThisEnd = prev(Distribution().end())->first;
    CDigFloat OtherStart = other.begin()->first;
    CDigFloat OtherEnd = prev(other.end())->first;
    
    switch(Operation)
    {
        case ProbDistOp::pdoPlus:
            TargetRangeStart = ThisStart + OtherStart;
            TargetRangeEnd = ThisEnd + OtherEnd;
            break;
            
        case ProbDistOp::pdoMinus:
            TargetRangeStart = ThisStart - OtherStart;
            TargetRangeEnd = ThisEnd + OtherEnd;
            break;
            
        case ProbDistOp::pdoMult:
            TargetRangeStart = min(min(ThisStart*OtherStart, ThisStart*OtherEnd),
                               min(ThisEnd*OtherStart,ThisEnd*OtherEnd));
            TargetRangeEnd = max(max(ThisStart*OtherStart, ThisStart*OtherEnd),
                               max(ThisEnd*OtherStart,ThisEnd*OtherEnd));
            break;
            
        case ProbDistOp::pdoDiv :
            
            // only set, if divisor is non-zero in every case:
            if(!(OtherStart <= 0 && OtherEnd >=0))
            {
                TargetRangeStart = min(min(ThisStart/OtherStart, ThisStart/OtherEnd),
                                min(ThisEnd/OtherStart,ThisEnd/OtherEnd));
                TargetRangeEnd = max(max(ThisStart/OtherStart, ThisStart/OtherEnd),
                                max(ThisEnd/OtherStart,ThisEnd/OtherEnd));
            }
            else
            {
                TargetRangeStart = 0;
                TargetRangeEnd = 0;
            }
            break;
            
    }   //endswitch(Operation)
}
