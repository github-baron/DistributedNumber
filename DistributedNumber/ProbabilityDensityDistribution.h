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

#ifndef CPROBABILITYDENSITYDISTRIBUTION_H
#define CPROBABILITYDENSITYDISTRIBUTION_H

// includes
#include "Distribution.h"

/**
 * @brief This class represents a probability density distribution.
 */
class CProbabilityDensityDistribution :  CDistribution
{
public:
    /**
     * Default constructor
     */
    CProbabilityDensityDistribution();

    /**
     * Copy constructor
     *
     * @param other TODO
     */
    CProbabilityDensityDistribution(const CProbabilityDensityDistribution& other);

    /**
     * Destructor
     */
    ~CProbabilityDensityDistribution();

    /**
     * Assignment operator
     *
     * @param other TODO
     * @return TODO
     */
    CProbabilityDensityDistribution& operator=(const CProbabilityDensityDistribution& other);

    /**
     * @todo write docs
     *
     * @param other TODO
     * @return TODO
     */
    bool operator==(const CProbabilityDensityDistribution& other) const;

    /**
     * @todo write docs
     *
     * @param other TODO
     * @return TODO
     */
    bool operator!=(const CProbabilityDensityDistribution& other) const;

    /**
     * Copy constructor
     *
     * @param other TODO
     */
    CProbabilityDensityDistribution(const CDistribution& other);
    
       /**
     * @brief adds new distribution element if distribution is not already normalized (see ::_Normalize).
     *        In case of adding the last element the surrounding elements will be set to zero (... if not done)
     *        and the whole distribution will be normalized (see ::_Normalize).
     *
     * @param variable CDigFloat& distribution variable
     * @param value CDigFloat&
     * @param bLastElement bool if set --> ::_Normalize is called and no further element can be added
     * @return bool
     */
    bool Add(const CDigFloat& variable, const CDigFloat value, bool bLastElement = false); 
  
protected:    
    /**
     * @brief normalize: divides by norm (integral over the distribution) 
     *
     * @param other CDistribution&
     * @return bool
     */
    void _Normalize();  
    /**
     * @brief initializes member
     *
     */
    void _Init();
    
    
    bool bNormalized;
    CDigFloat dfAbsIntegral;

};

#endif // CPROBABILITYDENSITYDISTRIBUTION_H
