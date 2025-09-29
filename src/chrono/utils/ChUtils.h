// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================

#ifndef CH_UTILS_H
#define CH_UTILS_H

#include <cstdio>
#include <stdexcept>
#include <algorithm>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/core/ChApiCE.h"
#include "chrono/utils/ChConstants.h"
#include "chrono/core/ChTypes.h"

namespace chrono {

#ifndef SOURCE_PATH_SIZE
    #define SOURCE_PATH_SIZE 0
#endif

#define __FILENAME__ ((__FILE__) + (SOURCE_PATH_SIZE))

#ifdef DEBUG_LOG
    #define ChDebugLog(x)                                                                       \
        do {                                                                                    \
            std::cerr << "[DBG " << __FILENAME__ << "::" << __func__ << "] " << x << std::endl; \
        } while (0)
#else
    #define ChDebugLog(x)
#endif

#define ChAssertAlways(exp)                                                                                    \
    {                                                                                                          \
        if (!(exp)) {                                                                                          \
            char msg[300];                                                                                     \
            std::sprintf(msg, "Expression '%s' returned false - file %s, line %d.", #exp, __FILE__, __LINE__); \
            std::cerr << msg << std::endl;                                                                     \
            throw std::runtime_error(msg);                                                                     \
        }                                                                                                      \
    }

/// Clamp and modify the specified value to lie within the given limits.
template <typename T>
void ChClampValue(T& value, T limitMin, T limitMax) {
    if (value < limitMin)
        value = limitMin;
    else if (value > limitMax)
        value = limitMax;
}

/// Clamp the specified value to lie within the given limits.
template <typename T>
T ChClamp(T value, T limitMin, T limitMax) {
    if (value < limitMin)
        return limitMin;
    if (value > limitMax)
        return limitMax;

    return value;
}

/// Signum function.
template <typename T>
int ChSignum(T x) {
    return (x > T(0)) - (x < T(0));
}

/// Wrap angle in range.
/// If symmetric, wrap in [-PI; PI).
/// If not symmetric, wrap in [0; 2PI).
template <typename T>
T ChWrapAngle(T angle, bool symmetric = true) {
    T wangle = angle;
    if (symmetric) {  // [-PI; +PI)
        wangle = std::fmod(wangle + CH_PI, CH_2PI);
        wangle < 0 ? wangle += CH_PI : wangle -= CH_PI;
    } else {  // [0; 2PI)
        wangle = std::fmod(wangle, CH_2PI);
        if (wangle < 0)
            wangle += CH_2PI;
    }
    return wangle;
}

namespace impl {
    /** @brief Implementation of the Newton-Raphson square root calculation
     *
     *  @return the square root value, as a double
     */
    CH_NODISCARD static constexpr inline double
        sqrt_newton_raphson(const double& x, const double& current, const double& previous) {
        return current == previous ? current : sqrt_newton_raphson(x, 0.5 * (current + x / current), current);
    }
}  // namespace impl

/** @brief Constexpr square root implementation using the Newton-Raphson method
 *
 *  @param x the value
 *
 *  @return the square root of the provided value
 */
CH_NODISCARD static constexpr inline double cesqrt(const double& x) {
    return x >= 0.0 && x < std::numeric_limits<double>::infinity() ? impl::sqrt_newton_raphson(x, x, 0.0) :
                                                                        std::numeric_limits<double>::quiet_NaN();
}

}  // end namespace chrono

#endif
