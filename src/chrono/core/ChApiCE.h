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

#ifndef CHAPICE_H
#define CHAPICE_H

#include <cstddef>
#include <cassert>
#include <iostream>
#include <array>

#include "chrono/ChVersion.h"
#include "chrono/core/ChPlatform.h"

// When compiling this library, remember to define CH_API_COMPILE
// (so that the symbols with 'ChApi' in front of them will be
// marked as exported). Otherwise, just do not define it if you
// link the library to your code, and the symbols will be imported.

#if defined(CH_API_COMPILE)
    #define ChApi ChApiEXPORT
#else
    #define ChApi ChApiIMPORT
#endif

// Define CH_NODISCARD macro for C++17 and later
#if __cplusplus >= 201703L
    #define CH_NODISCARD [[nodiscard]]
#else
    #define CH_NODISCARD
#endif


// Helper to detect if a type is constexpr constructible
template <typename T, int = (T{}, 0)>
constexpr bool is_constexpr_constructible (int)
{ return true; }

template <typename>
constexpr bool is_constexpr_constructible (long)
{ return false; }

// Type aliases
template<class T, unsigned L>
using CH_ARRAY = std::array<T, L>;

#endif
