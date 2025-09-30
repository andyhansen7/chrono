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
// Authors: Alessandro Tasora
// =============================================================================

#ifndef CHQUADRATURE
#define CHQUADRATURE

#include <vector>

#include "chrono/core/ChApiCE.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChTypes.h"
#include "chrono/utils/ChConstants.h"

namespace chrono {

/// Polynomial roots and weights for the Gauss-Legendre quadrature.
/// These quadrature tables are automatically managed by ChQuadrature.
class ChApi ChQuadratureTables {
  public:
    ChQuadratureTables(const int& order_from = 1, const int& order_to = 10);

    std::vector<std::vector<double>> Weight;
    std::vector<std::vector<double>> Lroots;

    void PrintTables() const;

  private:
    void glege_roots(ChMatrixDynamic<>& lcoef, const int& N, const int& ntable);
};

/// Polynomial roots and weights for quadrature over a triangle.
/// Assumes 2 natural 'area' coordinates u and v, with the 3rd assumed to be 1-u-v.
/// D. A. Dunavant, "High degree efficient symmetrical Gaussian quadrature rules for the triangle",
/// Int. J. Num. Meth. Engng, 21(1985):1129-1148.
class ChApi ChQuadratureTablesTriangle {
  public:
    ChQuadratureTablesTriangle();

    std::vector<std::vector<double>> Weight;
    std::vector<std::vector<double>> LrootsU;
    std::vector<std::vector<double>> LrootsV;
};

/// Polynomial roots and weights for quadrature over a tetrahedron.
/// Assumes 3 natural 'volume' coordinates u, v, and w, with the 4th assumed to be 1-u-v-w.
class ChApi ChQuadratureTablesTetrahedron {
  public:
    ChQuadratureTablesTetrahedron();

    std::vector<std::vector<double>> Weight;
    std::vector<std::vector<double>> LrootsU;
    std::vector<std::vector<double>> LrootsV;
    std::vector<std::vector<double>> LrootsW;
};

// -----------------------------------------------------------------------------

/// Base class for 1D integrand T = f(x) to be used in ChQuadrature.
/// The class is templated so that the computed valued can be either a scalar or a more complex object (e.g., a matrix).
template <class T = double>
class ChIntegrand1D {
  public:
    virtual ~ChIntegrand1D() {}

    /// Evaluate the function at point x, that is result T = f(x).
    virtual void Evaluate(T& result, const double& x) = 0;
};

/// Base class for 2D integrand T = f(x,y) to be used in ChQuadrature.
template <class T = double>
class ChIntegrand2D {
  public:
    virtual ~ChIntegrand2D() {}

    /// Evaluate the function at point x,y , that is result T = f(x,y).
    virtual void Evaluate(T& result, const double& x, const double& y) = 0;
};

/// Base class for 3D integrand T = f(x,y,z) to be used in ChQuadrature.
template <class T = double>
class ChIntegrand3D {
  public:
    virtual ~ChIntegrand3D() {}

    /// Evaluate the function at point x,y,z , that is result T = f(x,y,z)
    virtual void Evaluate(T& result, const double& x, const double& y, const double& z) = 0;
};

// -----------------------------------------------------------------------------

/// Gauss-Legendre quadrature in 1D, 2D, or 3D.
/// This class provides methods to integrate a function (the integrand) on a domain using the Gauss quadrature.
/// This method is most useful for polynomial integrands, in which case the result is exact if the order of quadrature
/// is greater or equal to the degree of the polynomial.
class ChApi ChQuadrature {
  public:
    /// Integrate the integrand T = f(x) over the 1D interval [xA, xB], with specified order of quadrature.
    /// Best if integrand is polynomial. For order in 1..10, precomputed polynomial coefficients are used.
    template <class T>
    static void Integrate1D(T& result,                    ///< result is returned here
                            ChIntegrand1D<T>& integrand,  ///< this is the integrand
                            const double& x_min,           ///< min limit for x domain
                            const double& x_max,           ///< min limit for x domain
                            const int& order               ///< order of integration
    ) {
        std::shared_ptr<ChQuadratureTables> mtables = nullptr;
        const std::vector<double>* lroots;
        const std::vector<double>* weight;

        if (static_cast<unsigned>(order) <= GetStaticTables()->Lroots.size()) {
            mtables = GetStaticTables();
            lroots = &mtables->Lroots[order - 1];
            weight = &mtables->Weight[order - 1];
        } else {
            mtables = chrono_types::make_shared<ChQuadratureTables>(order, order);
            mtables->PrintTables();
            lroots = &mtables->Lroots[0];
            weight = &mtables->Weight[0];
        }

        const double c1 = (x_max - x_min) / 2;
        const double c2 = (x_max + x_min) / 2;

        result *= 0;  // as result = 0, but works also for matrices.
        T val;        // temporary value for loop

        for (unsigned int i = 0; i < lroots->size(); i++) {
            integrand.Evaluate(val, (c1 * lroots->at(i) + c2));
            val *= weight->at(i);
            result += val;
        }
        result *= c1;  // result = c1 * sum;
    }

    /// Integrate the integrand T = f(x,y) over the 2D interval [xA, xB][yA, yB], with desired order of quadrature.
    /// Best if integrand is polynomial. For order in 1..10, precomputed polynomial coefficients are used.
    template <class T>
    static void Integrate2D(T& result,                    ///< result is returned here
                            ChIntegrand2D<T>& integrand,  ///< this is the integrand
                            const double& x_min,           ///< min limit for x domain
                            const double& x_max,           ///< min limit for x domain
                            const double& y_min,           ///< min limit for y domain
                            const double& y_max,           ///< min limit for y domain
                            const int& order               ///< order of integration
    ) {
        std::shared_ptr<ChQuadratureTables> mtables = nullptr;
        const std::vector<double>* lroots;
        const std::vector<double>* weight;

        if (static_cast<unsigned int>(order) <= GetStaticTables()->Lroots.size()) {
            mtables = GetStaticTables();
            lroots = &mtables->Lroots[order - 1];
            weight = &mtables->Weight[order - 1];
        } else {
            mtables = chrono_types::make_shared<ChQuadratureTables>(order, order);
            mtables->PrintTables();
            lroots = &mtables->Lroots[0];
            weight = &mtables->Weight[0];
        }

        const double Xc1 = (x_max - x_min) / 2;
        const double Xc2 = (x_max + x_min) / 2;
        const double Yc1 = (y_max - y_min) / 2;
        const double Yc2 = (y_max + y_min) / 2;

        result *= 0;  // as result = 0, but works also for matrices.
        T val;        // temporary value for loop

        for (unsigned int ix = 0; ix < lroots->size(); ix++)
            for (unsigned int iy = 0; iy < lroots->size(); iy++) {
                integrand.Evaluate(val, (Xc1 * lroots->at(ix) + Xc2), (Yc1 * lroots->at(iy) + Yc2));
                val *= (weight->at(ix) * weight->at(iy));
                result += val;
            }
        result *= (Xc1 * Yc1);
    }

    /// Integrate the integrand T = f(x,y,z) over the 3D interval [xA, xB][yA, yB][zA, zB], with desired order of
    /// quadrature. Best if integrand is polynomial. For order in 1..10, precomputed polynomial coefficients are used.
    template <class T>
    static void Integrate3D(T& result,                    ///< result is returned here
                            ChIntegrand3D<T>& integrand,  ///< this is the integrand
                            const double x_min,           ///< min limit for x domain
                            const double x_max,           ///< min limit for x domain
                            const double y_min,           ///< min limit for y domain
                            const double y_max,           ///< min limit for y domain
                            const double z_min,           ///< min limit for z domain
                            const double z_max,           ///< min limit for z domain
                            const int order               ///< order of integration
    ) {
        std::shared_ptr<ChQuadratureTables> mtables = nullptr;
        const std::vector<double>* lroots;
        const std::vector<double>* weight;

        if (static_cast<unsigned int>(order) <= GetStaticTables()->Lroots.size()) {
            mtables = GetStaticTables();
            lroots = &mtables->Lroots[order - 1];
            weight = &mtables->Weight[order - 1];
        } else {
            mtables = chrono_types::make_shared<ChQuadratureTables>(order, order);
            mtables->PrintTables();
            lroots = &mtables->Lroots[0];
            weight = &mtables->Weight[0];
        }

        const double Xc1 = (x_max - x_min) / 2;
        const double Xc2 = (x_max + x_min) / 2;
        const double Yc1 = (y_max - y_min) / 2;
        const double Yc2 = (y_max + y_min) / 2;
        const double Zc1 = (z_max - z_min) / 2;
        const double Zc2 = (z_max + z_min) / 2;

        result *= 0;  // as result = 0, but works also for matrices.
        T val;        // temporary value for loop

        for (unsigned int ix = 0; ix < lroots->size(); ix++)
            for (unsigned int iy = 0; iy < lroots->size(); iy++)
                for (unsigned int iz = 0; iz < lroots->size(); iz++) {
                    integrand.Evaluate(val, (Xc1 * lroots->at(ix) + Xc2), (Yc1 * lroots->at(iy) + Yc2),
                                       (Zc1 * lroots->at(iz) + Zc2));
                    val *= (weight->at(ix) * weight->at(iy) * weight->at(iz));
                    result += val;
                }
        result *= (Xc1 * Yc1 * Zc1);
    }

    /// Integrate the 2D integrand T = f(u,v) over a triangle, with desired order of quadrature.
    /// Best if integrand is polynomial. Two triangle coordinates are assumed to be 'area' coordinates u and v in
    /// [0...1], with the 3rd assumed to be 1-u-v.   For order between 1 and 5, use precomputed polynomial coefficients.
    template <class T>
    static void Integrate2Dtriangle(T& result,                    ///< result is returned here
                                    ChIntegrand2D<T>& integrand,  ///< this is the integrand
                                    const int order               ///< order of integration
    ) {
        if (static_cast<unsigned int>(order) > GetStaticTablesTriangle()->Weight.size())
            throw std::invalid_argument("Too high order of quadrature for triangle. Use lower order.");

        const auto mtables = GetStaticTablesTriangle();
        const std::vector<double>* lrootsU = &mtables->LrootsU[order - 1];
        const std::vector<double>* lrootsV = &mtables->LrootsV[order - 1];
        const std::vector<double>* weight = &mtables->Weight[order - 1];

        result *= 0;  // as result = 0, but works also for matrices.
        T val;        // temporary value for loop

        for (unsigned int i = 0; i < weight->size(); i++) {
            integrand.Evaluate(val, lrootsU->at(i), lrootsV->at(i));
            val *= weight->at(i) * 0.5;  // the 1/2 coefficient is not in the table
            result += val;
        }
    }

    /// Integrate the 3D integrand T = f(u,v,w) over a tetrahedron, with desired order of quadrature.
    /// Best if integrand is polynomial. Three tetrahedron coordinates are assumed to be 'volume' coordinates u, v, and
    /// w in [0...1], with the 4th assumed to be1-u-v-w.   For order between 1 and 5 use precomputed polynomial
    /// coefficients.
    template <class T>
    static void Integrate3Dtetrahedron(T& result,                    ///< result is returned here
                                       ChIntegrand3D<T>& integrand,  ///< this is the integrand
                                       const int order               ///< order of integration
    ) {
        if (static_cast<unsigned int>(order) > GetStaticTablesTetrahedron()->Weight.size())
            throw std::invalid_argument("Too high order of quadrature for tetrahedron. Use lower order.");

        const auto mtables = GetStaticTablesTetrahedron();
        const std::vector<double>* lrootsU = &mtables->LrootsU[order - 1];
        const std::vector<double>* lrootsV = &mtables->LrootsV[order - 1];
        const std::vector<double>* lrootsW = &mtables->LrootsW[order - 1];
        const std::vector<double>* weight = &mtables->Weight[order - 1];

        result *= 0;  // as result = 0, but works also for matrices.
        T val;        // temporary value for loop

        for (unsigned int i = 0; i < weight->size(); i++) {
            integrand.Evaluate(val, lrootsU->at(i), lrootsV->at(i), lrootsW->at(i));
            val *= weight->at(i) * CH_1_6;  // the 1/6 coefficient is not in the table
            result += val;
        }
    }

    /// Access a statically-allocated set of tables, from 0 to 10th order, with precomputed tables.
    CH_NODISCARD static inline std::shared_ptr<ChQuadratureTables> GetStaticTables();

    /// Access a statically-allocated set of tables for tetrahedron quadrature, with 5 precomputed tables.
    CH_NODISCARD static inline std::shared_ptr<ChQuadratureTablesTriangle> GetStaticTablesTriangle();

    /// Access a statically-allocated set of tables for tetrahedron quadrature, with 5 precomputed tables.
    /// Use Dunavant theory.
    CH_NODISCARD static inline std::shared_ptr<ChQuadratureTablesTetrahedron> GetStaticTablesTetrahedron();
};

// -----------------------------------------------------------------------------

#define CH_QUADRATURE_STATIC_TABLES 10

static std::shared_ptr<ChQuadratureTables> static_tables;
static std::shared_ptr<ChQuadratureTablesTriangle> static_tables_triangle;  // only 5 tables
static std::shared_ptr<ChQuadratureTablesTetrahedron> static_tables_tetrahedron;  // only 5 tables

// -----------------------------------------------------------------------------

CH_NODISCARD inline std::shared_ptr<ChQuadratureTables> ChQuadrature::GetStaticTables() {
    if (!static_tables) {
        static_tables = chrono_types::make_shared<ChQuadratureTables>(1, CH_QUADRATURE_STATIC_TABLES);
    }
    return static_tables;
}

CH_NODISCARD inline std::shared_ptr<ChQuadratureTablesTriangle> ChQuadrature::GetStaticTablesTriangle() {
    if (!static_tables_triangle) {
        static_tables_triangle = chrono_types::make_shared<ChQuadratureTablesTriangle>();
    }
    return static_tables_triangle;
}

CH_NODISCARD inline std::shared_ptr<ChQuadratureTablesTetrahedron> ChQuadrature::GetStaticTablesTetrahedron() {
    if (!static_tables_tetrahedron) {
        static_tables_tetrahedron = chrono_types::make_shared<ChQuadratureTablesTetrahedron>();
    }
    return static_tables_tetrahedron;
}

}  // end namespace chrono

#endif
