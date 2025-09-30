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
// Authors: Alessandro Tasora, Radu Serban
// =============================================================================

#ifndef CH_QUATERNION_H
#define CH_QUATERNION_H

#include <cmath>
#include <algorithm>
#include <limits>

#include "chrono/core/ChApiCE.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChVector3.h"
#include "chrono/utils/ChConstants.h"
#include "chrono/core/ChTypes.h"
#include "chrono/utils/ChUtils.h"

namespace chrono {

/// Class defining quaternion objects, that is four-dimensional numbers.
/// Unit quaternions can be used to represent rotations in 3d (in which case, they are also known as "Euler
/// parameters").
/// For more details, see the @ref manual_ChQuaternion manual page.
template <class Real = double>
class ChQuaternion {
  using ARRAY_TYPE = chrono_types::ChArray<Real, 4>;
  public:
    /// Default constructor.
    /// Note that this constructs a null quaternion {0,0,0,0}, not a {1,0,0,0} unit quaternion.
    constexpr inline ChQuaternion();

    /// Constructor from four scalars.
    /// The first is the real part, others are i,j,k imaginary parts
    constexpr inline ChQuaternion(const Real& e0, const Real& e1, const Real& e2, const Real& e3);

    /// Constructor from real part, and vector with i,j,k imaginary part.
    constexpr inline ChQuaternion(const Real& s, const ChVector3<Real>& v);

    /// Copy constructor.
    constexpr inline ChQuaternion(const ChQuaternion<Real>& other);

    /// Copy constructor with type change.
    template <class RealB>
    constexpr inline ChQuaternion(const ChQuaternion<RealB>& other);

    /// Access to components
    CH_NODISCARD constexpr inline Real& e0() { return m_data[0]; }
    CH_NODISCARD constexpr inline Real& e1() { return m_data[1]; }
    CH_NODISCARD constexpr inline Real& e2() { return m_data[2]; }
    CH_NODISCARD constexpr inline Real& e3() { return m_data[3]; }
    CH_NODISCARD constexpr inline const Real& e0() const { return m_data[0]; }
    CH_NODISCARD constexpr inline const Real& e1() const { return m_data[1]; }
    CH_NODISCARD constexpr inline const Real& e2() const { return m_data[2]; }
    CH_NODISCARD constexpr inline const Real& e3() const { return m_data[3]; }

    /// Access to underlying array storage.
    CH_NODISCARD constexpr inline ARRAY_TYPE& data() { return m_data; }
    CH_NODISCARD constexpr inline const ARRAY_TYPE& data() const { return m_data; }

    // EIGEN INTER-OPERABILITY

    /// Construct a quaternion from an Eigen vector expression.
    template <typename Derived>
    inline ChQuaternion(const Eigen::MatrixBase<Derived>& vec,
                 typename std::enable_if<(Derived::MaxRowsAtCompileTime == 1 || Derived::MaxColsAtCompileTime == 1),
                                         Derived>::type* = 0) {
        m_data[0] = vec(0);
        m_data[1] = vec(1);
        m_data[2] = vec(2);
        m_data[3] = vec(3);
    }

    /// View this quaternion as an Eigen vector.
    CH_NODISCARD inline Eigen::Map<Eigen::Matrix<Real, 4, 1>> eigen() { return Eigen::Map<Eigen::Matrix<Real, 4, 1>>(m_data.data()); }
    CH_NODISCARD inline Eigen::Map<const Eigen::Matrix<Real, 4, 1>> eigen() const {
        return Eigen::Map<const Eigen::Matrix<Real, 4, 1>>(m_data.data());
    }

    /// Assign an Eigen vector expression to this quaternion.
    template <typename Derived>
    inline ChQuaternion& operator=(const Eigen::MatrixBase<Derived>& vec) {
        m_data[0] = vec(0);
        m_data[1] = vec(1);
        m_data[2] = vec(2);
        m_data[3] = vec(3);
        return *this;
    }

    // SET & GET FUNCTIONS

    /// Sets the four values of the quaternion at once
    inline void Set(const Real& e0, const Real& e1, const Real& e2, const Real& e3);

    /// Sets the quaternion as a copy of another quaternion
    inline void Set(const ChQuaternion<Real>& q);

    /// Sets the quaternion with four components as a sample scalar
    inline void Set(const Real& s);

    /// Sets the quaternion as a null quaternion
    inline void SetNull();

    /// Sets the quaternion as a unit quaternion
    inline void SetUnit();

    /// Sets the scalar part only
    inline void SetScalar(const Real& s);

    /// Sets the vectorial part only
    inline void SetVector(const ChVector3<Real>& v);

    /// Return true if this quaternion is the null quaternion.
    CH_NODISCARD constexpr inline bool IsNull() const;

    /// Return true if this quaternion is the identity quaternion.
    CH_NODISCARD constexpr inline bool IsIdentity() const;

    /// Return true if quaternion is identical to other quaternion
    CH_NODISCARD constexpr inline bool Equals(const ChQuaternion<Real>& other) const;

    /// Return true if quaternion equals another quaternion, within a tolerance 'tol'
    CH_NODISCARD constexpr inline bool Equals(const ChQuaternion<Real>& other, const Real& tol) const;

    /// Gets the vectorial part only
    CH_NODISCARD constexpr inline ChVector3<Real> GetVector() const;

    /// Get the X axis of a coordinate system, given the quaternion which its alignment.
    /// It is assumed that the quaternion is normalized (and thus represents a rotation).
    CH_NODISCARD constexpr inline ChVector3<Real> GetAxisX() const;

    /// Get the Y axis of a coordinate system, given the quaternion which its alignment.
    /// It is assumed that the quaternion is normalized (and thus represents a rotation).
    CH_NODISCARD constexpr inline ChVector3<Real> GetAxisY() const;

    /// Get the Z axis of a coordinate system, given the quaternion which its alignment.
    /// It is assumed that the quaternion is normalized (and thus represents a rotation).
    CH_NODISCARD constexpr inline ChVector3<Real> GetAxisZ() const;

    // QUATERNION NORMS

    /// Compute the euclidean norm of the quaternion, that is its length or magnitude.
    CH_NODISCARD constexpr inline Real Length() const;

    /// Compute the squared euclidean norm of the quaternion.
    CH_NODISCARD constexpr inline Real Length2() const;

    /// Compute the infinity norm of the quaternion, that is the maximum absolute value of one of its elements.
    CH_NODISCARD constexpr inline Real LengthInf() const;

    // OPERATORS OVERLOADING

    /// Subscript operator.
    CH_NODISCARD constexpr inline Real& operator[](const unsigned& index);
    CH_NODISCARD constexpr inline const Real& operator[](const unsigned& index) const;

    /// Assignment operator: copy from another quaternion.
    constexpr inline ChQuaternion<Real>& operator=(const ChQuaternion<Real>& other);

    /// Operator for sign change.
    CH_NODISCARD constexpr inline ChQuaternion<Real> operator+() const;
    CH_NODISCARD constexpr inline ChQuaternion<Real> operator-() const;

    /// Operator for making a conjugate quaternion (the original is not changed).
    /// A conjugate quaternion has the vectorial part with changed sign.
    CH_NODISCARD constexpr inline ChQuaternion<Real> operator!() const;

    /// Operator for quaternion sum.
    CH_NODISCARD constexpr inline ChQuaternion<Real> operator+(const ChQuaternion<Real>& other) const;
    constexpr inline ChQuaternion<Real>& operator+=(const ChQuaternion<Real>& other);

    /// Operator for quaternion difference.
    CH_NODISCARD constexpr inline ChQuaternion<Real> operator-(const ChQuaternion<Real>& other) const;
    constexpr inline ChQuaternion<Real>& operator-=(const ChQuaternion<Real>& other);

    // NOTE
    // The following * and *= operators had a different behaviour prior to 13/9/2014,
    // but we assume no one used * *= in that previous form (element-by-element product).

    /// Operator for quaternion product: A*B means the typical quaternion product.
    /// Notes:
    /// - since unit quaternions can represent rotations, the product can represent a
    ///   concatenation of rotations as:
    ///        frame_rotation_2to0 = frame_rotation_1to0 * frame_rotation_2to1
    /// - pay attention to operator low precedence (see C++ precedence rules!)
    /// - quaternion product is not commutative.
    CH_NODISCARD constexpr inline ChQuaternion<Real> operator*(const ChQuaternion<Real>& other) const;

    /// Operator for quaternion product and assignment:
    /// A*=B means A'=A*B, with typical quaternion product.
    /// Notes:
    /// - since unit quaternions can represent rotations, the product can represent a
    ///   post-concatenation of a rotation in a kinematic chain.
    /// - quaternion product is not commutative.
    constexpr inline ChQuaternion<Real>& operator*=(const ChQuaternion<Real>& other);

    /// Operator for 'specular' quaternion product: A>>B = B*A.
    /// Notes:
    /// - since unit quaternions can represent rotations, the product can represent a
    ///   concatenation of rotations as:
    ///       frame_rotation_2to0 = frame_rotation_2to1 >> frame_rotation_1to0
    /// - pay attention to operator low precedence (see C++ precedence rules!)
    /// - quaternion product is not commutative.
    CH_NODISCARD constexpr inline ChQuaternion<Real> operator>>(const ChQuaternion<Real>& other) const;

    /// Operator for quaternion 'specular' product and assignment:
    /// A>>=B means A'=A>>B, or A'=B*A with typical quaternion product.
    /// Notes:
    /// - since unit quaternions can represent rotations, the product can represent a
    ///   pre-concatenation of a rotation in a kinematic chain.
    /// - quaternion product is not commutative.
    constexpr inline ChQuaternion<Real>& operator>>=(const ChQuaternion<Real>& other);

    // Operator for scaling the quaternion by a scalar value, as q*s.
    CH_NODISCARD constexpr inline ChQuaternion<Real> operator*(Real s) const;
    constexpr inline ChQuaternion<Real>& operator*=(Real s);

    /// Operator for element-wise division.
    /// Note that this is NOT the quaternion division operation.
    CH_NODISCARD constexpr inline ChQuaternion<Real> operator/(const ChQuaternion<Real>& other) const;
    constexpr inline ChQuaternion<Real>& operator/=(const ChQuaternion<Real>& other);

    /// Operator for scaling the quaternion by inverse of a scalar value, as q/s.
    CH_NODISCARD constexpr inline ChQuaternion<Real> operator/(Real s) const;
    constexpr inline ChQuaternion<Real>& operator/=(Real s);

    /// Operator for dot product: A^B means the scalar dot-product A*B.
    /// Note: pay attention to operator low precedence (see C++ precedence rules!)
    CH_NODISCARD constexpr inline Real operator^(const ChQuaternion<Real>& other) const;

    /// Component-wise comparison operators.
    CH_NODISCARD constexpr inline bool operator<=(const ChQuaternion<Real>& other) const noexcept;
    CH_NODISCARD constexpr inline bool operator>=(const ChQuaternion<Real>& other) const noexcept;
    CH_NODISCARD constexpr inline bool operator<(const ChQuaternion<Real>& other) const noexcept;
    CH_NODISCARD constexpr inline bool operator>(const ChQuaternion<Real>& other) const noexcept;
    CH_NODISCARD constexpr inline bool operator==(const ChQuaternion<Real>& other) const noexcept;
    CH_NODISCARD constexpr inline bool operator!=(const ChQuaternion<Real>& other) const noexcept;

    // FUNCTIONS

    /// Set this quaternion to the sum of A and B: this = A + B.
    constexpr inline void Add(const ChQuaternion<Real>& A, const ChQuaternion<Real>& B);

    /// Set this quaternion to the difference of A and B: this = A - B.
    constexpr inline void Sub(const ChQuaternion<Real>& A, const ChQuaternion<Real>& B);

    /// Set this quaternion to the quaternion product of the two quaternions A and B,
    /// following the classic Hamilton rule:  this = AxB.
    /// This is the true, typical quaternion product. It is NOT commutative.
    constexpr inline void Cross(const ChQuaternion<Real>& qa, const ChQuaternion<Real>& qb);

    /// Return the dot product with another quaternion: result = this ^ B.
    constexpr inline Real Dot(const ChQuaternion<Real>& B) const;

    /// Set this quaternion to the product of a quaternion A and scalar s: this = A * s.
    constexpr inline void Mul(const ChQuaternion<Real>& A, Real s);

    /// Scale this quaternion by a scalar: this *= s.
    constexpr inline void Scale(const Real& s);

    /// Normalize this quaternion in place, so that its euclidean length is 1.
    /// Return false if the original quaternion had zero length (in which case the quaternion
    /// is set to [1,0,0,0]) and return true otherwise.
    constexpr inline bool Normalize();

    /// Return a normalized copy of this quaternion, with euclidean length = 1.
    /// Not to be confused with Normalize() which normalizes in place.
    CH_NODISCARD constexpr inline ChQuaternion<Real> GetNormalized() const;

    /// Set this quaternion to the conjugate of the A quaternion.
    constexpr inline void Conjugate(const ChQuaternion<Real>& A);

    /// Conjugate this quaternion in place (its vectorial part changes sign).
    constexpr inline void Conjugate();

    /// Return a conjugated version of this quaternion.
    CH_NODISCARD constexpr inline ChQuaternion<Real> GetConjugate() const;

    /// Return the inverse of this quaternion.
    CH_NODISCARD constexpr inline ChQuaternion<Real> GetInverse() const;

    // TRANSFORMATIONS

    /// Rotate the vector A on the basis of this quaternion: res=p*[0,A]*p'
    /// (speed-optimized version). Endomorphism assumes p is already normalized.
    CH_NODISCARD constexpr inline ChVector3<Real> Rotate(const ChVector3<Real>& A) const;

    /// Rotate the vector A on the basis of conjugate of this quaternion: res=p'*[0,A]*p
    /// (speed-optimized version).  Endomorphism assumes p is already normalized.
    CH_NODISCARD constexpr inline ChVector3<Real> RotateBack(const ChVector3<Real>& A) const;

    // CONVERSIONS

    /// Set the quaternion from a rotation vector (ie. a 3D axis of rotation with length as angle of rotation)
    /// defined in absolute coords.
    /// If you need distinct axis and angle, use SetFromAngleAxis().
    constexpr inline void SetFromRotVec(const ChVector3<Real>& angle_axis);

    /// Get the rotation vector (ie. a 3D axis of rotation with length as angle of rotation) from a quaternion.
    /// If you need distinct axis and angle, use GetAngleAxis().
    CH_NODISCARD constexpr inline ChVector3<Real> GetRotVec() const;

    /// Set the quaternion from an angle of rotation and an axis, defined in absolute coords.
    /// The axis is supposed to be fixed, i.e. it is constant during rotation!
    /// NOTE, axis must be normalized!
    /// If you need directly the rotation vector=axis * angle, use SetFromRotVec().
    constexpr inline void SetFromAngleAxis(const Real& angle, const ChVector3<Real>& axis);

    /// Set the quaternion from an angle of rotation about X axis.
    constexpr inline void SetFromAngleX(const Real& angleX) { SetFromAngleAxis(angleX, ChVector3<Real>(1, 0, 0)); }

    /// Set the quaternion from an angle of rotation about Y axis.
    constexpr inline void SetFromAngleY(const Real& angleY) { SetFromAngleAxis(angleY, ChVector3<Real>(0, 1, 0)); }

    /// Set the quaternion from an angle of rotation about Z axis.
    constexpr inline void SetFromAngleZ(const Real& angleZ) { SetFromAngleAxis(angleZ, ChVector3<Real>(0, 0, 1)); }

    /// Convert the quaternion to an angle of rotation and an axis, defined in absolute coords.
    /// Resulting angle and axis must be passed as parameters.
    /// Note that angle is in [-PI....+PI] range. Also remember  (angle, axis) is the same of (-angle,-axis).
    /// If you need directly the rotation vector=axis * angle, use GetRotVec().
    constexpr inline void GetAngleAxis(Real& angle, ChVector3<Real>& axis) const;

    /// Set the quaternion from Cardan angles ZYX (Tait-Bryan sequence Z-Y'-X'', intrinsic rotations).
    /// NASA angles (heading, bank, attitude).
    constexpr inline void SetFromCardanAnglesZYX(const ChVector3<Real>& ang);

    /// Convert the quaternion to Cardan angles ZYX (Tait-Bryan sequence Z-Y'-X'', intrinsic rotations).
    /// NASA angles (heading, bank, attitude).
    CH_NODISCARD constexpr inline ChVector3<Real> GetCardanAnglesZYX() const;

    /// Set the quaternion from Cardan angles XYZ (Tait-Bryan sequence X-Y'-Z'', intrinsic rotations).
    constexpr inline void SetFromCardanAnglesXYZ(const ChVector3<Real>& ang);

    /// Convert the quaternion to Cardan angles XYZ (Tait-Bryan sequence X-Y'-Z'', intrinsic rotations).
    CH_NODISCARD constexpr inline ChVector3<Real> GetCardanAnglesXYZ() const;

    /// Set the quaternion dq/dt. Inputs: the vector of angular speed w specified in absolute coords,
    /// and the rotation expressed as a quaternion q.
    constexpr inline void SetDtFromAngVelAbs(const ChVector3<Real>& w, const ChQuaternion<Real>& q);

    /// Set the quaternion dq/dt. Inputs: the vector of angular speed w specified in relative coords,
    /// and the rotation expressed as a quaternion q.
    constexpr inline void SetDtFromAngVelRel(const ChVector3<Real>& w, const ChQuaternion<Real>& q);

    /// Compute the vector of angular speed 'w' specified in absolute coords,
    /// from this quaternion dq/dt and the rotation expressed as a quaternion q.
    constexpr inline void GetAngVelAbs(ChVector3<Real>& w, const ChQuaternion<Real>& q);

    /// Compute the vector of angular speed 'w' specified in relative coords,
    /// from this quaternion dq/dt and the rotation expressed as a quaternion q.
    constexpr inline void GetAngVelRel(ChVector3<Real>& w, const ChQuaternion<Real>& q);

    /// Set the quaternion ddq/dtdt. Inputs: the vector of angular acceleration 'a' specified
    /// in absolute coords, the rotation expressed as a quaternion q, the rotation speed
    /// as a quaternion 'q_dt'.
    constexpr inline void SetDt2FromAngAccAbs(const ChVector3<Real>& a, const ChQuaternion<Real>& q, const ChQuaternion<Real>& q_dt);

    /// Set the quaternion ddq/dtdt. Inputs: the vector of angular acceleration 'a' specified
    /// in relative coords, the rotation expressed as a quaternion q, the rotation speed
    /// as a quaternion 'q_dt'.
    constexpr inline void SetDt2FromAngAccRel(const ChVector3<Real>& a, const ChQuaternion<Real>& q, const ChQuaternion<Real>& q_dt);

    /// Set the quaternion dq/dt. Inputs:  the axis of rotation 'axis' (assuming it is already normalized
    /// and expressed in absolute coords), the angular speed 'angle_dt' (scalar value), and the
    /// rotation expressed as a quaternion 'q'.
    constexpr inline void SetDtFromAngleAxis(const ChQuaternion<Real>& q, Real angle_dt, const ChVector3<Real>& axis);

    /// Set the quaternion ddq/dtdt. Inputs: the axis of ang. acceleration 'axis' (assuming it is already
    /// normalized and expressed in absolute coords), the angular acceleration 'angle_dtdt' (scalar value),
    /// the rotation expressed as a quaternion 'quat' and th rotation speed 'q_dt'.
    constexpr inline void SetDt2FromAngleAxis(const ChQuaternion<Real>& q,
                             const ChQuaternion<Real>& q_dt,
                             const Real& angle_dtdt,
                             const ChVector3<Real>& axis);

    // SERIALIZATION

    /// Method to allow serialization of transient m_data to archives.
    constexpr inline void ArchiveOut(ChArchiveOut& archive_out);

    /// Method to allow de-serialization of transient m_data from archives.
    constexpr inline void ArchiveIn(ChArchiveIn& archive_in);

  private:
    ARRAY_TYPE m_data;  ///< quaternion data (e0, e1, e2, e3)

    template <typename RealB>
    friend class ChQuaternion;
};

CH_CLASS_VERSION(ChQuaternion<double>, 0)

// -----------------------------------------------------------------------------

/// Alias for double-precision quaternions.
/// <pre>
/// Instead of writing
///    ChQuaternion<double> q;
/// or
///    ChQuaternion<> q;
/// you can use:
///    ChQuaterniond q;
/// </pre>
typedef ChQuaternion<double> ChQuaterniond;

/// Alias for single-precision quaternions.
/// <pre>
/// Instead of writing
///    ChQuaternion<float> q;
/// you can use:
///    ChQuaternionf q;
/// </pre>
typedef ChQuaternion<float> ChQuaternionf;

// -----------------------------------------------------------------------------
// STATIC QUATERNION MATH OPERATIONS
//
// These functions are here for people which prefer to use static functions
// instead of ChQuaternion class' member functions.
// NOTE: sometimes a wise adoption of the following functions may give faster
// results than using overloaded operators +/-/* in the quaternion class.

/// Check if two quaternions are equal.
ChApi static constexpr inline bool Qequal(const ChQuaterniond& qa, const ChQuaterniond& qb);

/// Check if quaternion is not null.
ChApi static constexpr inline bool Qnotnull(const ChQuaterniond& qa);

ChApi static constexpr inline double Qlength(const ChQuaterniond& q);

ChApi static constexpr inline ChQuaterniond Qadd(const ChQuaterniond& qa, const ChQuaterniond& qb);

ChApi static constexpr inline ChQuaterniond Qsub(const ChQuaterniond& qa, const ChQuaterniond& qb);

ChApi static constexpr inline ChQuaterniond Qscale(const ChQuaterniond& q, const double& fact);

/// Return the norm two of the quaternion. Euler's parameters have norm = 1
ChApi static constexpr inline ChQuaterniond Qnorm(const ChQuaterniond& q);

/// Return the conjugate of the quaternion [s,v1,v2,v3] is [s,-v1,-v2,-v3]
ChApi static constexpr inline ChQuaterniond Qconjugate(const ChQuaterniond& q);

/// Return the product of two quaternions. It is non-commutative (like cross product in vectors).
ChApi static constexpr inline ChQuaterniond Qcross(const ChQuaterniond& qa, const ChQuaterniond& qb);

/// Get the quaternion time derivative from the vector of angular speed, with w specified in _local_ coords.
ChApi static constexpr inline ChQuaterniond QuatDtFromAngVelRel(const ChVector3d& w, const ChQuaterniond& q);

/// Get the quaternion time derivative from the vector of angular speed, with w specified in _absolute_ coords.
ChApi static constexpr inline ChQuaterniond QuatDtFromAngVelAbs(const ChVector3d& w, const ChQuaterniond& q);

/// Get the quaternion first derivative from the vector of angular acceleration with a specified in _absolute_ coords.
ChApi static constexpr inline ChQuaterniond QuatDt2FromAngAccAbs(const ChVector3d& a, const ChQuaterniond& q, const ChQuaterniond& q_dt);

///	Get the quaternion second derivative from the vector of angular acceleration with a specified in _relative_ coords.
ChApi static constexpr inline ChQuaterniond QuatDt2FromAngAccRel(const ChVector3d& a, const ChQuaterniond& q, const ChQuaterniond& q_dt);

/// Given the imaginary (vectorial) {e1 e2 e3} part of a quaternion,
/// find the entire quaternion q = {e0, e1, e2, e3}.
/// Note: singularities are possible.
ChApi static constexpr inline ChQuaterniond QuatFromImaginary(const ChVector3d& imm);

/// Given the imaginary (vectorial) {e1 e2 e3} part of a quaternion time derivative,
/// find the entire quaternion q = {e0, e1, e2, e3}.
/// Note: singularities are possible.
ChApi static constexpr inline ChQuaterniond QuatDtFromImaginary(const ChVector3d& im_dt, const ChQuaterniond& q);

/// Given the imaginary (vectorial) {e1 e2 e3} part of a quaternion second time derivative,
/// find the entire quaternion q = {e0, e1, e2, e3}.
/// Note: singularities are possible.
ChApi static constexpr inline ChQuaterniond QuatDt2FromImaginary(const ChVector3d& im_dtdt, const ChQuaterniond& q, const ChQuaterniond& q_dt);

ChApi static constexpr inline ChVector3d VaxisXfromQuat(const ChQuaterniond& quat);

/// Insertion of quaternion to output stream.
template <typename Real>
inline std::ostream& operator<<(std::ostream& out, const ChQuaternion<Real>& q) {
    out << q.e0() << "  " << q.e1() << "  " << q.e2() << "  " << q.e3();
    return out;
}

// =============================================================================
// IMPLEMENTATION OF ChQuaternion<Real> methods
// =============================================================================

// -----------------------------------------------------------------------------
// Constructors

template <class Real>
constexpr inline ChQuaternion<Real>::ChQuaternion() : m_data({0, 0, 0, 0}) {}

template <class Real>
constexpr inline ChQuaternion<Real>::ChQuaternion(const Real& e0, const Real& e1, const Real& e2, const Real& e3) : m_data({e0, e1, e2, e3}) {}

template <class Real>
constexpr inline ChQuaternion<Real>::ChQuaternion(const Real& s, const ChVector3<Real>& v) : m_data({s, v.x(), v.y(), v.z()}) {}

template <class Real>
constexpr inline ChQuaternion<Real>::ChQuaternion(const ChQuaternion<Real>& other) : m_data(other.m_data) {}

template <class Real>
template <class RealB>
constexpr inline ChQuaternion<Real>::ChQuaternion(const ChQuaternion<RealB>& other) : m_data({static_cast<Real>(other.m_data[0]),
                                                                                                   static_cast<Real>(other.m_data[1]),
                                                                                                   static_cast<Real>(other.m_data[2]),
                                                                                                   static_cast<Real>(other.m_data[3])}) {} 

// -----------------------------------------------------------------------------
// Subscript operators

template <class Real>
constexpr inline Real& ChQuaternion<Real>::operator[](const unsigned& index) {
    assert(index < 4);
    return m_data[index];
}

template <class Real>
constexpr inline const Real& ChQuaternion<Real>::operator[](const unsigned& index) const {
    assert(index < 4);
    return m_data[index];
}

// -----------------------------------------------------------------------------
// Assignments

template <class Real>
constexpr inline ChQuaternion<Real>& ChQuaternion<Real>::operator=(const ChQuaternion<Real>& other) {
    if (&other == this)
        return *this;
    m_data[0] = other.m_data[0];
    m_data[1] = other.m_data[1];
    m_data[2] = other.m_data[2];
    m_data[3] = other.m_data[3];
    return *this;
}

// -----------------------------------------------------------------------------
// Sign operators

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::operator+() const {
    return *this;
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::operator-() const {
    return ChQuaternion<Real>(-m_data[0], -m_data[1], -m_data[2], -m_data[3]);
}

// -----------------------------------------------------------------------------
// Arithmetic & quaternion operations

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::operator!() const {
    return ChQuaternion<Real>(m_data[0], -m_data[1], -m_data[2], -m_data[3]);
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::operator+(const ChQuaternion<Real>& other) const {
    return ChQuaternion<Real>(m_data[0] + other.m_data[0], m_data[1] + other.m_data[1], m_data[2] + other.m_data[2],
                              m_data[3] + other.m_data[3]);
}

template <class Real>
constexpr inline ChQuaternion<Real>& ChQuaternion<Real>::operator+=(const ChQuaternion<Real>& other) {
    m_data[0] += other.m_data[0];
    m_data[1] += other.m_data[1];
    m_data[2] += other.m_data[2];
    m_data[3] += other.m_data[3];
    return *this;
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::operator-(const ChQuaternion<Real>& other) const {
    return ChQuaternion<Real>(m_data[0] - other.m_data[0], m_data[1] - other.m_data[1], m_data[2] - other.m_data[2],
                              m_data[3] - other.m_data[3]);
}

template <class Real>
constexpr inline ChQuaternion<Real>& ChQuaternion<Real>::operator-=(const ChQuaternion<Real>& other) {
    m_data[0] -= other.m_data[0];
    m_data[1] -= other.m_data[1];
    m_data[2] -= other.m_data[2];
    m_data[3] -= other.m_data[3];
    return *this;
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::operator*(const ChQuaternion<Real>& other) const {
    ChQuaternion<Real> q;
    q.Cross(*this, other);
    return q;
}
template <class Real>
constexpr inline ChQuaternion<Real>& ChQuaternion<Real>::operator*=(const ChQuaternion<Real>& other) {
    this->Cross(*this, other);
    return *this;
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::operator>>(const ChQuaternion<Real>& other) const {
    ChQuaternion<Real> q;
    q.Cross(other, *this);
    return q;
}

template <class Real>
constexpr inline ChQuaternion<Real>& ChQuaternion<Real>::operator>>=(const ChQuaternion<Real>& other) {
    this->Cross(other, *this);
    return *this;
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::operator*(Real s) const {
    return ChQuaternion<Real>(m_data[0] * s, m_data[1] * s, m_data[2] * s, m_data[3] * s);
}

template <class Real>
constexpr inline ChQuaternion<Real>& ChQuaternion<Real>::operator*=(Real s) {
    m_data[0] *= s;
    m_data[1] *= s;
    m_data[2] *= s;
    m_data[3] *= s;
    return *this;
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::operator/(const ChQuaternion<Real>& other) const {
    return ChQuaternion<Real>(m_data[0] / other.m_data[0], m_data[1] / other.m_data[1], m_data[2] / other.m_data[2],
                              m_data[3] / other.m_data[3]);
}

template <class Real>
constexpr inline ChQuaternion<Real>& ChQuaternion<Real>::operator/=(const ChQuaternion<Real>& other) {
    m_data[0] /= other.m_data[0];
    m_data[1] /= other.m_data[1];
    m_data[2] /= other.m_data[2];
    m_data[3] /= other.m_data[3];
    return *this;
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::operator/(Real s) const {
    Real oos = 1 / s;
    return ChQuaternion<Real>(m_data[0] * oos, m_data[1] * oos, m_data[2] * oos, m_data[3] * oos);
}

template <class Real>
constexpr inline ChQuaternion<Real>& ChQuaternion<Real>::operator/=(Real s) {
    Real oos = 1 / s;
    m_data[0] *= oos;
    m_data[1] *= oos;
    m_data[2] *= oos;
    m_data[3] *= oos;
    return *this;
}

template <class Real>
constexpr inline Real ChQuaternion<Real>::operator^(const ChQuaternion<Real>& other) const {
    return this->Dot(other);
}

// -----------------------------------------------------------------------------
// Comparison operations

template <class Real>
constexpr inline bool ChQuaternion<Real>::operator<=(const ChQuaternion<Real>& other) const noexcept {
    return m_data[0] <= other.m_data[0] && m_data[1] <= other.m_data[1] && m_data[2] <= other.m_data[2] &&
           m_data[3] <= other.m_data[3];
}

template <class Real>
constexpr inline bool ChQuaternion<Real>::operator>=(const ChQuaternion<Real>& other) const noexcept {
    return m_data[0] >= other.m_data[0] && m_data[1] >= other.m_data[1] && m_data[2] >= other.m_data[2] &&
           m_data[3] >= other.m_data[3];
}

template <class Real>
constexpr inline bool ChQuaternion<Real>::operator<(const ChQuaternion<Real>& other) const noexcept {
    return m_data[0] < other.m_data[0] && m_data[1] < other.m_data[1] && m_data[2] < other.m_data[2] &&
           m_data[3] < other.m_data[3];
}

template <class Real>
constexpr inline bool ChQuaternion<Real>::operator>(const ChQuaternion<Real>& other) const noexcept {
    return m_data[0] > other.m_data[0] && m_data[1] > other.m_data[1] && m_data[2] > other.m_data[2] &&
           m_data[3] > other.m_data[3];
}

template <class Real>
constexpr inline bool ChQuaternion<Real>::operator==(const ChQuaternion<Real>& other) const noexcept {
    return other.m_data[0] == m_data[0] && other.m_data[1] == m_data[1] && other.m_data[2] == m_data[2] &&
           other.m_data[3] == m_data[3];
}

template <class Real>
constexpr inline bool ChQuaternion<Real>::operator!=(const ChQuaternion<Real>& other) const noexcept {
    return !(*this == other);
}

// -----------------------------------------------------------------------------
// Functions

template <class Real>
inline void ChQuaternion<Real>::Set(const Real& e0, const Real& e1, const Real& e2, const Real& e3) {
    m_data[0] = e0;
    m_data[1] = e1;
    m_data[2] = e2;
    m_data[3] = e3;
}

template <class Real>
inline void ChQuaternion<Real>::Set(const ChQuaternion<Real>& q) {
    m_data[0] = q.m_data[0];
    m_data[1] = q.m_data[1];
    m_data[2] = q.m_data[2];
    m_data[3] = q.m_data[3];
}

template <class Real>
inline void ChQuaternion<Real>::Set(const Real& s) {
    m_data[0] = s;
    m_data[1] = s;
    m_data[2] = s;
    m_data[3] = s;
}

template <class Real>
inline void ChQuaternion<Real>::SetNull() {
    m_data[0] = 0;
    m_data[1] = 0;
    m_data[2] = 0;
    m_data[3] = 0;
}

template <class Real>
inline void ChQuaternion<Real>::SetUnit() {
    m_data[0] = 1;
    m_data[1] = 0;
    m_data[2] = 0;
    m_data[3] = 0;
}

template <class Real>
inline void ChQuaternion<Real>::SetScalar(const Real& s) {
    m_data[0] = s;
}

template <class Real>
inline void ChQuaternion<Real>::SetVector(const ChVector3<Real>& v) {
    m_data[1] = v.x();
    m_data[2] = v.y();
    m_data[3] = v.z();
}

template <class Real>
constexpr inline bool ChQuaternion<Real>::IsNull() const {
    return m_data[0] == 0 && m_data[1] == 0 && m_data[2] == 0 && m_data[3] == 0;
}

template <class Real>
constexpr inline bool ChQuaternion<Real>::IsIdentity() const {
    return m_data[0] == 1 && m_data[1] == 0 && m_data[2] == 0 && m_data[3] == 0;
}

template <class Real>
constexpr inline bool ChQuaternion<Real>::Equals(const ChQuaternion<Real>& other) const {
    return (other.m_data[0] == m_data[0]) && (other.m_data[1] == m_data[1]) && (other.m_data[2] == m_data[2]) &&
           (other.m_data[3] == m_data[3]);
}

template <class Real>
constexpr inline bool ChQuaternion<Real>::Equals(const ChQuaternion<Real>& other, const Real& tol) const {
    return (std::abs(other.m_data[0] - m_data[0]) < tol) && (std::abs(other.m_data[1] - m_data[1]) < tol) &&
           (std::abs(other.m_data[2] - m_data[2]) < tol) && (std::abs(other.m_data[3] - m_data[3]) < tol);
}

template <class Real>
constexpr inline ChVector3<Real> ChQuaternion<Real>::GetVector() const {
    return ChVector3<Real>(m_data[1], m_data[2], m_data[3]);
}

template <class Real>
constexpr inline ChVector3<Real> ChQuaternion<Real>::GetAxisX() const {
    return ChVector3<Real>((m_data[0] * m_data[0] + m_data[1] * m_data[1]) * 2 - 1,
                           (m_data[1] * m_data[2] + m_data[0] * m_data[3]) * 2,
                           (m_data[1] * m_data[3] - m_data[0] * m_data[2]) * 2);
}

template <class Real>
constexpr inline ChVector3<Real> ChQuaternion<Real>::GetAxisY() const {
    return ChVector3<Real>((m_data[1] * m_data[2] - m_data[0] * m_data[3]) * 2,
                           (m_data[0] * m_data[0] + m_data[2] * m_data[2]) * 2 - 1,
                           (m_data[2] * m_data[3] + m_data[0] * m_data[1]) * 2);
}

template <class Real>
constexpr inline ChVector3<Real> ChQuaternion<Real>::GetAxisZ() const {
    return ChVector3<Real>((m_data[1] * m_data[3] + m_data[0] * m_data[2]) * 2,
                           (m_data[2] * m_data[3] - m_data[0] * m_data[1]) * 2,
                           (m_data[0] * m_data[0] + m_data[3] * m_data[3]) * 2 - 1);
}

template <class Real>
constexpr inline Real ChQuaternion<Real>::Length() const {
    return std::sqrt(Length2());
}

template <class Real>
constexpr inline Real ChQuaternion<Real>::Length2() const {
    return this->Dot(*this);
}

template <class Real>
constexpr inline Real ChQuaternion<Real>::LengthInf() const {
    Real e0e1 = std::max(std::abs(m_data[0]), std::abs(m_data[1]));
    Real e0e1e2 = std::max(e0e1, std::abs(m_data[2]));
    return std::max(e0e1e2, std::abs(m_data[3]));
}

template <class Real>
constexpr inline void ChQuaternion<Real>::Add(const ChQuaternion<Real>& A, const ChQuaternion<Real>& B) {
    m_data[0] = A.m_data[0] + B.m_data[0];
    m_data[1] = A.m_data[1] + B.m_data[1];
    m_data[2] = A.m_data[2] + B.m_data[2];
    m_data[3] = A.m_data[3] + B.m_data[3];
}

template <class Real>
constexpr inline void ChQuaternion<Real>::Sub(const ChQuaternion<Real>& A, const ChQuaternion<Real>& B) {
    m_data[0] = A.m_data[0] - B.m_data[0];
    m_data[1] = A.m_data[1] - B.m_data[1];
    m_data[2] = A.m_data[2] - B.m_data[2];
    m_data[3] = A.m_data[3] - B.m_data[3];
}

template <class Real>
constexpr inline void ChQuaternion<Real>::Cross(const ChQuaternion<Real>& qa, const ChQuaternion<Real>& qb) {
    Real w = qa.m_data[0] * qb.m_data[0] - qa.m_data[1] * qb.m_data[1] - qa.m_data[2] * qb.m_data[2] -
             qa.m_data[3] * qb.m_data[3];
    Real x = qa.m_data[0] * qb.m_data[1] + qa.m_data[1] * qb.m_data[0] - qa.m_data[3] * qb.m_data[2] +
             qa.m_data[2] * qb.m_data[3];
    Real y = qa.m_data[0] * qb.m_data[2] + qa.m_data[2] * qb.m_data[0] + qa.m_data[3] * qb.m_data[1] -
             qa.m_data[1] * qb.m_data[3];
    Real z = qa.m_data[0] * qb.m_data[3] + qa.m_data[3] * qb.m_data[0] - qa.m_data[2] * qb.m_data[1] +
             qa.m_data[1] * qb.m_data[2];
    m_data[0] = w;
    m_data[1] = x;
    m_data[2] = y;
    m_data[3] = z;
}

template <class Real>
constexpr inline Real ChQuaternion<Real>::Dot(const ChQuaternion<Real>& B) const {
    return (m_data[0] * B.m_data[0]) + (m_data[1] * B.m_data[1]) + (m_data[2] * B.m_data[2]) +
           (m_data[3] * B.m_data[3]);
}

template <class Real>
constexpr inline void ChQuaternion<Real>::Mul(const ChQuaternion<Real>& A, Real s) {
    m_data[0] = A.m_data[0] * s;
    m_data[1] = A.m_data[1] * s;
    m_data[2] = A.m_data[2] * s;
    m_data[3] = A.m_data[3] * s;
}

template <class Real>
constexpr inline void ChQuaternion<Real>::Scale(const Real& s) {
    m_data[0] *= s;
    m_data[1] *= s;
    m_data[2] *= s;
    m_data[3] *= s;
}

template <class Real>
constexpr inline bool ChQuaternion<Real>::Normalize() {
    Real length = this->Length();
    if (length < std::numeric_limits<Real>::min()) {
        m_data[0] = 1;
        m_data[1] = 0;
        m_data[2] = 0;
        m_data[3] = 0;
        return false;
    }
    this->Scale(1 / length);
    return true;
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::GetNormalized() const {
    ChQuaternion<Real> q(*this);
    q.Normalize();
    return q;
}

template <class Real>
constexpr inline void ChQuaternion<Real>::Conjugate(const ChQuaternion<Real>& A) {
    m_data[0] = +A.m_data[0];
    m_data[1] = -A.m_data[1];
    m_data[2] = -A.m_data[2];
    m_data[3] = -A.m_data[3];
}

template <class Real>
constexpr inline void ChQuaternion<Real>::Conjugate() {
    Conjugate(*this);
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::GetConjugate() const {
    return ChQuaternion<Real>(m_data[0], -m_data[1], -m_data[2], -m_data[3]);
}

template <class Real>
constexpr inline ChQuaternion<Real> ChQuaternion<Real>::GetInverse() const {
    ChQuaternion<Real> invq = this->GetConjugate();
    invq.Scale(1 / this->Length2());
    return invq;
}

template <class Real>
constexpr inline ChVector3<Real> ChQuaternion<Real>::Rotate(const ChVector3<Real>& A) const {
    Real e0e0 = m_data[0] * m_data[0];
    Real e1e1 = m_data[1] * m_data[1];
    Real e2e2 = m_data[2] * m_data[2];
    Real e3e3 = m_data[3] * m_data[3];
    Real e0e1 = m_data[0] * m_data[1];
    Real e0e2 = m_data[0] * m_data[2];
    Real e0e3 = m_data[0] * m_data[3];
    Real e1e2 = m_data[1] * m_data[2];
    Real e1e3 = m_data[1] * m_data[3];
    Real e2e3 = m_data[2] * m_data[3];
    return ChVector3<Real>(((e0e0 + e1e1) * 2 - 1) * A.x() + ((e1e2 - e0e3) * 2) * A.y() + ((e1e3 + e0e2) * 2) * A.z(),
                           ((e1e2 + e0e3) * 2) * A.x() + ((e0e0 + e2e2) * 2 - 1) * A.y() + ((e2e3 - e0e1) * 2) * A.z(),
                           ((e1e3 - e0e2) * 2) * A.x() + ((e2e3 + e0e1) * 2) * A.y() + ((e0e0 + e3e3) * 2 - 1) * A.z());
}

template <class Real>
constexpr inline ChVector3<Real> ChQuaternion<Real>::RotateBack(const ChVector3<Real>& A) const {
    Real e0e0 = +m_data[0] * m_data[0];
    Real e1e1 = +m_data[1] * m_data[1];
    Real e2e2 = +m_data[2] * m_data[2];
    Real e3e3 = +m_data[3] * m_data[3];
    Real e0e1 = -m_data[0] * m_data[1];
    Real e0e2 = -m_data[0] * m_data[2];
    Real e0e3 = -m_data[0] * m_data[3];
    Real e1e2 = +m_data[1] * m_data[2];
    Real e1e3 = +m_data[1] * m_data[3];
    Real e2e3 = +m_data[2] * m_data[3];
    return ChVector3<Real>(((e0e0 + e1e1) * 2 - 1) * A.x() + ((e1e2 - e0e3) * 2) * A.y() + ((e1e3 + e0e2) * 2) * A.z(),
                           ((e1e2 + e0e3) * 2) * A.x() + ((e0e0 + e2e2) * 2 - 1) * A.y() + ((e2e3 - e0e1) * 2) * A.z(),
                           ((e1e3 - e0e2) * 2) * A.x() + ((e2e3 + e0e1) * 2) * A.y() + ((e0e0 + e3e3) * 2 - 1) * A.z());
}

template <class Real>
constexpr inline void ChQuaternion<Real>::SetFromRotVec(const ChVector3<Real>& angle_axis) {
    Real theta_squared = angle_axis.Length2();
    // For non-zero rotation:
    if (theta_squared > 1e-30) {
        Real theta = std::sqrt(theta_squared);
        Real half_theta = theta / 2;
        Real k = std::sin(half_theta) / theta;
        m_data[0] = std::cos(half_theta);
        m_data[1] = angle_axis.x() * k;
        m_data[2] = angle_axis.y() * k;
        m_data[3] = angle_axis.z() * k;
    } else {
        // For almost zero rotation:
        Real k(0.5);
        m_data[0] = Real(1.0);
        m_data[1] = angle_axis.x() * k;
        m_data[2] = angle_axis.y() * k;
        m_data[3] = angle_axis.z() * k;
    }
}

template <class Real>
constexpr inline ChVector3<Real> ChQuaternion<Real>::GetRotVec() const {
    ChVector3<Real> angle_axis;
    Real sin_squared = m_data[1] * m_data[1] + m_data[2] * m_data[2] + m_data[3] * m_data[3];
    // For non-zero rotation
    if (sin_squared > 1e-30) {
        Real sin_theta = std::sqrt(sin_squared);
        Real k = 2 * std::atan2(sin_theta, m_data[0]) / sin_theta;
        angle_axis.x() = m_data[1] * k;
        angle_axis.y() = m_data[2] * k;
        angle_axis.z() = m_data[3] * k;
    } else {
        // For almost zero rotation
        Real k(2.0);
        angle_axis.x() = m_data[1] * k;
        angle_axis.y() = m_data[2] * k;
        angle_axis.z() = m_data[3] * k;
    }
    return angle_axis;
}

template <class Real>
constexpr inline void ChQuaternion<Real>::SetFromAngleAxis(const Real& angle, const ChVector3<Real>& axis) {
    Real halfang = (angle / 2);
    Real sinhalf = std::sin(halfang);
    m_data[0] = std::cos(halfang);
    m_data[1] = axis.x() * sinhalf;
    m_data[2] = axis.y() * sinhalf;
    m_data[3] = axis.z() * sinhalf;
}

template <class Real>
constexpr inline void ChQuaternion<Real>::GetAngleAxis(Real& angle, ChVector3<Real>& axis) const {
    Real sin_squared = m_data[1] * m_data[1] + m_data[2] * m_data[2] + m_data[3] * m_data[3];
    // For non-zero rotation
    if (sin_squared > 0) {
        Real sin_theta = std::sqrt(sin_squared);
        angle = 2 * std::atan2(sin_theta, m_data[0]);
        Real k = 1 / sin_theta;
        axis.x() = m_data[1] * k;
        axis.y() = m_data[2] * k;
        axis.z() = m_data[3] * k;
        axis.Normalize();
    } else {
        // For almost zero rotation
        angle = 0.0;
        axis.x() = 1;  // m_data[1] * 2.0;
        axis.y() = 0;  // m_data[2] * 2.0;
        axis.z() = 0;  // m_data[3] * 2.0;
    }
    // Ensure that angle is always in  [-PI...PI] range
    auto PI = static_cast<Real>(CH_PI);
    if (angle > PI) {
        angle -= 2 * PI;
    } else if (angle < -PI) {
        angle += 2 * PI;
    }
}

template <class Real>
constexpr inline void ChQuaternion<Real>::SetFromCardanAnglesZYX(const ChVector3<Real>& ang) {
    Real c1 = std::cos(ang.z() / 2);
    Real s1 = std::sin(ang.z() / 2);
    Real c2 = std::cos(ang.x() / 2);
    Real s2 = std::sin(ang.x() / 2);
    Real c3 = std::cos(ang.y() / 2);
    Real s3 = std::sin(ang.y() / 2);

    Real c1c2 = c1 * c2;
    Real s1s2 = s1 * s2;

    m_data[0] = c1c2 * c3 + s1s2 * s3;
    m_data[1] = c1c2 * s3 - s1s2 * c3;
    m_data[2] = c1 * s2 * c3 + s1 * c2 * s3;
    m_data[3] = s1 * c2 * c3 - c1 * s2 * s3;
}

template <class Real>
constexpr inline ChVector3<Real> ChQuaternion<Real>::GetCardanAnglesZYX() const {
    ChVector3<Real> nasa;
    Real sqw = m_data[0] * m_data[0];
    Real sqx = m_data[1] * m_data[1];
    Real sqy = m_data[2] * m_data[2];
    Real sqz = m_data[3] * m_data[3];
    // heading
    nasa.z() = std::atan2(2 * (m_data[1] * m_data[2] + m_data[3] * m_data[0]), (sqx - sqy - sqz + sqw));
    // bank
    nasa.y() = std::atan2(2 * (m_data[2] * m_data[3] + m_data[1] * m_data[0]), (-sqx - sqy + sqz + sqw));
    // attitude
    nasa.x() = std::asin(-2 * (m_data[1] * m_data[3] - m_data[2] * m_data[0]));
    return nasa;
}

template <class Real>
constexpr inline void ChQuaternion<Real>::SetFromCardanAnglesXYZ(const ChVector3<Real>& ang) {
    // Angles {phi;theta;psi} aka {roll;pitch;yaw}
    Real t0 = std::cos(ang.z() * Real(0.5));
    Real t1 = std::sin(ang.z() * Real(0.5));
    Real t2 = std::cos(ang.x() * Real(0.5));
    Real t3 = std::sin(ang.x() * Real(0.5));
    Real t4 = std::cos(ang.y() * Real(0.5));
    Real t5 = std::sin(ang.y() * Real(0.5));

    m_data[0] = t0 * t2 * t4 + t1 * t3 * t5;
    m_data[1] = t0 * t3 * t4 - t1 * t2 * t5;
    m_data[2] = t0 * t2 * t5 + t1 * t3 * t4;
    m_data[3] = t1 * t2 * t4 - t0 * t3 * t5;
}

template <class Real>
constexpr inline ChVector3<Real> ChQuaternion<Real>::GetCardanAnglesXYZ() const {
    // Angles {phi;theta;psi} aka {roll;pitch;yaw} rotation XYZ
    ChVector3<Real> euler;
    Real sq0 = m_data[0] * m_data[0];
    Real sq1 = m_data[1] * m_data[1];
    Real sq2 = m_data[2] * m_data[2];
    Real sq3 = m_data[3] * m_data[3];
    // roll
    euler.x() = std::atan2(2 * (m_data[2] * m_data[3] + m_data[0] * m_data[1]), sq3 - sq2 - sq1 + sq0);
    // pitch
    euler.y() = -std::asin(2 * (m_data[1] * m_data[3] - m_data[0] * m_data[2]));
    // yaw
    euler.z() = std::atan2(2 * (m_data[1] * m_data[2] + m_data[3] * m_data[0]), sq1 + sq0 - sq3 - sq2);
    return euler;
}

template <class Real>
constexpr inline void ChQuaternion<Real>::SetDtFromAngVelAbs(const ChVector3<Real>& w, const ChQuaternion<Real>& q) {
    ChQuaternion<Real> qwo(0, w);
    this->Cross(qwo, q);
    this->Scale((Real)0.5);  // {q_dt} = 1/2 {0,w}*{q}
}

template <class Real>
constexpr inline void ChQuaternion<Real>::SetDtFromAngVelRel(const ChVector3<Real>& w, const ChQuaternion<Real>& q) {
    ChQuaternion<Real> qwl(0, w);
    this->Cross(q, qwl);
    this->Scale((Real)0.5);  // {q_dt} = 1/2 {q}*{0,w_rel}
}

template <class Real>
constexpr inline void ChQuaternion<Real>::GetAngVelAbs(ChVector3<Real>& w, const ChQuaternion<Real>& q) {
    ChQuaternion<Real> qwo;
    qwo.Cross(*this, q.GetConjugate());
    w = qwo.GetVector();
    w.Scale(2);
}

template <class Real>
constexpr inline void ChQuaternion<Real>::GetAngVelRel(ChVector3<Real>& w, const ChQuaternion<Real>& q) {
    ChQuaternion<Real> qwl;
    qwl.Cross(q.GetConjugate(), *this);
    w = qwl.GetVector();
    w.Scale(2);
}

template <class Real>
constexpr inline void ChQuaternion<Real>::SetDt2FromAngAccAbs(const ChVector3<Real>& a,
                                                    const ChQuaternion<Real>& q,
                                                    const ChQuaternion<Real>& q_dt) {
    ChQuaternion<Real> qao(0, a);
    ChQuaternion<Real> qwo;
    ChQuaternion<Real> qtmpa;
    ChQuaternion<Real> qtmpb;
    qwo.Cross(q_dt, q.GetConjugate());
    qtmpb.Cross(qwo, q_dt);
    qtmpa.Cross(qao, q);
    qtmpa.Scale((Real)0.5);
    this->Add(qtmpa, qtmpb);
}

template <class Real>
constexpr inline void ChQuaternion<Real>::SetDt2FromAngAccRel(const ChVector3<Real>& a,
                                                    const ChQuaternion<Real>& q,
                                                    const ChQuaternion<Real>& q_dt) {
    ChQuaternion<Real> qal(0, a);
    ChQuaternion<Real> qwl;
    ChQuaternion<Real> qtmpa;
    ChQuaternion<Real> qtmpb;
    qwl.Cross(q.GetConjugate(), q_dt);
    qtmpb.Cross(q_dt, qwl);
    qtmpa.Cross(q, qal);
    qtmpa.Scale((Real)0.5);
    this->Add(qtmpa, qtmpb);
}

template <class Real>
constexpr inline void ChQuaternion<Real>::SetDtFromAngleAxis(const ChQuaternion<Real>& q,
                                                   Real angle_dt,
                                                   const ChVector3<Real>& axis) {
    this->SetDtFromAngVelAbs(angle_dt * axis, q);
}

template <class Real>
constexpr inline void ChQuaternion<Real>::SetDt2FromAngleAxis(const ChQuaternion<Real>& q,
                                                    const ChQuaternion<Real>& q_dt,
                                                    const Real& angle_dtdt,
                                                    const ChVector3<Real>& axis) {
    this->SetDt2FromAngAccAbs(angle_dtdt * axis, q, q_dt);
}

// -----------------------------------------------------------------------------
// Streaming operations

template <class Real>
constexpr inline void ChQuaternion<Real>::ArchiveOut(ChArchiveOut& archive_out) {
    // version number
    archive_out.VersionWrite<ChQuaterniond>();  // must use specialized template (any)
    // stream out all member m_data
    archive_out << CHNVP(m_data[0], "e0");
    archive_out << CHNVP(m_data[1], "e1");
    archive_out << CHNVP(m_data[2], "e2");
    archive_out << CHNVP(m_data[3], "e3");
}

template <class Real>
constexpr inline void ChQuaternion<Real>::ArchiveIn(ChArchiveIn& archive_in) {
    // version number
    /*int version =*/archive_in.VersionRead<ChQuaterniond>();  // must use specialized template (any)
    // stream in all member m_data
    archive_in >> CHNVP(m_data[0], "e0");
    archive_in >> CHNVP(m_data[1], "e1");
    archive_in >> CHNVP(m_data[2], "e2");
    archive_in >> CHNVP(m_data[3], "e3");
}

// -----------------------------------------------------------------------------
// CONSTANTS

/// Constant null quaternion: {0, 0, 0, 0}
ChApi static constexpr ChQuaterniond QNULL(0., 0., 0., 0.);


/// Constant unit quaternion: {1, 0, 0, 0} ,
/// corresponds to no rotation (diagonal rotation matrix)
ChApi static constexpr ChQuaterniond QUNIT(1., 0., 0., 0.);

// Constants for rotations of 90 degrees:
ChApi static constexpr ChQuaterniond Q_ROTATE_Y_TO_X(1 / CH_SQRT_2, 0, 0, -1 / CH_SQRT_2);
ChApi static constexpr ChQuaterniond Q_ROTATE_Y_TO_Z(1 / CH_SQRT_2, 1 / CH_SQRT_2, 0, 0);
ChApi static constexpr ChQuaterniond Q_ROTATE_X_TO_Y(1 / CH_SQRT_2, 0, 0, 1 / CH_SQRT_2);
ChApi static constexpr ChQuaterniond Q_ROTATE_X_TO_Z(1 / CH_SQRT_2, 0, -1 / CH_SQRT_2, 0);
ChApi static constexpr ChQuaterniond Q_ROTATE_Z_TO_Y(1 / CH_SQRT_2, -1 / CH_SQRT_2, 0, 0);
ChApi static constexpr ChQuaterniond Q_ROTATE_Z_TO_X(1 / CH_SQRT_2, 0, 1 / CH_SQRT_2, 0);

// Constants for rotations of 180 degrees:
ChApi static constexpr ChQuaterniond Q_FLIP_AROUND_X(0., 1., 0., 0.);
ChApi static constexpr ChQuaterniond Q_FLIP_AROUND_Y(0., 0., 1., 0.);
ChApi static constexpr ChQuaterniond Q_FLIP_AROUND_Z(0., 0., 0., 1.);

// -----------------------------------------------------------------------------

// Check if two quaternions are equal
CH_NODISCARD static constexpr inline bool Qequal(const ChQuaterniond& qa, const ChQuaterniond& qb) {
    return qa == qb;
}

// Check if quaternion is not null
CH_NODISCARD static constexpr inline bool Qnotnull(const ChQuaterniond& qa) {
    return (qa.e0() != 0) || (qa.e1() != 0) || (qa.e2() != 0) || (qa.e3() != 0);
}

CH_NODISCARD static constexpr inline double Qlength(const ChQuaterniond& q) {
    return (std::sqrt(std::pow(q.e0(), 2) + std::pow(q.e1(), 2) + std::pow(q.e2(), 2) + std::pow(q.e3(), 2)));
}

CH_NODISCARD static constexpr inline ChQuaterniond Qscale(const ChQuaterniond& q, const double& fact) {
    ChQuaterniond result;
    result.e0() = q.e0() * fact;
    result.e1() = q.e1() * fact;
    result.e2() = q.e2() * fact;
    result.e3() = q.e3() * fact;
    return result;
}

CH_NODISCARD static constexpr inline ChQuaterniond Qadd(const ChQuaterniond& qa, const ChQuaterniond& qb) {
    ChQuaterniond result;
    result.e0() = qa.e0() + qb.e0();
    result.e1() = qa.e1() + qb.e1();
    result.e2() = qa.e2() + qb.e2();
    result.e3() = qa.e3() + qb.e3();
    return result;
}

CH_NODISCARD static constexpr inline ChQuaterniond Qsub(const ChQuaterniond& qa, const ChQuaterniond& qb) {
    ChQuaterniond result;
    result.e0() = qa.e0() - qb.e0();
    result.e1() = qa.e1() - qb.e1();
    result.e2() = qa.e2() - qb.e2();
    result.e3() = qa.e3() - qb.e3();
    return result;
}

// Return the norm two of the quaternion. Euler's parameters have norm = 1
CH_NODISCARD static constexpr inline ChQuaterniond Qnorm(const ChQuaterniond& q) {
    const auto invlength = 1 / (Qlength(q));
    return Qscale(q, invlength);
}

// Return the conjugate of the quaternion [s,v1,v2,v3] is [s,-v1,-v2,-v3]
CH_NODISCARD static constexpr inline ChQuaterniond Qconjugate(const ChQuaterniond& q) {
    ChQuaterniond res;
    res.e0() = q.e0();
    res.e1() = -q.e1();
    res.e2() = -q.e2();
    res.e3() = -q.e3();
    return (res);
}

// Return the product of two quaternions. It is non-commutative (like cross product in vectors).
CH_NODISCARD static constexpr inline ChQuaterniond Qcross(const ChQuaterniond& qa, const ChQuaterniond& qb) {
    ChQuaterniond res;
    res.e0() = qa.e0() * qb.e0() - qa.e1() * qb.e1() - qa.e2() * qb.e2() - qa.e3() * qb.e3();
    res.e1() = qa.e0() * qb.e1() + qa.e1() * qb.e0() - qa.e3() * qb.e2() + qa.e2() * qb.e3();
    res.e2() = qa.e0() * qb.e2() + qa.e2() * qb.e0() + qa.e3() * qb.e1() - qa.e1() * qb.e3();
    res.e3() = qa.e0() * qb.e3() + qa.e3() * qb.e0() - qa.e2() * qb.e1() + qa.e1() * qb.e2();
    return (res);
}

// -----------------------------------------------------------------------------

// Get the quaternion time derivative from the vector of angular speed, with w specified in _absolute_ coords.
CH_NODISCARD static constexpr inline ChQuaterniond QuatDtFromAngVelAbs(const ChVector3d& w, const ChQuaterniond& q) {
    ChQuaterniond qw;
    const auto half = 0.5;

    qw.e0() = 0;
    qw.e1() = w.x();
    qw.e2() = w.y();
    qw.e3() = w.z();

    return Qscale(Qcross(qw, q), half);  // {q_dt} = 1/2 {0,w}*{q}
}

// Get the quaternion time derivative from the vector of angular speed, with w specified in _local_ coords.
CH_NODISCARD static constexpr inline ChQuaterniond QuatDtFromAngVelRel(const ChVector3d& w, const ChQuaterniond& q) {
    ChQuaterniond qw;
    const auto half = 0.5;

    qw.e0() = 0;
    qw.e1() = w.x();
    qw.e2() = w.y();
    qw.e3() = w.z();

    return Qscale(Qcross(q, qw), half);  // {q_dt} = 1/2 {q}*{0,w_rel}
}

// Get the quaternion first derivative from the vector of angular acceleration with a specified in _absolute_ coords.
CH_NODISCARD static constexpr inline ChQuaterniond QuatDt2FromAngAccAbs(const ChVector3d& a, const ChQuaterniond& q, const ChQuaterniond& q_dt) {
    ChQuaterniond ret;
    ret.SetDt2FromAngAccAbs(a, q, q_dt);
    return ret;
}

//	Get the quaternion second derivative from the vector of angular acceleration with a specified in _relative_ coords.
CH_NODISCARD static constexpr inline ChQuaterniond QuatDt2FromAngAccRel(const ChVector3d& a, const ChQuaterniond& q, const ChQuaterniond& q_dt) {
    ChQuaterniond ret;
    ret.SetDt2FromAngAccRel(a, q, q_dt);
    return ret;
}

// -----------------------------------------------------------------------------

// Given the imaginary (vectorial) {e1 e2 e3} part of a quaternion,
// find the entire quaternion q = {e0, e1, e2, e3}.
// Note: singularities are possible.
CH_NODISCARD static constexpr inline ChQuaterniond QuatFromImaginary(const ChVector3d& im) {
    ChQuaterniond q;
    q.e1() = im.x();
    q.e2() = im.y();
    q.e3() = im.z();
    q.e0() = std::sqrt(1 - q.e1() * q.e1() - q.e2() * q.e2() - q.e3() * q.e3());
    return q;
}

// Given the imaginary (vectorial) {e1 e2 e3} part of a quaternion time derivative,
// find the entire quaternion q = {e0, e1, e2, e3}.
// Note: singularities are possible.
CH_NODISCARD static constexpr inline ChQuaterniond QuatDtFromImaginary(const ChVector3d& im_dt, const ChQuaterniond& q) {
    ChQuaterniond q_dt;
    q_dt.e1() = im_dt.x();
    q_dt.e2() = im_dt.y();
    q_dt.e3() = im_dt.z();
    q_dt.e0() = (-q.e1() * q_dt.e1() - q.e2() * q_dt.e2() - q.e3() * q_dt.e3()) / q.e0();
    return q_dt;
}

// Given the imaginary (vectorial) {e1 e2 e3} part of a quaternion second time derivative,
// find the entire quaternion q = {e0, e1, e2, e3}.
// Note: singularities are possible.
CH_NODISCARD static constexpr inline ChQuaterniond QuatDt2FromImaginary(const ChVector3d& im_dtdt, const ChQuaterniond& q, const ChQuaterniond& q_dt) {
    ChQuaterniond q_dtdt;
    q_dtdt.e1() = im_dtdt.x();
    q_dtdt.e2() = im_dtdt.y();
    q_dtdt.e3() = im_dtdt.z();
    q_dtdt.e0() = (-q.e1() * q_dtdt.e1() - q.e2() * q_dtdt.e2() - q.e3() * q_dtdt.e3()                               //
                   - q_dt.e0() * q_dt.e0() - q_dt.e1() * q_dt.e1() - q_dt.e2() * q_dt.e2() - q_dt.e3() * q_dt.e3())  //
                  / q.e0();                                                                                          //
    return q_dtdt;
}

// -----------------------------------------------------------------------------

// Get the X axis of a coordsystem, given the quaternion which
// represents the alignment of the coordsystem.
CH_NODISCARD static constexpr inline ChVector3d VaxisXfromQuat(const ChQuaterniond& quat) {
    ChVector3d res;
    res.x() = (std::pow(quat.e0(), 2) + std::pow(quat.e1(), 2)) * 2 - 1;
    res.y() = ((quat.e1() * quat.e2()) + (quat.e0() * quat.e3())) * 2;
    res.z() = ((quat.e1() * quat.e3()) - (quat.e0() * quat.e2())) * 2;
    return res;
}


}  // end namespace chrono

#endif
