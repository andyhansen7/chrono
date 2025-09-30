// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2023 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================

#ifndef CH_ROTATION_H
#define CH_ROTATION_H

#include "chrono/core/ChApiCE.h"
#include "chrono/core/ChVector3.h"
#include "chrono/core/ChQuaternion.h"
#include "chrono/core/ChMatrix33.h"
#include "chrono/utils/ChConstants.h"
#include "chrono/core/ChRotation.h"
#include "chrono/utils/ChUtils.h"
#include "chrono/core/ChGlobal.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace chrono {

/// Definitions of various rotation representations for conversions.
enum class RotRepresentation {
    ANGLE_AXIS,         ///< angle-axis
    ROTATION_VECTOR,    ///< rotation vector (vector parallel with axis and length equal to angle)
    EULER_ANGLES_ZXZ,   ///< Euler angle sequence: Z - X' - Z'' (intrinsic rotations)
    CARDAN_ANGLES_ZXY,  ///< Cardan (Tait-Bryan) angle sequence: Z - X' - Y'' (intrinsic rotations)
    CARDAN_ANGLES_ZYX,  ///< Cardan (Tait-Bryan) angle sequence: Z - Y' - X'' (intrinsic rotations)
    CARDAN_ANGLES_XYZ,  ///< Cardan (Tait-Bryan) angle sequence: X - Y' - Z'' (intrinsic rotations)
    RODRIGUES           ///< Rodrigues parameters
};

/// Representation of an Euler/Cardan angle set.
struct AngleSet {
    RotRepresentation seq;  ///< angle sequence
    ChVector3d angles;      ///< Euler/Cardan angles
};

/// Representation of a rotation as angle-axis.
struct AngleAxis {
    double angle;     ///< rotation angle
    ChVector3d axis;  ///< rotation axis
};

// --------------------------

/// Convert from one set of Euler or Cardan angles to another.
ChApi CH_NODISCARD static inline AngleSet AngleSetFromAngleSet(const RotRepresentation& to_seq, const AngleSet& set);

/// Convert from a set of Euler angles to Rodrigues parameters.
ChApi CH_NODISCARD static inline ChVector3d RodriguesFromAngleSet(const AngleSet& set);

/// Convert from Rodrigues parameters to a set of Euler angles.
ChApi CH_NODISCARD static inline AngleSet AngleSetFromRodrigues(const RotRepresentation& to_seq, const ChVector3d& params);

// --------------------------

/// Convert from a quaternion to an angle-axis pair.
ChApi CH_NODISCARD static constexpr inline AngleAxis AngleAxisFromQuat(const ChQuaterniond& q);

/// Convert from an angle-axis pair to a quaternion.
/// The axis is supposed to be fixed, i.e. it is constant during rotation.
/// The 'axis' vector must be normalized.
ChApi CH_NODISCARD static constexpr inline ChQuaterniond QuatFromAngleAxis(const AngleAxis& angle_axis);

/// Convert from an angle and an axis to a quaternion.
/// The axis is supposed to be fixed, i.e. it is constant during rotation.
/// The 'axis' vector must be normalized.
ChApi CH_NODISCARD static constexpr inline ChQuaterniond QuatFromAngleAxis(const double& angle, const ChVector3d& axis);

/// Convert from a speed of rotation and an axis to a quaternion derivative.
/// The rotation axis is assumed to be represented in absolute coordinates.
ChApi CH_NODISCARD static constexpr inline ChQuaterniond QuatDtFromAngleAxis(const ChQuaterniond& quat, const double& angle_dt, const ChVector3d& axis);

/// Convert from a rotation acceleration and an axis to a quaternion second derivative.
/// The rotation axis is assumed to be represented in absolute coordinates.
ChApi CH_NODISCARD static constexpr inline ChQuaterniond QuatDt2FromAngleAxis(
    const double& angle_dtdt,
    const ChVector3d& axis,
    const ChQuaterniond& q,
    const ChQuaterniond& q_dt
);

/// Convert from a rotation about X axis to a quaternion.
ChApi CH_NODISCARD static constexpr inline ChQuaterniond QuatFromAngleX(const double& angle);

/// Convert from a rotation about Y axis to a quaternion.
ChApi CH_NODISCARD static constexpr inline ChQuaterniond QuatFromAngleY(const double& angle);

/// Convert from a rotation about Z axis to a quaternion.
ChApi CH_NODISCARD static constexpr inline ChQuaterniond QuatFromAngleZ(const double& angle);

// --------------------------

/// Convert from a quaternion to a rotation vector.
ChApi CH_NODISCARD static constexpr inline ChVector3d RotVecFromQuat(const ChQuaterniond& q);

/// Convert from a rotation vector to a quaternion.
ChApi CH_NODISCARD static constexpr inline ChQuaterniond QuatFromRotVec(const ChVector3d& vec);

// --------------------------

/// Convert from a quaternion to Rodrigues parameters.
ChApi CH_NODISCARD static inline ChVector3d RodriguesFromQuat(const ChQuaterniond& q);

/// Convert from Rodrigues parameters to a quaternion.
ChApi CH_NODISCARD static inline ChQuaterniond QuatFromRodrigues(const ChVector3d& params);

/// Convert from a set of Rodrigues parameter derivatives to a quaternion derivative.
ChApi CH_NODISCARD static inline ChQuaterniond QuatDtFromRodrigues(const ChVector3d& params, const ChQuaterniond& q);

/// Convert a set of Rodrigues parameter second derivatives to a quaternion second derivative.
ChApi CH_NODISCARD static inline ChQuaterniond QuatDt2FromRodrigues(const ChVector3d& params, const ChQuaterniond& q);

// --------------------------

/// Convert from a quaternion to a set of Euler angles.
ChApi CH_NODISCARD static inline AngleSet AngleSetFromQuat(const RotRepresentation& to_seq, const ChQuaterniond& q);

/// Convert from a set of Euler angles to a quaternion.
ChApi CH_NODISCARD static inline ChQuaterniond QuatFromAngleSet(const AngleSet& set);

/// Convert from a set of Euler angle derivatives to a quaternion derivative.
ChApi CH_NODISCARD static inline ChQuaterniond QuatDtFromAngleSet(const AngleSet& set, const ChQuaterniond& q);

/// Convert from a set of Euler angle second derivatives to a quaternion second derivative.
ChApi CH_NODISCARD static inline ChQuaterniond QuatDt2FromAngleSet(const AngleSet& set, const ChQuaterniond& q);

// --------------------------

/// Convert from a vector-to-vector rotation to a quaternion.
/// This quaternion represents the rotation that rotates the source vector to be aligned with the destination vector.
/// The vectors do not need to be normalized.
ChApi CH_NODISCARD static constexpr inline ChQuaterniond QuatFromVec2Vec(const ChVector3d& start, const ChVector3d& end);

// --------------------------

inline AngleSet AngleSetFromAngleSet(const RotRepresentation& to_seq, const AngleSet& set) {
    ChMatrix33<> R;

    switch (set.seq) {
        case RotRepresentation::EULER_ANGLES_ZXZ:
            R.SetFromEulerAnglesZXZ(set.angles);
            break;
        case RotRepresentation::CARDAN_ANGLES_ZXY:
            R.SetFromCardanAnglesZXY(set.angles);
            break;
        case RotRepresentation::CARDAN_ANGLES_ZYX:
            R.SetFromCardanAnglesZYX(set.angles);
            break;
        case RotRepresentation::CARDAN_ANGLES_XYZ:
            R.SetFromCardanAnglesXYZ(set.angles);
            break;
        default:
            std::cerr << "Unknown input angle set representation" << std::endl;
            throw std::runtime_error("Unknown input angle set representation");
            break;
    }

    ChVector3d to_angles;
    switch (to_seq) {
        case RotRepresentation::EULER_ANGLES_ZXZ:
            to_angles = R.GetEulerAnglesZXZ();
            break;
        case RotRepresentation::CARDAN_ANGLES_ZXY:
            to_angles = R.GetCardanAnglesZXY();
            break;
        case RotRepresentation::CARDAN_ANGLES_ZYX:
            to_angles = R.GetCardanAnglesZYX();
            break;
        case RotRepresentation::CARDAN_ANGLES_XYZ:
            to_angles = R.GetCardanAnglesXYZ();
            break;
        default:
            std::cerr << "Unknown output angle set representation" << std::endl;
            throw std::runtime_error("Unknown output angle set representation");
            break;
    }

    return AngleSet({to_seq, to_angles});
}

inline ChVector3d RodriguesFromAngleSet(const AngleSet& set) {
    ChMatrix33<> R;

    switch (set.seq) {
        case RotRepresentation::EULER_ANGLES_ZXZ:
            R.SetFromEulerAnglesZXZ(set.angles);
            break;
        case RotRepresentation::CARDAN_ANGLES_ZXY:
            R.SetFromCardanAnglesZXY(set.angles);
            break;
        case RotRepresentation::CARDAN_ANGLES_ZYX:
            R.SetFromCardanAnglesZYX(set.angles);
            break;
        case RotRepresentation::CARDAN_ANGLES_XYZ:
            R.SetFromCardanAnglesXYZ(set.angles);
            break;
        default:
            std::cerr << "Unknown input angle set representation" << std::endl;
            throw std::runtime_error("Unknown input angle set representation");
            break;
    }

    return R.GetRodriguesParameters();
}

inline AngleSet AngleSetFromRodrigues(const RotRepresentation& to_seq, const ChVector3d& params) {
    ChMatrix33<> R;
    R.SetFromRodriguesParameters(params);

    ChVector3d to_angles;
    switch (to_seq) {
        case RotRepresentation::EULER_ANGLES_ZXZ:
            to_angles = R.GetEulerAnglesZXZ();
            break;
        case RotRepresentation::CARDAN_ANGLES_ZXY:
            to_angles = R.GetCardanAnglesZXY();
            break;
        case RotRepresentation::CARDAN_ANGLES_ZYX:
            to_angles = R.GetCardanAnglesZYX();
            break;
        case RotRepresentation::CARDAN_ANGLES_XYZ:
            to_angles = R.GetCardanAnglesXYZ();
            break;
        default:
            std::cerr << "Unknown output angle set representation" << std::endl;
            throw std::runtime_error("Unknown output angle set representation");
            break;
    }

    return AngleSet({to_seq, to_angles});
}

// --------------------------

constexpr inline AngleAxis AngleAxisFromQuat(const ChQuaterniond& q) {
    double angle = 0.0;
    ChVector3d axis;

    if (std::abs(q.e0()) < 0.99999999) {
        const double arg = std::acos(q.e0());
        const double invsine = 1 / std::sin(arg);
        ChVector3d vtemp;
        vtemp.x() = invsine * q.e1();
        vtemp.y() = invsine * q.e2();
        vtemp.z() = invsine * q.e3();
        angle = 2 * arg;
        axis = Vnorm(vtemp);
    } else {
        axis.x() = 1;
        axis.y() = 0;
        axis.z() = 0;
        angle = 0;
    }

    return AngleAxis({angle, axis});
}

constexpr inline ChQuaterniond QuatFromAngleAxis(const double& angle, const ChVector3d& axis) {
    ChQuaterniond q;
    const auto halfang = angle / 2;
    const auto sinhalf = std::sin(halfang);

    q.e0() = std::cos(halfang);
    q.e1() = axis.x() * sinhalf;
    q.e2() = axis.y() * sinhalf;
    q.e3() = axis.z() * sinhalf;

    return q;
}

constexpr inline ChQuaterniond QuatFromAngleAxis(const AngleAxis& angle_axis) {
    return QuatFromAngleAxis(angle_axis.angle, angle_axis.axis);
}

constexpr inline ChQuaterniond QuatDtFromAngleAxis(const ChQuaterniond& quat, const double& angle_dt, const ChVector3d& axis) {
    return QuatDtFromAngVelAbs(Vmul(axis, angle_dt), quat);
}

constexpr inline ChQuaterniond QuatDt2FromAngleAxis(
    const double& angle_dtdt,
    const ChVector3d& axis,
    const ChQuaterniond& q,
    const ChQuaterniond& q_dt
) {
    return QuatDt2FromAngAccAbs(Vmul(axis, angle_dtdt), q, q_dt);
}

constexpr inline ChQuaterniond QuatFromAngleX(const double& angle) {
    return QuatFromAngleAxis({angle, ChVector3d(1, 0, 0)});
}

constexpr inline ChQuaterniond QuatFromAngleY(const double& angle) {
    return QuatFromAngleAxis({angle, ChVector3d(0, 1, 0)});
}

constexpr inline ChQuaterniond QuatFromAngleZ(const double& angle) {
    return QuatFromAngleAxis({angle, ChVector3d(0, 0, 1)});
}

// --------------------------

constexpr inline ChVector3d RotVecFromQuat(const ChQuaterniond& q) {
    return q.GetRotVec();
}

constexpr inline ChQuaterniond QuatFromRotVec(const ChVector3d& vec) {
    ChQuaterniond q;
    q.SetFromRotVec(vec);
    return q;
}

// --------------------------

inline ChVector3d RodriguesFromQuat(const ChQuaterniond& q) {
    ChMatrix33<> R(q);
    return R.GetRodriguesParameters();
}

inline ChQuaterniond QuatFromRodrigues(const ChVector3d& params) {
    ChMatrix33<> R;
    R.SetFromRodriguesParameters(params);
    return R.GetQuaternion();
}

inline ChQuaterniond QuatDtFromRodrigues(const ChVector3d& params, const ChQuaterniond& q) {
    const auto params1 = RodriguesFromQuat(q);
    const auto params2 = Vadd(params1, Vmul(params, FD_STEP));
    const auto q2 = QuatFromRodrigues(params2);

    return Qscale(Qsub(q2, q), 1 / FD_STEP);
}

inline ChQuaterniond QuatDt2FromRodrigues(const ChVector3d& params, const ChQuaterniond& q) {
    const auto params0 = RodriguesFromQuat(q);
    const auto paramsA = Vsub(params0, Vmul(params, FD_STEP));
    const auto paramsB = Vadd(params0, Vmul(params, FD_STEP));
    const auto qa = QuatFromRodrigues(paramsA);
    const auto qb = QuatFromRodrigues(paramsB);

    return Qscale(Qadd(Qadd(qa, qb), Qscale(q, -2)), 1 / FD_STEP);
}

// --------------------------

inline AngleSet AngleSetFromQuat(const RotRepresentation& to_seq, const ChQuaterniond& q) {
    ChMatrix33<> R(q);

    ChVector3d to_angles;
    switch (to_seq) {
        case RotRepresentation::EULER_ANGLES_ZXZ:
            to_angles = R.GetEulerAnglesZXZ();
            break;
        case RotRepresentation::CARDAN_ANGLES_ZXY:
            to_angles = R.GetCardanAnglesZXY();
            break;
        case RotRepresentation::CARDAN_ANGLES_ZYX:
            to_angles = R.GetCardanAnglesZYX();
            break;
        case RotRepresentation::CARDAN_ANGLES_XYZ:
            to_angles = R.GetCardanAnglesXYZ();
            break;
        default:
            std::cerr << "Unknown output angle set representation" << std::endl;
            throw std::runtime_error("Unknown output angle set representation");
            break;
    }

    return AngleSet({to_seq, to_angles});
}

inline ChQuaterniond QuatFromAngleSet(const AngleSet& set) {
    ChMatrix33<> R;

    switch (set.seq) {
        case RotRepresentation::EULER_ANGLES_ZXZ:
            R.SetFromEulerAnglesZXZ(set.angles);
            break;
        case RotRepresentation::CARDAN_ANGLES_ZXY:
            R.SetFromCardanAnglesZXY(set.angles);
            break;
        case RotRepresentation::CARDAN_ANGLES_ZYX:
            R.SetFromCardanAnglesZYX(set.angles);
            break;
        case RotRepresentation::CARDAN_ANGLES_XYZ:
            R.SetFromCardanAnglesXYZ(set.angles);
            break;
        default:
            std::cerr << "Unknown input angle set representation" << std::endl;
            throw std::runtime_error("Unknown input angle set representation");
            break;
    }

    return R.GetQuaternion();
}

inline ChQuaterniond QuatDtFromAngleSet(const AngleSet& set, const ChQuaterniond& q) {
    const auto set1 = AngleSetFromQuat(set.seq, q);
    const auto ang2 = Vadd(set1.angles, Vmul(set.angles, FD_STEP));
    const auto q2 = QuatFromAngleSet({set.seq, ang2});

    return Qscale(Qsub(q2, q), 1 / FD_STEP);
}

inline ChQuaterniond QuatDt2FromAngleSet(const AngleSet& set, const ChQuaterniond& q) {
    const auto set0 = AngleSetFromQuat(set.seq, q);
    const auto angA = Vsub(set0.angles, Vmul(set.angles, FD_STEP));
    const auto angB = Vadd(set0.angles, Vmul(set.angles, FD_STEP));
    const auto qa = QuatFromAngleSet({set.seq, angA});
    const auto qb = QuatFromAngleSet({set.seq, angB});

    return Qscale(Qadd(Qadd(qa, qb), Qscale(q, -2)), 1 / FD_STEP);
}

// --------------------------

constexpr inline ChQuaterniond QuatFromVec2Vec(const ChVector3d& start, const ChVector3d& end) {
    const double ANGLE_TOLERANCE = 1e-6;
    ChQuaterniond quat;
    double halfang = 0.0;
    double sinhalf = 0.0;
    ChVector3d axis;

    const double lenXlen = start.Length() * end.Length();
    axis = start % end;
    const double sinangle = ChClamp(axis.Length() / lenXlen, -1.0, +1.0);
    const double cosangle = ChClamp(start ^ end / lenXlen, -1.0, +1.0);

    // Consider three cases: Parallel, Opposite, non-collinear
    if (std::abs(sinangle) == 0.0 && cosangle > 0) {
        // fr_vect & to_vect are parallel
        quat.e0() = 1.0;
        quat.e1() = 0.0;
        quat.e2() = 0.0;
        quat.e3() = 0.0;
    } else if (std::abs(sinangle) < ANGLE_TOLERANCE && cosangle < 0) {
        // fr_vect & to_vect are opposite, i.e. ~180 deg apart
        axis = start.GetOrthogonalVector() + (-end).GetOrthogonalVector();
        axis.Normalize();
        quat.e0() = 0.0;
        quat.e1() = ChClamp(axis.x(), -1.0, +1.0);
        quat.e2() = ChClamp(axis.y(), -1.0, +1.0);
        quat.e3() = ChClamp(axis.z(), -1.0, +1.0);
    } else {
        // fr_vect & to_vect are not co-linear case
        axis.Normalize();
        halfang = 0.5 * std::atan2(sinangle, cosangle);
        sinhalf = std::sin(halfang);

        quat.e0() = std::cos(halfang);
        quat.e1() = sinhalf * axis.x();
        quat.e2() = sinhalf * axis.y();
        quat.e3() = sinhalf * axis.z();
    }
    return (quat);
}

}  // end namespace chrono

#endif
