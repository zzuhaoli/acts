// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/TrackStateProxy.hpp"
#include "Acts/EventData/TrackStateProxyConcept.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/Concepts.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <bitset>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <optional>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

namespace Acts {

// forward declarations
template <typename derived_t>
class MultiTrajectory;
class Surface;

namespace detail_lt {

/// Helper type that wraps two iterators
template <bool reverse, typename trajectory_t, std::size_t M, bool ReadOnly>
class TrackStateRange {
  using ProxyType = TrackStateProxy<trajectory_t, M, ReadOnly>;
  using IndexType = typename ProxyType::IndexType;
  static constexpr IndexType kInvalid = ProxyType::kInvalid;

 public:
  /// Iterator that wraps a track state proxy. The nullopt case signifies the
  /// end of the range, i.e. the "past-the-end" iterator
  struct Iterator {
    std::optional<ProxyType> proxy;

    using iterator_category = std::forward_iterator_tag;
    using value_type = ProxyType;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = void;

    Iterator& operator++() {
      if (!proxy) {
        return *this;
      }
      if constexpr (reverse) {
        if (proxy->hasPrevious()) {
          proxy = proxy->trajectory().getTrackState(proxy->previous());
          return *this;
        } else {
          proxy = std::nullopt;
          return *this;
        }
      } else {
        IndexType next =
            proxy->template component<IndexType, hashString("next")>();
        if (next != kInvalid) {
          proxy = proxy->trajectory().getTrackState(next);
          return *this;
        } else {
          proxy = std::nullopt;
          return *this;
        }
      }
    }

    bool operator==(const Iterator& other) const {
      if (!proxy && !other.proxy) {
        return true;
      }
      if (proxy && other.proxy) {
        return proxy->index() == other.proxy->index();
      }
      return false;
    }

    bool operator!=(const Iterator& other) const { return !(*this == other); }

    ProxyType operator*() const { return *proxy; }
    ProxyType operator*() { return *proxy; }
  };

  TrackStateRange(ProxyType _begin) : m_begin{_begin} {}
  TrackStateRange() : m_begin{std::nullopt} {}

  Iterator begin() { return m_begin; }
  Iterator end() { return Iterator{std::nullopt}; }

 private:
  Iterator m_begin;
};

// implement track state visitor concept
template <typename T, typename TS>
using call_operator_t = decltype(std::declval<T>()(std::declval<TS>()));

template <typename T, typename TS>
constexpr bool VisitorConcept = Concepts ::require<
    Concepts ::either<Concepts ::identical_to<bool, call_operator_t, T, TS>,
                      Concepts ::identical_to<void, call_operator_t, T, TS>>>;

}  // namespace detail_lt

/// This namespace contains typedefs and constant values that are used by
/// other parts of the @c MultiTrajectory implementation. It extracts these
/// from @c TrackStateTraits using the default maximum measurement dimension.
namespace MultiTrajectoryTraits {
constexpr unsigned int MeasurementSizeMax = eBoundSize;
using IndexType = TrackIndexType;
constexpr IndexType kInvalid = kTrackIndexInvalid;
}  // namespace MultiTrajectoryTraits

template <typename T>
struct IsReadOnlyMultiTrajectory;

/// Store a trajectory of track states with multiple components.
///
/// This container supports both simple, sequential trajectories as well
/// as combinatorial or multi-component trajectories. Each point can store
/// a parent point such that the trajectory forms a directed, acyclic graph
/// of sub-trajectories. From a set of endpoints, all possible sub-components
/// can be easily identified. Some functionality is provided to simplify
/// iterating over specific sub-components.
template <typename derived_t>
class MultiTrajectory {
 public:
  using Derived = derived_t;

  static constexpr bool ReadOnly = IsReadOnlyMultiTrajectory<Derived>::value;

  // Pull out type alias and re-expose them for ease of use.
  //
  static constexpr unsigned int MeasurementSizeMax =
      MultiTrajectoryTraits::MeasurementSizeMax;

 private:
  friend class TrackStateProxy<Derived, MeasurementSizeMax, true>;
  friend class TrackStateProxy<Derived, MeasurementSizeMax, false>;
  template <typename T>
  friend class MultiTrajectory;

 public:
  using ConstTrackStateProxy =
      Acts::TrackStateProxy<Derived, MeasurementSizeMax, true>;
  using TrackStateProxy =
      Acts::TrackStateProxy<Derived, MeasurementSizeMax, false>;

  using IndexType = typename TrackStateProxy::IndexType;
  static constexpr IndexType kInvalid = TrackStateProxy::kInvalid;

 protected:
  MultiTrajectory() = default;  // pseudo abstract base class

 private:
  /// Helper to static cast this to the Derived class for CRTP
  constexpr Derived& self() { return static_cast<Derived&>(*this); }
  /// Helper to static cast this to the Derived class for CRTP. Const version.
  constexpr const Derived& self() const {
    return static_cast<const Derived&>(*this);
  }

  /// Helper function to check if a component exists IF it is an optional one.
  /// Used in assertions
  constexpr bool checkOptional(HashedString key, IndexType istate) const {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
      case "filtered"_hash:
      case "smoothed"_hash:
      case "calibrated"_hash:
      case "jacobian"_hash:
      case "projector"_hash:
        return self().has_impl(key, istate);
      default:
        return true;
    }
  }

 public:
  /// Access a read-only point on the trajectory by index.
  /// @param istate The index to access
  /// @return Read only proxy to the stored track state
  ConstTrackStateProxy getTrackState(IndexType istate) const {
    return {*this, istate};
  }

  /// Access a writable point on the trajectory by index.
  /// @param istate The index to access
  /// @return Read-write proxy to the stored track state
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackStateProxy getTrackState(IndexType istate) {
    return {*this, istate};
  }

  /// Add a track state to the container and return a track state proxy to it
  /// This effectively calls @c addTrackState and @c getTrackState
  /// @note Only available if the track state container is not read-only
  /// @return a track state proxy to the newly added track state
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackStateProxy makeTrackState(
      TrackStatePropMask mask = TrackStatePropMask::All,
      IndexType iprevious = kInvalid) {
    return getTrackState(addTrackState(mask, iprevious));
  }

  /// Visit all previous states starting at a given endpoint.
  ///
  /// @param iendpoint  index of the last state
  /// @param callable   non-modifying functor to be called with each point
  template <typename F>
  void visitBackwards(IndexType iendpoint, F&& callable) const;

  /// Range for the track states from @p iendpoint to the trajectory start
  /// @param iendpoint Trajectory entry point to start from
  /// @return Iterator pair to iterate over
  /// @note Const version
  auto reverseTrackStateRange(IndexType iendpoint) const {
    using range_t =
        detail_lt::TrackStateRange<true, Derived, MeasurementSizeMax, true>;
    if (iendpoint == kInvalid) {
      return range_t{};
    }

    return range_t{getTrackState(iendpoint)};
  }

  /// Range for the track states from @p iendpoint to the trajectory start,
  /// i.e from the outside in.
  /// @param iendpoint Trajectory entry point to start from
  /// @return Iterator pair to iterate over
  /// @note Mutable version
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto reverseTrackStateRange(IndexType iendpoint) {
    using range_t =
        detail_lt::TrackStateRange<true, Derived, MeasurementSizeMax, false>;
    if (iendpoint == kInvalid) {
      return range_t{};
    }

    return range_t{getTrackState(iendpoint)};
  }

  /// Range for the track states from @p istartpoint to the trajectory end,
  /// i.e from inside out
  /// @param istartpoint Trajectory state index for the innermost track
  ///        state to start from
  /// @return Iterator pair to iterate over
  /// @note Const version
  auto forwardTrackStateRange(IndexType istartpoint) const {
    using range_t =
        detail_lt::TrackStateRange<false, Derived, MeasurementSizeMax, true>;
    if (istartpoint == kInvalid) {
      return range_t{};
    }

    return range_t{getTrackState(istartpoint)};
  }

  /// Range for the track states from @p istartpoint to the trajectory end,
  /// i.e from inside out
  /// @param istartpoint Trajectory state index for the innermost track
  ///        state to start from
  /// @return Iterator pair to iterate over
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto forwardTrackStateRange(IndexType istartpoint) {
    using range_t =
        detail_lt::TrackStateRange<false, Derived, MeasurementSizeMax, false>;
    if (istartpoint == kInvalid) {
      return range_t{};
    }

    return range_t{getTrackState(istartpoint)};
  }

  /// Apply a function to all previous states starting at a given endpoint.
  ///
  /// @param iendpoint  index of the last state
  /// @param callable   modifying functor to be called with each point
  ///
  /// @warning If the trajectory contains multiple components with common
  ///          points, this can have an impact on the other components.
  template <typename F, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void applyBackwards(IndexType iendpoint, F&& callable) {
    static_assert(detail_lt::VisitorConcept<F, TrackStateProxy>,
                  "Callable needs to satisfy VisitorConcept");

    if (iendpoint == MultiTrajectoryTraits::kInvalid) {
      throw std::runtime_error(
          "Cannot apply backwards with kInvalid as endpoint");
    }

    while (true) {
      auto ts = getTrackState(iendpoint);
      if constexpr (std::is_same_v<std::invoke_result_t<F, TrackStateProxy>,
                                   bool>) {
        bool proceed = callable(ts);
        // this point has no parent and ends the trajectory, or a break was
        // requested
        if (!proceed || !ts.hasPrevious()) {
          break;
        }
      } else {
        callable(ts);
        // this point has no parent and ends the trajectory
        if (!ts.hasPrevious()) {
          break;
        }
      }
      iendpoint = ts.previous();
    }
  }

  auto&& convertToReadOnly() const {
    auto&& cv = self().convertToReadOnly_impl();
    static_assert(
        IsReadOnlyMultiTrajectory<decltype(cv)>::value,
        "convertToReadOnly_impl does not return something that reports "
        "being ReadOnly.");
    return cv;
  }

  /// Clear the @c MultiTrajectory. Leaves the underlying storage untouched
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void clear() {
    self().clear_impl();
  }

  /// Returns the number of track states contained
  constexpr IndexType size() const { return self().size_impl(); }

  /// Add a track state without providing explicit information. Which components
  /// of the track state are initialized/allocated can be controlled via @p mask
  /// @param mask The bitmask that instructs which components to allocate and
  /// which to leave invalid
  /// @param iprevious index of the previous state, kInvalid if first
  /// @return Index of the newly added track state
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr IndexType addTrackState(
      TrackStatePropMask mask = TrackStatePropMask::All,
      IndexType iprevious = kInvalid) {
    return self().addTrackState_impl(mask, iprevious);
  }

  /// Add a column to the @c MultiTrajectory
  /// @tparam T Type of the column values to add
  /// @note This takes a string argument rather than a hashed string to maintain
  ///       compatibility with backends.
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void addColumn(const std::string& key) {
    self().template addColumn_impl<T>(key);
  }

  /// Check if a column with a key @p key exists.
  /// @param key Key to check for a column with
  /// @return True if the column exists, false if not.
  constexpr bool hasColumn(HashedString key) const {
    return self().hasColumn_impl(key);
  }

 protected:
  // These are internal helper functions which the @c TrackStateProxy class talks to

  /// Check for component existence of @p key in track satet @p istate
  /// @param key The key for which to check
  /// @param istate The track state index to check
  /// @return True if the component exists, false if not
  constexpr bool has(HashedString key, IndexType istate) const {
    return self().has_impl(key, istate);
  }

  /// Check for component existence of @p key in track satet @p istate
  /// @tparam key The key for which to check
  /// @param istate The track state index to check
  /// @return True if the component exists, false if not
  template <HashedString key>
  constexpr bool has(IndexType istate) const {
    return self().has_impl(key, istate);
  }

  /// Retrieve a parameter proxy instance for parameters at a given index
  /// @param parIdx Index into the parameter column
  /// @return Mutable proxy
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackStateProxy::Parameters parameters(IndexType parIdx) {
    return self().parameters_impl(parIdx);
  }

  /// Retrieve a parameter proxy instance for parameters at a given index
  /// @param parIdx Index into the parameter column
  /// @return Const proxy
  constexpr typename ConstTrackStateProxy::Parameters parameters(
      IndexType parIdx) const {
    return self().parameters_impl(parIdx);
  }

  /// Retrieve a covariance proxy instance for a covariance at a given index
  /// @param covIdx Index into the covariance column
  /// @return Mutable proxy
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackStateProxy::Covariance covariance(IndexType covIdx) {
    return self().covariance_impl(covIdx);
  }

  /// Retrieve a covariance proxy instance for a covariance at a given index
  /// @param covIdx Index into the covariance column
  /// @return Const proxy
  constexpr typename ConstTrackStateProxy::Covariance covariance(
      IndexType covIdx) const {
    return self().covariance_impl(covIdx);
  }

  /// Retrieve a jacobian proxy instance for a jacobian at a given index
  /// @param jacIdx Index into the jacobian column
  /// @return Mutable proxy
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackStateProxy::Covariance jacobian(IndexType jacIdx) {
    return self().jacobian_impl(jacIdx);
  }

  /// Retrieve a jacobian proxy instance for a jacobian at a given index
  /// @param jacIdx Index into the jacobian column
  /// @return Const proxy
  constexpr typename ConstTrackStateProxy::Covariance jacobian(
      IndexType jacIdx) const {
    return self().jacobian_impl(jacIdx);
  }

  /// Retrieve a measurement proxy instance for a measurement at a given index
  /// @tparam measdim the measurement dimension
  /// @param measIdx Index into the measurement column
  /// @return Mutable proxy
  template <std::size_t measdim, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr typename TrackStateProxy::template Measurement<measdim> measurement(
      IndexType measIdx) {
    return self().template measurement_impl<measdim>(measIdx);
  }

  /// Retrieve a measurement proxy instance for a measurement at a given index
  /// @tparam measdim the measurement dimension
  /// @param measIdx Index into the measurement column
  /// @return Const proxy
  template <std::size_t measdim>
  constexpr typename ConstTrackStateProxy::template Measurement<measdim>
  measurement(IndexType measIdx) const {
    return self().template measurement_impl<measdim>(measIdx);
  }

  /// Retrieve a measurement covariance proxy instance for a measurement at a
  /// given index
  /// @tparam measdim the measurement dimension
  /// @param covIdx Index into the measurement covariance column
  /// @return Mutable proxy
  template <std::size_t measdim, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr typename TrackStateProxy::template MeasurementCovariance<measdim>
  measurementCovariance(IndexType covIdx) {
    return self().template measurementCovariance_impl<measdim>(covIdx);
  }

  /// Retrieve a measurement covariance proxy instance for a measurement at a
  /// given index
  /// @param covIdx Index into the measurement covariance column
  /// @return Const proxy
  template <std::size_t measdim>
  constexpr
      typename ConstTrackStateProxy::template MeasurementCovariance<measdim>
      measurementCovariance(IndexType covIdx) const {
    return self().template measurementCovariance_impl<measdim>(covIdx);
  }

  /// Get the calibrated measurement size for a track state
  /// @param istate The track state
  /// @return the calibrated size
  IndexType calibratedSize(IndexType istate) const {
    return self().calibratedSize_impl(istate);
  }

  /// Share a shareable component from between track state.
  /// @param iself The track state index to share "into"
  /// @param iother The track state index to share from
  /// @param shareSource Which component to share from
  /// @param shareTarget Which component to share as. This doesn't have to be the same
  ///                    as @p shareSource, e.g. predicted can be shared as filtered.
  /// @note Shareable components are predicted, filtered, smoothed, calibrated, jacobian,
  ///       or projector. See @c TrackStatePropMask.
  /// @note The track states both need to be stored in the
  ///       same @c MultiTrajectory instance
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void shareFrom(IndexType iself, IndexType iother,
                           TrackStatePropMask shareSource,
                           TrackStatePropMask shareTarget) {
    self().shareFrom_impl(iself, iother, shareSource, shareTarget);
  }

  /// Unset an optional track state component
  /// @param target The component to unset
  /// @param istate The track state index to operate on
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void unset(TrackStatePropMask target, IndexType istate) {
    self().unset_impl(target, istate);
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @param istate The track state index to operate on
  /// @return Mutable reference to the component given by @p key
  template <typename T, HashedString key, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr T& component(IndexType istate) {
    assert(checkOptional(key, istate));
    return *std::any_cast<T*>(self().component_impl(key, istate));
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @param istate The track state index to operate on
  /// @return Mutable reference to the component given by @p key
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(HashedString key, IndexType istate) {
    assert(checkOptional(key, istate));
    return *std::any_cast<T*>(self().component_impl(key, istate));
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @param istate The track state index to operate on
  /// @return Const reference to the component given by @p key
  template <typename T, HashedString key>
  constexpr const T& component(IndexType istate) const {
    assert(checkOptional(key, istate));
    return *std::any_cast<const T*>(self().component_impl(key, istate));
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @param istate The track state index to operate on
  /// @return Const reference to the component given by @p key
  template <typename T>
  constexpr const T& component(HashedString key, IndexType istate) const {
    assert(checkOptional(key, istate));
    return *std::any_cast<const T*>(self().component_impl(key, istate));
  }

  /// Allocate storage for a calibrated measurement of specified dimension
  /// @param istate The track state to store for
  /// @param measdim the dimension of the measurement to store
  /// @note Is a noop if the track state already has an allocation
  ///       an the dimension is the same.
  void allocateCalibrated(IndexType istate, std::size_t measdim) {
    throw_assert(measdim > 0 && measdim <= eBoundSize,
                 "Invalid measurement dimension detected");

    self().allocateCalibrated_impl(istate, measdim);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setUncalibratedSourceLink(IndexType istate, SourceLink sourceLink) {
    self().setUncalibratedSourceLink_impl(istate, std::move(sourceLink));
  }

  SourceLink getUncalibratedSourceLink(IndexType istate) const {
    return self().getUncalibratedSourceLink_impl(istate);
  }

  const Surface* referenceSurface(IndexType istate) const {
    return self().referenceSurface_impl(istate);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setReferenceSurface(IndexType istate,
                           std::shared_ptr<const Surface> surface) {
    self().setReferenceSurface_impl(istate, std::move(surface));
  }

 private:
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void copyDynamicFrom(IndexType dstIdx, const T& src, IndexType srcIdx) {
    const auto& dynamicKeys = src.self().dynamicKeys_impl();
    for (const auto key : dynamicKeys) {
      std::any srcPtr = src.self().component_impl(key, srcIdx);
      self().copyDynamicFrom_impl(dstIdx, key, srcPtr);
    }
  }
};

}  // namespace Acts

#include "Acts/EventData/MultiTrajectory.ipp"
