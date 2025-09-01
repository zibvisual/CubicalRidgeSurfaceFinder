#pragma once

#include <parallel_hashmap/phmap.h>
template <typename T, typename S>
using BigHashMap = phmap::parallel_flat_hash_map<T,S>;
template <typename S>
using BigHashSet = phmap::parallel_flat_hash_set<S>;
template <typename T, typename S>
using SmallHashMap = phmap::flat_hash_map<T,S>;
template <typename S>
using SmallHashSet = phmap::flat_hash_set<S>;

// #include <unordered_map>
// #include <unordered_set>

// template <typename T, typename S>
// using BigHashMap = std::unordered_map<T, S>;
// template <typename S>
// using BigHashSet = std::unordered_set<S>;
// template <typename T, typename S>
// using SmallHashMap = std::unordered_map<T, S>;
// template <typename S>
// using SmallHashSet = std::unordered_set<S>;