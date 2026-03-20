#pragma once

#include <tuple>
#include <variant>
#include <vector>

#include <utils/Vec.hpp>

// TODO: instead of using std::tuple, we should create LineSet in the same way as tuple itself --> LineSet<T1,T2,T3> is nicer

// helper struct to check if generic is a tuple
template<typename T>
struct is_tuple_impl : std::false_type {};

template<typename... Ts>
struct is_tuple_impl<std::tuple<Ts...>> : std::true_type {};

template<typename T>
struct is_tuple : is_tuple_impl<std::decay_t<T>> {};

template <typename T>
struct LineSetPoint {
    VecFloat point;
    T data;

    /**
     * Print the inner data.
     */
    template <class N = T, typename std::enable_if<std::is_same<N, std::monostate>::value, int>::type = 0>
    void print_data(std::ostream& output) const {
        /* noop */
    }

    /**
     * Print the inner data.
     */
    template <class S = T, typename std::enable_if<std::is_same<S, std::string>::value, int>::type = 0>
    void print_data(std::ostream& output) const {
        output << data << std::endl;
    }

    /**
     * Print the inner data.
     */
    template <class I = T, typename std::enable_if<std::is_arithmetic<I>::value, int>::type = 0>
    void print_data(std::ostream& output) const {
        output << data << std::endl;
    }

    /**
     * Print the inner data.
     */
    template <class Tuple = T, typename std::enable_if<is_tuple<Tuple>::value, int>::type = 0>
    void print_data(std::ostream& output) const {
        std::apply([&output](auto const&... args){((output << std::to_string(args) << std::endl), ...);}, data);
    }
};

template <typename T>
class LineSet {
public:

    LineSet() : m_lines() {}

    void save(std::string output_path) const {
        write_line_data(output_path, *this);
    }

    void reserve(std::size_t size) {
        m_lines.reserve(size);
    }

    void add(std::vector<LineSetPoint<T>> line){
        m_lines.push_back(line);
    }

    // template <class Tuple = T, typename std::enable_if<std::is_same<typename Tuple::Distance, DistanceEnabled>::value, int>::type = 0>
    std::size_t size() const {
        return m_lines.size();
    }

    // void
    template <class N = T, typename std::enable_if<std::is_same<N, std::monostate>::value, int>::type = 0>
    std::size_t data_size() const {
        return 0;
    }

    // std::string
    template <class S = T, typename std::enable_if<std::is_same<S, std::string>::value, int>::type = 0>
    std::size_t data_size() const {
        return 1;
    }

    // primitive type
    template <class I = T, typename std::enable_if<std::is_arithmetic<I>::value, int>::type = 0>
    std::size_t data_size() const {
        return 1;
    }
    
    // tuple
    template <class Tuple = T, typename std::enable_if<is_tuple<Tuple>::value, int>::type = 0>
    std::size_t data_size() const {
        return std::tuple_size_v<T>;
    }

    const std::vector<std::vector<LineSetPoint<T>>>& lines() const {
        return m_lines;
    }
protected:
    std::vector<std::vector<LineSetPoint<T>>> m_lines;
};

template <typename T>
class LineSetBuilder {
public:
    LineSetBuilder() : m_set(), m_current_line(){}

    void reserve_lines(std::size_t size) {
        m_set.reserve(size);
    }

    void push_line() {
        m_set.add(std::move(m_current_line));
        m_current_line = std::vector<LineSetPoint<T>>();
    }

    void add_point(VecFloat point, T data){
        m_current_line.push_back(LineSetPoint{point, data});
    }

    LineSet<T> build() {
        auto result = std::move(m_set);
        m_set = LineSet<T>();
        return result;
    }

protected:
    LineSet<T> m_set;
    std::vector<LineSetPoint<T>> m_current_line;
};