#pragma once

#include <stdexcept>

namespace MatrExcepts {

class matrix_error : public std::runtime_error {
public:
    matrix_error(const std::string &what_arg) : std::runtime_error(what_arg) {}
    matrix_error(const matrix_error &) = default;
};

class wrong_init_list : public matrix_error {
public:
    wrong_init_list(const std::string &what_arg) : matrix_error(what_arg) {}
    wrong_init_list(const wrong_init_list &) = default;
};

class no_det : public matrix_error {
public:
    no_det(const std::string &what_arg) : matrix_error(what_arg) {}
    no_det(const no_det &) = default;
};

} // namespace MatrExcepts