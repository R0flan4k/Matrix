#pragma once

#include <stdexcept>

namespace MatrExcepts {

class matrix_error : public std::runtime_error {
public:
    matrix_error(const std::string &what_arg) : std::runtime_error(what_arg) {}
};

class wrong_init_list : public matrix_error {
public:
    wrong_init_list(const std::string &what_arg) : matrix_error(what_arg) {}
};

class no_det : public matrix_error {
public:
    no_det(const std::string &what_arg) : matrix_error(what_arg) {}
};

} // namespace MatrExcepts