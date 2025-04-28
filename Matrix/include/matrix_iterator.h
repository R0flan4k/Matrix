#pragma once

namespace Matrices {

namespace internal {

template <typename T> concept pointer = std::is_pointer_v<T>;

template <typename T, typename Difference = int, typename Pointer = T *,
          typename Reference = T &>
class random_access_iterator final {
    Pointer ptr_;

public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = Difference;
    using value_type = T;
    using reference = Reference;
    using pointer = Pointer;
    using size_type = std::size_t;

    random_access_iterator() noexcept : ptr_(Pointer()) {}
    random_access_iterator(const Pointer &ptr) noexcept : ptr_(ptr) {}

    random_access_iterator<T, Difference, Pointer, Reference> &
    operator+=(size_type n)
    {
        ptr_ += n;
        return *this;
    }
    random_access_iterator<T, Difference, Pointer, Reference> &
    operator-=(size_type n)
    {
        ptr_ -= n;
        return *this;
    }
    random_access_iterator<T, Difference, Pointer, Reference> operator++(int)
    {
        return random_access_iterator<T, Difference, Pointer, Reference>{
            ptr_++};
    }
    random_access_iterator<T, Difference, Pointer, Reference> &operator++()
    {
        return *this += 1;
    }
    random_access_iterator<T, Difference, Pointer, Reference> operator--(int)
    {
        return random_access_iterator<T, Difference, Pointer, Reference>{
            ptr_--};
    }
    random_access_iterator<T, Difference, Pointer, Reference> &operator--()
    {
        return *this -= 1;
    }
    reference operator*() const { return *ptr_; }
    pointer operator->() const { return ptr_; }

    friend difference_type operator-(
        const random_access_iterator<T, Difference, Pointer, Reference> &it1,
        const random_access_iterator<T, Difference, Pointer, Reference> &it2)
    {
        return it1.ptr_ - it2.ptr_;
    }

    friend bool operator==(
        const random_access_iterator<T, Difference, Pointer, Reference> &it1,
        const random_access_iterator<T, Difference, Pointer, Reference> &it2)
    {
        return it1.ptr_ == it2.ptr_;
    }

    reference operator[](size_type n) const { return ptr_[n]; }
};

template <typename T, typename Difference = std::ptrdiff_t,
          typename Pointer = T *, typename Reference = T &>
bool operator!=(
    const random_access_iterator<T, Difference, Pointer, Reference> &it1,
    const random_access_iterator<T, Difference, Pointer, Reference> &it2)
{
    return !(it1 == it2);
}

template <typename T, typename Difference = std::ptrdiff_t,
          typename Pointer = T *, typename Reference = T &>
bool operator>(
    const random_access_iterator<T, Difference, Pointer, Reference> &it1,
    const random_access_iterator<T, Difference, Pointer, Reference> &it2)
{
    return it1 - it2 > 0;
}

template <typename T, typename Difference = std::ptrdiff_t,
          typename Pointer = T *, typename Reference = T &>
bool operator<(
    const random_access_iterator<T, Difference, Pointer, Reference> &it1,
    const random_access_iterator<T, Difference, Pointer, Reference> &it2)
{
    return it1 - it2 < 0;
}

template <typename T, typename Difference = std::ptrdiff_t,
          typename Pointer = T *, typename Reference = T &>
bool operator>=(
    const random_access_iterator<T, Difference, Pointer, Reference> &it1,
    const random_access_iterator<T, Difference, Pointer, Reference> &it2)
{
    return it1 - it2 >= 0;
}

template <typename T, typename Difference = std::ptrdiff_t,
          typename Pointer = T *, typename Reference = T &>
bool operator<=(
    const random_access_iterator<T, Difference, Pointer, Reference> &it1,
    const random_access_iterator<T, Difference, Pointer, Reference> &it2)
{
    return it1 - it2 <= 0;
}

template <typename T, typename Difference = std::ptrdiff_t,
          typename Pointer = T *, typename Reference = T &>
random_access_iterator<T, Difference, Pointer, Reference>
operator+(random_access_iterator<T, Difference, Pointer, Reference> it,
          std::size_t n)
{
    return it += n;
}

template <typename T, typename Difference = std::ptrdiff_t,
          typename Pointer = T *, typename Reference = T &>
random_access_iterator<T, Difference, Pointer, Reference>
operator+(std::size_t n,
          random_access_iterator<T, Difference, Pointer, Reference> it)
{
    return it += n;
}

template <typename T, typename Difference = std::ptrdiff_t,
          typename Pointer = T *, typename Reference = T &>
random_access_iterator<T, Difference, Pointer, Reference>
operator-(random_access_iterator<T, Difference, Pointer, Reference> it,
          std::size_t n)
{
    return it -= n;
}

} // namespace internal

} // namespace Matrices