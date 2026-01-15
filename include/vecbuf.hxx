#ifndef HPP_STRING_VECBUF
#define HPP_STRING_VECBUF

#include <streambuf>
#include <string>
#include <vector>
#include <iostream>

/**
 * @brief A vector-backed stream buffer intended to support stable, contiguous storage for
 *      serialization workflows (notably BMI-related `GetValuePtr` access)
 *
 * Motivation:
 *      This class wraps a std::vector<char> and exposes it through the std::basic_streambuf put-area interface
 *      (pbase, pptr, epptr). The intent is to provide a contiguous buffer whose storage location can be treated
 *      as the "source of truth" and referenced by reliable pointers for downstream consumers (e.g. BMI restart
 *      serialization that wants to build strings and hand out stable pointers to their underlying bytes).
 *
 *  Reality:
 *      This implementation is incomplete. Most specifically, overflow() is not implemented and will throw
 *      an exception if called. Implementation is non-trivial as initial attempts result in partial functionality
 *      at best with many values never being added. Any code path that relies on overflow for single-character writes
 *      (common when used with iostreams) will guarantee exceptions
 *
 *      std::basic_streambuf does not guarantee that pointer juggling performed will remain correct under
 *      all iostream usage patterns.
 *
 *      Pointer stability is only as good as std::vector's stability. Any operation that re-allocates (reserve/growth)
 *      invalidates all previously obtained pointers into the vector. This class tries to pre-reserve to avoid
 *      reallocation, but cannot guarantee it in general.
 *
 *  Intended use pattern:
 *      Used to feed data in and out of boost::archive::binary_iarchive and boost::archive::binary_oarchive objects.
 */
template<class CharT = char, class Traits = std::char_traits<CharT>>
class vecbuf : public std::basic_streambuf<CharT> {
public:
    using streambuf = std::basic_streambuf<CharT, Traits>;
    using char_type = typename streambuf::char_type;
    using int_type = typename streambuf::int_type;
    using traits_type = typename streambuf::traits_type;
    using vector = std::vector<char_type>;
    using value_type = typename vector::value_type;
    using size_type = typename vector::size_type;

    /**
     * @brief Construct a vecbuf with an expected capacity
     *
     * This attempts to reserve storage up-front to reduce reallocations and therefore reduce pointer invalidation risk
     *
     * WARNING: This class is known-incomplete: overflow() will throw an exception if used until properly implemented
     *
     * @param capacity The initial capacity of the underlying vector
     */
    vecbuf(size_type capacity = 0) : vector_() {
        reserve(capacity);
        std::cerr << "vecbuf::overflow is not properly implemented. Expect runtime issues when used." << std::endl;
    }

    /**
     * @brief Reduce the capacity of the underlying vector
     *
     * WARNING: May invalidate pointers to the underlying data
     */
    constexpr void shrink_to_fit() {
        std::cerr <<
            "WARNING: vecbuf is being shrunk to fit content. " <<
                "Pointer position is not guaranteed to remain valid for outside consumers via BMI"
        << std::endl;
        vector_.shrink_to_fit();
    }

    /**
     * @brief Clears the underlying vector
     *
     * Note: capacity is not necessarily freed, so pointers may still be valid but content may no longer be
     * accessed through standard means
     */
    constexpr void clear() {
        std::cerr <<
            "WARNING: vecbuf is being cleared. " <<
            "Pointer position is not guaranteed to remain valid for outside consumers via BMI"
        << std::endl;
        vector_.clear();
    }

    /**
     * Reserve capacity for the underlying vector and rebind the streambuf put points to the vector storage
     *
     * This is the core operation that tries to keep pbase and epptr pointing at the underlying vectoring. Failing to
     * reserve correctly will result in the invalidation of pointers.
     *
     * WARNING:
     *  This operation may reallocate, resulting in the invalidation of previously shared pointers
     *
     * @param capacity The intended total length of the underlying vector
     */
    constexpr void reserve(size_type capacity) {
        std::cerr <<
            "WARNING: vecbuf is having memory reserved. "
            "Pointer position is not guaranteed to remain valid for outside consumers via BMI."
        << std::endl;
        vector_.reserve(capacity);
        setp_from_vector();
    }

    /**
     * Reserve a capacity in addition to the current conditions of the underlying vector
     * @param additional_capacity The amount of needed extra capacity
     */
    constexpr void reserve_additional(size_type additional_capacity) { reserve(size() + additional_capacity); }

    /**
     *
     * @return The data within the internal vector
     */
    constexpr const value_type* data() const { return vector_.data(); }

    /**
     *
     * @return The size of the internal vector
     */
    constexpr size_type size() const { return vector_.size(); }

    /**
     *
     * @return The capacity of the internal vector
     */
    constexpr size_type capacity() const { return vector_.capacity(); }

    /**
     * Writes `count` characters to the output sequence from the character array whose first element is pointed to by s.
     * The characters are written as if by repeated calls to sputc(). Writing stops when either count characters are
     * written or a call to sputc() would have returned Traits::eof().
     *
     *  If the put area becomes full (pptr() == epptr()), it is unspecified whether overflow() is actually called or
     *  its effect is achieved by other means.
     *
     * @param s The data to add to the vector buffer
     * @param count The number of items to add to the buffer
     * @return The number of characters written
     */
    std::streamsize xsputn(const char_type* s, std::streamsize count) override {
        std::cerr <<
            "vecbuf::xsputn is being used.  This function should theoretically handle multiple character insertion, " <<
                "but pointer juggling is not guaranteed to work."
        << std::endl;

        try {
            reserve_additional(count);
        }
        catch (const std::bad_alloc& error) {
            // reserve did not work, use slow algorithm
            return xsputn_slow(s, count);
        }
        // reserve worked, use fast algorithm
        return xsputn_fast(s, count);
    }

    /**
     * 	The intent of this function is to transmit characters from the put area of the stream buffer to the
     * 	associated character sequence.
     *
     *  Formally, this function ensures that there is space at the put area for at least one character.
     *  The base class version always fails, and a possibly-succeeding implementation can only be provided in derived
     *  classes (see implementation requirements). The standard library provides
     *  std::strstreambuf::overflow(), std::basic_stringbuf::overflow() (until C++26), and std::basic_filebuf::overflow().
     *
     *  NOTE:
     *  The sputc() and sputn() call this function in case of an overflow (pptr() == nullptr or pptr() >= epptr()).
     *
     * @param ch The singular character to add
     * @return Traits::eof or the character added
     * @throws std::runtime_error Throws a runtime_error since it is not properly implemented
     */
    int_type overflow(int_type ch) override {
        throw std::runtime_error(
            std::string("vecbuf::overflow has not been implemented. '") +
            std::to_string(ch) +
            std::string("' cannot be added to the buffer.")
        );
    }

protected:
    /**
     *
     * @return The pointer to the beginning of the data
     */
    constexpr value_type* pbase_from_vector() const {
        return const_cast<value_type*>(vector_.data());
    }

    /**
     *
     * @return The pointer to the last value in the underlying data
     */
    constexpr value_type* pptr_from_vector() const {
        return const_cast<value_type*>(vector_.data() + vector_.size());
    }

    /**
     *
     * @return The pointer to the last reserved position in the underlying data
     */
    constexpr value_type* epptr_from_vector() const { return const_cast<value_type*>(vector_.data()) + vector_.capacity(); }

    /**
     * Updates the internal pointers to the beginning and end of the underlying data
     *
     * Needed to be called whenever the underlying data changes in capacity or allocation
     */
    constexpr void setp_from_vector() {
        std::cerr <<
            "WARNING: Pointer positions are being manipulated. " <<
                "This may adversely affect outside structures that have retrieved access via BMI."
        << std::endl;
        streambuf::setp(pbase_from_vector(), epptr_from_vector()); streambuf::pbump(size());
    }

private:
    /**
     * The underlying storage structure, primarily accessed via the pointer interface
     *
     * Contiguous, though not stable
     */
    vector vector_;

    /**
     * A fast bulk write to the buffer using pointer logic that works when a prior `reserve_additional(count)` call
     * succeeded
     *
     * Requires two assumptions to be true to avoid undefined behavior:
     *  1. vector_.resize() does not reallocate
     *  2. pptr points to vector storage prior to a resize
     *
     * @param s The series of characters to add
     * @param count The number of characters to add
     * @return The number of characters inserted
     */
    std::streamsize xsputn_fast(const char_type* s, std::streamsize count) {
        // store current pptr (end of vector location)
        auto* old_pptr = pptr_from_vector();
        // resize the vector, does not move since space already reserved
        vector_.resize(vector_.size() + count);
        // directly memcpy new content to old pptr (end of vector before it was resized)
        traits_type::copy(old_pptr, s, count);
        // reserve() already calls setp_from_vector(), only adjust pptr to new epptr
        streambuf::pbump(count);

        return count;
    }

    /**
     * @brief Slow bulk write fallback - writes one character at a time
     *
     * @param s The series of characters to add
     * @param count The number of characters to add
     * @return The number of characters to add
     */
    std::streamsize xsputn_slow(const char_type* s, const std::streamsize count) {
        // reserving entire vector failed, emplace char for char
        std::streamsize written = 0;
        while (written < count) {
            try {
                // copy one char, should throw eventually std::bad_alloc
                vector_.emplace_back(s[written]);
            }
            catch (const std::bad_alloc& error) {
                // try overflow(), if eof return, else continue
                int_type c = this->overflow(traits_type::to_int_type(s[written]));
                if (traits_type::eq_int_type(c, traits_type::eof())) {
                    return written;
                }
            }
            // update pbase, pptr and epptr
            setp_from_vector();
            written++;
        }
        return written;
    }

};

#endif
