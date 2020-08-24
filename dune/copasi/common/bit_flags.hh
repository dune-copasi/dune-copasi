#ifndef DUNE_COPASI_COMMON_BITFLAGS_HH
#define DUNE_COPASI_COMMON_BITFLAGS_HH

#include <bitset>
#include <climits>

namespace Dune::Copasi {

//! Bitflag indicator
// This class might be specialized when a certain enum is wanted to do not have
// automatic convertion to bitflags
template<class T>
struct is_bit_flag : public std::is_enum<T>
{};

//! Alias for Bitflag indicator
template<class T>
constexpr bool is_bit_flag_v = is_bit_flag<T>::value;

/**
 * @brief     Bit flags for enumerators
 * @details   This class allows that an enum can be operated as a bit set of
 * flags. When an enum is a bit flag most operators expected for bits are
 * enabled to operate directly on that enum. This object is different than a
 * std::bitset because it is fully constexpr and because operators on the
 * underlying enum are automatic.
 * @code
 *   // Notice that each option correspond to a different bit
 *   enum class MyEnum {option_a = 0b01, option_b = 0b10};
 *   // `my_enum` is deduced as `BitFlags<MyEnum>` type because operator | on
 * `MyEnum` was overloaded to do so auto my_enum = MyEnum::option_a |
 * MyEnum::option_b; assert(my_enum.test(MyEnum::option_a));
 * @endcode
 *
 * @tparam     Enum     Enumerator to be converted into a bit flag
 */
//! Bit flags for enumerators
template<typename Enum>
class BitFlags
{
  // Ensure that Enum is a valid template
  static_assert(
    is_bit_flag_v<Enum>,
    "Enum is not a bit flag and it should not be instatiated in BitFlag class");
  static_assert(std::is_enum_v<Enum>, "Emum type must be a enummeration type");

  //! Bit set representation of the bit flag
  using BitSet = std::bitset<CHAR_BIT * sizeof(Enum)>;

  //! Unterlying type to operate with
  using UnderlyingType = std::underlying_type_t<Enum>;

  //! Value full of 0 bits
  static constexpr UnderlyingType null_value =
    UnderlyingType{} ^ UnderlyingType {};

public:
  //! Bitflag with all flags turned on
  static constexpr BitFlags<Enum> All =
    BitFlags<Enum>{ static_cast<Enum>(~null_value) };

  //! Bitflag with all flags turned off
  static constexpr BitFlags<Enum> None =
    BitFlags<Enum>{ static_cast<Enum>(null_value) };

  //! Default constructor
  // All flags are default initialized to false
  constexpr inline BitFlags()
    : _value{ null_value }
  {}

  //! Enum constructor
  // Initialization with an undelying enumerator
  constexpr inline BitFlags(const Enum& value)
    : _value(static_cast<UnderlyingType>(value))
  {}

  //! Bit set constructor
  // Initialization with a standard library bitset
  inline explicit BitFlags(const BitSet& bit_set)
    : BitFlags(static_cast<Enum>(bit_set.to_ullong()))
  {
    static_assert(sizeof(UnderlyingType) <= sizeof(unsigned long long),
                  "Underlying type is longer than maximum supported cast");
  }

  //! Return bitflag as bitset
  constexpr inline BitSet as_bitset() const
  {
    return BitSet{ static_cast<unsigned long long>(_value) };
  }
  //! Return bitflag as its underlying type
  constexpr inline const UnderlyingType as_underlying() const { return _value; }

  //! Return bitflag as its underlying type
  constexpr inline UnderlyingType& as_underlying() { return _value; }

  //! Return bitflag as its underlying enum
  constexpr inline const Enum as_enum() const
  {
    return static_cast<const Enum>(_value);
  }

  //! Return bitflag as its underlying type
  constexpr inline Enum& as_enum() { return static_cast<Enum&>(_value); }

  //! Implicit conversion to underlying enum
  constexpr inline operator Enum() const { return as_enum(); }

  //! Bitwise or operator with another bitflag
  constexpr inline BitFlags operator|(const BitFlags& rhs) const
  {
    return static_cast<Enum>(_value | rhs._value);
  }

  //! Bitwise and operator with another bitflag
  constexpr inline BitFlags operator&(const BitFlags& rhs) const
  {
    return static_cast<Enum>(_value & rhs._value);
  }

  //! Bitwise xor operator with another bitflag
  constexpr inline BitFlags operator^(const BitFlags& rhs) const
  {
    return static_cast<Enum>(_value ^ rhs._value);
  }

  //! Bitwise left shift operator
  constexpr inline BitFlags operator<<(std::size_t pos) const
  {
    return static_cast<Enum>(_value << pos);
  }

  //! Bitwise right shift operator
  constexpr inline BitFlags operator>>(std::size_t pos) const
  {
    return static_cast<Enum>(_value >> pos);
  }

  //! Bitwise or operator with this bitflag
  inline BitFlags& operator|=(const BitFlags& rhs)
  {
    _value |= rhs._value;
    return *this;
  }

  //! Bitwise and operator with this bitflag
  inline BitFlags& operator&=(const BitFlags& rhs)
  {
    _value &= rhs._value;
    return *this;
  }

  //! Bitwise xor operator with this bitflag
  inline BitFlags& operator^=(const BitFlags& rhs)
  {
    _value ^= rhs._value;
    return *this;
  }

  //! Bitwise left shift operator with this bitflag
  inline BitFlags& operator<<=(std::size_t pos)
  {
    _value <<= pos;
    return *this;
  }

  //! Bitwise right shift operator with this bitflag
  inline BitFlags& operator>>=(std::size_t pos)
  {
    _value >>= pos;
    return *this;
  }

  //! Bitwise not operator of this bitflag
  constexpr inline BitFlags operator~() const
  {
    return BitFlags{ static_cast<Enum>(~_value) };
  }

  //! Test if the required flag is active in the bitflag
  constexpr inline bool test(const BitFlags& flag) const
  {
    return (_value & flag._value) == flag._value;
  }

  //! Reset the required flags in the bitflag
  inline void reset(const BitFlags& flag) { _value &= ~flag._value; }

  //! Set the required flags to true
  inline void set(const BitFlags& flag, bool value = true)
  {
    value ? (_value |= flag._value) : reset(flag);
  }

  //! Flip the required flags
  inline void flip(const BitFlags& flag) { _value ^= flag; }

  //! Checks if all flags are set to true
  constexpr inline bool all() const { return as_bitset().all(); }

  //! Checks if any flags are set to true
  constexpr inline bool any() const { return as_bitset().any(); }

  //! Checks if none flags are set to true
  constexpr inline bool none() const { return as_bitset().none(); }

private:
  //! Underlying data
  UnderlyingType _value;
};

//! Bitwise or between two enums that may be bitflags
template<class Enum>
constexpr inline std::enable_if_t<is_bit_flag_v<Enum>, BitFlags<Enum>>
operator|(Enum lhs, Enum rhs)
{
  return BitFlags<Enum>{ lhs } | BitFlags<Enum>{ rhs };
}

//! Bitwise and between two enums that may be bitflags
template<class Enum>
constexpr inline std::enable_if_t<is_bit_flag_v<Enum>, BitFlags<Enum>>
operator&(Enum lhs, Enum rhs)
{
  return BitFlags<Enum>{ lhs } & BitFlags<Enum>{ rhs };
}

//! Bitwise xor between two enums that may be bitflags
template<class Enum>
constexpr inline std::enable_if_t<is_bit_flag_v<Enum>, BitFlags<Enum>>
operator^(Enum lhs, Enum rhs)
{
  return BitFlags<Enum>{ lhs } ^ BitFlags<Enum>{ rhs };
}

//! Return true if all flags of both enums (allowed to be bitflags) are equal
template<class Enum>
constexpr inline std::enable_if_t<is_bit_flag_v<Enum>, bool>
operator==(Enum lhs, Enum rhs)
{
  return BitFlags<Enum>{ lhs } == BitFlags<Enum>{ rhs };
}

//! Return true if any flags of both enums (allowed to be bitflags) are
//! different
template<class Enum>
constexpr inline std::enable_if_t<is_bit_flag_v<Enum>, bool>
operator!=(Enum lhs, Enum rhs)
{
  return BitFlags<Enum>{ lhs } != BitFlags<Enum>{ rhs };
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_COMMON_BITFLAGS_HH
