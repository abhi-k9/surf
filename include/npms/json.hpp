#pragma once
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>
namespace npms {
// Bounded, strict JSON. Duplicate keys, nonfinite numbers and malformed escapes
// fail.
class Json {
public:
  using Array = std::vector<Json>;
  using Object = std::map<std::string, Json>;
  std::variant<std::nullptr_t, bool, double, std::string, Array, Object> value = nullptr;
  Json() = default;
  Json(std::nullptr_t) : value(nullptr) {}
  Json(bool v) : value(v) {}
  Json(double v) : value(v) {
    if (!std::isfinite(v))
      throw std::runtime_error("nonfinite JSON number");
  }
  Json(int v) : Json(double(v)) {}
  Json(std::size_t v) : Json(double(v)) {}
  Json(const char *v) : value(std::string(v)) {}
  Json(std::string v) : value(std::move(v)) {}
  Json(Array v) : value(std::move(v)) {}
  Json(Object v) : value(std::move(v)) {}
  double number() const { return std::get<double>(value); }
  const std::string &string() const { return std::get<std::string>(value); }
  const Array &array() const { return std::get<Array>(value); }
  const Object &object() const { return std::get<Object>(value); }
  const Json &at(const std::string &k) const { return object().at(k); }
  const Json &at(std::size_t i) const { return array().at(i); }
  bool contains(const std::string &k) const {
    return std::holds_alternative<Object>(value) && object().contains(k);
  }
  std::size_t size() const { return array().size(); }
  static Json parse(const std::string &s) {
    if (s.size() > 2 * 1024 * 1024)
      throw std::runtime_error("JSON exceeds 2 MiB");
    struct Parser {
      const std::string &s;
      std::size_t p = 0;
      [[noreturn]] void fail() const {
        throw std::runtime_error("invalid JSON at byte " + std::to_string(p));
      }
      void ws() {
        while (p < s.size() && (s[p] == ' ' || s[p] == '\r' || s[p] == '\n' || s[p] == '\t'))
          ++p;
      }
      bool eat(char c) {
        ws();
        if (p < s.size() && s[p] == c) {
          ++p;
          return true;
        }
        return false;
      }
      unsigned hex() {
        unsigned v = 0;
        for (int i = 0; i < 4; ++i) {
          if (p == s.size())
            fail();
          char c = s[p++];
          v <<= 4;
          if (c >= '0' && c <= '9')
            v += c - '0';
          else if (c >= 'a' && c <= 'f')
            v += c - 'a' + 10;
          else if (c >= 'A' && c <= 'F')
            v += c - 'A' + 10;
          else
            fail();
        }
        return v;
      }
      std::string str() {
        if (!eat('"'))
          fail();
        std::string o;
        while (p < s.size()) {
          unsigned char c = static_cast<unsigned char>(s[p++]);
          if (c == '"')
            return o;
          if (c < 32)
            fail();
          if (c != '\\') {
            o += char(c);
            continue;
          }
          if (p == s.size())
            fail();
          char e = s[p++];
          switch (e) {
          case '"':
          case '\\':
          case '/':
            o += e;
            break;
          case 'b':
            o += '\b';
            break;
          case 'f':
            o += '\f';
            break;
          case 'n':
            o += '\n';
            break;
          case 'r':
            o += '\r';
            break;
          case 't':
            o += '\t';
            break;
          case 'u': {
            unsigned u = hex();
            if (u >= 0xd800 && u <= 0xdbff) {
              if (p + 2 > s.size() || s[p++] != '\\' || s[p++] != 'u')
                fail();
              unsigned lo = hex();
              if (lo < 0xdc00 || lo > 0xdfff)
                fail();
              u = 0x10000 + ((u - 0xd800) << 10) + (lo - 0xdc00);
            } else if (u >= 0xdc00 && u <= 0xdfff)
              fail();
            if (u < 128)
              o += char(u);
            else if (u < 2048) {
              o += char(192 | (u >> 6));
              o += char(128 | (u & 63));
            } else if (u < 65536) {
              o += char(224 | (u >> 12));
              o += char(128 | ((u >> 6) & 63));
              o += char(128 | (u & 63));
            } else {
              o += char(240 | (u >> 18));
              o += char(128 | ((u >> 12) & 63));
              o += char(128 | ((u >> 6) & 63));
              o += char(128 | (u & 63));
            }
            break;
          }
          default:
            fail();
          }
        }
        fail();
      }
      Json read(unsigned d) {
        if (d > 32)
          fail();
        ws();
        if (p == s.size())
          fail();
        char c = s[p];
        if (c == '"')
          return str();
        if (eat('{')) {
          Object o;
          if (eat('}'))
            return o;
          do {
            auto k = str();
            if (!eat(':'))
              fail();
            if (!o.emplace(k, read(d + 1)).second)
              fail();
          } while (eat(','));
          if (!eat('}'))
            fail();
          return o;
        }
        if (eat('[')) {
          Array a;
          if (eat(']'))
            return a;
          do {
            a.push_back(read(d + 1));
            if (a.size() > 100000)
              fail();
          } while (eat(','));
          if (!eat(']'))
            fail();
          return a;
        }
        for (auto lit : {std::string("null"), std::string("true"), std::string("false")})
          if (s.compare(p, lit.size(), lit) == 0) {
            p += lit.size();
            if (lit == "null")
              return nullptr;
            return lit == "true";
          }
        const std::size_t start = p;
        if (s[p] == '-')
          ++p;
        if (p == s.size())
          fail();
        if (s[p] == '0')
          ++p;
        else {
          if (s[p] < '1' || s[p] > '9')
            fail();
          while (p < s.size() && s[p] >= '0' && s[p] <= '9')
            ++p;
        }
        if (p < s.size() && s[p] == '.') {
          ++p;
          auto b = p;
          while (p < s.size() && s[p] >= '0' && s[p] <= '9')
            ++p;
          if (b == p)
            fail();
        }
        if (p < s.size() && (s[p] == 'e' || s[p] == 'E')) {
          ++p;
          if (p < s.size() && (s[p] == '+' || s[p] == '-'))
            ++p;
          auto b = p;
          while (p < s.size() && s[p] >= '0' && s[p] <= '9')
            ++p;
          if (b == p)
            fail();
        }
        double v = 0;
        auto r = std::from_chars(s.data() + start, s.data() + p, v);
        if (r.ec != std::errc{} || r.ptr != s.data() + p || !std::isfinite(v))
          fail();
        return v;
      }
    } p{s};
    auto v = p.read(0);
    p.ws();
    if (p.p != s.size())
      p.fail();
    return v;
  }
  static Json load(const std::string &path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f || f.tellg() < 0 || f.tellg() > 2 * 1024 * 1024)
      throw std::runtime_error("cannot read bounded JSON: " + path);
    std::string s(static_cast<std::size_t>(f.tellg()), '\0');
    f.seekg(0);
    if (!f.read(s.data(), static_cast<std::streamsize>(s.size())))
      throw std::runtime_error("JSON read failed");
    return parse(s);
  }
  std::string dump() const {
    std::ostringstream o;
    o << std::setprecision(17);
    write(o);
    return o.str();
  }
  void write(std::ostream &o) const {
    if (std::holds_alternative<std::nullptr_t>(value)) {
      o << "null";
      return;
    }
    if (auto v = std::get_if<bool>(&value)) {
      o << (*v ? "true" : "false");
      return;
    }
    if (auto v = std::get_if<double>(&value)) {
      o << *v;
      return;
    }
    if (auto v = std::get_if<std::string>(&value)) {
      o << '"';
      for (unsigned char c : *v) {
        if (c == '"' || c == '\\')
          o << '\\' << c;
        else if (c < 32)
          o << "\\u00" << "0123456789abcdef"[c >> 4] << "0123456789abcdef"[c & 15];
        else
          o << c;
      }
      o << '"';
      return;
    }
    if (auto v = std::get_if<Array>(&value)) {
      o << '[';
      bool first = true;
      for (auto &x : *v) {
        if (!first)
          o << ',';
        first = false;
        x.write(o);
      }
      o << ']';
      return;
    }
    o << '{';
    bool first = true;
    for (auto &[k, v] : object()) {
      if (!first)
        o << ',';
      first = false;
      Json(k).write(o);
      o << ':';
      v.write(o);
    }
    o << '}';
  }
};
inline std::size_t integer(const Json &j, std::size_t max) {
  double v = j.number();
  if (v < 0 || v > double(max) || std::floor(v) != v)
    throw std::runtime_error("integer outside allowed range");
  return static_cast<std::size_t>(v);
}
} // namespace npms
