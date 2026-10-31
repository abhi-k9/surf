#include "npms/core.hpp"
#include "npms/sha256.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <tuple>
namespace npms {
namespace {
constexpr std::array<unsigned, 6> fact = {1, 1, 2, 6, 24, 120};
Jet mul(const Jet &a, const Jet &b) {
  return {a[0] * b[0],
          a[1] * b[0] + a[0] * b[1],
          a[2] * b[0] + a[0] * b[2],
          a[3] * b[0] + 2 * a[1] * b[1] + a[0] * b[3],
          a[4] * b[0] + a[1] * b[2] + a[2] * b[1] + a[0] * b[4],
          a[5] * b[0] + 2 * a[2] * b[2] + a[0] * b[5]};
}
double dot(Vec3 a, Vec3 b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(Vec3 a, Vec3 b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
std::array<Jet, 21> polynomials(const Fixture &f, std::size_t cell, Vec2 xy) {
  auto lam = f.barycentric(cell, xy);
  auto t = f.cells.at(cell);
  double ax = t[1][0] - t[0][0], ay = t[1][1] - t[0][1], bx = t[2][0] - t[0][0],
         by = t[2][1] - t[0][1], det = ax * by - ay * bx;
  Jet u{lam[1], by / det, -bx / det, 0, 0, 0}, v{lam[2], -ay / det, ax / det, 0, 0, 0},
      w{lam[0], -u[1] - v[1], -u[2] - v[2], 0, 0, 0};
  std::array<std::array<Jet, 6>, 3> powers{};
  std::array<Jet, 3> coords{w, u, v};
  for (unsigned a = 0; a < 3; ++a) {
    powers[a][0] = {1, 0, 0, 0, 0, 0};
    for (unsigned i = 1; i <= 5; ++i)
      powers[a][i] = mul(powers[a][i - 1], coords[a]);
  }
  std::array<Jet, 21> out{};
  unsigned n = 0;
  for (unsigned i = 0; i <= 5; ++i)
    for (unsigned j = 0; j <= 5 - i; ++j) {
      auto k = 5 - i - j;
      auto q = mul(mul(powers[0][i], powers[1][j]), powers[2][k]);
      double scale = 120.0 / double(fact[i] * fact[j] * fact[k]);
      for (auto &z : q)
        z *= scale;
      out[n++] = q;
    }
  return out;
}
} // namespace
Fixture Fixture::load(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary | std::ios::ate);
  if (!in)
    throw std::runtime_error("cannot open fixture");
  auto size = in.tellg();
  if (size != 1194888 && size != 3674664)
    throw std::runtime_error("fixture byte length is not a certified layout");
  std::vector<unsigned char> bytes(static_cast<std::size_t>(size));
  in.seekg(0);
  if (!in.read(reinterpret_cast<char *>(bytes.data()), size))
    throw std::runtime_error("fixture read failed");
  if (std::string(bytes.begin(), bytes.begin() + 8) != "NPMSFAN1")
    throw std::runtime_error("bad fixture magic");
  std::size_t offset = 8;
  auto u32 = [&]() {
    std::uint32_t v = 0;
    for (unsigned j = 0; j < 4; ++j)
      v |= std::uint32_t(bytes.at(offset++)) << (8 * j);
    return v;
  };
  auto version = u32(), degree = u32(), n = u32(), macros = u32(), cells = u32(), bern = u32(),
       fine = u32(), coarse = u32();
  if (version != 1 || degree != 5 || (n != 5 && n != 9) || macros != n || cells != 12 ||
      bern != 21 || fine != 18 * n + 10 || coarse != 3 * (n + 1))
    throw std::runtime_error("unsupported fixture header/counts");
  if (static_cast<std::size_t>(size) != (n == 5 ? 1194888u : 3674664u))
    throw std::runtime_error("header/length mismatch");
  Fixture f;
  f.digest = sha256(bytes);
  const std::string expected =
      n == 5 ? "3831b9c5aa314da4391522ea99e46ff340fe5ddaa00326ad9db39f53884ee68d"
             : "a3db27d57bff39c4a2e49c268eb50ca91981bc1b4c6247876d3b444a38faad18";
  if (f.digest != expected)
    throw std::runtime_error("fixture SHA-256 mismatch");
  f.valence = n;
  f.fine_count = fine;
  f.coarse_count = coarse;
  auto real = [&]() {
    std::uint64_t bits = 0;
    for (unsigned j = 0; j < 8; ++j)
      bits |= std::uint64_t(bytes.at(offset++)) << (8 * j);
    double v = std::bit_cast<double>(bits);
    if (!std::isfinite(v))
      throw std::runtime_error("nonfinite fixture");
    return v;
  };
  f.macros.resize(n);
  f.cells.resize(n * 12);
  f.sites.resize(fine);
  f.parents.resize(fine);
  f.coarse_sites.resize(coarse);
  for (auto &t : f.macros)
    for (auto &p : t)
      for (auto &x : p)
        x = real();
  for (auto &t : f.cells)
    for (auto &p : t)
      for (auto &x : p)
        x = real();
  for (auto &p : f.sites)
    for (auto &x : p)
      x = real();
  for (auto &p : f.parents) {
    p = u32();
    if (p >= coarse)
      throw std::runtime_error("bad parent index");
  }
  f.coefficients.resize(fine * n * 12 * 21);
  for (auto &x : f.coefficients) {
    x = real();
    if (x < 0)
      throw std::runtime_error("negative certified coefficient");
  }
  for (auto &p : f.coarse_sites)
    for (auto &x : p)
      x = real();
  f.coarse_coefficients.resize(coarse * n * 12 * 21);
  for (auto &x : f.coarse_coefficients)
    x = real();
  if (offset != bytes.size())
    throw std::runtime_error("trailing binary payload");
  for (const auto &t : f.cells) {
    double det =
        (t[1][0] - t[0][0]) * (t[2][1] - t[0][1]) - (t[1][1] - t[0][1]) * (t[2][0] - t[0][0]);
    if (!std::isfinite(det) || std::abs(det) < 1e-14)
      throw std::runtime_error("degenerate parameter triangle");
  }
  return f;
}
std::array<double, 3> Fixture::barycentric(std::size_t cell, Vec2 xy) const {
  if (!std::isfinite(xy[0]) || !std::isfinite(xy[1]))
    throw std::runtime_error("nonfinite parameter point");
  auto t = cells.at(cell);
  double ax = t[1][0] - t[0][0], ay = t[1][1] - t[0][1], bx = t[2][0] - t[0][0],
         by = t[2][1] - t[0][1], det = ax * by - ay * bx;
  if (!std::isfinite(det) || std::abs(det) < 1e-14)
    throw std::runtime_error("degenerate parameter triangle");
  double x = xy[0] - t[0][0], y = xy[1] - t[0][1], u = (by * x - bx * y) / det,
         v = (-ay * x + ax * y) / det;
  return {1 - u - v, u, v};
}
std::size_t Fixture::locate(Vec2 p) const {
  for (std::size_t c = 0; c < cells.size(); ++c) {
    auto l = barycentric(c, p);
    if (*std::min_element(l.begin(), l.end()) >= -1e-12)
      return c;
  }
  throw std::runtime_error("OutsideDomain");
}
std::vector<Jet> Fixture::basis_on_cell(std::size_t cell, Vec2 p, bool coarse) const {
  auto q = polynomials(*this, cell, p);
  std::size_t count = coarse ? coarse_count : fine_count;
  const auto &coeff = coarse ? coarse_coefficients : coefficients;
  std::vector<Jet> out(count);
  for (std::size_t i = 0; i < count; ++i) {
    auto base = (i * cells.size() + cell) * 21;
    for (std::size_t k = 0; k < 21; ++k) {
      double c = coeff.at(base + k);
      if (c != 0)
        for (std::size_t d = 0; d < 6; ++d)
          out[i][d] += c * q[k][d];
    }
  }
  return out;
}
void Scene::validate(std::size_t n) const {
  if (n == 0 || positions.size() != n || weights.size() != n)
    throw std::runtime_error("scene control count mismatch");
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(weights[i]) || weights[i] <= 0)
      throw std::runtime_error("weights must be finite and strictly positive");
    for (double x : positions[i])
      if (!std::isfinite(x))
        throw std::runtime_error("control coordinates must be finite");
  }
}
Scene Scene::from_json(const Json &j, const Fixture &f) {
  if (j.at("schema").string() != "npms-fixed-fan-scene-v1" ||
      integer(j.at("valence"), 9) != f.valence)
    throw std::runtime_error("incompatible scene schema or valence");
  Scene s;
  if (j.at("control_xyz").size() != f.fine_count || j.at("weights").size() != f.fine_count)
    throw std::runtime_error("scene count mismatch");
  for (auto &p : j.at("control_xyz").array()) {
    if (p.size() != 3)
      throw std::runtime_error("control must have three coordinates");
    s.positions.push_back({p.at(0).number(), p.at(1).number(), p.at(2).number()});
  }
  for (auto &w : j.at("weights").array())
    s.weights.push_back(w.number());
  s.validate(f.fine_count);
  return s;
}
Json Scene::json(std::size_t n) const {
  Json::Array ps, ws;
  for (auto p : positions)
    ps.push_back(Json::Array{p[0], p[1], p[2]});
  for (double w : weights)
    ws.emplace_back(w);
  return Json::Object{
      {"schema", "npms-fixed-fan-scene-v1"}, {"valence", n}, {"control_xyz", ps}, {"weights", ws}};
}
Surface evaluate_basis(const std::vector<Jet> &b, const Scene &s) {
  s.validate(b.size());
  Surface r;
  r.minimum = std::numeric_limits<double>::max();
  std::array<Vec3, 6> a{};
  Jet w{};
  // Center homogeneous coordinates to avoid cancellation under translations.
  // Scale weights to preserve projective invariance even near DBL_MAX.
  const double scale = *std::max_element(s.weights.begin(), s.weights.end());
  for (std::size_t i = 0; i < b.size(); ++i) {
    r.partition += b[i][0];
    r.minimum = std::min(r.minimum, b[i][0]);
    for (std::size_t d = 0; d < 6; ++d) {
      double q = b[i][d] * (s.weights[i] / scale);
      w[d] += q;
      for (std::size_t k = 0; k < 3; ++k)
        a[d][k] += q * (s.positions[i][k] - s.positions[0][k]);
    }
  }
  if (!std::isfinite(w[0]) || w[0] <= 0)
    throw std::runtime_error("invalid rational denominator");
  r.denominator = w[0] * scale;
  for (std::size_t k = 0; k < 3; ++k) {
    auto &t = r.jet;
    t[0][k] = a[0][k] / w[0];
    t[1][k] = (a[1][k] - t[0][k] * w[1]) / w[0];
    t[2][k] = (a[2][k] - t[0][k] * w[2]) / w[0];
    t[3][k] = (a[3][k] - t[0][k] * w[3] - 2 * t[1][k] * w[1]) / w[0];
    t[4][k] = (a[4][k] - t[0][k] * w[4] - t[1][k] * w[2] - t[2][k] * w[1]) / w[0];
    t[5][k] = (a[5][k] - t[0][k] * w[5] - 2 * t[2][k] * w[2]) / w[0];
  }
  // Translate only the value back; derivatives of a constant translation vanish.
  for (std::size_t k = 0; k < 3; ++k)
    r.jet[0][k] += s.positions[0][k];
  for (auto row : r.jet)
    for (double x : row)
      if (!std::isfinite(x))
        throw std::runtime_error("surface arithmetic overflow");
  auto n = cross(r.jet[1], r.jet[2]);
  double E = dot(r.jet[1], r.jet[1]), F = dot(r.jet[1], r.jet[2]), G = dot(r.jet[2], r.jet[2]),
         det = dot(n, n);
  if (!std::isfinite(det) || !std::isfinite(E * G) || E * G <= 0 || det <= 1e-24 * E * G)
    return r;
  double length = std::sqrt(det);
  for (unsigned k = 0; k < 3; ++k)
    r.normal[k] = n[k] / length;
  double e = dot(r.normal, r.jet[3]), f = dot(r.normal, r.jet[4]), g = dot(r.normal, r.jet[5]);
  r.mean = (e * G - 2 * f * F + g * E) / (2 * det);
  r.gaussian = (e * g - f * f) / det;
  double disc = r.mean * r.mean - r.gaussian,
         tol = 1e-12 * std::max({1.0, r.mean * r.mean, std::abs(r.gaussian)});
  if (!std::isfinite(r.mean) || !std::isfinite(r.gaussian) || !std::isfinite(disc) || disc < -tol) {
    r.mean = 0;
    r.gaussian = 0;
    r.normal = {};
    return r;
  }
  // Only a roundoff-sized negative discriminant is replaced by zero.
  r.kmin = r.mean - std::sqrt(std::max(0.0, disc));
  r.kmax = r.mean + std::sqrt(std::max(0.0, disc));
  r.regular = true;
  return r;
}
Surface evaluate(const Fixture &f, const Scene &s, Vec2 xy, bool coarse) {
  return evaluate_basis(f.basis(xy, coarse), s);
}
Mesh tessellate(const Fixture &f, const Scene &s, unsigned res, std::size_t selected) {
  if (res < 1 || res > 32 || selected >= f.fine_count)
    throw std::runtime_error("resolution must be 1..32; selected control must exist");
  s.validate(f.fine_count);
  Mesh mesh;
  std::map<std::pair<long long, long long>, std::uint32_t> ids;
  std::set<std::array<std::uint32_t, 2>> ps, macro;
  auto edge = [](std::uint32_t a, std::uint32_t b) {
    return std::array<std::uint32_t, 2>{std::min(a, b), std::max(a, b)};
  };
  for (std::size_t c = 0; c < f.cells.size(); ++c) {
    auto t = f.cells[c];
    std::vector<std::vector<std::uint32_t>> grid(res + 1);
    for (unsigned i = 0; i <= res; ++i)
      for (unsigned j = 0; j <= res - i; ++j) {
        double u = double(i) / res, v = double(j) / res;
        Vec2 xy{(1 - u - v) * t[0][0] + u * t[1][0] + v * t[2][0],
                (1 - u - v) * t[0][1] + u * t[1][1] + v * t[2][1]};
        auto key = std::make_pair(std::llround(xy[0] * 1e12), std::llround(xy[1] * 1e12));
        auto [it, inserted] = ids.emplace(key, static_cast<std::uint32_t>(mesh.vertices.size()));
        if (inserted) {
          auto b = f.basis_on_cell(c, xy);
          auto surf = evaluate_basis(b, s);
          double maxw = *std::max_element(s.weights.begin(), s.weights.end()), den = 0;
          for (std::size_t k = 0; k < b.size(); ++k)
            den += b[k][0] * (s.weights[k] / maxw);
          mesh.vertices.push_back({xy, surf, b[selected][0] * (s.weights[selected] / maxw) / den});
        }
        grid[i].push_back(it->second);
      }
    auto add = [&](std::uint32_t a, std::uint32_t b, std::uint32_t d) {
      auto x = mesh.vertices[a].uv, y = mesh.vertices[b].uv, z = mesh.vertices[d].uv;
      if ((y[0] - x[0]) * (z[1] - x[1]) - (y[1] - x[1]) * (z[0] - x[0]) < 0)
        std::swap(b, d);
      mesh.triangles.push_back({a, b, d});
    };
    for (unsigned i = 0; i < res; ++i)
      for (unsigned j = 0; j < res - i; ++j) {
        add(grid[i][j], grid[i + 1][j], grid[i][j + 1]);
        if (i + j + 1 < res)
          add(grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1]);
      }
    auto mark = [&](std::uint32_t a, std::uint32_t b) {
      auto e = edge(a, b);
      ps.insert(e);
      auto x = mesh.vertices[a].uv, y = mesh.vertices[b].uv;
      auto mt = f.macros[c / 12];
      for (unsigned k = 0; k < 3; ++k) {
        auto v0 = mt[k], v1 = mt[(k + 1) % 3];
        auto line = [&](Vec2 p) {
          return std::abs((v1[0] - v0[0]) * (p[1] - v0[1]) - (v1[1] - v0[1]) * (p[0] - v0[0])) <
                 1e-10;
        };
        if (line(x) && line(y))
          macro.insert(e);
      }
    };
    for (unsigned k = 0; k < res; ++k) {
      mark(grid[0][k], grid[0][k + 1]);
      mark(grid[k][0], grid[k + 1][0]);
      mark(grid[k][res - k], grid[k + 1][res - k - 1]);
    }
  }
  mesh.ps_edges.assign(ps.begin(), ps.end());
  mesh.macro_edges.assign(macro.begin(), macro.end());
  return mesh;
}
void export_obj(const Mesh &m, const std::filesystem::path &path) {
  std::ofstream o(path);
  if (!o)
    throw std::runtime_error("cannot write OBJ");
  o << std::setprecision(17)
    << "# NPMS fixed-fan demo; analytic normals. Invalid normals are zero.\n";
  for (auto &v : m.vertices)
    o << "v " << v.surface.jet[0][0] << ' ' << v.surface.jet[0][1] << ' ' << v.surface.jet[0][2]
      << '\n';
  for (auto &v : m.vertices)
    o << "vn " << v.surface.normal[0] << ' ' << v.surface.normal[1] << ' ' << v.surface.normal[2]
      << '\n';
  for (auto t : m.triangles) {
    o << 'f';
    for (auto i : t)
      o << ' ' << i + 1 << "//" << i + 1;
    o << '\n';
  }
  if (!o)
    throw std::runtime_error("OBJ write failed");
}
} // namespace npms
