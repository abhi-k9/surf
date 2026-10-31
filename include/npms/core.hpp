#pragma once
#include "json.hpp"
#include <array>
#include <cstdint>
#include <filesystem>
#include <vector>
namespace npms {
using Vec2 = std::array<double, 2>;
using Vec3 = std::array<double, 3>;
using Jet = std::array<double, 6>;
struct Fixture {
  std::size_t valence{}, fine_count{}, coarse_count{};
  std::vector<std::array<Vec2, 3>> macros, cells;
  std::vector<Vec2> sites, coarse_sites;
  std::vector<std::uint32_t> parents;
  std::vector<double> coefficients, coarse_coefficients;
  std::string digest;
  static Fixture load(const std::filesystem::path &);
  std::array<double, 3> barycentric(std::size_t, Vec2) const;
  std::size_t locate(Vec2) const;
  std::vector<Jet> basis_on_cell(std::size_t, Vec2, bool coarse = false) const;
  std::vector<Jet> basis(Vec2 p, bool coarse = false) const {
    return basis_on_cell(locate(p), p, coarse);
  }
};
struct Scene {
  std::vector<Vec3> positions;
  std::vector<double> weights;
  void validate(std::size_t) const;
  static Scene from_json(const Json &, const Fixture &);
  Json json(std::size_t valence) const;
};
struct Surface {
  std::array<Vec3, 6> jet{};
  Vec3 normal{};
  bool regular = false;
  double mean = 0, gaussian = 0, kmin = 0, kmax = 0, denominator = 0, partition = 0, minimum = 0;
};
Surface evaluate_basis(const std::vector<Jet> &, const Scene &);
Surface evaluate(const Fixture &, const Scene &, Vec2, bool coarse = false);
struct Vertex {
  Vec2 uv;
  Surface surface;
  double influence;
};
struct Mesh {
  std::vector<Vertex> vertices;
  std::vector<std::array<std::uint32_t, 3>> triangles;
  std::vector<std::array<std::uint32_t, 2>> ps_edges, macro_edges;
};
Mesh tessellate(const Fixture &, const Scene &, unsigned resolution, std::size_t selected = 0);
void export_obj(const Mesh &, const std::filesystem::path &);
} // namespace npms
