#include "npms/core.hpp"
#include <chrono>
#include <iostream>
#include <numeric>
using namespace npms;
namespace {
Json surface_json(const Surface &s) {
  Json::Array rows;
  for (auto p : s.jet)
    rows.push_back(Json::Array{p[0], p[1], p[2]});
  return Json::Object{{"jet", rows},
                      {"normal", Json::Array{s.normal[0], s.normal[1], s.normal[2]}},
                      {"regular", s.regular},
                      {"mean", s.regular ? Json(s.mean) : Json()},
                      {"gaussian", s.regular ? Json(s.gaussian) : Json()},
                      {"kmin", s.regular ? Json(s.kmin) : Json()},
                      {"kmax", s.regular ? Json(s.kmax) : Json()},
                      {"partition", s.partition},
                      {"minimum", s.minimum}};
}
Json mesh_json(const Mesh &m, double ms) {
  Json::Array vs, ts, ps, mac;
  double residual = 0, minbasis = 1;
  std::size_t invalid = 0;
  for (auto &v : m.vertices) {
    auto &s = v.surface;
    vs.push_back(Json::Array{s.jet[0][0], s.jet[0][1], s.jet[0][2], s.normal[0], s.normal[1],
                             s.normal[2], v.influence, s.regular ? Json(s.mean) : Json(),
                             s.regular ? Json(s.gaussian) : Json(), s.partition - 1, s.minimum,
                             s.regular});
    residual = std::max(residual, std::abs(s.partition - 1));
    minbasis = std::min(minbasis, s.minimum);
    if (!s.regular)
      ++invalid;
  }
  for (auto t : m.triangles)
    ts.push_back(Json::Array{std::size_t(t[0]), std::size_t(t[1]), std::size_t(t[2])});
  for (auto e : m.ps_edges)
    ps.push_back(Json::Array{std::size_t(e[0]), std::size_t(e[1])});
  for (auto e : m.macro_edges)
    mac.push_back(Json::Array{std::size_t(e[0]), std::size_t(e[1])});
  return Json::Object{
      {"vertices", vs},
      {"triangles", ts},
      {"ps_edges", ps},
      {"macro_edges", mac},
      {"ms", ms},
      {"partition_residual", residual},
      {"minimum_basis", minbasis},
      {"invalid", invalid},
      {"buffer_bytes", m.vertices.size() * sizeof(Vertex) +
                           m.triangles.size() * sizeof(std::array<std::uint32_t, 3>)}};
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc < 4) {
      std::cerr << "Usage: npms_cli HANDOFF_DIR {5|9} "
                   "inventory|worker|sample|obj|benchmark [scene.json] "
                   "[resolution]\n";
      return 2;
    }
    std::filesystem::path root = argv[1];
    int n = std::stoi(argv[2]);
    if (n != 5 && n != 9)
      throw std::runtime_error("valence must be 5 or 9");
    auto f = Fixture::load(root / "fixtures" / ("fan" + std::to_string(n) + "_full42_demo_v1.bin"));
    std::string cmd = argv[3];
    if (cmd == "inventory") {
      Json::Array groups, supports;
      for (std::size_t i = 0; i < f.fine_count; ++i) {
        Json::Array support;
        for (std::size_t m = 0; m < f.valence; ++m) {
          bool active = false;
          for (std::size_t k = 0; k < 252; ++k)
            if (f.coefficients[(i * f.valence + m) * 252 + k] != 0)
              active = true;
          if (active)
            support.emplace_back(m);
        }
        supports.emplace_back(support);
      }
      for (std::size_t p = 0; p < f.coarse_count; ++p)
        groups.emplace_back(std::size_t(std::count(f.parents.begin(), f.parents.end(), p)));
      auto [lo, hi] = std::minmax_element(f.coefficients.begin(), f.coefficients.end());
      std::cout << Json(Json::Object{{"valence", n},
                                     {"controls", f.fine_count},
                                     {"parents", f.coarse_count},
                                     {"microcells", f.cells.size()},
                                     {"coeff_min", *lo},
                                     {"coeff_max", *hi},
                                     {"parent_group_sizes", groups},
                                     {"supports", supports},
                                     {"sha256", f.digest}})
                       .dump()
                << '\n';
      return 0;
    }
    auto defaultScene = Scene::from_json(
        Json::load((root / "scenes" / ("fan" + std::to_string(n) + "_rational.json")).string()), f);
    if (cmd == "worker") {
      // One JSON request per line; no filesystem paths or executable commands
      // accepted.
      for (;;) {
        std::string line;
        char c = 0;
        while (std::cin.get(c) && c != '\n') {
          if (line.size() >= 2 * 1024 * 1024)
            throw std::runtime_error("worker line exceeds 2 MiB");
          line += c;
        }
        if (line.empty() && !std::cin)
          break;
        try {
          auto j = Json::parse(line);
          auto s = j.contains("scene") ? Scene::from_json(j.at("scene"), f) : defaultScene;
          if (j.contains("xy")) {
            auto xy = j.at("xy");
            if (xy.size() != 2)
              throw std::runtime_error("expected xy pair");
            auto basis = f.basis({xy.at(0).number(), xy.at(1).number()});
            auto result = surface_json(evaluate_basis(basis, s));
            if (j.contains("basis")) {
              Json::Array rows;
              for (auto b : basis) {
                Json::Array row;
                for (double v : b)
                  row.emplace_back(v);
                rows.emplace_back(row);
              }
              std::get<Json::Object>(result.value).emplace("basis", rows);
            }
            std::cout << result.dump() << '\n';
          } else {
            unsigned res = static_cast<unsigned>(integer(j.at("resolution"), 32));
            auto selected = integer(j.at("selected"), f.fine_count - 1);
            auto start = std::chrono::steady_clock::now();
            auto m = tessellate(f, s, res, selected);
            double ms =
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                    .count();
            std::cout << mesh_json(m, ms).dump() << '\n';
          }
        } catch (const std::exception &e) {
          std::cout << Json(Json::Object{{"error", e.what()}}).dump() << '\n';
        }
        std::cout.flush();
      }
      return 0;
    }
    auto s = argc > 4 ? Scene::from_json(Json::load(argv[4]), f) : defaultScene;
    unsigned res = argc > 5 ? static_cast<unsigned>(std::stoul(argv[5])) : 8;
    if (cmd == "sample") {
      auto golden = Json::load(
          (root / "fixtures" / ("fan" + std::to_string(n) + "_golden_samples.json")).string());
      for (const auto &row : golden.at("samples").array()) {
        auto xy = row.at("xy");
        std::cout << surface_json(evaluate(f, s, {xy.at(0).number(), xy.at(1).number()})).dump()
                  << '\n';
      }
      return 0;
    }
    if (cmd == "obj") {
      auto m = tessellate(f, s, res);
      auto path = "fan" + std::to_string(n) + ".obj";
      export_obj(m, path);
      std::cout << path << ' ' << m.vertices.size() << " vertices " << m.triangles.size()
                << " triangles\n";
      return 0;
    }
    if (cmd == "benchmark") {
      Json::Array evals, meshes;
      double sink = 0;
      for (int count : {1000, 10000, 100000}) {
        std::vector<double> runs;
        for (int repeat = 0; repeat < 3; ++repeat) {
          auto start = std::chrono::steady_clock::now();
          for (int i = 0; i < count; ++i) {
            auto c = static_cast<std::size_t>(i) % f.cells.size();
            auto t = f.cells[c];
            Vec2 p{.2 * t[0][0] + .3 * t[1][0] + .5 * t[2][0],
                   .2 * t[0][1] + .3 * t[1][1] + .5 * t[2][1]};
            sink += evaluate(f, s, p).jet[0][2];
          }
          runs.push_back(
              std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                  .count());
        }
        std::sort(runs.begin(), runs.end());
        evals.push_back(Json::Object{{"count", count},
                                     {"min_ms", runs.front()},
                                     {"median_ms", runs[1]},
                                     {"max_ms", runs.back()}});
      }
      for (unsigned r : {4u, 8u, 16u, 32u}) {
        auto start = std::chrono::steady_clock::now();
        auto m = tessellate(f, s, r);
        double ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                .count();
        meshes.push_back(Json::Object{
            {"resolution", std::size_t(r)},
            {"vertices", m.vertices.size()},
            {"triangles", m.triangles.size()},
            {"ms", ms},
            {"buffer_bytes", m.vertices.size() * sizeof(Vertex) +
                                 m.triangles.size() * sizeof(std::array<std::uint32_t, 3>)}});
      }
      std::cout << Json(Json::Object{{"valence", n},
                                     {"evaluations", evals},
                                     {"tessellation", meshes},
                                     {"checksum", sink}})
                       .dump()
                << '\n';
      return 0;
    }
    throw std::runtime_error("unknown command");
  } catch (const std::exception &e) {
    std::cerr << "NPMS: " << e.what() << '\n';
    return 1;
  }
}
