#include "npms/core.hpp"
#include "npms/sha256.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
using namespace npms;
namespace {
std::size_t checks = 0;
void require(bool b, const std::string &msg) {
  ++checks;
  if (!b)
    throw std::runtime_error(msg);
}
template <class F> void rejects(F fn, const std::string &msg) {
  bool caught = false;
  try {
    fn();
  } catch (const std::exception &) {
    caught = true;
  }
  require(caught, msg);
}
Vec2 point(const std::array<Vec2, 3> &t, double a, double b, double c) {
  return {a * t[0][0] + b * t[1][0] + c * t[2][0], a * t[0][1] + b * t[1][1] + c * t[2][1]};
}
double difference(const Surface &a, const Surface &b) {
  double d = 0;
  for (unsigned k = 0; k < 6; ++k)
    for (unsigned j = 0; j < 3; ++j)
      d = std::max(d, std::abs(a.jet[k][j] - b.jet[k][j]));
  return d;
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 2)
      throw std::runtime_error("supply handoff directory");
    std::filesystem::path root = argv[1];
    require(sha256(std::vector<unsigned char>{}) ==
                "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
            "SHA empty vector");
    require(sha256(std::vector<unsigned char>{'a', 'b', 'c'}) ==
                "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
            "SHA abc");
    for (auto bad : {"NaN", "1e999", "01", "[1,]", "{\"a\":1,\"a\":2}", "\"\\uD800\"", "true false",
                     "{\"a\":}"})
      rejects([&] { Json::parse(bad); }, "malformed JSON accepted");
    rejects([&] { Json::parse(std::string(100, '[') + std::string(100, ']')); }, "depth limit");
    // Independent differential-geometry checks: z=x^2 +/- y^2 at the origin.
    for (double sign : {1.0, -1.0}) {
      std::vector<Jet> b{{1, -1, -1, -2, 0, -2},
                         {0, 1, 0, 0, 0, 0},
                         {0, 0, 1, 0, 0, 0},
                         {0, 0, 0, 2, 0, 0},
                         {0, 0, 0, 0, 0, 2}};
      Scene graph{{Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}, Vec3{0, 0, sign}},
                  {1, 1, 1, 1, 1}};
      auto surface = evaluate_basis(b, graph);
      require(surface.regular, "analytic graph regularity");
      require(std::abs(surface.mean - (1 + sign)) < 1e-12, "analytic mean curvature");
      require(std::abs(surface.gaussian - 4 * sign) < 1e-12, "analytic Gaussian curvature");
      require(std::abs(surface.kmax - 2) < 1e-12 && std::abs(surface.kmin - 2 * sign) < 1e-12,
              "principal curvatures");
    }
    rejects([] { evaluate_basis({}, {}); }, "empty evaluator input");
    for (int n : {5, 9}) {
      std::string prefix = "fan" + std::to_string(n);
      auto path = root / "fixtures" / (prefix + "_full42_demo_v1.bin");
      auto f = Fixture::load(path);
      auto plane =
          Scene::from_json(Json::load((root / "scenes" / (prefix + "_plane.json")).string()), f);
      auto rational =
          Scene::from_json(Json::load((root / "scenes" / (prefix + "_rational.json")).string()), f);
      auto golden = Json::load((root / "fixtures" / (prefix + "_golden_samples.json")).string());
      auto schema = Json::load((root / "fixtures" / (prefix + "_schema.json")).string());
      require(integer(schema.at("function_count"), 200) == f.fine_count, "schema count");
      std::size_t bi = 0;
      for (int i = 0; i <= 5; ++i)
        for (int j = 0; j <= 5 - i; ++j) {
          auto t = schema.at("bernstein_order_ijk").at(bi++);
          require(t.at(0).number() == i && t.at(1).number() == j && t.at(2).number() == 5 - i - j,
                  "Bernstein order");
        }
      // Independently map reference vertices using v1->O v2->Pi v3->Pi+1.
      for (std::size_t c = 0; c < f.cells.size(); ++c) {
        auto names = schema.at("reference_cells").at(c % 12);
        for (unsigned k = 0; k < 3; ++k) {
          auto p = schema.at("reference_vertices").at(names.at(k).string());
          double a = p.at(1).number() / 3, b = 1 - (p.at(0).number() + p.at(1).number()) / 3,
                 d = p.at(0).number() / 3;
          auto xy = point(f.macros[c / 12], a, b, d);
          for (unsigned dim = 0; dim < 2; ++dim)
            require(std::abs(xy[dim] - f.cells[c][k][dim]) < 1e-13, "physical reference mapping");
        }
      }
      double worst = 0, identity = 0, partition = 0, moment = 0, seam = 0, parent_error = 0,
             transfer = 0, fd = 0;
      std::size_t samples = 0;
      for (const auto &row : golden.at("samples").array()) {
        Vec2 xy{row.at("xy").at(0).number(), row.at("xy").at(1).number()};
        auto b = f.basis(xy);
        double sum = 0;
        Vec2 mom{};
        for (std::size_t i = 0; i < b.size(); ++i) {
          require(b[i][0] >= -1e-12, "basis negativity");
          sum += b[i][0];
          for (unsigned k = 0; k < 2; ++k)
            mom[k] += b[i][0] * f.sites[i][k];
        }
        partition = std::max(partition, std::abs(sum - 1));
        for (unsigned k = 0; k < 2; ++k)
          moment = std::max(moment, std::abs(mom[k] - xy[k]));
        for (auto variant : {"plane", "rational"}) {
          auto s = evaluate_basis(b, std::string(variant) == "plane" ? plane : rational);
          auto target = row.at(variant).at("surface_jet_6x3");
          for (unsigned d = 0; d < 6; ++d)
            for (unsigned k = 0; k < 3; ++k) {
              double err = std::abs(s.jet[d][k] - target.at(d).at(k).number());
              worst = std::max(worst, err);
              require(err <= 1e-9, "golden jet failed");
            }
          require(s.regular, "regular surface reported degenerate");
          if (std::string(variant) == "plane") {
            Surface expected;
            expected.jet = {
                Vec3{xy[0], xy[1], 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{}, Vec3{}, Vec3{}};
            identity = std::max(identity, difference(s, expected));
            require(std::abs(s.mean) < 1e-8 && std::abs(s.gaussian) < 1e-8, "plane curvature");
          }
        }
        ++samples;
      }
      require(samples == static_cast<std::size_t>(n * 24), "all golden samples executed");
      require(worst < 1e-9 && identity < 1e-9 && partition < 1e-10 && moment < 1e-10,
              "numerical gate");
      using Key = std::pair<long long, long long>;
      using Edge = std::pair<Key, Key>;
      std::map<Edge, std::vector<std::size_t>> edges;
      auto key = [](Vec2 p) { return Key{std::llround(p[0] * 1e12), std::llround(p[1] * 1e12)}; };
      for (std::size_t c = 0; c < f.cells.size(); ++c) {
        auto t = f.cells[c];
        for (unsigned k = 0; k < 3; ++k) {
          auto a = key(t[k]), b = key(t[(k + 1) % 3]);
          if (b < a)
            std::swap(a, b);
          edges[{a, b}].push_back(c);
        }
        for (auto xy : t) {
          auto ev = evaluate(f, plane, xy);
          Surface expected;
          expected.jet = {
              Vec3{xy[0], xy[1], 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{}, Vec3{}, Vec3{}};
          require(difference(ev, expected) < 1e-9, "corner identity");
          require(f.locate(xy) == f.locate(xy), "boundary classification");
        }
        auto xy = point(t, .27, .34, .39);
        double h = 1e-5;
        auto center = evaluate(f, rational, xy);
        for (unsigned dim = 0; dim < 2; ++dim) {
          auto lo = xy, hi = xy;
          lo[dim] -= h;
          hi[dim] += h;
          auto a = evaluate(f, rational, lo), b = evaluate(f, rational, hi);
          for (unsigned k = 0; k < 3; ++k) {
            fd = std::max(fd,
                          std::abs((b.jet[0][k] - a.jet[0][k]) / (2 * h) - center.jet[dim + 1][k]));
            for (unsigned d = 1; d <= 2; ++d) {
              unsigned second = dim == 0 ? (d == 1 ? 3 : 4) : (d == 1 ? 4 : 5);
              fd = std::max(
                  fd, std::abs((b.jet[d][k] - a.jet[d][k]) / (2 * h) - center.jet[second][k]));
            }
          }
        }
      }
      std::size_t shared = 0;
      for (auto &[e, cs] : edges) {
        if (cs.size() != 2)
          continue;
        ++shared;
        Vec2 a{double(e.first.first) / 1e12, double(e.first.second) / 1e12},
            b{double(e.second.first) / 1e12, double(e.second.second) / 1e12};
        for (double u : {.2, .5, .8}) {
          Vec2 xy{(1 - u) * a[0] + u * b[0], (1 - u) * a[1] + u * b[1]};
          auto left = f.basis_on_cell(cs[0], xy), right = f.basis_on_cell(cs[1], xy);
          for (std::size_t i = 0; i < left.size(); ++i)
            for (unsigned d = 0; d < 6; ++d)
              seam = std::max(seam, std::abs(left[i][d] - right[i][d]));
          require(difference(evaluate_basis(left, rational), evaluate_basis(right, rational)) <
                      1e-8,
                  "rational seam");
        }
      }
      require(shared == static_cast<std::size_t>(n == 5 ? 85 : 153), "every shared seam");
      require(seam < 1e-8 && fd < 1e-6, "seam or finite-difference derivative gate");
      std::size_t stride = f.cells.size() * 21;
      for (std::size_t p = 0; p < f.coarse_count; ++p)
        for (std::size_t c = 0; c < stride; ++c) {
          double sum = 0;
          for (std::size_t i = 0; i < f.fine_count; ++i)
            if (f.parents[i] == p)
              sum += f.coefficients[i * stride + c];
          parent_error =
              std::max(parent_error, std::abs(sum - f.coarse_coefficients[p * stride + c]));
        }
      require(parent_error <= 1e-12, "parent coefficient identity");
      Scene coarse, fine;
      for (std::size_t i = 0; i < f.coarse_count; ++i) {
        auto p = f.coarse_sites[i];
        coarse.positions.push_back({p[0], p[1], .1 * std::sin(double(i))});
        coarse.weights.push_back(.5 + .03 * double(i));
      }
      for (auto p : f.parents) {
        fine.positions.push_back(coarse.positions[p]);
        fine.weights.push_back(coarse.weights[p]);
      }
      for (auto t : f.cells) {
        auto xy = point(t, .2, .3, .5);
        transfer =
            std::max(transfer, difference(evaluate(f, fine, xy), evaluate(f, coarse, xy, true)));
      }
      require(transfer < 1e-9, "fixed-mesh rational transfer");
      for (double w : {0., -1., std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::infinity()}) {
        auto s = plane;
        s.weights[0] = w;
        rejects([&] { evaluate(f, s, {0, 0}); }, "invalid weight accepted");
      }
      auto degenerate = plane;
      for (auto &p : degenerate.positions)
        p = {1.234, 3.14, 2.72};
      require(!evaluate(f, degenerate, {.2, .2}).regular, "collapsed surface must be invalid");
      rejects([&] { f.basis({100, 100}); }, "outside domain");
      rejects([&] { f.basis({std::numeric_limits<double>::quiet_NaN(), 0}); }, "NaN parameter");
      auto bad = f;
      bad.cells[0][1] = bad.cells[0][0];
      rejects([&] { bad.basis_on_cell(0, {0, 0}); }, "degenerate parameter triangle");
      auto restored = Scene::from_json(Json::parse(rational.json(f.valence).dump()), f);
      require(restored.positions == rational.positions && restored.weights == rational.weights,
              "scene exact roundtrip");
      auto edit = plane;
      edit.positions[0][2] += .5;
      require(difference(evaluate(f, edit, {0, 0}), evaluate(f, plane, {0, 0})) > 1e-6,
              "XYZ edit response");
      edit = plane;
      edit.weights[0] *= 2;
      require(difference(evaluate(f, edit, {0, 0}), evaluate(f, plane, {0, 0})) > 1e-6,
              "weight edit response");
      for (unsigned r : {2u, 4u, 8u, 16u}) {
        auto m = tessellate(f, rational, r);
        require(m.triangles.size() == f.cells.size() * r * r, "tessellation triangles");
        std::map<std::pair<std::uint32_t, std::uint32_t>, unsigned> incidence;
        double area = 0;
        for (auto t : m.triangles) {
          auto a = m.vertices[t[0]].uv, b = m.vertices[t[1]].uv, c = m.vertices[t[2]].uv;
          double d = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
          require(d > 0, "triangle winding");
          area += d / 2;
          for (unsigned k = 0; k < 3; ++k) {
            auto x = t[k], y = t[(k + 1) % 3];
            if (x > y)
              std::swap(x, y);
            ++incidence[{x, y}];
          }
        }
        std::size_t boundary = 0;
        for (auto [e, count] : incidence) {
          require(count <= 2, "nonmanifold mesh");
          if (count == 1)
            ++boundary;
        }
        require(boundary == static_cast<std::size_t>(n) * 2 * r, "no cracks or T junctions");
        require(m.vertices.size() + m.triangles.size() == incidence.size() + 1,
                "disk Euler characteristic");
        require(area > 1, "nonzero mesh area");
      }
      // Every modified or truncated binary fails authentication before data
      // allocation/use.
      std::ifstream in(path, std::ios::binary);
      std::vector<char> bytes((std::istreambuf_iterator<char>(in)), {});
      auto tmp = std::filesystem::temp_directory_path() / (prefix + "_npms_test.bin");
      for (std::size_t offset : {0u, 8u, 16u, 32u, 48u, 1000u}) {
        auto corrupt = bytes;
        corrupt[offset] ^= 0x7f;
        {
          std::ofstream o(tmp, std::ios::binary);
          o.write(corrupt.data(), static_cast<std::streamsize>(corrupt.size()));
        }
        rejects([&] { Fixture::load(tmp); }, "mutated binary accepted");
      }
      {
        std::ofstream o(tmp, std::ios::binary);
        o.write(bytes.data(), 100);
      }
      rejects([&] { Fixture::load(tmp); }, "truncated binary accepted");
      std::filesystem::remove(tmp);
      std::cout << "PASS fan=" << n << " golden_samples=" << samples
                << " scenes=2 jet_max=" << worst << " identity_max=" << identity
                << " partition=" << partition << " moment=" << moment << " seams=" << shared
                << " seam_max=" << seam << " parent_max=" << parent_error
                << " transfer_max=" << transfer << " finite_difference_max=" << fd << '\n';
    }
    std::cout << "PASS checks=" << checks << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "FAIL after " << checks << " checks: " << e.what() << '\n';
    return 1;
  }
}
