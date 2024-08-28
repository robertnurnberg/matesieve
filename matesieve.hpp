#include <zlib.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "external/json.hpp"

struct TestMetaData {
  std::optional<std::string> book, new_tc, resolved_base, resolved_new, tc;
  std::optional<int> threads;
  std::optional<bool> sprt;
  std::optional<std::vector<int>> pentanomial;
};

template <typename T = std::string>
std::optional<T> get_optional(const nlohmann::json &j, const char *name) {
  const auto it = j.find(name);
  if (it != j.end()) {
    return std::optional<T>(j[name]);
  } else {
    return std::nullopt;
  }
}

void from_json(const nlohmann::json &nlohmann_json_j,
               TestMetaData &nlohmann_json_t) {
  auto &j = nlohmann_json_j["args"];

  nlohmann_json_t.sprt =
      j.contains("sprt") ? std::optional<bool>(true) : std::nullopt;

  nlohmann_json_t.book = get_optional(j, "book");
  nlohmann_json_t.new_tc = get_optional(j, "new_tc");
  nlohmann_json_t.resolved_base = get_optional(j, "resolved_base");
  nlohmann_json_t.resolved_new = get_optional(j, "resolved_new");
  nlohmann_json_t.tc = get_optional(j, "tc");
  nlohmann_json_t.threads = get_optional<int>(j, "threads");

  auto &jr = nlohmann_json_j["results"];
  nlohmann_json_t.pentanomial =
      get_optional<std::vector<int>>(jr, "pentanomial");
}

[[nodiscard]] inline std::vector<std::string>
get_files(const std::string &path, bool recursive = false) {
  std::vector<std::string> files;

  for (const auto &entry : std::filesystem::directory_iterator(path)) {
    if (std::filesystem::is_regular_file(entry)) {
      std::string stem = entry.path().stem().string();
      std::string extension = entry.path().extension().string();
      if (extension == ".gz") {
        if (stem.size() >= 4 && stem.substr(stem.size() - 4) == ".pgn") {
          files.push_back(entry.path().string());
        }
      } else if (extension == ".pgn") {
        files.push_back(entry.path().string());
      }
    } else if (recursive && std::filesystem::is_directory(entry)) {
      auto subdir_files = get_files(entry.path().string(), true);
      files.insert(files.end(), subdir_files.begin(), subdir_files.end());
    }
  }

  return files;
}

[[nodiscard]] inline std::vector<std::vector<std::string>>
split_chunks(const std::vector<std::string> &pgns, int target_chunks) {
  const int chunks_size = (pgns.size() + target_chunks - 1) / target_chunks;

  auto begin = pgns.begin();
  auto end = pgns.end();

  std::vector<std::vector<std::string>> chunks;

  while (begin != end) {
    auto next =
        std::next(begin, std::min(chunks_size,
                                  static_cast<int>(std::distance(begin, end))));
    chunks.push_back(std::vector<std::string>(begin, next));
    begin = next;
  }

  return chunks;
}

class CommandLine {
public:
  CommandLine(int argc, char const *argv[]) {
    for (int i = 1; i < argc; ++i) {
      args.emplace_back(argv[i]);
    }
  }

  bool has_argument(const std::string &arg,
                    bool without_parameter = false) const {
    const auto pos = std::find(args.begin(), args.end(), arg);
    return pos != args.end() &&
           (without_parameter || std::next(pos) != args.end());
  }

  std::string get_argument(const std::string &arg,
                           std::string default_value = "") const {
    auto it = std::find(args.begin(), args.end(), arg);

    if (it != args.end() && std::next(it) != args.end()) {
      return *std::next(it);
    }

    return default_value;
  }

private:
  std::vector<std::string> args;
};
