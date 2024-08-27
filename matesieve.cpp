#include "matesieve.hpp"

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "external/chess.hpp"
#include "external/gzip/gzstream.h"
#include "external/parallel_hashmap/phmap.h"
#include "external/threadpool.hpp"

namespace fs = std::filesystem;
using json = nlohmann::json;

using namespace chess;

using PackedBoard = std::array<std::uint8_t, 24>;

namespace std {
template <> struct hash<PackedBoard> {
  size_t operator()(const PackedBoard pbfen) const {
    std::string_view sv(reinterpret_cast<const char *>(pbfen.data()),
                        pbfen.size());
    return std::hash<std::string_view>{}(sv);
  }
};
} // namespace std

// unordered map to count positions (using compressed format)
using fen_map_t = phmap::parallel_flat_hash_map<
    PackedBoard, std::int64_t, std::hash<PackedBoard>,
    std::equal_to<PackedBoard>,
    std::allocator<std::pair<PackedBoard, std::int64_t>>, 8, std::mutex>;

fen_map_t fen_map;

// map to collect metadata for tests
using map_meta = std::unordered_map<std::string, TestMetaData>;

std::atomic<std::size_t> total_files = 0;
std::atomic<std::size_t> total_games = 0;

namespace analysis {

class Analyze : public pgn::Visitor {
public:
  Analyze(const std::string &regex_engine, std::mutex &progress_output)
      : regex_engine(regex_engine), progress_output(progress_output) {}

  virtual ~Analyze() {}

  void startPgn() override {}

  void header(std::string_view key, std::string_view value) override {

    if (key == "FEN") {
      board.setFen(value);
    }

    if (key == "Variant" && value == "fischerandom") {
      board.set960(true);
    }

    if (key == "White") {
      white = value;
    }

    if (key == "Black") {
      black = value;
    }

    if (key == "WhiteElo") {
      whiteElo = std::atoi(value.data());
    }
    if (key == "BlackElo") {
      blackElo = std::atoi(value.data());
    }
  }

  void startMoves() override {
    do_filter = !regex_engine.empty();

    if (do_filter) {
      if (white.empty() || black.empty()) {
        this->skipPgn(true);
        return;
      }

      std::regex regex(regex_engine);

      if (std::regex_match(white, regex)) {
        filter_side = Color::WHITE;
      }

      if (std::regex_match(black, regex)) {
        if (filter_side == Color::NONE) {
          filter_side = Color::BLACK;
        } else {
          do_filter = false;
        }
      }
    }
    total_games++;
  }

  void move(std::string_view move, std::string_view comment) override {

    // openbench uses Nf3 {+0.57 17/28 583 363004}, fishtest Nf3 {+0.57/17}
    const size_t delimiter_pos = comment.find_first_of(" /");

    if (!do_filter || filter_side == board.sideToMove()) {
      if (delimiter_pos != std::string::npos && comment != "book") {
        const auto match_eval = comment.substr(0, delimiter_pos);

        if (match_eval[1] == 'M') {
          auto key = Board::Compact::encode(board);
          std::int64_t mate =
              std::stoll(std::string(comment.substr(2, delimiter_pos)));
          if (match_eval[0] == '-') {
            mate = -mate;
          }

          fen_map.lazy_emplace_l(
              std::move(key),
              [&](fen_map_t::value_type &p) {
                if (std::abs(mate) < std::abs(p.second))
                  p.second = mate;
              },
              [&](const fen_map_t::constructor &ctor) {
                ctor(std::move(key), mate);
              });
        }
      }
    }
    Move m;

    m = uci::parseSan(board, move, moves);

    // chess-lib may call move() with empty strings for move
    if (m == Move::NO_MOVE) {
      this->skipPgn(true);
      return;
    }

    board.makeMove<true>(m);
  }

  void endPgn() override {
    board.set960(false);
    board.setFen(constants::STARTPOS);

    filter_side = Color::NONE;

    white.clear();
    black.clear();

    whiteElo = blackElo = 0;
  }

private:
  const std::string &regex_engine;
  std::mutex &progress_output;

  Board board;
  Movelist moves;

  bool skip = false;

  bool do_filter = false;
  Color filter_side = Color::NONE;

  std::string white;
  std::string black;

  int whiteElo = 0, blackElo = 0;
};

void ana_files(const std::vector<std::string> &files,
               const std::string &regex_engine, std::mutex &progress_output) {

  for (const auto &file : files) {
    std::string move_counter;
    const auto pgn_iterator = [&](std::istream &iss) {
      auto vis = std::make_unique<Analyze>(regex_engine, progress_output);

      pgn::StreamParser parser(iss);

      try {
        parser.readGames(*vis);
      } catch (const std::exception &e) {
        std::cout << "Error when parsing: " << file << std::endl;
        std::cerr << e.what() << '\n';
      }
    };

    if (file.size() >= 3 && file.substr(file.size() - 3) == ".gz") {
      igzstream input(file.c_str());
      pgn_iterator(input);
    } else {
      std::ifstream pgn_stream(file);
      pgn_iterator(pgn_stream);
      pgn_stream.close();
    }

    ++total_files;

    // Limit the scope of the lock
    {
      const std::lock_guard<std::mutex> lock(progress_output);
      std::cout << "\rProcessed " << total_files << " files" << std::flush;
    }
  }
}

} // namespace analysis

[[nodiscard]] map_meta get_metadata(const std::vector<std::string> &file_list,
                                    bool allow_duplicates) {
  map_meta meta_map;
  std::unordered_map<std::string, std::string>
      test_map; // map to check for duplicate tests
  std::set<std::string> test_warned;
  for (const auto &pathname : file_list) {
    fs::path path(pathname);
    std::string filename = path.filename().string();
    std::string test_id = filename.substr(0, filename.find_first_of("-."));
    std::string test_filename = (path.parent_path() / test_id).string();

    if (test_map.find(test_id) == test_map.end()) {
      test_map[test_id] = test_filename;
    } else if (test_map[test_id] != test_filename) {
      if (test_warned.find(test_filename) == test_warned.end()) {
        std::cout << (allow_duplicates ? "Warning" : "Error")
                  << ": Detected a duplicate of test " << test_id
                  << " in directory " << path.parent_path().string()
                  << std::endl;
        test_warned.insert(test_filename);

        if (!allow_duplicates) {
          std::cout << "Use --allowDuplicates to continue nonetheless."
                    << std::endl;
          std::exit(1);
        }
      }
    }

    // load the JSON data from disk, only once for each test
    if (meta_map.find(test_filename) == meta_map.end()) {
      std::ifstream json_file(test_filename + ".json");

      if (!json_file.is_open())
        continue;

      json metadata = json::parse(json_file);

      meta_map[test_filename] = metadata.get<TestMetaData>();
    }
  }
  return meta_map;
}

void filter_files_book(std::vector<std::string> &file_list,
                       const map_meta &meta_map, const std::regex &regex_book,
                       bool invert) {
  const auto pred = [&regex_book, invert,
                     &meta_map](const std::string &pathname) {
    fs::path path(pathname);
    std::string filename = path.filename().string();
    std::string test_id = filename.substr(0, filename.find_first_of("-."));
    std::string test_filename = (path.parent_path() / test_id).string();

    // check if metadata and "book" entry exist
    if (meta_map.find(test_filename) != meta_map.end() &&
        meta_map.at(test_filename).book.has_value()) {
      bool match =
          std::regex_match(meta_map.at(test_filename).book.value(), regex_book);

      return invert ? match : !match;
    }

    // missing metadata or "book" entry can never match
    return true;
  };

  file_list.erase(std::remove_if(file_list.begin(), file_list.end(), pred),
                  file_list.end());
}

void process(const std::vector<std::string> &files_pgn,
             const std::string &regex_engine, int concurrency) {
  // Create more chunks than threads to prevent threads from idling.
  int target_chunks = 4 * concurrency;

  auto files_chunked = split_chunks(files_pgn, target_chunks);

  std::cout << "Found " << files_pgn.size() << " .pgn(.gz) files, creating "
            << files_chunked.size() << " chunks for processing." << std::endl;

  // Mutex for progress output
  std::mutex progress_output;

  // Create a thread pool
  ThreadPool pool(concurrency);

  for (const auto &files : files_chunked) {

    pool.enqueue([&files, &regex_engine, &progress_output, &files_chunked]() {
      analysis::ana_files(files, regex_engine, progress_output);
    });
  }

  // Wait for all threads to finish
  pool.wait();
}

void print_usage(char const *program_name) {
  std::stringstream ss;

  // clang-format off
    ss << "Usage: " << program_name << " [options]" << "\n";
    ss << "Options:" << "\n";
    ss << "  --file <path>         Path to .pgn(.gz) file" << "\n";
    ss << "  --dir <path>          Path to directory containing .pgn(.gz) files (default: pgns)" << "\n";
    ss << "  -r                    Search for .pgn(.gz) files recursively in subdirectories" << "\n";
    ss << "  --allowDuplicates     Allow duplicate directories for test pgns" << "\n";
    ss << "  --concurrency <N>     Number of concurrent threads to use (default: maximum)" << "\n";
    ss << "  --matchEngine <regex> Filter data based on engine name" << "\n";
    ss << "  --matchBook <regex>   Filter data based on book name" << "\n";
    ss << "  --matchBookInvert     Invert the filter" << "\n";
    ss << "  -o <path>             Path to output epd file (default: matesieve.epd)" << "\n";
    ss << "  --help                Print this help message" << "\n";
  // clang-format on

  std::cout << ss.str();
}

int main(int argc, char const *argv[]) {
  const std::vector<std::string> args(argv + 1, argv + argc);

  std::vector<std::string> files_pgn;
  std::string regex_engine, regex_book, filename = "matesieve.epd";
  int concurrency = std::max(1, int(std::thread::hardware_concurrency()));

  std::vector<std::string>::const_iterator pos;

  if (std::find(args.begin(), args.end(), "--help") != args.end()) {
    print_usage(argv[0]);
    return 0;
  }

  if (find_argument(args, pos, "--concurrency")) {
    concurrency = std::stoi(*std::next(pos));
  }

  if (find_argument(args, pos, "--file")) {
    files_pgn = {*std::next(pos)};
  } else {
    std::string path = "./pgns";

    if (find_argument(args, pos, "--dir")) {
      path = *std::next(pos);
    }

    bool recursive = find_argument(args, pos, "-r", true);
    std::cout << "Looking " << (recursive ? "(recursively) " : "")
              << "for pgn files in " << path << std::endl;

    files_pgn = get_files(path, recursive);
  }

  // sort to easily check for "duplicate" files, i.e. "foo.pgn.gz" and "foo.pgn"
  std::sort(files_pgn.begin(), files_pgn.end());

  for (size_t i = 1; i < files_pgn.size(); ++i) {
    if (files_pgn[i].find(files_pgn[i - 1]) == 0) {
      std::cout << "Error: \"Duplicate\" files: " << files_pgn[i - 1] << " and "
                << files_pgn[i] << std::endl;
      std::exit(1);
    }
  }

  bool allow_duplicates = find_argument(args, pos, "--allowDuplicates", true);
  auto meta_map = get_metadata(files_pgn, allow_duplicates);

  if (find_argument(args, pos, "--matchBook")) {
    regex_book = *std::next(pos);

    if (!regex_book.empty()) {
      bool invert = find_argument(args, pos, "--matchBookInvert", true);
      std::cout << "Filtering pgn files " << (invert ? "not " : "")
                << "matching the book name " << regex_book << std::endl;
      std::regex regex(regex_book);
      filter_files_book(files_pgn, meta_map, regex, invert);
    }
  }

  if (find_argument(args, pos, "--matchEngine")) {
    regex_engine = *std::next(pos);
  }

  if (find_argument(args, pos, "-o")) {
    filename = *std::next(pos);
  }

  std::ofstream out_file(filename);

  const auto t0 = std::chrono::high_resolution_clock::now();

  process(files_pgn, regex_engine, concurrency);

  for (const auto &pair : fen_map) {
    auto board = Board::Compact::decode(pair.first);
    std::string fen = board.getFen(false);
    out_file << fen << " bm #" << pair.second << ";\n";
  }
  out_file.close();

  const auto t1 = std::chrono::high_resolution_clock::now();

  std::cout << "\nSaved " << fen_map.size() << " unique mates from "
            << total_games << " games to " << filename << "."
            << "\nTotal time for processing: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                       .count() /
                   1000.0
            << " s" << std::endl;

  return 0;
}
