# Find mates in (fishtest) pgn data

Given a set of .pgn(.gz) files with engine evaluations, find all the positions 
with a mate score and store them in the format used by [matetools](https://github.com/robertnurnberg/matetools).

Example usage:

```
> ./matesieve --dir ./pgns -r
Looking (recursively) for pgn files in ./pgns
Found 292 .pgn(.gz) files, creating 30 chunks for processing.
Processed 292 files
Saved 1734136 unique mates from 3804577 games to matesieve.epd.
Total time for processing: 68.72 s
```

```
Usage: ./matesieve [options]
Options:
  --file <path>         Path to .pgn(.gz) file
  --dir <path>          Path to directory containing .pgn(.gz) files (default: pgns)
  -r                    Search for .pgn(.gz) files recursively in subdirectories
  --allowDuplicates     Allow duplicate directories for test pgns
  --concurrency <N>     Number of concurrent threads to use (default: maximum)
  --matchEngine <regex> Filter data based on engine name
  --matchBook <regex>   Filter data based on book name
  --matchBookInvert     Invert the filter
  -o <path>             Path to output epd file (default: matesieve.epd)
  --help                Print this help message
```

The code is based on a [WDL model](https://github.com/official-stockfish/WDL_model) and [fastpopular](https://github.com/vondele/fastpopular).
