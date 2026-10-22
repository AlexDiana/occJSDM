#!/usr/bin/env python3
"""Make an instrumented package copy; never alter production sampler code.

Usage: python3 dev/simstudy/audit_rng_threads.py SOURCE EMPTY_DESTINATION
Then run validate_rng_safety.R with --source=DESTINATION.
The audit fails immediately if any instrumented RNG leaf leaves R's thread.
"""
from pathlib import Path
import re
import shutil
import sys

source, destination = map(lambda p: Path(p).resolve(), sys.argv[1:])
destination.mkdir(parents=True, exist_ok=False)
for filename in ("DESCRIPTION", "NAMESPACE"):
    shutil.copy2(source / filename, destination / filename)
for directory in ("R", "src"):
    shutil.copytree(source / directory, destination / directory,
                    ignore=shutil.ignore_patterns("*.o", "*.so", "*.dll", "*.dylib"))

header = destination / "src" / "rng.h"
text = header.read_text()
guard = r'''
// Validation-only instrumentation, installed in a separate package copy.
#include <thread>
#include <map>
#include <string>
#include <stdexcept>
inline std::map<std::string, unsigned long>& occjsdm_audit_counts() {
  static std::map<std::string, unsigned long> counts;
  return counts;
}
inline void occjsdm_audit_touch(const char* leaf) {
  // setOccJSDMSeed() establishes this on R's main thread before sampling.
  static const std::thread::id owner = std::this_thread::get_id();
  if (std::this_thread::get_id() != owner)
    throw std::runtime_error(std::string("RNG outside main thread: ") + leaf);
  ++occjsdm_audit_counts()[leaf];
}
'''
text = text.replace("#define OCCJSDM_RNG_H", "#define OCCJSDM_RNG_H\n" + guard)
for name in ("runif", "rnorm"):
    signature = f"inline double {name}() {{"
    assert text.count(signature) == 1
    text = text.replace(signature, signature + f'\n  occjsdm_audit_touch("cpp {name}");')
header.write_text(text)

for filename in ("functions.cpp", "jsdm.cpp"):
    path = destination / "src" / filename
    result = []
    for line in path.read_text().splitlines():
        # Every actual R/Armadillo RNG statement is a standalone statement
        # in these two sources. Comments deliberately earn no instrumentation.
        code = line.split("//", 1)[0]
        matches = re.findall(r"R::r[a-z]+|arma::randn", code)
        if matches:
            indent = line[:len(line) - len(line.lstrip())]
            result.append(indent + f'occjsdm_audit_touch("{matches[0]}");')
        result.append(line)
        if "void setOccJSDMSeed(unsigned int seed) {" in code:
            result.append('  occjsdm_audit_touch("seed");')
    path.write_text("\n".join(result) + "\n")

with (destination / "src" / "functions.cpp").open("a") as out:
    out.write('''
// [[Rcpp::export]]
Rcpp::List rngAuditCounts() {
  Rcpp::List counts;
  for (const auto& entry : occjsdm_audit_counts())
    counts[entry.first] = static_cast<double>(entry.second);
  return counts;
}
''')
print(destination)
