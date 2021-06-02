// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <fuzzer/FuzzedDataProvider.h>

#include <memory>
#include <sstream>
#include <string>

#include "GraphMol/FileParsers/FileParsers.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  FuzzedDataProvider fdp(data, size);

  const bool sanitize = fdp.ConsumeIntegralInRange(0, 1);
  const bool remove_hs = fdp.ConsumeIntegralInRange(0, 1);
  const bool strict_parsing = fdp.ConsumeIntegralInRange(0, 1);
  std::istringstream data_stream(fdp.ConsumeRemainingBytesAsString());
  unsigned int num_lines = 0;  // output parameter.

  std::unique_ptr<RDKit::RWMol> result = nullptr;

  try {
    result.reset(RDKit::MolDataStreamToMol(&data_stream, num_lines, sanitize,
                                           remove_hs, strict_parsing));
  } catch (...) {
  }

  return 0;
}
