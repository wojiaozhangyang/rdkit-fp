
#include <RDGeneral/export.h>
#ifndef __RD_PubchemFPS_H__
#define __RD_PubchemFPS_H__
#include <string>

class ExplicitBitVect;
namespace RDKit {
class ROMol;
namespace PubchemFingerprints {
const std::string pubchemFingerprintVersion = "2.0.0";


RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *getFingerprintAsBitVect(
    const ROMol &mol);
}  // namespace PubchemFingerprints
}  // namespace RDKit

#endif
