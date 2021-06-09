// $Id$
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Contribution from Roger Sayle
#include <vector>
#include "Pubchem.h"
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolOps.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace {
struct Patterns {

};

boost::flyweight<std::unique_ptr<Patterns>, boost::flyweights::no_tracking>
    gpats;
void GenerateFP(const RDKit::ROMol &mol, ExplicitBitVect &fp) {
  if (!gpats.get()) {
    gpats = std::unique_ptr<Patterns>(new Patterns());
  }
  const Patterns &pats = *(gpats.get());
  PRECONDITION(fp.size() == 882, "bad fingerprint");
  fp.clearBits();

  if (!mol.getNumAtoms()) {
    return;
  }

  std::vector<RDKit::MatchVectType> matches;
  RDKit::RWMol::ConstAtomIterator atom;
  RDKit::MatchVectType match;
  unsigned int count;

  // for (atom = mol.beginAtoms(); atom != mol.endAtoms(); ++atom) {
  //   switch ((*atom)->getAtomicNum()) {
  //   }
  // }


  /* BIT 125 ================》Special comment this is modified to consider the code useless */
  // RDKit::RingInfo *info = mol.getRingInfo();
  // unsigned int ringcount = info->numRings();
  // unsigned int nArom = 0;
  // for (unsigned int i = 0; i < ringcount; i++) {
  //   bool isArom = true;
  //   const std::vector<int> *ring = &info->bondRings()[i];
  //   std::vector<int>::const_iterator iter;
  //   for (iter = ring->begin(); iter != ring->end(); ++iter) {
  //     if (!mol.getBondWithIdx(*iter)->getIsAromatic()) {
  //       isArom = false;
  //       break;
  //     }
  //   }
  //   if (isArom) {
  //     if (nArom) {
  //       fp.setBit(125);
  //       break;
  //     } else {
  //       nArom++;
  //     }
  //   }
  // }

  /* BIT 881 */
  std::vector<int> mapping;
  if (RDKit::MolOps::getMolFrags(mol, mapping) > 1) {
    fp.setBit(881);
  }
}
}  // namespace

namespace RDKit {
namespace PubchemFingerprints {
ExplicitBitVect *getFingerprintAsBitVect(const ROMol &mol) {
  auto *fp = new ExplicitBitVect(882);
  GenerateFP(mol, *fp);
  return fp;
}
}  // namespace PubchemFingerprints
}  // namespace RDKit
