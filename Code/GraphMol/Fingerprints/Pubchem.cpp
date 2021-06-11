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
  std::unique_ptr<RDKit::ROMol> bit_1 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[H]"));
  std::unique_ptr<RDKit::ROMol> bit_2 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[H]"));
  std::unique_ptr<RDKit::ROMol> bit_3 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[H]"));
  std::unique_ptr<RDKit::ROMol> bit_4 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[H]"));
  std::unique_ptr<RDKit::ROMol> bit_5 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li]"));
  std::unique_ptr<RDKit::ROMol> bit_6 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li]"));
  std::unique_ptr<RDKit::ROMol> bit_7 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]"));
  std::unique_ptr<RDKit::ROMol> bit_8 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]"));
  std::unique_ptr<RDKit::ROMol> bit_9 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]"));
  std::unique_ptr<RDKit::ROMol> bit_10 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_11 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_12 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_13 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_14 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_15 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_16 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_17 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_18 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_19 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_20 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_21 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_22 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_23 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_24 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[F]"));
  std::unique_ptr<RDKit::ROMol> bit_25 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[F]"));
  std::unique_ptr<RDKit::ROMol> bit_26 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[F]"));
  std::unique_ptr<RDKit::ROMol> bit_27 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Na]"));
  std::unique_ptr<RDKit::ROMol> bit_28 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Na]"));
  std::unique_ptr<RDKit::ROMol> bit_29 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#14]"));
  std::unique_ptr<RDKit::ROMol> bit_30 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#14]"));
  std::unique_ptr<RDKit::ROMol> bit_31 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_32 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_33 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_34 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]"));
  std::unique_ptr<RDKit::ROMol> bit_35 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]"));
  std::unique_ptr<RDKit::ROMol> bit_36 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]"));
  std::unique_ptr<RDKit::ROMol> bit_37 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]"));
  std::unique_ptr<RDKit::ROMol> bit_38 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Cl]"));
  std::unique_ptr<RDKit::ROMol> bit_39 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Cl]"));
  std::unique_ptr<RDKit::ROMol> bit_40 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Cl]"));
  std::unique_ptr<RDKit::ROMol> bit_41 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Cl]"));
  std::unique_ptr<RDKit::ROMol> bit_42 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[K]"));
  std::unique_ptr<RDKit::ROMol> bit_43 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[K]"));
  std::unique_ptr<RDKit::ROMol> bit_44 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Br]"));
  std::unique_ptr<RDKit::ROMol> bit_45 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Br]"));
  std::unique_ptr<RDKit::ROMol> bit_46 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Br]"));
  std::unique_ptr<RDKit::ROMol> bit_47 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[I]"));
  std::unique_ptr<RDKit::ROMol> bit_48 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[I]"));
  std::unique_ptr<RDKit::ROMol> bit_49 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[I]"));
  std::unique_ptr<RDKit::ROMol> bit_50 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Be]"));
  std::unique_ptr<RDKit::ROMol> bit_51 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Mg]"));
  std::unique_ptr<RDKit::ROMol> bit_52 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Al]"));
  std::unique_ptr<RDKit::ROMol> bit_53 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ca]"));
  std::unique_ptr<RDKit::ROMol> bit_54 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Sc]"));
  std::unique_ptr<RDKit::ROMol> bit_55 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ti]"));
  std::unique_ptr<RDKit::ROMol> bit_56 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[V]"));
  std::unique_ptr<RDKit::ROMol> bit_57 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Cr]"));
  std::unique_ptr<RDKit::ROMol> bit_58 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Mn]"));
  std::unique_ptr<RDKit::ROMol> bit_59 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Fe]"));
  std::unique_ptr<RDKit::ROMol> bit_60 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Co]"));
  std::unique_ptr<RDKit::ROMol> bit_61 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ni]"));
  std::unique_ptr<RDKit::ROMol> bit_62 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Cu]"));
  std::unique_ptr<RDKit::ROMol> bit_63 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Zn]"));
  std::unique_ptr<RDKit::ROMol> bit_64 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ga]"));
  std::unique_ptr<RDKit::ROMol> bit_65 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#32]"));
  std::unique_ptr<RDKit::ROMol> bit_66 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#33]"));
  std::unique_ptr<RDKit::ROMol> bit_67 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#34]"));
  std::unique_ptr<RDKit::ROMol> bit_68 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Kr]"));
  std::unique_ptr<RDKit::ROMol> bit_69 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Rb]"));
  std::unique_ptr<RDKit::ROMol> bit_70 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Sr]"));
  std::unique_ptr<RDKit::ROMol> bit_71 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Y]"));
  std::unique_ptr<RDKit::ROMol> bit_72 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Zr]"));
  std::unique_ptr<RDKit::ROMol> bit_73 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Nb]"));
  std::unique_ptr<RDKit::ROMol> bit_74 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Mo]"));
  std::unique_ptr<RDKit::ROMol> bit_75 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ru]"));
  std::unique_ptr<RDKit::ROMol> bit_76 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Rh]"));
  std::unique_ptr<RDKit::ROMol> bit_77 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Pd]"));
  std::unique_ptr<RDKit::ROMol> bit_78 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ag]"));
  std::unique_ptr<RDKit::ROMol> bit_79 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Cd]"));
  std::unique_ptr<RDKit::ROMol> bit_80 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[In]"));
  std::unique_ptr<RDKit::ROMol> bit_81 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Sn]"));
  std::unique_ptr<RDKit::ROMol> bit_82 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Sb]"));
  std::unique_ptr<RDKit::ROMol> bit_83 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#52]"));
  std::unique_ptr<RDKit::ROMol> bit_84 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Xe]"));
  std::unique_ptr<RDKit::ROMol> bit_85 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Cs]"));
  std::unique_ptr<RDKit::ROMol> bit_86 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ba]"));
  std::unique_ptr<RDKit::ROMol> bit_87 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Lu]"));
  std::unique_ptr<RDKit::ROMol> bit_88 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Hf]"));
  std::unique_ptr<RDKit::ROMol> bit_89 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ta]"));
  std::unique_ptr<RDKit::ROMol> bit_90 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[W]"));
  std::unique_ptr<RDKit::ROMol> bit_91 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Re]"));
  std::unique_ptr<RDKit::ROMol> bit_92 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Os]"));
  std::unique_ptr<RDKit::ROMol> bit_93 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ir]"));
  std::unique_ptr<RDKit::ROMol> bit_94 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Pt]"));
  std::unique_ptr<RDKit::ROMol> bit_95 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Au]"));
  std::unique_ptr<RDKit::ROMol> bit_96 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Hg]"));
  std::unique_ptr<RDKit::ROMol> bit_97 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Tl]"));
  std::unique_ptr<RDKit::ROMol> bit_98 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Pb]"));
  std::unique_ptr<RDKit::ROMol> bit_99 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Bi]"));
  std::unique_ptr<RDKit::ROMol> bit_100 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[La]"));
  std::unique_ptr<RDKit::ROMol> bit_101 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ce]"));
  std::unique_ptr<RDKit::ROMol> bit_102 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Pr]"));
  std::unique_ptr<RDKit::ROMol> bit_103 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Nd]"));
  std::unique_ptr<RDKit::ROMol> bit_104 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Pm]"));
  std::unique_ptr<RDKit::ROMol> bit_105 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Sm]"));
  std::unique_ptr<RDKit::ROMol> bit_106 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Eu]"));
  std::unique_ptr<RDKit::ROMol> bit_107 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Gd]"));
  std::unique_ptr<RDKit::ROMol> bit_108 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Tb]"));
  std::unique_ptr<RDKit::ROMol> bit_109 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Dy]"));
  std::unique_ptr<RDKit::ROMol> bit_110 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Ho]"));
  std::unique_ptr<RDKit::ROMol> bit_111 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Er]"));
  std::unique_ptr<RDKit::ROMol> bit_112 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Tm]"));
  std::unique_ptr<RDKit::ROMol> bit_113 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Yb]"));
  std::unique_ptr<RDKit::ROMol> bit_114 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Tc]"));
  std::unique_ptr<RDKit::ROMol> bit_115 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[U]"));
  std::unique_ptr<RDKit::ROMol> bit_116 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_117 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1cc1"));
  std::unique_ptr<RDKit::ROMol> bit_118 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aa1"));
  std::unique_ptr<RDKit::ROMol> bit_119 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aa1"));
  std::unique_ptr<RDKit::ROMol> bit_120 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_121 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_122 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_123 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_124 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1cc1"));
  std::unique_ptr<RDKit::ROMol> bit_125 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aa1"));
  std::unique_ptr<RDKit::ROMol> bit_126 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aa1"));
  std::unique_ptr<RDKit::ROMol> bit_127 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_128 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_129 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_130 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_131 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_132 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaa1"));
  std::unique_ptr<RDKit::ROMol> bit_133 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaa1"));
  std::unique_ptr<RDKit::ROMol> bit_134 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_135 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_136 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_137 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_138 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_139 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaa1"));
  std::unique_ptr<RDKit::ROMol> bit_140 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaa1"));
  std::unique_ptr<RDKit::ROMol> bit_141 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_142 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_143 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_144 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_145 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_146 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_147 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_148 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_149 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_150 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_151 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_152 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_153 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_154 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_155 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_156 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_157 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_158 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_159 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_160 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_161 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_162 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_163 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_164 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_165 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_166 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_167 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_168 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_169 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_170 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_171 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_172 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_173 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_174 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_175 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_176 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_177 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_178 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_179 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_180 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1ccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_181 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_182 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_183 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_184 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_185 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_186 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_187 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1ccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_188 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_189 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_190 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_191 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_192 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_193 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_194 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1ccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_195 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_196 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_197 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_198 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_199 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_200 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_201 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1ccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_202 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_203 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_204 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_205 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_206 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_207 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_208 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1ccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_209 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_210 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_211 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_212 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_213 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_214 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_215 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1cccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_216 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_217 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_218 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_219 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_220 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_221 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_222 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1cccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_223 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_224 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_225 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_226 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_227 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_228 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_229 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1ccccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_230 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_231 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_232 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_233 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_234 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_235 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_236 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1ccccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_237 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_238 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_239 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_240 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_241 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_242 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~1~*~*~*~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_243 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1cccccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_244 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_245 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_246 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("C~1~C~C~C~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_247 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N~1~A~A~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_248 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_249 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("*~1~*~*~*~*~*~*~*~*~*1"));
  std::unique_ptr<RDKit::ROMol> bit_250 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c1ccccccccc1"));
  std::unique_ptr<RDKit::ROMol> bit_251 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n1aaaaaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_252 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[a;!#6]1aaaaaaaaa1"));
  std::unique_ptr<RDKit::ROMol> bit_253 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("C~1~C~C~C~C~C~C~C~C~C1"));
  std::unique_ptr<RDKit::ROMol> bit_254 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("N~1~A~A~A~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_255 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[A;!#6]~1~A~A~A~A~A~A~A~A~A1"));
  std::unique_ptr<RDKit::ROMol> bit_256 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[aR]"));
  std::unique_ptr<RDKit::ROMol> bit_257 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[aR;!#6]"));
  std::unique_ptr<RDKit::ROMol> bit_258 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[aR]"));
  std::unique_ptr<RDKit::ROMol> bit_259 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[aR;!#6]"));
  std::unique_ptr<RDKit::ROMol> bit_260 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[aR]"));
  std::unique_ptr<RDKit::ROMol> bit_261 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[aR;!#6]"));
  std::unique_ptr<RDKit::ROMol> bit_262 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[aR]"));
  std::unique_ptr<RDKit::ROMol> bit_263 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[aR;!#6]"));
  std::unique_ptr<RDKit::ROMol> bit_264 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_265 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li]~[Li]"));
  std::unique_ptr<RDKit::ROMol> bit_266 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li]~[#5]"));
  std::unique_ptr<RDKit::ROMol> bit_267 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li]~[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_268 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li]~[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_269 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li]~[F]"));
  std::unique_ptr<RDKit::ROMol> bit_270 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li]~[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_271 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li]~[#16]"));
  std::unique_ptr<RDKit::ROMol> bit_272 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Li]~[Cl]"));
  std::unique_ptr<RDKit::ROMol> bit_273 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_274 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]~[#5]"));
  std::unique_ptr<RDKit::ROMol> bit_275 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]~[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_276 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]~[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_277 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]~[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_278 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]~[F]"));
  std::unique_ptr<RDKit::ROMol> bit_279 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]~[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_280 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]~[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_281 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]~[#16]"));
  std::unique_ptr<RDKit::ROMol> bit_282 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]~[Cl]"));
  std::unique_ptr<RDKit::ROMol> bit_283 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#5]~[Br]"));
  std::unique_ptr<RDKit::ROMol> bit_284 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_285 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_286 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_287 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_288 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[F]"));
  std::unique_ptr<RDKit::ROMol> bit_289 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[Na]"));
  std::unique_ptr<RDKit::ROMol> bit_290 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[Mg]"));
  std::unique_ptr<RDKit::ROMol> bit_291 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[Al]"));
  std::unique_ptr<RDKit::ROMol> bit_292 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[#14]"));
  std::unique_ptr<RDKit::ROMol> bit_293 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_294 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[#16]"));
  std::unique_ptr<RDKit::ROMol> bit_295 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[Cl]"));
  std::unique_ptr<RDKit::ROMol> bit_296 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[#33]"));
  std::unique_ptr<RDKit::ROMol> bit_297 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[#34]"));
  std::unique_ptr<RDKit::ROMol> bit_298 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[Br]"));
  std::unique_ptr<RDKit::ROMol> bit_299 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[I]"));
  std::unique_ptr<RDKit::ROMol> bit_300 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_301 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_302 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_303 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[F]"));
  std::unique_ptr<RDKit::ROMol> bit_304 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[#14]"));
  std::unique_ptr<RDKit::ROMol> bit_305 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_306 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[#16]"));
  std::unique_ptr<RDKit::ROMol> bit_307 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[Cl]"));
  std::unique_ptr<RDKit::ROMol> bit_308 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[Br]"));
  std::unique_ptr<RDKit::ROMol> bit_309 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_310 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_311 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~[Mg]"));
  std::unique_ptr<RDKit::ROMol> bit_312 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~[Na]"));
  std::unique_ptr<RDKit::ROMol> bit_313 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~[Al]"));
  std::unique_ptr<RDKit::ROMol> bit_314 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~[#14]"));
  std::unique_ptr<RDKit::ROMol> bit_315 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_316 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~[K]"));
  std::unique_ptr<RDKit::ROMol> bit_317 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[F]~[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_318 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[F]~[#16]"));
  std::unique_ptr<RDKit::ROMol> bit_319 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Al;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_320 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[Al]~[Cl]"));
  std::unique_ptr<RDKit::ROMol> bit_321 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#14;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_322 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#14]~[#14]"));
  std::unique_ptr<RDKit::ROMol> bit_323 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#14]~[Cl]"));
  std::unique_ptr<RDKit::ROMol> bit_324 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#15;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_325 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#15]~[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_326 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#33;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_327 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#33]~[#33]"));
  std::unique_ptr<RDKit::ROMol> bit_328 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~Br)(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_329 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](~Br)(~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_330 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](~Br)"));
  std::unique_ptr<RDKit::ROMol> bit_331 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~Br)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_332 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~Br)(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_333 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_334 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](~[#6])(~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_335 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](~[#6])(~[#6])(~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_336 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0](~[#6])(~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_337 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](~[#6])(~[#6])(~[#6])(~[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_338 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](~[#6])(~[#6])(~[#6])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_339 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0](~[#6])(~[#6])(~[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_340 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0](~[#6])(~[#6])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_341 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](~[#6])(~[#6])(~[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_342 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](~[#6])(~[#6])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_343 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~[#6])(~Cl)"));
  std::unique_ptr<RDKit::ROMol> bit_344 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](~[#6])(~Cl)"));
  std::unique_ptr<RDKit::ROMol> bit_345 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_346 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0](~[#6])(~[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_347 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0](~[#6])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_348 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0](~[#6])(~[#8])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_349 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0](~[#6])(~[#15])"));
  std::unique_ptr<RDKit::ROMol> bit_350 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0](~[#6])(~[#16])"));
  std::unique_ptr<RDKit::ROMol> bit_351 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~[#6])(~I)"));
  std::unique_ptr<RDKit::ROMol> bit_352 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~[#6])(~[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_353 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~[#6])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_354 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~[#6])(~[#16])"));
  std::unique_ptr<RDKit::ROMol> bit_355 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~[#6])(~[#14])"));
  std::unique_ptr<RDKit::ROMol> bit_356 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#6])(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_357 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#6])(:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_358 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#6])(:c)(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_359 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#6])(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_360 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#6])(:n)(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_361 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~Cl)(~Cl)"));
  std::unique_ptr<RDKit::ROMol> bit_362 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](~Cl)"));
  std::unique_ptr<RDKit::ROMol> bit_363 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~Cl)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_364 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~F)(~F)"));
  std::unique_ptr<RDKit::ROMol> bit_365 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~F)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_366 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](~[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_367 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_368 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0](~[#8])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_369 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](~[#16])"));
  std::unique_ptr<RDKit::ROMol> bit_370 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](~[#14])"));
  std::unique_ptr<RDKit::ROMol> bit_371 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[c;!H0](:c)"));
  std::unique_ptr<RDKit::ROMol> bit_372 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[c;!H0](:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_373 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[c;!H0](:c)(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_374 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[c;!H0](:n)"));
  std::unique_ptr<RDKit::ROMol> bit_375 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0;!H1;!H2]"));
  std::unique_ptr<RDKit::ROMol> bit_376 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~[#7])(~[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_377 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#7])(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_378 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#7])(:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_379 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#7])(:c)(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_380 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#7])(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_381 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](~[#8])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_382 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#8])(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_383 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#8])(:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_384 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(~[#16])(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_385 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_386 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(:c)(:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_387 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(:c)(:c)(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_388 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(:c)(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_389 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(:c)(:n)(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_390 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c(:n)(:n)"));
  std::unique_ptr<RDKit::ROMol> bit_391 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7](~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_392 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7](~[#6])(~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_393 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7;!H0](~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_394 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7;!H0](~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_395 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7;!H0](~[#6])(~[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_396 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7](~[#6])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_397 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n(~[#6])(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_398 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n(~[#6])(:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_399 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7;!H0](~[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_400 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[n;!H0](:c)"));
  std::unique_ptr<RDKit::ROMol> bit_401 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[n;!H0](:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_402 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7](~[#8])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_403 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n(~[#8])(:o)"));
  std::unique_ptr<RDKit::ROMol> bit_404 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n(:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_405 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n(:c)(:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_406 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8](~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_407 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8;!H0](~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_408 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8](~[#6])(~[#15])"));
  std::unique_ptr<RDKit::ROMol> bit_409 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8;!H0](~[#16])"));
  std::unique_ptr<RDKit::ROMol> bit_410 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("o(:c)(:c)"));
  std::unique_ptr<RDKit::ROMol> bit_411 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#15](~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_412 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#15](~[#8])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_413 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16](~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_414 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16;!H0](~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_415 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16](~[#6])(~[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_416 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#14](~[#6])(~[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_417 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_418 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]#[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_419 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_420 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]#[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_421 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_422 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=,:[s,S]"));
  std::unique_ptr<RDKit::ROMol> bit_423 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]=,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_424 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_425 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]=,:[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_426 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#15]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_427 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#15]=,:[#15]"));
  std::unique_ptr<RDKit::ROMol> bit_428 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]#[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_429 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0]#[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_430 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]#[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_431 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(-,:[#6])(=,:[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_432 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(-,:[#6])(=,:[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_433 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(-,:[#6])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_434 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(Cl)(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_435 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6&!H0](-,:[#6])(=,:[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_436 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6&!H0](-,:[#6])(=,:[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_437 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6&!H0](-,:[#6])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_438 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(-,:[#7])(=,:[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_439 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(-,:[#7])(=,:[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_440 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(-,:[#7])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_441 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(-,:[#8])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_442 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(=,:[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_443 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(=,:[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_444 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#6])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_445 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6](Cl)(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_446 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0](-,:[#7])(=,:[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_447 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](=,:[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_448 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](=,:[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_449 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0](=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_450 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#7])(=,:[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_451 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#7])(=,:[#7])"));
  std::unique_ptr<RDKit::ROMol> bit_452 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#7])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_453 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6](-,:[#8])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_454 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7](-,:[#6])(=,:[#6])"));
  std::unique_ptr<RDKit::ROMol> bit_455 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7](-,:[#6])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_456 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7](-,:[#8])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_457 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#15](-,:[#8])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_458 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16](-,:[#6])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_459 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16](-,:[#8])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_460 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16](=,:[#8])(=,:[#8])"));
  std::unique_ptr<RDKit::ROMol> bit_461 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]-,:[#6]-,:[#6]#C"));
  std::unique_ptr<RDKit::ROMol> bit_462 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]=,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_463 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_464 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n:c-,:[#16;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_465 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_466 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#16]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_467 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N#[#6]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_468 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]=,:[#7]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_469 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#16]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_470 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]-,:[#16]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_471 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:c-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_472 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("s:c:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_473 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:n:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_474 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]-,:c:n:c"));
  std::unique_ptr<RDKit::ROMol> bit_475 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("s:c:c:n"));
  std::unique_ptr<RDKit::ROMol> bit_476 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16]-,:[#6]=,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_477 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#8]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_478 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:[#7]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_479 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16]-,:[#6]=,:[#7;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_480 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16]-,:[#6]-,:[#16]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_481 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:s:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_482 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]-,:[#16]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_483 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:n-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_484 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:[#16]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_485 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:c:n:c"));
  std::unique_ptr<RDKit::ROMol> bit_486 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n:c:c:n"));
  std::unique_ptr<RDKit::ROMol> bit_487 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:c:n:n"));
  std::unique_ptr<RDKit::ROMol> bit_488 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]=,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_489 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]=,:[#7;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_490 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#16]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_491 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_492 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]-,:n:[c;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_493 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:c:o:c"));
  std::unique_ptr<RDKit::ROMol> bit_494 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]=,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_495 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]=,:[#6]-,:c:n"));
  std::unique_ptr<RDKit::ROMol> bit_496 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]-,:[#7]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_497 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n:n-,:[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_498 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]-,:c:c:n"));
  std::unique_ptr<RDKit::ROMol> bit_499 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]=,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_500 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:c:c:n"));
  std::unique_ptr<RDKit::ROMol> bit_501 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]-,:[#16]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_502 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Cl-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_503 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]=,:[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_504 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Cl-,:c:[c;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_505 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n:c:n-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_506 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Cl-,:c:c-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_507 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]-,:c:n:c"));
  std::unique_ptr<RDKit::ROMol> bit_508 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#16]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_509 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16]=,:[#6]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_510 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Br-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_511 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7;!H0]-,:[#7;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_512 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16]=,:[#6]-,:[#7;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_513 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#33]-,:[#8;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_514 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("s:c:[c;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_515 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#7]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_516 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#7]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_517 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0]=,:[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_518 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#7]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_519 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#7]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_520 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]=,:[#6]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_521 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_522 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:n-,:[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_523 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#7]-,:[#7;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_524 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_525 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]=,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_526 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#33]-,:c:[c;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_527 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Cl-,:c:c-,:Cl"));
  std::unique_ptr<RDKit::ROMol> bit_528 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:c:[n;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_529 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7;!H0]-,:[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_530 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl-,:[#6]-,:[#6]-,:Cl"));
  std::unique_ptr<RDKit::ROMol> bit_531 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n:c-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_532 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_533 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]-,:c:[c;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_534 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]-,:c:c-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_535 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]-,:c:c-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_536 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_537 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_538 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_539 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]=,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_540 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]=,:[#6]-,:[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_541 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#7]-,:[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_542 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_543 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]-,:c:[c;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_544 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]-,:c:c-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_545 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]-,:c:c-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_546 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_547 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:c:[c;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_548 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:c:c-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_549 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]-,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_550 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_551 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_552 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl-,:[#6]-,:[#6]-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_553 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:c-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_554 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_555 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Br-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_556 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]=,:[#6]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_557 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]=,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_558 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n:c-,:[#8;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_559 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]=,:[#7]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_560 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#7;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_561 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_562 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl-,:[#6]-,:[#6]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_563 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Br-,:[#6]-,:[#6]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_564 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_565 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]=,:[#6]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_566 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:c-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_567 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_568 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_569 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("N#[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_570 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_571 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:c-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_572 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6;!H0]-,:[#8;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_573 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("n:c:n:c"));
  std::unique_ptr<RDKit::ROMol> bit_574 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_575 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_576 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:c:c-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_577 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]=,:[#6]-,:c:[c;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_578 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:c-,:[#7]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_579 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]-,:c:c-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_580 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_581 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_582 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_583 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_584 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl-,:c:c-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_585 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("c:c-,:[#6]=,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_586 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:c:c-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_587 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#16]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_588 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:c:c-,:[#8;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_589 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_590 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:c:c-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_591 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:c:c-,:[#8;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_592 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_593 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_594 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#6]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_595 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#8]-,:[#6]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_596 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("c:c-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_597 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]=,:[#6]-,:[#7]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_598 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_599 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Cl-,:c:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_600 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0]-,:[#6]=,:[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_601 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:c:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_602 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:c:c:c-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_603 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_604 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]-,:c:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_605 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#8]-,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_606 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_607 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:c:c-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_608 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_609 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_610 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl-,:[#6]-,:[#6]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_611 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#8]-,:[#6]-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_612 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#6]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_613 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#8]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_614 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#7]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_615 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#8]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_616 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#6]-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_617 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:c:n:n:c"));
  std::unique_ptr<RDKit::ROMol> bit_618 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#6]-,:[#8;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_619 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:c-,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_620 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]=,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_621 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("c:c-,:[#8]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_622 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-,:c:c:c:n"));
  std::unique_ptr<RDKit::ROMol> bit_623 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#8]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_624 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_625 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:c:c-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_626 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:c:c-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_627 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#8]-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_628 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]=,:[#33]-,:c:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_629 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#7]-,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_630 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]-,:c:c:c-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_631 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:c:c-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_632 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:c:c-,:[#8;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_633 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#8]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_634 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_635 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_636 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#7]-,:[#6]-,:[#7;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_637 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#7]-,:[#6]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_638 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_639 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_640 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_641 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]=,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_642 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_643 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_644 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6;!H0]-,:[#6]-,:[#7;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_645 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]=,:[#7]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_646 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#7]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_647 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#7]-,:[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_648 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#7]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_649 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#7]-,:c:c-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_650 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#7]-,:c:c-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_651 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#7]-,:[#6]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_652 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]-,:c:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_653 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]-,:c:c:c-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_654 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]-,:c:c:c-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_655 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#7]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_656 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_657 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#7]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_658 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#7]-,:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_659 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#16]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_660 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_661 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]=,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_662 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#8]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_663 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_664 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#8;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_665 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]=,:[#6]-,:[#6]=,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_666 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:c:c-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_667 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]=,:[#6]-,:[#6]-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_668 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]=,:[#6]-,:[#6]-,:[#8;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_669 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:c:c-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_670 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl-,:c:c-,:[#6]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_671 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Br-,:c:c:c-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_672 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]=,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_673 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]=,:[#6;!H0]"));
  std::unique_ptr<RDKit::ROMol> bit_674 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]=,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_675 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#7]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_676 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Br-,:[#6]-,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_677 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("N#[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_678 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]=,:[#6]-,:c:c"));
  std::unique_ptr<RDKit::ROMol> bit_679 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#6]=,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_680 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_681 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_682 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_683 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_684 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_685 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_686 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_687 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_688 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_689 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_690 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_691 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_692 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_693 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_694 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_695 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]=,:[#8]"));
  std::unique_ptr<RDKit::ROMol> bit_696 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#7]"));
  std::unique_ptr<RDKit::ROMol> bit_697 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_698 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6](-,:[#6])-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_699 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_700 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6](-,:[#6])-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_701 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#8]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_702 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6](-,:[#8])-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_703 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#7]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_704 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#8]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6](-,:[#7])-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_705 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_706 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6](-,:[#8])-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_707 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6](=,:[#8])-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_708 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
          "[#8]=,:[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6](-,:[#7])-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_709 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6](-,:[#6])-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_710 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6](-,:[#6])-,:[#6]-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_711 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6]-,:[#6](-,:[#6])-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_712 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6](-,:[#6])(-,:[#6])-,:[#6]-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_713 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6]-,:[#6](-,:[#6])-,:[#6](-,:[#6])-,:[#6]"));
  std::unique_ptr<RDKit::ROMol> bit_714 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1ccc([#6])cc1"));
  std::unique_ptr<RDKit::ROMol> bit_715 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1ccc([#8])cc1"));
  std::unique_ptr<RDKit::ROMol> bit_716 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1ccc([#16])cc1"));
  std::unique_ptr<RDKit::ROMol> bit_717 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1ccc([#7])cc1"));
  std::unique_ptr<RDKit::ROMol> bit_718 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1ccc(Cl)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_719 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1ccc(Br)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_720 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1ccc([#8])cc1"));
  std::unique_ptr<RDKit::ROMol> bit_721 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1ccc([#16])cc1"));
  std::unique_ptr<RDKit::ROMol> bit_722 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1ccc([#7])cc1"));
  std::unique_ptr<RDKit::ROMol> bit_723 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1ccc(Cl)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_724 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1ccc(Br)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_725 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1ccc([#16])cc1"));
  std::unique_ptr<RDKit::ROMol> bit_726 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1ccc([#7])cc1"));
  std::unique_ptr<RDKit::ROMol> bit_727 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1ccc(Cl)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_728 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1ccc(Br)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_729 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]c1ccc([#7])cc1"));
  std::unique_ptr<RDKit::ROMol> bit_730 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]c1ccc(Cl)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_731 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]c1ccc(Br)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_732 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Clc1ccc(Cl)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_733 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Clc1ccc(Br)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_734 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Brc1ccc(Br)cc1"));
  std::unique_ptr<RDKit::ROMol> bit_735 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1cc([#6])ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_736 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1cc([#8])ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_737 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1cc([#16])ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_738 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1cc([#7])ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_739 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1cc(Cl)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_740 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1cc(Br)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_741 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1cc([#8])ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_742 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1cc([#16])ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_743 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1cc([#7])ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_744 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1cc(Cl)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_745 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1cc(Br)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_746 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1cc([#16])ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_747 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1cc([#7])ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_748 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1cc(Cl)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_749 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1cc(Br)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_750 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]c1cc([#7])ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_751 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]c1cc(Cl)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_752 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]c1cc(Br)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_753 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Clc1cc(Cl)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_754 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Clc1cc(Br)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_755 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Brc1cc(Br)ccc1"));
  std::unique_ptr<RDKit::ROMol> bit_756 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1c([#6])cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_757 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1c([#8])cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_758 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1c([#16])cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_759 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1c([#7])cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_760 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1c(Cl)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_761 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]c1c(Br)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_762 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1c([#8])cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_763 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1c([#16])cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_764 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1c([#7])cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_765 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1c(Cl)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_766 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]c1c(Br)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_767 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1c([#16])cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_768 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1c([#7])cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_769 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1c(Cl)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_770 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]c1c(Br)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_771 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]c1c([#7])cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_772 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]c1c(Cl)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_773 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]c1c(Br)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_774 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Clc1c(Cl)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_775 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Clc1c(Br)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_776 =
      std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("Brc1c(Br)cccc1"));
  std::unique_ptr<RDKit::ROMol> bit_777 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6][#6]([#6])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_778 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6][#6]([#8])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_779 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6][#6]([#16])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_780 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6][#6]([#7])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_781 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6][#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_782 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_783 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6][#6]([#8])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_784 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6][#6]([#16])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_785 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6][#6]([#7])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_786 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6][#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_787 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_788 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6][#6]([#16])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_789 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6][#6]([#7])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_790 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6][#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_791 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_792 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6][#6][#6]([#7])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_793 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6][#6][#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_794 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6][#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_795 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl[#6]1[#6][#6][#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_796 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl[#6]1[#6][#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_797 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Br[#6]1[#6][#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_798 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6]([#6])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_799 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6]([#8])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_800 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6]([#16])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_801 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6]([#7])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_802 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6](Cl)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_803 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_804 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6]([#8])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_805 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6]([#16])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_806 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6]([#7])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_807 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6](Cl)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_808 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_809 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6]([#16])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_810 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6]([#7])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_811 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6](Cl)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_812 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_813 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6][#6]([#7])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_814 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6][#6](Cl)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_815 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6][#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_816 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl[#6]1[#6][#6](Cl)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_817 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl[#6]1[#6][#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_818 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Br[#6]1[#6][#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_819 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6]([#6])[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_820 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6]([#8])[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_821 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6]([#16])[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_822 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6]([#7])[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_823 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6](Cl)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_824 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6](Br)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_825 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6]([#8])[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_826 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6]([#16])[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_827 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6]([#7])[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_828 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6](Cl)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_829 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6](Br)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_830 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6]([#16])[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_831 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6]([#7])[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_832 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6](Cl)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_833 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6](Br)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_834 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6]([#7])[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_835 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6](Cl)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_836 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6](Br)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_837 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl[#6]1[#6](Cl)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_838 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl[#6]1[#6](Br)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_839 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Br[#6]1[#6](Br)[#6][#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_840 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6]([#6])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_841 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6]([#8])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_842 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6]([#16])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_843 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6]([#7])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_844 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_845 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_846 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6]([#8])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_847 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6]([#16])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_848 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6]([#7])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_849 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_850 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_851 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6]([#16])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_852 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6]([#7])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_853 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_854 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_855 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6][#6]([#7])[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_856 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6][#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_857 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_858 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl[#6]1[#6][#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_859 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl[#6]1[#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_860 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Br[#6]1[#6][#6](Br)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_861 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6]([#6])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_862 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6]([#8])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_863 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6]([#16])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_864 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6]([#7])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_865 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6](Cl)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_866 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#6][#6]1[#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_867 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6]([#8])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_868 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6]([#16])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_869 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6]([#7])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_870 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6](Cl)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_871 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#8][#6]1[#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_872 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6]([#16])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_873 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6]([#7])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_874 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6](Cl)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_875 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#16][#6]1[#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_876 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6]([#7])[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_877 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6](Cl)[#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_878 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("[#7][#6]1[#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_879 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl[#6]1[#6](Cl)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_880 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Cl[#6]1[#6](Br)[#6][#6][#6]1"));
  std::unique_ptr<RDKit::ROMol> bit_881 = std::unique_ptr<RDKit::ROMol>(
      RDKit::SmartsToMol("Br[#6]1[#6](Br)[#6][#6][#6]1"));
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
  count = RDKit::SubstructMatch(mol, *pats.bit_1, matches, true, true);
  if (count > 4) {
    fp.setBit(1);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_2, matches, true, true);
  if (count > 8) {
    fp.setBit(2);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_3, matches, true, true);
  if (count > 16) {
    fp.setBit(3);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_4, matches, true, true);
  if (count > 32) {
    fp.setBit(4);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_5, matches, true, true);
  if (count > 1) {
    fp.setBit(5);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_6, matches, true, true);
  if (count > 2) {
    fp.setBit(6);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_7, matches, true, true);
  if (count > 1) {
    fp.setBit(7);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_8, matches, true, true);
  if (count > 2) {
    fp.setBit(8);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_9, matches, true, true);
  if (count > 4) {
    fp.setBit(9);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_10, matches, true, true);
  if (count > 2) {
    fp.setBit(10);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_11, matches, true, true);
  if (count > 4) {
    fp.setBit(11);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_12, matches, true, true);
  if (count > 8) {
    fp.setBit(12);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_13, matches, true, true);
  if (count > 16) {
    fp.setBit(13);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_14, matches, true, true);
  if (count > 32) {
    fp.setBit(14);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_15, matches, true, true);
  if (count > 1) {
    fp.setBit(15);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_16, matches, true, true);
  if (count > 2) {
    fp.setBit(16);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_17, matches, true, true);
  if (count > 4) {
    fp.setBit(17);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_18, matches, true, true);
  if (count > 8) {
    fp.setBit(18);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_19, matches, true, true);
  if (count > 1) {
    fp.setBit(19);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_20, matches, true, true);
  if (count > 2) {
    fp.setBit(20);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_21, matches, true, true);
  if (count > 4) {
    fp.setBit(21);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_22, matches, true, true);
  if (count > 8) {
    fp.setBit(22);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_23, matches, true, true);
  if (count > 16) {
    fp.setBit(23);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_24, matches, true, true);
  if (count > 1) {
    fp.setBit(24);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_25, matches, true, true);
  if (count > 2) {
    fp.setBit(25);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_26, matches, true, true);
  if (count > 4) {
    fp.setBit(26);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_27, matches, true, true);
  if (count > 1) {
    fp.setBit(27);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_28, matches, true, true);
  if (count > 2) {
    fp.setBit(28);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_29, matches, true, true);
  if (count > 1) {
    fp.setBit(29);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_30, matches, true, true);
  if (count > 2) {
    fp.setBit(30);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_31, matches, true, true);
  if (count > 1) {
    fp.setBit(31);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_32, matches, true, true);
  if (count > 2) {
    fp.setBit(32);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_33, matches, true, true);
  if (count > 4) {
    fp.setBit(33);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_34, matches, true, true);
  if (count > 1) {
    fp.setBit(34);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_35, matches, true, true);
  if (count > 2) {
    fp.setBit(35);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_36, matches, true, true);
  if (count > 4) {
    fp.setBit(36);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_37, matches, true, true);
  if (count > 8) {
    fp.setBit(37);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_38, matches, true, true);
  if (count > 1) {
    fp.setBit(38);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_39, matches, true, true);
  if (count > 2) {
    fp.setBit(39);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_40, matches, true, true);
  if (count > 4) {
    fp.setBit(40);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_41, matches, true, true);
  if (count > 8) {
    fp.setBit(41);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_42, matches, true, true);
  if (count > 1) {
    fp.setBit(42);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_43, matches, true, true);
  if (count > 2) {
    fp.setBit(43);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_44, matches, true, true);
  if (count > 1) {
    fp.setBit(44);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_45, matches, true, true);
  if (count > 2) {
    fp.setBit(45);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_46, matches, true, true);
  if (count > 4) {
    fp.setBit(46);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_47, matches, true, true);
  if (count > 1) {
    fp.setBit(47);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_48, matches, true, true);
  if (count > 2) {
    fp.setBit(48);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_49, matches, true, true);
  if (count > 4) {
    fp.setBit(49);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_50, matches, true, true);
  if (count > 1) {
    fp.setBit(50);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_51, matches, true, true);
  if (count > 1) {
    fp.setBit(51);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_52, matches, true, true);
  if (count > 1) {
    fp.setBit(52);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_53, matches, true, true);
  if (count > 1) {
    fp.setBit(53);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_54, matches, true, true);
  if (count > 1) {
    fp.setBit(54);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_55, matches, true, true);
  if (count > 1) {
    fp.setBit(55);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_56, matches, true, true);
  if (count > 1) {
    fp.setBit(56);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_57, matches, true, true);
  if (count > 1) {
    fp.setBit(57);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_58, matches, true, true);
  if (count > 1) {
    fp.setBit(58);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_59, matches, true, true);
  if (count > 1) {
    fp.setBit(59);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_60, matches, true, true);
  if (count > 1) {
    fp.setBit(60);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_61, matches, true, true);
  if (count > 1) {
    fp.setBit(61);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_62, matches, true, true);
  if (count > 1) {
    fp.setBit(62);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_63, matches, true, true);
  if (count > 1) {
    fp.setBit(63);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_64, matches, true, true);
  if (count > 1) {
    fp.setBit(64);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_65, matches, true, true);
  if (count > 1) {
    fp.setBit(65);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_66, matches, true, true);
  if (count > 1) {
    fp.setBit(66);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_67, matches, true, true);
  if (count > 1) {
    fp.setBit(67);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_68, matches, true, true);
  if (count > 1) {
    fp.setBit(68);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_69, matches, true, true);
  if (count > 1) {
    fp.setBit(69);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_70, matches, true, true);
  if (count > 1) {
    fp.setBit(70);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_71, matches, true, true);
  if (count > 1) {
    fp.setBit(71);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_72, matches, true, true);
  if (count > 1) {
    fp.setBit(72);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_73, matches, true, true);
  if (count > 1) {
    fp.setBit(73);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_74, matches, true, true);
  if (count > 1) {
    fp.setBit(74);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_75, matches, true, true);
  if (count > 1) {
    fp.setBit(75);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_76, matches, true, true);
  if (count > 1) {
    fp.setBit(76);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_77, matches, true, true);
  if (count > 1) {
    fp.setBit(77);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_78, matches, true, true);
  if (count > 1) {
    fp.setBit(78);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_79, matches, true, true);
  if (count > 1) {
    fp.setBit(79);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_80, matches, true, true);
  if (count > 1) {
    fp.setBit(80);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_81, matches, true, true);
  if (count > 1) {
    fp.setBit(81);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_82, matches, true, true);
  if (count > 1) {
    fp.setBit(82);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_83, matches, true, true);
  if (count > 1) {
    fp.setBit(83);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_84, matches, true, true);
  if (count > 1) {
    fp.setBit(84);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_85, matches, true, true);
  if (count > 1) {
    fp.setBit(85);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_86, matches, true, true);
  if (count > 1) {
    fp.setBit(86);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_87, matches, true, true);
  if (count > 1) {
    fp.setBit(87);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_88, matches, true, true);
  if (count > 1) {
    fp.setBit(88);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_89, matches, true, true);
  if (count > 1) {
    fp.setBit(89);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_90, matches, true, true);
  if (count > 1) {
    fp.setBit(90);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_91, matches, true, true);
  if (count > 1) {
    fp.setBit(91);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_92, matches, true, true);
  if (count > 1) {
    fp.setBit(92);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_93, matches, true, true);
  if (count > 1) {
    fp.setBit(93);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_94, matches, true, true);
  if (count > 1) {
    fp.setBit(94);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_95, matches, true, true);
  if (count > 1) {
    fp.setBit(95);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_96, matches, true, true);
  if (count > 1) {
    fp.setBit(96);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_97, matches, true, true);
  if (count > 1) {
    fp.setBit(97);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_98, matches, true, true);
  if (count > 1) {
    fp.setBit(98);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_99, matches, true, true);
  if (count > 1) {
    fp.setBit(99);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_100, matches, true, true);
  if (count > 1) {
    fp.setBit(100);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_101, matches, true, true);
  if (count > 1) {
    fp.setBit(101);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_102, matches, true, true);
  if (count > 1) {
    fp.setBit(102);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_103, matches, true, true);
  if (count > 1) {
    fp.setBit(103);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_104, matches, true, true);
  if (count > 1) {
    fp.setBit(104);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_105, matches, true, true);
  if (count > 1) {
    fp.setBit(105);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_106, matches, true, true);
  if (count > 1) {
    fp.setBit(106);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_107, matches, true, true);
  if (count > 1) {
    fp.setBit(107);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_108, matches, true, true);
  if (count > 1) {
    fp.setBit(108);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_109, matches, true, true);
  if (count > 1) {
    fp.setBit(109);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_110, matches, true, true);
  if (count > 1) {
    fp.setBit(110);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_111, matches, true, true);
  if (count > 1) {
    fp.setBit(111);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_112, matches, true, true);
  if (count > 1) {
    fp.setBit(112);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_113, matches, true, true);
  if (count > 1) {
    fp.setBit(113);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_114, matches, true, true);
  if (count > 1) {
    fp.setBit(114);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_115, matches, true, true);
  if (count > 1) {
    fp.setBit(115);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_116, matches, true, true);
  if (count > 1) {
    fp.setBit(116);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_117, matches, true, true);
  if (count > 1) {
    fp.setBit(117);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_118, matches, true, true);
  if (count > 1) {
    fp.setBit(118);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_119, matches, true, true);
  if (count > 1) {
    fp.setBit(119);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_120, matches, true, true);
  if (count > 1) {
    fp.setBit(120);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_121, matches, true, true);
  if (count > 1) {
    fp.setBit(121);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_122, matches, true, true);
  if (count > 1) {
    fp.setBit(122);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_123, matches, true, true);
  if (count > 2) {
    fp.setBit(123);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_124, matches, true, true);
  if (count > 2) {
    fp.setBit(124);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_125, matches, true, true);
  if (count > 2) {
    fp.setBit(125);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_126, matches, true, true);
  if (count > 2) {
    fp.setBit(126);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_127, matches, true, true);
  if (count > 2) {
    fp.setBit(127);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_128, matches, true, true);
  if (count > 2) {
    fp.setBit(128);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_129, matches, true, true);
  if (count > 2) {
    fp.setBit(129);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_130, matches, true, true);
  if (count > 1) {
    fp.setBit(130);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_131, matches, true, true);
  if (count > 1) {
    fp.setBit(131);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_132, matches, true, true);
  if (count > 1) {
    fp.setBit(132);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_133, matches, true, true);
  if (count > 1) {
    fp.setBit(133);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_134, matches, true, true);
  if (count > 1) {
    fp.setBit(134);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_135, matches, true, true);
  if (count > 1) {
    fp.setBit(135);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_136, matches, true, true);
  if (count > 1) {
    fp.setBit(136);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_137, matches, true, true);
  if (count > 2) {
    fp.setBit(137);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_138, matches, true, true);
  if (count > 2) {
    fp.setBit(138);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_139, matches, true, true);
  if (count > 2) {
    fp.setBit(139);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_140, matches, true, true);
  if (count > 2) {
    fp.setBit(140);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_141, matches, true, true);
  if (count > 2) {
    fp.setBit(141);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_142, matches, true, true);
  if (count > 2) {
    fp.setBit(142);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_143, matches, true, true);
  if (count > 2) {
    fp.setBit(143);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_144, matches, true, true);
  if (count > 1) {
    fp.setBit(144);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_145, matches, true, true);
  if (count > 1) {
    fp.setBit(145);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_146, matches, true, true);
  if (count > 1) {
    fp.setBit(146);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_147, matches, true, true);
  if (count > 1) {
    fp.setBit(147);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_148, matches, true, true);
  if (count > 1) {
    fp.setBit(148);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_149, matches, true, true);
  if (count > 1) {
    fp.setBit(149);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_150, matches, true, true);
  if (count > 1) {
    fp.setBit(150);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_151, matches, true, true);
  if (count > 2) {
    fp.setBit(151);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_152, matches, true, true);
  if (count > 2) {
    fp.setBit(152);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_153, matches, true, true);
  if (count > 2) {
    fp.setBit(153);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_154, matches, true, true);
  if (count > 2) {
    fp.setBit(154);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_155, matches, true, true);
  if (count > 2) {
    fp.setBit(155);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_156, matches, true, true);
  if (count > 2) {
    fp.setBit(156);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_157, matches, true, true);
  if (count > 2) {
    fp.setBit(157);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_158, matches, true, true);
  if (count > 3) {
    fp.setBit(158);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_159, matches, true, true);
  if (count > 3) {
    fp.setBit(159);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_160, matches, true, true);
  if (count > 3) {
    fp.setBit(160);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_161, matches, true, true);
  if (count > 3) {
    fp.setBit(161);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_162, matches, true, true);
  if (count > 3) {
    fp.setBit(162);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_163, matches, true, true);
  if (count > 3) {
    fp.setBit(163);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_164, matches, true, true);
  if (count > 3) {
    fp.setBit(164);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_165, matches, true, true);
  if (count > 4) {
    fp.setBit(165);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_166, matches, true, true);
  if (count > 4) {
    fp.setBit(166);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_167, matches, true, true);
  if (count > 4) {
    fp.setBit(167);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_168, matches, true, true);
  if (count > 4) {
    fp.setBit(168);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_169, matches, true, true);
  if (count > 4) {
    fp.setBit(169);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_170, matches, true, true);
  if (count > 4) {
    fp.setBit(170);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_171, matches, true, true);
  if (count > 4) {
    fp.setBit(171);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_172, matches, true, true);
  if (count > 5) {
    fp.setBit(172);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_173, matches, true, true);
  if (count > 5) {
    fp.setBit(173);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_174, matches, true, true);
  if (count > 5) {
    fp.setBit(174);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_175, matches, true, true);
  if (count > 5) {
    fp.setBit(175);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_176, matches, true, true);
  if (count > 5) {
    fp.setBit(176);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_177, matches, true, true);
  if (count > 5) {
    fp.setBit(177);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_178, matches, true, true);
  if (count > 5) {
    fp.setBit(178);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_179, matches, true, true);
  if (count > 1) {
    fp.setBit(179);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_180, matches, true, true);
  if (count > 1) {
    fp.setBit(180);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_181, matches, true, true);
  if (count > 1) {
    fp.setBit(181);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_182, matches, true, true);
  if (count > 1) {
    fp.setBit(182);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_183, matches, true, true);
  if (count > 1) {
    fp.setBit(183);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_184, matches, true, true);
  if (count > 1) {
    fp.setBit(184);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_185, matches, true, true);
  if (count > 1) {
    fp.setBit(185);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_186, matches, true, true);
  if (count > 2) {
    fp.setBit(186);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_187, matches, true, true);
  if (count > 2) {
    fp.setBit(187);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_188, matches, true, true);
  if (count > 2) {
    fp.setBit(188);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_189, matches, true, true);
  if (count > 2) {
    fp.setBit(189);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_190, matches, true, true);
  if (count > 2) {
    fp.setBit(190);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_191, matches, true, true);
  if (count > 2) {
    fp.setBit(191);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_192, matches, true, true);
  if (count > 2) {
    fp.setBit(192);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_193, matches, true, true);
  if (count > 3) {
    fp.setBit(193);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_194, matches, true, true);
  if (count > 3) {
    fp.setBit(194);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_195, matches, true, true);
  if (count > 3) {
    fp.setBit(195);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_196, matches, true, true);
  if (count > 3) {
    fp.setBit(196);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_197, matches, true, true);
  if (count > 3) {
    fp.setBit(197);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_198, matches, true, true);
  if (count > 3) {
    fp.setBit(198);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_199, matches, true, true);
  if (count > 3) {
    fp.setBit(199);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_200, matches, true, true);
  if (count > 4) {
    fp.setBit(200);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_201, matches, true, true);
  if (count > 4) {
    fp.setBit(201);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_202, matches, true, true);
  if (count > 4) {
    fp.setBit(202);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_203, matches, true, true);
  if (count > 4) {
    fp.setBit(203);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_204, matches, true, true);
  if (count > 4) {
    fp.setBit(204);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_205, matches, true, true);
  if (count > 4) {
    fp.setBit(205);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_206, matches, true, true);
  if (count > 4) {
    fp.setBit(206);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_207, matches, true, true);
  if (count > 5) {
    fp.setBit(207);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_208, matches, true, true);
  if (count > 5) {
    fp.setBit(208);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_209, matches, true, true);
  if (count > 5) {
    fp.setBit(209);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_210, matches, true, true);
  if (count > 5) {
    fp.setBit(210);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_211, matches, true, true);
  if (count > 5) {
    fp.setBit(211);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_212, matches, true, true);
  if (count > 5) {
    fp.setBit(212);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_213, matches, true, true);
  if (count > 5) {
    fp.setBit(213);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_214, matches, true, true);
  if (count > 1) {
    fp.setBit(214);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_215, matches, true, true);
  if (count > 1) {
    fp.setBit(215);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_216, matches, true, true);
  if (count > 1) {
    fp.setBit(216);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_217, matches, true, true);
  if (count > 1) {
    fp.setBit(217);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_218, matches, true, true);
  if (count > 1) {
    fp.setBit(218);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_219, matches, true, true);
  if (count > 1) {
    fp.setBit(219);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_220, matches, true, true);
  if (count > 1) {
    fp.setBit(220);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_221, matches, true, true);
  if (count > 2) {
    fp.setBit(221);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_222, matches, true, true);
  if (count > 2) {
    fp.setBit(222);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_223, matches, true, true);
  if (count > 2) {
    fp.setBit(223);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_224, matches, true, true);
  if (count > 2) {
    fp.setBit(224);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_225, matches, true, true);
  if (count > 2) {
    fp.setBit(225);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_226, matches, true, true);
  if (count > 2) {
    fp.setBit(226);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_227, matches, true, true);
  if (count > 2) {
    fp.setBit(227);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_228, matches, true, true);
  if (count > 1) {
    fp.setBit(228);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_229, matches, true, true);
  if (count > 1) {
    fp.setBit(229);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_230, matches, true, true);
  if (count > 1) {
    fp.setBit(230);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_231, matches, true, true);
  if (count > 1) {
    fp.setBit(231);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_232, matches, true, true);
  if (count > 1) {
    fp.setBit(232);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_233, matches, true, true);
  if (count > 1) {
    fp.setBit(233);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_234, matches, true, true);
  if (count > 1) {
    fp.setBit(234);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_235, matches, true, true);
  if (count > 2) {
    fp.setBit(235);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_236, matches, true, true);
  if (count > 2) {
    fp.setBit(236);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_237, matches, true, true);
  if (count > 2) {
    fp.setBit(237);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_238, matches, true, true);
  if (count > 2) {
    fp.setBit(238);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_239, matches, true, true);
  if (count > 2) {
    fp.setBit(239);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_240, matches, true, true);
  if (count > 2) {
    fp.setBit(240);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_241, matches, true, true);
  if (count > 2) {
    fp.setBit(241);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_242, matches, true, true);
  if (count > 1) {
    fp.setBit(242);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_243, matches, true, true);
  if (count > 1) {
    fp.setBit(243);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_244, matches, true, true);
  if (count > 1) {
    fp.setBit(244);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_245, matches, true, true);
  if (count > 1) {
    fp.setBit(245);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_246, matches, true, true);
  if (count > 1) {
    fp.setBit(246);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_247, matches, true, true);
  if (count > 1) {
    fp.setBit(247);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_248, matches, true, true);
  if (count > 1) {
    fp.setBit(248);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_249, matches, true, true);
  if (count > 1) {
    fp.setBit(249);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_250, matches, true, true);
  if (count > 1) {
    fp.setBit(250);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_251, matches, true, true);
  if (count > 1) {
    fp.setBit(251);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_252, matches, true, true);
  if (count > 1) {
    fp.setBit(252);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_253, matches, true, true);
  if (count > 1) {
    fp.setBit(253);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_254, matches, true, true);
  if (count > 1) {
    fp.setBit(254);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_255, matches, true, true);
  if (count > 1) {
    fp.setBit(255);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_256, matches, true, true);
  if (count > 1) {
    fp.setBit(256);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_257, matches, true, true);
  if (count > 1) {
    fp.setBit(257);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_258, matches, true, true);
  if (count > 2) {
    fp.setBit(258);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_259, matches, true, true);
  if (count > 2) {
    fp.setBit(259);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_260, matches, true, true);
  if (count > 3) {
    fp.setBit(260);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_261, matches, true, true);
  if (count > 3) {
    fp.setBit(261);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_262, matches, true, true);
  if (count > 4) {
    fp.setBit(262);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_263, matches, true, true);
  if (count > 4) {
    fp.setBit(263);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_264, matches, true, true);
  if (count > 1) {
    fp.setBit(264);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_265, matches, true, true);
  if (count > 1) {
    fp.setBit(265);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_266, matches, true, true);
  if (count > 1) {
    fp.setBit(266);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_267, matches, true, true);
  if (count > 1) {
    fp.setBit(267);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_268, matches, true, true);
  if (count > 1) {
    fp.setBit(268);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_269, matches, true, true);
  if (count > 1) {
    fp.setBit(269);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_270, matches, true, true);
  if (count > 1) {
    fp.setBit(270);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_271, matches, true, true);
  if (count > 1) {
    fp.setBit(271);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_272, matches, true, true);
  if (count > 1) {
    fp.setBit(272);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_273, matches, true, true);
  if (count > 1) {
    fp.setBit(273);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_274, matches, true, true);
  if (count > 1) {
    fp.setBit(274);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_275, matches, true, true);
  if (count > 1) {
    fp.setBit(275);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_276, matches, true, true);
  if (count > 1) {
    fp.setBit(276);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_277, matches, true, true);
  if (count > 1) {
    fp.setBit(277);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_278, matches, true, true);
  if (count > 1) {
    fp.setBit(278);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_279, matches, true, true);
  if (count > 1) {
    fp.setBit(279);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_280, matches, true, true);
  if (count > 1) {
    fp.setBit(280);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_281, matches, true, true);
  if (count > 1) {
    fp.setBit(281);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_282, matches, true, true);
  if (count > 1) {
    fp.setBit(282);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_283, matches, true, true);
  if (count > 1) {
    fp.setBit(283);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_284, matches, true, true);
  if (count > 1) {
    fp.setBit(284);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_285, matches, true, true);
  if (count > 1) {
    fp.setBit(285);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_286, matches, true, true);
  if (count > 1) {
    fp.setBit(286);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_287, matches, true, true);
  if (count > 1) {
    fp.setBit(287);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_288, matches, true, true);
  if (count > 1) {
    fp.setBit(288);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_289, matches, true, true);
  if (count > 1) {
    fp.setBit(289);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_290, matches, true, true);
  if (count > 1) {
    fp.setBit(290);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_291, matches, true, true);
  if (count > 1) {
    fp.setBit(291);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_292, matches, true, true);
  if (count > 1) {
    fp.setBit(292);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_293, matches, true, true);
  if (count > 1) {
    fp.setBit(293);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_294, matches, true, true);
  if (count > 1) {
    fp.setBit(294);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_295, matches, true, true);
  if (count > 1) {
    fp.setBit(295);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_296, matches, true, true);
  if (count > 1) {
    fp.setBit(296);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_297, matches, true, true);
  if (count > 1) {
    fp.setBit(297);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_298, matches, true, true);
  if (count > 1) {
    fp.setBit(298);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_299, matches, true, true);
  if (count > 1) {
    fp.setBit(299);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_300, matches, true, true);
  if (count > 1) {
    fp.setBit(300);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_301, matches, true, true);
  if (count > 1) {
    fp.setBit(301);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_302, matches, true, true);
  if (count > 1) {
    fp.setBit(302);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_303, matches, true, true);
  if (count > 1) {
    fp.setBit(303);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_304, matches, true, true);
  if (count > 1) {
    fp.setBit(304);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_305, matches, true, true);
  if (count > 1) {
    fp.setBit(305);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_306, matches, true, true);
  if (count > 1) {
    fp.setBit(306);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_307, matches, true, true);
  if (count > 1) {
    fp.setBit(307);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_308, matches, true, true);
  if (count > 1) {
    fp.setBit(308);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_309, matches, true, true);
  if (count > 1) {
    fp.setBit(309);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_310, matches, true, true);
  if (count > 1) {
    fp.setBit(310);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_311, matches, true, true);
  if (count > 1) {
    fp.setBit(311);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_312, matches, true, true);
  if (count > 1) {
    fp.setBit(312);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_313, matches, true, true);
  if (count > 1) {
    fp.setBit(313);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_314, matches, true, true);
  if (count > 1) {
    fp.setBit(314);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_315, matches, true, true);
  if (count > 1) {
    fp.setBit(315);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_316, matches, true, true);
  if (count > 1) {
    fp.setBit(316);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_317, matches, true, true);
  if (count > 1) {
    fp.setBit(317);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_318, matches, true, true);
  if (count > 1) {
    fp.setBit(318);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_319, matches, true, true);
  if (count > 1) {
    fp.setBit(319);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_320, matches, true, true);
  if (count > 1) {
    fp.setBit(320);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_321, matches, true, true);
  if (count > 1) {
    fp.setBit(321);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_322, matches, true, true);
  if (count > 1) {
    fp.setBit(322);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_323, matches, true, true);
  if (count > 1) {
    fp.setBit(323);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_324, matches, true, true);
  if (count > 1) {
    fp.setBit(324);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_325, matches, true, true);
  if (count > 1) {
    fp.setBit(325);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_326, matches, true, true);
  if (count > 1) {
    fp.setBit(326);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_327, matches, true, true);
  if (count > 1) {
    fp.setBit(327);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_328, matches, true, true);
  if (count > 1) {
    fp.setBit(328);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_329, matches, true, true);
  if (count > 1) {
    fp.setBit(329);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_330, matches, true, true);
  if (count > 1) {
    fp.setBit(330);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_331, matches, true, true);
  if (count > 1) {
    fp.setBit(331);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_332, matches, true, true);
  if (count > 1) {
    fp.setBit(332);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_333, matches, true, true);
  if (count > 1) {
    fp.setBit(333);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_334, matches, true, true);
  if (count > 1) {
    fp.setBit(334);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_335, matches, true, true);
  if (count > 1) {
    fp.setBit(335);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_336, matches, true, true);
  if (count > 1) {
    fp.setBit(336);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_337, matches, true, true);
  if (count > 1) {
    fp.setBit(337);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_338, matches, true, true);
  if (count > 1) {
    fp.setBit(338);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_339, matches, true, true);
  if (count > 1) {
    fp.setBit(339);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_340, matches, true, true);
  if (count > 1) {
    fp.setBit(340);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_341, matches, true, true);
  if (count > 1) {
    fp.setBit(341);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_342, matches, true, true);
  if (count > 1) {
    fp.setBit(342);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_343, matches, true, true);
  if (count > 1) {
    fp.setBit(343);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_344, matches, true, true);
  if (count > 1) {
    fp.setBit(344);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_345, matches, true, true);
  if (count > 1) {
    fp.setBit(345);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_346, matches, true, true);
  if (count > 1) {
    fp.setBit(346);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_347, matches, true, true);
  if (count > 1) {
    fp.setBit(347);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_348, matches, true, true);
  if (count > 1) {
    fp.setBit(348);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_349, matches, true, true);
  if (count > 1) {
    fp.setBit(349);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_350, matches, true, true);
  if (count > 1) {
    fp.setBit(350);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_351, matches, true, true);
  if (count > 1) {
    fp.setBit(351);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_352, matches, true, true);
  if (count > 1) {
    fp.setBit(352);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_353, matches, true, true);
  if (count > 1) {
    fp.setBit(353);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_354, matches, true, true);
  if (count > 1) {
    fp.setBit(354);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_355, matches, true, true);
  if (count > 1) {
    fp.setBit(355);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_356, matches, true, true);
  if (count > 1) {
    fp.setBit(356);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_357, matches, true, true);
  if (count > 1) {
    fp.setBit(357);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_358, matches, true, true);
  if (count > 1) {
    fp.setBit(358);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_359, matches, true, true);
  if (count > 1) {
    fp.setBit(359);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_360, matches, true, true);
  if (count > 1) {
    fp.setBit(360);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_361, matches, true, true);
  if (count > 1) {
    fp.setBit(361);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_362, matches, true, true);
  if (count > 1) {
    fp.setBit(362);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_363, matches, true, true);
  if (count > 1) {
    fp.setBit(363);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_364, matches, true, true);
  if (count > 1) {
    fp.setBit(364);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_365, matches, true, true);
  if (count > 1) {
    fp.setBit(365);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_366, matches, true, true);
  if (count > 1) {
    fp.setBit(366);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_367, matches, true, true);
  if (count > 1) {
    fp.setBit(367);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_368, matches, true, true);
  if (count > 1) {
    fp.setBit(368);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_369, matches, true, true);
  if (count > 1) {
    fp.setBit(369);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_370, matches, true, true);
  if (count > 1) {
    fp.setBit(370);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_371, matches, true, true);
  if (count > 1) {
    fp.setBit(371);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_372, matches, true, true);
  if (count > 1) {
    fp.setBit(372);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_373, matches, true, true);
  if (count > 1) {
    fp.setBit(373);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_374, matches, true, true);
  if (count > 1) {
    fp.setBit(374);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_375, matches, true, true);
  if (count > 1) {
    fp.setBit(375);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_376, matches, true, true);
  if (count > 1) {
    fp.setBit(376);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_377, matches, true, true);
  if (count > 1) {
    fp.setBit(377);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_378, matches, true, true);
  if (count > 1) {
    fp.setBit(378);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_379, matches, true, true);
  if (count > 1) {
    fp.setBit(379);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_380, matches, true, true);
  if (count > 1) {
    fp.setBit(380);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_381, matches, true, true);
  if (count > 1) {
    fp.setBit(381);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_382, matches, true, true);
  if (count > 1) {
    fp.setBit(382);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_383, matches, true, true);
  if (count > 1) {
    fp.setBit(383);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_384, matches, true, true);
  if (count > 1) {
    fp.setBit(384);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_385, matches, true, true);
  if (count > 1) {
    fp.setBit(385);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_386, matches, true, true);
  if (count > 1) {
    fp.setBit(386);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_387, matches, true, true);
  if (count > 1) {
    fp.setBit(387);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_388, matches, true, true);
  if (count > 1) {
    fp.setBit(388);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_389, matches, true, true);
  if (count > 1) {
    fp.setBit(389);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_390, matches, true, true);
  if (count > 1) {
    fp.setBit(390);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_391, matches, true, true);
  if (count > 1) {
    fp.setBit(391);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_392, matches, true, true);
  if (count > 1) {
    fp.setBit(392);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_393, matches, true, true);
  if (count > 1) {
    fp.setBit(393);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_394, matches, true, true);
  if (count > 1) {
    fp.setBit(394);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_395, matches, true, true);
  if (count > 1) {
    fp.setBit(395);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_396, matches, true, true);
  if (count > 1) {
    fp.setBit(396);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_397, matches, true, true);
  if (count > 1) {
    fp.setBit(397);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_398, matches, true, true);
  if (count > 1) {
    fp.setBit(398);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_399, matches, true, true);
  if (count > 1) {
    fp.setBit(399);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_400, matches, true, true);
  if (count > 1) {
    fp.setBit(400);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_401, matches, true, true);
  if (count > 1) {
    fp.setBit(401);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_402, matches, true, true);
  if (count > 1) {
    fp.setBit(402);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_403, matches, true, true);
  if (count > 1) {
    fp.setBit(403);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_404, matches, true, true);
  if (count > 1) {
    fp.setBit(404);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_405, matches, true, true);
  if (count > 1) {
    fp.setBit(405);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_406, matches, true, true);
  if (count > 1) {
    fp.setBit(406);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_407, matches, true, true);
  if (count > 1) {
    fp.setBit(407);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_408, matches, true, true);
  if (count > 1) {
    fp.setBit(408);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_409, matches, true, true);
  if (count > 1) {
    fp.setBit(409);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_410, matches, true, true);
  if (count > 1) {
    fp.setBit(410);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_411, matches, true, true);
  if (count > 1) {
    fp.setBit(411);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_412, matches, true, true);
  if (count > 1) {
    fp.setBit(412);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_413, matches, true, true);
  if (count > 1) {
    fp.setBit(413);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_414, matches, true, true);
  if (count > 1) {
    fp.setBit(414);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_415, matches, true, true);
  if (count > 1) {
    fp.setBit(415);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_416, matches, true, true);
  if (count > 1) {
    fp.setBit(416);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_417, matches, true, true);
  if (count > 1) {
    fp.setBit(417);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_418, matches, true, true);
  if (count > 1) {
    fp.setBit(418);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_419, matches, true, true);
  if (count > 1) {
    fp.setBit(419);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_420, matches, true, true);
  if (count > 1) {
    fp.setBit(420);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_421, matches, true, true);
  if (count > 1) {
    fp.setBit(421);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_422, matches, true, true);
  if (count > 1) {
    fp.setBit(422);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_423, matches, true, true);
  if (count > 1) {
    fp.setBit(423);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_424, matches, true, true);
  if (count > 1) {
    fp.setBit(424);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_425, matches, true, true);
  if (count > 1) {
    fp.setBit(425);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_426, matches, true, true);
  if (count > 1) {
    fp.setBit(426);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_427, matches, true, true);
  if (count > 1) {
    fp.setBit(427);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_428, matches, true, true);
  if (count > 1) {
    fp.setBit(428);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_429, matches, true, true);
  if (count > 1) {
    fp.setBit(429);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_430, matches, true, true);
  if (count > 1) {
    fp.setBit(430);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_431, matches, true, true);
  if (count > 1) {
    fp.setBit(431);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_432, matches, true, true);
  if (count > 1) {
    fp.setBit(432);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_433, matches, true, true);
  if (count > 1) {
    fp.setBit(433);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_434, matches, true, true);
  if (count > 1) {
    fp.setBit(434);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_435, matches, true, true);
  if (count > 1) {
    fp.setBit(435);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_436, matches, true, true);
  if (count > 1) {
    fp.setBit(436);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_437, matches, true, true);
  if (count > 1) {
    fp.setBit(437);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_438, matches, true, true);
  if (count > 1) {
    fp.setBit(438);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_439, matches, true, true);
  if (count > 1) {
    fp.setBit(439);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_440, matches, true, true);
  if (count > 1) {
    fp.setBit(440);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_441, matches, true, true);
  if (count > 1) {
    fp.setBit(441);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_442, matches, true, true);
  if (count > 1) {
    fp.setBit(442);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_443, matches, true, true);
  if (count > 1) {
    fp.setBit(443);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_444, matches, true, true);
  if (count > 1) {
    fp.setBit(444);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_445, matches, true, true);
  if (count > 1) {
    fp.setBit(445);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_446, matches, true, true);
  if (count > 1) {
    fp.setBit(446);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_447, matches, true, true);
  if (count > 1) {
    fp.setBit(447);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_448, matches, true, true);
  if (count > 1) {
    fp.setBit(448);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_449, matches, true, true);
  if (count > 1) {
    fp.setBit(449);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_450, matches, true, true);
  if (count > 1) {
    fp.setBit(450);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_451, matches, true, true);
  if (count > 1) {
    fp.setBit(451);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_452, matches, true, true);
  if (count > 1) {
    fp.setBit(452);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_453, matches, true, true);
  if (count > 1) {
    fp.setBit(453);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_454, matches, true, true);
  if (count > 1) {
    fp.setBit(454);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_455, matches, true, true);
  if (count > 1) {
    fp.setBit(455);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_456, matches, true, true);
  if (count > 1) {
    fp.setBit(456);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_457, matches, true, true);
  if (count > 1) {
    fp.setBit(457);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_458, matches, true, true);
  if (count > 1) {
    fp.setBit(458);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_459, matches, true, true);
  if (count > 1) {
    fp.setBit(459);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_460, matches, true, true);
  if (count > 1) {
    fp.setBit(460);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_461, matches, true, true);
  if (count > 1) {
    fp.setBit(461);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_462, matches, true, true);
  if (count > 1) {
    fp.setBit(462);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_463, matches, true, true);
  if (count > 1) {
    fp.setBit(463);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_464, matches, true, true);
  if (count > 1) {
    fp.setBit(464);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_465, matches, true, true);
  if (count > 1) {
    fp.setBit(465);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_466, matches, true, true);
  if (count > 1) {
    fp.setBit(466);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_467, matches, true, true);
  if (count > 1) {
    fp.setBit(467);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_468, matches, true, true);
  if (count > 1) {
    fp.setBit(468);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_469, matches, true, true);
  if (count > 1) {
    fp.setBit(469);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_470, matches, true, true);
  if (count > 1) {
    fp.setBit(470);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_471, matches, true, true);
  if (count > 1) {
    fp.setBit(471);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_472, matches, true, true);
  if (count > 1) {
    fp.setBit(472);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_473, matches, true, true);
  if (count > 1) {
    fp.setBit(473);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_474, matches, true, true);
  if (count > 1) {
    fp.setBit(474);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_475, matches, true, true);
  if (count > 1) {
    fp.setBit(475);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_476, matches, true, true);
  if (count > 1) {
    fp.setBit(476);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_477, matches, true, true);
  if (count > 1) {
    fp.setBit(477);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_478, matches, true, true);
  if (count > 1) {
    fp.setBit(478);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_479, matches, true, true);
  if (count > 1) {
    fp.setBit(479);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_480, matches, true, true);
  if (count > 1) {
    fp.setBit(480);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_481, matches, true, true);
  if (count > 1) {
    fp.setBit(481);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_482, matches, true, true);
  if (count > 1) {
    fp.setBit(482);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_483, matches, true, true);
  if (count > 1) {
    fp.setBit(483);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_484, matches, true, true);
  if (count > 1) {
    fp.setBit(484);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_485, matches, true, true);
  if (count > 1) {
    fp.setBit(485);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_486, matches, true, true);
  if (count > 1) {
    fp.setBit(486);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_487, matches, true, true);
  if (count > 1) {
    fp.setBit(487);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_488, matches, true, true);
  if (count > 1) {
    fp.setBit(488);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_489, matches, true, true);
  if (count > 1) {
    fp.setBit(489);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_490, matches, true, true);
  if (count > 1) {
    fp.setBit(490);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_491, matches, true, true);
  if (count > 1) {
    fp.setBit(491);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_492, matches, true, true);
  if (count > 1) {
    fp.setBit(492);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_493, matches, true, true);
  if (count > 1) {
    fp.setBit(493);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_494, matches, true, true);
  if (count > 1) {
    fp.setBit(494);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_495, matches, true, true);
  if (count > 1) {
    fp.setBit(495);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_496, matches, true, true);
  if (count > 1) {
    fp.setBit(496);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_497, matches, true, true);
  if (count > 1) {
    fp.setBit(497);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_498, matches, true, true);
  if (count > 1) {
    fp.setBit(498);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_499, matches, true, true);
  if (count > 1) {
    fp.setBit(499);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_500, matches, true, true);
  if (count > 1) {
    fp.setBit(500);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_501, matches, true, true);
  if (count > 1) {
    fp.setBit(501);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_502, matches, true, true);
  if (count > 1) {
    fp.setBit(502);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_503, matches, true, true);
  if (count > 1) {
    fp.setBit(503);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_504, matches, true, true);
  if (count > 1) {
    fp.setBit(504);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_505, matches, true, true);
  if (count > 1) {
    fp.setBit(505);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_506, matches, true, true);
  if (count > 1) {
    fp.setBit(506);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_507, matches, true, true);
  if (count > 1) {
    fp.setBit(507);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_508, matches, true, true);
  if (count > 1) {
    fp.setBit(508);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_509, matches, true, true);
  if (count > 1) {
    fp.setBit(509);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_510, matches, true, true);
  if (count > 1) {
    fp.setBit(510);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_511, matches, true, true);
  if (count > 1) {
    fp.setBit(511);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_512, matches, true, true);
  if (count > 1) {
    fp.setBit(512);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_513, matches, true, true);
  if (count > 1) {
    fp.setBit(513);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_514, matches, true, true);
  if (count > 1) {
    fp.setBit(514);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_515, matches, true, true);
  if (count > 1) {
    fp.setBit(515);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_516, matches, true, true);
  if (count > 1) {
    fp.setBit(516);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_517, matches, true, true);
  if (count > 1) {
    fp.setBit(517);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_518, matches, true, true);
  if (count > 1) {
    fp.setBit(518);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_519, matches, true, true);
  if (count > 1) {
    fp.setBit(519);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_520, matches, true, true);
  if (count > 1) {
    fp.setBit(520);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_521, matches, true, true);
  if (count > 1) {
    fp.setBit(521);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_522, matches, true, true);
  if (count > 1) {
    fp.setBit(522);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_523, matches, true, true);
  if (count > 1) {
    fp.setBit(523);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_524, matches, true, true);
  if (count > 1) {
    fp.setBit(524);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_525, matches, true, true);
  if (count > 1) {
    fp.setBit(525);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_526, matches, true, true);
  if (count > 1) {
    fp.setBit(526);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_527, matches, true, true);
  if (count > 1) {
    fp.setBit(527);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_528, matches, true, true);
  if (count > 1) {
    fp.setBit(528);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_529, matches, true, true);
  if (count > 1) {
    fp.setBit(529);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_530, matches, true, true);
  if (count > 1) {
    fp.setBit(530);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_531, matches, true, true);
  if (count > 1) {
    fp.setBit(531);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_532, matches, true, true);
  if (count > 1) {
    fp.setBit(532);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_533, matches, true, true);
  if (count > 1) {
    fp.setBit(533);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_534, matches, true, true);
  if (count > 1) {
    fp.setBit(534);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_535, matches, true, true);
  if (count > 1) {
    fp.setBit(535);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_536, matches, true, true);
  if (count > 1) {
    fp.setBit(536);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_537, matches, true, true);
  if (count > 1) {
    fp.setBit(537);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_538, matches, true, true);
  if (count > 1) {
    fp.setBit(538);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_539, matches, true, true);
  if (count > 1) {
    fp.setBit(539);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_540, matches, true, true);
  if (count > 1) {
    fp.setBit(540);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_541, matches, true, true);
  if (count > 1) {
    fp.setBit(541);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_542, matches, true, true);
  if (count > 1) {
    fp.setBit(542);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_543, matches, true, true);
  if (count > 1) {
    fp.setBit(543);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_544, matches, true, true);
  if (count > 1) {
    fp.setBit(544);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_545, matches, true, true);
  if (count > 1) {
    fp.setBit(545);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_546, matches, true, true);
  if (count > 1) {
    fp.setBit(546);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_547, matches, true, true);
  if (count > 1) {
    fp.setBit(547);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_548, matches, true, true);
  if (count > 1) {
    fp.setBit(548);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_549, matches, true, true);
  if (count > 1) {
    fp.setBit(549);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_550, matches, true, true);
  if (count > 1) {
    fp.setBit(550);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_551, matches, true, true);
  if (count > 1) {
    fp.setBit(551);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_552, matches, true, true);
  if (count > 1) {
    fp.setBit(552);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_553, matches, true, true);
  if (count > 1) {
    fp.setBit(553);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_554, matches, true, true);
  if (count > 1) {
    fp.setBit(554);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_555, matches, true, true);
  if (count > 1) {
    fp.setBit(555);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_556, matches, true, true);
  if (count > 1) {
    fp.setBit(556);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_557, matches, true, true);
  if (count > 1) {
    fp.setBit(557);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_558, matches, true, true);
  if (count > 1) {
    fp.setBit(558);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_559, matches, true, true);
  if (count > 1) {
    fp.setBit(559);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_560, matches, true, true);
  if (count > 1) {
    fp.setBit(560);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_561, matches, true, true);
  if (count > 1) {
    fp.setBit(561);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_562, matches, true, true);
  if (count > 1) {
    fp.setBit(562);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_563, matches, true, true);
  if (count > 1) {
    fp.setBit(563);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_564, matches, true, true);
  if (count > 1) {
    fp.setBit(564);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_565, matches, true, true);
  if (count > 1) {
    fp.setBit(565);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_566, matches, true, true);
  if (count > 1) {
    fp.setBit(566);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_567, matches, true, true);
  if (count > 1) {
    fp.setBit(567);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_568, matches, true, true);
  if (count > 1) {
    fp.setBit(568);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_569, matches, true, true);
  if (count > 1) {
    fp.setBit(569);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_570, matches, true, true);
  if (count > 1) {
    fp.setBit(570);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_571, matches, true, true);
  if (count > 1) {
    fp.setBit(571);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_572, matches, true, true);
  if (count > 1) {
    fp.setBit(572);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_573, matches, true, true);
  if (count > 1) {
    fp.setBit(573);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_574, matches, true, true);
  if (count > 1) {
    fp.setBit(574);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_575, matches, true, true);
  if (count > 1) {
    fp.setBit(575);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_576, matches, true, true);
  if (count > 1) {
    fp.setBit(576);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_577, matches, true, true);
  if (count > 1) {
    fp.setBit(577);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_578, matches, true, true);
  if (count > 1) {
    fp.setBit(578);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_579, matches, true, true);
  if (count > 1) {
    fp.setBit(579);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_580, matches, true, true);
  if (count > 1) {
    fp.setBit(580);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_581, matches, true, true);
  if (count > 1) {
    fp.setBit(581);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_582, matches, true, true);
  if (count > 1) {
    fp.setBit(582);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_583, matches, true, true);
  if (count > 1) {
    fp.setBit(583);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_584, matches, true, true);
  if (count > 1) {
    fp.setBit(584);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_585, matches, true, true);
  if (count > 1) {
    fp.setBit(585);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_586, matches, true, true);
  if (count > 1) {
    fp.setBit(586);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_587, matches, true, true);
  if (count > 1) {
    fp.setBit(587);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_588, matches, true, true);
  if (count > 1) {
    fp.setBit(588);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_589, matches, true, true);
  if (count > 1) {
    fp.setBit(589);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_590, matches, true, true);
  if (count > 1) {
    fp.setBit(590);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_591, matches, true, true);
  if (count > 1) {
    fp.setBit(591);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_592, matches, true, true);
  if (count > 1) {
    fp.setBit(592);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_593, matches, true, true);
  if (count > 1) {
    fp.setBit(593);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_594, matches, true, true);
  if (count > 1) {
    fp.setBit(594);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_595, matches, true, true);
  if (count > 1) {
    fp.setBit(595);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_596, matches, true, true);
  if (count > 1) {
    fp.setBit(596);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_597, matches, true, true);
  if (count > 1) {
    fp.setBit(597);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_598, matches, true, true);
  if (count > 1) {
    fp.setBit(598);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_599, matches, true, true);
  if (count > 1) {
    fp.setBit(599);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_600, matches, true, true);
  if (count > 1) {
    fp.setBit(600);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_601, matches, true, true);
  if (count > 1) {
    fp.setBit(601);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_602, matches, true, true);
  if (count > 1) {
    fp.setBit(602);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_603, matches, true, true);
  if (count > 1) {
    fp.setBit(603);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_604, matches, true, true);
  if (count > 1) {
    fp.setBit(604);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_605, matches, true, true);
  if (count > 1) {
    fp.setBit(605);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_606, matches, true, true);
  if (count > 1) {
    fp.setBit(606);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_607, matches, true, true);
  if (count > 1) {
    fp.setBit(607);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_608, matches, true, true);
  if (count > 1) {
    fp.setBit(608);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_609, matches, true, true);
  if (count > 1) {
    fp.setBit(609);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_610, matches, true, true);
  if (count > 1) {
    fp.setBit(610);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_611, matches, true, true);
  if (count > 1) {
    fp.setBit(611);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_612, matches, true, true);
  if (count > 1) {
    fp.setBit(612);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_613, matches, true, true);
  if (count > 1) {
    fp.setBit(613);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_614, matches, true, true);
  if (count > 1) {
    fp.setBit(614);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_615, matches, true, true);
  if (count > 1) {
    fp.setBit(615);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_616, matches, true, true);
  if (count > 1) {
    fp.setBit(616);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_617, matches, true, true);
  if (count > 1) {
    fp.setBit(617);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_618, matches, true, true);
  if (count > 1) {
    fp.setBit(618);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_619, matches, true, true);
  if (count > 1) {
    fp.setBit(619);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_620, matches, true, true);
  if (count > 1) {
    fp.setBit(620);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_621, matches, true, true);
  if (count > 1) {
    fp.setBit(621);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_622, matches, true, true);
  if (count > 1) {
    fp.setBit(622);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_623, matches, true, true);
  if (count > 1) {
    fp.setBit(623);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_624, matches, true, true);
  if (count > 1) {
    fp.setBit(624);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_625, matches, true, true);
  if (count > 1) {
    fp.setBit(625);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_626, matches, true, true);
  if (count > 1) {
    fp.setBit(626);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_627, matches, true, true);
  if (count > 1) {
    fp.setBit(627);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_628, matches, true, true);
  if (count > 1) {
    fp.setBit(628);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_629, matches, true, true);
  if (count > 1) {
    fp.setBit(629);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_630, matches, true, true);
  if (count > 1) {
    fp.setBit(630);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_631, matches, true, true);
  if (count > 1) {
    fp.setBit(631);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_632, matches, true, true);
  if (count > 1) {
    fp.setBit(632);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_633, matches, true, true);
  if (count > 1) {
    fp.setBit(633);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_634, matches, true, true);
  if (count > 1) {
    fp.setBit(634);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_635, matches, true, true);
  if (count > 1) {
    fp.setBit(635);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_636, matches, true, true);
  if (count > 1) {
    fp.setBit(636);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_637, matches, true, true);
  if (count > 1) {
    fp.setBit(637);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_638, matches, true, true);
  if (count > 1) {
    fp.setBit(638);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_639, matches, true, true);
  if (count > 1) {
    fp.setBit(639);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_640, matches, true, true);
  if (count > 1) {
    fp.setBit(640);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_641, matches, true, true);
  if (count > 1) {
    fp.setBit(641);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_642, matches, true, true);
  if (count > 1) {
    fp.setBit(642);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_643, matches, true, true);
  if (count > 1) {
    fp.setBit(643);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_644, matches, true, true);
  if (count > 1) {
    fp.setBit(644);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_645, matches, true, true);
  if (count > 1) {
    fp.setBit(645);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_646, matches, true, true);
  if (count > 1) {
    fp.setBit(646);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_647, matches, true, true);
  if (count > 1) {
    fp.setBit(647);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_648, matches, true, true);
  if (count > 1) {
    fp.setBit(648);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_649, matches, true, true);
  if (count > 1) {
    fp.setBit(649);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_650, matches, true, true);
  if (count > 1) {
    fp.setBit(650);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_651, matches, true, true);
  if (count > 1) {
    fp.setBit(651);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_652, matches, true, true);
  if (count > 1) {
    fp.setBit(652);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_653, matches, true, true);
  if (count > 1) {
    fp.setBit(653);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_654, matches, true, true);
  if (count > 1) {
    fp.setBit(654);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_655, matches, true, true);
  if (count > 1) {
    fp.setBit(655);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_656, matches, true, true);
  if (count > 1) {
    fp.setBit(656);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_657, matches, true, true);
  if (count > 1) {
    fp.setBit(657);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_658, matches, true, true);
  if (count > 1) {
    fp.setBit(658);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_659, matches, true, true);
  if (count > 1) {
    fp.setBit(659);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_660, matches, true, true);
  if (count > 1) {
    fp.setBit(660);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_661, matches, true, true);
  if (count > 1) {
    fp.setBit(661);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_662, matches, true, true);
  if (count > 1) {
    fp.setBit(662);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_663, matches, true, true);
  if (count > 1) {
    fp.setBit(663);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_664, matches, true, true);
  if (count > 1) {
    fp.setBit(664);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_665, matches, true, true);
  if (count > 1) {
    fp.setBit(665);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_666, matches, true, true);
  if (count > 1) {
    fp.setBit(666);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_667, matches, true, true);
  if (count > 1) {
    fp.setBit(667);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_668, matches, true, true);
  if (count > 1) {
    fp.setBit(668);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_669, matches, true, true);
  if (count > 1) {
    fp.setBit(669);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_670, matches, true, true);
  if (count > 1) {
    fp.setBit(670);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_671, matches, true, true);
  if (count > 1) {
    fp.setBit(671);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_672, matches, true, true);
  if (count > 1) {
    fp.setBit(672);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_673, matches, true, true);
  if (count > 1) {
    fp.setBit(673);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_674, matches, true, true);
  if (count > 1) {
    fp.setBit(674);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_675, matches, true, true);
  if (count > 1) {
    fp.setBit(675);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_676, matches, true, true);
  if (count > 1) {
    fp.setBit(676);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_677, matches, true, true);
  if (count > 1) {
    fp.setBit(677);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_678, matches, true, true);
  if (count > 1) {
    fp.setBit(678);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_679, matches, true, true);
  if (count > 1) {
    fp.setBit(679);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_680, matches, true, true);
  if (count > 1) {
    fp.setBit(680);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_681, matches, true, true);
  if (count > 1) {
    fp.setBit(681);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_682, matches, true, true);
  if (count > 1) {
    fp.setBit(682);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_683, matches, true, true);
  if (count > 1) {
    fp.setBit(683);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_684, matches, true, true);
  if (count > 1) {
    fp.setBit(684);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_685, matches, true, true);
  if (count > 1) {
    fp.setBit(685);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_686, matches, true, true);
  if (count > 1) {
    fp.setBit(686);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_687, matches, true, true);
  if (count > 1) {
    fp.setBit(687);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_688, matches, true, true);
  if (count > 1) {
    fp.setBit(688);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_689, matches, true, true);
  if (count > 1) {
    fp.setBit(689);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_690, matches, true, true);
  if (count > 1) {
    fp.setBit(690);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_691, matches, true, true);
  if (count > 1) {
    fp.setBit(691);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_692, matches, true, true);
  if (count > 1) {
    fp.setBit(692);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_693, matches, true, true);
  if (count > 1) {
    fp.setBit(693);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_694, matches, true, true);
  if (count > 1) {
    fp.setBit(694);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_695, matches, true, true);
  if (count > 1) {
    fp.setBit(695);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_696, matches, true, true);
  if (count > 1) {
    fp.setBit(696);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_697, matches, true, true);
  if (count > 1) {
    fp.setBit(697);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_698, matches, true, true);
  if (count > 1) {
    fp.setBit(698);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_699, matches, true, true);
  if (count > 1) {
    fp.setBit(699);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_700, matches, true, true);
  if (count > 1) {
    fp.setBit(700);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_701, matches, true, true);
  if (count > 1) {
    fp.setBit(701);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_702, matches, true, true);
  if (count > 1) {
    fp.setBit(702);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_703, matches, true, true);
  if (count > 1) {
    fp.setBit(703);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_704, matches, true, true);
  if (count > 1) {
    fp.setBit(704);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_705, matches, true, true);
  if (count > 1) {
    fp.setBit(705);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_706, matches, true, true);
  if (count > 1) {
    fp.setBit(706);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_707, matches, true, true);
  if (count > 1) {
    fp.setBit(707);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_708, matches, true, true);
  if (count > 1) {
    fp.setBit(708);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_709, matches, true, true);
  if (count > 1) {
    fp.setBit(709);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_710, matches, true, true);
  if (count > 1) {
    fp.setBit(710);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_711, matches, true, true);
  if (count > 1) {
    fp.setBit(711);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_712, matches, true, true);
  if (count > 1) {
    fp.setBit(712);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_713, matches, true, true);
  if (count > 1) {
    fp.setBit(713);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_714, matches, true, true);
  if (count > 1) {
    fp.setBit(714);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_715, matches, true, true);
  if (count > 1) {
    fp.setBit(715);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_716, matches, true, true);
  if (count > 1) {
    fp.setBit(716);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_717, matches, true, true);
  if (count > 1) {
    fp.setBit(717);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_718, matches, true, true);
  if (count > 1) {
    fp.setBit(718);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_719, matches, true, true);
  if (count > 1) {
    fp.setBit(719);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_720, matches, true, true);
  if (count > 1) {
    fp.setBit(720);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_721, matches, true, true);
  if (count > 1) {
    fp.setBit(721);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_722, matches, true, true);
  if (count > 1) {
    fp.setBit(722);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_723, matches, true, true);
  if (count > 1) {
    fp.setBit(723);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_724, matches, true, true);
  if (count > 1) {
    fp.setBit(724);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_725, matches, true, true);
  if (count > 1) {
    fp.setBit(725);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_726, matches, true, true);
  if (count > 1) {
    fp.setBit(726);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_727, matches, true, true);
  if (count > 1) {
    fp.setBit(727);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_728, matches, true, true);
  if (count > 1) {
    fp.setBit(728);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_729, matches, true, true);
  if (count > 1) {
    fp.setBit(729);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_730, matches, true, true);
  if (count > 1) {
    fp.setBit(730);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_731, matches, true, true);
  if (count > 1) {
    fp.setBit(731);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_732, matches, true, true);
  if (count > 1) {
    fp.setBit(732);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_733, matches, true, true);
  if (count > 1) {
    fp.setBit(733);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_734, matches, true, true);
  if (count > 1) {
    fp.setBit(734);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_735, matches, true, true);
  if (count > 1) {
    fp.setBit(735);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_736, matches, true, true);
  if (count > 1) {
    fp.setBit(736);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_737, matches, true, true);
  if (count > 1) {
    fp.setBit(737);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_738, matches, true, true);
  if (count > 1) {
    fp.setBit(738);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_739, matches, true, true);
  if (count > 1) {
    fp.setBit(739);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_740, matches, true, true);
  if (count > 1) {
    fp.setBit(740);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_741, matches, true, true);
  if (count > 1) {
    fp.setBit(741);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_742, matches, true, true);
  if (count > 1) {
    fp.setBit(742);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_743, matches, true, true);
  if (count > 1) {
    fp.setBit(743);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_744, matches, true, true);
  if (count > 1) {
    fp.setBit(744);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_745, matches, true, true);
  if (count > 1) {
    fp.setBit(745);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_746, matches, true, true);
  if (count > 1) {
    fp.setBit(746);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_747, matches, true, true);
  if (count > 1) {
    fp.setBit(747);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_748, matches, true, true);
  if (count > 1) {
    fp.setBit(748);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_749, matches, true, true);
  if (count > 1) {
    fp.setBit(749);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_750, matches, true, true);
  if (count > 1) {
    fp.setBit(750);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_751, matches, true, true);
  if (count > 1) {
    fp.setBit(751);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_752, matches, true, true);
  if (count > 1) {
    fp.setBit(752);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_753, matches, true, true);
  if (count > 1) {
    fp.setBit(753);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_754, matches, true, true);
  if (count > 1) {
    fp.setBit(754);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_755, matches, true, true);
  if (count > 1) {
    fp.setBit(755);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_756, matches, true, true);
  if (count > 1) {
    fp.setBit(756);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_757, matches, true, true);
  if (count > 1) {
    fp.setBit(757);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_758, matches, true, true);
  if (count > 1) {
    fp.setBit(758);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_759, matches, true, true);
  if (count > 1) {
    fp.setBit(759);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_760, matches, true, true);
  if (count > 1) {
    fp.setBit(760);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_761, matches, true, true);
  if (count > 1) {
    fp.setBit(761);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_762, matches, true, true);
  if (count > 1) {
    fp.setBit(762);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_763, matches, true, true);
  if (count > 1) {
    fp.setBit(763);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_764, matches, true, true);
  if (count > 1) {
    fp.setBit(764);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_765, matches, true, true);
  if (count > 1) {
    fp.setBit(765);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_766, matches, true, true);
  if (count > 1) {
    fp.setBit(766);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_767, matches, true, true);
  if (count > 1) {
    fp.setBit(767);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_768, matches, true, true);
  if (count > 1) {
    fp.setBit(768);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_769, matches, true, true);
  if (count > 1) {
    fp.setBit(769);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_770, matches, true, true);
  if (count > 1) {
    fp.setBit(770);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_771, matches, true, true);
  if (count > 1) {
    fp.setBit(771);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_772, matches, true, true);
  if (count > 1) {
    fp.setBit(772);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_773, matches, true, true);
  if (count > 1) {
    fp.setBit(773);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_774, matches, true, true);
  if (count > 1) {
    fp.setBit(774);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_775, matches, true, true);
  if (count > 1) {
    fp.setBit(775);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_776, matches, true, true);
  if (count > 1) {
    fp.setBit(776);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_777, matches, true, true);
  if (count > 1) {
    fp.setBit(777);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_778, matches, true, true);
  if (count > 1) {
    fp.setBit(778);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_779, matches, true, true);
  if (count > 1) {
    fp.setBit(779);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_780, matches, true, true);
  if (count > 1) {
    fp.setBit(780);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_781, matches, true, true);
  if (count > 1) {
    fp.setBit(781);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_782, matches, true, true);
  if (count > 1) {
    fp.setBit(782);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_783, matches, true, true);
  if (count > 1) {
    fp.setBit(783);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_784, matches, true, true);
  if (count > 1) {
    fp.setBit(784);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_785, matches, true, true);
  if (count > 1) {
    fp.setBit(785);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_786, matches, true, true);
  if (count > 1) {
    fp.setBit(786);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_787, matches, true, true);
  if (count > 1) {
    fp.setBit(787);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_788, matches, true, true);
  if (count > 1) {
    fp.setBit(788);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_789, matches, true, true);
  if (count > 1) {
    fp.setBit(789);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_790, matches, true, true);
  if (count > 1) {
    fp.setBit(790);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_791, matches, true, true);
  if (count > 1) {
    fp.setBit(791);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_792, matches, true, true);
  if (count > 1) {
    fp.setBit(792);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_793, matches, true, true);
  if (count > 1) {
    fp.setBit(793);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_794, matches, true, true);
  if (count > 1) {
    fp.setBit(794);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_795, matches, true, true);
  if (count > 1) {
    fp.setBit(795);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_796, matches, true, true);
  if (count > 1) {
    fp.setBit(796);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_797, matches, true, true);
  if (count > 1) {
    fp.setBit(797);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_798, matches, true, true);
  if (count > 1) {
    fp.setBit(798);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_799, matches, true, true);
  if (count > 1) {
    fp.setBit(799);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_800, matches, true, true);
  if (count > 1) {
    fp.setBit(800);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_801, matches, true, true);
  if (count > 1) {
    fp.setBit(801);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_802, matches, true, true);
  if (count > 1) {
    fp.setBit(802);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_803, matches, true, true);
  if (count > 1) {
    fp.setBit(803);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_804, matches, true, true);
  if (count > 1) {
    fp.setBit(804);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_805, matches, true, true);
  if (count > 1) {
    fp.setBit(805);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_806, matches, true, true);
  if (count > 1) {
    fp.setBit(806);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_807, matches, true, true);
  if (count > 1) {
    fp.setBit(807);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_808, matches, true, true);
  if (count > 1) {
    fp.setBit(808);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_809, matches, true, true);
  if (count > 1) {
    fp.setBit(809);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_810, matches, true, true);
  if (count > 1) {
    fp.setBit(810);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_811, matches, true, true);
  if (count > 1) {
    fp.setBit(811);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_812, matches, true, true);
  if (count > 1) {
    fp.setBit(812);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_813, matches, true, true);
  if (count > 1) {
    fp.setBit(813);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_814, matches, true, true);
  if (count > 1) {
    fp.setBit(814);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_815, matches, true, true);
  if (count > 1) {
    fp.setBit(815);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_816, matches, true, true);
  if (count > 1) {
    fp.setBit(816);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_817, matches, true, true);
  if (count > 1) {
    fp.setBit(817);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_818, matches, true, true);
  if (count > 1) {
    fp.setBit(818);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_819, matches, true, true);
  if (count > 1) {
    fp.setBit(819);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_820, matches, true, true);
  if (count > 1) {
    fp.setBit(820);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_821, matches, true, true);
  if (count > 1) {
    fp.setBit(821);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_822, matches, true, true);
  if (count > 1) {
    fp.setBit(822);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_823, matches, true, true);
  if (count > 1) {
    fp.setBit(823);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_824, matches, true, true);
  if (count > 1) {
    fp.setBit(824);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_825, matches, true, true);
  if (count > 1) {
    fp.setBit(825);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_826, matches, true, true);
  if (count > 1) {
    fp.setBit(826);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_827, matches, true, true);
  if (count > 1) {
    fp.setBit(827);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_828, matches, true, true);
  if (count > 1) {
    fp.setBit(828);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_829, matches, true, true);
  if (count > 1) {
    fp.setBit(829);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_830, matches, true, true);
  if (count > 1) {
    fp.setBit(830);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_831, matches, true, true);
  if (count > 1) {
    fp.setBit(831);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_832, matches, true, true);
  if (count > 1) {
    fp.setBit(832);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_833, matches, true, true);
  if (count > 1) {
    fp.setBit(833);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_834, matches, true, true);
  if (count > 1) {
    fp.setBit(834);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_835, matches, true, true);
  if (count > 1) {
    fp.setBit(835);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_836, matches, true, true);
  if (count > 1) {
    fp.setBit(836);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_837, matches, true, true);
  if (count > 1) {
    fp.setBit(837);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_838, matches, true, true);
  if (count > 1) {
    fp.setBit(838);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_839, matches, true, true);
  if (count > 1) {
    fp.setBit(839);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_840, matches, true, true);
  if (count > 1) {
    fp.setBit(840);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_841, matches, true, true);
  if (count > 1) {
    fp.setBit(841);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_842, matches, true, true);
  if (count > 1) {
    fp.setBit(842);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_843, matches, true, true);
  if (count > 1) {
    fp.setBit(843);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_844, matches, true, true);
  if (count > 1) {
    fp.setBit(844);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_845, matches, true, true);
  if (count > 1) {
    fp.setBit(845);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_846, matches, true, true);
  if (count > 1) {
    fp.setBit(846);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_847, matches, true, true);
  if (count > 1) {
    fp.setBit(847);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_848, matches, true, true);
  if (count > 1) {
    fp.setBit(848);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_849, matches, true, true);
  if (count > 1) {
    fp.setBit(849);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_850, matches, true, true);
  if (count > 1) {
    fp.setBit(850);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_851, matches, true, true);
  if (count > 1) {
    fp.setBit(851);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_852, matches, true, true);
  if (count > 1) {
    fp.setBit(852);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_853, matches, true, true);
  if (count > 1) {
    fp.setBit(853);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_854, matches, true, true);
  if (count > 1) {
    fp.setBit(854);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_855, matches, true, true);
  if (count > 1) {
    fp.setBit(855);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_856, matches, true, true);
  if (count > 1) {
    fp.setBit(856);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_857, matches, true, true);
  if (count > 1) {
    fp.setBit(857);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_858, matches, true, true);
  if (count > 1) {
    fp.setBit(858);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_859, matches, true, true);
  if (count > 1) {
    fp.setBit(859);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_860, matches, true, true);
  if (count > 1) {
    fp.setBit(860);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_861, matches, true, true);
  if (count > 1) {
    fp.setBit(861);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_862, matches, true, true);
  if (count > 1) {
    fp.setBit(862);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_863, matches, true, true);
  if (count > 1) {
    fp.setBit(863);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_864, matches, true, true);
  if (count > 1) {
    fp.setBit(864);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_865, matches, true, true);
  if (count > 1) {
    fp.setBit(865);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_866, matches, true, true);
  if (count > 1) {
    fp.setBit(866);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_867, matches, true, true);
  if (count > 1) {
    fp.setBit(867);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_868, matches, true, true);
  if (count > 1) {
    fp.setBit(868);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_869, matches, true, true);
  if (count > 1) {
    fp.setBit(869);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_870, matches, true, true);
  if (count > 1) {
    fp.setBit(870);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_871, matches, true, true);
  if (count > 1) {
    fp.setBit(871);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_872, matches, true, true);
  if (count > 1) {
    fp.setBit(872);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_873, matches, true, true);
  if (count > 1) {
    fp.setBit(873);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_874, matches, true, true);
  if (count > 1) {
    fp.setBit(874);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_875, matches, true, true);
  if (count > 1) {
    fp.setBit(875);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_876, matches, true, true);
  if (count > 1) {
    fp.setBit(876);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_877, matches, true, true);
  if (count > 1) {
    fp.setBit(877);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_878, matches, true, true);
  if (count > 1) {
    fp.setBit(878);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_879, matches, true, true);
  if (count > 1) {
    fp.setBit(879);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_880, matches, true, true);
  if (count > 1) {
    fp.setBit(880);
  }
  count = RDKit::SubstructMatch(mol, *pats.bit_881, matches, true, true);
  if (count > 1) {
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
