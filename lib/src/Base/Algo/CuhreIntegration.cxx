//                                               -*- C++ -*-
/**
 *  @brief Implement the Cuhre adaptive cubature method for multidimensional integration
 *
 *  Copyright 2005-2027 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "openturns/CuhreIntegration.hxx"
#include "openturns/Exception.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/ResourceMap.hxx"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

BEGIN_NAMESPACE_OPENTURNS

/**
 * Native serial implementation of the Cuhre algorithm: deterministic adaptive
 * subdivision of the unit hypercube driven by Genz-Malik cubature rules.
 * Only the serial code path is implemented: no parallel evaluation of the
 * integrand, no state file. Each subregion is evaluated with a single
 * vectorized Function call.
 *
 * References:
 * - Thomas Hahn, "Cuba: a library for multidimensional numerical
 *   integration", Computer Physics Communications, volume 168, issue 2,
 *   pages 78-95, 2005. The Cuhre driver reimplemented here (adaptive
 *   subdivision, weighted combination of successive iterations, chi-square
 *   probability of the error) follows the Cuhre routine of the Cuba library.
 * - J. Berntsen, T. O. Espelid and A. Genz, "An adaptive algorithm for the
 *   approximate calculation of multiple integrals", ACM Transactions on
 *   Mathematical Software, volume 17, issue 4, pages 437-451, 1991
 *   (algorithm DCUHRE). The degree-7, 9, 11 and 13 cubature rules tabulated
 *   below are the DCUHRE rules as distributed with the Cuba library.
 */
namespace
{

/* Number of rules combined in each cubature rule: the approximation rule
 * plus the null rules used for the error estimate. */
const UnsignedInteger cuhreRuleNumber = 5;

/* Tiny value guarding divisions by a vanishing error, matches Cuba NOTZERO. */
const Scalar cuhreNotZero = std::ldexp(1.0, -104);

/* One fully symmetric point set of a cubature rule: pointNumber_ points are
 * generated from generator_ by central symmetry and permutation, sharing
 * one weight per rule. scale_ and norm_ describe, for each null rule, the
 * null combination of largest normalized value (see sampleCuhreRegion). */
struct CuhreSet
{
  UnsignedInteger pointNumber_ = 0;
  Scalar weight_[cuhreRuleNumber];
  Scalar scale_[cuhreRuleNumber];
  Scalar norm_[cuhreRuleNumber];
  Point generator_;
};

/* Cubature rule: collection of point sets plus error estimate coefficients. */
struct CuhreRule
{
  std::vector<CuhreSet> sets_;
  Scalar errorCoefficient_[3];
  UnsignedInteger pointNumber_ = 0;
};

/* Subregion of the unit hypercube with its local integral estimates. */
struct CuhreRegion
{
  UnsignedInteger depth_ = 0;
  Point lower_;
  Point upper_;
  Point average_;
  Point error_;
  Indices bisectDimension_;
};

/* Integral estimates of one component accumulated over all subregions. */
struct CuhreTotals
{
  Scalar average_ = 0.0;
  Scalar error_ = 0.0;
  Scalar lastAverage_ = 0.0;
  Scalar lastError_ = 0.0;
  Scalar weightSum_ = 0.0;
  Scalar averageSum_ = 0.0;
  Scalar guess_ = 0.0;
  Scalar chiSum_ = 0.0;
  Scalar chiSquareSum_ = 0.0;
  Scalar chiSquare_ = 0.0;
};

void appendCuhreSet(CuhreRule & rule,
                    const UnsignedInteger pointNumber,
                    const Scalar weights[cuhreRuleNumber],
                    const Point & generator)
{
  CuhreSet set;
  set.pointNumber_ = pointNumber;
  std::copy(weights, weights + cuhreRuleNumber, set.weight_);
  set.generator_ = generator;
  rule.sets_.push_back(set);
  rule.pointNumber_ += pointNumber;
}

/* Precompute the scale/norm factors of the null-rule error search. */
void finalizeCuhreRule(CuhreRule & rule)
{
  for (UnsignedInteger r = 1; r < cuhreRuleNumber - 1; ++r)
  {
    for (CuhreSet & set : rule.sets_)
    {
      const Scalar scale = (set.weight_[r] == 0.0) ? 100.0 : -set.weight_[r + 1] / set.weight_[r];
      Scalar sum = 0.0;
      for (const CuhreSet & other : rule.sets_)
        sum += 1.0 * other.pointNumber_ * std::fabs(other.weight_[r + 1] + scale * other.weight_[r]);
      set.scale_[r] = scale;
      set.norm_[r] = 1.0 / sum;
    }
  }
}

/* Degree-13 rule, 2D only. */
CuhreRule buildCuhreRule13()
{
  static const Scalar weights[][cuhreRuleNumber] =
  {
    { 0.00844923090033615,      0.3213775489050763,     0.3372900883288987,
      -0.8264123822525677,      0.6539094339575232 },
    { 0.023771474018994404,    -0.1767341636743844,    -0.1644903060344491,
      0.306583861409436,       -0.2041614154422632},
    { 0.02940016170142405,      0.07347600537466073,     0.07707849911634623,
      0.002389292538329435,    -0.174698151579499 },
    { 0.006644436465817374,    -0.03638022004364754,    -0.03804478358506311,
      -0.1343024157997222,     0.03937939671417803 },
    { 0.0042536044255016,       0.021252979220987123,    0.02223559940380806,
      0.08833366840533902,     0.006974520545933992 },
    { 0,                        0.1460984204026913,      0.1480693879765931,
      0,                        0 },
    { 0.0040664827465935255,    0.017476132861520992,    4.467143702185815e-6,
      0.0009786283074168292,   0.0066677021717782585 },
    { 0.03362231646315497,      0.1444954045641582,      0.150894476707413,
      -0.1319227889147519,     0.05512960621544304 },
    { 0.033200804136503725,     0.0001307687976001325,   3.6472001075162155e-5,
      0.00799001220015063,     0.05443846381278608 },
    { 0.014093686924979677,     0.0005380992313941161,   0.000577719899901388,
      0.0033917470797606257,   0.02310903863953934 },
    { 0.000977069770327625,     0.0001042259576889814,   0.0001041757313688177,
      0.0022949157182832643,   0.01506937747477189 },
    { 0.007531996943580376,    -0.001401152865045733,   -0.001452822267047819,
      -0.01358584986119197,    -0.060570216489018905 },
    { 0.02577183086722915,      0.008041788181514763,    0.008338339968783704,
      0.04025866859057809,     0.04225737654686337},
    { 0.015625,                -0.1420416552759383,     -0.147279632923196,
      0.003760268580063992,    0.02561989142113099 }
  };
  static const Scalar generators[] =
  {
    0.12585646717265545,      0.3506966822267133,
    0.4795480315800981,       0.4978005239276064,
    0.25,                     0.07972723291487795,
    0.1904495567970094,       0.3291384627633596,
    0.43807365825146577,      0.499121592026599,
    0.4895111329084231,       0.32461421628226944,
    0.43637106005656195,      0.1791307322940614,
    0.2833333333333333,       0.1038888888888889
  };
  const UnsignedInteger dimension = 2;
  CuhreRule rule;
  appendCuhreSet(rule, 1, weights[0], Point(dimension, 0.0));
  for (UnsignedInteger k = 0; k < 5; ++k)
  {
    Point generator(dimension, 0.0);
    generator[0] = generators[k];
    appendCuhreSet(rule, 2 * dimension, weights[1 + k], generator);
  }
  for (UnsignedInteger k = 0; k < 5; ++k)
  {
    Point generator(dimension, generators[5 + k]);
    appendCuhreSet(rule, 2 * dimension * (dimension - 1), weights[6 + k], generator);
  }
  for (UnsignedInteger k = 0; k < 3; ++k)
  {
    Point generator(dimension, 0.0);
    generator[0] = generators[10 + 2 * k];
    generator[1] = generators[11 + 2 * k];
    appendCuhreSet(rule, 4 * dimension * (dimension - 1), weights[11 + k], generator);
  }
  rule.errorCoefficient_[0] = 10.0;
  rule.errorCoefficient_[1] = 1.0;
  rule.errorCoefficient_[2] = 5.0;
  finalizeCuhreRule(rule);
  return rule;
}

/* Degree-11 rule, 3D only. */
CuhreRule buildCuhreRule11()
{
  static const Scalar weights[][cuhreRuleNumber] =
  {
    { 0.0009903847688882167,   1.715006248224684,       1.936014978949526,
      0.517082819560576,       2.05440450381852 },
    { 0.0084964717409851,     -0.3755893815889209,     -0.3673449403754268,
      0.01445269144914044,     0.013777599884901202 },
    { 0.00013587331735072814,  0.1488632145140549,      0.02929778657898176,
      -0.3601489663995932,    -0.576806291790441 },
    { 0.022982920777660364,   -0.2497046640620823,     -0.1151883520260315,
      0.3628307003418485,      0.03726835047700328 },
    { 0.004202649722286289,    0.1792501419135204,      0.05086658220872218,
      0.007148802650872729,    0.0068148789397772195 },
    { 0.0012671889041675774,   0.0034461267589738897,   0.04453911087786469,
      -0.09222852896022966,    0.057231697338518496 },
    { 0.0002109560854981544,  -0.005140483185555825,   -0.022878282571259,
      0.01719339732471725,    -0.044930187438112855 },
    { 0.016830857056410086,    0.006536017839876424,    0.02908926216345833,
      -0.102141653746035,      0.027292365738663484 },
    { 0.00021876823557504823, -0.00065134549392297,    -0.002898884350669207,
      -0.007504397861080493,   0.000354747395055699 },
    { 0.009690420479796819,   -0.006304672433547204,   -0.028059634133074954,
      0.01648362537726711,     0.01571366799739551 },
    { 0.030773311284628138,    0.01266959399788263,     0.05638741361145884,
      0.05234610158469334,     0.049900992192785674 },
    { 0.0084974310856038,     -0.005454241018647931,   -0.02427469611942451,
      0.014454323316130661,    0.0137791555266677 },
    { 0.0017749535291258914,   0.004826995274768427,    0.021483070341828422,
      0.003019236275367777,    0.0028782064230998723 }
  };
  static const Scalar generators[] =
  {
    0.095,                    0.25,
    0.375,                    0.4,
    0.4975,                   0.49936724991757,
    0.38968518428362114,      0.49998494965443835,
    0.3951318612385894,       0.22016983438253684,
    0.4774686911397297,       0.2189239229503431,
    0.4830546566815374,       0.2288552938881567
  };
  const UnsignedInteger dimension = 3;
  CuhreRule rule;
  appendCuhreSet(rule, 1, weights[0], Point(dimension, 0.0));
  for (UnsignedInteger k = 0; k < 5; ++k)
  {
    Point generator(dimension, 0.0);
    generator[0] = generators[k];
    appendCuhreSet(rule, 2 * dimension, weights[1 + k], generator);
  }
  for (UnsignedInteger k = 0; k < 2; ++k)
  {
    Point generator(dimension, 0.0);
    generator[0] = generators[5 + k];
    generator[1] = generators[5 + k];
    appendCuhreSet(rule, 2 * dimension * (dimension - 1), weights[6 + k], generator);
  }
  for (UnsignedInteger k = 0; k < 3; ++k)
  {
    Point generator(dimension, 0.0);
    generator[0] = generators[7 + k];
    generator[1] = generators[7 + k];
    generator[2] = generators[7 + k];
    appendCuhreSet(rule, 4 * dimension * (dimension - 1) * (dimension - 2) / 3, weights[8 + k], generator);
  }
  for (UnsignedInteger k = 0; k < 2; ++k)
  {
    Point generator(dimension, 0.0);
    generator[0] = generators[10 + 2 * k];
    generator[1] = generators[11 + 2 * k];
    generator[2] = generators[11 + 2 * k];
    if (k == 1) generator[1] = generators[10 + 2 * k];
    appendCuhreSet(rule, 4 * dimension * (dimension - 1) * (dimension - 2), weights[11 + k], generator);
  }
  rule.errorCoefficient_[0] = 4.0;
  rule.errorCoefficient_[1] = 0.5;
  rule.errorCoefficient_[2] = 3.0;
  finalizeCuhreRule(rule);
  return rule;
}

/* Degree-9 rule, any dimension (default rule). */
CuhreRule buildCuhreRule9(const UnsignedInteger dimension)
{
  if (dimension > 20) throw InvalidArgumentException(HERE) << "Cuhre degree-9 rule needs 2^dimension corner points, which is infeasible with dimension=" << dimension;
  static const Scalar weights[] =
  {
    -0.002361170967885511788400941242259231309691,
    0.1141539002385732526821323741697655347686,
    -0.6383392007670238909386026193674701393074,
    0.7484998850468520800423030047583803945205,
    -0.001432401703339912514196154599769007103671,
    0.05747150786448972594860897296200006759892,
    -0.1422510457143424323449521620935950679394,
    -0.06287502873828697998942424881040490136987,
    0.2545911332489590890011611142429070613156,
    -1.207328566678236261002219995185143356737,
    0.8956736576416067650809467826488567200939,
    -0.3647935698604914666100134551377381205297,
    0.003541756451678267682601411863388846964536,
    -0.07260936739589367960492815865074633743652,
    0.1055749162521899101218622863269817454540,
    0.002148602555009868771294231899653510655506,
    -0.03226856389295394999786630399875134318006,
    0.01063678399023121748083624225818915724455,
    0.01468910249614349017540783437728097691502,
    0.5113470834646759143109387357149329909126,
    0.4597644812080634464633352781605214342691,
    0.1823967849302457333050067275688690602649,
    -0.04508628929435784075980562738240804429658,
    0.2141588352435279340097929526588394300172,
    -0.02735154652654564472203690086290223507436,
    0.05494106704871123410060080562462135546101,
    0.1193759620257077529708962121565290178730,
    0.6508951939192025059314756320878023215278,
    0.1474493982943446016775696826942585013243,
    0.05769338449097348357291272856892627722165,
    0.03499962660214358382244159694487155861542,
    -1.386862771927828143599782668709014266770,
    -0.2386668732575008878964134721962088068396,
    0.01553241727660705326386197156586357005224,
    0.003532809960709087023561817517751309380604,
    0.09231719987444221619017126187763868745587,
    0.02254314464717892037990281369120402214829,
    0.01367577326327282236101845043145111753718,
    -0.3254475969596012529657378160439011607639,
    0.001770878225839133841300705931694423482268,
    0.001074301277504934385647115949826755327753,
    0.2515001149531479199576969952416196054795
  };
  static const Scalar generators[] =
  {
    0.4779536579022695061928604197171830064732,
    0.2030285873691198677998034402373279133258,
    0.4476273546261781288207704806530998539285,
    0.125,
    0.3430378987808781457001426145164678603407
  };
  const Scalar dim = 1.0 * dimension;
  const UnsignedInteger twoToDimension = 1u << dimension;
  CuhreRule rule;
  {
    const Scalar w[cuhreRuleNumber] =
    {
      dim * (dim * (dim * weights[0] + weights[1]) + weights[2]) + weights[3],
      dim * (dim * (dim * weights[4] + weights[5]) + weights[6]) - weights[7],
      dim * weights[8] - (dim * (dim * (dim * weights[4] + weights[5]) + weights[6]) - weights[7]),
      dim * (dim * weights[9] + weights[10]) - 1.0 + (dim * (dim * (dim * weights[0] + weights[1]) + weights[2]) + weights[3]),
      dim * weights[11] + 1.0 - (dim * (dim * (dim * weights[0] + weights[1]) + weights[2]) + weights[3])
    };
    appendCuhreSet(rule, 1, w, Point(dimension, 0.0));
  }
  {
    const Scalar w0 = dim * (dim * weights[12] + weights[13]) + weights[14];
    const Scalar w1 = dim * (dim * weights[15] + weights[16]) + weights[17];
    const Scalar w[cuhreRuleNumber] =
    {
      w0,
      w1,
      weights[18] - w1,
      dim * weights[19] + weights[20] + w0,
      weights[21] - w0
    };
    Point generator(dimension, 0.0);
    generator[0] = generators[0];
    appendCuhreSet(rule, 2 * dimension, w, generator);
  }
  {
    const Scalar w0 = dim * weights[22] + weights[23];
    const Scalar w1 = dim * weights[24] + weights[25];
    const Scalar w[cuhreRuleNumber] =
    {
      w0,
      w1,
      weights[26] - w1,
      dim * weights[27] + weights[28],
      -w0
    };
    Point generator(dimension, 0.0);
    generator[0] = generators[1];
    appendCuhreSet(rule, 2 * dimension, w, generator);
  }
  {
    const Scalar w[cuhreRuleNumber] = { weights[29], weights[30], -weights[29], weights[31], -weights[29] };
    Point generator(dimension, 0.0);
    generator[0] = generators[2];
    appendCuhreSet(rule, 2 * dimension, w, generator);
  }
  {
    const Scalar w[cuhreRuleNumber] = { 0.0, 0.0, weights[32], 0.0, 0.0 };
    Point generator(dimension, 0.0);
    generator[0] = generators[3];
    appendCuhreSet(rule, 2 * dimension, w, generator);
  }
  {
    const Scalar w0 = weights[33] - dim * weights[12];
    const Scalar w1 = weights[34] - dim * weights[15];
    const Scalar w[cuhreRuleNumber] = { w0, w1, -w1, weights[35] + w0, -w0 };
    Point generator(dimension, 0.0);
    generator[0] = generators[0];
    generator[1] = generators[0];
    appendCuhreSet(rule, 2 * dimension * (dimension - 1), w, generator);
  }
  {
    const Scalar w[cuhreRuleNumber] = { weights[36], weights[37], -weights[37], weights[38], -weights[36] };
    Point generator(dimension, 0.0);
    generator[0] = generators[0];
    generator[1] = generators[1];
    appendCuhreSet(rule, 4 * dimension * (dimension - 1), w, generator);
  }
  {
    const Scalar w[cuhreRuleNumber] = { weights[39], weights[40], -weights[40], weights[39], -weights[39] };
    Point generator(dimension, 0.0);
    generator[0] = generators[0];
    generator[1] = generators[0];
    generator[2] = generators[0];
    appendCuhreSet(rule, 4 * dimension * (dimension - 1) * (dimension - 2) / 3, w, generator);
  }
  {
    const Scalar w0 = weights[41] / (1.0 * twoToDimension);
    const Scalar w1 = weights[7] / (1.0 * twoToDimension);
    const Scalar w[cuhreRuleNumber] = { w0, w1, -w1, w0, -w0 };
    Point generator(dimension, generators[4]);
    appendCuhreSet(rule, twoToDimension, w, generator);
  }
  rule.errorCoefficient_[0] = 5.0;
  rule.errorCoefficient_[1] = 1.0;
  rule.errorCoefficient_[2] = 5.0;
  finalizeCuhreRule(rule);
  return rule;
}

/* Degree-7 rule, any dimension. */
CuhreRule buildCuhreRule7(const UnsignedInteger dimension)
{
  if (dimension > 20) throw InvalidArgumentException(HERE) << "Cuhre degree-7 rule needs 2^dimension corner points, which is infeasible with dimension=" << dimension;
  static const Scalar weights[] =
  {
    0.01941786667474838842844534313920462333850,
    -0.4038525770115018254611834753723880293161,
    0.6448566876746598222277360730193089551024,
    0.01177982690775806141012214458820955067854,
    -0.1804131874073360901182293138710989490609,
    -0.08878582808133504444306598174517276122439,
    0.05632864580828594137378124255408286479947,
    -0.009708933337374194214222671569602311669249,
    -0.9912917677958235813775106862002319060386,
    -0.1775716561626700888861319634903455224488,
    0.1235939803204323357183625846672135876752,
    0.07497814870203369068087999555157339703666,
    0.5548914705142355977605994477355651401434,
    0.08804124152269277122645182458858273865209,
    0.02111835845551338508329573367808085283304,
    -0.009930220323965333308685820460105538586058,
    -0.06410005328501090417895544042025034295870,
    0.03038172903822100765927778829870429682489,
    0.005889913453879030705061072294104775339268,
    -0.004854466668687097107111335784801155834624,
    0.3551433123253401777722639269806910448976
  };
  static const Scalar generators[] =
  {
    0.4779536579022695061928604197171830064732,
    0.2030285873691198677998034402373279133258,
    0.375,
    0.3430378987808781457001426145164678603407
  };
  const Scalar dim = 1.0 * dimension;
  const UnsignedInteger twoToDimension = 1u << dimension;
  CuhreRule rule;
  {
    const Scalar w0 = dim * (dim * weights[0] + weights[1]) + weights[2];
    const Scalar w1 = dim * (dim * weights[3] + weights[4]) - weights[5];
    const Scalar w[cuhreRuleNumber] =
    {
      w0,
      w1,
      dim * weights[6] - w1,
      dim * (dim * weights[7] + weights[8]) - weights[9],
      1.0 - w0
    };
    appendCuhreSet(rule, 1, w, Point(dimension, 0.0));
  }
  {
    const Scalar w[cuhreRuleNumber] = { weights[10], weights[11], -weights[10], weights[12], -weights[10] };
    Point generator(dimension, 0.0);
    generator[0] = generators[1];
    appendCuhreSet(rule, 2 * dimension, w, generator);
  }
  {
    const Scalar w0 = weights[13] - dim * weights[0];
    const Scalar w1 = weights[14] - dim * weights[3];
    const Scalar w[cuhreRuleNumber] =
    {
      w0,
      w1,
      weights[15] - w1,
      weights[16] - dim * weights[7],
      -w0
    };
    Point generator(dimension, 0.0);
    generator[0] = generators[0];
    appendCuhreSet(rule, 2 * dimension, w, generator);
  }
  {
    const Scalar w[cuhreRuleNumber] = { 0.0, 0.0, weights[17], 0.0, 0.0 };
    Point generator(dimension, 0.0);
    generator[0] = generators[2];
    appendCuhreSet(rule, 2 * dimension, w, generator);
  }
  {
    const Scalar w[cuhreRuleNumber] = { -weights[7], weights[18], -weights[18], weights[19], weights[7] };
    Point generator(dimension, 0.0);
    generator[0] = generators[0];
    generator[1] = generators[0];
    appendCuhreSet(rule, 2 * dimension * (dimension - 1), w, generator);
  }
  {
    const Scalar w0 = weights[20] / (1.0 * twoToDimension);
    const Scalar w1 = weights[5] / (1.0 * twoToDimension);
    const Scalar w[cuhreRuleNumber] = { w0, w1, -w1, weights[9] / (1.0 * twoToDimension), -w0 };
    Point generator(dimension, generators[3]);
    appendCuhreSet(rule, twoToDimension, w, generator);
  }
  rule.errorCoefficient_[0] = 5.0;
  rule.errorCoefficient_[1] = 1.0;
  rule.errorCoefficient_[2] = 5.0;
  finalizeCuhreRule(rule);
  return rule;
}

/* Select the cubature rule: key 13 (2D only), 11 (3D only), 9, 7, else default. */
CuhreRule buildCuhreRule(const UnsignedInteger dimension,
                         const UnsignedInteger key)
{
  if ((key == 13) && (dimension == 2)) return buildCuhreRule13();
  if ((key == 11) && (dimension == 3)) return buildCuhreRule11();
  if (key == 9) return buildCuhreRule9(dimension);
  if (key == 7) return buildCuhreRule7(dimension);
  if (dimension == 2) return buildCuhreRule13();
  if (dimension == 3) return buildCuhreRule11();
  return buildCuhreRule9(dimension);
}

/* Append to buffer the centrally symmetric sums of all the distinct
 * permutations of generator over the box [lower, upper]. The generator is
 * passed by value as it is permuted in place and restored on exit. */
void appendCuhreGeneratorPoints(const Point & lower,
                               const Point & upper,
                               std::vector<Scalar> generator,
                               std::vector<Scalar> & buffer)
{
  const UnsignedInteger dimension = lower.getDimension();
  for (;;)
  {
    for (UnsignedInteger dim = 0; dim < dimension; ++dim)
      buffer.push_back(0.5 * (lower[dim] + upper[dim]) + generator[dim] * (lower[dim] - upper[dim]));
    UnsignedInteger dim = 0;
    while (dim < dimension)
    {
      generator[dim] = -generator[dim];
      if (generator[dim] < 0.0) break;
      ++dim;
    }
    if (dim < dimension) continue;
    int signedDim = 0;
    Bool permuted = false;
    for (signedDim = 1; signedDim < static_cast<int>(dimension); ++signedDim)
    {
      const Scalar guard = generator[signedDim];
      if (generator[signedDim - 1] > guard)
      {
        int index = signedDim;
        int down = signedDim - 1;
        int other = signedDim;
        for (int i = 0; i < --other; ++i)
        {
          const Scalar tmp = generator[i];
          generator[i] = generator[other];
          generator[other] = tmp;
          if (tmp <= guard) --down;
          if (generator[i] > guard) index = i;
        }
        if (generator[down] <= guard) down = index;
        generator[signedDim] = generator[down];
        generator[down] = guard;
        permuted = true;
        break;
      }
    }
    if (permuted) continue;
    std::reverse(generator.begin(), generator.end());
    return;
  }
}

/* Evaluate the cubature rule over a subregion: single vectorized Function
 * call over all the rule points mapped onto the integration interval, scaled
 * by the interval volume. Returns the number of integrand evaluations. */
UnsignedInteger sampleCuhreRegion(const CuhreRule & rule,
                                 const Function & function,
                                 const Point & intervalLower,
                                 const Point & intervalUpper,
                                 const Scalar intervalVolume,
                                 CuhreRegion & region)
{
  const UnsignedInteger inputDimension = intervalLower.getDimension();
  const UnsignedInteger outputDimension = region.average_.getDimension();
  std::vector<Scalar> unitPoints;
  unitPoints.reserve(rule.pointNumber_ * inputDimension);
  for (const CuhreSet & set : rule.sets_)
  {
    if (set.pointNumber_ == 0) continue;
    const std::vector<Scalar> generator(set.generator_.begin(), set.generator_.end());
    appendCuhreGeneratorPoints(region.lower_, region.upper_, generator, unitPoints);
  }
  if (unitPoints.size() != rule.pointNumber_ * inputDimension)
    throw InternalException(HERE) << "Cuhre point generation mismatch: got " << unitPoints.size() << " coordinates, expected " << rule.pointNumber_ * inputDimension;
  Sample inputSample(rule.pointNumber_, inputDimension);
  for (UnsignedInteger p = 0; p < rule.pointNumber_; ++p)
    for (UnsignedInteger d = 0; d < inputDimension; ++d)
      inputSample(p, d) = intervalLower[d] + (intervalUpper[d] - intervalLower[d]) * unitPoints[p * inputDimension + d];
  Sample values(function(inputSample));
  if ((values.getSize() != rule.pointNumber_) || (values.getDimension() != outputDimension))
    throw InternalException(HERE) << "Cuhre integrand evaluation mismatch: got size=" << values.getSize() << " dimension=" << values.getDimension();
  for (UnsignedInteger p = 0; p < rule.pointNumber_; ++p)
    for (UnsignedInteger c = 0; c < outputDimension; ++c)
      values(p, c) *= intervalVolume;
  /* Volume of the subregion within the unit hypercube. */
  const Scalar regionVolume = std::ldexp(1.0, -static_cast<int>(region.depth_));
  /* Ratio of the generators driving the fourth-difference bisection choice. */
  const Scalar generatorRatio = rule.sets_[2].generator_[0] / rule.sets_[1].generator_[0];
  const Scalar ratio = generatorRatio * generatorRatio;
  /* Points 1 to 2*ndim hold the first single-generator set, the next
   * 2*ndim points the second one. */
  const UnsignedInteger pairOffset = 2 * inputDimension;
  const UnsignedInteger setOffset = 1 + pairOffset;
  for (UnsignedInteger c = 0; c < outputDimension; ++c)
  {
    const Scalar base = values(0, c) * 2.0 * (1.0 - ratio);
    UnsignedInteger bisectDim = 0;
    Scalar maxRange = 0.0;
    for (UnsignedInteger d = 0; d < inputDimension; ++d)
    {
      const Scalar range = region.upper_[d] - region.lower_[d];
      if (range > maxRange)
      {
        maxRange = range;
        bisectDim = d;
      }
    }
    Scalar maxDiff = 0.0;
    for (UnsignedInteger d = 0; d < inputDimension; ++d)
    {
      const Scalar fourthDiff = std::fabs(base + ratio * (values(1 + 2 * d, c) + values(2 + 2 * d, c)) - (values(setOffset + 2 * d, c) + values(setOffset + 1 + 2 * d, c)));
      if (fourthDiff > maxDiff)
      {
        maxDiff = fourthDiff;
        bisectDim = d;
      }
    }
    Scalar sums[cuhreRuleNumber] = {};
    UnsignedInteger offset = 0;
    for (const CuhreSet & set : rule.sets_)
    {
      for (UnsignedInteger k = 0; k < set.pointNumber_; ++k)
      {
        const Scalar value = values(offset + k, c);
        for (UnsignedInteger r = 0; r < cuhreRuleNumber; ++r)
          sums[r] += value * set.weight_[r];
      }
      offset += set.pointNumber_;
    }
    /* Search for the null rule of largest normalized value. */
    for (UnsignedInteger r = 1; r < cuhreRuleNumber - 1; ++r)
    {
      Scalar maxError = 0.0;
      for (const CuhreSet & set : rule.sets_)
        maxError = std::max(maxError, std::fabs(sums[r + 1] + set.scale_[r] * sums[r]) * set.norm_[r]);
      sums[r] = maxError;
    }
    region.average_[c] = regionVolume * sums[0];
    const Bool smallErrors = (rule.errorCoefficient_[0] * sums[1] <= sums[2]) && (rule.errorCoefficient_[0] * sums[2] <= sums[3]);
    region.error_[c] = regionVolume * (smallErrors ? rule.errorCoefficient_[1] * sums[1] : rule.errorCoefficient_[2] * std::max(sums[1], std::max(sums[2], sums[3])));
    region.bisectDimension_[c] = bisectDim;
  }
  return rule.pointNumber_;
}

/* Adaptive subdivision driver. Stops once every component meets its tolerance
 * and minimumEvaluations have been spent; throws once maximumEvaluations is
 * exceeded. useLastSample selects the weighting of successive iterations. */
Point integrateCuhre(const Function & function,
                     const Point & intervalLower,
                     const Point & intervalUpper,
                     const Scalar intervalVolume,
                     const UnsignedInteger key,
                     const Scalar relativeError,
                     const Scalar absoluteError,
                     const UnsignedInteger minimumEvaluations,
                     const UnsignedInteger maximumEvaluations,
                     const Bool useLastSample)
{
  const UnsignedInteger inputDimension = intervalLower.getDimension();
  const UnsignedInteger outputDimension = function.getOutputDimension();
  const CuhreRule rule(buildCuhreRule(inputDimension, key));
  const Scalar absoluteTolerance = std::max(absoluteError, cuhreNotZero);
  const UnsignedInteger minimumRequired = std::max(minimumEvaluations, rule.pointNumber_ + 1);
  std::vector<CuhreRegion> regions;
  CuhreRegion first;
  first.lower_ = Point(inputDimension, 0.0);
  first.upper_ = Point(inputDimension, 1.0);
  first.average_ = Point(outputDimension, 0.0);
  first.error_ = Point(outputDimension, 0.0);
  first.bisectDimension_ = Indices(outputDimension, 0);
  regions.push_back(first);
  std::vector<CuhreTotals> totals(outputDimension);
  UnsignedInteger evaluations = sampleCuhreRegion(rule, function, intervalLower, intervalUpper, intervalVolume, regions[0]);
  for (UnsignedInteger c = 0; c < outputDimension; ++c)
  {
    CuhreTotals & total = totals[c];
    total.average_ = total.lastAverage_ = total.guess_ = regions[0].average_[c];
    total.error_ = total.lastError_ = regions[0].error_[c];
    total.weightSum_ = 1.0 / std::max(total.error_ * total.error_, cuhreNotZero);
    total.averageSum_ = total.weightSum_ * total.average_;
  }
  UnsignedInteger regionNumber = 1;
  for (;;)
  {
    Scalar maxRatio = -std::numeric_limits<Scalar>::max();
    UnsignedInteger worstComponent = 0;
    for (UnsignedInteger c = 0; c < outputDimension; ++c)
    {
      const Scalar tolerance = std::max(relativeError * std::fabs(totals[c].average_), absoluteTolerance);
      const Scalar ratio = totals[c].error_ / tolerance;
      if (ratio > maxRatio)
      {
        maxRatio = ratio;
        worstComponent = c;
      }
    }
    if ((maxRatio <= 1.0) && (evaluations >= minimumRequired)) break;
    if (evaluations >= maximumEvaluations)
      throw InternalException(HERE) << "Cuhre integration failed to reach the requested accuracy within " << maximumEvaluations << " integrand evaluations on " << regionNumber << " subregions";
    UnsignedInteger worstRegion = 0;
    Scalar maxError = -std::numeric_limits<Scalar>::max();
    const UnsignedInteger regionCount = regions.size();
    for (UnsignedInteger i = 0; i < regionCount; ++i)
    {
      if (regions[i].error_[worstComponent] > maxError)
      {
        maxError = regions[i].error_[worstComponent];
        worstRegion = i;
      }
    }
    CuhreRegion left(regions[worstRegion]);
    const Point oldAverage(left.average_);
    const Point oldError(left.error_);
    const UnsignedInteger bisectDim = left.bisectDimension_[worstComponent];
    const Scalar middle = 0.5 * (left.lower_[bisectDim] + left.upper_[bisectDim]);
    CuhreRegion right(left);
    left.depth_ += 1;
    right.depth_ += 1;
    left.upper_[bisectDim] = middle;
    right.lower_[bisectDim] = middle;
    evaluations += sampleCuhreRegion(rule, function, intervalLower, intervalUpper, intervalVolume, left);
    evaluations += sampleCuhreRegion(rule, function, intervalLower, intervalUpper, intervalVolume, right);
    for (UnsignedInteger c = 0; c < outputDimension; ++c)
    {
      CuhreTotals & total = totals[c];
      Scalar difference = left.average_[c] + right.average_[c] - oldAverage[c];
      total.lastAverage_ += difference;
      difference = std::fabs(0.25 * difference);
      Scalar error = left.error_[c] + right.error_[c];
      if (error > 0.0)
      {
        const Scalar factor = 1.0 + 2.0 * difference / error;
        left.error_[c] *= factor;
        right.error_[c] *= factor;
      }
      left.error_[c] += difference;
      right.error_[c] += difference;
      total.lastError_ += left.error_[c] + right.error_[c] - oldError[c];
      Scalar weight = 1.0 / std::max(total.lastError_ * total.lastError_, cuhreNotZero);
      total.weightSum_ += weight;
      const Scalar inverseWeightSum = 1.0 / total.weightSum_;
      total.averageSum_ += weight * total.lastAverage_;
      const Scalar average = inverseWeightSum * total.averageSum_;
      weight *= total.lastAverage_ - total.guess_;
      total.chiSum_ += weight;
      total.chiSquareSum_ += weight * total.lastAverage_;
      total.chiSquare_ = total.chiSquareSum_ - average * total.chiSum_;
      if (useLastSample)
      {
        total.average_ = total.lastAverage_;
        total.error_ = total.lastError_;
      }
      else
      {
        total.average_ = average;
        total.error_ = std::sqrt(inverseWeightSum);
      }
    }
    /* Store the bisected subregions with their inflated errors, as the
     * largest-error search below reads them back. */
    regions[worstRegion] = left;
    regions.push_back(right);
    ++regionNumber;
  }
  Point integral(outputDimension);
  for (UnsignedInteger c = 0; c < outputDimension; ++c)
    integral[c] = totals[c].average_;
  return integral;
}

} // anonymous namespace

CLASSNAMEINIT(CuhreIntegration)

static const Factory<CuhreIntegration> Factory_CuhreIntegration;

/* Parameters constructor */
CuhreIntegration::CuhreIntegration(const String algorithmName)
  : IntegrationAlgorithmImplementation()
  , maximumRelativeError_(ResourceMap::GetAsScalar("CuhreIntegration-DefaultMaximumRelativeError"))
  , maximumAbsoluteError_(ResourceMap::GetAsScalar("CuhreIntegration-DefaultMaximumAbsoluteError"))
  , maximumCallsNumber_(ResourceMap::GetAsUnsignedInteger("CuhreIntegration-DefaultMaximumCallsNumber"))
{
  setAlgorithmName(algorithmName);
}

/* Virtual constructor */
CuhreIntegration * CuhreIntegration::clone() const
{
  return new CuhreIntegration(*this);
}

/* Compute an approximation of \int_{I}f(x)dx, where I
 * is an multidimensional interval and f a multidimensional function
 */
Point CuhreIntegration::integrate(const Function & function,
                                 const Interval & interval) const
{
  if (function.getInputDimension() != interval.getDimension())
    throw InvalidArgumentException(HERE) << "Error: interval dimension and input dimension have to match, here interval dimension=" << interval.getDimension() << " and input dimension=" << function.getInputDimension();
  if (algorithmName_ != "cuhre")
    throw NotYetImplementedException(HERE) << "Error: only the cuhre algorithm is supported, here algorithm=" << algorithmName_;
  const UnsignedInteger inputDimension = interval.getDimension();
  const UnsignedInteger outputDimension = function.getOutputDimension();
  if (inputDimension < 2)
    throw InvalidArgumentException(HERE) << "Cuhre minimum input dimension is 2, got " << inputDimension;
  if (inputDimension > 1024)
    throw InvalidArgumentException(HERE) << "Cuhre maximum input dimension is 1024, got " << inputDimension;
  if (outputDimension > 1024)
    throw InvalidArgumentException(HERE) << "Cuhre maximum output dimension is 1024, got " << outputDimension;
  const UnsignedInteger key = ResourceMap::GetAsUnsignedInteger("CuhreIntegration-key");
  const UnsignedInteger minimumEvaluations = ResourceMap::GetAsUnsignedInteger("CuhreIntegration-mineval");
  const UnsignedInteger flags = ResourceMap::GetAsUnsignedInteger("CuhreIntegration-flags");
  const Bool useLastSample = (flags & 4) != 0;
  return integrateCuhre(function, interval.getLowerBound(), interval.getUpperBound(), interval.getVolume(), key, maximumRelativeError_, maximumAbsoluteError_, minimumEvaluations, maximumCallsNumber_, useLastSample);
}

/* maximumRelativeError accessor */
Scalar CuhreIntegration::getMaximumRelativeError() const
{
  return maximumRelativeError_;
}

void CuhreIntegration::setMaximumRelativeError(const Scalar maximumRelativeError)
{
  if (!(maximumRelativeError > 0)) throw InvalidArgumentException(HERE) << "Error: maximumRelativeError must be positive, here maximumRelativeError=" << maximumRelativeError;
  maximumRelativeError_ = maximumRelativeError;
}

/* maximumAbsoluteError accessor */
Scalar CuhreIntegration::getMaximumAbsoluteError() const
{
  return maximumAbsoluteError_;
}

void CuhreIntegration::setMaximumAbsoluteError(const Scalar maximumAbsoluteError)
{
  if (!(maximumAbsoluteError > 0)) throw InvalidArgumentException(HERE) << "Error: maximumAbsoluteError must be positive, here maximumAbsoluteError=" << maximumAbsoluteError;
  maximumAbsoluteError_ = maximumAbsoluteError;
}

/* algorithmName accessor */
String CuhreIntegration::getAlgorithmName() const
{
  return algorithmName_;
}

void CuhreIntegration::setAlgorithmName(const String algorithmName)
{
  if (!GetAlgorithmNames().contains(algorithmName)) throw InvalidArgumentException(HERE) << "Error: the name of the integration routine should be \"cuhre\", here algorithm=" << algorithmName;
  algorithmName_ = algorithmName;
}

/* maximumCallsNumber accessor */
UnsignedInteger CuhreIntegration::getMaximumCallsNumber() const
{
  return maximumCallsNumber_;
}

void CuhreIntegration::setMaximumCallsNumber(const UnsignedInteger maximumCallsNumber)
{
  maximumCallsNumber_ = maximumCallsNumber;
}

/* Return the names of the available integration algorithms */
Description CuhreIntegration::GetAlgorithmNames()
{
  return {"cuhre"};
}

/* String converter */
String CuhreIntegration::__repr__() const
{
  OSS oss(true);
  oss << "class=" << CuhreIntegration::GetClassName()
      << ", maximumRelativeError=" << maximumRelativeError_
      << ", maximumAbsoluteError=" << maximumAbsoluteError_
      << ", algorithmName=" << algorithmName_
      << ", maximumCallsNumber=" << maximumCallsNumber_;
  return oss;
}

/* String converter */
String CuhreIntegration::__str__(const String & ) const
{
  OSS oss(false);
  oss << CuhreIntegration::GetClassName()
      << "(algorithmName=" << algorithmName_
      << ", maximumRelativeError=" << maximumRelativeError_
      << ", maximumAbsoluteError=" << maximumAbsoluteError_
      << ", maximumCallsNumber=" << maximumCallsNumber_
      << ")";
  return oss;
}

/* Method save() stores the object through the StorageManager */
void CuhreIntegration::save(Advocate & adv) const
{
  IntegrationAlgorithmImplementation::save(adv);
  adv.saveAttribute("maximumRelativeError_", maximumRelativeError_);
  adv.saveAttribute("maximumAbsoluteError_", maximumAbsoluteError_);
  adv.saveAttribute("algorithmName_", algorithmName_);
  adv.saveAttribute("maximumCallsNumber_", maximumCallsNumber_);
}

/* Method load() reloads the object from the StorageManager */
void CuhreIntegration::load(Advocate & adv)
{
  IntegrationAlgorithmImplementation::load(adv);
  adv.loadAttribute("maximumRelativeError_", maximumRelativeError_);
  adv.loadAttribute("maximumAbsoluteError_", maximumAbsoluteError_);
  adv.loadAttribute("algorithmName_", algorithmName_);
  adv.loadAttribute("maximumCallsNumber_", maximumCallsNumber_);
}

END_NAMESPACE_OPENTURNS
