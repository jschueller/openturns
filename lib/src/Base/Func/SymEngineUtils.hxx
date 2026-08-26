//                                               -*- C++ -*-
/**
 *  @brief Utilities for SymEngine symbolic differentiation integration.
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef OPENTURNS_SYMENGINEUTILS_HXX
#define OPENTURNS_SYMENGINEUTILS_HXX

#include "openturns/OTconfig.hxx"

#ifdef OPENTURNS_HAVE_SYMENGINE

#include <symengine/constants.h>
#include <symengine/symbol.h>
#include <string>
#include <set>

BEGIN_NAMESPACE_OPENTURNS

/**
 * Return the set of SymEngine built-in constant names that can
 * shadow user variables if used as variable names.
 *
 * SymEngine's parser natively recognizes these as mathematical
 * constants. If a user formula has variables with these names,
 * the parser would treat them as constants instead of variables,
 * producing wrong derivatives.
 */
inline std::set<std::string> GetSymEngineBuiltinNames()
{
  return {"E", "pi", "I", "EulerGamma", "Catalan", "GoldenRatio"};
}

/**
 * Rename variables that clash with SymEngine built-in constants
 * to avoid the parser treating them as constants.
 *
 * For each variable name that matches a SymEngine built-in,
 * the formula is rewritten with a prefixed name (e.g. "E" -> "_OT_E_").
 * The renaming map is returned so the caller can reverse it.
 */
inline std::map<std::string, std::string> RenameClashingVariables(
    String & formula,
    const Description & variableNames)
{
  const std::set<std::string> builtins = GetSymEngineBuiltinNames();
  std::map<std::string, std::string> renameMap;

  for (UnsignedInteger i = 0; i < variableNames.getSize(); ++i)
  {
    const std::string varName = variableNames[i];
    if (builtins.find(varName) != builtins.end())
    {
      const std::string safeName = "_OT_" + varName + "_";
      renameMap[varName] = safeName;

      // Replace standalone occurrences of varName in the formula
      std::string::size_type pos = 0;
      while (pos < formula.size())
      {
        pos = formula.find(varName, pos);
        if (pos == std::string::npos) break;

        const bool atStart = (pos == 0);
        const bool prevIsIdent = !atStart && (std::isalnum(formula[pos - 1]) || formula[pos - 1] == '_');
        const bool nextIsIdent = (pos + varName.size() < formula.size()) && (std::isalnum(formula[pos + varName.size()]) || formula[pos + varName.size()] == '_');

        if (!prevIsIdent && !nextIsIdent)
        {
          formula.replace(pos, varName.size(), safeName);
          pos += safeName.size();
        }
        else
        {
          pos += varName.size();
        }
      }
    }
  }
  return renameMap;
}

/**
 * Reverse the variable renaming applied by RenameClashingVariables.
 *
 * After SymEngine differentiation and conversion to ExprTk format,
 * this restores original variable names.
 */
inline void RestoreVariableNames(
    String & formula,
    const std::map<std::string, std::string> & renameMap)
{
  for (const auto & entry : renameMap)
  {
    const std::string safeName = entry.second;
    const std::string originalName = entry.first;
    std::string::size_type pos = 0;
    while ((pos = formula.find(safeName, pos)) != std::string::npos)
    {
      formula.replace(pos, safeName.size(), originalName);
      pos += originalName.size();
    }
  }
}

/**
 * Replace ExprTk function names that SymEngine cannot differentiate
 * with mathematically equivalent expressions that SymEngine can handle.
 *
 * SymEngine lacks analytical derivatives for: log2, log10, log1p, expm1,
 * cbrt, lngamma. This function rewrites these to equivalent forms using
 * functions SymEngine knows (log, exp, pow, gamma).
 */
inline void ReplaceUnsupportedFunctions(String & formula)
{
  // Helper: find the matching close-paren for the open-paren at pos.
  // pos must point to '('.
  auto findMatchingParen = [&formula](std::string::size_type openPos) -> std::string::size_type
  {
    UnsignedInteger depth = 0;
    for (std::string::size_type i = openPos; i < formula.size(); ++i)
    {
      if (formula[i] == '(') ++depth;
      else if (formula[i] == ')')
      {
        --depth;
        if (depth == 0) return i;
      }
    }
    return std::string::npos;
  };

  // Each rule: function name to find, replacement for the function-name+open-paren,
  // and text appended after the matching close-paren.
  struct Rule {
    std::string funcName;   // e.g. "log2"
    std::string replacement; // replaces "funcName(" in the formula
    std::string afterParen; // appended after the matching ')'
  };

  const std::vector<Rule> rules = {
    {"log10",   "log(10)^(-1)*log(",  ""},
    {"log1p",   "log(1+",             ""},
    {"log2",    "log(2)^(-1)*log(",   ""},
    {"expm1",   "(exp(",              "-1)"},
    {"cbrt",    "(",                  "^(1.0/3.0)"},
  };

  for (const auto & rule : rules)
  {
    const std::string pattern = rule.funcName + "(";
    std::string::size_type searchPos = 0;
    while (true)
    {
      std::string::size_type pos = formula.find(pattern, searchPos);
      if (pos == std::string::npos) break;

      const std::string::size_type openParen = pos + rule.funcName.size();
      const std::string::size_type closeParen = findMatchingParen(openParen);
      if (closeParen == std::string::npos) break;

      const std::ptrdiff_t lenDiff = static_cast<std::ptrdiff_t>(rule.replacement.size()) - static_cast<std::ptrdiff_t>(pattern.size());
      formula.replace(pos, pattern.size(), rule.replacement);
      const std::string::size_type newCloseParen = closeParen + lenDiff;
      if (!rule.afterParen.empty())
        formula.insert(newCloseParen + 1, rule.afterParen);
      searchPos = newCloseParen + 1 + rule.afterParen.size();
    }
  }
}

/**
 * Return a constants map that registers ExprTk constant names
 * as SymEngine mathematical constants.
 *
 * ExprTk registers "e_" and "pi_" as constants. ev3 recognizes both
 * "e_"/"_e" and "pi_"/"_pi". SymEngine's parser only recognizes "E"
 * and "pi" natively. This map bridges the gap so that formulas using
 * ExprTk/ev3 constant names can be parsed by SymEngine.
 */
inline std::map<const std::string, const SymEngine::RCP<const SymEngine::Basic>> GetExprTkConstantsMap()
{
  std::map<const std::string, const SymEngine::RCP<const SymEngine::Basic>> constants;
  const SymEngine::RCP<const SymEngine::Basic> euler = SymEngine::E;
  const SymEngine::RCP<const SymEngine::Basic> pi = SymEngine::pi;
  constants.insert({"e_", euler});
  constants.insert({"pi_", pi});
  return constants;
}

/**
 * Convert a formula string from SymEngine's __str__() format
 * to ExprTk-compatible syntax.
 *
 * SymEngine differences vs ExprTk:
 *   - Power operator: ** vs ^
 *   - ceiling(x) vs ceil(x)
 *   - truncate(x) vs trunc(x)
 *   - Euler's number: E vs e_
 *   - Pi: pi vs pi_
 *
 * ev3 differences vs ExprTk (for reference):
 *   - ev3 outputs ^ for power (compatible with ExprTk)
 *   - ev3 outputs ceil, trunc (compatible with ExprTk)
 *   - ev3 outputs numeric values for constants (compatible with ExprTk)
 */
inline String ConvertSymEngineToExprTk(const String & symengineFormula)
{
  String result = symengineFormula;

  // Replace ** with ^ for power operator
  {
    std::string::size_type pos = result.find("**");
    while (pos != std::string::npos)
    {
      result.replace(pos, 2, "^");
      pos = result.find("**", pos);
    }
  }

  // Replace ceiling( with ceil(
  {
    const std::string from = "ceiling(";
    const std::string to = "ceil(";
    std::string::size_type pos = result.find(from);
    while (pos != std::string::npos)
    {
      result.replace(pos, from.size(), to);
      pos = result.find(from, pos + to.size());
    }
  }

  // Replace truncate( with trunc(
  {
    const std::string from = "truncate(";
    const std::string to = "trunc(";
    std::string::size_type pos = result.find(from);
    while (pos != std::string::npos)
    {
      result.replace(pos, from.size(), to);
      pos = result.find(from, pos + to.size());
    }
  }

  // Replace standalone E (Euler's number) with e_
  // Must not replace E inside function names (erf, erfc, exp, etc.)
  {
    std::string::size_type pos = 0;
    while (pos < result.size())
    {
      pos = result.find('E', pos);
      if (pos == std::string::npos) break;

      const bool atStart = (pos == 0);
      const bool prevIsIdent = !atStart && (std::isalnum(result[pos - 1]) || result[pos - 1] == '_');
      const bool nextIsIdent = (pos + 1 < result.size()) && (std::isalnum(result[pos + 1]) || result[pos + 1] == '_');

      if (!prevIsIdent && !nextIsIdent)
      {
        result.replace(pos, 1, "e_");
      }
      else
      {
        ++pos;
      }
    }
  }

  // Replace standalone pi with pi_
  // SymEngine outputs "pi" for the constant, but ExprTk uses "pi_"
  {
    const std::string from = "pi";
    const std::string to = "pi_";
    std::string::size_type pos = 0;
    while (pos < result.size())
    {
      pos = result.find(from, pos);
      if (pos == std::string::npos) break;

      const bool atStart = (pos == 0);
      const bool prevIsIdent = !atStart && (std::isalnum(result[pos - 1]) || result[pos - 1] == '_');
      const bool nextIsIdent = (pos + from.size() < result.size()) && (std::isalnum(result[pos + from.size()]) || result[pos + from.size()] == '_');

      if (!prevIsIdent && !nextIsIdent)
      {
        result.replace(pos, from.size(), to);
        pos += to.size();
      }
      else
      {
        pos += from.size();
      }
    }
  }

  return result;
}

/**
 * Check if a formula string produced by SymEngine contains unevaluated
 * Derivative or Subs expressions, which cannot be evaluated by ExprTk.
 *
 * This happens when SymEngine cannot compute an analytical derivative
 * for a function (e.g. abs, cbrt, expm1, min, max) and leaves the
 * derivative unevaluated as Derivative(func, var) or
 * Subs(Derivative(...), ...).
 */
inline Bool HasUnevaluatedDerivatives(const String & formula)
{
  return formula.find("Derivative(") != std::string::npos
      || formula.find("Subs(") != std::string::npos;
}

END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_HAVE_SYMENGINE */

#endif /* OPENTURNS_SYMENGINEUTILS_HXX */
